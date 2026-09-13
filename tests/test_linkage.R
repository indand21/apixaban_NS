# test_linkage.R
# Unit tests for PBPK-QSP linkage

source("R/00_packages.R")
source("R/01_utils.R")
source("R/10_pbpk_parameters.R")
source("R/11_pbpk_odes.R")
source("R/12_pbpk_simulate.R")
source("R/20_coag_parameters.R")
source("R/21_coag_odes.R")
source("R/22_coag_ckd_modifiers.R")
source("R/23_coag_simulate.R")
source("R/24_platelet_model.R")
source("R/25_platelet_simulate.R")
source("R/30_linkage.R")
source("R/40_tga_analysis.R")
source("R/41_pk_analysis.R")
source("R/42_bleeding_risk.R")

cat("=== Linkage Module Tests ===\n\n")
pass <- 0; fail <- 0

# Test 1: Unit conversion ng/mL -> nM
cat("Test 1: Unit conversion ng/mL to nM... ")
conc_nM <- ng_ml_to_nM(100)  # 100 ng/mL
expected <- 100 * 1000 / 459.5  # ~217.6 nM
if (abs(conc_nM - expected) < 0.1) {
  cat(sprintf("PASS (%.1f nM)\n", conc_nM)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (got %.1f, expected %.1f)\n", conc_nM, expected)); fail <- fail + 1
}

# Test 2: Round-trip conversion
cat("Test 2: Round-trip nM -> ng/mL -> nM... ")
original <- 200  # nM
round_trip <- ng_ml_to_nM(nM_to_ng_ml(original))
if (abs(round_trip - original) < 0.001) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 3: extract_trough_nM returns positive value
cat("Test 3: extract_trough_nM returns positive value... ")
sim <- simulate_pbpk(dose_mg = 5, n_doses = 14, egfr = 120)
trough_nM <- extract_trough_nM(sim)
if (trough_nM > 0) {
  cat(sprintf("PASS (%.2f nM)\n", trough_nM)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 4: NS modifiers applied correctly
cat("Test 4: NS modifiers increase FVIII... ")
ic_normal <- get_coag_initial_conditions()
ic_ns <- apply_ckd_modifiers(ic_normal, "NS_Severe")
if (ic_ns["VIII"] > ic_normal["VIII"]) {
  cat(sprintf("PASS (%.3f -> %.3f nM)\n", ic_normal["VIII"], ic_ns["VIII"])); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 5: NS modifiers decrease Protein C
cat("Test 5: NS modifiers decrease Protein C... ")
if (ic_ns["PC"] < ic_normal["PC"]) {
  cat(sprintf("PASS (%.1f -> %.1f nM)\n", ic_normal["PC"], ic_ns["PC"])); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 6: Full linked simulation runs
cat("Test 6: run_linked_simulation completes... ")
tryCatch({
  res <- run_linked_simulation(dose_mg = 5, egfr = 120, ckd_stage = "Normal")
  if (!is.null(res$tga_metrics) && !is.null(res$pk_metrics)) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (missing metrics)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 7: NS alters TGA metrics
cat("Test 7: NS severity alters TGA metrics vs Normal... ")
tryCatch({
  res_normal <- run_linked_simulation(dose_mg = 5, egfr = 100, ckd_stage = "Normal")
  res_ns <- run_linked_simulation(dose_mg = 5, egfr = 90, ckd_stage = "NS_Severe")
  etp_diff <- abs(res_ns$tga_metrics$ETP - res_normal$tga_metrics$ETP)
  if (etp_diff > 0) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (identical ETP)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 8: NS PK modifiers are applied in linked simulation, and applied with
# restrictive-clearance semantics.
#
# eGFR is held constant across stages: the original version of this test
# compared Normal at eGFR 100 with NS_Severe at eGFR 90, confounding renal
# function with disease severity.
#
# It also asserted free_ratio > 1.1, which encoded the artefact in which
# CL_hepatic and CL_renal were applied to TOTAL concentration and the fu
# multiplier passed straight through to unbound exposure. Because apixaban is
# restrictively cleared, fu alone cannot move unbound exposure at all: only the
# CLint_hepatic modifier can. Expectations are therefore derived from
# data/ns_pk_modifiers.csv so the test survives revision of those values.
cat("Test 8: NS modifiers propagate to linked PK with restrictive semantics... ")
tryCatch({
  clear_baseline_cache()
  res_normal <- run_linked_simulation(dose_mg = 5, egfr = 100, ckd_stage = "Normal")
  clear_baseline_cache()
  res_ns <- run_linked_simulation(dose_mg = 5, egfr = 100, ckd_stage = "NS_Severe")
  mods <- get_ckd_pk_modifiers("NS_Severe")

  free_ratio    <- res_ns$Cp_free_nM / res_normal$Cp_free_nM
  total_ratio   <- res_ns$pk_metrics$Ctrough_total / res_normal$pk_metrics$Ctrough_total
  freeauc_ratio <- res_ns$pk_metrics$AUCtau_free / res_normal$pk_metrics$AUCtau_free
  fu_mult       <- mods$fu_plasma
  clint_mult    <- mods$CLint_hepatic

  # AUC_unbound = F x Dose / CLint_total. Scaling only the hepatic component by
  # clint_mult bounds the rise in unbound AUC at 1/clint_mult, and renal
  # intrinsic clearance keeps the realised value strictly below that bound.
  ok <- free_ratio < fu_mult &&                 # no pass-through of the fu multiplier
        total_ratio < free_ratio &&             # total falls further than unbound
        freeauc_ratio <= 1 / clint_mult &&      # bounded by the CLint modifier
        if (clint_mult < 1) freeauc_ratio > 1 else freeauc_ratio <= 1.02

  if (ok) {
    cat(sprintf("PASS (trough free %.2fx, total %.2fx, AUC free %.2fx; fu mult %.2f, CLint mult %.2f)\n",
                free_ratio, total_ratio, freeauc_ratio, fu_mult, clint_mult)); pass <- pass + 1
  } else {
    cat(sprintf("FAIL (trough free %.2fx, total %.2fx, AUC free %.2fx; fu mult %.2f, CLint mult %.2f)\n",
                free_ratio, total_ratio, freeauc_ratio, fu_mult, clint_mult)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 9: run_linked_simulation returns hemostatic metrics and TI
cat("Test 9: run_linked_simulation returns HC and TI... ")
tryCatch({
  clear_baseline_cache()
  res <- run_linked_simulation(dose_mg = 5, egfr = 120, ckd_stage = "Normal",
                                compute_ti = TRUE)
  if (!is.null(res$hemostatic_metrics) &&
      res$hemostatic_metrics$HC > 0 &&
      !is.null(res$therapeutic_index) &&
      res$therapeutic_index$TI > 0) {
    cat(sprintf("PASS (HC=%.0f, TI=%.3f)\n",
                res$hemostatic_metrics$HC, res$therapeutic_index$TI))
    pass <- pass + 1
  } else {
    cat("FAIL\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
if (fail > 0) quit(status = 1)
