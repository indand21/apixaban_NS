# test_platelet.R
# Unit tests for platelet model and bleeding risk module

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

cat("=== Platelet Model Tests ===\n\n")
pass <- 0; fail <- 0

# Test 1: Platelet parameters load correctly
cat("Test 1: Platelet parameters load... ")
plt_params <- get_platelet_params()
if (length(plt_params) >= 14 && plt_params$k_thr > 0 && plt_params$PLT_rest_0 == 3000) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 2: NS platelet modifiers load
cat("Test 2: NS platelet modifier table loads... ")
mod_table <- get_ckd_platelet_modifier_table()
if (nrow(mod_table) >= 5 && "NS_Severe" %in% names(mod_table)) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 3: NS modifiers increase platelet parameters
cat("Test 3: NS_Severe increases k_thr... ")
p_normal <- get_platelet_params()
p_ns <- apply_ckd_platelet_modifiers(p_normal, "NS_Severe")
if (p_ns$k_thr > p_normal$k_thr) {
  cat(sprintf("PASS (%.4f -> %.4f)\n", p_normal$k_thr, p_ns$k_thr)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 4: Platelet initial conditions
cat("Test 4: Platelet ICs have 6 species... ")
ic <- get_platelet_initial_conditions()
if (length(ic) == 6 && ic["PLT_rest"] == 3000 && ic["PLT_act"] == 0) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 5: Platelet ODE runs with no-drug coag result
cat("Test 5: Platelet simulation runs (no drug)... ")
tryCatch({
  coag_nodrug <- simulate_coagulation(Cp_free_nM = 0)
  plt_res <- simulate_platelets(coag_nodrug, ckd_stage = "Normal")
  if (nrow(plt_res) > 10 && "PLT_agg" %in% names(plt_res) && max(plt_res$PLT_agg) > 0) {
    cat(sprintf("PASS (peak PLT_agg=%.1f nM)\n", max(plt_res$PLT_agg))); pass <- pass + 1
  } else {
    cat("FAIL (no aggregation)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 6: HC is positive for no-drug
cat("Test 6: HC is positive for no-drug baseline... ")
tryCatch({
  coag_nodrug <- simulate_coagulation(Cp_free_nM = 0)
  plt_nodrug <- simulate_platelets(coag_nodrug, ckd_stage = "Normal")
  hemo <- compute_hemostatic_metrics(plt_nodrug)
  if (hemo$HC > 0 && hemo$peak_PLT_agg > 0) {
    cat(sprintf("PASS (HC=%.0f nM*s, peak=%.1f nM)\n", hemo$HC, hemo$peak_PLT_agg))
    pass <- pass + 1
  } else {
    cat("FAIL\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 7: Drug reduces HC (more thrombin suppression = less platelet activation)
cat("Test 7: Drug reduces HC compared to no-drug... ")
tryCatch({
  coag_nodrug <- simulate_coagulation(Cp_free_nM = 0)
  coag_drug <- simulate_coagulation(Cp_free_nM = 100)  # High drug concentration
  plt_nodrug <- simulate_platelets(coag_nodrug, ckd_stage = "Normal")
  plt_drug <- simulate_platelets(coag_drug, ckd_stage = "Normal")
  hc_nodrug <- compute_hemostatic_metrics(plt_nodrug)$HC
  hc_drug <- compute_hemostatic_metrics(plt_drug)$HC
  if (hc_drug < hc_nodrug) {
    cat(sprintf("PASS (no-drug HC=%.0f, drug HC=%.0f)\n", hc_nodrug, hc_drug))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (no-drug=%.0f, drug=%.0f)\n", hc_nodrug, hc_drug)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 8: NS increases HC under same thrombin forcing
cat("Test 8: NS_Severe increases HC under same thrombin forcing... ")
tryCatch({
  coag_nodrug <- simulate_coagulation(Cp_free_nM = 0)
  plt_normal <- simulate_platelets(coag_nodrug, ckd_stage = "Normal")
  plt_ns <- simulate_platelets(coag_nodrug, ckd_stage = "NS_Severe")
  hc_normal <- compute_hemostatic_metrics(plt_normal)$HC
  hc_ns <- compute_hemostatic_metrics(plt_ns)$HC
  if (hc_ns > hc_normal) {
    cat(sprintf("PASS (Normal HC=%.0f, NS_Severe HC=%.0f)\n", hc_normal, hc_ns))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (Normal=%.0f, NS_Severe=%.0f)\n", hc_normal, hc_ns)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 9: Therapeutic index computation
cat("Test 9: TI = efficacy * safety, correct bounds... ")
tryCatch({
  # Primary endpoint is peak-thrombin suppression; the ETP variant is retained
  # alongside, so both are asserted here.
  tga_nodrug <- list(ETP = 10000, peak_thrombin = 600)
  tga_drug <- list(ETP = 5000, peak_thrombin = 300)
  hemo_nodrug <- list(HC = 100000)
  hemo_drug <- list(HC = 80000)
  ti <- compute_therapeutic_index(tga_drug, tga_nodrug, hemo_drug, hemo_nodrug)
  expected_eff <- 1 - 300/600      # 0.5 from peak thrombin
  expected_saf <- 80000/100000     # 0.8
  expected_ti <- 0.5 * 0.8         # 0.4
  if (abs(ti$efficacy_score - expected_eff) < 1e-6 &&
      abs(ti$safety_score - expected_saf) < 1e-6 &&
      abs(ti$TI - expected_ti) < 1e-6 &&
      identical(ti$objective, "peak") &&
      abs(ti$efficacy_score_etp - (1 - 5000/10000)) < 1e-6 &&
      abs(ti$TI_etp - 0.4) < 1e-6) {
    cat(sprintf("PASS (TI=%.3f)\n", ti$TI)); pass <- pass + 1
  } else {
    cat(sprintf("FAIL (got TI=%.3f, expected %.3f)\n", ti$TI, expected_ti)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 10: TI is 0 when dose is 0 (no efficacy)
cat("Test 10: TI boundaries — no drug means TI=0... ")
ti_zero <- compute_therapeutic_index(
  tga_drug = list(ETP = 10000, peak_thrombin = 600),
  tga_nodrug = list(ETP = 10000, peak_thrombin = 600),
  hemo_drug = list(HC = 100000), hemo_nodrug = list(HC = 100000)
)
if (abs(ti_zero$TI) < 1e-6) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat(sprintf("FAIL (TI=%.6f, expected 0)\n", ti_zero$TI)); fail <- fail + 1
}

# Test 11: Full linked simulation with TI
cat("Test 11: run_linked_simulation returns TI... ")
tryCatch({
  clear_baseline_cache()
  res <- run_linked_simulation(dose_mg = 5, egfr = 120, ckd_stage = "Normal",
                                compute_ti = TRUE)
  if (!is.null(res$therapeutic_index) &&
      res$therapeutic_index$TI > 0 &&
      res$therapeutic_index$TI < 1 &&
      !is.null(res$hemostatic_metrics) &&
      res$hemostatic_metrics$HC > 0) {
    cat(sprintf("PASS (TI=%.3f)\n", res$therapeutic_index$TI)); pass <- pass + 1
  } else {
    cat("FAIL (bad TI or HC values)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 12: NS_Severe produces a different TI than Normal at same dose
cat("Test 12: NS_Severe TI differs from Normal at 5 mg... ")
tryCatch({
  clear_baseline_cache()
  res_normal <- run_linked_simulation(dose_mg = 5, egfr = 120, ckd_stage = "Normal",
                                       compute_ti = TRUE)
  res_ns <- run_linked_simulation(dose_mg = 5, egfr = 90, ckd_stage = "NS_Severe",
                                  compute_ti = TRUE)
  if (abs(res_ns$therapeutic_index$TI - res_normal$therapeutic_index$TI) > 1e-6) {
    cat(sprintf("PASS (Normal=%.3f, NS_Severe=%.3f)\n",
                res_normal$therapeutic_index$TI, res_ns$therapeutic_index$TI))
    pass <- pass + 1
  } else {
    cat("FAIL\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
