# test_coagulation.R
# Unit tests for Hockin-Mann coagulation module

source("R/00_packages.R")
source("R/01_utils.R")
source("R/20_coag_parameters.R")
source("R/21_coag_odes.R")
source("R/22_coag_ckd_modifiers.R")
source("R/23_coag_simulate.R")
source("R/40_tga_analysis.R")

cat("=== Coagulation Module Tests ===\n\n")
pass <- 0; fail <- 0

# Test 1: Rate constants loaded
cat("Test 1: Rate constants are complete (k1-k41)... ")
rc <- get_coag_rate_constants()
expected_keys <- paste0("k", setdiff(1:41,17))
if (all(expected_keys %in% names(rc))) {
  cat("PASS\n"); pass <- pass + 1
} else {
  missing <- setdiff(expected_keys, names(rc))
  cat(sprintf("FAIL (missing: %s)\n", paste(missing, collapse = ", "))); fail <- fail + 1
}

# Test 2: Initial conditions have 34 species
cat("Test 2: Initial conditions have 34 species... ")
ic <- get_coag_initial_conditions()
if (length(ic) == 34) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat(sprintf("FAIL (got %d species)\n", length(ic))); fail <- fail + 1
}

# Test 3: All rate constants are positive
cat("Test 3: Active rates positive and implicit TM pathway disabled... ")
all_pos <- all(unlist(rc[names(rc)!='k39']) > 0) && rc$k39 == 0
if (all_pos) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 4: Baseline simulation runs
cat("Test 4: Baseline simulation runs without error... ")
tryCatch({
  sim <- simulate_coagulation(Cp_free_nM = 0, TF_pM = 25, t_end = 600)
  if (nrow(sim) > 0) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (empty result)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 5: Thrombin is generated (peak > 10 nM)
cat("Test 5: Thrombin is generated (peak > 10 nM)... ")
sim <- simulate_coagulation(Cp_free_nM = 0, TF_pM = 25, t_end = 1200)
peak <- max(sim$IIa)
if (peak > 10) {
  cat(sprintf("PASS (peak=%.1f nM)\n", peak)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (peak=%.1f nM)\n", peak)); fail <- fail + 1
}

# Test 6: No negative species
cat("Test 6: No negative species concentrations... ")
# Check all species columns (exclude time and auxiliary)
species_cols <- names(sim)[2:35]  # First 34 after time
has_neg <- any(sim[, species_cols] < -1e-6, na.rm = TRUE)
if (!has_neg) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 7: Mass conservation - total factor II (II + IIa + mIIa + bound forms)
cat("Test 7: Prothrombin mass conservation (approximate)... ")
II_total_start <- ic["II"]
II_total_end <- sim$II[nrow(sim)] + sim$IIa[nrow(sim)] + sim$mIIa[nrow(sim)] +
                sim$Va_Xa_II[nrow(sim)] + sim$ATIII_IIa[nrow(sim)]
ratio <- II_total_end / II_total_start
if (abs(ratio - 1) < 0.05) {
  cat(sprintf("PASS (ratio=%.3f)\n", ratio)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (ratio=%.3f)\n", ratio)); fail <- fail + 1
}

# Test 8: Drug reduces thrombin generation
cat("Test 8: Apixaban reduces peak thrombin... ")
sim_drug <- simulate_coagulation(Cp_free_nM = 100, TF_pM = 25, t_end = 1200)
peak_drug <- max(sim_drug$IIa)
if (peak_drug < peak) {
  cat(sprintf("PASS (%.1f < %.1f nM)\n", peak_drug, peak)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 9: exclude_species = "TF" leaves TF unmodified
cat("Test 9: exclude_species='TF' leaves TF unmodified... ")
ic_base <- get_coag_initial_conditions(TF_pM = 5)
tf_before <- ic_base["TF"]
ic_excl <- apply_ckd_modifiers(ic_base, ckd_stage = "NS_Severe", exclude_species = "TF")
tf_after <- ic_excl["TF"]
viii_after <- ic_excl["VIII"]
# TF should be unchanged, but VIII should be modified in severe NS
if (tf_after == tf_before && viii_after != ic_base["VIII"]) {
  cat(sprintf("PASS (TF=%.4f unchanged, VIII modified)\n", tf_after)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (TF: %.4f->%.4f, VIII: %.4f->%.4f)\n",
              tf_before, tf_after, ic_base["VIII"], viii_after)); fail <- fail + 1
}

# Test 10: exclude_species = NULL gives same result as original
cat("Test 10: Default fixed TF differs from explicit endogenous TF scaling... ")
ic_default <- apply_ckd_modifiers(ic_base, ckd_stage = "NS_Severe")
ic_null    <- apply_ckd_modifiers(ic_base, ckd_stage = "NS_Severe", exclude_species = NULL)
if (ic_default['TF']==ic_base['TF'] && ic_null['TF']>ic_default['TF'] && identical(ic_default[names(ic_default)!='TF'],ic_null[names(ic_null)!='TF'])) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL (results differ)\n"); fail <- fail + 1
}

# Test 11: Age coag modifiers — age=80 increases VIII, decreases PC
cat("Test 11: Age=80 increases VIII and decreases PC vs baseline... ")
ic_base80 <- get_coag_initial_conditions(TF_pM = 25)
ic_aged   <- apply_age_coag_modifiers(ic_base80, age = 80)
viii_ratio <- ic_aged["VIII"] / ic_base80["VIII"]
pc_ratio   <- ic_aged["PC"] / ic_base80["PC"]
# VIII: 1 + 0.05*(80-40)/10 = 1.20; PC: 1 - 0.03*(80-50)/10 = 0.91
if (viii_ratio > 1.0 && pc_ratio < 1.0 && abs(viii_ratio - 1.20) < 0.01 &&
    abs(pc_ratio - 0.91) < 0.01) {
  cat(sprintf("PASS (VIII=%.2fx, PC=%.2fx)\n", viii_ratio, pc_ratio)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (VIII=%.2fx, PC=%.2fx)\n", viii_ratio, pc_ratio)); fail <- fail + 1
}

# Test 12: age=NULL returns ic unchanged
cat("Test 12: apply_age_coag_modifiers with age=NULL returns ic unchanged... ")
ic_null <- apply_age_coag_modifiers(ic_base80, age = NULL)
if (identical(ic_base80, ic_null)) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
if (fail > 0) quit(status=1)
