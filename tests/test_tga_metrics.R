# test_tga_metrics.R
# Unit tests for TGA metric extraction

source("R/00_packages.R")
source("R/01_utils.R")
source("R/20_coag_parameters.R")
source("R/21_coag_odes.R")
source("R/23_coag_simulate.R")
source("R/40_tga_analysis.R")

cat("=== TGA Metrics Tests ===\n\n")
pass <- 0; fail <- 0

# Run a baseline simulation for testing
sim <- simulate_coagulation(Cp_free_nM = 0, TF_pM = 25, t_end = 1200)

# Test 1: ETP is positive
cat("Test 1: ETP is positive... ")
metrics <- compute_tga_metrics(sim)
if (metrics$ETP > 0) {
  cat(sprintf("PASS (ETP=%.0f nM*s)\n", metrics$ETP)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 2: Peak thrombin matches max of IIa column
cat("Test 2: Peak thrombin matches max(IIa)... ")
expected_peak <- max(sim$IIa)
if (abs(metrics$peak_thrombin - expected_peak) < 0.001) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 3: Lag time < time-to-peak
cat("Test 3: Lag time < time-to-peak... ")
if (!is.na(metrics$lag_time) && metrics$lag_time < metrics$time_to_peak) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 4: ETP_nM_min = ETP / 60
cat("Test 4: ETP_nM_min = ETP / 60... ")
if (abs(metrics$ETP_nM_min - metrics$ETP / 60) < 0.01) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 5: Velocity index is positive
cat("Test 5: Velocity index is positive... ")
if (!is.na(metrics$velocity_index) && metrics$velocity_index > 0) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 6: Drug monotonically reduces ETP
cat("Test 6: ETP monotonically decreases with drug concentration... ")
concs <- c(0, 10, 50, 200)
etps <- sapply(concs, function(c) {
  s <- simulate_coagulation(Cp_free_nM = c, TF_pM = 25, t_end = 1200)
  compute_tga_metrics(s)$ETP
})
monotonic <- all(diff(etps) <= 0)
if (monotonic) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat(sprintf("FAIL (ETPs: %s)\n", paste(round(etps), collapse = ", "))); fail <- fail + 1
}

# Test 7: compute_tga_metrics_df returns one-row data frame
cat("Test 7: compute_tga_metrics_df returns data frame... ")
df <- compute_tga_metrics_df(sim, ckd_stage = "Normal", dose_mg = 5)
if (is.data.frame(df) && nrow(df) == 1 && "ckd_stage" %in% names(df)) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
