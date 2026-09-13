stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_coag_validation.R
# Validate Hockin-Mann coagulation model against published TGA data

rm(list = ls())

source("R/00_packages.R")
source("R/01_utils.R")
source("R/20_coag_parameters.R")
source("R/21_coag_odes.R")
source("R/23_coag_simulate.R")
source("R/40_tga_analysis.R")

cat("=== Hockin-Mann Coagulation Model Validation ===\n\n")

# --- Baseline: 25 pM TF, no drug ---
cat("Simulating baseline thrombin generation (25 pM TF, no drug)...\n")
baseline <- simulate_coagulation(Cp_free_nM = 0, TF_pM = 25, t_end = 1200)
baseline_metrics <- compute_tga_metrics(baseline)

cat("\n--- Baseline TGA Metrics ---\n")
cat(sprintf("  Peak thrombin: %.1f nM (expected: 300-400 nM)\n", baseline_metrics$peak_thrombin))
cat(sprintf("  Time-to-peak:  %.0f s (expected: 300-500 s)\n", baseline_metrics$time_to_peak))
cat(sprintf("  Lag time:      %.0f s (expected: 60-100 s)\n", baseline_metrics$lag_time))
cat(sprintf("  ETP:           %.0f nM*s (%.0f nM*min)\n", baseline_metrics$ETP, baseline_metrics$ETP_nM_min))

# --- Acceptance criteria ---
cat("\n--- Acceptance Criteria ---\n")
checks <- list(
  peak_thrombin = list(val = baseline_metrics$peak_thrombin, lo = 100, hi = 600, name = "Peak thrombin (nM)"),
  lag_time      = list(val = baseline_metrics$lag_time, lo = 30, hi = 200, name = "Lag time (s)"),
  time_to_peak  = list(val = baseline_metrics$time_to_peak, lo = 200, hi = 600, name = "Time-to-peak (s)")
)

for (ch in checks) {
  pass <- !is.na(ch$val) & ch$val >= ch$lo & ch$val <= ch$hi
  cat(sprintf("  %s: %.1f [%.0f-%.0f] %s\n",
              ch$name, ch$val, ch$lo, ch$hi, ifelse(pass, "PASS", "CHECK")))
}

# --- TF dose-response ---
cat("\n--- TF Dose-Response ---\n")
tf_levels <- c(1, 5, 10, 25, 50, 100)  # pM
tf_results <- purrr::map_dfr(tf_levels, function(tf) {
  sim <- simulate_coagulation(Cp_free_nM = 0, TF_pM = tf, t_end = 1200)
  metrics <- compute_tga_metrics(sim)
  data.frame(
    TF_pM = tf,
    peak_thrombin = metrics$peak_thrombin,
    lag_time = metrics$lag_time,
    ETP = metrics$ETP
  )
})
print(tf_results, digits = 3)
cat("Expected: Peak thrombin increases and lag time decreases with increasing TF\n")

# --- Drug dose-response (apixaban) ---
cat("\n--- Apixaban Dose-Response ---\n")
drug_concs <- c(0, 1, 5, 10, 50, 100, 500)  # nM
drug_results <- purrr::map_dfr(drug_concs, function(conc) {
  sim <- simulate_coagulation(Cp_free_nM = conc, TF_pM = 25, t_end = 1200)
  metrics <- compute_tga_metrics(sim)
  data.frame(
    Cp_free_nM    = conc,
    peak_thrombin = metrics$peak_thrombin,
    lag_time      = metrics$lag_time,
    ETP           = metrics$ETP,
    ETP_pct       = metrics$ETP / baseline_metrics$ETP * 100
  )
})
print(drug_results, digits = 3)
cat("Expected: Monotonic decrease in ETP with increasing drug concentration\n")

etp_monotonic <- all(diff(drug_results$ETP) <= 0)
cat(sprintf("  Monotonic ETP decrease: %s\n", ifelse(etp_monotonic, "PASS", "FAIL")))

# --- Plots ---
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Baseline thrombin generation curve
p1 <- ggplot(baseline, aes(x = time / 60, y = IIa)) +
  geom_line(linewidth = 1, color = "darkred") +
  labs(x = "Time (min)", y = "Thrombin (nM)",
       title = "Baseline thrombin generation (25 pM TF, no drug)") +
  theme_bw()

# Drug dose-response curves
drug_sims <- purrr::map_dfr(drug_concs, function(conc) {
  sim <- simulate_coagulation(Cp_free_nM = conc, TF_pM = 25, t_end = 1200)
  data.frame(
    time_min   = sim$time / 60,
    IIa        = sim$IIa,
    conc_label = paste0(conc, " nM")
  )
})
drug_sims$conc_label <- factor(drug_sims$conc_label,
                               levels = paste0(drug_concs, " nM"))

p2 <- ggplot(drug_sims, aes(x = time_min, y = IIa, color = conc_label)) +
  geom_line(linewidth = 0.7) +
  labs(x = "Time (min)", y = "Thrombin (nM)",
       color = "Apixaban (free)",
       title = "Thrombin generation with apixaban dose-response") +
  theme_bw() +
  theme(legend.position = "right")

library(patchwork)
p_combined <- p1 / p2
ggsave("output/figures/coag_validation.png", p_combined, width = 10, height = 10, dpi = 300)
cat("\nValidation plots saved to output/figures/coag_validation.png\n")
