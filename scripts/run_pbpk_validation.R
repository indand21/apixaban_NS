stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_pbpk_validation.R
# Validate PBPK model against published apixaban PK data

rm(list = ls())

source("R/00_packages.R")
source("R/01_utils.R")
source("R/10_pbpk_parameters.R")
source("R/11_pbpk_odes.R")
source("R/12_pbpk_simulate.R")
source("R/41_pk_analysis.R")
source("R/51_tables.R")

cat("=== PBPK Model Validation ===\n\n")

# --- Published reference values (5 mg BID, steady-state, normal renal function) ---
# Sources: FDA Clinical Pharmacology Review, Byon et al. 2019
observed <- list(
  Cmax_total    = 171,   # ng/mL (range ~130-230)
  Ctrough_total = 51,    # ng/mL (range ~30-80)
  AUCtau_total  = 1100   # ng*h/mL (approximate)
)

# --- Simulate 5 mg BID for normal renal function ---
cat("Simulating 5 mg BID, eGFR=120 (Normal)...\n")
sim <- simulate_pbpk(dose_mg = 5, tau = 12, n_doses = 14, egfr = 120)
pk <- compute_pk_metrics(sim, tau = 12)

cat("\n--- Predicted vs Observed (Normal, 5 mg BID) ---\n")
cat(sprintf("  Cmax (total):    Predicted = %.1f ng/mL | Observed = %d ng/mL | Ratio = %.2f\n",
            pk$Cmax_total, observed$Cmax_total, pk$Cmax_total / observed$Cmax_total))
cat(sprintf("  Ctrough (total): Predicted = %.1f ng/mL | Observed = %d ng/mL | Ratio = %.2f\n",
            pk$Ctrough_total, observed$Ctrough_total, pk$Ctrough_total / observed$Ctrough_total))
cat(sprintf("  AUCtau (total):  Predicted = %.0f ng*h/mL | Observed = %d ng*h/mL | Ratio = %.2f\n",
            pk$AUCtau_total, observed$AUCtau_total, pk$AUCtau_total / observed$AUCtau_total))

# --- Acceptance criteria ---
cat("\n--- Acceptance Criteria (within 0.5-2.0 fold) ---\n")
ratios <- c(
  Cmax    = pk$Cmax_total / observed$Cmax_total,
  Ctrough = pk$Ctrough_total / observed$Ctrough_total,
  AUCtau  = pk$AUCtau_total / observed$AUCtau_total
)

for (nm in names(ratios)) {
  pass <- ratios[nm] >= 0.5 & ratios[nm] <= 2.0
  cat(sprintf("  %s: ratio=%.2f %s\n", nm, ratios[nm], ifelse(pass, "PASS", "FAIL")))
}

# --- CKD stage validation ---
cat("\n\n--- PK across CKD stages (5 mg BID) ---\n")
# Published: CKD3-4 show ~30-40% increase in AUC (FDA label)
pk_ckd <- compute_pk_across_ckd(dose_mg = 5)
print(pk_ckd[, c("ckd_stage", "egfr", "Cmax_total", "Ctrough_total", "AUCtau_total")], digits = 3)

# Check directional: AUC should increase with decreasing eGFR
cat("\nDirectional check: AUCtau increases with decreasing eGFR\n")
auc_increasing <- all(diff(pk_ckd$AUCtau_total) >= 0)
cat(sprintf("  Monotonic increase: %s\n", ifelse(auc_increasing, "PASS", "FAIL")))

# --- Plot validation ---
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Steady-state profile
ss <- extract_steady_state(sim)
p <- ggplot(ss, aes(x = time_in_interval, y = C_total_ng_mL)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_hline(yintercept = 171, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 51, linetype = "dashed", color = "orange") +
  annotate("text", x = 1, y = 180, label = "Observed Cmax (171)", color = "red", size = 3) +
  annotate("text", x = 8, y = 60, label = "Observed Ctrough (51)", color = "orange", size = 3) +
  labs(x = "Time in dosing interval (h)",
       y = "Total plasma concentration (ng/mL)",
       title = "PBPK Validation: Apixaban 5 mg BID (Normal renal function)") +
  theme_bw()

ggsave("output/figures/pbpk_validation.png", p, width = 8, height = 5, dpi = 300)
cat("\nValidation plot saved to output/figures/pbpk_validation.png\n")
