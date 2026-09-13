stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_dose_finding_TF5.R
# Nephrotic syndrome dose finding: TF = 5 pM (CAT assay protocol).
# Identical to run_dose_finding.R except every run_linked_simulation() and
# run_coupled_linked_simulation() call explicitly passes TF_pM = 5.
# Outputs overwrite the canonical optimal_doses.csv and coupled_optimal_doses.csv;
# the previous TF=25 versions are preserved as *_TF25_default.csv.

rm(list = ls())

TF_pM_RUN <- 5

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
source("R/26_coupled_odes.R")
source("R/27_coupled_simulate.R")
source("R/30_linkage.R")
source("R/40_tga_analysis.R")
source("R/41_pk_analysis.R")
source("R/42_bleeding_risk.R")
source("R/50_visualization.R")

ckd_stages  <- c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")
egfr_values <- c(Normal = 100, NS_Mild = 90, NS_Moderate = 90, NS_Severe = 90)

clear_baseline_cache()

ti_objective <- function(dose_mg, egfr, ckd_stage) {
  res <- run_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
    TF_pM = TF_pM_RUN,
    compute_ti = TRUE
  )
  res$therapeutic_index$TI
}

cat("Finding TI-optimal dose per nephrotic syndrome severity (TF =", TF_pM_RUN, "pM)...\n")
cat("==========================================\n\n")

dose_results <- data.frame(
  ckd_stage      = character(),
  egfr           = numeric(),
  optimal_dose   = numeric(),
  TI             = numeric(),
  efficacy_score = numeric(),
  safety_score   = numeric(),
  ETP            = numeric(),
  HC             = numeric(),
  stringsAsFactors = FALSE
)

for (stage in ckd_stages) {
  egfr <- egfr_values[stage]
  cat(sprintf("Optimizing dose for %s (eGFR=%d)...", stage, egfr))

  tryCatch({
    opt <- optimize(
      ti_objective,
      interval = c(0.1, 20),
      egfr = egfr, ckd_stage = stage,
      maximum = TRUE, tol = 0.05
    )
    opt_dose <- opt$maximum

    verify <- run_linked_simulation(
      dose_mg = opt_dose, egfr = egfr, ckd_stage = stage,
      TF_pM = TF_pM_RUN,
      compute_ti = TRUE
    )

    dose_results <- rbind(dose_results, data.frame(
      ckd_stage      = stage,
      egfr           = egfr,
      optimal_dose   = round(opt_dose, 2),
      TI             = verify$therapeutic_index$TI,
      efficacy_score = verify$therapeutic_index$efficacy_score,
      safety_score   = verify$therapeutic_index$safety_score,
      ETP            = verify$tga_metrics$ETP,
      HC             = verify$hemostatic_metrics$HC,
      stringsAsFactors = FALSE
    ))

    cat(sprintf(" %.2f mg BID (TI=%.3f, Eff=%.2f, Saf=%.2f)\n",
                opt_dose, verify$therapeutic_index$TI,
                verify$therapeutic_index$efficacy_score,
                verify$therapeutic_index$safety_score))
  }, error = function(e) {
    cat(sprintf(" FAILED: %s\n", e$message))
    dose_results <<- rbind(dose_results, data.frame(
      ckd_stage = stage, egfr = egfr,
      optimal_dose = NA, TI = NA,
      efficacy_score = NA, safety_score = NA,
      ETP = NA, HC = NA, stringsAsFactors = FALSE
    ))
  })
}

# --- PD-equivalent doses (ETP-matching) ---
cat("\n--- PD-Equivalent Doses (ETP-matching) ---\n")
ref_result <- run_linked_simulation(
  dose_mg = 5, egfr = 120, ckd_stage = "Normal",
  TF_pM = TF_pM_RUN
)
ref_ETP <- ref_result$tga_metrics$ETP

etp_objective <- function(dose_mg, egfr, ckd_stage, target_etp) {
  res <- run_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
    TF_pM = TF_pM_RUN,
    compute_ti = FALSE
  )
  res$tga_metrics$ETP - target_etp
}

dose_results$pd_equiv_dose <- NA
dose_results$pd_equiv_dose[dose_results$ckd_stage == "Normal"] <- 5.0

for (stage in ckd_stages[ckd_stages != "Normal"]) {
  egfr <- egfr_values[stage]
  tryCatch({
    root <- uniroot(
      etp_objective,
      interval = c(0.5, 20),
      egfr = egfr, ckd_stage = stage, target_etp = ref_ETP,
      tol = 0.01
    )
    dose_results$pd_equiv_dose[dose_results$ckd_stage == stage] <- round(root$root, 2)
    cat(sprintf("  %s: PD-equiv = %.2f mg\n", stage, root$root))
  }, error = function(e) {
    cat(sprintf("  %s: PD-equiv FAILED (%s)\n", stage, e$message))
  })
}

cat("\n=== UNCOUPLED DOSE OPTIMIZATION RESULTS (TF =", TF_pM_RUN, "pM) ===\n")
print(dose_results, digits = 3)

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(dose_results, "output/tables/optimal_doses.csv", row.names = FALSE)

# =================================================================
# COUPLED Dose-Finding
# =================================================================
cat("\n=================================================================\n")
cat("COUPLED dose-finding (TF =", TF_pM_RUN, "pM)\n")
cat("=================================================================\n\n")

clear_coupled_baseline_cache()

ti_objective_coupled <- function(dose_mg, egfr, ckd_stage) {
  res <- run_coupled_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
    TF_pM = TF_pM_RUN,
    compute_ti = TRUE
  )
  res$therapeutic_index$TI
}

coupled_dose_results <- data.frame(
  ckd_stage      = character(),
  egfr           = numeric(),
  optimal_dose   = numeric(),
  TI             = numeric(),
  efficacy_score = numeric(),
  safety_score   = numeric(),
  ETP            = numeric(),
  HC             = numeric(),
  stringsAsFactors = FALSE
)

for (stage in ckd_stages) {
  egfr <- egfr_values[stage]
  cat(sprintf("COUPLED optimizing dose for %s (eGFR=%d)...", stage, egfr))

  tryCatch({
    opt <- optimize(
      ti_objective_coupled,
      interval = c(0.1, 20),
      egfr = egfr, ckd_stage = stage,
      maximum = TRUE, tol = 0.05
    )
    opt_dose <- opt$maximum

    verify <- run_coupled_linked_simulation(
      dose_mg = opt_dose, egfr = egfr, ckd_stage = stage,
      TF_pM = TF_pM_RUN,
      compute_ti = TRUE
    )

    coupled_dose_results <- rbind(coupled_dose_results, data.frame(
      ckd_stage      = stage,
      egfr           = egfr,
      optimal_dose   = round(opt_dose, 2),
      TI             = verify$therapeutic_index$TI,
      efficacy_score = verify$therapeutic_index$efficacy_score,
      safety_score   = verify$therapeutic_index$safety_score,
      ETP            = verify$tga_metrics$ETP,
      HC             = verify$hemostatic_metrics$HC,
      stringsAsFactors = FALSE
    ))

    cat(sprintf(" %.2f mg BID (TI=%.3f, Eff=%.2f, Saf=%.2f)\n",
                opt_dose, verify$therapeutic_index$TI,
                verify$therapeutic_index$efficacy_score,
                verify$therapeutic_index$safety_score))
  }, error = function(e) {
    cat(sprintf(" FAILED: %s\n", e$message))
    coupled_dose_results <<- rbind(coupled_dose_results, data.frame(
      ckd_stage = stage, egfr = egfr,
      optimal_dose = NA, TI = NA,
      efficacy_score = NA, safety_score = NA,
      ETP = NA, HC = NA, stringsAsFactors = FALSE
    ))
  })
}

cat("\n=== COUPLED DOSE OPTIMIZATION RESULTS (TF =", TF_pM_RUN, "pM) ===\n")
print(coupled_dose_results, digits = 3)

write.csv(coupled_dose_results, "output/tables/coupled_optimal_doses.csv", row.names = FALSE)

# --- Comparison ---
cat("\n=== UNCOUPLED vs COUPLED OPTIMAL DOSES (TF =", TF_pM_RUN, "pM) ===\n")
compare <- merge(
  dose_results[, c("ckd_stage", "optimal_dose", "TI")],
  coupled_dose_results[, c("ckd_stage", "optimal_dose", "TI")],
  by = "ckd_stage",
  suffixes = c("_uncoupled", "_coupled")
)
compare$dose_shift <- compare$optimal_dose_coupled - compare$optimal_dose_uncoupled
print(compare, digits = 3)

cat("\nAll results saved to output/tables/\n")
