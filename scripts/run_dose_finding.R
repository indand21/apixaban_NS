stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_dose_finding.R
# Find optimal dose per nephrotic syndrome severity maximizing therapeutic index (TI)
# TI = Efficacy_score * Safety_score balances thrombin suppression vs bleeding risk

rm(list = ls())

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

# --- Configuration ---
ckd_stages  <- c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")
egfr_values <- c(Normal = 100, NS_Mild = 90, NS_Moderate = 90, NS_Severe = 90)

# Clear baseline cache to start fresh
clear_baseline_cache()

# --- TI objective: negative TI for minimize (optimize maximizes by default) ---
ti_objective <- function(dose_mg, egfr, ckd_stage) {
  res <- run_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
    compute_ti = TRUE
  )
  res$therapeutic_index$TI
}

# --- Find optimal dose for each disease stage ---
cat("Finding TI-optimal dose per nephrotic syndrome severity...\n")
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
    # Maximize TI over dose range 0.5-20 mg
    opt <- optimize(
      ti_objective,
      interval = c(0.5, 20),
      egfr = egfr,
      ckd_stage = stage,
      maximum = TRUE,
      tol = 0.05
    )

    opt_dose <- opt$maximum
    opt_TI   <- opt$objective

    # Run at optimal dose to get full metrics
    verify <- run_linked_simulation(
      dose_mg = opt_dose, egfr = egfr, ckd_stage = stage,
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
      ETP = NA, HC = NA,
      stringsAsFactors = FALSE
    ))
  })
}

# --- Also compute PD-equivalent doses (legacy) for comparison ---
cat("\n--- PD-Equivalent Doses (ETP-matching, for comparison) ---\n")
ref_result <- run_linked_simulation(dose_mg = 5, egfr = 120, ckd_stage = "Normal")
ref_ETP <- ref_result$tga_metrics$ETP

etp_objective <- function(dose_mg, egfr, ckd_stage, target_etp) {
  res <- run_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
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

# --- Output ---
cat("\n=== DOSE OPTIMIZATION RESULTS ===\n")
print(dose_results, digits = 3)

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(dose_results, "output/tables/optimal_doses.csv", row.names = FALSE)

# --- Plots ---
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# TI vs dose curve for each nephrotic syndrome severity
cat("\nGenerating TI-vs-dose curves...\n")
dose_sweep <- seq(0.5, 15, by = 0.5)
ti_curves <- list()
for (stage in ckd_stages) {
  egfr <- egfr_values[stage]
  for (d in dose_sweep) {
    tryCatch({
      res <- run_linked_simulation(
        dose_mg = d, egfr = egfr, ckd_stage = stage, compute_ti = TRUE
      )
      ti_curves[[length(ti_curves) + 1]] <- data.frame(
        dose_mg        = d,
        ckd_stage      = stage,
        TI             = res$therapeutic_index$TI,
        efficacy_score = res$therapeutic_index$efficacy_score,
        safety_score   = res$therapeutic_index$safety_score,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {})
  }
}
ti_curve_df <- do.call(rbind, ti_curves)
ti_curve_df$ckd_stage <- factor(ti_curve_df$ckd_stage, levels = ckd_stages)

p_ti <- ggplot(ti_curve_df, aes(x = dose_mg, y = TI, color = ckd_stage)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = dose_results, aes(x = optimal_dose, y = TI),
             size = 3, shape = 18) +
  labs(x = "Dose (mg BID)", y = "Therapeutic Index",
       color = "NS Severity",
       title = "Therapeutic index vs dose by nephrotic syndrome severity",
       subtitle = "Diamonds = optimal dose; TI = Efficacy x Safety") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/figures/ti_vs_dose.png", p_ti, width = 10, height = 6, dpi = 300)

# Comparison: TI-optimal vs PD-equivalent doses
p_compare <- plot_dose_comparison(dose_results)
ggsave("output/figures/dose_comparison.png", p_compare, width = 10, height = 6, dpi = 300)

cat("\nUncoupled dose-finding complete.\n")

# =================================================================
# COUPLED (Two-Way) Dose-Finding
# =================================================================
cat("\n=================================================================\n")
cat("COUPLED dose-finding: PS_exp -> prothrombinase feedback\n")
cat("=================================================================\n\n")

clear_coupled_baseline_cache()

# Coupled TI objective
ti_objective_coupled <- function(dose_mg, egfr, ckd_stage) {
  res <- run_coupled_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
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
      interval = c(0.5, 20),
      egfr = egfr,
      ckd_stage = stage,
      maximum = TRUE,
      tol = 0.05
    )

    opt_dose <- opt$maximum
    verify <- run_coupled_linked_simulation(
      dose_mg = opt_dose, egfr = egfr, ckd_stage = stage,
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
      ETP = NA, HC = NA,
      stringsAsFactors = FALSE
    ))
  })
}

cat("\n=== COUPLED DOSE OPTIMIZATION RESULTS ===\n")
print(coupled_dose_results, digits = 3)

write.csv(coupled_dose_results, "output/tables/coupled_optimal_doses.csv", row.names = FALSE)

# --- Comparison table ---
cat("\n=== UNCOUPLED vs COUPLED OPTIMAL DOSES ===\n")
compare <- merge(
  dose_results[, c("ckd_stage", "optimal_dose", "TI")],
  coupled_dose_results[, c("ckd_stage", "optimal_dose", "TI")],
  by = "ckd_stage",
  suffixes = c("_uncoupled", "_coupled")
)
compare$dose_shift <- compare$optimal_dose_coupled - compare$optimal_dose_uncoupled
print(compare, digits = 3)

# --- Coupled TI-vs-dose curves ---
cat("\nGenerating coupled TI-vs-dose curves...\n")
coupled_dose_sweep <- seq(0.5, 15, by = 0.5)
coupled_ti_curves <- list()
for (stage in ckd_stages) {
  egfr <- egfr_values[stage]
  for (d in coupled_dose_sweep) {
    tryCatch({
      res <- run_coupled_linked_simulation(
        dose_mg = d, egfr = egfr, ckd_stage = stage, compute_ti = TRUE
      )
      coupled_ti_curves[[length(coupled_ti_curves) + 1]] <- data.frame(
        dose_mg        = d,
        ckd_stage      = stage,
        TI             = res$therapeutic_index$TI,
        efficacy_score = res$therapeutic_index$efficacy_score,
        safety_score   = res$therapeutic_index$safety_score,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {})
  }
}
coupled_ti_curve_df <- do.call(rbind, coupled_ti_curves)
coupled_ti_curve_df$ckd_stage <- factor(coupled_ti_curve_df$ckd_stage, levels = ckd_stages)

p_ti_coupled <- ggplot(coupled_ti_curve_df, aes(x = dose_mg, y = TI, color = ckd_stage)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = coupled_dose_results, aes(x = optimal_dose, y = TI),
             size = 3, shape = 18) +
  labs(x = "Dose (mg BID)", y = "Therapeutic Index",
       color = "NS Severity",
       title = "COUPLED: Therapeutic index vs dose by nephrotic syndrome severity",
       subtitle = "Diamonds = optimal dose; TI = Efficacy x Safety; PS_exp feedback active") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("output/figures/coupled_ti_vs_dose.png", p_ti_coupled, width = 10, height = 6, dpi = 300)

cat("\nAll results saved to output/\n")
