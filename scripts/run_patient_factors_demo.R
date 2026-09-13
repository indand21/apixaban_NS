stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_patient_factors_demo.R
# Explore body weight x age x nephrotic syndrome severity effects on PK and optimal dosing
# Produces heatmaps and tables for FDA 2-of-3 dose-reduction criteria analysis
# Updated: TF=5 pM (CAT-standard), WP1-7 model improvements

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
source("R/43_clinical_outcome_mapping.R")  # WP4: ETP->stroke, HC->bleeding mapping
source("R/50_visualization.R")

# Canonical TF concentration (Hemker CAT-standard, manuscript §2.3.3)
TF_pM_RUN <- 5

# --- Configuration ---
weights   <- c(50, 70, 90)       # kg
ages      <- c(40, 65, 80)       # years
ckd_stages <- c("Normal", "NS_Moderate", "NS_Severe")
egfr_map   <- c(Normal = 100, NS_Moderate = 90, NS_Severe = 90)

dose_ref <- 5  # mg BID reference dose

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Clear caches
clear_baseline_cache()

# =====================================================================
# Part 1: PK metrics at 5 mg across BW x age x NS grid
# =====================================================================
cat("=== Patient Factors: PK Analysis ===\n")
cat("Grid: BW =", paste(weights, collapse = ", "), "kg\n")
cat("      Age =", paste(ages, collapse = ", "), "yr\n")
cat("      NS severity =", paste(ckd_stages, collapse = ", "), "\n\n")

pk_results <- data.frame(
  body_weight = numeric(), age = numeric(), ckd_stage = character(),
  egfr = numeric(), dose_mg = numeric(),
  Cmax_total = numeric(), Ctrough_total = numeric(),
  AUCtau_total = numeric(), Cmax_free = numeric(),
  Ctrough_free = numeric(), AUCtau_free = numeric(),
  stringsAsFactors = FALSE
)

for (bw in weights) {
  for (ag in ages) {
    for (stage in ckd_stages) {
      egfr <- egfr_map[stage]
      cat(sprintf("  BW=%d kg, Age=%d yr, %s (eGFR=%d)... ", bw, ag, stage, egfr))

      tryCatch({
        params <- get_pbpk_params(egfr = egfr, body_weight = bw,
                                  ckd_stage = stage, age = ag)
        sim <- simulate_pbpk(dose_mg = dose_ref, tau = 12, n_doses = 14,
                             egfr = egfr, params = params)
        pk <- compute_pk_metrics(sim, tau = 12)

        pk_results <- rbind(pk_results, data.frame(
          body_weight  = bw, age = ag, ckd_stage = stage,
          egfr = egfr, dose_mg = dose_ref,
          Cmax_total = pk$Cmax_total, Ctrough_total = pk$Ctrough_total,
          AUCtau_total = pk$AUCtau_total, Cmax_free = pk$Cmax_free,
          Ctrough_free = pk$Ctrough_free, AUCtau_free = pk$AUCtau_free,
          stringsAsFactors = FALSE
        ))
        cat(sprintf("AUC=%.1f ng*h/mL\n", pk$AUCtau_total))
      }, error = function(e) {
        cat(sprintf("FAILED: %s\n", e$message))
      })
    }
  }
}

write.csv(pk_results, "output/tables/patient_factors_pk.csv", row.names = FALSE)
cat("\nPK results saved to output/tables/patient_factors_pk.csv\n")

# =====================================================================
# Part 2: Dose optimization across grid (uncoupled TI)
# =====================================================================
cat("\n=== Patient Factors: Dose Optimization ===\n")

ti_objective_pf <- function(dose_mg, egfr, ckd_stage, body_weight, age) {
  res <- run_linked_simulation(
    dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
    TF_pM = TF_pM_RUN,
    body_weight = body_weight, age = age, compute_ti = TRUE
  )
  res$therapeutic_index$TI
}

opt_results <- data.frame(
  body_weight = numeric(), age = numeric(), ckd_stage = character(),
  egfr = numeric(), optimal_dose = numeric(), TI = numeric(),
  efficacy_score = numeric(), safety_score = numeric(),
  ETP = numeric(), HC = numeric(),
  stringsAsFactors = FALSE
)

for (bw in weights) {
  for (ag in ages) {
    for (stage in ckd_stages) {
      egfr <- egfr_map[stage]
      cat(sprintf("  Optimizing BW=%d, Age=%d, %s... ", bw, ag, stage))

      tryCatch({
        opt <- optimize(
          ti_objective_pf,
          interval = c(0.1, 20),
          egfr = egfr, ckd_stage = stage,
          body_weight = bw, age = ag,
          maximum = TRUE, tol = 0.05
        )

        # Verify at optimal dose
        verify <- run_linked_simulation(
          dose_mg = opt$maximum, egfr = egfr, ckd_stage = stage,
          TF_pM = TF_pM_RUN,
          body_weight = bw, age = ag, compute_ti = TRUE
        )

        opt_results <- rbind(opt_results, data.frame(
          body_weight    = bw, age = ag, ckd_stage = stage,
          egfr           = egfr,
          optimal_dose   = round(opt$maximum, 2),
          TI             = verify$therapeutic_index$TI,
          efficacy_score = verify$therapeutic_index$efficacy_score,
          safety_score   = verify$therapeutic_index$safety_score,
          ETP            = verify$tga_metrics$ETP,
          HC             = verify$hemostatic_metrics$HC,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("%.2f mg (TI=%.3f)\n", opt$maximum, verify$therapeutic_index$TI))
      }, error = function(e) {
        cat(sprintf("FAILED: %s\n", e$message))
        opt_results <<- rbind(opt_results, data.frame(
          body_weight = bw, age = ag, ckd_stage = stage, egfr = egfr,
          optimal_dose = NA, TI = NA, efficacy_score = NA,
          safety_score = NA, ETP = NA, HC = NA,
          stringsAsFactors = FALSE
        ))
      })
    }
  }
}

write.csv(opt_results, "output/tables/patient_factors_optimal_doses.csv", row.names = FALSE)
cat("\nOptimal doses saved to output/tables/patient_factors_optimal_doses.csv\n")

# =====================================================================
# Part 3: FDA 2-of-3 criteria flagging
# =====================================================================
cat("\n=== FDA Dose-Reduction Criteria (2 of 3) ===\n")
cat("Criteria: Age >= 80, BW <= 60 kg, SCr >= 1.5 mg/dL\n")
cat("(SCr >= 1.5 approximated as eGFR <= 25; NS defaults preserve eGFR)\n\n")

opt_results$flag_age  <- opt_results$age >= 80
opt_results$flag_bw   <- opt_results$body_weight <= 60
opt_results$flag_scr  <- opt_results$egfr <= 25
opt_results$criteria_met <- opt_results$flag_age + opt_results$flag_bw + opt_results$flag_scr
opt_results$fda_reduce <- opt_results$criteria_met >= 2

flagged <- opt_results[opt_results$fda_reduce, ]
if (nrow(flagged) > 0) {
  cat("Patients meeting 2-of-3 criteria (FDA recommends 2.5 mg BID):\n")
  for (i in seq_len(nrow(flagged))) {
    cat(sprintf("  BW=%d kg, Age=%d yr, %s -> model optimal=%.2f mg\n",
                flagged$body_weight[i], flagged$age[i], flagged$ckd_stage[i],
                flagged$optimal_dose[i]))
  }
} else {
  cat("No patients in current grid meet 2-of-3 criteria.\n")
}

# =====================================================================
# Part 4: Figures
# =====================================================================
cat("\n=== Generating figures ===\n")

# 4a: Dose heatmap (BW x Age, faceted by NS severity)
opt_results$bw_label  <- paste0(opt_results$body_weight, " kg")
opt_results$age_label <- paste0(opt_results$age, " yr")

p_heatmap <- ggplot(opt_results,
       aes(x = age_label, y = bw_label, fill = optimal_dose)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", optimal_dose)), size = 4) +
  facet_wrap(~ ckd_stage, ncol = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 1.0, name = "Optimal\nDose (mg)") +
  labs(x = "Age", y = "Body Weight",
       title = sprintf("TI-Optimal Dose by Body Weight, Age, and NS Severity (TF=%d pM)", TF_pM_RUN),
       subtitle = "Uncoupled model; blue = lower dose, red = higher dose") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        panel.grid = element_blank())

ggsave("output/figures/patient_factors_dose_heatmap.png", p_heatmap,
       width = 12, height = 5, dpi = 300)
cat("  Saved: patient_factors_dose_heatmap.png\n")

# 4b: AUC comparison bar chart
pk_results$label <- paste0(pk_results$body_weight, "kg/", pk_results$age, "yr")
pk_results$ckd_stage <- factor(pk_results$ckd_stage, levels = ckd_stages)

p_auc <- ggplot(pk_results,
       aes(x = label, y = AUCtau_total, fill = ckd_stage)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(x = "Patient (BW/Age)", y = "AUCtau (ng*h/mL)",
       fill = "NS Severity",
       title = sprintf("Steady-State AUC at %d mg BID by Patient Factors", dose_ref),
       subtitle = sprintf("TF=%d pM | Higher free exposure in hypoalbuminemic/proteinuric NS", TF_pM_RUN)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("output/figures/patient_factors_auc_comparison.png", p_auc,
       width = 10, height = 6, dpi = 300)
cat("  Saved: patient_factors_auc_comparison.png\n")

# =====================================================================
# Part 5: WP4 Clinical Outcome Mapping
# =====================================================================
cat("\n=== WP4: Clinical Outcome Mapping ===\n")
cat("Computing stroke/bleeding hazard rates at each patient's optimal dose...\n\n")

# Get healthy baseline HC for bleeding hazard calibration
clear_baseline_cache()
ic_healthy <- get_coag_initial_conditions(TF_pM = TF_pM_RUN)
coag_healthy <- simulate_coagulation(ic = ic_healthy, Cp_free_nM = 0, TF_pM = TF_pM_RUN)
plt_healthy <- simulate_platelets(coag_healthy)
hc_healthy <- compute_hemostatic_metrics(plt_healthy)$HC

opt_results$stroke_rate_yr <- NA_real_
opt_results$bleed_rate_yr  <- NA_real_
opt_results$utility_eq     <- NA_real_

for (i in seq_len(nrow(opt_results))) {
  if (!is.na(opt_results$ETP[i]) && !is.na(opt_results$HC[i])) {
    opt_results$stroke_rate_yr[i] <- etp_to_stroke_hazard(opt_results$ETP[i])
    opt_results$bleed_rate_yr[i]  <- hc_to_bleeding_hazard(opt_results$HC[i], hc_healthy)
    opt_results$utility_eq[i]     <- clinical_utility(opt_results$ETP[i], opt_results$HC[i],
                                                       hc_healthy_nodrug = hc_healthy)
  }
}

write.csv(opt_results, "output/tables/patient_factors_optimal_doses.csv", row.names = FALSE)
cat("Updated optimal doses with clinical outcome columns\n")

# 5a: Outcome trade-off heatmap
if (any(!is.na(opt_results$utility_eq))) {
  p_utility <- ggplot(opt_results,
         aes(x = age_label, y = bw_label, fill = utility_eq * 100)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", utility_eq * 100)), size = 3.5) +
    facet_wrap(~ ckd_stage, ncol = 3) +
    scale_fill_gradient(low = "#F7FCB9", high = "#D7301F",
                        name = "Combined\nEvent Rate\n(%/yr)") +
    labs(x = "Age", y = "Body Weight",
         title = sprintf("Annual Stroke+Bleeding Rate at Optimal Dose (TF=%d pM)", TF_pM_RUN),
         subtitle = "WP4: ETP->stroke (Besser 2008), HC->bleeding (Tantry 2013)") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey90"),
          panel.grid = element_blank())

  ggsave("output/figures/patient_factors_outcome_heatmap.png", p_utility,
         width = 12, height = 5, dpi = 300)
  cat("  Saved: patient_factors_outcome_heatmap.png\n")
}

# =====================================================================
# Summary
# =====================================================================
cat("\n=== Summary ===\n")
cat(sprintf("Model: TF=%d pM, dose interval [0.1, 20] mg\n", TF_pM_RUN))
cat(sprintf("Total scenarios: %d\n", nrow(opt_results)))
cat(sprintf("Dose range: %.2f - %.2f mg\n",
            min(opt_results$optimal_dose, na.rm = TRUE),
            max(opt_results$optimal_dose, na.rm = TRUE)))
cat(sprintf("FDA 2-of-3 flagged: %d / %d\n",
            sum(opt_results$fda_reduce, na.rm = TRUE), nrow(opt_results)))
if (any(!is.na(opt_results$utility_eq))) {
  cat(sprintf("Combined event rate range: %.2f%% - %.2f%% per year\n",
              min(opt_results$utility_eq, na.rm = TRUE) * 100,
              max(opt_results$utility_eq, na.rm = TRUE) * 100))
}
cat("\nAll outputs saved to output/tables/ and output/figures/\n")
