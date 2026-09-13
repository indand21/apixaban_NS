stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_linked_simulation.R
# Main script: Run linked PBPK-QSP across nephrotic syndrome severity groups and doses
# Outputs combined TGA metrics and PK metrics tables

# --- Setup ---
rm(list = ls())
setwd(here::here())  # Or set to project root manually

# Source all model files
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
source("R/51_tables.R")

# --- Scenario Definition ---
ckd_stages <- c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")
egfr_values <- c(Normal = 100, NS_Mild = 90, NS_Moderate = 90, NS_Severe = 90)
doses <- c(2.5, 5.0)  # mg BID

# --- Run All Scenarios ---
cat("Running linked PBPK-QSP simulations across nephrotic syndrome severity groups and doses...\n")
cat("=================================================================\n\n")

results_list <- list()
summary_rows <- list()
counter <- 0

for (dose_mg in doses) {
  for (stage in ckd_stages) {
    counter <- counter + 1
    egfr <- egfr_values[stage]

    cat(sprintf("[%d/%d] Dose: %.1f mg BID | NS: %s (eGFR=%d)...",
                counter, length(doses) * length(ckd_stages),
                dose_mg, stage, egfr))

    # Run linked simulation
    res <- run_linked_simulation(
      dose_mg   = dose_mg,
      egfr      = egfr,
      ckd_stage = stage,
      use_trough = TRUE,
      TF_pM     = 25,
      tau       = 12,
      n_doses   = 14
    )

    results_list[[paste0(dose_mg, "_", stage)]] <- res

    # Compile summary row
    ti <- res$therapeutic_index
    hm <- res$hemostatic_metrics
    summary_rows[[counter]] <- data.frame(
      dose_mg        = dose_mg,
      ckd_stage      = stage,
      egfr           = egfr,
      Cmax_total     = res$pk_metrics$Cmax_total,
      Ctrough_total  = res$pk_metrics$Ctrough_total,
      AUCtau_total   = res$pk_metrics$AUCtau_total,
      Cp_free_nM     = res$Cp_free_nM,
      ETP            = res$tga_metrics$ETP,
      ETP_nM_min     = res$tga_metrics$ETP_nM_min,
      peak_thrombin  = res$tga_metrics$peak_thrombin,
      lag_time       = res$tga_metrics$lag_time,
      time_to_peak   = res$tga_metrics$time_to_peak,
      velocity_index = res$tga_metrics$velocity_index,
      HC             = hm$HC,
      HC_nM_min      = hm$HC_nM_min,
      peak_PLT_agg   = hm$peak_PLT_agg,
      hemostatic_lag = hm$hemostatic_lag,
      efficacy_score = if (!is.null(ti)) ti$efficacy_score else NA,
      safety_score   = if (!is.null(ti)) ti$safety_score else NA,
      TI             = if (!is.null(ti)) ti$TI else NA,
      stringsAsFactors = FALSE
    )

    cat(" done\n")
  }
}

# --- Combine Results ---
summary_df <- do.call(rbind, summary_rows)
summary_df$ckd_stage <- factor(summary_df$ckd_stage, levels = ckd_stages)

cat("\n=== SUMMARY TABLE ===\n\n")
print(summary_df, digits = 3)

# --- Save outputs ---
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

write.csv(summary_df, "output/tables/linked_simulation_summary.csv", row.names = FALSE)
cat("\nSummary saved to output/tables/linked_simulation_summary.csv\n")

# --- Generate Plots ---
cat("\nGenerating plots...\n")

# PK profiles
pk_plot <- plot_pk_profiles(results_list, ckd_stages, doses)
ggsave("output/figures/pk_profiles_by_ns.png", pk_plot, width = 12, height = 8, dpi = 300)

# TGA curves
tga_plot <- plot_tga_curves(results_list, ckd_stages, doses)
ggsave("output/figures/tga_curves_by_ns.png", tga_plot, width = 12, height = 8, dpi = 300)

# TGA metrics heatmap
heatmap_plot <- plot_tga_heatmap(summary_df)
ggsave("output/figures/tga_metrics_heatmap.png", heatmap_plot, width = 10, height = 8, dpi = 300)

# Benefit-risk plot
if ("TI" %in% names(summary_df)) {
  br_plot <- plot_benefit_risk(summary_df)
  ggsave("output/figures/benefit_risk.png", br_plot, width = 12, height = 6, dpi = 300)
}

# Summary table
make_summary_table(summary_df, "output/tables/formatted_summary.csv")

cat("\nAll uncoupled simulations complete.\n")

# =================================================================
# COUPLED (Two-Way) Simulations: PS_exp -> Prothrombinase Feedback
# =================================================================
cat("\n=================================================================\n")
cat("Running COUPLED (two-way) simulations: PS_exp -> prothrombinase...\n")
cat("=================================================================\n\n")

coupled_results_list <- list()
coupled_summary_rows <- list()
counter <- 0

for (dose_mg in doses) {
  for (stage in ckd_stages) {
    counter <- counter + 1
    egfr <- egfr_values[stage]

    cat(sprintf("[%d/%d] COUPLED | Dose: %.1f mg BID | NS: %s (eGFR=%d)...",
                counter, length(doses) * length(ckd_stages),
                dose_mg, stage, egfr))

    res_c <- run_coupled_linked_simulation(
      dose_mg    = dose_mg,
      egfr       = egfr,
      ckd_stage  = stage,
      use_trough = TRUE,
      TF_pM      = 25,
      tau        = 12,
      n_doses    = 14
    )

    coupled_results_list[[paste0(dose_mg, "_", stage)]] <- res_c

    ti <- res_c$therapeutic_index
    hm <- res_c$hemostatic_metrics
    coupled_summary_rows[[counter]] <- data.frame(
      mode           = "coupled",
      dose_mg        = dose_mg,
      ckd_stage      = stage,
      egfr           = egfr,
      Cp_free_nM     = res_c$Cp_free_nM,
      ETP            = res_c$tga_metrics$ETP,
      ETP_nM_min     = res_c$tga_metrics$ETP_nM_min,
      peak_thrombin  = res_c$tga_metrics$peak_thrombin,
      lag_time       = res_c$tga_metrics$lag_time,
      time_to_peak   = res_c$tga_metrics$time_to_peak,
      HC             = hm$HC,
      HC_nM_min      = hm$HC_nM_min,
      peak_PLT_agg   = hm$peak_PLT_agg,
      efficacy_score = if (!is.null(ti)) ti$efficacy_score else NA,
      safety_score   = if (!is.null(ti)) ti$safety_score else NA,
      TI             = if (!is.null(ti)) ti$TI else NA,
      stringsAsFactors = FALSE
    )

    cat(" done\n")
  }
}

coupled_summary_df <- do.call(rbind, coupled_summary_rows)
coupled_summary_df$ckd_stage <- factor(coupled_summary_df$ckd_stage, levels = ckd_stages)

cat("\n=== COUPLED SUMMARY TABLE ===\n\n")
print(coupled_summary_df, digits = 3)

write.csv(coupled_summary_df, "output/tables/coupled_simulation_summary.csv", row.names = FALSE)
cat("\nCoupled summary saved to output/tables/coupled_simulation_summary.csv\n")

# --- Comparison: Coupled vs Uncoupled ---
cat("\n=== COUPLED vs UNCOUPLED COMPARISON ===\n")
cat("(Showing peak thrombin ratios: coupled/uncoupled)\n\n")

for (dose_mg in doses) {
  cat(sprintf("Dose: %.1f mg BID\n", dose_mg))
  for (stage in ckd_stages) {
    key <- paste0(dose_mg, "_", stage)
    peak_unc <- results_list[[key]]$tga_metrics$peak_thrombin
    peak_cpl <- coupled_results_list[[key]]$tga_metrics$peak_thrombin
    ratio <- peak_cpl / peak_unc
    cat(sprintf("  %6s: uncoupled=%.1f nM, coupled=%.1f nM, ratio=%.3f\n",
                stage, peak_unc, peak_cpl, ratio))
  }
}

cat("\nAll simulations complete. Results in output/\n")
