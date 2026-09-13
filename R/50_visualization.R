# 50_visualization.R
# Plotting functions for PBPK-QSP results

#' Plot PK profiles by disease stage
#' @param results_list Named list from run_linked_simulation
#' @param ckd_stages Character vector of CKD stage names
#' @param doses Numeric vector of doses
#' @return ggplot object (patchwork)
plot_pk_profiles <- function(results_list, ckd_stages, doses) {
  pk_data <- purrr::map_dfr(names(results_list), function(key) {
    res <- results_list[[key]]
    ss <- extract_steady_state(res$pbpk_result)
    data.frame(
      time_h    = ss$time_in_interval,
      Cp_total  = ss$C_total_ng_mL,
      Cp_free   = ss$C_free_ng_mL,
      dose_mg   = res$dose_mg,
      ckd_stage = res$ckd_stage,
      stringsAsFactors = FALSE
    )
  })

  pk_data$ckd_stage <- factor(pk_data$ckd_stage, levels = ckd_stages)
  pk_data$dose_label <- paste0(pk_data$dose_mg, " mg BID")

  p1 <- ggplot(pk_data, aes(x = time_h, y = Cp_total, color = ckd_stage)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ dose_label) +
    labs(x = "Time in dosing interval (h)",
         y = "Total plasma concentration (ng/mL)",
         color = "NS Severity",
         title = "Steady-state apixaban PK by nephrotic syndrome severity") +
    theme_bw() +
    theme(legend.position = "bottom")

  p2 <- ggplot(pk_data, aes(x = time_h, y = Cp_free, color = ckd_stage)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ dose_label) +
    labs(x = "Time in dosing interval (h)",
         y = "Free plasma concentration (ng/mL)",
         color = "NS Severity",
         title = "Free (unbound) apixaban by nephrotic syndrome severity") +
    theme_bw() +
    theme(legend.position = "bottom")

  p1 / p2
}

#' Plot thrombin generation curves by disease stage
#' @param results_list Named list from run_linked_simulation
#' @param ckd_stages Character vector of CKD stage names
#' @param doses Numeric vector of doses
#' @return ggplot object (patchwork)
plot_tga_curves <- function(results_list, ckd_stages, doses) {
  tga_data <- purrr::map_dfr(names(results_list), function(key) {
    res <- results_list[[key]]
    data.frame(
      time_s    = res$coag_result$time,
      IIa       = res$coag_result$IIa,
      dose_mg   = res$dose_mg,
      ckd_stage = res$ckd_stage,
      stringsAsFactors = FALSE
    )
  })

  tga_data$ckd_stage <- factor(tga_data$ckd_stage, levels = ckd_stages)
  tga_data$dose_label <- paste0(tga_data$dose_mg, " mg BID")
  tga_data$time_min <- tga_data$time_s / 60

  ggplot(tga_data, aes(x = time_min, y = IIa, color = ckd_stage)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ dose_label) +
    labs(x = "Time (min)",
         y = "Thrombin (nM)",
         color = "NS Severity",
         title = "Thrombin generation curves by NS severity and dose") +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Plot TGA metrics heatmap (dose x disease stage)
#' @param summary_df Summary data frame from run_linked_simulation
#' @return ggplot object (patchwork)
plot_tga_heatmap <- function(summary_df) {
  # Normalize metrics to Normal/5mg reference
  ref <- summary_df[summary_df$ckd_stage == "Normal" & summary_df$dose_mg == 5, ]

  summary_df$ETP_ratio <- summary_df$ETP / ref$ETP
  summary_df$peak_ratio <- summary_df$peak_thrombin / ref$peak_thrombin
  summary_df$lag_ratio <- summary_df$lag_time / ref$lag_time

  summary_df$dose_label <- paste0(summary_df$dose_mg, " mg")

  p_etp <- ggplot(summary_df, aes(x = ckd_stage, y = dose_label, fill = ETP_ratio)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", ETP_ratio)), size = 3.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
    labs(title = "ETP (ratio to Normal/5mg)", x = "", y = "", fill = "Ratio") +
    theme_minimal()

  p_peak <- ggplot(summary_df, aes(x = ckd_stage, y = dose_label, fill = peak_ratio)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", peak_ratio)), size = 3.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
    labs(title = "Peak thrombin (ratio to Normal/5mg)", x = "", y = "", fill = "Ratio") +
    theme_minimal()

  p_lag <- ggplot(summary_df, aes(x = ckd_stage, y = dose_label, fill = lag_ratio)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", lag_ratio)), size = 3.5) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1) +
    labs(title = "Lag time (ratio to Normal/5mg)", x = "", y = "", fill = "Ratio") +
    theme_minimal()

  p_etp / p_peak / p_lag
}

#' Plot PD-equivalent doses by disease stage
#' @param dose_df Data frame with columns: ckd_stage, pd_equiv_dose
#' @return ggplot object
plot_pd_equivalent_doses <- function(dose_df) {
  dose_df$ckd_stage <- factor(dose_df$ckd_stage,
                              levels = unique(dose_df$ckd_stage))

  ggplot(dose_df, aes(x = ckd_stage, y = pd_equiv_dose)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = 5.0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "orange") +
    annotate("text", x = 0.5, y = 5.2, label = "Standard 5 mg", hjust = 0, color = "red", size = 3) +
    annotate("text", x = 0.5, y = 2.7, label = "Reduced 2.5 mg", hjust = 0, color = "orange", size = 3) +
    labs(x = "NS Severity", y = "PD-equivalent dose (mg BID)",
         title = "PD-optimized apixaban dose by nephrotic syndrome severity",
         subtitle = "Dose matching ETP of 5 mg BID in normal renal function") +
    theme_bw()
}

#' Plot benefit-risk: efficacy vs safety trade-off + TI
#' @param summary_df Summary data frame with efficacy_score, safety_score, TI columns
#' @return ggplot object (patchwork)
plot_benefit_risk <- function(summary_df) {
  summary_df$dose_label <- paste0(summary_df$dose_mg, " mg")

  # Panel 1: Efficacy vs Safety scatter
  p1 <- ggplot(summary_df, aes(x = safety_score, y = efficacy_score,
                                 color = ckd_stage, shape = dose_label)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Safety score (HC retained)",
         y = "Efficacy score (ETP suppressed)",
         color = "NS Severity", shape = "Dose",
         title = "Benefit-risk trade-off") +
    theme_bw() +
    theme(legend.position = "bottom")

  # Panel 2: TI by CKD stage and dose
  p2 <- ggplot(summary_df, aes(x = ckd_stage, y = TI, fill = dose_label)) +
    geom_col(position = "dodge", alpha = 0.8) +
    labs(x = "NS Severity", y = "Therapeutic Index",
         fill = "Dose",
         title = "Therapeutic index by NS severity and dose") +
    theme_bw() +
    theme(legend.position = "bottom")

  p1 | p2
}

#' Plot TI-optimal vs PD-equivalent dose comparison
#' @param dose_df Data frame with optimal_dose and pd_equiv_dose columns
#' @return ggplot object
plot_dose_comparison <- function(dose_df) {
  dose_df$ckd_stage <- factor(dose_df$ckd_stage,
                              levels = unique(dose_df$ckd_stage))

  plot_data <- tidyr::pivot_longer(
    dose_df[, c("ckd_stage", "optimal_dose", "pd_equiv_dose")],
    cols = c("optimal_dose", "pd_equiv_dose"),
    names_to = "method",
    values_to = "dose_mg"
  )
  plot_data$method <- ifelse(plot_data$method == "optimal_dose",
                             "TI-optimal", "PD-equivalent (ETP)")

  ggplot(plot_data, aes(x = ckd_stage, y = dose_mg, fill = method)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_hline(yintercept = 5.0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "orange") +
    annotate("text", x = 0.5, y = 5.2, label = "Standard 5 mg",
             hjust = 0, color = "red", size = 3) +
    annotate("text", x = 0.5, y = 2.7, label = "Reduced 2.5 mg",
             hjust = 0, color = "orange", size = 3) +
    scale_fill_manual(values = c("TI-optimal" = "steelblue",
                                  "PD-equivalent (ETP)" = "coral")) +
    labs(x = "NS Severity", y = "Dose (mg BID)", fill = "Method",
         title = "Dose optimization: TI-optimal vs PD-equivalent",
         subtitle = "TI-optimal accounts for bleeding risk; PD-equiv matches ETP only") +
    theme_bw() +
    theme(legend.position = "bottom")
}
