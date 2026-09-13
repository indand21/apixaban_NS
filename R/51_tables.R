# 51_tables.R
# Summary table generation

#' Create formatted summary table
#' @param summary_df Data frame from run_linked_simulation
#' @param output_path File path for CSV output
#' @return Data frame (invisibly)
make_summary_table <- function(summary_df, output_path = NULL) {
  fmt <- summary_df %>%
    mutate(
      Cmax_total    = round(Cmax_total, 1),
      Ctrough_total = round(Ctrough_total, 1),
      AUCtau_total  = round(AUCtau_total, 0),
      Cp_free_nM    = round(Cp_free_nM, 2),
      ETP_nM_min    = round(ETP_nM_min, 0),
      peak_thrombin = round(peak_thrombin, 1),
      lag_time_s    = round(lag_time, 1),
      time_to_peak_s = round(time_to_peak, 1),
      velocity_index = round(velocity_index, 2)
    )

  # Add HC and TI columns if present
  if ("HC" %in% names(summary_df)) {
    fmt <- fmt %>%
      mutate(
        HC_nM_min      = round(HC_nM_min, 0),
        peak_PLT_agg   = round(peak_PLT_agg, 1),
        hemostatic_lag_s = round(hemostatic_lag, 1)
      )
  }
  if ("TI" %in% names(summary_df)) {
    fmt <- fmt %>%
      mutate(
        efficacy_score = round(efficacy_score, 3),
        safety_score   = round(safety_score, 3),
        TI             = round(TI, 3)
      )
  }

  # Select columns based on availability
  base_cols <- c("dose_mg", "ckd_stage", "egfr",
                 "Cmax_total", "Ctrough_total", "AUCtau_total",
                 "Cp_free_nM", "ETP_nM_min", "peak_thrombin",
                 "lag_time_s", "time_to_peak_s", "velocity_index")
  hc_cols <- c("HC_nM_min", "peak_PLT_agg", "hemostatic_lag_s")
  ti_cols <- c("efficacy_score", "safety_score", "TI")

  sel_cols <- base_cols
  if ("HC" %in% names(summary_df)) sel_cols <- c(sel_cols, hc_cols)
  if ("TI" %in% names(summary_df)) sel_cols <- c(sel_cols, ti_cols)

  fmt <- fmt %>% select(all_of(sel_cols))

  if (!is.null(output_path)) {
    write.csv(fmt, output_path, row.names = FALSE)
    cat("Formatted table saved to:", output_path, "\n")
  }

  invisible(fmt)
}

#' Create PK comparison table (observed vs predicted)
#' @param pk_df Data frame with predicted PK metrics
#' @param observed Named list of observed values (e.g., list(Cmax=171, Ctrough=51))
#' @return Data frame with comparison
make_pk_comparison_table <- function(pk_df, observed = NULL) {
  if (is.null(observed)) {
    # Published apixaban 5 mg BID steady-state PK (normal renal function)
    observed <- data.frame(
      metric   = c("Cmax_total", "Ctrough_total", "AUCtau_total"),
      observed = c(171, 51, 1100),  # ng/mL, ng/mL, ng*h/mL
      source   = c("FDA label", "FDA label", "FDA label"),
      stringsAsFactors = FALSE
    )
  }

  comparison <- pk_df %>%
    filter(ckd_stage == "Normal", dose_mg == 5) %>%
    select(Cmax_total, Ctrough_total, AUCtau_total) %>%
    tidyr::pivot_longer(everything(), names_to = "metric", values_to = "predicted") %>%
    left_join(observed, by = "metric") %>%
    mutate(
      ratio = round(predicted / observed, 3),
      pct_error = round((predicted - observed) / observed * 100, 1)
    )

  comparison
}
