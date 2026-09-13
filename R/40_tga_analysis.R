# 40_tga_analysis.R
# Thrombin Generation Assay (TGA) metric extraction

#' Compute TGA metrics from coagulation simulation
#' @param coag_result Data frame from simulate_coagulation
#' @param thrombin_col Column name for thrombin concentration (default "IIa")
#' @param lag_threshold Thrombin threshold for lag time (nM, default 10)
#' @return Named list: ETP, peak_thrombin, lag_time, time_to_peak
compute_tga_metrics <- function(coag_result,
                                thrombin_col = "IIa",
                                lag_threshold = 10) {

  time_s  <- coag_result$time
  thrombin <- coag_result[[thrombin_col]]

  # ETP: Endogenous Thrombin Potential (area under thrombin curve, nM*s)
  ETP <- compute_auc(time_s, thrombin)

  # Convert ETP to nM*min for clinical comparison
  ETP_nM_min <- ETP / 60

  # Peak thrombin (nM)
  peak_idx <- which.max(thrombin)
  peak_thrombin <- thrombin[peak_idx]

  # Time-to-peak (seconds)
  time_to_peak <- time_s[peak_idx]

  # Lag time: time to first exceed threshold
  above_thresh <- which(thrombin >= lag_threshold)
  if (length(above_thresh) > 0) {
    # Linear interpolation for more precise lag time
    idx <- above_thresh[1]
    if (idx > 1) {
      t1 <- time_s[idx - 1]; t2 <- time_s[idx]
      c1 <- thrombin[idx - 1]; c2 <- thrombin[idx]
      lag_time <- t1 + (lag_threshold - c1) * (t2 - t1) / (c2 - c1)
    } else {
      lag_time <- time_s[1]
    }
  } else {
    lag_time <- NA  # Threshold never reached
  }

  # Velocity index: peak / (time_to_peak - lag_time)
  if (!is.na(lag_time) && time_to_peak > lag_time) {
    velocity_index <- peak_thrombin / (time_to_peak - lag_time)
  } else {
    velocity_index <- NA
  }

  list(
    ETP            = ETP,
    ETP_nM_min     = ETP_nM_min,
    peak_thrombin  = peak_thrombin,
    time_to_peak   = time_to_peak,
    lag_time       = lag_time,
    velocity_index = velocity_index
  )
}

#' Compute TGA metrics as a one-row data frame (for binding)
#' @param coag_result Data frame from simulate_coagulation
#' @param ... Additional columns to include (e.g., ckd_stage, dose_mg)
#' @return One-row data frame of TGA metrics
compute_tga_metrics_df <- function(coag_result, ...) {
  metrics <- compute_tga_metrics(coag_result)
  extra <- list(...)
  data.frame(c(metrics, extra), stringsAsFactors = FALSE)
}
