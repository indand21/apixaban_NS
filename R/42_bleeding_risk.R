# 42_bleeding_risk.R
# Exploratory platelet aggregate metrics and bounded composite score

#' Compute hemostatic capacity metrics from platelet simulation
#' @param plt_result Data frame from simulate_platelets
#' @param agg_threshold Threshold in model-equivalent platelet units (default 10)
#' @return Named list: HC, peak_PLT_agg, hemostatic_lag
compute_hemostatic_metrics <- function(plt_result, agg_threshold = 10) {
  time_s  <- plt_result$time
  PLT_agg <- plt_result$PLT_agg

  # HC: integral of PLT_agg(t) in model-equivalent units times seconds
  HC <- compute_auc(time_s, PLT_agg)

  # Legacy HC_nM_min field is actually model-equivalent units times minutes
  HC_nM_min <- HC / 60

  # Peak aggregate state (model-equivalent units)
  peak_idx <- which.max(PLT_agg)
  peak_PLT_agg <- PLT_agg[peak_idx]
  time_to_peak_agg <- time_s[peak_idx]

  # Hemostatic lag: time to PLT_agg > threshold
  above_thresh <- which(PLT_agg >= agg_threshold)
  if (length(above_thresh) > 0) {
    idx <- above_thresh[1]
    if (idx > 1) {
      t1 <- time_s[idx - 1]; t2 <- time_s[idx]
      c1 <- PLT_agg[idx - 1]; c2 <- PLT_agg[idx]
      hemostatic_lag <- t1 + (agg_threshold - c1) * (t2 - t1) / (c2 - c1)
    } else {
      hemostatic_lag <- time_s[1]
    }
  } else {
    hemostatic_lag <- NA
  }

  list(
    HC               = HC,
    HC_nM_min        = HC_nM_min,
    peak_PLT_agg     = peak_PLT_agg,
    time_to_peak_agg = time_to_peak_agg,
    hemostatic_lag   = hemostatic_lag
  )
}

#' Compute exploratory composite score from drug vs no-drug metrics
#'
#' Both denominators use the matched individual and feedback mode without drug.
#' Historical field names efficacy_score, safety_score and TI are compatibility
#' aliases only. These quantities are not validated clinical efficacy or safety.
#'
#' @section Choice of efficacy endpoint:
#' The primary efficacy endpoint is PEAK THROMBIN suppression, not ETP
#' suppression. At therapeutic trough exposure the ETP of this model changes by
#' about 0.01 percent, because ETP integrates a thrombin curve whose area is
#' dominated by the slow termination phase rather than by the drug-sensitive
#' propagation burst. An objective built on that quantity is numerically flat,
#' so the dose search returned search-boundary solutions rather than interior
#' maxima. Peak thrombin moves by tens of percent over the same exposure range
#' and is the drug-sensitive summary of the same simulation.
#' The ETP-based score is retained and reported alongside as a sensitivity
#' (fields efficacy_score_etp and TI_etp) so the change is auditable and the
#' previous behaviour remains inspectable.
#'
#' @param tga_drug TGA metrics with drug (from compute_tga_metrics)
#' @param tga_nodrug TGA metrics without drug (matched no-drug baseline)
#' @param hemo_drug Hemostatic metrics with drug (from compute_hemostatic_metrics)
#' @param hemo_nodrug Hemostatic metrics without drug (matched no-drug baseline)
#' @param hemo_healthy_nodrug Optional healthy no-drug reference. Retained for
#'        interface compatibility; not used for the reported score.
#' @param objective Efficacy endpoint driving the reported score: "peak"
#'        (default, peak-thrombin suppression) or "etp" (legacy ETP suppression).
#' @return Named list: efficacy_score, safety_score, TI, plus the endpoint
#'         variants efficacy_score_peak, efficacy_score_etp, TI_peak, TI_etp
#'         and the label objective.
compute_therapeutic_index <- function(tga_drug, tga_nodrug,
                                      hemo_drug, hemo_nodrug,
                                      hemo_healthy_nodrug = NULL,
                                      objective = c("peak", "etp")) {

  objective <- match.arg(objective)

  # Bounded suppression fraction: 1 - drug/no-drug, clamped to [0, 1].
  # A missing or non-finite endpoint gives NA (the caller decides whether that
  # is fatal); a no-drug baseline of zero means there is nothing to suppress,
  # which is a well-defined score of zero rather than an error.
  suppression <- function(drug, nodrug) {
    if (is.null(drug) || is.null(nodrug) ||
        !is.finite(drug) || !is.finite(nodrug)) return(NA_real_)
    if (nodrug <= 0) return(0)
    max(0, min(1, 1 - drug / nodrug))
  }

  efficacy_etp  <- suppression(tga_drug$ETP,           tga_nodrug$ETP)
  efficacy_peak <- suppression(tga_drug$peak_thrombin, tga_nodrug$peak_thrombin)

  if (is.na(efficacy_etp)) efficacy_etp <- 0

  if (objective == "peak" && is.na(efficacy_peak)) {
    stop("objective = 'peak' requires peak_thrombin in both TGA metric lists. ",
         "Use compute_tga_metrics(), or pass objective = 'etp'.")
  }

  # Matched platelet aggregate retention; not an observed bleeding risk
  ref_HC <- hemo_nodrug$HC  # Matched individual no-drug reference
  if (is.finite(ref_HC) && ref_HC > 0) {
    safety_score <- hemo_drug$HC / ref_HC
  } else {
    safety_score <- 0
  }
  safety_score <- max(0, min(1, safety_score))

  # Exploratory composite scores, not therapeutic indices or validated safety metrics
  TI_peak <- if (is.na(efficacy_peak)) NA_real_ else efficacy_peak * safety_score
  TI_etp  <- efficacy_etp * safety_score

  efficacy_score <- if (objective == "peak") efficacy_peak else efficacy_etp
  TI             <- if (objective == "peak") TI_peak       else TI_etp

  list(
    efficacy_score      = efficacy_score,
    safety_score        = safety_score,
    TI                  = TI,
    objective           = objective,
    efficacy_score_peak = efficacy_peak,
    efficacy_score_etp  = efficacy_etp,
    TI_peak             = TI_peak,
    TI_etp              = TI_etp
  )
}
