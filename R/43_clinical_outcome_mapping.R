# 43_clinical_outcome_mapping.R
# WP4: Map ETP and HC surrogates to clinical stroke/bleeding hazard rates
# Based on literature-derived associations from data/etp_outcome_associations.csv

#' Compute annual stroke hazard rate from ETP
#' Based on Besser 2008: HR 1.30 per 100 nM*min ETP increase
#' and baseline AF stroke rate from ARISTOTLE: ~1.5%/yr on apixaban
#' @param etp_nM_min ETP in nM*min
#' @param baseline_rate Annual baseline stroke rate (default 0.015)
#' @param hr_per_100 Hazard ratio per 100 nM*min ETP increase (default 1.30)
#' @param etp_reference Reference ETP for HR=1 (default 400 nM*min, approximate healthy)
#' @return Annual stroke hazard rate
etp_to_stroke_hazard <- function(etp_nM_min, baseline_rate = 0.015,
                                  hr_per_100 = 1.30, etp_reference = 400) {
  stop("ETP-to-stroke mapping is unvalidated and disabled. Besser 2008 studied VTE recurrence, not stroke.")
  # Log-linear model: hazard = baseline * HR^(delta_ETP / 100)
  delta_etp <- etp_nM_min - etp_reference
  hr <- hr_per_100 ^ (delta_etp / 100)
  baseline_rate * hr
}

#' Compute annual major bleeding hazard rate from hemostatic capacity
#' Inverse mapping: lower HC → higher bleeding risk
#' Based on platelet reactivity data (Tantry 2013) adapted to HC metric:
#' HR ~2.0 when HC falls below 50% of healthy baseline
#' @param hc Hemostatic capacity (nM*s)
#' @param hc_healthy_nodrug HC for healthy patient without drug (~494,000 nM*s)
#' @param baseline_rate Annual baseline major bleeding rate (default 0.025 for AF on DOAC)
#' @param hr_at_50pct HR when HC is at 50% of healthy (default 2.0)
#' @return Annual major bleeding hazard rate
hc_to_bleeding_hazard <- function(hc, hc_healthy_nodrug = 494000,
                                   baseline_rate = 0.025,
                                   hr_at_50pct = 2.0) {
  warning("Hypothetical mathematical mapping only; not a validated clinical bleeding hazard.")
  # HC fraction retained
  hc_fraction <- hc / hc_healthy_nodrug
  hc_fraction <- max(0.01, min(1.0, hc_fraction))  # clamp

  # Log-linear model calibrated to: at 50% HC, HR = hr_at_50pct
  # log(HR) = -k * log(hc_fraction), where k = log(hr_at_50pct) / log(0.5)
  k <- log(hr_at_50pct) / log(0.5)
  hr <- hc_fraction ^ k

  baseline_rate * hr
}

#' Compute clinical utility score (lower = better)
#' U(dose) = w_stroke * stroke_rate + w_bleed * bleed_rate
#' @param etp_nM_min ETP at the given dose
#' @param hc HC at the given dose
#' @param w_stroke Weight for stroke events (default 1.0)
#' @param w_bleed Weight for bleeding events (default 1.0)
#' @param hc_healthy_nodrug Healthy baseline HC for bleeding calibration
#' @return Clinical utility (annual combined event rate, weighted)
clinical_utility <- function(etp_nM_min, hc,
                              w_stroke = 1.0, w_bleed = 1.0,
                              hc_healthy_nodrug = 494000) {
  stroke_rate  <- etp_to_stroke_hazard(etp_nM_min)
  bleed_rate   <- hc_to_bleeding_hazard(hc, hc_healthy_nodrug)
  w_stroke * stroke_rate + w_bleed * bleed_rate
}
