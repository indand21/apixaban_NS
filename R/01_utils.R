# 01_utils.R
# Shared utility functions

#' Compute AUC using trapezoidal rule
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @return Scalar AUC value
compute_auc <- function(time, conc) {
  pracma::trapz(time, conc)
}

#' Extract Cmax and Tmax from a concentration-time profile
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @return Named list with Cmax and Tmax
compute_cmax <- function(time, conc) {
  idx <- which.max(conc)
  list(Cmax = conc[idx], Tmax = time[idx])
}

#' Extract Ctrough (concentration at end of dosing interval)
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @param tau Dosing interval in hours (default 12 for BID)
#' @return Ctrough value
compute_ctrough <- function(time, conc, tau = 12) {
  # Find closest time to end of last complete interval
  n_intervals <- floor(max(time) / tau)
  t_trough <- n_intervals * tau
  idx <- which.min(abs(time - t_trough))
  conc[idx]
}

#' Convert eGFR to CKD stage
#' @param egfr eGFR in mL/min/1.73m2
#' @return Character CKD stage label
egfr_to_ckd_stage <- function(egfr) {
  if (egfr >= 90) return("Normal")
  if (egfr >= 60) return("CKD2")
  if (egfr >= 30) return("CKD3")
  if (egfr >= 15) return("CKD4")
  if (egfr > 0)   return("CKD5")
  return("CKD5D")
}

#' Create dosing events table for deSolve
#' @param dose_mg Dose in mg
#' @param tau Dosing interval in hours (default 12 for BID)
#' @param n_doses Number of doses
#' @param bioavailability Oral bioavailability (F)
#' @param lag_h Absorption lag time in hours. Drug enters the gut compartment
#'   this long after ingestion, so the dosing interval still begins at
#'   ingestion and Tmax is measured from ingestion, as it is reported
#'   clinically. Must be shorter than tau.
#' @return Data frame with columns: var, time, value, method
make_dosing_events <- function(dose_mg, tau = 12, n_doses = 14,
                               bioavailability = 0.50, lag_h = 0) {
  if (lag_h < 0 || lag_h >= tau) {
    stop("lag_h must be in [0, tau); got ", lag_h)
  }
  dose_times <- seq(0, by = tau, length.out = n_doses) + lag_h
  data.frame(
    var    = rep("A_gut", n_doses),
    time   = dose_times,
    value  = rep(dose_mg * bioavailability, n_doses),
    method = rep("add", n_doses),
    stringsAsFactors = FALSE
  )
}

#' Convert apixaban concentration from ng/mL to nM
#' @param conc_ng_mL Concentration in ng/mL
#' @param MW Molecular weight of apixaban (459.5 g/mol)
#' @return Concentration in nM
ng_ml_to_nM <- function(conc_ng_mL, MW = 459.5) {
  conc_ng_mL * 1000 / MW
}

#' Convert apixaban concentration from nM to ng/mL
#' @param conc_nM Concentration in nM
#' @param MW Molecular weight of apixaban (459.5 g/mol)
#' @return Concentration in ng/mL
nM_to_ng_ml <- function(conc_nM, MW = 459.5) {
  conc_nM * MW / 1000
}
