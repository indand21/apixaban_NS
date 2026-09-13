# 41_pk_analysis.R
# PK metric extraction from PBPK simulation results

#' Compute steady-state PK metrics
#' @param pbpk_result Data frame from simulate_pbpk
#' @param tau Dosing interval in hours (default 12)
#' @return Named list of PK metrics
compute_pk_metrics <- function(pbpk_result, tau = 12) {
  ss <- extract_steady_state(pbpk_result, tau = tau)

  # Total plasma concentration (ng/mL)
  cmax_info <- compute_cmax(ss$time_in_interval, ss$C_total_ng_mL)
  ctrough   <- ss$C_total_ng_mL[nrow(ss)]
  auc_tau   <- compute_auc(ss$time_in_interval, ss$C_total_ng_mL)

  # Free plasma concentration (ng/mL)
  cmax_free_info <- compute_cmax(ss$time_in_interval, ss$C_free_ng_mL)
  ctrough_free   <- ss$C_free_ng_mL[nrow(ss)]
  auc_tau_free   <- compute_auc(ss$time_in_interval, ss$C_free_ng_mL)

  list(
    Cmax_total    = cmax_info$Cmax,
    Tmax          = cmax_info$Tmax,
    Ctrough_total = ctrough,
    AUCtau_total  = auc_tau,
    Cmax_free     = cmax_free_info$Cmax,
    Ctrough_free  = ctrough_free,
    AUCtau_free   = auc_tau_free
  )
}

#' Compute PK metrics across multiple disease stages
#' @param dose_mg Dose in mg
#' @param egfr_values Named vector of eGFR values (names = disease stage labels)
#' @param tau Dosing interval in hours
#' @param n_doses Number of doses for steady state
#' @return Data frame with PK metrics per disease stage
compute_pk_across_ckd <- function(dose_mg = 5,
                                  egfr_values = c(Normal = 100, NS_Mild = 90,
                                                  NS_Moderate = 90,
                                                  NS_Severe = 90),
                                  tau = 12, n_doses = 14,
                                  body_weight = 70, age = NULL) {

  results <- purrr::map_dfr(names(egfr_values), function(stage) {
    egfr <- egfr_values[stage]
    params <- get_pbpk_params(egfr = egfr, body_weight = body_weight,
                              ckd_stage = stage, age = age)
    sim <- simulate_pbpk(dose_mg = dose_mg, tau = tau, n_doses = n_doses,
                         egfr = egfr, params = params)
    pk <- compute_pk_metrics(sim, tau = tau)
    data.frame(
      ckd_stage     = stage,
      egfr          = egfr,
      dose_mg       = dose_mg,
      Cmax_total    = pk$Cmax_total,
      Tmax          = pk$Tmax,
      Ctrough_total = pk$Ctrough_total,
      AUCtau_total  = pk$AUCtau_total,
      Cmax_free     = pk$Cmax_free,
      Ctrough_free  = pk$Ctrough_free,
      AUCtau_free   = pk$AUCtau_free,
      stringsAsFactors = FALSE
    )
  })

  results$ckd_stage <- factor(results$ckd_stage,
                              levels = names(egfr_values))
  results
}
