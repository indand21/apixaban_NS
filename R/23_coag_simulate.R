# 23_coag_simulate.R
# Coagulation QSP simulation wrapper

#' Run Hockin-Mann coagulation simulation
#' @param ic Named vector of initial conditions (nM). If NULL, uses defaults.
#' @param rate_constants Named list of rate constants. If NULL, uses defaults.
#' @param Cp_free_nM Free apixaban concentration in nM (scalar or function of t)
#' @param Ki_apixaban Ki for FXa inhibition in nM (default 0.08)
#' @param TF_pM TF trigger concentration in pM (default 25)
#' @param t_end End time in seconds (default 1200)
#' @param dt Time step in seconds (default 0.5)
#' @return Data frame with time and all species concentrations
simulate_coagulation <- function(ic = NULL,
                                 rate_constants = NULL,
                                 Cp_free_nM = 0,
                                 Ki_apixaban = 0.25,
                                 TF_pM = 5,
                                 t_end = 1200,
                                 dt = 0.5) {

  # Default parameters
  if (is.null(rate_constants)) {
    rate_constants <- get_coag_rate_constants()
  }
  if (is.null(ic)) {
    ic <- get_coag_initial_conditions(TF_pM = TF_pM)
  }

  # Combine rate constants with drug parameters
  params <- c(
    rate_constants,
    list(
      Cp_free_nM  = Cp_free_nM,
      Ki_apixaban = Ki_apixaban
    )
  )

  # Time vector
  times <- seq(0, t_end, by = dt)

  # Solve with lsoda (stiff solver), tight tolerances
  out <- deSolve::ode(
    y      = ic,
    times  = times,
    func   = coag_odes,
    parms  = params,
    method = "lsoda",
    atol   = 1e-10,
    rtol   = 1e-8,
    maxsteps = 50000
  )

  result <- as.data.frame(out)

  # Store metadata
  attr(result, "Cp_free_nM")  <- Cp_free_nM
  attr(result, "Ki_apixaban") <- Ki_apixaban
  attr(result, "TF_pM")       <- TF_pM

  result
}

#' Run coagulation simulation at multiple drug concentrations
#' @param conc_nM_vec Vector of free drug concentrations in nM
#' @param ... Additional arguments passed to simulate_coagulation
#' @return List of simulation results, named by concentration
simulate_coag_dose_response <- function(conc_nM_vec, ...) {
  results <- lapply(conc_nM_vec, function(conc) {
    simulate_coagulation(Cp_free_nM = conc, ...)
  })
  names(results) <- paste0(round(conc_nM_vec, 2), "_nM")
  results
}
