# 25_platelet_simulate.R
# Platelet model simulation wrapper
# Takes IIa(t) time-course from Hockin-Mann as forcing function

#' Create thrombin forcing function from coagulation result
#' @param coag_result Data frame from simulate_coagulation (must have time, IIa cols)
#' @return Interpolation function: t_seconds -> IIa (nM)
create_thrombin_forcing <- function(coag_result) {
  approxfun(coag_result$time, coag_result$IIa, rule = 2)
}

#' Run platelet model simulation
#' @param coag_result Data frame from simulate_coagulation (provides IIa time-course)
#' @param ckd_stage CKD stage for platelet dysfunction modifiers
#' @param plt_params Platelet parameters (if NULL, loads defaults and applies CKD mods)
#' @param t_end End time in seconds (default: match coag simulation)
#' @param dt Time step in seconds (default 0.5)
#' @param data_dir Path to data directory
#' @return Data frame with time and 6 platelet species
simulate_platelets <- function(coag_result,
                               ckd_stage = "Normal",
                               plt_params = NULL,
                               t_end = NULL,
                               dt = 0.5,
                               data_dir = "data/", age = NULL) {

  # Default: match coag simulation time span

  if (is.null(t_end)) {
    t_end <- max(coag_result$time)
  }

  # Get and modify platelet parameters for CKD
  if (is.null(plt_params)) {
    plt_params <- get_platelet_params(data_dir)
    plt_params <- apply_ckd_platelet_modifiers(plt_params, ckd_stage, data_dir)
    plt_params <- apply_age_platelet_modifiers(plt_params, age)
  }

  # Create thrombin forcing function from coag output
  IIa_forcing <- create_thrombin_forcing(coag_result)

  # Combine parameters with forcing
  ode_params <- c(plt_params, list(IIa_forcing = IIa_forcing))

  # Initial conditions
  ic <- get_platelet_initial_conditions(plt_params)

  # Time vector
  times <- seq(0, t_end, by = dt)

  # Solve ODE
  out <- deSolve::ode(
    y      = ic,
    times  = times,
    func   = platelet_odes,
    parms  = ode_params,
    method = "lsoda",
    atol   = 1e-8,
    rtol   = 1e-6,
    maxsteps = 50000
  )

  result <- as.data.frame(out)

  # Store metadata
  attr(result, "ckd_stage") <- ckd_stage

  result
}
