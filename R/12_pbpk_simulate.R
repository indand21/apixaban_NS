# 12_pbpk_simulate.R
# PBPK simulation wrapper

#' Run PBPK simulation for apixaban
#' @param dose_mg Dose in mg (default 5)
#' @param tau Dosing interval in hours (default 12 for BID)
#' @param n_doses Number of doses (default 14, i.e., 7 days BID)
#' @param egfr eGFR in mL/min/1.73m2
#' @param body_weight Body weight in kg
#' @param dt Output time step in hours (default 0.1)
#' @param params Optional pre-computed parameter list
#' @return Data frame with time and all state/auxiliary variables
simulate_pbpk <- function(dose_mg = 5,
                          tau = 12,
                          n_doses = 14,
                          egfr = 120,
                          body_weight = 70,
                          dt = 0.1,
                          params = NULL) {

  # Get parameters
  if (is.null(params)) {
    params <- get_pbpk_params(egfr = egfr, body_weight = body_weight)
  }

  # Initial state: all compartments empty
  state <- c(
    A_gut        = 0,
    A_liver      = 0,
    A_kidney     = 0,
    A_plasma     = 0,
    A_peripheral = 0
  )

  # Dosing events: add dose * F_oral to gut at each dose time
  events <- make_dosing_events(
    dose_mg = dose_mg,
    tau = tau,
    n_doses = n_doses,
    bioavailability = params$F_oral,
    lag_h = if (is.null(params$Tlag)) 0 else params$Tlag
  )

  # Time vector
  t_end <- (n_doses - 1) * tau + tau  # One interval after last dose
  times <- seq(0, t_end, by = dt)

  # Solve ODE system
  out <- deSolve::ode(
    y      = state,
    times  = times,
    func   = pbpk_odes,
    parms  = params,
    method = "lsoda",
    events = list(data = events),
    atol   = 1e-8,
    rtol   = 1e-6
  )

  # Convert to data frame
  result <- as.data.frame(out)

  # Add metadata
  attr(result, "dose_mg")  <- dose_mg

  attr(result, "tau")      <- tau
  attr(result, "n_doses")  <- n_doses
  attr(result, "egfr")     <- egfr
  attr(result, "params")   <- params

  result
}

#' Extract steady-state interval from PBPK result
#' @param pbpk_result Data frame from simulate_pbpk
#' @param tau Dosing interval in hours
#' @return Data frame for last complete dosing interval
extract_steady_state <- function(pbpk_result, tau = 12) {
  t_max <- max(pbpk_result$time)
  t_start <- t_max - tau
  ss <- pbpk_result[pbpk_result$time >= t_start & pbpk_result$time <= t_max, ]
  ss$time_in_interval <- ss$time - t_start
  ss
}
