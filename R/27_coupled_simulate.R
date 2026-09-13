# 27_coupled_simulate.R
# Simulation wrapper for the coupled 40-species ODE system

# Species name vectors for splitting results
.coag_species_names <- c(
  "TF", "VII", "VIIa", "TF_VII", "TF_VIIa",
  "X", "Xa", "TF_VIIa_X", "TF_VIIa_Xa",
  "IX", "IXa", "TF_VIIa_IX",
  "II", "IIa",
  "VIII", "VIIIa", "IXa_VIIIa", "IXa_VIIIa_X",
  "V", "Va", "Va_Xa", "Va_Xa_II",
  "mIIa",
  "TFPI", "TFPI_Xa", "TFPI_Xa_TF_VIIa",
  "ATIII", "ATIII_IIa", "ATIII_Xa", "ATIII_IXa",
  "PC", "APC",
  "VIIIa_i", "Va_i"
)

.plt_species_names <- c(
  "PLT_rest", "PLT_act", "PLT_agg", "ADP", "TxA2", "PS_exp"
)

#' Run coupled coagulation-platelet simulation
#'
#' Merges 34 coagulation + 6 platelet species into a single 40-species ODE.
#' PS_exp feeds back to enhance prothrombinase (R13, R15) via Michaelis-Menten.
#' IIa drives platelet activation directly within the same timestep.
#'
#' @param ic_coag Named vector of 34 coagulation initial conditions (nM)
#' @param ic_plt Named vector of 6 platelet initial conditions (nM)
#' @param rate_constants Named list of coagulation rate constants
#' @param plt_params Named list of platelet parameters
#' @param coupling_params Named list: Km_PS, n_ps, alpha_assembly, alpha_catalysis
#' @param Cp_free_nM Free apixaban concentration (nM, scalar or function of t)
#' @param Ki_apixaban Ki for FXa inhibition (nM, default 0.08)
#' @param t_end End time in seconds (default 1200)
#' @param dt Time step in seconds (default 0.5)
#' @return List with coag_result, plt_result, full_result data frames
simulate_coupled <- function(ic_coag,
                             ic_plt,
                             rate_constants,
                             plt_params,
                             coupling_params,
                             Cp_free_nM = 0,
                             Ki_apixaban = 0.25,
                             t_end = 1200,
                             dt = 0.5) {

  # Combine initial conditions: 40-element named vector
  ic_combined <- c(ic_coag, ic_plt)
  stopifnot(length(ic_combined) == 40)

  # Combine all parameters into one list
  all_params <- c(
    rate_constants,
    plt_params,
    coupling_params,
    list(
      Cp_free_nM  = Cp_free_nM,
      Ki_apixaban = Ki_apixaban
    )
  )

  # Per-species tolerances: tight for coag (1e-10), relaxed for platelet (1e-8)
  atol_vec <- rep(1e-14, 40)

  # Time vector
  times <- seq(0, t_end, by = dt)

  # Solve combined ODE
  out <- deSolve::ode(
    y        = ic_combined,
    times    = times,
    func     = coupled_odes,
    parms    = all_params,
    method   = "lsoda",
    atol     = atol_vec,
    rtol     = 1e-12,
    maxsteps = 100000
  )

  result <- as.data.frame(out)

  # Split into coag and platelet data frames
  coag_result <- result[, c("time", .coag_species_names)]
  plt_result  <- result[, c("time", .plt_species_names)]

  # Copy auxiliary columns if present
  aux_cols <- intersect(
    c("total_thrombin", "inhib_factor", "drug_nM",
      "ps_enhance_asm", "ps_enhance_cat", "IIa_input"),
    names(result)
  )
  if (length(aux_cols) > 0) {
    coag_result <- cbind(coag_result, result[, aux_cols, drop = FALSE])
    plt_result  <- cbind(plt_result,  result[, aux_cols, drop = FALSE])
  }

  # Store metadata
  attr(coag_result, "Cp_free_nM")  <- Cp_free_nM
  attr(coag_result, "Ki_apixaban") <- Ki_apixaban
  attr(coag_result, "coupled")     <- TRUE
  attr(plt_result, "coupled")      <- TRUE

  list(
    coag_result = coag_result,
    plt_result  = plt_result,
    full_result = result
  )
}
