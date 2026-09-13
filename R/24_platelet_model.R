# 24_platelet_model.R
# Simplified mechanistic platelet activation/aggregation model
# Six phenomenological states and eight rates; no validated bleeding endpoint
# Inherited phenomenological parameterization; not a reproduction of a published platelet model

#' Get platelet model parameters
#' @param data_dir Path to data directory
#' @return Named list of rate constants and Km values
get_platelet_params <- function(data_dir = "data/") {
  csv_path <- file.path(data_dir, "platelet_parameters.csv")
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    params <- as.list(setNames(df$value, df$parameter))
    return(params)
  }

  # Fallback hardcoded
  list(
    k_thr     = 0.015,    # 1/s, thrombin-induced activation
    Km_thr    = 100.0,    # nM, PAR half-max (full aggregation EC50)
    n_hill    = 2.0,      # Hill coefficient for cooperative PAR signaling
    k_adp     = 0.005,    # 1/s, ADP-induced activation
    Km_adp    = 1000.0,   # nM
    k_txa     = 0.004,    # 1/s, TxA2-induced activation
    Km_txa    = 200.0,    # nM
    k_agg     = 0.003,    # 1/s, aggregation rate
    Km_fbg    = 500.0,    # nM, fibrinogen Km
    FBG       = 7500.0,   # nM, plasma fibrinogen
    k_rel_adp = 100.0,    # nM/s per nM PLT_act
    k_txa_syn = 25.0,     # nM/s per nM PLT_act
    k_ps      = 0.003,    # 1/s, PS exposure
    k_adp_deg = 0.05,     # 1/s, ADP degradation
    k_txa_deg = 0.023,    # 1/s, TxA2 degradation
    k_deact   = 0.008,    # 1/s, platelet deactivation/desensitization
    PLT_rest_0 = 3000.0   # model-equivalent units, initial resting platelet pool
  )
}

#' Get disease platelet modifier table
#' @param data_dir Path to data directory
#' @return Data frame with multipliers by CKD stage
get_ckd_platelet_modifier_table <- function(data_dir = "data/") {
  ns_path <- file.path(data_dir, "ns_platelet_modifiers.csv")
  if (file.exists(ns_path)) {
    return(read.csv(ns_path, stringsAsFactors = FALSE))
  }

  stop("Missing NS platelet modifier table; no legacy CKD fallback is permitted.")
}

#' Apply disease platelet modifiers to platelet parameters
#' @param params Named list of platelet parameters (from get_platelet_params)
#' @param ckd_stage Character CKD stage
#' @param data_dir Path to data directory
#' @return Modified parameter list
apply_ckd_platelet_modifiers <- function(params, ckd_stage = "Normal",
                                          data_dir = "data/") {
  if (ckd_stage == "Normal") return(params)

  mod_table <- get_ckd_platelet_modifier_table(data_dir)

  if (!ckd_stage %in% names(mod_table)) {
    valid_stages <- setdiff(names(mod_table), c("parameter", "mechanism", "source"))
    stop("Unknown disease stage: ", ckd_stage,
         ". Must be one of: ", paste(valid_stages, collapse = ", "))
  }

  for (i in seq_len(nrow(mod_table))) {
    pname <- mod_table$parameter[i]
    if (pname %in% names(params)) {
      params[[pname]] <- params[[pname]] * mod_table[[ckd_stage]][i]
    }
  }

  params
}

#' Apply age-dependent platelet modifiers
#' Platelet count declines with age (>60 yr: -10% per decade, floor 0.5x)
#' @param params Named list of platelet parameters
#' @param age Age in years (NULL = no modification)
#' @return Modified parameter list
apply_age_platelet_modifiers <- function(params, age = NULL) {
  if (is.null(age)) return(params)

  # PLT_rest_0: -10% per decade over 60, floor at 0.5x
  plt_factor <- max(0.5, 1 - 0.10 * max(0, age - 60) / 10)
  params$PLT_rest_0 <- params$PLT_rest_0 * plt_factor

  params
}

#' Get platelet model initial conditions
#' @param params Platelet parameters (for PLT_rest_0)
#' @return Named vector of initial concentrations (nM)
get_platelet_initial_conditions <- function(params = NULL) {
  if (is.null(params)) params <- get_platelet_params()
  c(
    PLT_rest = params$PLT_rest_0,
    PLT_act  = 0.0,
    PLT_agg  = 0.0,
    ADP      = 0.0,
    TxA2     = 0.0,
    PS_exp   = 0.0
  )
}

#' Platelet ODE function for deSolve
#' @param t Time (seconds)
#' @param state Named vector of 6 species concentrations (nM)
#' @param params Named list containing rate constants and IIa_forcing function
#' @return List of derivatives
platelet_odes <- function(t, state, params) {
  with(as.list(c(state, params)), {

    # Get thrombin concentration from forcing function
    if (is.function(IIa_forcing)) {
      IIa <- IIa_forcing(t)
    } else {
      IIa <- IIa_forcing
    }

    # Clamp negative values to zero (numerical artifact protection)
    PLT_rest <- max(PLT_rest, 0)
    PLT_act  <- max(PLT_act, 0)
    ADP      <- max(ADP, 0)
    TxA2     <- max(TxA2, 0)

    # R1: PLT_rest -> PLT_act (thrombin + ADP + TxA2 driven activation)
    # Thrombin uses Hill kinetics for cooperative PAR1/PAR4 signaling
    IIa_n <- max(0, IIa)^n_hill
    Km_n  <- Km_thr^n_hill
    activation_rate <- PLT_rest * (
      k_thr * IIa_n / (IIa_n + Km_n) +
      k_adp * ADP / (ADP + Km_adp) +
      k_txa * TxA2 / (TxA2 + Km_txa)
    )

    # R2: PLT_act -> PLT_agg (GPIIb/IIIa + fibrinogen)
    # Thrombin-dependent aggregation gate is phenomenological.
    # No explicit fibrin formation or crosslinking state is represented.
    aggregation_rate <- k_agg * PLT_act * FBG / (FBG + Km_fbg) * IIa_n / (IIa_n + Km_n)

    # R3: PLT_act -> ADP release (dense granule secretion)
    adp_release_rate <- k_rel_adp * PLT_act

    # R4: PLT_act -> TxA2 synthesis (COX-1)
    txa2_syn_rate <- k_txa_syn * PLT_act

    # R5: PLT_act -> PS_exp (phosphatidylserine exposure)
    ps_rate <- k_ps * PLT_act

    # R6: ADP decay (ectonucleotidases)
    adp_decay_rate <- k_adp_deg * ADP

    # R7: TxA2 decay (short half-life ~30s)
    txa2_decay_rate <- k_txa_deg * TxA2

    # R8: PLT_act deactivation/desensitization (receptor internalization)
    deactivation_rate <- k_deact * PLT_act

    # ODEs
    dPLT_rest <- -activation_rate
    dPLT_act  <- activation_rate - aggregation_rate - deactivation_rate
    dPLT_agg  <- aggregation_rate
    dADP      <- adp_release_rate - adp_decay_rate
    dTxA2     <- txa2_syn_rate - txa2_decay_rate
    dPS_exp   <- ps_rate

    list(c(dPLT_rest, dPLT_act, dPLT_agg, dADP, dTxA2, dPS_exp),
         IIa_input = IIa)
  })
}
