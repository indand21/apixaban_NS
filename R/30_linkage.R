# 30_linkage.R
# PBPK -> QSP bridge: unit conversion and drug forcing functions

#' Create time-varying drug forcing function from PBPK result
#' @param pbpk_result Data frame from simulate_pbpk
#' @param MW Molecular weight of apixaban (g/mol)
#' @return Function of time (in seconds) returning Cp_free in nM
create_drug_forcing <- function(pbpk_result, MW = 459.5) {
  # Extract free plasma concentration in ng/mL
  ss <- extract_steady_state(pbpk_result, tau = attr(pbpk_result, "tau"))
  time_h   <- ss$time_in_interval
  Cp_free  <- ss$C_free_ng_mL

  # Convert ng/mL to nM: C_nM = C_ng_mL * 1000 / MW
  Cp_nM <- Cp_free * 1000 / MW

  # Create interpolation function (time in hours -> nM)
  interp_fn <- approxfun(time_h, Cp_nM, rule = 2)

  # Return function of time in SECONDS (coag model uses seconds)
  function(t_sec) {
    t_h <- (t_sec / 3600) %% max(time_h)
    interp_fn(t_h)
  }
}

#' Create constant drug forcing (for quasi-steady-state TGA)
#' @param Cp_free_nM Scalar free drug concentration in nM
#' @return Scalar value (not a function, for efficiency)
create_constant_forcing <- function(Cp_free_nM) {
  Cp_free_nM
}

#' Convert apixaban Ctrough from PBPK to constant nM forcing
#' @param pbpk_result Data frame from simulate_pbpk
#' @param tau Dosing interval (hours)
#' @param MW Molecular weight
#' @return Cp_free at trough in nM
extract_trough_nM <- function(pbpk_result, tau = 12, MW = 459.5) {
  ss <- extract_steady_state(pbpk_result, tau = tau)
  ctrough_free_ng_mL <- ss$C_free_ng_mL[nrow(ss)]
  ctrough_free_ng_mL * 1000 / MW
}

#' Convert apixaban Cmax from PBPK to nM
#' @param pbpk_result Data frame from simulate_pbpk
#' @param tau Dosing interval (hours)
#' @param MW Molecular weight
#' @return Cp_free at Cmax in nM
extract_cmax_nM <- function(pbpk_result, tau = 12, MW = 459.5) {
  ss <- extract_steady_state(pbpk_result, tau = tau)
  cmax_free_ng_mL <- max(ss$C_free_ng_mL)
  cmax_free_ng_mL * 1000 / MW
}

# --- Baseline cache for therapeutic index computation ---
# Caches no-drug baseline results per CKD stage to avoid redundant computation
.baseline_cache <- new.env(parent = emptyenv())

#' Clear the baseline cache (e.g., between test runs)
clear_baseline_cache <- function() {
  rm(list = ls(.baseline_cache), envir = .baseline_cache)
}

#' Get or compute no-drug baseline for a CKD stage
#' @param ckd_stage CKD stage label
#' @param egfr eGFR value
#' @param TF_pM TF trigger in pM
#' @param t_coag_end Coagulation simulation end time (seconds)
#' @param body_weight Body weight in kg (default 70)
#' @param age Age in years (default NULL)
#' @return List with tga_metrics and hemostatic_metrics for no-drug baseline
get_nodrug_baseline <- function(ckd_stage, egfr, TF_pM = 5, t_coag_end = 1200,
                                body_weight = 70, age = NULL) {
  cache_key <- paste0(ckd_stage, "_", TF_pM, "_end", t_coag_end, "_bw", body_weight,
                      "_age", ifelse(is.null(age), "NA", age))

  if (exists(cache_key, envir = .baseline_cache)) {
    return(get(cache_key, envir = .baseline_cache))
  }

  # Run coagulation without drug
  ic <- get_coag_initial_conditions(TF_pM = TF_pM)
  ic <- apply_ckd_modifiers(ic, ckd_stage = ckd_stage)
  ic <- apply_age_coag_modifiers(ic, age = age)

  coag_nodrug <- simulate_coagulation(
    ic = ic, Cp_free_nM = 0, TF_pM = TF_pM, t_end = t_coag_end
  )
  tga_nodrug <- compute_tga_metrics(coag_nodrug)

  # Run platelet model on no-drug thrombin
  plt_nodrug <- simulate_platelets(coag_nodrug, ckd_stage = ckd_stage, age = age)
  hemo_nodrug <- compute_hemostatic_metrics(plt_nodrug)

  baseline <- list(
    tga_metrics        = tga_nodrug,
    hemostatic_metrics = hemo_nodrug
  )

  assign(cache_key, baseline, envir = .baseline_cache)
  baseline
}

# --- Coupled baseline cache ---
.coupled_baseline_cache <- new.env(parent = emptyenv())

#' Clear the coupled baseline cache
clear_coupled_baseline_cache <- function() {
  rm(list = ls(.coupled_baseline_cache), envir = .coupled_baseline_cache)
}

#' Get or compute no-drug coupled baseline for a CKD stage
#' @param ckd_stage CKD stage label
#' @param egfr eGFR value
#' @param TF_pM TF trigger in pM
#' @param t_coag_end Simulation end time (seconds)
#' @param body_weight Body weight in kg (default 70)
#' @param age Age in years (default NULL)
#' @return List with tga_metrics and hemostatic_metrics for coupled no-drug baseline
get_coupled_nodrug_baseline <- function(ckd_stage, egfr, TF_pM = 5,
                                         t_coag_end = 1200,
                                         body_weight = 70, age = NULL) {
  cache_key <- paste0(ckd_stage, "_", TF_pM, "_end", t_coag_end, "_bw", body_weight,
                      "_age", ifelse(is.null(age), "NA", age))

  if (exists(cache_key, envir = .coupled_baseline_cache)) {
    return(get(cache_key, envir = .coupled_baseline_cache))
  }

  # Coag ICs with CKD + age modifiers
  ic_coag <- get_coag_initial_conditions(TF_pM = TF_pM)
  ic_coag <- apply_ckd_modifiers(ic_coag, ckd_stage = ckd_stage)
  ic_coag <- apply_age_coag_modifiers(ic_coag, age = age)

  # Platelet ICs with CKD + age modifiers
  plt_params <- get_platelet_params()
  plt_params <- apply_ckd_platelet_modifiers(plt_params, ckd_stage)
  plt_params <- apply_age_platelet_modifiers(plt_params, age = age)
  ic_plt <- get_platelet_initial_conditions(plt_params)

  # Coupling and coag rate constants
  coupling_params <- get_coupling_params()
  rate_constants  <- get_coag_rate_constants()

  # Run coupled simulation with no drug
  coupled_out <- simulate_coupled(
    ic_coag         = ic_coag,
    ic_plt          = ic_plt,
    rate_constants  = rate_constants,
    plt_params      = plt_params,
    coupling_params = coupling_params,
    Cp_free_nM      = 0,
    t_end           = t_coag_end
  )

  tga_nodrug  <- compute_tga_metrics(coupled_out$coag_result)
  hemo_nodrug <- compute_hemostatic_metrics(coupled_out$plt_result)

  baseline <- list(
    tga_metrics        = tga_nodrug,
    hemostatic_metrics = hemo_nodrug
  )

  assign(cache_key, baseline, envir = .coupled_baseline_cache)
  baseline
}

#' Run linked PBPK-QSP simulation for a single scenario
#' @param dose_mg Apixaban dose in mg
#' @param egfr eGFR value
#' @param ckd_stage CKD stage label
#' @param use_trough If TRUE, use Ctrough; if FALSE, use time-varying forcing
#' @param TF_pM TF trigger in pM
#' @param tau Dosing interval (hours)
#' @param n_doses Number of PBPK doses for steady state
#' @param t_coag_end Coagulation simulation end time (seconds)
#' @param compute_ti If TRUE, compute therapeutic index (requires no-drug baseline)
#' @return List with pbpk_result, coag_result, pk_metrics, tga_metrics,
#'         plt_result, hemostatic_metrics, therapeutic_index (if compute_ti=TRUE)
run_linked_simulation <- function(dose_mg = 5,
                                  egfr = 120,
                                  ckd_stage = "Normal",
                                  use_trough = TRUE,
                                  TF_pM = 5,
                                  tau = 12,
                                  n_doses = 14,
                                  t_coag_end = 1200,
                                  compute_ti = TRUE,
                                  body_weight = 70,
                                  age = NULL) {

  # Step 1: Run PBPK (with CKD-dependent PK modifiers + BW/age scaling)
  params <- get_pbpk_params(egfr = egfr, body_weight = body_weight,
                            ckd_stage = ckd_stage, age = age)
  pbpk_result <- simulate_pbpk(
    dose_mg = dose_mg, tau = tau, n_doses = n_doses,
    egfr = egfr, params = params
  )

  # Step 2: Extract PK metrics
  pk_metrics <- compute_pk_metrics(pbpk_result, tau = tau)

  # Step 3: Create drug forcing for QSP
  if (use_trough) {
    Cp_free_nM <- extract_trough_nM(pbpk_result, tau = tau)
  } else {
    Cp_free_nM <- create_drug_forcing(pbpk_result)
  }

  # Step 4: Get CKD- and age-modified initial conditions
  ic <- get_coag_initial_conditions(TF_pM = TF_pM)
  ic <- apply_ckd_modifiers(ic, ckd_stage = ckd_stage)
  ic <- apply_age_coag_modifiers(ic, age = age)

  # Step 5: Run coagulation QSP
  coag_result <- simulate_coagulation(
    ic = ic,
    Cp_free_nM = Cp_free_nM,
    TF_pM = TF_pM,
    t_end = t_coag_end
  )

  # Step 6: Compute TGA metrics
  tga_metrics <- compute_tga_metrics(coag_result)

  # Step 7: Run platelet model with IIa(t) from coagulation
  plt_result <- simulate_platelets(coag_result, ckd_stage = ckd_stage, age = age)

  # Step 8: Compute hemostatic capacity metrics
  hemostatic_metrics <- compute_hemostatic_metrics(plt_result)

  # Step 9: Compute therapeutic index (vs no-drug baseline)
  # Safety score uses healthy (Normal) no-drug HC as reference to capture
  # both CKD platelet dysfunction AND drug-induced HC reduction
  therapeutic_index <- NULL
  if (compute_ti) {
    baseline <- get_nodrug_baseline(ckd_stage, egfr, TF_pM, t_coag_end,
                                    body_weight = body_weight, age = age)
    # Get healthy reference for absolute safety scoring
    healthy_baseline <- get_nodrug_baseline("Normal", 120, TF_pM, t_coag_end,
                                            body_weight = body_weight, age = age)
    therapeutic_index <- compute_therapeutic_index(
      tga_drug            = tga_metrics,
      tga_nodrug          = baseline$tga_metrics,
      hemo_drug           = hemostatic_metrics,
      hemo_nodrug         = baseline$hemostatic_metrics,
      hemo_healthy_nodrug = healthy_baseline$hemostatic_metrics
    )
  }

  list(
    pbpk_result        = pbpk_result,
    coag_result        = coag_result,
    plt_result         = plt_result,
    pk_metrics         = pk_metrics,
    tga_metrics        = tga_metrics,
    hemostatic_metrics = hemostatic_metrics,
    therapeutic_index  = therapeutic_index,
    Cp_free_nM         = if (is.numeric(Cp_free_nM)) Cp_free_nM else NA,
    dose_mg            = dose_mg,
    egfr               = egfr,
    ckd_stage          = ckd_stage,
    body_weight        = body_weight,
    age                = age
  )
}

#' Run coupled (two-way) linked PBPK-QSP simulation
#'
#' Same interface as run_linked_simulation but uses the merged 40-species ODE
#' with PS_exp -> prothrombinase feedback. Coagulation and platelet models are
#' solved simultaneously in one deSolve call.
#'
#' @inheritParams run_linked_simulation
#' @return List with same structure as run_linked_simulation, plus coupled=TRUE flag
run_coupled_linked_simulation <- function(dose_mg = 5,
                                           egfr = 120,
                                           ckd_stage = "Normal",
                                           use_trough = TRUE,
                                           TF_pM = 5,
                                           tau = 12,
                                           n_doses = 14,
                                           t_coag_end = 1200,
                                           compute_ti = TRUE,
                                           body_weight = 70,
                                           age = NULL) {

  # Step 1: Run PBPK (with BW/age scaling)
  params <- get_pbpk_params(egfr = egfr, body_weight = body_weight,
                            ckd_stage = ckd_stage, age = age)
  pbpk_result <- simulate_pbpk(
    dose_mg = dose_mg, tau = tau, n_doses = n_doses,
    egfr = egfr, params = params
  )

  # Step 2: PK metrics
  pk_metrics <- compute_pk_metrics(pbpk_result, tau = tau)

  # Step 3: Drug forcing
  if (use_trough) {
    Cp_free_nM <- extract_trough_nM(pbpk_result, tau = tau)
  } else {
    Cp_free_nM <- create_drug_forcing(pbpk_result)
  }

  # Step 4: CKD- and age-modified coag ICs
  ic_coag <- get_coag_initial_conditions(TF_pM = TF_pM)
  ic_coag <- apply_ckd_modifiers(ic_coag, ckd_stage = ckd_stage)
  ic_coag <- apply_age_coag_modifiers(ic_coag, age = age)

  # Step 5: CKD- and age-modified platelet params and ICs
  plt_params <- get_platelet_params()
  plt_params <- apply_ckd_platelet_modifiers(plt_params, ckd_stage)
  plt_params <- apply_age_platelet_modifiers(plt_params, age = age)
  ic_plt <- get_platelet_initial_conditions(plt_params)

  # Step 6: Get coupling and rate constant parameters
  coupling_params <- get_coupling_params()
  rate_constants  <- get_coag_rate_constants()

  # Step 7: Run coupled 40-species ODE
  coupled_out <- simulate_coupled(
    ic_coag         = ic_coag,
    ic_plt          = ic_plt,
    rate_constants  = rate_constants,
    plt_params      = plt_params,
    coupling_params = coupling_params,
    Cp_free_nM      = Cp_free_nM,
    t_end           = t_coag_end
  )

  coag_result <- coupled_out$coag_result
  plt_result  <- coupled_out$plt_result

  # Step 8: Compute TGA metrics
  tga_metrics <- compute_tga_metrics(coag_result)

  # Step 9: Compute hemostatic capacity
  hemostatic_metrics <- compute_hemostatic_metrics(plt_result)

  # Step 10: Therapeutic index (coupled baselines)
  therapeutic_index <- NULL
  if (compute_ti) {
    baseline <- get_coupled_nodrug_baseline(ckd_stage, egfr, TF_pM, t_coag_end,
                                            body_weight = body_weight, age = age)
    healthy_baseline <- get_coupled_nodrug_baseline("Normal", 120, TF_pM, t_coag_end,
                                                    body_weight = body_weight, age = age)
    therapeutic_index <- compute_therapeutic_index(
      tga_drug            = tga_metrics,
      tga_nodrug          = baseline$tga_metrics,
      hemo_drug           = hemostatic_metrics,
      hemo_nodrug         = baseline$hemostatic_metrics,
      hemo_healthy_nodrug = healthy_baseline$hemostatic_metrics
    )
  }

  list(
    pbpk_result        = pbpk_result,
    coag_result        = coag_result,
    plt_result         = plt_result,
    pk_metrics         = pk_metrics,
    tga_metrics        = tga_metrics,
    hemostatic_metrics = hemostatic_metrics,
    therapeutic_index  = therapeutic_index,
    Cp_free_nM         = if (is.numeric(Cp_free_nM)) Cp_free_nM else NA,
    dose_mg            = dose_mg,
    egfr               = egfr,
    ckd_stage          = ckd_stage,
    body_weight        = body_weight,
    age                = age,
    coupled            = TRUE
  )
}
