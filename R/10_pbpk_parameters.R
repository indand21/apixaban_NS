# 10_pbpk_parameters.R
# 5-compartment PBPK parameters for apixaban
# Sources: Xu et al. 2021, FDA Clinical Pharmacology Review, Byon et al. 2019

#' Load disease-dependent PK modifiers from CSV
#' @param ckd_stage Disease stage label. In this NS project valid stages are
#'   "Normal", "NS_Mild", "NS_Moderate", "NS_Severe".
#' @return Named list with multipliers/terms: fu_plasma, CLint_hepatic,
#'   CL_proteinuria
#' @details The hepatic modifier scales the INTRINSIC clearance of unbound drug
#'   (`CLint_hepatic`), not the observed total-drug clearance. It was renamed from
#'   `CL_hepatic` because that ambiguity is what allowed binding to be varied
#'   without any compensating change in elimination.
get_ckd_pk_modifiers <- function(ckd_stage = "Normal") {
  ns_path <- file.path("data", "ns_pk_modifiers.csv")
  csv_path <- ns_path # No silent fallback to a different disease
  if (!file.exists(csv_path)) {
    stop("Disease PK modifiers file not found: ", csv_path)
  }
  mod_df <- read.csv(csv_path, stringsAsFactors = FALSE)

  valid_stages <- setdiff(names(mod_df), c("parameter", "source"))
  if (!(ckd_stage %in% valid_stages)) {
    stop("Unknown disease stage '", ckd_stage, "'. Valid stages: ",
         paste(valid_stages, collapse = ", "))
  }

  required <- c("fu_plasma", "CLint_hepatic", "CL_proteinuria")
  missing <- setdiff(required, mod_df$parameter)
  if (length(missing)) {
    stop("Modifier file ", csv_path, " is missing required parameter row(s): ",
         paste(missing, collapse = ", "),
         ". Note that 'CL_hepatic' was renamed to 'CLint_hepatic'.")
  }

  mods <- setNames(mod_df[[ckd_stage]], mod_df$parameter)
  as.list(mods)
}

#' Load calibrated absorption/distribution parameters
#'
#' These three quantities are NOT independently measured values. They are fitted
#' to the four aggregate steady-state moments reported by Frost et al. 2013 for
#' 5 mg BID on day 7 (Cmax, Ctrough, AUCtau, Tmax) by
#' scripts/calibrate_pk_absorption_distribution.R. Systemic clearance is held at
#' the label values (3.3 L/h total, 0.9 L/h renal) and is not fitted.
#'
#' Note on F_oral: oral data identify only the ratio CL/F, never F and CL
#' separately. Holding CL at the label systemic value, the observed AUCtau
#' implies CL/F of about 4.75 L/h and therefore F of about 0.69, which exceeds
#' the label's approximately 50 percent. That gap is a real and unresolved
#' conflict between the label's intravenous clearance and the observed oral
#' exposure; it is absorbed into F here and reported, not hidden.
#'
#' @param path Path to the calibration table
#' @return Named numeric vector of calibrated parameters
get_pk_calibration <- function(path = file.path("data", "pk_calibration.csv")) {
  if (!file.exists(path)) {
    stop("PK calibration table not found: ", path)
  }
  d <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("F_oral", "ka", "Tlag")
  missing <- setdiff(required, d$parameter)
  if (length(missing)) {
    stop("PK calibration table is missing required row(s): ",
         paste(missing, collapse = ", "))
  }
  setNames(d$value, d$parameter)
}

#' Recompute every derived quantity that depends on fu or intrinsic clearance
#'
#' Two families of derived quantities depend on the unbound fraction:
#'   - effective TOTAL-drug clearances: CL = fu x CLint (restrictive clearance),
#'     except CL_proteinuria, which removes BOUND drug and is weighted (1 - fu);
#'   - tissue:plasma partition coefficients: Kp = Kp_ref x fu / fu_normal.
#'     Tissue accumulation is driven by the UNBOUND plasma concentration, so
#'     with tissue binding unchanged Kp scales with the plasma unbound fraction.
#'
#' Together these make volume of distribution rise with fu and, for a drug that
#' is restrictively cleared AND restrictively distributed, leave the terminal
#' half-life less sensitive to binding; it need not be exactly invariant.
#'
#' Call this after perturbing `CLint_*`, `fu_plasma` or `CL_proteinuria` on a
#' parameter list, so that the quantities the ODEs use cannot drift out of step
#' with the reported ones.
#'
#' @param params Full parameter list from get_pbpk_params()
#' @return The same list with CL_renal, CL_hepatic, CL_total, Kp_liver,
#'   Kp_kidney and Kp_peripheral (re)set
recompute_derived_params <- function(params) {
  params$CL_hepatic <- params$CLint_hepatic * params$fu_plasma
  params$CL_renal   <- params$CLint_renal   * params$fu_plasma
  params$CL_total   <- params$CL_renal + params$CL_hepatic +
                       params$CL_proteinuria * (1 - params$fu_plasma)
  fu_rel <- params$fu_plasma / params$fu_normal
  params$Kp_liver      <- params$Kp_liver_ref      * fu_rel
  params$Kp_kidney     <- params$Kp_kidney_ref     * fu_rel
  params$Kp_peripheral <- params$Kp_peripheral_ref * fu_rel
  params
}

#' Get PBPK parameters for apixaban
#' @param egfr eGFR in mL/min/1.73m2 (default 120 = normal)
#' @param body_weight Body weight in kg (default 70)
#' @param ckd_stage Disease stage label for PK modifiers (default "Normal")
#' @param age Age in years (default NULL = no age scaling)
#' @return Named list of all PBPK parameters
get_pbpk_params <- function(egfr = 120, body_weight = 70, ckd_stage = "Normal",
                            age = NULL) {

  # --- Drug properties ---
  MW         <- 459.5    # g/mol
  .cal       <- get_pk_calibration()
  F_oral_label <- 0.50   # Eliquis label value, retained for reporting only
  F_oral     <- unname(.cal["F_oral"])  # Calibrated; see get_pk_calibration()
  fu_normal  <- 0.13     # Fraction unbound in NORMAL plasma (87% protein bound)
  fu_plasma  <- fu_normal
  ka         <- unname(.cal["ka"])      # Calibrated against observed Tmax
  Tlag       <- unname(.cal["Tlag"])    # Calibrated absorption lag (h)

  # --- Clearance ---
  # IV systemic reference is 3.3 L/h, not apparent oral clearance.
  # fu*CLint is a low-extraction local coefficient; finite flow makes exact
  # whole-body systemic clearance slightly smaller. No fitted IV model is claimed.
  # Bioavailability is applied only at the dose input.
  # Renal fraction ~27% of total clearance
  CL_total_ref   <- 3.3  # Systemic IV clearance; do not multiply by F
  CL_renal_ref   <- 0.9  # Systemic renal clearance
  CL_hepatic_ref <- CL_total_ref - CL_renal_ref  # Lumped nonrenal clearance

  # Apixaban is a low-extraction drug, so hepatic metabolism, glomerular
  # filtration and tubular secretion all act on the UNBOUND fraction only
  # (restrictive clearance). The ODEs must therefore be driven by intrinsic
  # clearances of unbound drug, obtained by dividing the observed total-drug
  # clearances by the normal fraction unbound:
  #     CL_total_observed = fu x CL_int   =>   CL_int = CL_total_observed / fu_normal
  # At fu_plasma = fu_normal the two parameterisations are numerically identical,
  # so the healthy-subject baseline is unchanged; only states that alter plasma
  # binding diverge. Driving the ODEs with total-drug clearance instead would
  # make unbound exposure scale linearly with fu, which is not physiological.
  CLint_hepatic_ref <- CL_hepatic_ref / fu_normal
  CLint_renal_ref   <- CL_renal_ref   / fu_normal

  # --- Disease-dependent PK modifiers ---
  ckd_mods   <- get_ckd_pk_modifiers(ckd_stage)
  fu_plasma  <- fu_plasma * ckd_mods$fu_plasma
  CLint_hepatic_ref <- CLint_hepatic_ref * ckd_mods$CLint_hepatic
  CL_proteinuria_ref <- if (!is.null(ckd_mods$CL_proteinuria)) ckd_mods$CL_proteinuria else 0

  # --- eGFR-dependent renal clearance ---
  CLint_renal <- CLint_renal_ref * (egfr / 120)
  CLint_hepatic <- CLint_hepatic_ref
  CL_proteinuria <- CL_proteinuria_ref

  # --- Allometric body weight scaling (reference = 70 kg) ---
  bw_ratio     <- body_weight / 70
  bw_allo      <- bw_ratio^0.75  # allometric exponent for CL and organ volumes

  CLint_hepatic <- CLint_hepatic * bw_allo
  CL_proteinuria <- CL_proteinuria * bw_allo

  # --- Age-dependent hepatic CL decline (>40 yr: -0.7% per year, floor 50%) ---
  if (!is.null(age)) {
    age_factor <- max(0.5, 1 - 0.007 * max(0, age - 40))
    CLint_hepatic <- CLint_hepatic * age_factor
  }

  # --- Compartment volumes (BW-scaled) ---
  V_gut        <- 1.0  * bw_allo   # L (transit compartment)
  V_liver      <- 1.5  * bw_allo   # L
  V_kidney     <- 0.5  * bw_allo   # L
  V_plasma     <- 4.5  * bw_ratio  # L (central; scales linearly with BW)
  V_peripheral <- 14.5 * bw_allo   # L (peripheral tissue volume)

  # --- Inter-compartmental clearances (BW-scaled) ---
  Q_liver      <- 90.0 * bw_allo   # L/h hepatic blood flow
  Q_kidney     <- 72.0 * bw_allo   # L/h renal blood flow
  Q_peripheral <- 4.0  * bw_allo   # L/h peripheral distribution

  # --- Tissue partition coefficients (Kp) ---
  # Reference values are the tissue:plasma ratios AT NORMAL plasma binding,
  # chosen so Vd_ss = V_plasma + sum(V_t x Kp_t) is approximately 28-30 L.
  # Tissue accumulation is driven by the UNBOUND plasma concentration, so with
  # tissue binding unchanged Kp scales with the plasma unbound fraction:
  #     Kp = Kp_ref x fu_plasma / fu_normal
  # At fu_plasma = fu_normal this is a no-op, so the healthy baseline and all
  # existing validation are unchanged. Holding Kp fixed while fu rises would
  # keep volume of distribution constant and make the terminal half-life fall
  # spuriously as clearance rises, which is not physiological for a drug that
  # is restrictively distributed as well as restrictively cleared. The derived
  # Kp values are set by recompute_derived_params().
  Kp_liver_ref      <- 2.5    # Liver-to-plasma partition at fu_normal
  Kp_kidney_ref     <- 2.0    # Kidney-to-plasma partition at fu_normal
  Kp_peripheral_ref <- 1.3    # Peripheral-to-plasma partition at fu_normal

  params <- list(
    MW             = MW,
    F_oral         = F_oral,
    F_oral_label   = F_oral_label,
    fu_normal      = fu_normal,
    fu_plasma      = fu_plasma,
    ka             = ka,
    Tlag           = Tlag,
    CLint_renal    = CLint_renal,
    CLint_hepatic  = CLint_hepatic,
    CL_proteinuria = CL_proteinuria,
    CL_renal_ref   = CL_renal_ref,
    CL_hepatic_ref = CL_hepatic_ref,
    V_gut          = V_gut,
    V_liver        = V_liver,
    V_kidney       = V_kidney,
    V_plasma       = V_plasma,
    V_peripheral   = V_peripheral,
    Q_liver        = Q_liver,
    Q_kidney       = Q_kidney,
    Q_peripheral   = Q_peripheral,
    Kp_liver_ref      = Kp_liver_ref,
    Kp_kidney_ref     = Kp_kidney_ref,
    Kp_peripheral_ref = Kp_peripheral_ref,
    egfr           = egfr,
    body_weight    = body_weight,
    age            = age
  )
  recompute_derived_params(params)
}
