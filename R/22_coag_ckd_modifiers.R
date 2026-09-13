# 22_coag_ckd_modifiers.R
# Disease-stage specific modifications to coagulation factor levels.
# In this nephrotic syndrome sibling project the historical function names are
# retained for compatibility, but the preferred input table is ns_factor_levels.csv.

#' Load the calibrated nephrotic-syndrome severity scale
#'
#' The per-species multipliers in data/ns_factor_levels.csv are provisional
#' literature-directed guesses, not fitted values. Applied at face value they
#' produced an NS_Severe:Normal ETP ratio of 2.93 against an observed 1.20
#' (Kelddal 2025: 1096 vs 910 nM.min), a 248 percent overprediction of the
#' severe-NS endogenous thrombin potential.
#'
#' This scalar shrinks every species multiplier towards 1 by a common factor:
#'     modifier' = 1 + s x (modifier - 1)
#' so s = 1 reproduces the table unchanged and s = 0 removes the NS coagulation
#' gradient entirely. It is calibrated by scripts/calibrate_ns_severity.R to the
#' single observed ETP ratio.
#'
#' IMPORTANT: one aggregate ratio identifies only this one scalar. It cannot
#' distinguish an over-aggressive factor table from a cascade whose ETP is
#' over-sensitive to that table, and it does not license the relative weighting
#' between species (ATIII, FVIII and the rest keep their assumed proportions).
#'
#' The option "ns_severity_scale", when set, overrides the table. That exists so
#' scripts/calibrate_ns_severity.R can sweep the scalar without rewriting its own
#' input file mid-search; production runs leave it unset.
#'
#' @param path Path to the calibration table
#' @return Numeric severity scale
get_ns_severity_scale <- function(path = file.path("data", "ns_severity_calibration.csv")) {
  o <- getOption("ns_severity_scale")
  if (!is.null(o)) {
    if (length(o) != 1 || !is.finite(o)) {
      stop("option 'ns_severity_scale' must be a single finite number")
    }
    return(o)
  }
  if (!file.exists(path)) {
    stop("NS severity calibration table not found: ", path)
  }
  d <- read.csv(path, stringsAsFactors = FALSE)
  v <- d$value[d$parameter == "severity_scale"]
  if (length(v) != 1 || !is.finite(v)) {
    stop("ns_severity_calibration.csv must contain exactly one finite severity_scale row")
  }
  v
}

#' Get CKD modifier table
#' @param data_dir Path to data directory
#' @return Data frame with multiplicative factors by CKD stage
get_ckd_modifier_table <- function(data_dir = "data/") {
  ns_path <- file.path(data_dir, "ns_factor_levels.csv")
  if (file.exists(ns_path)) {
    return(read.csv(ns_path, stringsAsFactors = FALSE))
  }

  stop("Missing NS coagulation modifier table; no legacy CKD fallback is permitted.")
}

#' Get disease modifiers for a specific stage
#' @param ckd_stage Character stage label. In this NS project valid stages are
#'   "Normal", "NS_Mild", "NS_Moderate", "NS_Severe".
#' @param data_dir Path to data directory
#' @param severity_scale Common shrinkage applied to every multiplier, as
#'   1 + s x (modifier - 1). Defaults to the calibrated value; pass 1 to read
#'   the raw table.
#' @return Named vector of multiplicative factors
get_ckd_modifiers <- function(ckd_stage = "Normal", data_dir = "data/",
                              severity_scale = get_ns_severity_scale()) {
  mod_table <- get_ckd_modifier_table(data_dir)

  if (!ckd_stage %in% names(mod_table)) {
    valid_stages <- setdiff(names(mod_table), c("species", "source"))
    stop("Unknown disease stage: ", ckd_stage,
         ". Must be one of: ", paste(valid_stages, collapse = ", "))
  }

  mods <- mod_table[[ckd_stage]]
  names(mods) <- mod_table$species
  # Shrink every multiplier towards 1 by the calibrated severity scale
  1 + severity_scale * (mods - 1)
}

#' Apply disease modifiers to coagulation initial conditions
#' @param ic Named vector of initial conditions (from get_coag_initial_conditions)
#' @param ckd_stage Character CKD stage
#' @param data_dir Path to data directory
#' @param exclude_species Character vector of species names to exclude from modification
#'   (e.g., "TF" to keep TF at standard trigger level for clinical CAT comparison)
#' @param severity_scale Common shrinkage of the NS factor multipliers; see
#'   get_ns_severity_scale()
#' @return Modified named vector of initial conditions
apply_ckd_modifiers <- function(ic, ckd_stage = "Normal", data_dir = "data/",
                                exclude_species = "TF",
                                severity_scale = get_ns_severity_scale()) {
  if (ckd_stage == "Normal") return(ic)

  mods <- get_ckd_modifiers(ckd_stage, data_dir, severity_scale = severity_scale)

  # Remove excluded species from modifiers
  if (!is.null(exclude_species)) {
    mods <- mods[!names(mods) %in% exclude_species]
  }

  # Map modifier species names to IC species names
  species_map <- c(
    TF   = "TF",
    VII  = "VII",
    VIII = "VIII",
    V    = "V",
    II   = "II",
    IX   = "IX",
    X    = "X",
    ATIII = "ATIII",
    PC   = "PC",
    TFPI = "TFPI"
  )

  for (mod_name in names(mods)) {
    ic_name <- species_map[mod_name]
    if (!is.na(ic_name) && ic_name %in% names(ic)) {
      ic[ic_name] <- ic[ic_name] * mods[mod_name]
    }
  }

  ic
}

#' Apply age-dependent coagulation factor modifiers
#' Multiplicative adjustments applied after CKD modifiers.
#' Literature: Favaloro 2014, Mari 2008.
#' @param ic Named vector of initial conditions
#' @param age Age in years (NULL = no modification)
#' @return Modified named vector of initial conditions
apply_age_coag_modifiers <- function(ic, age = NULL) {
  if (is.null(age)) return(ic)

  # Factor VIII: +5% per decade over 40
  if ("VIII" %in% names(ic)) {
    viii_factor <- 1 + 0.05 * max(0, age - 40) / 10
    ic["VIII"] <- ic["VIII"] * viii_factor
  }

  # Protein C: -3% per decade over 50, floor at 0.7x
  if ("PC" %in% names(ic)) {
    pc_factor <- max(0.7, 1 - 0.03 * max(0, age - 50) / 10)
    ic["PC"] <- ic["PC"] * pc_factor
  }

  # Factor II (prothrombin): +3% per decade over 40
  if ("II" %in% names(ic)) {
    ii_factor <- 1 + 0.03 * max(0, age - 40) / 10
    ic["II"] <- ic["II"] * ii_factor
  }

  ic
}
