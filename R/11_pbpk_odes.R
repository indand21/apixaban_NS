# 11_pbpk_odes.R
# 5-compartment PBPK ODE system for apixaban
# Compartments: gut, liver, kidney, plasma, peripheral
#
# A deep peripheral compartment was added and then removed during the September
# 2026 Tier-1 revision. It was intended to fix the within-interval decay, which
# was too steep (peak-to-trough 4.2 against an observed 2.6, trough
# underpredicted by 42 percent). Calibration drove its flow to 0.031 L/h and its
# partition coefficient to 0.039, contributing 0.39 L of a 28.5 L steady-state
# distribution volume; removing it changed peak by 0.05 percent and trough by
# 0.20 percent. It also split the terminal phase into two near-degenerate modes
# (7.95 and 9.37 h), making the apparent terminal half-life depend on the
# regression window (8.78 h over 48-72 h, 9.34 h over 250-400 h) where the model
# without it gives a clean 8.51 h. It was therefore dropped as unsupported
# structure that introduced an artifact. The within-interval decay is instead
# explained by the calibrated absorption lag and slower absorption rate.
# State variables are AMOUNTS (mg) in each compartment

#' PBPK ODE function for deSolve
#' @param t Time (hours)
#' @param state Named vector of state variables (amounts in mg)
#' @param params Named list of PBPK parameters
#' @return List of derivatives and auxiliary variables
pbpk_odes <- function(t, state, params) {
  with(as.list(c(state, params)), {

    # --- Concentrations (mg/L) from amounts ---
    C_plasma     <- A_plasma / V_plasma
    C_liver      <- A_liver / V_liver
    C_kidney     <- A_kidney / V_kidney
    C_peripheral <- A_peripheral / V_peripheral

    # --- Free (venous) concentrations leaving tissues ---
    # Using well-stirred model: C_venous = C_tissue / Kp
    C_liver_venous      <- C_liver / Kp_liver
    C_kidney_venous     <- C_kidney / Kp_kidney
    C_peripheral_venous <- C_peripheral / Kp_peripheral

    # --- ODEs (amounts, mg) ---

    # Gut: first-order absorption
    dA_gut <- -ka * A_gut

    # Liver: portal absorption + arterial inflow - venous outflow - elimination
    # Apixaban is a low-extraction drug (ER ~0.04), so hepatic elimination is
    # RESTRICTED to the unbound fraction: rate = CLint_hepatic x fu x C_venous.
    # Applying CL_hepatic to total concentration instead would decouple
    # elimination from binding and make unbound exposure scale linearly with fu,
    # which is not physiological (AUC_unbound = F x Dose / CLint, independent of fu).
    dA_liver <- (ka * A_gut) +                  # absorption from gut (portal)
                Q_liver * C_plasma -             # arterial inflow
                Q_liver * C_liver_venous -       # venous outflow
                CLint_hepatic * fu_plasma * C_liver_venous   # hepatic elimination

    CL_proteinuria_eff <- if (exists("CL_proteinuria")) CL_proteinuria else 0

    # Kidney: arterial inflow - venous outflow - renal elimination.
    # Glomerular filtration and tubular secretion also act on unbound drug only.
    # CL_proteinuria is the exception: it represents urinary loss of apixaban
    # bound to albumin that crosses the damaged glomerular barrier, so it is the
    # one elimination route that removes BOUND drug and is weighted by (1 - fu).
    dA_kidney <- Q_kidney * C_plasma -           # arterial inflow
                 Q_kidney * C_kidney_venous -    # venous outflow
                 CLint_renal * fu_plasma * C_kidney_venous -  # renal elimination
                 CL_proteinuria_eff * (1 - fu_plasma) * C_kidney_venous

    # Peripheral: distribution
    dA_peripheral <- Q_peripheral * C_plasma -   # arterial inflow
                     Q_peripheral * C_peripheral_venous  # venous outflow

    # Plasma: venous return from all tissues minus arterial outflow
    dA_plasma <- Q_liver * C_liver_venous +      # venous return from liver
                 Q_kidney * C_kidney_venous +    # venous return from kidney
                 Q_peripheral * C_peripheral_venous -  # venous return from peripheral
                 (Q_liver + Q_kidney + Q_peripheral) * C_plasma  # total arterial outflow

    # --- Auxiliary outputs ---
    C_plasma_total <- C_plasma             # Total plasma concentration (mg/L = ug/mL)
    C_plasma_free  <- C_plasma * fu_plasma # Free plasma concentration (mg/L)
    # Convert to ng/mL: 1 mg/L = 1000 ng/mL
    C_total_ng_mL  <- C_plasma_total * 1000
    C_free_ng_mL   <- C_plasma_free * 1000

    list(
      c(dA_gut, dA_liver, dA_kidney, dA_plasma, dA_peripheral),
      C_plasma_total = C_plasma_total,
      C_plasma_free  = C_plasma_free,
      C_total_ng_mL  = C_total_ng_mL,
      C_free_ng_mL   = C_free_ng_mL
    )
  })
}
