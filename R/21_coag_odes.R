# 21_coag_odes.R
# Reduced Hockin-Mann-inspired cascade; not the original published model.
# 34 species, 27+ reactions, mass-action kinetics
# FXa inhibition by apixaban applied via competitive inhibition factor

#' Coagulation ODE function for deSolve
#' @param t Time (seconds)
#' @param state Named vector of 34 species concentrations (nM)
#' @param params Named list containing:
#'   - rate constants (k1-k41)
#'   - Cp_free_nM: free apixaban concentration (nM), scalar or function of t
#'   - Ki_apixaban: free FXa inhibition constant at 37 C (nM), default 0.25
#' @return List of derivatives
coag_odes <- function(t, state, params) {
  with(as.list(c(state, params)), {

    # --- FXa inhibition by apixaban ---
    # Get drug concentration (constant or time-varying)
    if (is.function(Cp_free_nM)) {
      drug_conc <- Cp_free_nM(t)
    } else {
      drug_conc <- Cp_free_nM
    }
    # Competitive inhibition factor for Xa-dependent reactions
    inhib_Xa <- 1.0 / (1.0 + drug_conc / Ki_apixaban)

    # Rapid-equilibrium mixed inhibition of total prothrombinase pools.
    # Luettgen et al. 2011, DOI 10.3109/14756366.2010.535793 (37 C).
    phi_e <- 1 / (1 + drug_conc / 0.62)
    phi_es <- 1 / (1 + drug_conc / 1.68)
    # =====================================================
    # REACTION RATES (all in nM/s)
    # =====================================================

    # --- INITIATION PHASE ---

    # R1: TF + VII <-> TF:VII
    r1_fwd <- k1 * TF * VII
    r1_rev <- k2 * TF_VII

    # R2: TF + VIIa <-> TF:VIIa
    r2_fwd <- k3 * TF * VIIa
    r2_rev <- k4 * TF_VIIa

    # R3: TF:VIIa activates TF:VII -> TF:VIIa (autoactivation)
    r3 <- k5 * TF_VIIa * TF_VII

    # R4: Xa activates TF:VII -> TF:VIIa (Xa-dependent, inhibited)
    r4 <- k6 * Xa * inhib_Xa * TF_VII

    # R5: IIa activates VII -> VIIa
    r5 <- k7 * IIa * VII

    # R6: TF:VIIa + X <-> TF:VIIa:X
    r6_fwd <- k8 * TF_VIIa * X
    r6_rev <- k9 * TF_VIIa_X

    # R7: TF:VIIa:X -> TF:VIIa + Xa (extrinsic Xa generation)
    r7 <- k10 * TF_VIIa_X

    # R8: TF:VIIa + Xa <-> TF:VIIa:Xa (product inhibition by Xa)
    r8_fwd <- k11 * TF_VIIa * Xa * inhib_Xa
    r8_rev <- k12 * TF_VIIa_Xa

    # R9: TF:VIIa + IX <-> TF:VIIa:IX
    r9_fwd <- k13 * TF_VIIa * IX
    r9_rev <- k14 * TF_VIIa_IX

    # R10: TF:VIIa:IX -> TF:VIIa + IXa
    r10 <- k15 * TF_VIIa_IX

    # --- PROPAGATION PHASE ---

    # R11: Free Xa + II (slow prothrombin activation without Va)
    # Xa-dependent, inhibited by apixaban
    r11_fwd <- k16 * Xa * inhib_Xa * II

    # R12: Xa:II -> mIIa (meizothrombin, slow pathway)
    # Modeled as: free Xa converts II to mIIa slowly
    # This rate is already captured in r11

    # R13: Va + Xa <-> Va:Xa (prothrombinase assembly)
    # Both bound and unbound Xa are represented in the total assembly pool
    r13_fwd <- k19 * Va * Xa
    r13_rev <- k20 * Va_Xa

    # R14: Va:Xa + II <-> Va:Xa:II (substrate binding)
    # Prothrombinase complex is Xa-dependent
    r14_fwd <- k21 * Va_Xa * phi_e * II
    r14_rev <- k22 * Va_Xa_II * phi_es

    # R15: Va:Xa:II -> Va:Xa + IIa (prothrombinase catalysis)
    # The catalytic step of prothrombinase
    r15 <- k23 * Va_Xa_II * phi_es

    # R16: IXa + VIIIa <-> IXa:VIIIa (tenase assembly)
    r16_fwd <- k24 * IXa * VIIIa
    r16_rev <- k25 * IXa_VIIIa

    # R17: IXa:VIIIa + X <-> IXa:VIIIa:X
    r17_fwd <- k26 * IXa_VIIIa * X
    r17_rev <- k27 * IXa_VIIIa_X

    # R18: IXa:VIIIa:X -> IXa:VIIIa + Xa (intrinsic tenase Xa generation)
    r18 <- k28 * IXa_VIIIa_X

    # --- POSITIVE FEEDBACK ---

    # R19: IIa activates VIII -> VIIIa
    r19 <- k29 * IIa * VIII

    # R20: VIIIa spontaneous inactivation
    r20 <- k30 * VIIIa

    # R21: IIa activates V -> Va
    r21 <- k31 * IIa * V

    # --- INHIBITORY PATHWAYS ---

    # R22: ATIII + IIa -> ATIII:IIa (irreversible)
    r22 <- k32 * ATIII * IIa

    # R23: ATIII + Xa -> ATIII:Xa (irreversible)
    # ATIII binds free Xa; apixaban competes for Xa active site
    r23 <- k33 * ATIII * Xa * inhib_Xa

    # R24: ATIII + IXa -> ATIII:IXa (irreversible)
    r24 <- k34 * ATIII * IXa

    # R25: TFPI + Xa <-> TFPI:Xa
    r25_fwd <- k35 * TFPI * Xa * inhib_Xa
    r25_rev <- k36 * TFPI_Xa

    # R26: TFPI:Xa + TF:VIIa <-> TFPI:Xa:TF:VIIa (quaternary complex)
    r26_fwd <- k37 * TFPI_Xa * TF_VIIa
    r26_rev <- k38 * TFPI_Xa_TF_VIIa

    # R27: optional implicit TM-dependent PC activation; k39=0 in primary assay
    r27 <- k39 * IIa * PC

    # R28: APC + Va -> APC + Va_i (APC inactivates Va)
    r28 <- k40 * APC * Va

    # R29: APC + VIIIa -> APC + VIIIa_i (APC inactivates VIIIa)
    r29 <- k41 * APC * VIIIa

    # =====================================================
    # ODEs for 34 species (nM/s)
    # =====================================================

    # 1. TF
    dTF <- -r1_fwd + r1_rev - r2_fwd + r2_rev

    # 2. VII
    dVII <- -r1_fwd + r1_rev - r5

    # 3. VIIa
    dVIIa <- -r2_fwd + r2_rev + r5

    # 4. TF:VII
    dTF_VII <- r1_fwd - r1_rev - r3 - r4

    # 5. TF:VIIa
    dTF_VIIa <- r2_fwd - r2_rev + r3 + r4 -
                r6_fwd + r6_rev + r7 -
                r8_fwd + r8_rev -
                r9_fwd + r9_rev + r10 -
                r26_fwd + r26_rev

    # 6. X (Factor X)
    dX <- -r6_fwd + r6_rev - r17_fwd + r17_rev

    # 7. Xa (free Factor Xa)
    dXa <- r7 + r18 -           # Generation from extrinsic/intrinsic tenase
           r8_fwd + r8_rev -    # Product inhibition complex
           r13_fwd + r13_rev -  # Prothrombinase assembly
           r23 -                # ATIII inhibition
           r25_fwd + r25_rev    # TFPI inhibition

    # 8. TF:VIIa:X
    dTF_VIIa_X <- r6_fwd - r6_rev - r7

    # 9. TF:VIIa:Xa
    dTF_VIIa_Xa <- r8_fwd - r8_rev

    # 10. IX
    dIX <- -r9_fwd + r9_rev

    # 11. IXa
    dIXa <- r10 - r16_fwd + r16_rev - r24

    # 12. TF:VIIa:IX
    dTF_VIIa_IX <- r9_fwd - r9_rev - r10

    # 13. II (Prothrombin)
    dII <- -r11_fwd - r14_fwd + r14_rev

    # 14. IIa (Thrombin)
    dIIa <- r15 + k18 * mIIa - r22  # From prothrombinase + mIIa conversion - ATIII

    # 15. VIII
    dVIII <- -r19

    # 16. VIIIa
    dVIIIa <- r19 - r20 - r16_fwd + r16_rev - r29

    # 17. IXa:VIIIa (tenase complex)
    dIXa_VIIIa <- r16_fwd - r16_rev - r17_fwd + r17_rev + r18

    # 18. IXa:VIIIa:X
    dIXa_VIIIa_X <- r17_fwd - r17_rev - r18

    # 19. V
    dV <- -r21

    # 20. Va
    dVa <- r21 - r13_fwd + r13_rev - r28

    # 21. Va:Xa (prothrombinase)
    dVa_Xa <- r13_fwd - r13_rev - r14_fwd + r14_rev + r15

    # 22. Va:Xa:II
    dVa_Xa_II <- r14_fwd - r14_rev - r15

    # 23. mIIa (meizothrombin)
    dmIIa <- r11_fwd - k18 * mIIa

    # 24. TFPI
    dTFPI <- -r25_fwd + r25_rev

    # 25. TFPI:Xa
    dTFPI_Xa <- r25_fwd - r25_rev - r26_fwd + r26_rev

    # 26. TFPI:Xa:TF:VIIa
    dTFPI_Xa_TF_VIIa <- r26_fwd - r26_rev

    # 27. ATIII
    dATIII <- -r22 - r23 - r24

    # 28. ATIII:IIa
    dATIII_IIa <- r22

    # 29. ATIII:Xa
    dATIII_Xa <- r23

    # 30. ATIII:IXa
    dATIII_IXa <- r24

    # 31. PC (Protein C)
    dPC <- -r27

    # 32. APC (Activated Protein C)
    dAPC <- r27

    # 33. VIIIa_i (inactivated)
    dVIIIa_i <- r20 + r29

    # 34. Va_i (inactivated)
    dVa_i <- r28

    list(c(
      dTF, dVII, dVIIa, dTF_VII, dTF_VIIa,
      dX, dXa, dTF_VIIa_X, dTF_VIIa_Xa,
      dIX, dIXa, dTF_VIIa_IX,
      dII, dIIa,
      dVIII, dVIIIa, dIXa_VIIIa, dIXa_VIIIa_X,
      dV, dVa, dVa_Xa, dVa_Xa_II,
      dmIIa,
      dTFPI, dTFPI_Xa, dTFPI_Xa_TF_VIIa,
      dATIII, dATIII_IIa, dATIII_Xa, dATIII_IXa,
      dPC, dAPC,
      dVIIIa_i, dVa_i
    ),
    # Auxiliary outputs
    total_thrombin = IIa + 1.2 * mIIa,  # mIIa has ~120% clotting activity
    inhib_factor   = inhib_Xa,
    drug_nM        = drug_conc
    )
  })
}
