# test_pbpk.R
# Unit tests for PBPK module

source("R/00_packages.R")
source("R/01_utils.R")
source("R/10_pbpk_parameters.R")
source("R/11_pbpk_odes.R")
source("R/12_pbpk_simulate.R")
source("R/41_pk_analysis.R")

cat("=== PBPK Module Tests ===\n\n")
pass <- 0; fail <- 0

# Test 1: Parameters return valid list
cat("Test 1: get_pbpk_params returns complete parameter list... ")
params <- get_pbpk_params(egfr = 120)
required_names <- c("MW", "F_oral", "fu_normal", "fu_plasma", "ka", "CL_total",
                    "CL_renal", "CL_hepatic", "CLint_renal", "CLint_hepatic",
                    "CL_proteinuria", "V_plasma", "V_peripheral")
if (all(required_names %in% names(params))) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 2: CL_renal scales with eGFR
cat("Test 2: CL_renal scales with eGFR... ")
p_normal <- get_pbpk_params(egfr = 120)
p_ckd3   <- get_pbpk_params(egfr = 45)
ratio <- p_ckd3$CL_renal / p_normal$CL_renal
expected <- 45 / 120
if (abs(ratio - expected) < 0.01) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat(sprintf("FAIL (ratio=%.3f, expected=%.3f)\n", ratio, expected)); fail <- fail + 1
}

# Test 3: Simulation runs without error
cat("Test 3: simulate_pbpk runs to completion... ")
tryCatch({
  sim <- simulate_pbpk(dose_mg = 5, n_doses = 4, egfr = 120)
  if (nrow(sim) > 0 && max(sim$time) > 0) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (empty result)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 4: No negative concentrations
cat("Test 4: No negative concentrations... ")
sim <- simulate_pbpk(dose_mg = 5, n_doses = 14, egfr = 120)
has_neg <- any(sim$C_total_ng_mL < -1e-6, na.rm = TRUE)
if (!has_neg) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL (negative concentrations found)\n"); fail <- fail + 1
}

# Test 5: Cmax in reasonable range for 5 mg BID
cat("Test 5: Cmax in reasonable range (50-500 ng/mL) for 5 mg BID... ")
pk <- compute_pk_metrics(sim)
if (pk$Cmax_total > 50 && pk$Cmax_total < 500) {
  cat(sprintf("PASS (Cmax=%.1f)\n", pk$Cmax_total)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (Cmax=%.1f)\n", pk$Cmax_total)); fail <- fail + 1
}

# Test 6: AUC increases with reduced CL
cat("Test 6: AUC increases with reduced renal clearance... ")
sim_ckd <- simulate_pbpk(dose_mg = 5, n_doses = 14, egfr = 30)
pk_ckd <- compute_pk_metrics(sim_ckd)
if (pk_ckd$AUCtau_total > pk$AUCtau_total) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 7: NS PK modifiers load correctly
cat("Test 7: NS PK modifiers load for all stages... ")
tryCatch({
  all_ok <- TRUE
  for (stg in c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")) {
    m <- get_ckd_pk_modifiers(stg)
    if (!all(c("fu_plasma", "CLint_hepatic", "CL_proteinuria") %in% names(m))) all_ok <- FALSE
  }
  if (all_ok) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (missing modifier fields)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 8: fu_plasma increases in severe NS
cat("Test 8: fu_plasma increases in NS_Severe vs Normal... ")
p_normal <- get_pbpk_params(egfr = 100, ckd_stage = "Normal")
p_ns_severe <- get_pbpk_params(egfr = 90, ckd_stage = "NS_Severe")
if (p_ns_severe$fu_plasma > p_normal$fu_plasma) {
  cat(sprintf("PASS (%.4f -> %.4f)\n", p_normal$fu_plasma, p_ns_severe$fu_plasma)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 9: proteinuria clearance increases in severe NS
cat("Test 9: CL_proteinuria increases in NS_Severe vs Normal... ")
if (p_ns_severe$CL_proteinuria > p_normal$CL_proteinuria) {
  cat(sprintf("PASS (%.4f -> %.4f)\n", p_normal$CL_proteinuria, p_ns_severe$CL_proteinuria)); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 10: Normal disease stage leaves params unchanged
cat("Test 10: Normal disease stage gives identity modifiers... ")
m_normal <- get_ckd_pk_modifiers("Normal")
if (m_normal$fu_plasma == 1.0 && m_normal$CLint_hepatic == 1.0 &&
    m_normal$CL_proteinuria == 0.0) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# Test 11: BW scaling — 50 kg gets higher AUC than 90 kg
cat("Test 11: 50 kg patient has higher AUC than 90 kg at same dose... ")
tryCatch({
  p50 <- get_pbpk_params(egfr = 120, body_weight = 50)
  p90 <- get_pbpk_params(egfr = 120, body_weight = 90)
  sim50 <- simulate_pbpk(dose_mg = 5, n_doses = 14, egfr = 120, params = p50)
  sim90 <- simulate_pbpk(dose_mg = 5, n_doses = 14, egfr = 120, params = p90)
  pk50 <- compute_pk_metrics(sim50)
  pk90 <- compute_pk_metrics(sim90)
  if (pk50$AUCtau_total > pk90$AUCtau_total) {
    cat(sprintf("PASS (AUC50=%.1f > AUC90=%.1f)\n", pk50$AUCtau_total, pk90$AUCtau_total))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (AUC50=%.1f, AUC90=%.1f)\n", pk50$AUCtau_total, pk90$AUCtau_total))
    fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 12: Age CL decline — age=80 CL_hepatic < age=40
cat("Test 12: Age=80 has lower CL_hepatic than age=40... ")
p40 <- get_pbpk_params(egfr = 120, age = 40)
p80 <- get_pbpk_params(egfr = 120, age = 80)
age_ratio <- p80$CL_hepatic / p40$CL_hepatic
# Expected: max(0.5, 1 - 0.007*40) = 0.72
if (p80$CL_hepatic < p40$CL_hepatic && abs(age_ratio - 0.72) < 0.01) {
  cat(sprintf("PASS (ratio=%.3f, expected ~0.72)\n", age_ratio)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (ratio=%.3f)\n", age_ratio)); fail <- fail + 1
}

# Test 13: Default backward compat — age=NULL skips age scaling
cat("Test 13: age=NULL gives same result as no age param... ")
p_default <- get_pbpk_params(egfr = 120)
p_null    <- get_pbpk_params(egfr = 120, age = NULL)
if (identical(p_default$CL_hepatic, p_null$CL_hepatic) &&
    identical(p_default$V_plasma, p_null$V_plasma)) {
  cat("PASS\n"); pass <- pass + 1
} else {
  cat("FAIL\n"); fail <- fail + 1
}

# ---------------------------------------------------------------------------
# Restrictive-clearance regression tests
#
# Apixaban is a low-extraction drug, so every elimination route except urinary
# loss of albumin-bound drug acts on the UNBOUND fraction only. Consequently
#     AUC_unbound = F x Dose / CLint
# which is INDEPENDENT of fu. Changing plasma protein binding alone must
# therefore leave unbound exposure unchanged while total exposure moves by
# 1/fu. These tests guard that invariant.
#
# They exist because the model previously applied CL_hepatic and CL_renal to
# TOTAL venous concentration, so the nephrotic syndrome fu multiplier
# (data/ns_pk_modifiers.csv) passed straight through to unbound exposure and
# produced a spurious +19% to +32% rise in unbound AUC across NS_Moderate to
# NS_Severe. That rise was the CSV constant, not a mechanism.
# ---------------------------------------------------------------------------

ss_auc <- function(params, dose_mg = 5, tau = 12, n_doses = 14) {
  sim <- simulate_pbpk(dose_mg = dose_mg, tau = tau, n_doses = n_doses,
                       egfr = 120, params = params)
  compute_pk_metrics(sim, tau = tau)
}

p_ref <- get_pbpk_params(egfr = 120, ckd_stage = "Normal")
pk_ref <- ss_auc(p_ref)

# Test 14: fu alone does NOT change unbound exposure (the core invariant)
cat("Test 14: scaling fu alone leaves unbound AUC invariant... ")
p_fu <- p_ref
p_fu$fu_plasma <- p_ref$fu_plasma * 1.45   # same multiplier as NS_Severe
p_fu <- recompute_derived_params(p_fu)
pk_fu <- ss_auc(p_fu)
free_ratio  <- pk_fu$AUCtau_free  / pk_ref$AUCtau_free
total_ratio <- pk_fu$AUCtau_total / pk_ref$AUCtau_total
# CL_proteinuria is 0 at the Normal stage, so both relations are exact
if (abs(free_ratio - 1) < 0.005 && abs(total_ratio - 1 / 1.45) < 0.01) {
  cat(sprintf("PASS (free x%.4f, total x%.4f, expected 1.000 and %.4f)\n",
              free_ratio, total_ratio, 1 / 1.45)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (free x%.4f expected 1.000; total x%.4f expected %.4f)\n",
              free_ratio, total_ratio, 1 / 1.45)); fail <- fail + 1
}

# Test 15: the reparameterisation is a no-op at normal binding
cat("Test 15: effective clearances at fu_normal match the reference values... ")
ok15 <- abs(p_ref$CL_hepatic - 2.40) < 1e-9 &&
        abs(p_ref$CL_renal   - 0.90) < 1e-9 &&
        abs(p_ref$CL_total   - 3.30) < 1e-9 &&
        abs(p_ref$CLint_hepatic - 2.40 / 0.13) < 1e-9
if (ok15) {
  cat(sprintf("PASS (CL_hep=%.3f CL_ren=%.3f CL_tot=%.3f)\n",
              p_ref$CL_hepatic, p_ref$CL_renal, p_ref$CL_total)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (CL_hep=%.6f CL_ren=%.6f CL_tot=%.6f CLint_hep=%.6f)\n",
              p_ref$CL_hepatic, p_ref$CL_renal, p_ref$CL_total,
              p_ref$CLint_hepatic)); fail <- fail + 1
}

# Test 16: unbound exposure responds to CLint_hepatic and NEVER to fu.
# Assertions are derived from data/ns_pk_modifiers.csv rather than hard-coded,
# so the test stays valid when the disease modifiers are revised.
cat("Test 16: unbound AUC responds to CLint, is invariant to fu... ")
tryCatch({
  mods <- get_ckd_pk_modifiers("NS_Severe")
  p_n  <- get_pbpk_params(egfr = 120, ckd_stage = "Normal")
  p_s  <- get_pbpk_params(egfr = 120, ckd_stage = "NS_Severe")

  free_auc <- function(p) {
    compute_pk_metrics(simulate_pbpk(dose_mg = 5, tau = 12, n_doses = 14,
                                     egfr = 120, params = p), tau = 12)$AUCtau_free
  }
  mk <- function(fu_mult, clint_mult, cprot) {
    q <- p_n
    q$fu_plasma      <- p_n$fu_plasma * fu_mult
    q$CLint_hepatic  <- p_n$CLint_hepatic * clint_mult
    q$CL_proteinuria <- cprot
    recompute_derived_params(q)
  }

  base   <- free_auc(mk(1, 1, 0))
  fu_on  <- free_auc(mk(mods$fu_plasma, 1, 0)) / base
  cl_on  <- free_auc(mk(1, mods$CLint_hepatic, 0)) / base
  both   <- free_auc(mk(mods$fu_plasma, mods$CLint_hepatic, 0)) / base
  actual <- free_auc(p_s) / free_auc(p_n)

  # Analytic expectation for scaling only the HEPATIC intrinsic clearance:
  #   AUC_unbound = F x Dose / CLint_total, and CLint_total = CLint_hep + CLint_ren
  expect_cl <- (p_n$CLint_hepatic + p_n$CLint_renal) /
               (p_n$CLint_hepatic * mods$CLint_hepatic + p_n$CLint_renal)

  ok <- abs(fu_on - 1) < 0.005 &&            # fu contributes nothing
        abs(both - cl_on) < 0.005 &&         # adding fu on top changes nothing
        abs(cl_on - expect_cl) < 0.01 &&     # matches restrictive-clearance algebra
        actual <= both + 1e-6                # bound-drug loss can only lower it
  if (ok) {
    cat(sprintf("PASS (fu x%.4f, CLint x%.4f [exp %.4f], both x%.4f, NS_Severe x%.4f)\n",
                fu_on, cl_on, expect_cl, both, actual)); pass <- pass + 1
  } else {
    cat(sprintf("FAIL (fu x%.4f, CLint x%.4f [exp %.4f], both x%.4f, NS_Severe x%.4f)\n",
                fu_on, cl_on, expect_cl, both, actual)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 17: reduced intrinsic clearance is the ONLY lever that raises unbound
# exposure. 0.82 is the magnitude Derebail et al. 2023 imply (free CL/F ~18%
# lower in NS); the resulting unbound AUC rise must be about +15%.
cat("Test 17: lowering CLint_hepatic raises unbound AUC... ")
p_lo <- p_ref
p_lo$CLint_hepatic <- p_lo$CLint_hepatic * 0.82
p_lo <- recompute_derived_params(p_lo)
pk_lo <- ss_auc(p_lo)
lo_free <- pk_lo$AUCtau_free / pk_ref$AUCtau_free
if (lo_free > 1.10 && lo_free < 1.21) {
  cat(sprintf("PASS (unbound %+.1f%%)\n", 100 * (lo_free - 1))); pass <- pass + 1
} else {
  cat(sprintf("FAIL (unbound %+.1f%%, expected about +15%%)\n",
              100 * (lo_free - 1))); fail <- fail + 1
}

# Test 18: urinary loss of BOUND drug is a genuine extra elimination route, so
# it must LOWER unbound exposure. It can never account for a rise.
cat("Test 18: CL_proteinuria lowers unbound AUC... ")
p_pr <- p_ref
p_pr$CL_proteinuria <- 0.20
p_pr <- recompute_derived_params(p_pr)
pk_pr <- ss_auc(p_pr)
pr_free <- pk_pr$AUCtau_free / pk_ref$AUCtau_free
if (pr_free < 0.99) {
  cat(sprintf("PASS (unbound %+.1f%%)\n", 100 * (pr_free - 1))); pass <- pass + 1
} else {
  cat(sprintf("FAIL (unbound %+.1f%%, expected a fall)\n",
              100 * (pr_free - 1))); fail <- fail + 1
}

# Test 19: derived clearances cannot drift from the quantities the ODEs use
cat("Test 19: recompute_derived_params keeps CL_total consistent... ")
p_d <- p_ref
p_d$fu_plasma <- 0.30
p_d$CLint_hepatic <- p_d$CLint_hepatic * 0.9
p_d$CL_proteinuria <- 0.15
p_d <- recompute_derived_params(p_d)
expected_total <- p_d$CLint_hepatic * p_d$fu_plasma +
                  p_d$CLint_renal * p_d$fu_plasma +
                  p_d$CL_proteinuria * (1 - p_d$fu_plasma)
if (abs(p_d$CL_total - expected_total) < 1e-12 &&
    abs(p_d$CL_hepatic - p_d$CLint_hepatic * p_d$fu_plasma) < 1e-12 &&
    abs(p_d$Kp_liver - p_d$Kp_liver_ref * p_d$fu_plasma / p_d$fu_normal) < 1e-12) {
  cat(sprintf("PASS (CL_total=%.4f, Kp_liver=%.4f)\n",
              p_d$CL_total, p_d$Kp_liver)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (CL_total=%.6f expected %.6f; Kp_liver=%.6f expected %.6f)\n",
              p_d$CL_total, expected_total, p_d$Kp_liver,
              p_d$Kp_liver_ref * p_d$fu_plasma / p_d$fu_normal)); fail <- fail + 1
}

vd_ss <- function(p) {
  p$V_plasma + p$V_liver * p$Kp_liver + p$V_kidney * p$Kp_kidney +
    p$V_peripheral * p$Kp_peripheral
}

# Test 20: Kp scales with fu, so volume of distribution rises with fu.
# Tissue accumulation is driven by unbound plasma concentration, so with tissue
# binding unchanged Kp = Kp_ref x fu / fu_normal.
cat("Test 20: Kp and Vd_ss scale with the unbound fraction... ")
p_sv <- get_pbpk_params(egfr = 120, ckd_stage = "NS_Severe")
fu_rel <- p_sv$fu_plasma / p_sv$fu_normal
ok20 <- abs(p_sv$Kp_liver / p_ref$Kp_liver - fu_rel) < 1e-9 &&
        abs(p_sv$Kp_kidney / p_ref$Kp_kidney - fu_rel) < 1e-9 &&
        abs(p_sv$Kp_peripheral / p_ref$Kp_peripheral - fu_rel) < 1e-9 &&
        vd_ss(p_sv) > vd_ss(p_ref) * 1.2
if (ok20) {
  cat(sprintf("PASS (Kp x%.4f = fu x%.4f; Vd_ss %.1f -> %.1f L)\n",
              p_sv$Kp_liver / p_ref$Kp_liver, fu_rel,
              vd_ss(p_ref), vd_ss(p_sv))); pass <- pass + 1
} else {
  cat(sprintf("FAIL (Kp ratios %.4f/%.4f/%.4f vs fu %.4f; Vd_ss %.1f -> %.1f)\n",
              p_sv$Kp_liver / p_ref$Kp_liver, p_sv$Kp_kidney / p_ref$Kp_kidney,
              p_sv$Kp_peripheral / p_ref$Kp_peripheral, fu_rel,
              vd_ss(p_ref), vd_ss(p_sv))); fail <- fail + 1
}

# Test 21: the Kp scaling is a no-op at normal binding, so the healthy baseline
# and all existing validation are unchanged.
cat("Test 21: Kp equals its reference value at fu_normal... ")
expected_vd <- p_ref$V_plasma + p_ref$V_liver * 2.5 + p_ref$V_kidney * 2.0 +
               p_ref$V_peripheral * 1.3
ok21 <- abs(p_ref$Kp_liver - 2.5) < 1e-9 &&
        abs(p_ref$Kp_kidney - 2.0) < 1e-9 &&
        abs(p_ref$Kp_peripheral - 1.3) < 1e-9 &&
        abs(vd_ss(p_ref) - expected_vd) < 1e-9
if (ok21) {
  cat(sprintf("PASS (Kp %.1f/%.1f/%.1f, Vd_ss %.2f L)\n",
              p_ref$Kp_liver, p_ref$Kp_kidney, p_ref$Kp_peripheral,
              vd_ss(p_ref))); pass <- pass + 1
} else {
  cat(sprintf("FAIL (Kp %.4f/%.4f/%.4f, Vd_ss %.4f)\n",
              p_ref$Kp_liver, p_ref$Kp_kidney, p_ref$Kp_peripheral,
              vd_ss(p_ref))); fail <- fail + 1
}

# Test 22: scaling Kp with fu removes the spurious half-life collapse.
# The terminal slope of a mammillary model is governed by
#     lambda_z ~ k10 x k21 / (k10 + k12 + k21)
# with k10 = CL/V_plasma, k12 = Q/V_plasma, k21 = Q/V_tissue_effective.
# Under joint fu-scaling of CL and tissue volumes the NUMERATOR is exactly
# invariant (k10 gains f, k21 loses f) and the denominator moves only modestly,
# so lambda_z - and hence t_half - stays near its normal value. Under the old
# fixed-Kp parameterisation k21 did not fall while k10 rose by f, so lambda_z
# rose by about f and t_half collapsed by about 1/f (-31% at fu_rel = 1.45).
cat("Test 22: terminal half-life is preserved, not collapsed, by fu scaling... ")
thalf <- function(p) {
  st <- c(A_gut = 0, A_liver = 0, A_kidney = 0, A_plasma = 0, A_peripheral = 0)
  ev <- data.frame(var = "A_gut", time = 0, value = 10 * p$F_oral, method = "add")
  o <- as.data.frame(deSolve::ode(y = st, times = seq(0, 72, by = 0.05),
        func = pbpk_odes, parms = p, method = "lsoda",
        events = list(data = ev), atol = 1e-10, rtol = 1e-8))
  tl <- o[o$time >= 48, ]
  lz <- unname(-coef(lm(log(C_total_ng_mL) ~ time, data = tl))[2])
  log(2) / lz
}
t_ref <- thalf(p_ref)
t_fu  <- thalf(p_fu)     # fu x1.45, CLint unchanged, CL_proteinuria = 0
fu_rel22  <- p_fu$fu_plasma / p_fu$fu_normal
fixed_kp  <- 1 / fu_rel22          # what the old parameterisation implied
obs22     <- t_fu / t_ref
# Compare with the exact linear-system terminal eigenvalue, not a false
# assumption that binding leaves half-life exactly unchanged.
eigen_half <- function(p) {
  z <- c(A_gut=0,A_liver=0,A_kidney=0,A_plasma=0,A_peripheral=0)
  mat <- sapply(seq_along(z),function(j){v<-z;v[j]<-1;pbpk_odes(0,v,p)[[1]]})
  log(2)/(-max(Re(eigen(mat)$values)))
}
expected_ratio <- eigen_half(p_fu)/eigen_half(p_ref)
if (abs(obs22 / expected_ratio - 1) < 0.005 && obs22 > fixed_kp) {
  cat(sprintf("PASS (t_half %.2f -> %.2f h, ratio %.4f; fixed-Kp would give %.4f)\n",
              t_ref, t_fu, obs22, fixed_kp)); pass <- pass + 1
} else {
  cat(sprintf("FAIL (t_half ratio %.4f; expected ~1.0, fixed-Kp %.4f)\n",
              obs22, fixed_kp)); fail <- fail + 1
}

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
if (fail > 0) quit(status = 1)
