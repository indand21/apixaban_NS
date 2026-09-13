# test_coupling.R
# Unit tests for two-way coupled coagulation-platelet system

source("R/00_packages.R")
source("R/01_utils.R")
source("R/10_pbpk_parameters.R")
source("R/11_pbpk_odes.R")
source("R/12_pbpk_simulate.R")
source("R/20_coag_parameters.R")
source("R/21_coag_odes.R")
source("R/22_coag_ckd_modifiers.R")
source("R/23_coag_simulate.R")
source("R/24_platelet_model.R")
source("R/25_platelet_simulate.R")
source("R/26_coupled_odes.R")
source("R/27_coupled_simulate.R")
source("R/30_linkage.R")
source("R/40_tga_analysis.R")
source("R/41_pk_analysis.R")
source("R/42_bleeding_risk.R")

cat("=== Coupled System Tests ===\n\n")
pass <- 0; fail <- 0

# --- Helper: run a quick coupled sim ---
run_coupled_test <- function(Cp_free_nM = 0, ckd_stage = "Normal", t_end = 1200) {
  ic_coag <- get_coag_initial_conditions(TF_pM = 25)
  ic_coag <- apply_ckd_modifiers(ic_coag, ckd_stage = ckd_stage)

  plt_params <- get_platelet_params()
  plt_params <- apply_ckd_platelet_modifiers(plt_params, ckd_stage)
  ic_plt <- get_platelet_initial_conditions(plt_params)

  coupling_params <- get_coupling_params()
  rate_constants  <- get_coag_rate_constants()

  simulate_coupled(
    ic_coag         = ic_coag,
    ic_plt          = ic_plt,
    rate_constants  = rate_constants,
    plt_params      = plt_params,
    coupling_params = coupling_params,
    Cp_free_nM      = Cp_free_nM,
    t_end           = t_end
  )
}

# Test 1: Coupled ODE runs to completion without NaN/Inf
cat("Test 1: Coupled ODE runs to completion (no NaN/Inf)... ")
tryCatch({
  out <- run_coupled_test()
  has_nan <- any(is.nan(unlist(out$full_result)))
  has_inf <- any(is.infinite(unlist(out$full_result[, -1])))  # skip time col
  if (!has_nan && !has_inf) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    cat("FAIL (NaN or Inf detected)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 2: PS_exp > 0 in coupled sim (platelets expose PS)
cat("Test 2: PS_exp > 0 in coupled simulation... ")
tryCatch({
  out <- run_coupled_test()
  max_ps <- max(out$plt_result$PS_exp)
  if (max_ps > 0) {
    cat(sprintf("PASS (max PS_exp=%.2f nM)\n", max_ps)); pass <- pass + 1
  } else {
    cat("FAIL (PS_exp never exceeds 0)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 3: Coupled IIa peak >= uncoupled IIa peak (PS feedback amplifies thrombin)
cat("Test 3: Coupled IIa peak >= uncoupled IIa peak... ")
tryCatch({
  # Coupled
  out_coupled <- run_coupled_test()
  peak_coupled <- max(out_coupled$coag_result$IIa)

  # Uncoupled
  ic_coag <- get_coag_initial_conditions(TF_pM = 25)
  out_uncoupled <- simulate_coagulation(ic = ic_coag, Cp_free_nM = 0, t_end = 1200)
  peak_uncoupled <- max(out_uncoupled$IIa)

  if (peak_coupled >= peak_uncoupled) {
    cat(sprintf("PASS (coupled=%.1f, uncoupled=%.1f nM)\n",
                peak_coupled, peak_uncoupled))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (coupled=%.1f < uncoupled=%.1f nM)\n",
                peak_coupled, peak_uncoupled))
    fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 4: Drug still suppresses IIa in coupled mode
cat("Test 4: Drug suppresses IIa in coupled mode... ")
tryCatch({
  out_nodrug <- run_coupled_test(Cp_free_nM = 0)
  out_drug   <- run_coupled_test(Cp_free_nM = 200)  # ~Ctrough at 5mg

  peak_nodrug <- max(out_nodrug$coag_result$IIa)
  peak_drug   <- max(out_drug$coag_result$IIa)

  if (peak_drug < peak_nodrug) {
    cat(sprintf("PASS (no-drug=%.1f, drug=%.1f nM)\n",
                peak_nodrug, peak_drug))
    pass <- pass + 1
  } else {
    cat("FAIL (drug did not suppress IIa)\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 5: Severe NS increases coupled thrombin generation while preserving platelet output
cat("Test 5: NS_Severe increases coupled IIa and has positive platelet output... ")
tryCatch({
  out_normal <- run_coupled_test(ckd_stage = "Normal")
  out_ns <- run_coupled_test(ckd_stage = "NS_Severe")

  ps_normal <- max(out_normal$plt_result$PS_exp)
  ps_ns <- max(out_ns$plt_result$PS_exp)
  iia_normal <- max(out_normal$coag_result$IIa)
  iia_ns <- max(out_ns$coag_result$IIa)
  hc_ns <- compute_hemostatic_metrics(out_ns$plt_result)$HC

  if (iia_ns > iia_normal && ps_ns > 0 && hc_ns > 0) {
    cat(sprintf("PASS (IIa %.1f -> %.1f nM, PS=%.1f nM, HC=%.0f)\n",
                iia_normal, iia_ns, ps_ns, hc_ns))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (PS %.1f -> %.1f, IIa %.1f -> %.1f, HC=%.0f)\n",
                ps_normal, ps_ns, iia_normal, iia_ns, hc_ns))
    fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 6: PS enhancement factor is bounded (Michaelis-Menten saturates)
cat("Test 6: PS enhancement is bounded (saturates at alpha+1)... ")
tryCatch({
  cp <- get_coupling_params()
  # At very high PS_exp, enhancement should approach 1 + alpha
  enh_high <- compute_ps_enhancement(1e6, cp$alpha_catalysis, cp$Km_PS, cp$n_ps)
  enh_zero <- compute_ps_enhancement(0, cp$alpha_catalysis, cp$Km_PS, cp$n_ps)
  max_enh <- 1 + cp$alpha_catalysis

  if (enh_zero == 1.0 && abs(enh_high - max_enh) < 0.01) {
    cat(sprintf("PASS (at 0: %.1f, at 1e6: %.1f, max: %.1f)\n",
                enh_zero, enh_high, max_enh))
    pass <- pass + 1
  } else {
    cat(sprintf("FAIL (at 0: %.3f, at 1e6: %.3f)\n", enh_zero, enh_high))
    fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 7: No negative species concentrations
cat("Test 7: No negative species concentrations... ")
tryCatch({
  out <- run_coupled_test()
  # Check all 40 species columns (skip time and aux columns)
  species_cols <- c(.coag_species_names, .plt_species_names)
  all_positive <- TRUE
  for (col in species_cols) {
    if (any(out$full_result[[col]] < -1e-8)) {
      cat(sprintf("FAIL (negative values in %s, min=%.2e)\n",
                  col, min(out$full_result[[col]])))
      all_positive <- FALSE
      break
    }
  }
  if (all_positive) {
    cat("PASS\n"); pass <- pass + 1
  } else {
    fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 8: HC and TI computable from coupled results
cat("Test 8: HC and TI computable from coupled results... ")
tryCatch({
  clear_coupled_baseline_cache()
  res <- run_coupled_linked_simulation(
    dose_mg = 5, egfr = 120, ckd_stage = "Normal", compute_ti = TRUE
  )
  if (!is.null(res$hemostatic_metrics) && res$hemostatic_metrics$HC > 0 &&
      !is.null(res$therapeutic_index) && res$therapeutic_index$TI > 0) {
    cat(sprintf("PASS (HC=%.0f, TI=%.3f)\n",
                res$hemostatic_metrics$HC, res$therapeutic_index$TI))
    pass <- pass + 1
  } else {
    cat("FAIL\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 9: Coupled result has 40 species
cat("Test 9: Coupled result has 40 species... ")
tryCatch({
  out <- run_coupled_test()
  n_species <- ncol(out$full_result) - 1  # subtract time column
  # 40 species + auxiliary columns
  if (n_species >= 40) {
    cat(sprintf("PASS (%d columns including aux)\n", n_species)); pass <- pass + 1
  } else {
    cat(sprintf("FAIL (%d species, expected >= 40)\n", n_species)); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

# Test 10: Coupled simulation for NS_Severe completes
cat("Test 10: Coupled simulation for NS_Severe completes... ")
tryCatch({
  out <- run_coupled_test(ckd_stage = "NS_Severe")
  if (nrow(out$coag_result) > 0 && max(out$coag_result$IIa) > 0) {
    cat(sprintf("PASS (peak IIa=%.1f nM)\n", max(out$coag_result$IIa)))
    pass <- pass + 1
  } else {
    cat("FAIL\n"); fail <- fail + 1
  }
}, error = function(e) {
  cat(sprintf("FAIL (%s)\n", e$message)); fail <<- fail + 1
})

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", pass, fail))
