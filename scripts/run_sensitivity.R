stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_sensitivity.R
# One-at-a-time (OAT) sensitivity analysis
# Perturb key parameters ±50% and measure impact on TGA metrics

rm(list = ls())

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

# --- Configuration ---
dose_mg    <- 5.0
egfr       <- 120
ckd_stage  <- "Normal"
perturb    <- 0.50  # ±50%

# Parameters to perturb and their module
sensitivity_params <- list(
  # PBPK parameters
  list(name = "Ki_apixaban", module = "coag",  default = 0.08),
  list(name = "ka",          module = "pbpk",  default = 1.0),
  list(name = "CLint_renal",   module = "pbpk",  default = NULL),  # Will use from params
  list(name = "CLint_hepatic", module = "pbpk",  default = NULL),
  list(name = "fu_plasma",   module = "pbpk",  default = NULL),
  # Key Hockin-Mann rate constants
  list(name = "k23",         module = "coag",  default = 1344.0),  # Prothrombinase kcat
  list(name = "k28",         module = "coag",  default = 8.2),     # Tenase kcat
  list(name = "k32",         module = "coag",  default = 7.1e-6),  # ATIII-IIa
  list(name = "k33",         module = "coag",  default = 2.3e-5),  # ATIII-Xa
  list(name = "k39",         module = "coag",  default = 6.2e-6),  # PC activation
  # Platelet model parameters
  list(name = "k_thr",       module = "platelet", default = 0.020),   # Thrombin-platelet activation
  list(name = "k_agg",       module = "platelet", default = 0.005),   # Aggregation rate
  list(name = "k_txa_syn",   module = "platelet", default = 50.0),    # TxA2 synthesis
  list(name = "k_rel_adp",   module = "platelet", default = 200.0)    # ADP release
)

# --- Baseline ---
cat("Running baseline simulation...\n")
clear_baseline_cache()
baseline <- run_linked_simulation(dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage,
                                   compute_ti = TRUE)
baseline_metrics <- baseline$tga_metrics
baseline_hemo <- baseline$hemostatic_metrics
baseline_ti <- baseline$therapeutic_index
cat(sprintf("Baseline ETP=%.0f nM*s, Peak=%.1f nM, Lag=%.1f s\n",
            baseline_metrics$ETP, baseline_metrics$peak_thrombin, baseline_metrics$lag_time))
cat(sprintf("Baseline HC=%.0f nM*s, TI=%.3f\n\n",
            baseline_hemo$HC, baseline_ti$TI))

# --- OAT Sensitivity ---
cat("Running sensitivity analysis (±50%)...\n")
cat("=========================================\n\n")

sa_results <- list()

for (sp in sensitivity_params) {
  param_name <- sp$name

  for (direction in c("low", "high")) {
    mult <- if (direction == "low") (1 - perturb) else (1 + perturb)
    label <- paste0(param_name, "_", direction)

    cat(sprintf("  %s (x%.2f)...", label, mult))

    tryCatch({
      coag_res <- NULL
      plt_params_perturbed <- NULL

      if (sp$module == "pbpk") {
        # Perturb PBPK parameter
        params <- get_pbpk_params(egfr = egfr, ckd_stage = ckd_stage)
        if (!(param_name %in% names(params))) {
          stop("Sensitivity parameter '", param_name, "' is not in the PBPK ",
               "parameter list, so the perturbation would be silently ignored. ",
               "Note that clearance parameters were renamed to CLint_renal / ",
               "CLint_hepatic when elimination was made fu-dependent.")
        }
        params[[param_name]] <- params[[param_name]] * mult
        # Both CLint_* and fu_plasma feed the derived total-drug clearances
        params <- recompute_derived_params(params)
        pbpk_res <- simulate_pbpk(dose_mg = dose_mg, egfr = egfr, params = params)
        Cp_free_nM <- extract_trough_nM(pbpk_res)

        ic <- get_coag_initial_conditions()
        coag_res <- simulate_coagulation(ic = ic, Cp_free_nM = Cp_free_nM)
        tga <- compute_tga_metrics(coag_res)

      } else if (sp$module == "coag") {
        # Perturb coagulation parameter
        pbpk_res <- simulate_pbpk(dose_mg = dose_mg, egfr = egfr)
        Cp_free_nM <- extract_trough_nM(pbpk_res)

        rc <- get_coag_rate_constants()
        Ki <- 0.08

        if (param_name == "Ki_apixaban") {
          Ki <- Ki * mult
        } else if (param_name %in% names(rc)) {
          rc[[param_name]] <- rc[[param_name]] * mult
        }

        ic <- get_coag_initial_conditions()
        coag_res <- simulate_coagulation(
          ic = ic, rate_constants = rc,
          Cp_free_nM = Cp_free_nM, Ki_apixaban = Ki
        )
        tga <- compute_tga_metrics(coag_res)

      } else if (sp$module == "platelet") {
        # Perturb platelet parameter — coag unchanged, platelet params modified
        pbpk_res <- simulate_pbpk(dose_mg = dose_mg, egfr = egfr)
        Cp_free_nM <- extract_trough_nM(pbpk_res)
        ic <- get_coag_initial_conditions()
        coag_res <- simulate_coagulation(ic = ic, Cp_free_nM = Cp_free_nM)
        tga <- compute_tga_metrics(coag_res)

        plt_params_perturbed <- get_platelet_params()
        if (param_name %in% names(plt_params_perturbed)) {
          plt_params_perturbed[[param_name]] <- plt_params_perturbed[[param_name]] * mult
        }
      }

      # Run platelet model on the coag result
      if (!is.null(coag_res)) {
        plt_res <- simulate_platelets(coag_res, ckd_stage = ckd_stage,
                                       plt_params = plt_params_perturbed)
        hemo <- compute_hemostatic_metrics(plt_res)
      } else {
        hemo <- list(HC = NA, peak_PLT_agg = NA)
      }

      sa_results[[label]] <- data.frame(
        parameter      = param_name,
        direction      = direction,
        multiplier     = mult,
        ETP            = tga$ETP,
        peak_thrombin  = tga$peak_thrombin,
        lag_time       = tga$lag_time,
        HC             = hemo$HC,
        ETP_pct_change = (tga$ETP - baseline_metrics$ETP) / baseline_metrics$ETP * 100,
        peak_pct_change = (tga$peak_thrombin - baseline_metrics$peak_thrombin) /
                          baseline_metrics$peak_thrombin * 100,
        lag_pct_change  = if (!is.na(tga$lag_time) & !is.na(baseline_metrics$lag_time))
                            (tga$lag_time - baseline_metrics$lag_time) /
                            baseline_metrics$lag_time * 100
                          else NA,
        HC_pct_change   = if (!is.na(hemo$HC) & baseline_hemo$HC > 0)
                            (hemo$HC - baseline_hemo$HC) / baseline_hemo$HC * 100
                          else NA,
        stringsAsFactors = FALSE
      )

      cat(" done\n")

    }, error = function(e) {
      cat(sprintf(" ERROR: %s\n", e$message))
    })
  }
}

# --- Compile and save ---
sa_df <- do.call(rbind, sa_results)
rownames(sa_df) <- NULL

cat("\n=== SENSITIVITY ANALYSIS RESULTS ===\n")
print(sa_df[, c("parameter", "direction", "ETP_pct_change", "peak_pct_change",
                 "lag_pct_change", "HC_pct_change")],
      digits = 2)

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(sa_df, "output/tables/sensitivity_analysis.csv", row.names = FALSE)

# --- Tornado plot ---
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Compute range for each parameter
sa_range <- sa_df %>%
  group_by(parameter) %>%
  summarize(
    ETP_range = diff(range(ETP_pct_change)),
    peak_range = diff(range(peak_pct_change)),
    .groups = "drop"
  ) %>%
  arrange(desc(ETP_range))

# Tornado plot for ETP
p_tornado <- sa_df %>%
  mutate(parameter = factor(parameter, levels = rev(sa_range$parameter))) %>%
  ggplot(aes(x = ETP_pct_change, y = parameter, fill = direction)) +
  geom_col(position = "identity", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid") +
  scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
  labs(x = "% Change in ETP", y = "Parameter",
       title = "Sensitivity of ETP to ±50% parameter perturbation",
       fill = "Direction") +
  theme_bw()

ggsave("output/figures/sensitivity_tornado_etp.png", p_tornado, width = 10, height = 6, dpi = 300)

# Tornado plot for HC
sa_hc <- sa_df[!is.na(sa_df$HC_pct_change), ]
if (nrow(sa_hc) > 0) {
  sa_hc_range <- sa_hc %>%
    group_by(parameter) %>%
    summarize(HC_range = diff(range(HC_pct_change)), .groups = "drop") %>%
    arrange(desc(HC_range))

  p_tornado_hc <- sa_hc %>%
    mutate(parameter = factor(parameter, levels = rev(sa_hc_range$parameter))) %>%
    ggplot(aes(x = HC_pct_change, y = parameter, fill = direction)) +
    geom_col(position = "identity", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid") +
    scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
    labs(x = "% Change in Hemostatic Capacity", y = "Parameter",
         title = "Sensitivity of HC to +/-50% parameter perturbation",
         fill = "Direction") +
    theme_bw()

  ggsave("output/figures/sensitivity_tornado_hc.png", p_tornado_hc, width = 10, height = 6, dpi = 300)
}

cat("\nUncoupled sensitivity complete.\n")

# =================================================================
# COUPLED Sensitivity Analysis
# =================================================================
cat("\n=================================================================\n")
cat("COUPLED sensitivity analysis: PS_exp -> prothrombinase feedback\n")
cat("=================================================================\n\n")

# Coupled sensitivity parameters: all uncoupled params + coupling params
coupled_sensitivity_params <- c(
  sensitivity_params,
  list(
    list(name = "alpha_assembly",  module = "coupling", default = NULL),
    list(name = "alpha_catalysis", module = "coupling", default = NULL),
    list(name = "Km_PS",           module = "coupling", default = NULL),
    list(name = "n_ps",            module = "coupling", default = NULL),
    list(name = "k_ps",            module = "platelet", default = NULL)  # PS exposure rate
  )
)

# --- Coupled baseline ---
cat("Running coupled baseline...\n")
clear_coupled_baseline_cache()
coupled_baseline <- run_coupled_linked_simulation(
  dose_mg = dose_mg, egfr = egfr, ckd_stage = ckd_stage, compute_ti = TRUE
)
cb_metrics <- coupled_baseline$tga_metrics
cb_hemo    <- coupled_baseline$hemostatic_metrics
cb_ti      <- coupled_baseline$therapeutic_index
cat(sprintf("Coupled baseline ETP=%.0f nM*s, Peak=%.1f nM, Lag=%.1f s\n",
            cb_metrics$ETP, cb_metrics$peak_thrombin, cb_metrics$lag_time))
cat(sprintf("Coupled baseline HC=%.0f nM*s, TI=%.3f\n\n",
            cb_hemo$HC, cb_ti$TI))

# --- Helper: run a single coupled sensitivity sim ---
run_coupled_sa <- function(param_name, module, mult) {
  # Get all default parameters
  pbpk_params     <- get_pbpk_params(egfr = egfr, ckd_stage = ckd_stage)
  rate_constants  <- get_coag_rate_constants()
  plt_params      <- get_platelet_params()
  plt_params      <- apply_ckd_platelet_modifiers(plt_params, ckd_stage)
  coupling_params <- get_coupling_params()
  Ki <- 0.08

  # Apply perturbation to the right module
  if (module == "pbpk") {
    if (param_name %in% names(pbpk_params)) {
      pbpk_params[[param_name]] <- pbpk_params[[param_name]] * mult
    }
  } else if (module == "coag") {
    if (param_name == "Ki_apixaban") {
      Ki <- Ki * mult
    } else if (param_name %in% names(rate_constants)) {
      rate_constants[[param_name]] <- rate_constants[[param_name]] * mult
    }
  } else if (module == "platelet") {
    if (param_name %in% names(plt_params)) {
      plt_params[[param_name]] <- plt_params[[param_name]] * mult
    }
  } else if (module == "coupling") {
    if (param_name %in% names(coupling_params)) {
      coupling_params[[param_name]] <- coupling_params[[param_name]] * mult
    }
  }

  # Run PBPK
  pbpk_res <- simulate_pbpk(dose_mg = dose_mg, egfr = egfr, params = pbpk_params)
  Cp_free_nM <- extract_trough_nM(pbpk_res)

  # Coag ICs
  ic_coag <- get_coag_initial_conditions()
  ic_coag <- apply_ckd_modifiers(ic_coag, ckd_stage = ckd_stage)

  # Platelet ICs
  ic_plt <- get_platelet_initial_conditions(plt_params)

  # Run coupled ODE
  coupled_out <- simulate_coupled(
    ic_coag         = ic_coag,
    ic_plt          = ic_plt,
    rate_constants  = rate_constants,
    plt_params      = plt_params,
    coupling_params = coupling_params,
    Cp_free_nM      = Cp_free_nM,
    Ki_apixaban     = Ki
  )

  tga  <- compute_tga_metrics(coupled_out$coag_result)
  hemo <- compute_hemostatic_metrics(coupled_out$plt_result)

  list(tga = tga, hemo = hemo)
}

# --- Run coupled OAT ---
cat("Running coupled OAT sensitivity (±50%)...\n\n")
coupled_sa_results <- list()

for (sp in coupled_sensitivity_params) {
  param_name <- sp$name

  for (direction in c("low", "high")) {
    mult <- if (direction == "low") (1 - perturb) else (1 + perturb)
    label <- paste0(param_name, "_", direction)

    cat(sprintf("  %s (x%.2f)...", label, mult))

    tryCatch({
      res <- run_coupled_sa(param_name, sp$module, mult)
      tga  <- res$tga
      hemo <- res$hemo

      coupled_sa_results[[label]] <- data.frame(
        parameter       = param_name,
        direction       = direction,
        multiplier      = mult,
        ETP             = tga$ETP,
        peak_thrombin   = tga$peak_thrombin,
        lag_time        = tga$lag_time,
        HC              = hemo$HC,
        ETP_pct_change  = (tga$ETP - cb_metrics$ETP) / cb_metrics$ETP * 100,
        peak_pct_change = (tga$peak_thrombin - cb_metrics$peak_thrombin) /
                          cb_metrics$peak_thrombin * 100,
        lag_pct_change  = if (!is.na(tga$lag_time) & !is.na(cb_metrics$lag_time))
                            (tga$lag_time - cb_metrics$lag_time) /
                            cb_metrics$lag_time * 100
                          else NA,
        HC_pct_change   = if (!is.na(hemo$HC) & cb_hemo$HC > 0)
                            (hemo$HC - cb_hemo$HC) / cb_hemo$HC * 100
                          else NA,
        stringsAsFactors = FALSE
      )

      cat(" done\n")
    }, error = function(e) {
      cat(sprintf(" ERROR: %s\n", e$message))
    })
  }
}

coupled_sa_df <- do.call(rbind, coupled_sa_results)
rownames(coupled_sa_df) <- NULL

cat("\n=== COUPLED SENSITIVITY ANALYSIS RESULTS ===\n")
print(coupled_sa_df[, c("parameter", "direction", "ETP_pct_change", "peak_pct_change",
                          "lag_pct_change", "HC_pct_change")],
      digits = 2)

write.csv(coupled_sa_df, "output/tables/coupled_sensitivity_analysis.csv", row.names = FALSE)

# --- Coupled tornado plots ---
coupled_sa_range <- coupled_sa_df %>%
  group_by(parameter) %>%
  summarize(
    ETP_range = diff(range(ETP_pct_change)),
    peak_range = diff(range(peak_pct_change)),
    .groups = "drop"
  ) %>%
  arrange(desc(ETP_range))

p_tornado_coupled <- coupled_sa_df %>%
  mutate(parameter = factor(parameter, levels = rev(coupled_sa_range$parameter))) %>%
  ggplot(aes(x = ETP_pct_change, y = parameter, fill = direction)) +
  geom_col(position = "identity", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid") +
  scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
  labs(x = "% Change in ETP", y = "Parameter",
       title = "COUPLED: Sensitivity of ETP to +/-50% parameter perturbation",
       fill = "Direction") +
  theme_bw()

ggsave("output/figures/coupled_sensitivity_tornado_etp.png", p_tornado_coupled,
       width = 10, height = 7, dpi = 300)

coupled_sa_hc <- coupled_sa_df[!is.na(coupled_sa_df$HC_pct_change), ]
if (nrow(coupled_sa_hc) > 0) {
  coupled_hc_range <- coupled_sa_hc %>%
    group_by(parameter) %>%
    summarize(HC_range = diff(range(HC_pct_change)), .groups = "drop") %>%
    arrange(desc(HC_range))

  p_tornado_hc_coupled <- coupled_sa_hc %>%
    mutate(parameter = factor(parameter, levels = rev(coupled_hc_range$parameter))) %>%
    ggplot(aes(x = HC_pct_change, y = parameter, fill = direction)) +
    geom_col(position = "identity", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid") +
    scale_fill_manual(values = c("low" = "steelblue", "high" = "coral")) +
    labs(x = "% Change in Hemostatic Capacity", y = "Parameter",
         title = "COUPLED: Sensitivity of HC to +/-50% parameter perturbation",
         fill = "Direction") +
    theme_bw()

  ggsave("output/figures/coupled_sensitivity_tornado_hc.png", p_tornado_hc_coupled,
         width = 10, height = 7, dpi = 300)
}

cat("\nAll results saved to output/\n")
