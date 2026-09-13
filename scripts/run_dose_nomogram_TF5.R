stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_dose_nomogram_TF5.R
# Nephrotic syndrome nomogram regeneration: TF = 5 pM (CAT assay protocol).
# Identical to run_dose_nomogram.R except every run_linked_simulation() call
# explicitly passes TF_pM = 5 (instead of relying on the R/20_coag_parameters.R
# default of 25 pM). Written to align Table 6 + §3.6 with Table 5 (also TF=5).
# Outputs overwrite the canonical nomogram files; the previous TF=25 versions
# have been preserved as *_TF25_default.csv.

rm(list = ls())

TF_pM_RUN <- 5

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
source("R/50_visualization.R")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

cat("=== Dose Nomogram (TF =", TF_pM_RUN, "pM) ===\n")

bw_seq     <- seq(45, 100, by = 10)
age_seq    <- seq(30, 90, by = 10)
ckd_stages <- c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")
egfr_map   <- c(Normal = 100, NS_Mild = 90, NS_Moderate = 90, NS_Severe = 90)

grid <- expand.grid(
  body_weight = bw_seq,
  age         = age_seq,
  ckd_stage   = ckd_stages,
  stringsAsFactors = FALSE
)
grid$egfr <- egfr_map[grid$ckd_stage]

cat(sprintf("Grid size: %d BW x %d ages x %d NS severity groups = %d scenarios\n",
            length(bw_seq), length(age_seq), length(ckd_stages), nrow(grid)))

n_cores <- min(parallel::detectCores() - 1, 12)
cat(sprintf("Using %d parallel workers\n", n_cores))

source_files <- c(
  "R/00_packages.R", "R/01_utils.R",
  "R/10_pbpk_parameters.R", "R/11_pbpk_odes.R", "R/12_pbpk_simulate.R",
  "R/20_coag_parameters.R", "R/21_coag_odes.R", "R/22_coag_ckd_modifiers.R",
  "R/23_coag_simulate.R", "R/24_platelet_model.R", "R/25_platelet_simulate.R",
  "R/26_coupled_odes.R", "R/27_coupled_simulate.R",
  "R/30_linkage.R", "R/40_tga_analysis.R", "R/41_pk_analysis.R",
  "R/42_bleeding_risk.R"
)

optimize_one <- function(row) {
  bw    <- row[["body_weight"]]
  ag    <- row[["age"]]
  egfr  <- row[["egfr"]]
  stage <- row[["ckd_stage"]]

  ti_obj <- function(dose_mg) {
    res <- run_linked_simulation(
      dose_mg = dose_mg, egfr = egfr, ckd_stage = stage,
      body_weight = bw, age = ag,
      TF_pM = TF_pM_RUN,
      compute_ti = TRUE
    )
    res$therapeutic_index$TI
  }

  tryCatch({
    opt <- optimize(ti_obj, interval = c(0.1, 20), maximum = TRUE, tol = 0.05)

    verify <- run_linked_simulation(
      dose_mg = opt$maximum, egfr = egfr, ckd_stage = stage,
      body_weight = bw, age = ag,
      TF_pM = TF_pM_RUN,
      compute_ti = TRUE
    )

    data.frame(
      optimal_dose   = round(opt$maximum, 3),
      TI             = verify$therapeutic_index$TI,
      efficacy_score = verify$therapeutic_index$efficacy_score,
      safety_score   = verify$therapeutic_index$safety_score,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      optimal_dose = NA_real_, TI = NA_real_,
      efficacy_score = NA_real_, safety_score = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}

t_start <- Sys.time()

cl <- parallel::makeCluster(n_cores)
parent_wd <- getwd()
parallel::clusterExport(cl, c("source_files", "parent_wd", "TF_pM_RUN"))
parallel::clusterEvalQ(cl, {
  setwd(parent_wd)
  for (f in source_files) source(f)
  clear_baseline_cache()
  TRUE
})

grid_list <- split(grid, seq_len(nrow(grid)))
cat(sprintf("Dispatching %d scenarios to %d workers...\n", nrow(grid), n_cores))
results <- parallel::parLapply(cl, grid_list, optimize_one)
parallel::stopCluster(cl)

elapsed <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("Parallel optimization complete: %.1f minutes\n", as.numeric(elapsed)))

results_df <- do.call(rbind, results)
grid$optimal_dose   <- results_df$optimal_dose
grid$TI             <- results_df$TI
grid$efficacy_score <- results_df$efficacy_score
grid$safety_score   <- results_df$safety_score

n_success <- sum(!is.na(grid$optimal_dose))
n_fail    <- sum(is.na(grid$optimal_dose))
cat(sprintf("Successful: %d / %d (failed: %d)\n", n_success, nrow(grid), n_fail))

for (stage in ckd_stages) {
  idx <- grid$ckd_stage == stage
  cat(sprintf("  %s: dose range %.2f - %.2f mg\n",
              stage,
              min(grid$optimal_dose[idx], na.rm = TRUE),
              max(grid$optimal_dose[idx], na.rm = TRUE)))
}

write.csv(grid, "output/tables/nomogram_training_data.csv", row.names = FALSE)
cat("Saved: output/tables/nomogram_training_data.csv\n")

# --- GAM fitting ---
cat("\n=== GAM Fitting ===\n")
grid$ckd_f <- factor(grid$ckd_stage, levels = ckd_stages)

gam_models <- list()
for (stage in ckd_stages) {
  sub <- grid[grid$ckd_stage == stage & !is.na(grid$optimal_dose), ]
  gam_models[[stage]] <- gam(
    optimal_dose ~ te(body_weight, age, k = c(5, 5)),
    data = sub, method = "REML"
  )
  s <- summary(gam_models[[stage]])
  cat(sprintf("  %s: R^2=%.4f, deviance explained=%.1f%%, edf=%.1f\n",
              stage, s$r.sq, s$dev.expl * 100, sum(s$edf)))
}

gam_pooled <- gam(
  optimal_dose ~ ckd_f + te(body_weight, age, by = ckd_f, k = c(5, 5)),
  data = grid[!is.na(grid$optimal_dose), ],
  method = "REML"
)
saveRDS(list(per_stage = gam_models, pooled = gam_pooled),
        "output/nomogram_gam_models.rds")
cat("Saved: output/nomogram_gam_models.rds\n")

# --- Linear coefficients ---
cat("\n=== Linear Dosing Equations ===\n")

lm_coefs <- data.frame(
  ckd_stage  = character(),
  intercept  = numeric(),
  coef_bw    = numeric(),
  coef_age   = numeric(),
  R2_linear  = numeric(),
  R2_gam     = numeric(),
  stringsAsFactors = FALSE
)

for (stage in ckd_stages) {
  sub <- grid[grid$ckd_stage == stage & !is.na(grid$optimal_dose), ]
  lm_fit <- lm(optimal_dose ~ body_weight + age, data = sub)
  lm_s   <- summary(lm_fit)
  gam_r2 <- summary(gam_models[[stage]])$r.sq
  coefs <- coef(lm_fit)

  lm_coefs <- rbind(lm_coefs, data.frame(
    ckd_stage  = stage,
    intercept  = round(coefs[1], 4),
    coef_bw    = round(coefs[2], 4),
    coef_age   = round(coefs[3], 4),
    R2_linear  = round(lm_s$r.squared, 4),
    R2_gam     = round(gam_r2, 4),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s: dose = %.3f + %.4f*BW + %.4f*Age  (R^2=%.4f, GAM R^2=%.4f)\n",
              stage, coefs[1], coefs[2], coefs[3], lm_s$r.squared, gam_r2))
}

write.csv(lm_coefs, "output/tables/nomogram_linear_coefficients.csv", row.names = FALSE)
cat("Saved: output/tables/nomogram_linear_coefficients.csv\n")

# --- Quick sanity check: BW=70 / Age=65 / Normal predicted dose ---
ref_normal <- lm_coefs[lm_coefs$ckd_stage == "Normal", ]
ref_dose <- ref_normal$intercept + ref_normal$coef_bw * 70 + ref_normal$coef_age * 65
cat(sprintf("\nReference patient (BW=70 kg, Age=65 yr, Normal): %.2f mg BID\n",
            ref_dose))
cat("(Manuscript §3.6 narrative claims approximately 1.0 mg BID)\n")

# --- Contour plots ---
cat("\n=== Generating Contour Plots ===\n")

pred_grid <- expand.grid(
  body_weight = seq(45, 100, by = 1),
  age         = seq(30, 90, by = 1),
  stringsAsFactors = FALSE
)

all_pred <- list()
for (stage in ckd_stages) {
  pg <- pred_grid
  pg$ckd_stage <- stage
  pg$ckd_f     <- factor(stage, levels = ckd_stages)
  pg$pred_dose <- predict(gam_models[[stage]], newdata = pg)
  all_pred[[stage]] <- pg
}
pred_df <- do.call(rbind, all_pred)
pred_df$ckd_stage <- factor(pred_df$ckd_stage, levels = ckd_stages)

fda_rect <- data.frame(xmin = 80, xmax = 90, ymin = 45, ymax = 60)

p_contour <- ggplot(pred_df, aes(x = age, y = body_weight, z = pred_dose)) +
  geom_tile(aes(fill = pred_dose)) +
  geom_contour(breaks = c(0.5, 1, 2.5), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_contour(breaks = seq(0.25, 5, by = 0.25), color = "grey30", linewidth = 0.3) +
  geom_rect(data = fda_rect, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "red", linewidth = 1, linetype = "solid") +
  annotate("text", x = 85, y = 47, label = "FDA\n2-of-3", color = "red",
           size = 2.5, fontface = "bold") +
  facet_wrap(~ ckd_stage, ncol = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 1, name = "Optimal\nDose (mg)") +
  labs(x = "Age (years)", y = "Body Weight (kg)",
       title = sprintf("Dose Nomogram: GAM-Predicted Optimal Dose (TF = %d pM)", TF_pM_RUN),
       subtitle = "Dashed lines = 0.5, 1, 2.5 mg reference; red box = FDA 2-of-3 region") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        panel.grid = element_blank(),
        legend.position = "right")

ggsave("output/figures/dose_nomogram_contours.png", p_contour,
       width = 18, height = 10, dpi = 300)
cat("Saved: output/figures/dose_nomogram_contours.png\n")

# --- Lookup table ---
cat("\n=== Clinical Lookup Table ===\n")

lookup_bw  <- c(50, 60, 70, 80, 90, 100)
lookup_age <- c(40, 50, 60, 70, 80)

lookup_grid <- expand.grid(body_weight = lookup_bw, age = lookup_age,
                            stringsAsFactors = FALSE)

lookup_results <- list()
for (stage in ckd_stages) {
  lg <- lookup_grid
  lg$ckd_f <- factor(stage, levels = ckd_stages)
  pred <- predict(gam_models[[stage]], newdata = lg)
  lookup_results[[stage]] <- round(pred * 2) / 2
}

lookup_df <- lookup_grid
for (stage in ckd_stages) lookup_df[[stage]] <- lookup_results[[stage]]
lookup_df <- lookup_df[order(lookup_df$body_weight, lookup_df$age), ]

write.csv(lookup_df, "output/tables/dose_nomogram_lookup.csv", row.names = FALSE)
cat("Saved: output/tables/dose_nomogram_lookup.csv\n")

cat("\nClinical Dose Nomogram (mg BID, rounded to 0.5 mg)\n")
cat("===================================================\n")
cat(sprintf("%-6s %-5s | %-7s %-7s %-7s %-7s %-7s %-7s\n",
            "BW", "Age", ckd_stages[1], ckd_stages[2], ckd_stages[3],
            ckd_stages[4], ckd_stages[5], ckd_stages[6]))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_len(nrow(lookup_df))) {
  cat(sprintf("%-6d %-5d | %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f\n",
              lookup_df$body_weight[i], lookup_df$age[i],
              lookup_df[[ckd_stages[1]]][i], lookup_df[[ckd_stages[2]]][i],
              lookup_df[[ckd_stages[3]]][i], lookup_df[[ckd_stages[4]]][i],
              lookup_df[[ckd_stages[5]]][i], lookup_df[[ckd_stages[6]]][i]))
}

cat("\n=== Nomogram Summary (TF =", TF_pM_RUN, "pM) ===\n")
cat(sprintf("Training data: %d scenarios (%d successful)\n",
            nrow(grid), sum(!is.na(grid$optimal_dose))))
cat(sprintf("GAM per-stage R^2: %.4f - %.4f\n",
            min(sapply(gam_models, function(m) summary(m)$r.sq)),
            max(sapply(gam_models, function(m) summary(m)$r.sq))))
cat(sprintf("Linear R^2: %.4f - %.4f\n",
            min(lm_coefs$R2_linear), max(lm_coefs$R2_linear)))
cat(sprintf("Lookup table dose range: %.1f - %.1f mg\n",
            min(unlist(lookup_results)), max(unlist(lookup_results))))
cat("Done.\n")
