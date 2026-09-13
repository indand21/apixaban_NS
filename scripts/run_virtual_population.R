stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_virtual_population.R
# WP5: Sample virtual patient populations per nephrotic syndrome severity,
# compute TI-optimal dose for each, report distributions.
# Uses PSOCK parallelism for speed.

rm(list = ls())
setwd(here::here())

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
source("R/60_population_sampler.R")

cat("=============================================================\n")
cat("   WP5: Nephrotic Syndrome Virtual Population Dose Optimization\n")
cat("=============================================================\n\n")

TF_pM_RUN  <- 5
N_PATIENTS <- as.integer(Sys.getenv("NS_N_PATIENTS", "200"))  # per NS severity group
ckd_stages <- c("Normal", "NS_Mild", "NS_Moderate", "NS_Severe")
egfr_values <- c(Normal = 100, NS_Mild = 90, NS_Moderate = 90, NS_Severe = 90)

cat(sprintf("Patients per stage: %d | TF = %d pM\n", N_PATIENTS, TF_pM_RUN))
cat(sprintf("Total simulations: %d dose-finding runs\n\n",
            N_PATIENTS * length(ckd_stages)))

# --- Set up parallel cluster ---
n_cores <- min(parallel::detectCores() - 1, 10)
cat(sprintf("Using %d parallel workers\n\n", n_cores))

source_files <- c(
  "R/00_packages.R", "R/01_utils.R",
  "R/10_pbpk_parameters.R", "R/11_pbpk_odes.R", "R/12_pbpk_simulate.R",
  "R/20_coag_parameters.R", "R/21_coag_odes.R", "R/22_coag_ckd_modifiers.R",
  "R/23_coag_simulate.R", "R/24_platelet_model.R", "R/25_platelet_simulate.R",
  "R/26_coupled_odes.R", "R/27_coupled_simulate.R",
  "R/30_linkage.R", "R/40_tga_analysis.R", "R/41_pk_analysis.R",
  "R/42_bleeding_risk.R", "R/60_population_sampler.R"
)

# Worker function
optimize_one_patient <- function(args) {
  stage <- args$stage
  egfr  <- args$egfr
  seed  <- args$seed

  patient <- sample_virtual_patient(ckd_stage = stage, egfr = egfr, seed = seed)
  opt <- find_patient_optimal_dose(patient, TF_pM = TF_pM_RUN,
                                    dose_interval = c(0.1, 15))

  data.frame(
    ckd_stage      = stage,
    egfr           = egfr,
    patient_id     = seed,
    body_weight    = patient$body_weight,
    age            = patient$age,
    optimal_dose   = opt$optimal_dose,
    max_TI         = opt$max_TI,
    efficacy_score = opt$efficacy_score,
    safety_score   = opt$safety_score,
    stringsAsFactors = FALSE
  )
}

# Build task list
tasks <- list()
counter <- 0
for (stage in ckd_stages) {
  egfr <- egfr_values[stage]
  for (i in seq_len(N_PATIENTS)) {
    counter <- counter + 1
    tasks[[counter]] <- list(
      stage = stage,
      egfr  = egfr,
      seed  = counter * 1000 + i  # reproducible seeds
    )
  }
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

cat(sprintf("Dispatching %d patients to %d workers...\n", length(tasks), n_cores))
results <- parallel::parLapply(cl, tasks, optimize_one_patient)
parallel::stopCluster(cl)

elapsed <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("Population optimization complete: %.1f minutes\n\n", as.numeric(elapsed)))

# Combine results
pop_df <- do.call(rbind, results)
pop_df$ckd_stage <- factor(pop_df$ckd_stage, levels = ckd_stages)

n_success <- sum(!is.na(pop_df$optimal_dose))
n_fail    <- sum(is.na(pop_df$optimal_dose))
cat(sprintf("Successful: %d / %d (failed: %d)\n\n", n_success, nrow(pop_df), n_fail))

# --- Per-stage summary statistics ---
cat("=== Population Dose Distribution Summary ===\n\n")
cat(sprintf("%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %8s\n",
            "Stage", "Median", "Mean", "SD", "Q5", "Q95", "IQR", "Coverage"))
cat(paste(rep("-", 70), collapse = ""), "\n")

pop_summary <- data.frame()

for (stage in ckd_stages) {
  sub <- pop_df[pop_df$ckd_stage == stage & !is.na(pop_df$optimal_dose), ]
  if (nrow(sub) == 0) next

  med  <- median(sub$optimal_dose)
  mn   <- mean(sub$optimal_dose)
  s    <- sd(sub$optimal_dose)
  q5   <- quantile(sub$optimal_dose, 0.05)
  q95  <- quantile(sub$optimal_dose, 0.95)
  iqr  <- IQR(sub$optimal_dose)

  # Coverage: fraction of patients within ±20% of median
  within_20pct <- sum(abs(sub$optimal_dose - med) / med <= 0.20) / nrow(sub) * 100

  cat(sprintf("%-8s  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %7.1f%%\n",
              stage, med, mn, s, q5, q95, iqr, within_20pct))

  pop_summary <- rbind(pop_summary, data.frame(
    ckd_stage = stage,
    n = nrow(sub),
    median_dose = round(med, 3),
    mean_dose = round(mn, 3),
    sd_dose = round(s, 3),
    q05 = round(q5, 3),
    q25 = round(quantile(sub$optimal_dose, 0.25), 3),
    q75 = round(quantile(sub$optimal_dose, 0.75), 3),
    q95 = round(q95, 3),
    coverage_20pct = round(within_20pct, 1),
    median_TI = round(median(sub$max_TI, na.rm = TRUE), 4),
    median_efficacy = round(median(sub$efficacy_score, na.rm = TRUE), 3),
    median_safety = round(median(sub$safety_score, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  ))
}

# --- Save outputs ---
dir.create("output/wp5_population", showWarnings = FALSE, recursive = TRUE)
write.csv(pop_df, "output/wp5_population/virtual_population_doses.csv",
          row.names = FALSE)
write.csv(pop_summary, "output/wp5_population/population_dose_summary.csv",
          row.names = FALSE)

# --- Visualization ---
p_violin <- ggplot(pop_df[!is.na(pop_df$optimal_dose), ],
                    aes(x = ckd_stage, y = optimal_dose, fill = ckd_stage)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.05, 0.5, 0.95)) +
  geom_hline(yintercept = 5.0, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "orange") +
  annotate("text", x = 6.3, y = 5.1, label = "5 mg label", color = "red", size = 3) +
  annotate("text", x = 6.3, y = 2.6, label = "2.5 mg label", color = "orange", size = 3) +
  labs(x = "NS Severity", y = "TI-Optimal Dose (mg BID)",
       title = "WP5: NS Virtual Population TI-Optimal Dose Distributions",
       subtitle = sprintf("N=%d per stage | TF=%d pM | Lines = 5th/50th/95th percentiles",
                          N_PATIENTS, TF_pM_RUN)) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("output/wp5_population/dose_distributions_violin.png", p_violin,
       width = 10, height = 7, dpi = 300)

p_bw_age <- ggplot(pop_df[!is.na(pop_df$optimal_dose) &
                            pop_df$ckd_stage %in% c("Normal", "NS_Moderate", "NS_Severe"), ],
                    aes(x = body_weight, y = optimal_dose, color = factor(round(age/10)*10))) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ ckd_stage, scales = "free_y") +
  labs(x = "Body Weight (kg)", y = "TI-Optimal Dose (mg BID)",
       color = "Age (decade)",
       title = "Dose vs Body Weight by Age Decade",
       subtitle = "Selected NS severity groups") +
  theme_bw()

ggsave("output/wp5_population/dose_vs_bw_by_age.png", p_bw_age,
       width = 14, height = 5, dpi = 300)

cat("\nSaved:\n")
cat("  output/wp5_population/virtual_population_doses.csv\n")
cat("  output/wp5_population/population_dose_summary.csv\n")
cat("  output/wp5_population/dose_distributions_violin.png\n")
cat("  output/wp5_population/dose_vs_bw_by_age.png\n")

cat("\n=============================================================\n")
cat("   WP5 COMPLETE\n")
cat("=============================================================\n")
