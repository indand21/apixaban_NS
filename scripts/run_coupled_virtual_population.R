stop('Retired legacy workflow. Use scripts/run_ns_pipeline.R; archived dose and validation claims are not current.')
# run_coupled_virtual_population.R
# WP5: Coupled model virtual population dose optimization
# Uses the 40-ODE coupled coagulation-platelet system with PS feedback

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
cat("   WP5: Coupled Virtual Population Dose Optimization\n")
cat("   40-ODE coupled system with PS-prothrombinase feedback\n")
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

# Worker function for coupled model
optimize_one_patient_coupled <- function(args) {
  stage <- args$stage
  egfr  <- args$egfr
  seed  <- args$seed
  
  patient <- sample_virtual_patient(ckd_stage = stage, egfr = egfr, seed = seed)
  
  # TI objective function using coupled model
  ti_obj <- function(dose_mg) {
    res <- tryCatch(
      run_coupled_patient_simulation(patient, dose_mg, TF_pM = TF_pM_RUN,
                                      compute_ti = TRUE),
      error = function(e) NULL
    )
    if (is.null(res) || is.null(res$therapeutic_index)) return(-Inf)
    res$therapeutic_index$TI
  }
  
  # Optimize
  opt <- tryCatch(
    optimize(ti_obj, interval = c(0.1, 20), maximum = TRUE, tol = 0.1),
    error = function(e) list(maximum = NA, objective = NA)
  )
  
  if (is.na(opt$maximum)) {
    return(data.frame(
      ckd_stage      = stage,
      egfr           = egfr,
      patient_id     = seed,
      body_weight    = patient$body_weight,
      age            = patient$age,
      optimal_dose   = NA,
      max_TI         = NA,
      efficacy_score = NA,
      safety_score   = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  # Verify optimal dose
  verify <- tryCatch(
    run_coupled_patient_simulation(patient, opt$maximum, TF_pM = TF_pM_RUN,
                                    compute_ti = TRUE),
    error = function(e) NULL
  )
  
  data.frame(
    ckd_stage      = stage,
    egfr           = egfr,
    patient_id     = seed,
    body_weight    = patient$body_weight,
    age            = patient$age,
    optimal_dose   = opt$maximum,
    max_TI         = if (!is.null(verify)) verify$therapeutic_index$TI else NA,
    efficacy_score = if (!is.null(verify)) verify$therapeutic_index$efficacy_score else NA,
    safety_score   = if (!is.null(verify)) verify$therapeutic_index$safety_score else NA,
    stringsAsFactors = FALSE
  )
}

# Build task list
tasks <- list()
for (stage in ckd_stages) {
  for (i in 1:N_PATIENTS) {
    tasks[[length(tasks) + 1]] <- list(
      stage = stage,
      egfr = egfr_values[stage],
      seed = (i - 1) * 1000 + which(ckd_stages == stage)
    )
  }
}

# Export functions and run in parallel
cluster <- parallel::makePSOCKcluster(n_cores)
parallel::clusterExport(cluster, c("tasks", "TF_pM_RUN"), envir = environment())
parallel::clusterEvalQ(cluster, {
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
})

cat("Starting coupled virtual population optimization...\n")
cat(sprintf("Total tasks: %d | Workers: %d\n\n", length(tasks), n_cores))

# Run optimization
results <- parallel::parLapply(cluster, tasks, function(task) {
  optimize_one_patient_coupled(task)
})

parallel::stopCluster(cluster)

# Combine results
all_results <- do.call(rbind, results)

# Compute summary statistics
summary_stats <- aggregate(
  optimal_dose ~ ckd_stage + egfr,
  data = all_results,
  FUN = function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    c(
      n = length(x),
      median = median(x),
      mean = mean(x),
      sd = sd(x),
      q05 = quantile(x, 0.05),
      q25 = quantile(x, 0.25),
      q75 = quantile(x, 0.75),
      q95 = quantile(x, 0.95)
    )
  }
)

# Flatten summary matrix
summary_df <- do.call(data.frame, summary_stats)

# Compute coverage metric
coverage_stats <- aggregate(
  optimal_dose ~ ckd_stage,
  data = all_results,
  FUN = function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA)
    med <- median(x)
    mean(abs(x - med) <= 0.2 * med) * 100
  }
)
names(coverage_stats)[2] <- "coverage_20pct"

summary_df <- merge(summary_df, coverage_stats, by = "ckd_stage")

# Compute median TI, efficacy, safety per stage
ti_summary <- aggregate(
  cbind(max_TI, efficacy_score, safety_score) ~ ckd_stage,
  data = all_results,
  FUN = function(x) median(x, na.rm = TRUE)
)

summary_df <- merge(summary_df, ti_summary, by = "ckd_stage")

# Print results
cat("\n=============================================================\n")
cat("   COUPLED VIRTUAL POPULATION RESULTS\n")
cat("=============================================================\n\n")

print(summary_df[, c("ckd_stage", "optimal_dose.n", "optimal_dose.median",
                     "optimal_dose.mean", "optimal_dose.sd",
                     "optimal_dose.q05", "optimal_dose.q95", "coverage_20pct")])

# Save results
output_dir <- "output/wp5_population"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(all_results,
          file.path(output_dir, "coupled_virtual_population_doses.csv"),
          row.names = FALSE)
write.csv(summary_df,
          file.path(output_dir, "coupled_population_dose_summary.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  output/wp5_population/coupled_virtual_population_doses.csv\n")
cat("  output/wp5_population/coupled_population_dose_summary.csv\n")

# Create visualization
library(ggplot2)

all_results$ckd_stage <- factor(all_results$ckd_stage, levels = ckd_stages)

p_violin <- ggplot(all_results[!is.na(all_results$optimal_dose), ],
                   aes(x = ckd_stage, y = optimal_dose, fill = ckd_stage)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = c("Normal" = "#1b9e77", "NS_Mild" = "#d95f02",
                               "NS_Moderate" = "#7570b3", "NS_Severe" = "#e7298a")) +
  labs(x = "NS Severity", y = "TI-Optimal Dose (mg BID)",
       title = "Coupled Model: TI-Optimal Dose Distributions",
       subtitle = "Virtual Population (n=200 per stage) | 40-ODE coupled system") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "coupled_dose_distributions_violin.png"),
       p_violin, width = 8, height = 6, dpi = 300)

cat("\nFigure saved: output/wp5_population/coupled_dose_distributions_violin.png\n")

cat("\n=============================================================\n")
cat("   COUPLED VIRTUAL POPULATION ANALYSIS COMPLETE\n")
cat("=============================================================\n")
