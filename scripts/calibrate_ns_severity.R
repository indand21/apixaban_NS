# Calibrate the nephrotic-syndrome coagulation severity scale to the one
# aggregate observation that constrains it: the NS:healthy ETP ratio reported by
# Kelddal et al. 2025 (1096 vs 910 nM.min in patients on 5 mg BID).
#
# The per-species multipliers in data/ns_factor_levels.csv are provisional
# literature-directed guesses. Applied unshrunk they gave an NS_Severe:Normal
# ETP ratio of 2.93 against the observed 1.20. This script finds the single
# common shrinkage s in modifier' = 1 + s x (modifier - 1) that reproduces the
# observed ratio.
#
# ONE ratio identifies ONE scalar. This calibration cannot distinguish an
# over-aggressive factor table from a cascade whose ETP is over-sensitive to
# that table, and it does not validate the relative weighting between species.
source('scripts/load_project.R')

obs <- read.csv('data/external_reference.csv', stringsAsFactors = FALSE)
k <- obs[obs$study == 'Kelddal2025' & obs$metric == 'ETP_nM_min', ]
etp_healthy <- k$value[k$cohort == 'Healthy']
etp_ns      <- k$value[k$cohort == 'NS']
target_ratio <- etp_ns / etp_healthy
cat(sprintf('Observed ETP: healthy %.0f, NS %.0f, ratio %.4f\n',
            etp_healthy, etp_ns, target_ratio))

# Conditions matching the Kelddal rows of the external comparison
etp_for <- function(stage, egfr) {
  prep <- prepare_patient(reference_patient(stage, egfr, 70, 51))
  evaluate_prepared(prep, 5, compute_ti = FALSE)$tga_metrics$ETP_nM_min
}

ratio_at <- function(s) {
  old <- getOption('ns_severity_scale'); on.exit(options(ns_severity_scale = old))
  options(ns_severity_scale = s)
  r <- etp_for('NS_Severe', 81) / etp_for('Normal', 107)
  cat(sprintf('  s = %.5f -> ETP ratio %.4f\n', s, r))
  r
}

r1 <- ratio_at(1)
cat(sprintf('\nUnshrunk (s=1) ratio %.4f vs observed %.4f\n', r1, target_ratio))
if (r1 < target_ratio) {
  stop('Unshrunk ratio already below target; shrinkage cannot reach it.')
}

fit <- uniroot(function(s) ratio_at(s) - target_ratio, interval = c(0, 1),
               tol = 1e-5)
s_hat <- fit$root
cat(sprintf('\nCalibrated severity_scale = %.5f\n', s_hat))

options(ns_severity_scale = s_hat)
shrunk <- get_ckd_modifiers('NS_Severe')
options(ns_severity_scale = NULL)
cat('\nNS_Severe multipliers after shrinkage:\n'); print(round(shrunk, 4))

write.csv(data.frame(parameter = 'severity_scale', value = s_hat,
  source = sprintf('Calibrated on %s to the Kelddal 2025 NS:healthy ETP ratio of %.4f (1096/910 nM.min) under the external-comparison conditions. One aggregate ratio identifies only this one scalar; the relative weighting between species in ns_factor_levels.csv remains an untested assumption.', Sys.Date(), target_ratio)),
  'data/ns_severity_calibration.csv', row.names = FALSE)

dir.create('output/revised', recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(quantity = c('observed_ratio','unshrunk_ratio','calibrated_scale','calibrated_ratio'),
                     value = c(target_ratio, r1, s_hat, ratio_at(s_hat))),
          'output/revised/ns_severity_fit.csv', row.names = FALSE)
cat('\nNS SEVERITY CALIBRATION COMPLETE\n')
