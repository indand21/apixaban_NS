# Calibrate absorption parameters to the Frost 2013 aggregate steady-state
# moments (5 mg BID, day 7, n=6).
#
# Systemic clearance is NOT fitted: CL_total stays at the label 3.3 L/h with a
# 0.9 L/h renal component, preserving the September 2026 audit correction.
#
# Only four aggregate moments are available (no digitised concentration-time
# profile), so this is moment matching, not profile fitting, and it yields no
# standard errors.
#
# F_oral is profiled out analytically rather than searched. At steady state
# AUCtau = F x Dose / CL exactly, and AUCtau is independent of the absorption
# parameters, so F is determined in closed form by the observed AUC once CL is
# fixed. That leaves a 2-parameter search (ka, Tlag) against the 3 remaining
# moments (Cmax, Ctrough, Tmax).
#
# ka and Tlag are only weakly separated by aggregate moments: both delay and
# broaden the peak, and with no concentration-time profile the data cannot tell
# a slow first-order input from a delayed faster one. A small ridge penalty
# pulls Tlag towards zero so the lag is only as large as the Tmax observation
# actually requires. The reported lag is therefore a fitted descriptor of an
# observed delay, not an independently identified physiological transit time.
#
# An earlier version of this script also fitted a deep peripheral compartment.
# It was removed after review found it contributed negligibly to the fit
# while making the terminal half-life window-dependent.
source('scripts/load_project.R')

obs <- read.csv('data/external_reference.csv', stringsAsFactors = FALSE)
obs <- obs[obs$study == 'Frost2013', ]
target <- setNames(obs$value, obs$metric)
cat('Targets (Frost 2013, 5 mg BID day 7):\n'); print(target)

F_REF <- 0.50   # F used inside the simulation; rescaled analytically afterwards
TLAG_RANGE  <- c(0, 3)    # hours; must stay well inside the 12 h interval
LAMBDA_TLAG <- 0.01       # ridge weight pulling Tlag towards zero

tlag_from <- function(z) TLAG_RANGE[1] + diff(TLAG_RANGE) * plogis(z)

# Moments at F = F_REF for a given absorption shape. Reference subject matches
# the Frost comparison row in run_revised_analysis.R (eGFR 120, 70 kg, age 40).
sim_shape <- function(ka, Tlag, dt = 0.05, n_doses = 28) {
  p <- get_pbpk_params(egfr = 120, body_weight = 70, ckd_stage = 'Normal', age = 40)
  p$F_oral <- F_REF; p$ka <- ka; p$Tlag <- Tlag
  p <- recompute_derived_params(p)
  s <- simulate_pbpk(dose_mg = 5, tau = 12, n_doses = n_doses, dt = dt, params = p)
  m <- compute_pk_metrics(s, tau = 12)
  c(Cmax_total = m$Cmax_total, Ctrough_total = m$Ctrough_total,
    AUCtau_total = m$AUCtau_total, Tmax = m$Tmax)
}

# Rescale the exposure moments by the F that matches the observed AUC exactly.
scale_to_F <- function(m) {
  F_oral <- F_REF * target['AUCtau_total'] / m['AUCtau_total']
  r <- F_oral / F_REF
  list(F_oral = unname(F_oral),
       moments = c(Cmax_total = unname(m['Cmax_total'] * r),
                   Ctrough_total = unname(m['Ctrough_total'] * r),
                   AUCtau_total = unname(target['AUCtau_total']),
                   Tmax = unname(m['Tmax'])))
}

# Tmax is a median of n=6 with a reported range of 2-4 h, not a precise
# statistic, so it is weighted down rather than treated as exact.
wts <- c(Cmax_total = 1, Ctrough_total = 1, Tmax = 0.35)

cost <- function(theta) {
  Tlag <- tlag_from(theta[2])
  m <- tryCatch(sim_shape(exp(theta[1]), Tlag), error = function(e) NULL)
  if (is.null(m) || any(!is.finite(m))) return(1e6)
  sc <- scale_to_F(m)$moments
  if (sc['Cmax_total'] <= 0 || sc['Ctrough_total'] <= 0) return(1e6)
  sum(wts * (log(sc[names(wts)] / target[names(wts)]))^2) + LAMBDA_TLAG * Tlag^2
}

set.seed(20260909)
start <- c(log(0.70), qlogis(0.5/3))
best <- NULL
for (i in 1:4) {
  s0 <- if (i == 1) start else start + rnorm(2, 0, 0.4)
  fit <- optim(s0, cost, method = 'Nelder-Mead',
               control = list(maxit = 400, reltol = 1e-10))
  cat(sprintf('  restart %d: cost %.6g  (ka=%.4f Tlag=%.3f h)\n',
              i, fit$value, exp(fit$par[1]), tlag_from(fit$par[2])))
  if (is.null(best) || fit$value < best$value) best <- fit
}

ka <- exp(best$par[1]); Tlag <- tlag_from(best$par[2])
# Final evaluation at production settings (dt = 0.1, 28 doses)
m_prod <- sim_shape(ka, Tlag, dt = 0.1, n_doses = 28)
sc <- scale_to_F(m_prod)
est <- c(F_oral = sc$F_oral, ka = ka, Tlag = Tlag)

# Verify steady state really is attained at the production dose count
m56 <- sim_shape(ka, Tlag, dt = 0.1, n_doses = 56)
ss_err <- abs(m_prod['AUCtau_total'] / m56['AUCtau_total'] - 1)
cat(sprintf('\nSteady-state check, 28 vs 56 doses: relative AUC difference %.3g\n', ss_err))
if (ss_err > 1e-5) stop('Steady state not attained at 28 doses.')

cat('\nCalibrated parameters:\n'); print(round(est, 5))
cmp <- data.frame(metric = names(target), observed = as.numeric(target),
                  fitted = as.numeric(sc$moments[names(target)]))
cmp$rel_error_pct <- 100 * (cmp$fitted / cmp$observed - 1)
cat('\nFitted vs observed:\n'); print(cmp, row.names = FALSE, digits = 4)
cat(sprintf('\nImplied CL/F = %.3f L/h (label systemic CL = 3.30 L/h, label F ~ 0.50)\n',
            3.3 / est['F_oral']))

src <- c(
  sprintf('Calibrated to the Frost 2013 5 mg BID day-7 AUCtau on %s, with CL held at the label 3.3 L/h. Oral data identify only CL/F; this value exceeds the label ~0.50 and the gap is an unresolved conflict between the label IV clearance and the observed oral exposure.', Sys.Date()),
  'Calibrated absorption rate constant; replaces the inherited 0.7 that gave Tmax 1.4 h against an observed median of 4 h.',
  'Calibrated absorption lag. Only weakly separated from ka by aggregate moments, and ridge-penalised towards zero; a fitted descriptor of the observed delay to peak, not an identified transit time.')
write.csv(data.frame(parameter = names(est), value = unname(est),
                     units = c('fraction','1/h','h'), source = src),
          'data/pk_calibration.csv', row.names = FALSE)
dir.create('output/revised', recursive = TRUE, showWarnings = FALSE)
write.csv(cmp, 'output/revised/pk_calibration_fit.csv', row.names = FALSE)
cat('\nPK CALIBRATION COMPLETE\n')
