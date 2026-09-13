# Apixaban NS QSP

Revised 9 September 2026 (Tier-1). This exploratory model links apixaban pharmacokinetics, a reduced coagulation cascade, and a phenomenological platelet module. Numerical verification does not establish clinical qualification. The model does not support clinical dose selection. Absorption and the nephrotic-syndrome coagulation gradient are now calibrated to published aggregate summaries, so agreement with those summaries is a fit rather than a test.

## Repository contents

This repository contains the model code, calibration inputs, tests, and
technical documentation only. Running the pipeline (below) regenerates
results locally under `output/revised` (full-precision tables, a source
manifest, and software-version information) and `output/_run_logs`
(verification logs); neither is checked in, since both are fully
reproducible from the code and data here.

- Calibration inputs: `data/pk_calibration.csv` and
  `data/ns_severity_calibration.csv`; refitting them regenerates fit
  summaries under `output/revised/pk_calibration_fit.csv` and
  `ns_severity_fit.csv`.

## Reproduction

Run from the project root:

```powershell
Rscript --vanilla scripts/run_ns_pipeline.R
```

The pipeline compiles the C solver if necessary, runs all tests, and regenerates results and figures under output/revised. To skip recompiling and rerunning the analysis and only rebuild figures/tables from existing results, add --figures-only. NS_N_PATIENTS sets synthetic individuals per scenario (default 200); NS_WORKERS sets parallel workers (default 6). TF propagation always uses 200 paired prior draws. Individual cases are atomically checkpointed under output/revised/checkpoints with input-hash namespaces.

Calibration is deliberately outside the pipeline. The pipeline consumes data/pk_calibration.csv and data/ns_severity_calibration.csv as authoritative inputs. To refit, run scripts/calibrate_pk_absorption_distribution.R and then scripts/calibrate_ns_severity.R, in that order, because the severity fit depends on the exposure the pharmacokinetic fit produces. The test suite fails if either table stops reproducing its target, so a stale table cannot pass silently. Do not edit a script while it is running; the interpreter rereads the file and the run aborts.

R 4.5 with compatible Rtools and packages listed in R/00_packages.R are required. No packages are silently installed.

For long computations, copy R, data, src, scripts, tests, and output into a local temporary working directory, run there, and copy verified outputs back. Avoid keeping long-running output streams on a streamed cloud drive.

## Scientific interpretation

Primary scenarios hold eGFR at 100 mL/min/1.73 m2, age at 40 years, and weight at 70 kg. Labels Normal, NS_Mild, NS_Moderate, and NS_Severe denote provisional model perturbations, not validated severity classes. Historical ckd_stage names are compatibility interfaces. Assay TF is fixed at 5 pM; implicit thrombomodulin-dependent protein C activation is disabled in the primary assay.

The score combines PEAK THROMBIN suppression with matched platelet aggregate retention. It previously used ETP suppression, which is about 1e-4 at therapeutic trough exposure and left the objective flat enough that the dose search returned its 20-mg search bound instead of an interior maximum. The ETP formulation is still computed and reported alongside as efficacy_score_etp, TI_etp and score_etp. Historical fields TI, efficacy_score, safety_score, and optimal_dose are compatibility labels, not clinical endpoints or recommended doses. Report boundary maxima and weak objectives. No prescribing nomogram, clinical stroke prediction, or validated bleeding probability is produced.

Of the eight external comparison rows, the four Frost 2013 moments and the Kelddal 2025 NS:healthy ETP ratio were used for calibration; the two Kelddal trough concentrations and the absolute ETP levels were not, and remain out-of-sample checks of a model calibrated on a different study. Core rate constants, platelet scaling, the species weighting inside the NS factor table, and the remaining disease modifiers require further experimental qualification.

## Scope of this repository

This repository ships the model, its calibration inputs, its test suite, and
its technical documentation. It does not include manuscript sources,
generated results, or superseded/archived material. Legacy analysis scripts
that were superseded during development stop explicitly (see their `stop()`
guard) rather than silently regenerating outdated outputs.
