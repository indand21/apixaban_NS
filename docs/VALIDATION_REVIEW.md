# Validation review

## Overall assessment

Share with caveats as an exploratory numerical modeling study. Not qualified for clinical dosing, thrombosis prediction, bleeding prediction, or regulatory decision-making.

Updated after the Tier-1 revision (docs/TIER1_REVISION_2026-09-09.md). Absorption and the NS coagulation gradient are now calibrated, so five of the eight external comparison rows are fitted and agreement with them is not evidence of validation. The status of each row is given below and in the manuscript.

## Methodology and corrections

The review covers the project source, parameter and modifier tables, numerical tests, reference data, analysis scripts, inherited outputs, manuscript, supplement, and figure provenance. The original manuscript described CKD rather than the current NS project and used stale numerical results. Those assets are archived and superseded.

The revised primary comparison controls age, weight, renal function, and assay TF across disease scenarios. Rate loaders use authoritative CSV values with unit conversion. Drug effects include formed prothrombinase; implicit protein C activation is disabled in the no-TM primary assay. All synthetic individual effects propagate to simulations, with matched individual and feedback-mode no-drug baselines. Prior propagation and descriptive percentile ranges are not called posterior estimates or credible intervals.

The diagnostic dose search includes zero and both endpoints, brackets local maxima, skips exactly flat local plateaus, and reports boundary solutions. It does not assume an interior maximum or turn a weak or boundary-constrained score into a clinical recommendation. Stroke mapping is disabled; the separate toy HC hazard has a corrected sign and an explicit warning but is excluded from the analysis.

## Calculation checks

The final test suite records its exact status in output/revised/test_suite_status.csv. Revision-specific checks cover R/C equation agreement, all-state trajectories, ten conserved coagulation moieties, platelet balance with the desensitized sink, PK mass balance, dose scaling, steady state, unit conversion, complete propagation of sampled effects, lognormal mean/CV specification, flat objective behavior, and CSV field counts. Source and result checks are written to output/revised/revision_tests.csv and result_integrity.csv. Numerical verification is not clinical validation.

The low-extraction PK approximation was distinguished from exact finite-flow behavior. A previous half-life invariance assertion was replaced by comparison with the linear-system eigenvalue. All reported ETP values are converted from nM.s to nM.min. HC retains model-equivalent units. Manuscript claims are assembled from full-precision outputs and recorded separately in manuscript_claims.json.

## Presentation review

Figures use explicit units and scenario definitions. Feedback and feedforward series have distinct line styles. Log score axes omit zeros explicitly, and sensitivity heatmaps label numeric percentages. No dose heatmap is presented as a prescribing nomogram. External observations are shown in an exact table with source-specific metrics and cohort definitions rather than a pooled pass/fail statistic. All 22 final Word pages were rendered with Microsoft Word and visually inspected. Tables, repeated table headers, native equations, figure captions, and page boundaries passed review. The final abstract contains 176 words. Automated checks found all eight native tables, five embedded figures, and no unresolved placeholders or out-of-bounds text. All 119 numerical/component checks and nine result-integrity checks passed.

## Required scientific caveats

The inherited cascade and platelet kinetics remain incompletely source-qualified. The model does not reproduce the original Hockin topology or purified prothrombinase substrate kinetics. Platelet and surface normalizations are arbitrary. NS severity gradients, age effects, tissue binding scaling, and feedback parameters are provisional. Aggregate clinical comparison uses an imperfect cohort proxy and incompletely matched assay conditions. Independent synthetic distributions do not reproduce a measured NS population, and no clinical outcome observation model exists.

External agreement after Tier-1 calibration is as follows. The four Frost 2013 moments were fitted: peak +2.6 percent, trough +7.6 percent, interval AUC 0.0 percent by construction, and time to peak 3.30 h against an observed median of 4 h with a reported range of 2 to 4 h. The Kelddal 2025 NS:healthy ETP ratio was fitted and is reproduced at 1.204.

Two rows were not fitted and are therefore genuine out-of-sample checks of a model calibrated on a different study. Predicted steady-state troughs are 60.6 ng/mL against an observed 51 (95 percent CI 39 to 64) in healthy participants and 41.2 against an observed 35 (CI 28 to 43) in nephrotic syndrome. Both fall inside the reported intervals, having previously been 35 and 37 percent low. This is the strongest external evidence the project holds, and it is still only two aggregate means.

One systematic discrepancy remains. Absolute ETP is overpredicted by 42.9 percent in both cohorts, outside both reported intervals. Because the calibration fixed only the ratio, what survives is a single common scale factor rather than a structural gradient error, and the unverified tissue-factor reagent concentration of the comparator assay is the obvious candidate. It has not been resolved and must not be described as agreement.

These limitations do not prevent reporting the exploratory study accurately. They do prevent claiming a clinically correct dose, validated therapeutic index, or patient-level predictive performance. Appropriate next work requires measured assay-matched inputs and independent clinical qualification, not simply more synthetic simulations.

## Handoff details

The earlier manuscripts, tables, figures, and cover letter remain recoverable in archive. No external files, clinical protocols outside this project, patient records, submissions, messages, or publications were changed. The cloud-drive interruption was handled using a local temporary working copy and input-hashed per-case checkpoints. Author declarations and full submission metadata were not invented and must be confirmed by the authors before journal submission.
