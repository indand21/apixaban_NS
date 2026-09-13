# Project revision audit

SUPERSEDED IN PART. This document records the revision that made the project numerically
correct and reproducible. It deliberately left three scientific problems in place: a
degenerate ETP-based objective, an uncalibrated absorption model, and an NS coagulation
gradient contradicted by the observed thrombin-potential ratio. Those were addressed
afterwards; see docs/TIER1_REVISION_2026-09-09.md, which takes precedence where the two
disagree. In particular, the statements below that the model underpredicts healthy trough
exposure and overpredicts severe-NS ETP describe the pre-Tier-1 model.

The September 2026 revision replaces the inherited CKD manuscript with a reproducible nephrotic-syndrome analysis. Clinical dose recommendations require evidence beyond this exploratory model.

## Scope

Audit all R modules, disease and kinetic inputs, tests, analysis scripts, outputs, and manuscript components. Preserve the pre-revision files in an archive before replacing generated results. Record corrections and distinguish numerical verification from clinical qualification.

## Confirmed issues

- Systemic clearance 3.3 L/h and renal clearance 0.9 L/h were incorrectly multiplied by oral bioavailability.
- The manuscript and outputs predate September binding/distribution corrections.
- The inherited coagulation model is a reduced, modified cascade, not a faithful Hockin-Mann reproduction; protein-C extension and unused parameters require explicit documentation.
- Formed prothrombinase was not directly inhibited by apixaban.
- The rate/species CSV files were read but ignored by parameter loaders.
- Sequential platelet simulations omitted age effects; virtual-patient code sampled but omitted volume, absorption and platelet effects.
- Virtual patients used mismatched population baselines, including uncoupled baselines for coupled simulations.
- The hypothetical bleeding-hazard equation had an inverted sign; the outcome transfer functions are not qualified for NS or stroke prediction.
- Baseline cache keys omitted simulation duration; time-varying forcing started at the first dose instead of the steady-state interval.
- Prior propagation was mislabeled as posterior inference, and surrogate optima were misrepresented as clinical recommendations.

## Implemented corrections

The active cascade and platelet feedback now share component equations. The compiled implementation is checked against the R equations. Numerical investigation identified tolerance-sensitive platelet amplification; primary coagulation-platelet integration now uses relative/absolute tolerances 1e-12/1e-14. Tightened endpoint comparisons pass a 0.01% relative-error criterion, and all ten coagulation moieties conserve to numerical precision.

Systemic reference clearance was restored to 3.3 L/h with a 0.9-L/h renal component. Disease nonrenal intrinsic-clearance reduction is a sensitivity hypothesis rather than an imposed calibrated fact. The primary bound-loss coefficients are 0/0.005/0.01/0.02 L/h. Fixed TF, the no-implicit-TM primary assay, and explicit inhibition of formed prothrombinase replace the prior inconsistent assay assumptions. Exact purified-enzyme agreement is not claimed for the inherited cascade rates.

Rate tables are authoritative and use explicit unit conversion. Unused k17 and k42 are removed; k39 is disabled in the primary assay and retained as a structural switch. Platelet and HC molar labels were corrected to model-equivalent units. The desensitized sink is included in balance checks. Feedback coefficients are exploratory assumptions, not fitted dose bounds.

Every sampled volume, absorption, binding, coagulation, and platelet effect is propagated. Lognormal draws use the correct mean-one CV parameterization. Identical draws are paired across disease stages and feedback modes. Baselines belong to the same individual, feedback mode, and assay duration. The search includes zero and both endpoints and avoids redundant refinement of exactly flat plateaus. Boundary frequencies and weak-objective flags are preserved.

Clinical stroke mapping is disabled. The separate hypothetical HC hazard has a corrected direction and an explicit warning, and is excluded from all revised results. No clinical event-rate model, therapeutic-index validation, or prescribing nomogram is reported. TF uncertainty is labeled prior propagation without Bayesian updating.

## Regenerated analysis

The canonical scripts regenerate 32 reference scenario/dose/mode rows, eight descriptive external comparisons, 216 structural-sensitivity rows, 88 local-sensitivity rows, 1600 records from 200 paired synthetic individuals, a 100-point patient-factor grid, and 800 TF-prior records. Full precision outputs, result-integrity checks, checksums and session information are in output/revised. Manuscript numerical claims are assembled from those files rather than copied from the historical text.

All 119 executed checks passed, comprising 47 revision-specific checks and 72 component checks across seven files. Exact executed pass/fail status is recorded in output/revised/test_suite_status.csv and the logs. All nine final result-integrity checks also passed. The finalizer refuses to mark analysis complete when tests or expected row/key checks fail.

## Source qualification and remaining limits

Primary source checks are documented in docs/PRIMARY_SOURCE_NOTES.md. The NS comparator has internally inconsistent unit typography and an abstract-versus-results CI/P-value discrepancy; the full results interval is retained, with the discrepancy disclosed. No external observation was invented to fill a calibration gap.

The model underpredicts healthy apixaban trough concentrations and overpredicts severe-scenario ETP relative to the NS aggregate comparator. It remains an exploratory reduced model with provisional kinetics, disease gradients, platelet normalization, age rules, and feedback. The manuscript explicitly explains these failed comparisons and does not convert numerical verification into a clinical-validation claim.

## Reproducibility and archival actions

An interruption of the Google Drive mount terminated uncheckpointed workers. The remaining analysis was restarted in a local temporary working copy. Per-case CSV checkpoints are written atomically and separated by input hashes, preventing a subsequent interruption from discarding completed individual calculations. The compiled solver no longer repeatedly probes the cloud filesystem once loaded.

Historical manuscripts, supplements, cover letters, tables, and figures were moved recoverably to archive/retired_legacy_assets_20260909; an earlier complete snapshot is in archive/pre_revision_20260909. Active legacy analysis scripts now stop with a pointer to the revised pipeline. README_NS.md, CLI_CONTEXT.md, and scientific notes no longer repeat stale dose recommendations.

The final manuscript and supplement are generated in the manuscript package folder with editable Markdown sources and five rebuilt figures. Document generation used the verified installed Pandoc/Python/Word toolchain because the managed document runtime was unavailable. The canonical renderer was attempted and could not find LibreOffice; installed Microsoft Word exported both documents successfully. Every final page was rendered and visually inspected: 13 manuscript pages and nine supplementary pages. Native table counts were two and six, embedded figure counts four and one, and native equation counts 19 and five. No unresolved placeholders or out-of-bounds text were detected. Complete authorship and author declarations must be verified before journal submission; missing personal or funding details were not fabricated.
