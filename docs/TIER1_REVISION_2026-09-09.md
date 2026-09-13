# Tier-1 revision audit

This revision follows docs/REVISION_AUDIT_2026-09-09.md. That earlier pass made the
project numerically correct and reproducible; it deliberately left three scientific
problems in place. This pass addresses those three. It changes what the analysis
concludes, so the earlier results, preserved in archive/pre_tier1_20260909, are
superseded rather than corrected.

Nothing here converts the model into a validated dosing tool. Two of the three changes
are calibrations to published aggregate summaries, which means agreement with those
summaries is now a fit rather than a test.

## 1. The efficacy endpoint was degenerate

compute_therapeutic_index scored efficacy as ETP suppression, 1 - ETP(D)/ETP(0). At
reference trough exposure that quantity is 1.1e-4. The composite score is the product of
that term and platelet retention, so the score inherited the near-zero magnitude, its
maximum was 8.99e-05, and the surface was flat enough that the dose search returned the
20-mg search bound instead of an interior maximum in the moderate and severe scenarios.
In the synthetic population 83.5 percent of severe-NS maximizers sat at a boundary, and
87.0 percent in the tissue-factor propagation. Those numbers described the objective,
not the pharmacology.

The same simulation already produced a drug-sensitive summary: in that same
pre-calibration model peak thrombin fell 38.7 percent at the same exposure and time to
peak nearly doubled, against the 1.1e-4 of the ETP term. The efficacy term is
now peak-thrombin suppression. The ETP formulation is computed and reported beside every
result as efficacy_score_etp, TI_etp and score_etp, so the change is auditable and the
previous behaviour stays inspectable.

This is a change of endpoint, not a discovery that the drug works better than thought.
Peak thrombin and ETP answer different questions, and the structural-sensitivity result is
still reported. Its values did change with the recalibration and are now more extreme, not
less: severe-NS ETP suppression is 2.48, 0.028 and 0.025 percent over 600, 1200 and 2400
seconds, against 7.22, 1.88 and 2.24 percent before. Peak-thrombin suppression over the
same windows is 42.47 percent and does not vary with the window at all, because a peak is a
point statistic rather than an integral. That contrast is the clearest single argument for
the endpoint change.

## 2. The pharmacokinetic model could not reproduce the observed profile shape

The inherited model matched steady-state peak concentration to within 6 percent but
underpredicted the trough by 42.4 percent, underpredicted interval AUC by 28.3 percent,
and put Tmax at 1.4 h against an observed median of 4 h. The peak-to-trough ratio was
4.2 against an observed 2.6.

Two separate defects were involved.

The AUC error was not a structural problem and could not be fixed by one. At steady
state AUCtau = F x Dose / CL exactly, independent of absorption and distribution.
Holding CL at the label systemic value of 3.3 L/h, the observed AUC forces
CL/F = 4.74 L/h and therefore F = 0.697, well above the label's approximately 0.50.
Oral data identify only CL/F, never F and CL separately. The earlier audit established
that CL must not be multiplied by F, and that correction is preserved: CL stays at 3.3
and 0.9 L/h and is not fitted. The conflict is therefore absorbed into F and reported.
The fitted F should not be read as an estimate of absolute bioavailability; it is the
value that reconciles a label intravenous clearance with an observed oral exposure that
disagree.

The shape error was structural. An absorption lag was added, and ka and Tlag were fitted
with F to the four Frost 2013 aggregate moments. Because no digitised concentration-time
profile was available, this is moment matching, not profile fitting, and yields no
standard errors. ka and Tlag are only weakly separated by aggregate moments, so a ridge
penalty keeps the lag as small as the Tmax observation requires. The fitted lag describes
an observed delay; it is not an identified transit time.

### A deep peripheral compartment was tried and rejected

The plan for this item was an absorption lag plus a second, slowly equilibrating
peripheral compartment. The compartment was implemented, calibrated, and then removed,
because the evidence did not support it. That negative result is recorded here rather
than quietly dropped.

An unconstrained first fit ran away to a deep pool with a 380-hour equilibration time
constant. That is a slow sink, not a distribution compartment: the 14-day simulation
never reaches steady state, the final-interval AUC is depressed, and the profiled F
inflates to 0.736 to compensate. Bounding the equilibration time constant to at most 18 h
fixed that, and steady-state attainment was then verified directly (28 versus 56 doses
agreeing to 7.4e-12).

With that bound the calibration drove the compartment's flow to 0.031 L/h and its
partition coefficient to 0.039, contributing 0.39 L of a 28.5 L steady-state distribution
volume. Setting its flow to zero changed peak concentration by 0.05 percent, trough by
0.20 percent, and AUC and Tmax not at all.

It was not merely inert. It split the terminal phase into two near-degenerate modes of
7.95 and 9.37 h, so the apparent terminal half-life depended on the regression window:
8.78 h over 48-72 h, 9.14 h over 120-200 h, 9.34 h over 250-400 h. The model without it
returns a clean 8.51 h in every window. A structure that contributes nothing to the data
it was fitted to, while making a previously well-defined quantity window-dependent, is an
artifact, so it was removed.

The within-interval decay problem is therefore explained by delayed and slower absorption
alone. The delivered model is the five-compartment one with a calibrated lag.

## 3. The nephrotic-syndrome coagulation gradient contradicted the one observation

Applied at face value, the multipliers in data/ns_factor_levels.csv produced a
severe-NS to reference ETP ratio of 2.92 against an observed 1.20 (1096 versus
910 nM.min). The healthy baseline was 43 percent high; the gradient was the larger
error, and antithrombin at 0.60 with factor VIII at 1.70 did most of the work.

A single severity scale s now shrinks every multiplier towards one as
m' = 1 + s(m - 1), calibrated by uniroot to that observed ratio. The fitted value is
s = 0.2380, which reproduces the ratio to 1.2044. Applied severe-NS multipliers become
1.17 for factor VIII and 0.90 for antithrombin.

The limits of this are important and are stated in the manuscript rather than buried.
One aggregate ratio identifies exactly one scalar. The calibration cannot distinguish an
over-aggressive factor table from a cascade whose thrombin integral is over-sensitive to
that table, and it does not test the assumed relative weighting between species, which
simply moves together. The shrunk multipliers are the gradient this model needs in order
to match one observation. They are not measured factor levels and must not be quoted as
such.

One external check was not used in the fit and is worth recording. A 9.5 percent
antithrombin reduction in a scenario labelled severe looks mild against the urinary
antithrombin loss classically described in nephrotic syndrome, which would suggest the
calibration had gone too far. It did not. A multicentre cross-sectional analysis of 47
adults with nephrotic syndrome (Kelddal et al., Kidney360 2025;6(11):1960-1969,
doi:10.34067/KID.0000000865) reports antithrombin at 0.94 and free protein S at 1.13,
both within the normal range and not different from healthy controls, with the
prothrombotic state instead driven by elevated thrombomodulin, syndecan-1 and von
Willebrand factor and by impaired fibrinolysis. The calibrated multiplier of 0.90 sits
close to that measured 0.94 despite having been fitted only to a thrombin-potential
ratio. This is a single independent comparison in a different cohort, not a validation,
but it supports the reading that the original 0.60 assumption was too aggressive and it
identifies the endothelial and fibrinolytic axes as what the model currently omits.

## What is now fitted rather than assumed

| Quantity | Status | Fitted to |
| --- | --- | --- |
| F_oral, ka, Tlag | Calibrated | Frost 2013 four aggregate moments |
| NS severity scale | Calibrated | Kelddal 2025 NS:healthy ETP ratio |
| CL_total, CL_renal | Held at label values, not fitted | not applicable |
| Species weighting in ns_factor_levels.csv | Assumed | not applicable |
| Platelet, feedback, age and binding gradients | Assumed | not applicable |

Refitting is a deliberate, separate step and is not part of the pipeline:

    Rscript --vanilla scripts/calibrate_pk_absorption_distribution.R
    Rscript --vanilla scripts/calibrate_ns_severity.R

The pipeline consumes data/pk_calibration.csv and data/ns_severity_calibration.csv as
authoritative inputs, and the test suite verifies that they still reproduce their
targets. Run the pharmacokinetic calibration first; the severity calibration depends on
the drug exposure it produces.

## Regenerated results

All 132 executed checks pass, comprising 60 revision-specific checks (was 47) and 72
component checks across seven files, and all nine result-integrity checks pass.

### The dose search now returns interior maxima

Every scenario and feedback mode has an interior optimum. Boundary solutions and
weak-objective flags are zero everywhere, in the reference scenarios, in all 1600 paired
synthetic records, and in all 800 tissue-factor draws.

| Quantity | Before | After |
| --- | --- | --- |
| Reference feedforward maximizer | 1.71 mg, score 8.99e-05 | 3.62 mg, score 0.278 |
| Severe-NS feedforward maximizer | 20.00 mg (at bound) | 4.55 mg |
| Severe-NS boundary rate, population | 83.5 percent | 0.0 percent |
| Severe-NS boundary rate, TF draws | 87.0 percent | 0.0 percent |
| Severe-NS TF 2.5-97.5 percentile | 7.21 to 20.00 mg | 0.87 to 8.03 mg |

The ETP-based score, recomputed at the same maximizers, remains between 6e-05 and 2e-04.
That is the direct evidence that the boundary solutions described the objective rather
than the pharmacology.

These are still exploratory diagnostics. The safety factor is an unanchored platelet
integral with no validated mapping to bleeding, and the percentile ranges are propagated
input distributions, not credible intervals. The maximizers must not be read as doses.

### External comparison

| Row | Status | Before | After |
| --- | --- | --- | --- |
| Frost healthy Cmax | fitted | -6.0 percent | +2.6 percent |
| Frost healthy trough | fitted | -42.4 percent | +7.6 percent |
| Frost healthy AUCtau | fitted | -28.3 percent | 0.0 percent |
| Frost healthy Tmax | fitted | 1.4 h | 3.30 h (observed median 4, range 2 to 4) |
| Kelddal healthy trough | HELD OUT | -35.0 percent | +18.8 percent, inside the 39 to 64 CI |
| Kelddal NS trough | HELD OUT | -36.5 percent | +17.7 percent, inside the 28 to 43 CI |
| Kelddal healthy ETP | level not fitted | +42.9 percent | +42.9 percent |
| Kelddal NS ETP | ratio fitted | +248.2 percent | +42.9 percent |

The two held-out rows carry the weight. Absorption was calibrated on Frost 2013 alone,
and the model then predicted a different study's steady-state troughs inside their
reported confidence intervals, having previously been about 35 percent low in both
cohorts. That is two aggregate means, not a validation.

The residual ETP discrepancy is now a single common scale factor of 42.9 percent shared
by both cohorts rather than a structural gradient error. The unverified tissue-factor
reagent concentration of the comparator assay is the obvious candidate and has not been
resolved.

### The central mechanistic result is unchanged

Total interval AUC still falls 31.4 percent from reference to severe NS while unbound AUC
falls 0.52 percent. Absolute exposures rose with the calibrated input fraction, but the
dissociation between total and unbound exposure, which is the manuscript's actual
contribution, is unaffected by all three Tier-1 changes.

## Consequence for validation claims

Before this revision the external comparison was weak but independent. It is now partly
a fit. Of the eight comparison rows, the four Frost moments and the Kelddal ETP ratio
were used for calibration. The two Kelddal trough concentrations were not, and the
absolute ETP levels were not, so those remain out-of-sample checks of a model calibrated
on a different study. The manuscript states which rows are fitted and which are not, and
does not describe agreement with fitted summaries as validation.
