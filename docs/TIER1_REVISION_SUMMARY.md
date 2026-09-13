# Tier-1 revision: summary

Trimmed summary of the Tier-1 revision that produced the current model. It
changed three things and, unlike a bugfix pass, changed what the analysis
concludes: earlier results are superseded, not merely corrected. Two of the
three changes are calibrations to published aggregate summaries, so
agreement with those summaries is now a fit rather than a test. Nothing here
converts the model into a validated dosing tool.

## 1. The efficacy endpoint was degenerate

The dose-optimization score originally scored efficacy as ETP (endogenous
thrombin potential) suppression, `1 - ETP(dose)/ETP(0)`. At reference trough
exposure that quantity is **1.1e-4**, so the composite score inherited this
near-zero magnitude (max 8.99e-05), and the surface was flat enough that the
dose search returned the 20-mg search bound instead of an interior maximum
in the moderate and severe scenarios. In the synthetic population, 83.5% of
severe-NS maximizers sat at that boundary (87.0% in the tissue-factor
propagation). Those numbers described the objective, not the pharmacology.

The same simulation already contained a drug-sensitive summary: peak
thrombin fell 38.7% at the same exposure (time to peak nearly doubled),
against ETP's 1.1e-4. **The efficacy term is now peak-thrombin suppression.**
The ETP formulation is still computed and reported beside every result
(`efficacy_score_etp`, `TI_etp`, `score_etp`) so the change stays auditable.

This is a change of endpoint, not a discovery that the drug works better
than thought: peak thrombin and ETP answer different questions. After
recalibration the contrast sharpened, not softened: severe-NS ETP
suppression is 2.48%, 0.028%, and 0.025% over 600/1200/2400-second assay
windows (versus 7.22%, 1.88%, 2.24% before), while peak-thrombin suppression
is 42.47% and does not vary with the window at all, because a peak is a
point statistic rather than an integral. That contrast is the single
clearest argument for the endpoint change.

**Effect on dose optimization:**

| Quantity | Before | After |
| --- | --- | --- |
| Reference feedforward maximizer | 1.71 mg, score 8.99e-05 | 3.62 mg, score 0.278 |
| Severe-NS feedforward maximizer | 20.00 mg (at search bound) | 4.55 mg |
| Severe-NS boundary rate, population | 83.5% | 0.0% |
| Severe-NS boundary rate, TF draws | 87.0% | 0.0% |

Every scenario and feedback mode now has an interior optimum; boundary
solutions and weak-objective flags are zero across all 1600 paired
synthetic records and all 800 tissue-factor draws. The ETP-based score,
recomputed at the same maximizers, still sits between 6e-05 and 2e-04, the
direct evidence that the boundary solutions described the objective, not
the pharmacology. These remain exploratory diagnostics: the safety factor
is an unanchored platelet integral, and maximizers must not be read as
doses.

## 2. The pharmacokinetic model could not reproduce the observed profile shape

The inherited model underpredicted trough by 42.4% and interval AUC by
28.3%, and put Tmax at 1.4 h against an observed median of 4 h. Two fixes:
an absorption lag was added, and bioavailability (F), absorption rate (ka),
and lag (Tlag) were calibrated to the four Frost 2013 aggregate
steady-state moments, which is moment matching, not profile fitting, since no
digitized concentration-time curve was available. Systemic clearance was
**not** fitted and stays at the label values (3.3 and 0.9 L/h); because oral
data identify only CL/F, the fitted F (~0.70) exceeds the label's ~0.50,
a reported conflict rather than a hidden one.

A deep peripheral compartment was also tried and rejected: it contributed
0.39 L of a 28.5 L steady-state distribution volume, changed the fitted
moments by at most 0.2%, and split the terminal phase into two
near-degenerate modes that made terminal half-life window-dependent
(8.78-9.34 h across regression windows, versus a clean 8.51 h without it).
A structure that explains none of the data it was fitted to, while making a
previously well-defined quantity window-dependent, is an artifact, so it
was removed. The delivered model is the five-compartment one with a
calibrated lag only.

## 3. The nephrotic-syndrome coagulation gradient contradicted the one observation

Applied at face value, the literature-directed factor multipliers produced
a severe-NS to reference ETP ratio of 2.92 against an observed 1.20. A
single severity scale `s` now shrinks every multiplier toward one as
`m' = 1 + s(m-1)`, calibrated to that one ratio (fitted `s = 0.238`).
One aggregate ratio identifies exactly one scalar; the calibration cannot
distinguish an over-aggressive factor table from an over-sensitive cascade,
and it does not test the assumed relative weighting between species. The
shrunk multipliers are the gradient this model needs to match one
observation; they are not measured factor levels.

## What's fitted vs. assumed

| Quantity | Status | Fitted to |
| --- | --- | --- |
| F_oral, ka, Tlag | Calibrated | Frost 2013 four aggregate moments |
| NS severity scale | Calibrated | Kelddal 2025 NS:healthy ETP ratio |
| CL_total, CL_renal | Held at label values | not applicable |
| Species weighting, platelet/feedback/age/binding gradients | Assumed | not applicable |

## What stayed independent

Of the eight external-comparison rows, the four Frost moments and the
Kelddal NS:healthy ETP ratio were calibration targets. The two Kelddal
trough concentrations and the absolute ETP levels were **not** fitted and
remain out-of-sample checks: the model, calibrated on Frost 2013 alone,
predicted a different study's steady-state troughs inside their reported
confidence intervals (previously ~35% low in both cohorts). That's two
aggregate means, not a validation.

The central mechanistic result is unaffected by all three changes: total
interval AUC still falls 31.4% from reference to severe NS while unbound
AUC falls only 0.52%, the dissociation between total and unbound exposure
that motivates the model in the first place.
