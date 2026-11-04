# Attempt S-9 findings: the survival-weighted YFB fit at 2,000 genes

Written 2026-11-04. Registration:
[`2026-11-03-survival-weight-p2000-plan.md`](../plans/2026-11-03-survival-weight-p2000-plan.md).
Outputs: `results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid/` (`gate_cells.csv`,
`gate.csv`, `summary.csv`, `cv_chosen_cells.csv`, `weight_curves.png`) and
`results/benchmark_sim/outputs/supervision_diagnostics/s9_rep_n400_b3.3/`. All 480 registered
fits ran without error. No real or external cohort was used.

## Summary

- **Gate S passes at w_S = 1.** With the z-scored input it was written for,
  `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` recovers the planted program in 10 of 10 datasets
  in every p = 2,000 cell, against 0, 0 and 1 of 10 for the S-7 unsupervised reference. Mean
  false-active is 0.9, 0.1 and 0.0, and all 30 fits converge. All three cells qualify. The
  verdict is the same with the canonical w_S = 0 fit as the reference.
- **Training CV picks this weight.** In 55 of 60 datasets the arm with the highest training CV C
  is w_S = 1. The CV-chosen fits at p = 2,000 recover 10, 9 and 10 of 10 with mean false-active
  0.9, 0.1 and 0.0. So the weight can be chosen without the planted truth.
- **Larger weights are worse.** At p = 2,000 every arm with w_S ≥ 40 stops at 300 sweeps without
  converging and loses the program. At p = 50 every weight recovers it, but false activity rises
  steadily with w_S.
- **Remaining problem: coefficient calibration at p = 2,000.** At w_S = 1 the validation
  calibration slope is 1.20, 0.42 and 0.18 (median, by planted variance), and no 95% interval
  for the matched coefficient contains β* in any p = 2,000 fit. The ranking of patients is good
  (validation C 0.94 to 0.99), but the size of the predicted risk differences is too large.
  This is a tempering question, not one that downweighting genomics addresses (see below).
- **Replication passes.** The October 15 bed gives 18 of 20 recovered at noise 0.25 and 19 of 20
  at noise 1, with the same recovered seeds as the October 30 run. Per-seed differences in
  sweeps and correlations are explained by an ELBO fix committed after that run.

## Gate S at p = 2,000

Recovery is out of 10 datasets; reference is the S-7 unsupervised fit on the same datasets.

| Arm | Variance 0.25% | Variance 1% | Variance 5% | Qualifying cells |
|---|---|---|---|---|
| w_S = 1 | 10 (ref 0), FA 0.9, conv 10 | 10 (ref 0), FA 0.1, conv 10 | 10 (ref 1), FA 0.0, conv 10 | **3** |
| w_S = 10 | 10, FA 2.2, conv 10 | 10, FA 1.8, conv 10 | 10, FA 0.8, conv 10 | 1 |
| w_S = 40 (p/50) | 3, FA 2.6, conv 0 | 4, FA 1.5, conv 0 | 10, FA 1.3, conv 0 | 0 |
| w_S = 100 | 0, FA 2.2, conv 0 | 0, FA 1.0, conv 0 | 8, FA 1.8, conv 0 | 0 |
| w_S = 1,000 | 1, FA 1.5, conv 0 | 0, FA 1.0, conv 0 | 0, FA 1.0, conv 0 | 0 |
| w_S = 10,000 | 1, FA 3.9, conv 0 | 1, FA 1.9, conv 0 | 2, FA 1.1, conv 0 | 0 |
| share 0.2 | 2, FA 0.0, conv 10 | 2, FA 0.0, conv 10 | 6, FA 0.0, conv 10 | 1 |

FA is mean false-active programs; conv is the number of converged fits. The gate reads the best
registered arm, which uses the planted truth; here CV picks the same arm (next section).

**Effect of centering.** S-7 ran the same fitter at w_S = 1 on the uncentered nonnegative input
and recovered 10 of 10 in each cell, but with mean false-active 1.0, 1.5 and 1.1, so it would have
failed the gate. With centered input, false activity falls to 0.9, 0.1 and 0.0. The uncentered
fit had to spend a program on the column means. That program could pick up survival signal.

## Training cross-validation as the weight selector

Each dataset's arm with the highest five-fold training CV C (w_S = 0 excluded):

| p | Variance | Arms chosen | Recovered | Mean FA | Converged |
|---|---|---|---|---|---|
| 50 | 0.25% | w1 ×10 | 10 | 1.3 | 9 |
| 50 | 1% | w1 ×10 | 10 | 0.5 | 10 |
| 50 | 5% | w1 ×10 | 10 | 0.3 | 10 |
| 2,000 | 0.25% | w1 ×10 | 10 | 0.9 | 10 |
| 2,000 | 1% | w1 ×8, share0.2 ×2 | 9 | 0.1 | 10 |
| 2,000 | 5% | w1 ×7, w10 ×2, share0.2 ×1 | 10 | 0.0 | 10 |

## Per-weight curves

`weight_curves.png` plots recovery, false activity, survival share, calibration slope, β*
coverage, convergence, validation C and CV C against w_S for both p.

**p = 50.** All weights from 1 to 10,000 recover 10 of 10 in every cell and converge, except one
w_S = 1 fit. Mean false-active rises with w_S: at 0.25% variance it is 1.3, 2.2, 3.2, 4.4 and
6.0 for w_S = 1, 10, 100, 1,000 and 10,000. The calibration slope is 1.04 to 1.07 at w_S = 1 and
about 0.85 at every larger weight. Validation C changes little (0.89 to 0.92).

**p = 2,000.** w_S = 1 and 10 recover everything and converge. From w_S = 40 upward no fit
converges within 300 sweeps, recovery collapses, and validation C falls (at 1% variance: 0.972 at
w_S = 1, 0.904 at 40, 0.815 at 10,000).

## Survival share from p = 50 to p = 2,000

The final in-fit median share of the matched program's F precision that comes from the Cox term,
at w_S = 1:

| Variance | p = 50 | p = 2,000 |
|---|---|---|
| 0.25% | 0.998 | 0.202 |
| 1% | 0.909 | 0.049 |
| 5% | 0.789 | 0.010 |

As expected, the share falls by a factor of 5 to 80 as p grows at a fixed weight. The fall does
not cost recovery: w_S = 1 recovers 10 of 10 at p = 2,000 even at a share of 0.01. Raising the
weight to restore the p = 50 share (w_S = 40 gives shares of 0.45, 0.06 and 0.74) breaks convergence instead.
So a p-scaled weight is not needed, and w_S = p/50 is harmful.

**share0.2 arm.** The calibrated weight came out between 0.002 and 0.02, close to unsupervised.
The w_S = 0 fit, which the registration uses to calibrate the weight, has coefficients of
essentially zero (max |E β| = 0 to 3×10⁻⁷). This is by construction: with `canonicalize = TRUE`,
w_S also weights the Cox term of the β update, so w_S = 0 removes survival from β as well. The
registration's description of this arm ("β is still fitted") was wrong. The survival term in the precision is therefore
near zero at that state, and the solver gives a weight that does not represent the supervised
fit. The arm is reported as registered. It recovers 2, 2 and 6 of 10 with no false activity.
A share-targeted weight would need to be calibrated at a supervised state.

## Non-convergence at large w_S

At p = 2,000, arms with w_S ≥ 40 reach the 300-sweep cap with max |E β| between 17 and 26
(against 5 to 8 at w_S = 1). Reconstruction stays bounded: relative loss is at most 0.985
against 0.979 for w_S = 0. The minimum noise precision τ_j is at least 0.92. No lost-variance
records occurred. This is the β mechanism found in the earlier attempts (overconfident
coefficients with bounded reconstruction), not the τ mechanism (falling τ_j with
reconstruction loss above 1). At p = 50 the same weights converge but carry more false-active
programs, which is the same mechanism in a milder form.

## Coefficient calibration at w_S = 1

| p | Variance | Calibration slope (median) | β* covered | Max \|E β\| | Validation C |
|---|---|---|---|---|---|
| 50 | 0.25% | 1.04 | 1/10 | 3.8 | 0.902 |
| 50 | 1% | 1.07 | 2/10 | 4.1 | 0.908 |
| 50 | 5% | 1.07 | 3/10 | 5.5 | 0.920 |
| 2,000 | 0.25% | 1.20 | 0/10 | 5.4 | 0.944 |
| 2,000 | 1% | 0.42 | 0/10 | 7.3 | 0.972 |
| 2,000 | 5% | 0.18 | 0/10 | 8.3 | 0.986 |

A calibration slope below 1 means the validation risk differences are smaller than the fitted
coefficients predict. At p = 2,000 and the two larger variances the fitted coefficient is about
2 to 5 times too large for out-of-sample risk. The intervals are also too narrow, since none
covers β*, and coverage is low even at p = 50. Discrimination is not affected.

**Why S-10 is not run.** The work order runs S-10 if S-9 fails Gate S, or if it passes but leaves
a problem S-10 is derived to fix. S-10 downweights the genomics term, which raises survival's
share of the F precision. In this grid, raising that share (by raising w_S) made calibration and
false activity worse, not better. The remaining problem is the scale of the Cox likelihood
relative to the posterior for β, which is what a tempering calibration addresses.

**Follow-up, not implemented here.** Calibrate the Cox-likelihood weight by SafeBayes
(Grünwald) or by the Holmes–Walker power-likelihood calibration (2017), both of which choose the
power on the likelihood from the data. A prior on w_S does not
substitute for this: it averages over weights, it does not calibrate the likelihood's scale.

## Replication cell

`run_supervision_diagnostics.R s9_rep 20 2 400 3.3` (p = 50, n = 400, β₃ = 3.3, w_S = 1, no
supervised start):

| Noise | Recovered | Converged | October 30 (A1) recovered | October 30 converged |
|---|---|---|---|---|
| 0.25 | 18/20 | 20/20 | 18/20 | 19/20 |
| 1 | 19/20 | 20/20 | 19/20 | 19/20 |

The criterion (at least 15 of 20 at noise 0.25) is met, and the recovered seeds are the same.
Sweep counts and correlations differ in 17 of 40 fits (largest |Δr| 0.034), and two fits that
hit the 300-sweep cap on October 30 now converge. A1 ran at 04:07 on October 30. The ELBO fix
that counts the Gaussian log-normalizer once per observation in `compute_ebnm_kl` was committed at
04:47 that day, and it changes where the relative-ELBO stopping rule fires. Refitting three seeds
with the pre-fix function restored reproduces the October 30 sweeps and correlations exactly
(300 sweeps not converged, r = 0.9312; 300 not converged, r = 0.8951; 71 sweeps, r = 0.8861).
The R and compiled sweep engines agree to the same stopping sweep in the seeds checked, so the
engine is not a source of the difference.

## Decision and next steps

S-9 passes Gate S at w_S = 1, the fit works in simulation, and training CV selects the same weight.
Under the 2026-11-04 decision the real-data steps follow, each registered before it runs:

1. Step S3 of the realignment plan: training CV on the single-log D4 data with
   `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` on z-scored input, reporting the identified
   programs, which are survival-active, and training CV C.
2. The feature-set-size curve (2,000 genes to all genes) and the survival-screening comparison,
   with screening inside folds.
3. The null-data calibration check at w_S = 1.
4. One external-cohort scoring of the fixed model, after training CV has fixed every choice.

The calibration problem above is carried as a known limitation. Real-data C is unaffected by it,
but coefficients and hazard ratios from these fits should not be read as calibrated.
