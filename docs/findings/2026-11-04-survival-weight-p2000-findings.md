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
- **Remaining problem: supervision too weak at p = 2,000 (corrected 2026-11-05).** With exact
  event times, the w_S = 1 validation calibration slope is 1.45, 2.16 and 3.23 (median, by
  planted variance), against about 1.0 for the true predictor. Slopes above 1 mean the fitted risk
  differences are too small: the coefficients are attenuated, not overconfident. No 95% interval
  for the matched coefficient contains β* in any p = 2,000 fit. Ranking is good (validation C 0.94
  to 0.99). The first version of this document reported slopes of 1.20, 0.42 and 0.18 and called
  the coefficients overconfident. That was an artifact of `coxph`'s default tie handling, which
  merged the extreme simulated times (see the calibration section).
- **Strength of the test bed.** The planted survival effect is far stronger than on the real data
  (validation C 0.94 to 0.99, SD of the true predictor 9 to 26). Gate S was passed under strong
  survival signal.
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
records occurred, so the τ mechanism (falling τ_j with reconstruction loss above 1) does not
apply. *Corrected 2026-11-05:* the first version called this the β mechanism (overconfident
coefficients). The max |E β| of 17 to 26 is close to the true β* of 9 to 26, so the coefficients
are not inflated relative to the truth. Two features are mechanical: under `canonicalize = TRUE`,
w_S multiplies the Cox precision in the β update, so the posterior SD of β shrinks by 1/√w_S,
and false-active counts at the 1.96 SD threshold rise with w_S for that reason alone. With exact
event times the validation calibration slope at w_S ≥ 100 is 0.07 to 0.17, so these fits do
overfit the training outcome, but the cause of non-convergence is not yet identified.

## Coefficient calibration (corrected 2026-11-05)

The registered calibration slope is the coefficient of a Cox regression of validation survival on
the fitted risk score. The first version computed it with `coxph`'s default `timefix = TRUE`. At
p = 2,000 the simulated times span about 10⁻³⁸ to 10¹⁰, and the default tie correction merges the
smallest ones (94 distinct times out of 120 in one checked dataset). That flattens every slope:
the true predictor itself gets a slope of 0.04 to 0.27 instead of about 1. The slopes below are
recomputed from the saved fits with `timefix = FALSE` (`grid/cal_slope_exact.csv`). C is
unaffected. The real-data runs are unaffected, because their times are on an ordinary scale.

Median validation calibration slope (1 = calibrated, above 1 = risk differences too small, below
1 = too large):

| p | Variance | True predictor | w_S = 1 | w_S = 10 | w_S = 100 | w_S = 10,000 | Risk SD, w_S = 1 / truth |
|---|---|---|---|---|---|---|---|
| 50 | 0.25% | 1.05 | 1.04 | 0.88 | 0.85 | 0.85 | 3.4 / 3.6 |
| 50 | 1% | 1.01 | 1.07 | 0.87 | 0.85 | 0.85 | 3.9 / 4.1 |
| 50 | 5% | 0.98 | 1.07 | 0.87 | 0.85 | 0.85 | 4.9 / 5.1 |
| 2,000 | 0.25% | 1.04 | 1.45 | 0.26 | 0.12 | 0.07 | 4.4 / 9.5 |
| 2,000 | 1% | 0.99 | 2.16 | 0.71 | 0.17 | 0.08 | 6.3 / 15.5 |
| 2,000 | 5% | 0.99 | 3.23 | 1.74 | 0.16 | 0.10 | 7.5 / 25.6 |

At p = 50, w_S = 1 is calibrated and larger weights are mildly overconfident. At p = 2,000,
w_S = 1 is attenuated by a factor of about 1.5 to 3: its risk score has a third to a half of the
true predictor's spread. w_S = 10 is closer to calibrated at the two larger variances, and
w_S ≥ 100 overfits. A weight between 1 and 10, or one chosen by a calibration criterion, is
where the supervision strength is right in this test bed.

| p | Variance | β* covered (w_S = 1) | Max \|E β\| | Validation C |
|---|---|---|---|---|
| 50 | 0.25% | 1/10 | 3.8 | 0.902 |
| 50 | 1% | 2/10 | 4.1 | 0.908 |
| 50 | 5% | 3/10 | 5.5 | 0.920 |
| 2,000 | 0.25% | 0/10 | 5.4 | 0.944 |
| 2,000 | 1% | 0/10 | 7.3 | 0.972 |
| 2,000 | 5% | 0/10 | 8.3 | 0.986 |

The β* coverage failure is real. At p = 2,000 the median matched E β is 6.9 and 7.7 at 1% and
5% variance, against a median β* of 16.2 and 26.0: the intervals miss because the coefficients
are shrunk toward zero.

**S-10 (corrected 2026-11-05).** The first version skipped S-10 on the grounds that raising
survival's share made calibration worse. With the corrected slopes that reasoning does not hold:
at p = 2,000 a moderately larger survival weight moves calibration toward 1 before overfitting
starts. The open question is how to set the supervision strength from the data while keeping the
fit stable. Downweighting genomics (S-10), a calibrated likelihood power (SafeBayes or
Holmes–Walker), and loss normalization as in DeSurv are all candidates.

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

The attenuation at p = 2,000 is carried as a known limitation: at w_S = 1 the supervision is
weaker than the data support. Real-data C is unaffected by the slope error.
