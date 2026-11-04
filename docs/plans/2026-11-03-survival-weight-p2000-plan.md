# Attempt S-9: the survival-weighted YFB fit at 2,000 genes

Registered 2026-11-04, before any S-9 fit was run. The work order is
[`2026-10-30-downweighted-genomics-prompt.md`](../prompts/2026-10-30-downweighted-genomics-prompt.md),
Part A. The attempts log carries the matching S-9 entry. Nothing here changes Gate S or any S-7
or S-8 record.

## Goals and deliverables (kept current)

Decision of 2026-11-04 (DECISIONS.md): Gate S is reported literally, but real-data work follows
if the fit works in simulation, without waiting for a literal pass.

1. **S-9 in simulation:** per-weight curves (recovery, false-active, survival share, β calibration,
   convergence, validation and CV C), the Gate S verdict, and whether training CV picks the same
   weight as the truth-based best arm.
2. **If it works, on the real data:**
   - the programs the fit identifies, and which are survival-active;
   - out-of-sample C from training cross-validation;
   - the feature-set-size curve, 2,000 genes to all genes, with numbers and an elbow plot;
   - the survival-screening comparison rerun with the survival-weighted fit.
   - after those choices are fixed by training CV, one external-cohort scoring of the fixed model
     (assessment only; external cohorts choose nothing).
3. **If it does not work:** S-10 (downweighted genomics), registered first.
4. **After the single-modality items succeed:** extend to the M-modality model.
5. **Reporting:** keep α_F = 0 results, labelled; present the supervised fit as preferred once it
   is carried forward. Update the October 30 chapter, tracker, DECISIONS.md and review handoff.

## Question

`fit_ssbmf_yfb(canonicalize = TRUE, w_S)` adds a Cox term to the gene-by-gene F update:

A_j = τ_j Σ_i E[l_ik²] + w_S E[β_k²] Σ_i w_i y_ij²,

and the matching term in B_j. On October 15, with w_S = 1 at p = 50 and n = 400, this fit
recovered a low-variance planted program in 17 of 20 datasets where the two-step methods recovered
none. The genomics term in A_j grows with n and its weight does not shrink with p, while the Cox
term has about one contribution per event, so at fixed w_S survival's share of the F precision is
expected to fall as p grows. S-9 asks whether the fit still recovers the program at p = 2,000,
and whether a larger or p-scaled weight is needed.

## Prior run found before registration (step A1)

The S-7 grid already contains this fitter as its `canonical_w1` arm (w_S = 1, same priors, K,
sweeps and α as below), at p = 50 and 2,000 on the same seeds. It was registered in S-1 as a
baseline but never reported. Its saved metrics at p = 2,000:

| Planted variance | Recovered | Mean false-active | Converged | Final survival share, matched program (median) |
|---|---|---|---|---|
| 0.25% | 10/10 | 1.0 | 10/10 | 0.27 |
| 1% | 10/10 | 1.5 | 10/10 | 0.16 |
| 5% | 10/10 | 1.1 | 10/10 | 0.03 |

The S-7 unsupervised reference recovered 0, 0 and 1 of 10 in the same cells. Read against Gate S,
w_S = 1 would qualify only in the 0.25% cell (the other two exceed mean false-active 1), so it
would fail the gate on false activity, not on recovery or convergence.

This does not settle the question, for one reason: S-7 gave this fitter the nonnegative,
uncentered X/sd input prepared for the new fitter, while `fit_ssbmf_yfb()` is documented for
column-centered input. Without centering, one program has to absorb the column means. S-9 reruns
the fitter on the input it was written for and adds the weight range.

## Data

- The S-7 corrected shared-loading generator, `yfb_decision_data(seed, p, v, truth_mode = "shared")`
  in `results/benchmark_sim/run_yfb_supervision_decision.R`: seeds 1–10, p ∈ {50, 2,000}, planted
  raw centered variance share v ∈ {0.25%, 1%, 5%}, n = 400 training and 120 validation patients,
  coefficient 3.3 per unit normalized prepared projection, noise 0.25, about 30% censoring.
  60 datasets, generated exactly as in S-7.
- **Preprocessing for this fitter:** the signed per-cohort z-score. The generator's training matrix
  is X/sd within the cohort, so centering it gives the z-score; the planted direction and the
  generated risk are unchanged. The validation cohort is z-scored within itself.
- **Why centered input for this fitter (checked 2026-11-04 against the realignment plan's Phase P
  and the two Gemini documents).** The preferred setup for the nonnegative model is X/sd within
  each cohort, kept nonnegative for reconstruction, with centering only inside the Cox projection
  (Phase P item 3; Gemini guide, Issue 2 of the reconciliation list). That setup belongs to the new
  fitter, whose point-exponential priors need nonnegative input. `fit_ssbmf_yfb()` has signed
  point-Laplace L and F priors, no feature intercept, and uses the same matrix for reconstruction
  and for the projection. Given uncentered input, it must spend a program on the column means.
  Given the centered matrix, it sees what the Gemini setup's Cox projection sees, and its signed
  reconstruction is well defined. The point-exponential prior on z-scored input collapsed in the
  October 15 test bed, so this fitter keeps point-Laplace. The S-7 uncentered `canonical_w1`
  result stays on record as the mismatched-input comparison.
- **Variance share under each preprocessing**, recorded per dataset: the planted component's share
  of the centered sum of squares (what the z-scored fitter sees, equal to the S-7 prepared share)
  and its share of the uncentered sum of squares (what the nonnegative uncentered input carries).

## Arms

All arms are the S-7 `canonical_w1` call with the z-scored input: point-Laplace L and F,
per-program Normal β, K = 7, at most 300 sweeps, tolerance 10⁻⁵ and α from `config/globals.yml`,
no sign correction, no supervised start.

1. w_S = 0: the same fitter without survival in the F update (β is still fitted). Same-fitter
   reference, and the source of the share0.2 weight.
2. w_S ∈ {1, 10, 100, 1,000, 10,000}.
3. w_S = p/50: 1 at p = 50 (the w_S = 1 fit is reused, not refitted) and 40 at p = 2,000.
4. share0.2: the w_S that makes the median survival share of the F precision equal 0.2 at the
   w_S = 0 fit. Shares are w A_S/(A_G + w A_S) over genes and live programs, with A_G and A_S
   computed at the w_S = 0 fit's final state and Cox working weights at its predictor; the S-2
   solver (`yfb_calibrate_gamma()`, no damping) finds w. The weight is fixed before the supervised
   fit and never revised. If it cannot be computed (no live program), the arm is a failed fit.

The S-7 unsupervised reference (new fitter, γ = 0) is read from its saved checkpoints.
A smoke check confirms that the generator and that reference reproduce the saved S-7 values.

## Replication cell

The A1 calibrated bed (`results/benchmark_sim/run_supervision_diagnostics.R`, arm set `s9_rep`:
p = 50, n = 400, β₃ = 3.3, seeds 20261016–20261035, noise 0.25 and 1, w_S = 1, supervised start off).
On October 30 this arm recovered 18 of 20 at noise 0.25 (October 15: 17 of 20 for the selected
supervised fit). **Replication criterion:** at least 15 of 20 recovered at noise 0.25, and every
per-seed difference from the saved A1 `results.csv` explained. If it fails, the grid result is not
interpreted until the difference is diagnosed.

## Metrics, per fit

- **Recovered:** best |r| between the planted loading direction and E[F_k] at least 0.7, the
  matched program survival-active (|E β_k| ≥ 1.96 posterior SD, SD > 0), coefficient of the correct
  sign. Same as S-7 and A1.
- **False-active:** survival-active programs other than the match.
- **Converged:** the fitter's own flag (relative ELBO change below tolerance before 300 sweeps).
- **Survival share:** the final in-fit median over genes of w_S A_S/A_F, for the matched program
  and the median over live programs.
- **β calibration:** (i) the validation Cox calibration slope of the risk score (1 is calibrated,
  below 1 is overconfident); (ii) whether E β_k ± 1.96 SD contains β*, the coefficient of the true
  training predictor projected on the matched Z_k.
- **Validation C** on the 120 validation patients, orientation fixed at fit time.
- **Training cross-validated C:** five folds stratified by event, fixed per dataset (seed + 7000);
  each part z-scored within itself; share0.2 recalibrated inside each fold from that fold's w_S = 0
  fit.
- **Mechanism diagnostics:** sweeps used, ELBO decreases, max |E β|, reconstruction loss
  ‖Y − E[L]E[F]ᵀ‖²/‖Y‖², min τ_j, and lost-variance records (E β² − (E β)² ≤ 0 on a live program).
- A fit that errors counts as not converged, not recovered and seven false-active programs.

## Gate S (unchanged) and how it is read

At p = 2,000, a single arm must qualify in at least two of the three variance cells. A cell
qualifies when recovery exceeds the S-7 unsupervised reference by at least 3 of 10 datasets,
mean false-active is at most 1 and at least 9 of 10 fits converge (the conservative reading
registered in S-1). The gate is read at the best registered arm, which uses the planted truth and
is allowed only in simulation. The same verdict is also reported with the canonical w_S = 0 fit as
the reference; the S-7 reference controls.

**How the weight would be chosen in practice.** In each dataset the arm with the highest training
CV C is recorded. Per cell, report how often each arm is chosen, the recovery, false-active and
convergence of the CV-chosen fits, and whether CV picks the arm that the gate reads. This is
reported, not gated.

## Reports

- Per-weight curves at p = 50 and 2,000 for recovery, false-active, survival share, calibration
  slope, β coverage, convergence, validation C and CV C, plotted against w_S.
- How the survival share changes from p = 50 to p = 2,000 at each weight, and whether a fall in
  share explains any loss of recovery.
- Non-convergence at large w_S, read as the τ mechanism (falling min τ_j with reconstruction loss
  above 1), the β mechanism (rising |E β| and false-active counts with bounded reconstruction) or
  something new.
- Findings: `docs/findings/2026-11-04-survival-weight-p2000-findings.md`. It also names, as a
  follow-up not implemented here, a calibration rule for the tempering weight (SafeBayes or the
  Holmes–Walker power-likelihood calibration); a prior on w_S is not a substitute.

## Next steps fixed in advance

- **If S-9 passes:** run step S3 of the realignment plan with the passing weight (real-data
  training CV on the single-log data, no external cohort), then the null-data calibration check of
  Part B step 3.2 at that weight. Both get their own registered entries first.
- **Feature-set size (requested 2026-11-04, following the advisor's suggestion).** If a weight
  passes, the S3 training CV is run over a range of feature-set sizes, from 2,000 genes to the full
  gene set. The smallest set whose training CV C and program correlations match the largest set
  within one standard error is preferred, as the more parsimonious model. This uses training data
  only. The earlier fixed-K feature-count curve (attempt Task5-F1) used the historical α_F = 0
  model and is not a substitute.
  Report the numbers per size (CV C ± SE, convergence, program correlations with the largest set)
  and an elbow plot of CV C against feature-set size.
- **Survival screening before the joint fit (requested 2026-11-04).** Rerun the October 16
  feature-selection comparison (`docs/plans/2026-10-16-feature-selection-experiment-plan.md`,
  runner `results/feature_selection/run_feature_selection_comparison.R`) with the passing
  survival-weighted fit in place of the α_F = 0 fit. Arms: the unsupervised mean/variance screen
  U(N) and the same genes followed by a fold-internal cohort-stratified Cox screen S(N), at the
  same N grid. The Cox screen runs on fold-training patients only, so held-out patients never
  inform it. The question is whether screening adds anything once survival already shapes the
  loadings; the October 26 α_F = 0 results (S(N) above U(N) by 0.01–0.04 in CV C) stay on record
  beside the new ones, labelled as the unsupervised-loadings model.
- **If S-9 fails:** run Part B (S-10, downweighted genomics) under its own registration, and
  say which failure it is meant to fix.

## Software and compute

`results/benchmark_sim/run_yfb_survival_weight_p2000.R`: PSOCK workers, one dataset per work item,
one BLAS thread each, checkpoints per dataset and arm under
`results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid/`. About 42 fits per dataset
including CV. No fitter or update code changes; the replication adds one arm set to the A1 runner.
No external or real cohort is used.
