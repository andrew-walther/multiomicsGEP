# Simulation design for the ssbmf w = 1 panel (registered before any run)

Written 2026-10-26, before the first panel run. This is the reduced local
panel of the rebuild plan (Task 10, run locally instead of on Longleaf): 5
datasets per cell instead of 20. The full panel (20 datasets, semi-synthetic
designs) runs on Longleaf after review.

## Data-generating model

`ssbmf_simulate()`: Y_m = 1 mu_m^T + L F_m^T + E_m with nonnegative
activities (exponential, 20% zeros) and nonnegative loadings on disjoint
sets of 10% of each block's features; noise in groups of 20 features that
share a nuisance factor unrelated to the programs (within-group noise
correlation 0.3, `feature_correlation = "block"`). Survival is Weibull
(shape 1.5) with log hazard eta = sum_k beta_k Z_k / sd(Z_k), Z_k the
projection of the training-centered data on the true loadings; 30%
censoring. n = 200 training and 200 validation samples. K_true = 4.

## Cells

| Cell | Blocks (features) | Noise SD | beta (per SD of true score) | Purpose |
|---|---|---|---|---|
| A  | expression 2,000 | 1 (lower) | (0.7, 0.35, 0, 0) | one block, O1–O3, O7 |
| A2 | expression 2,000 | 2 (higher) | (0.7, 0.35, 0, 0) | one block, higher noise |
| B  | expression 2,000 + methylation 10,000 | 1 | (0.7, 0.35, 0, 0) | two blocks, O1–O3, O7 |
| C  | expression 2,000 | 1 | (0, 0, 0, 1.0); program 4 explains about 1% of variance | O4, one block |
| C2 | expression 2,000 + methylation 10,000 | 1 | (0, 0, 0, 1.0); program 4 about 1% in each block | O4, two blocks |
| C3 | expression 2,000 | 1 | (0, 0, 0, 1.0); program 4 explains about 0.25% of variance | O4, harder |

Cell C3 was added after one pilot fit (cell C, seed 1, K_init = 6), in
which all three arms recovered the 1% program (|r| = 0.97–0.98), so that
cell cannot show whether supervision helps. Nothing about the method or its
settings was changed after the pilot.

K_init in {4, 6, 8} (K_true, K_true + 2, 2 K_true). Datasets: seeds 1–5 per
cell (the same 5 datasets for every K_init and arm).

## Arms

1. `ssbmf_fit` at w_S = 1, F prior point-Laplace (provisional default).
2. `ssbmf_fit` at w_S = 1, F prior point-exponential.
3. Two-step comparator: `ssbmf_fit` with w_S = 0 (unsupervised EBMF with
   the same model, priors, initialization and intercept as arm 1), then
   `survival::coxph` on the projection scores Z = sum_m Ytilde_m F_m of the
   programs that are not inactive. Validation risk = Cox linear predictor of
   the validation projection scores.

## Metrics per fit

- O1: number of steps that lowered the objective by more than
  1e-8 x |objective| (`step_decreases`), worst relative step.
- O2: converged (relative objective change < 1e-6 within 1,000 sweeps).
- O3: for each true program with pve >= 2%, best match over fitted programs
  by |r| of the loadings concatenated over blocks; recovered if >= 0.9.
  Criterion: >= 90% of (fit, program) pairs.
- O4 (cells C, C2): the low-variance program is recovered if its best match
  has |r| >= 0.7, that fitted program is survival-active, and the sign of its
  coefficient (after aligning the loading sign with the truth) is positive.
  Survival-active: lfsr <= 0.05 for ssbmf (one-sided 95%); Wald p <= 0.10
  for the two-step Cox fit (the same two-sided 90% interval excluding 0).
  Criterion: >= 80% of datasets and more often than the two-step arm.
- O7: for every true program matched at |r| >= 0.9 (and the low-variance
  program when matched at >= 0.7), whether the 90% posterior interval of
  beta, on the per-SD-of-score scale and sign-aligned, covers the true value.
  Criterion: coverage between 80% and 95%.
- Also: validation C (Harrell, risk = eta, `reverse = TRUE`), sweeps,
  runtime, rejected steps, survival share of the F precision for the
  matched prognostic programs.

## Rules

No parameter is changed after seeing results. If O1–O3 fail, the work
returns to the derivation. If O4 fails at w_S = 1, the balance options of
rebuild Part 4 are tried in the plan's order, and any survival weight
chosen by training cross-validation is reported as the last resort.

## Balance panel (rebuild Part 4), added after the w = 1 panel

O4 was not met at w_S = 1, so the balance options were run on the O4 cells
(C, C2, C3) with the same datasets and metrics, in the rebuild plan's order:

1. `ssbmf_eff`: effective-feature block weight (no tuning parameter), all
   K_init.
2. `ssbmf_cvw`: w_S chosen by 5-fold training CV with the one-standard-error
   rule (smallest w_S, then smallest K), the last resort. Grid {1, 10, 100}.
   The grid {1, 1e2, 1e3, 1e4, 1e5} first written in the runner was reduced
   after one pilot dataset (cell C3, seed 2, K_init = 6): fits at
   w_S >= 1,000 took 2.5-5 minutes each with 400-2,300 rejected steps and
   gave no better loading recovery (|r| 0.41-0.42) or validation C than
   w_S = 100. One-block cells use all K_init; the two-block cell C2 uses
   K_init = 6 only, for run time.
