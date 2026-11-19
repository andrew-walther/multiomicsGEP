# Near-threshold simulation at K = 5 with four survival-signal settings

Registered 2026-11-19, before any fit. Author: Andrew Walther.

## Why

The November 20 progress-book figure compares validation C of the supervised fit
(w_S = 1) and the two-step fit in the near-threshold simulation of November 5
([registration](2026-11-05-near-threshold-plan.md)). That run used K = 7, the
working rank at the time, and two survival settings: very strong (true-risk C 0.95)
and moderate (0.74). The working model now uses K = 5. Neither setting resembles the
real data, where fitted models reach C of about 0.63 to 0.66. This rerun shows where
supervision at w_S = 1 helps and where it has no effect, including a setting close to
the real data.

This is a descriptive figure rerun, not a model-selection study. It adopts nothing
and changes no working-model setting.

## Design

- **Generator:** unchanged, `yfb_decision_data(seed, p = 2,000, share, truth_mode =
  "shared")`. Three true programs, each loading on 15% of genes. Program 3 is
  prognostic and scaled to a planted share of 0.05%, 0.1% or 0.2% of expression
  variance. 400 training and 120 validation patients, Weibull survival with shape 1.5
  and about 30% censoring. Seeds 1 to 10.
- **Survival settings.** The true log hazard is the patient's Program 3 score times a
  coefficient. Three settings rescale the training true risk score to a fixed SD,
  the same way the November 5 moderate setting does. The fourth keeps the original
  coefficient of 3.3.

  | Setting | Runner name | True risk score SD | True-risk C (calibration) |
  |---|---|---|---|
  | Real-data-like | `real_like` | 0.8 | 0.67 |
  | Moderate | `realistic` (November 5 name) | 1.1 | 0.74 |
  | Middle | `middle` | 2.5 | 0.86 |
  | Very strong | `strong` | about 6 to 10 (unscaled) | 0.95 |

  The SDs for the new settings were chosen on 2026-11-19 from generated data only, with
  no model fits: true-risk C averaged over seeds 1 to 5 and the three shares. SD 0.8
  is the closest setting whose ceiling (0.67) is consistent with the C that fitted
  models reach on the real data. SD 2.5 puts the middle setting near halfway between
  0.74 and 0.95.
- **Arms:** two-step (w_S = 0 loadings, then Cox on all live program scores) and
  supervised YFB (w_S = 1, survival-block weight, mean-field β). Both use
  point-Laplace L and F, Normal β, canonicalization and K = 5.
- **Size:** 10 seeds × 3 shares × 4 settings × 2 arms = 240 fits.

## Outcomes (same definitions as November 5)

- Validation C of each arm's full risk score Y~val~Fβ on the 120 validation patients.
- Recovery: the fitted program with the largest absolute loading correlation with the
  true Program 3 has |r| ≥ 0.7, |β/SD(β)| ≥ 1.96 and the correct survival direction.
- Convergence, sweeps and false-active programs, reported as before.
- **New: loading similarity between the two arms.** For each dataset and each live
  supervised program (column norm above 10⁻⁸), the largest absolute correlation of its
  loadings with any two-step program. The dataset value is the mean over live programs,
  and its minimum is also kept. It measures whether supervision changed the
  factorization.

## Code and outputs

`results/benchmark_sim/run_yfb_near_threshold.R` gains optional environment settings
for the output folder, K, survival settings and arms, plus the two new settings.
Unset, the runner reproduces the November 5 design. Outputs go to
`results/benchmark_sim/outputs/yfb_near_threshold_k5_2026-11-19/`: `metrics.csv`,
`loading_similarity.csv`, `task_status.txt` and `provenance.txt` are committed, and
per-fit checkpoints stay ignored. The November 5 results are not touched.

Any task error, nonconverged fit or undefined C is reported with the figure, not
dropped.

## Result (2026-11-19)

All 240 fits converged within 29 sweeps, with no task errors and no undefined C.
Supervised fits had 4 or 5 live programs and two-step fits the same number in every
dataset. The observed true-risk C on the validation patients was 0.69, 0.74, 0.86 and
0.95 for the four settings. The real-data-like value is slightly above the 0.67 from the
calibration, which used seeds 1 to 5 only. Per-cell summaries are in
`cell_summary.csv`, written by `docs/progress_book/figures/plot_2026-11-20.R`.

| Setting | Share | Mean C gain | Recovered, supervised vs two-step | Median loading r |
|---|---|---|---|---|
| Real-data-like | 0.05% / 0.1% / 0.2% | +0.017 / +0.005 / +0.003 | 0 vs 0 / 8 vs 8 / 10 vs 10 | 0.981 / 0.998 / 0.999 |
| Moderate | 0.05% / 0.1% / 0.2% | +0.030 / +0.008 / +0.002 | 0 vs 0 / 9 vs 8 / 10 vs 10 | 0.978 / 0.995 / 0.998 |
| Strong | 0.05% / 0.1% / 0.2% | +0.077 / +0.020 / +0.001 | 3 vs 0 / 10 vs 8 / 10 vs 10 | 0.926 / 0.982 / 0.996 |
| Very strong | 0.05% / 0.1% / 0.2% | +0.112 / +0.039 / +0.008 | 4 vs 0 / 10 vs 8 / 10 vs 10 | 0.873 / 0.975 / 0.994 |

Supervision changes the loadings and helps prediction only as the survival signal
strengthens. In the real-data-like setting the two fits find nearly the same programs.
The figure and text are in the November 20 progress-book chapter. Nothing is adopted.
