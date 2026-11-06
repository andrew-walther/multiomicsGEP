# Attempt S-13: survival-block weight and joint q(β), weight chosen by CV — findings

Registration: [`docs/plans/2026-11-05-supervision-strength-plan.md`](../plans/2026-11-05-supervision-strength-plan.md)
(section "Registration"). Fitter commit `01095a2`. Runner
`results/benchmark_sim/run_yfb_supervision_strength.R`. Summary
`results/benchmark_sim/summarize_yfb_supervision_strength.R`. Output
`results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05/`. Synthetic data only.

The grid has 480 signal fits: 2 arms × 4 cells × 10 seeds × 6 weights, each with five CV folds.
There are 360 null fits on permuted survival. One signal fit errored (weak cell, seed 9, A2,
w_S = 50), and 7 of the 80 w_S = 50 fits lost a CV fold. Neither affects the selected weights,
because the CV never chose a weight above 10.

The grid driver's log ends with a parse error. It was written after the task-status file, when
the R process read further into the runner script, which had been edited while the grid was
running (`bbafa4c` added the two-step mode). All 480 tasks had finished and checkpointed. The
summary reads only the checkpoints.

## Arms and cells

- **A1:** w_S multiplies the Cox term and KL(q(β)‖g_β) together. β is untempered, and F sees
  survival at weight w_S.
- **A2:** A1 with a joint Gaussian q(β) and per-program prior variances by type-II maximum
  likelihood.
- **Weights:** {1, 2, 5, 10, 20, 50}. The weight is chosen per dataset by five-fold CV partial
  log-likelihood, with ties going to the smaller weight.
- **Cells:** the S-9 p = 2,000 cells (planted variance 0.0025, 0.01, 0.05), a weak-survival
  cell (oracle C about 0.68), and permuted-survival nulls of the S-9 cells.
- **Matched reference:** the two-step fit, which takes the w_S = 0 loadings on the same centered
  input and then fits a Cox model on the live program scores. Its Cox fit was rerun with exact
  training times (see "Correction to the two-step reference").

## Verdict against the registered criteria

| Criterion | A1 | A2 |
|---|---|---|
| Gate S, S-9 definition (≥ 3 above the S-7 reference, FA ≤ 1, ≥ 9/10 converged, in ≥ 2 of 3 cells) | Passes (3 of 3) | Passes (3 of 3) |
| Median calibration slope in [0.8, 1.25] in ≥ 2 of 3 S-9 cells | **Fails** (1 of 3) | **Fails** (1 of 3) |
| ≥ 9/10 converged per cell at the selected weight | Passes | Passes |
| Null mean false-active ≤ 1 at the modal selected weight | **Fails** (1.1 at w_S = 5) | Passes (0.9 at w_S = 2) |

Neither arm passes every registered criterion.

**Gate S does not measure a supervision effect here.** Its reference is the S-7 unsupervised fit,
which was fitted to uncentered input and recovers 0 or 1 of 10. The matched two-step fit recovers
10 of 10 in every cell, as do both supervised arms. Against the matched reference, the recovery
gain is 0 in every cell. This is the same point registered as motivation for S-14.

## Results at the CV-selected weight

| Arm | Cell | Selected w_S (count) | Recovered | Two-step recovered | FA | Two-step FA | C | Two-step C | Slope | Two-step slope |
|---|---|---|---|---|---|---|---|---|---|---|
| A1 | 0.0025 | 1 (1), 2 (9) | 10 | 10 | 0.9 | 1.7 | 0.946 | 0.939 | 1.24 | 1.16 |
| A2 | 0.0025 | 2 (10) | 10 | 10 | 0.6 | 1.7 | 0.946 | 0.939 | 1.22 | 1.16 |
| A1 | 0.01 | 5 (10) | 10 | 10 | 0.0 | 0.9 | 0.972 | 0.972 | 1.27 | 1.10 |
| A2 | 0.01 | 5 (10) | 10 | 10 | 0.0 | 0.9 | 0.972 | 0.972 | 1.28 | 1.10 |
| A1 | 0.05 | 10 (10) | 10 | 10 | 0.0 | 0.2 | 0.986 | 0.986 | 1.45 | 0.96 |
| A2 | 0.05 | 10 (10) | 10 | 10 | 0.0 | 0.2 | 0.986 | 0.986 | 1.45 | 0.96 |
| A1 | weak | 1 (8), 2 (2) | 10 | 10 | 1.2 | 0.5 | 0.671 | 0.668 | 1.05 | 1.04 |
| A2 | weak | 1 (10) | 10 | 10 | 0.9 | 0.5 | 0.672 | 0.668 | 1.05 | 1.04 |

FA is the mean number of false-active programs, C is the mean validation C, and the slope is the
median exact-times calibration slope. The mean paired C difference from the two-step fit is
+0.007 in the 0.0025 cell (supervised higher in 10 of 10), +0.003 to +0.004 in the weak cell, and
0.000 in the 0.01 and 0.05 cells.

Per-weight results for every cell are in `weight_curves.csv` and `weight_curves.png`.

## What this shows

1. **CV partial log-likelihood picks a sensible weight without manual tuning.** The selected weight
   grows with the strength of the planted program: 1 or 2 in the weak and 0.0025 cells, 5 at
   0.01, 10 at 0.05. It never chose w_S = 20 or 50. At those weights the full fits stop converging
   (0 or 1 of 10 per cell), the survival share of the matched program approaches 1, and the
   calibration slope falls well below 1.
2. **At the selected weight, supervision matches the two-step fit and does not beat it.** Recovery
   is the same (10 of 10 everywhere). Validation C is equal in the stronger cells and higher by
   0.003 to 0.007 in the two weakest. In the two stronger S-9 cells, supervision has fewer false
   actives (0.0 against 0.9 and 0.2).
3. **Supervised coefficients stay attenuated in the stronger cells.** The selected-weight slopes
   are 1.24, 1.27 and 1.45, and the two-step slopes are 1.16, 1.10 and 0.96. The CV objective
   prefers a weight below the one that would give a slope of 1. In the 0.05 cell, w_S = 20 gives
   a slope of 1.03 but does not converge.
4. **The A1 null failure is narrow.** On permuted survival, w_S = 1 gives 0.2 to 0.4 false-active
   programs and w_S = 5 gives 0.9 to 1.1. The registered check uses each arm's modal selected
   weight from the signal cells (5 for A1, 2 for A2), not a weight chosen by CV on the null data.
5. **The joint q(β) arm matches the mean-field arm.** The difference is at most 0.004 in C, with
   the same recoveries and slightly fewer false actives in two cells.
6. **Survival-only directions.** In the weak cell, the 1 to 2 active programs with PVE below 1% are
   the false actives. In the 0.0025 cell the matched program itself is labelled a survival-only
   direction, because the planted program's own PVE is about 0.34%, below the 1% cut. That label
   reflects the planted share, not a detached direction.

## Correction to the two-step reference

The two-step arm first used `coxph()` with its default `timefix = TRUE` on the training times.
In the strong-survival S-9 cells, the simulated times span many orders of magnitude, so the
default merges near-equal times. That shrinks the two-step coefficients, giving validation
slopes of 1.80, 6.70 and 12.74. It is the same artifact found earlier in the S-9 validation
slopes. With exact times (`timefix = FALSE`), the two-step slopes are 1.16, 1.10 and 0.96, and
its validation C rises (for example 0.973 to 0.986 at 0.05). The default-timefix table is kept
as `twostep_timefix_default.csv`. All results above use the exact-times rerun, and the same fix
was applied to the S-14 two-step arm.

## Limits

- Ten datasets per cell and one generator.
- The S-9 cells are not a regime where supervision is needed for recovery. S-14 tests weaker
  planted programs.
- The two-step false-active count uses Cox Wald tests. The supervised arms use 1.96 posterior SDs.

## Added 2026-11-06: how much the false-active programs contribute

This analysis is descriptive and post hoc, from `results/benchmark_sim/summarize_false_active_contribution.R`
(output `false_active_contribution.csv`). It changes no registered verdict. For each CV-selected
fit, it splits the validation risk score into the matched program's part and the rest, and it
reports standardized coefficients (log hazard ratio per training SD of the program score).

| Arm | Cell | Mean FA | Risk-score variance from other actives (median) | Std. \|β\| matched (median) | Largest std. \|β\| among false actives (median) | C change from removing them |
|---|---|---|---|---|---|---|
| A1 | 0.0025 | 0.9 | < 0.1% | 6.0 | 0.13 | −0.001 |
| A2 | 0.0025 | 0.6 | < 0.1% | 6.1 | 0.07 | 0.000 |
| A1 | weak | 1.2 | 0.8% | 0.74 | 0.13 | 0.000 |
| A2 | weak | 0.9 | 0.7% | 0.74 | 0.11 | +0.002 |

The 0.01 and 0.05 cells have no false actives. Where they occur, the false-active programs have
coefficients about one fifth to one fiftieth of the matched program's, carry under 1% of the
risk-score variance, and change validation C by at most 0.002. For prediction they are
negligible. The registered false-active count does not measure this, so future registrations
should bound the contribution (risk-variance share or standardized coefficient) as well as the
count.
