# Attempt S-14: supervision near the detection threshold — findings

Registration: [`docs/plans/2026-11-05-near-threshold-plan.md`](../plans/2026-11-05-near-threshold-plan.md).
Runner: `results/benchmark_sim/run_yfb_near_threshold.R`, commit `bbafa4c`. Output:
`results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05/metrics.csv` (240 fits, no errors,
all converged). Synthetic data only. The two-step arm was rerun with exact training times
(`coxph(..., timefix = FALSE)`), commit noted in the attempts log; the default `timefix = TRUE`
merges near-equal extreme times and attenuates the two-step Cox coefficients. The rerun changed
only the strong-survival two-step slopes and false-active counts (0.1%: slope 1.25 → 1.22;
0.2%: slope 1.53 → 1.14, false-active 1.5 → 1.8). The default-timefix checkpoints are kept in
`fits_twostep_timefix_default/` (not in Git).

## Verdict against the registered criteria

| Criterion | Result |
|---|---|
| Primary: A1 at w_S = 1 recovers ≥ 3 more datasets than the two-step fit in ≥ 2 of 3 strong cells, with mean false-active ≤ 1 and ≥ 9/10 converged | **Fails** (0 of 3 cells qualify) |
| Secondary: mean paired validation C difference (A1 w_S = 1 minus two-step) positive in ≥ 2 of 3 strong cells | **Passes** (3 of 3; A1 has the higher C in 29 of 30 strong datasets) |

The 0.05% cell has the recovery gain (5 vs 0) but a mean false-active count of 2.5. The 0.1% cell
gains 2 (10 vs 8). The 0.2% cell gains 0 (10 vs 10).

## Results by cell

Recovered counts are out of 10. "FA" is the mean number of false-active programs. C is the mean
validation C, and the slope is the median exact-times calibration slope. All 240 fits converged.

| Survival | Share | Arm | Recovered | Mean \|r\| | FA | C | Oracle C | Slope |
|---|---|---|---|---|---|---|---|---|
| strong | 0.05% | two-step | 0 | 0.20 | 2.8 | 0.644 | 0.936 | 1.19 |
| strong | 0.05% | A1 w_S = 1 | 5 | 0.62 | 2.5 | 0.800 | 0.936 | 1.58 |
| strong | 0.05% | A1 w_S = 5 | 9 | 0.76 | 1.0 | 0.844 | 0.936 | 0.65 |
| strong | 0.05% | A2 w_S = 1 | 5 | 0.63 | 2.4 | 0.801 | 0.936 | 1.58 |
| strong | 0.1% | two-step | 8 | 0.73 | 2.4 | 0.861 | 0.948 | 1.22 |
| strong | 0.1% | A1 w_S = 1 | 10 | 0.87 | 1.6 | 0.901 | 0.948 | 1.08 |
| strong | 0.1% | A1 w_S = 5 | 10 | 0.77 | 1.1 | 0.873 | 0.948 | 0.27 |
| strong | 0.1% | A2 w_S = 1 | 10 | 0.86 | 1.6 | 0.897 | 0.948 | 1.12 |
| strong | 0.2% | two-step | 10 | 0.89 | 1.8 | 0.928 | 0.958 | 1.14 |
| strong | 0.2% | A1 w_S = 1 | 10 | 0.92 | 1.0 | 0.936 | 0.958 | 1.34 |
| strong | 0.2% | A1 w_S = 5 | 10 | 0.85 | 0.9 | 0.920 | 0.958 | 0.43 |
| strong | 0.2% | A2 w_S = 1 | 10 | 0.92 | 1.0 | 0.936 | 0.958 | 1.31 |
| realistic | 0.05% | two-step | 0 | 0.20 | 1.6 | 0.592 | 0.740 | 1.19 |
| realistic | 0.05% | A1 w_S = 1 | 0 | 0.34 | 1.9 | 0.630 | 0.740 | 1.36 |
| realistic | 0.05% | A1 w_S = 5 | 0 | 0.40 | 0.6 | 0.650 | 0.740 | 0.24 |
| realistic | 0.05% | A2 w_S = 1 | 0 | 0.34 | 1.7 | 0.630 | 0.740 | 1.31 |
| realistic | 0.1% | two-step | 8 | 0.73 | 1.1 | 0.705 | 0.736 | 1.24 |
| realistic | 0.1% | A1 w_S = 1 | 9 | 0.77 | 1.0 | 0.712 | 0.736 | 1.16 |
| realistic | 0.1% | A1 w_S = 5 | 0 | 0.30 | 0.4 | 0.634 | 0.736 | 0.04 |
| realistic | 0.1% | A2 w_S = 1 | 9 | 0.77 | 1.0 | 0.712 | 0.736 | 1.16 |
| realistic | 0.2% | two-step | 10 | 0.89 | 0.8 | 0.727 | 0.735 | 1.11 |
| realistic | 0.2% | A1 w_S = 1 | 10 | 0.90 | 0.5 | 0.728 | 0.735 | 1.12 |
| realistic | 0.2% | A1 w_S = 5 | 7 | 0.69 | 0.8 | 0.671 | 0.735 | 0.23 |
| realistic | 0.2% | A2 w_S = 1 | 10 | 0.90 | 0.5 | 0.728 | 0.735 | 1.12 |

Mean paired validation C difference from the two-step fit, with the number of the 10 datasets in
which A1 at w_S = 1 has the higher C:

| Survival | Share | A1 w_S = 1 | A1 w_S = 5 | A2 w_S = 1 | A1 w_S = 1 higher |
|---|---|---|---|---|---|
| strong | 0.05% | +0.156 | +0.200 | +0.157 | 10 |
| strong | 0.1% | +0.041 | +0.013 | +0.036 | 10 |
| strong | 0.2% | +0.008 | −0.008 | +0.008 | 9 |
| realistic | 0.05% | +0.038 | +0.058 | +0.038 | 9 |
| realistic | 0.1% | +0.008 | −0.070 | +0.008 | 8 |
| realistic | 0.2% | +0.002 | −0.055 | +0.002 | 8 |

## What this shows

1. **Supervision helps prediction near the threshold, and the gain grows as the program weakens.**
   With strong survival, A1 at w_S = 1 has the higher validation C in 29 of the 30 paired strong
   datasets. The gain is +0.156 at a share of 0.05%, where the two-step fit does not find the
   program at all (mean |r| 0.20), and +0.008 at 0.2%, where both fits find it.
2. **The registered recovery criterion fails on false activity, not recovery.** In the 0.05%
   strong cell the supervised fit recovers 5 programs to the two-step fit's 0, but averages 2.5
   false-active programs. Every arm has this problem in the near-threshold cells. The two-step fit
   averages 2.8 false actives in the same cell and 2.3 across the strong cells. When the planted
   direction is only partly resolved, other programs carry part of it and become significant. The
   false-active count then measures how the direction is split across programs rather than
   spurious survival signal. The criterion was registered as written and is reported as failed.
3. **A1 at w_S = 5 trades calibration for recovery.** It recovers 9 of 10 in the 0.05% strong cell
   with 1.0 false actives, but its calibration slopes are 0.04 to 0.65 in every cell, and with
   realistic survival it loses recovery (0 of 10 at 0.1%, against 9 for w_S = 1). Under realistic
   survival a fixed w_S = 5 overfits. This repeats S-9: a weight large enough to move the loadings
   in a strong-survival cell overfits a weak one, so the weight has to be chosen from the data,
   which is what S-13 tests.
4. **With realistic survival (oracle C about 0.74) supervision adds little.** At 0.1% and 0.2% all
   w_S = 1 fits and the two-step fit reach C within 0.01 to 0.03 of the oracle. At 0.05% no arm
   recovers the program. The supervised fits gain +0.04 in C but stay 0.11 below the oracle.
5. **The joint q(β) arm (A2) matches A1 to within 0.005 in C** and recovers the same datasets. On
   these data the mean-field β posterior is not the limiting factor.
6. **Calibration at w_S = 1 is close to 1 under realistic survival (1.12 to 1.36).** With strong
   survival the slopes are 1.08 to 1.58: the coefficients are attenuated, as in the corrected S-9.

## Limits

- Ten datasets per cell, one generator (`truth_mode = "shared"`, K = 7, p = 2,000). The planted
  shares are below the shares of the S-9 cells.
- The two-step false-active count uses Cox Wald tests on the live program scores. The supervised
  arms use 1.96 posterior SDs. The two rules are not the same test.
- The realistic-survival rescaling was calibrated from the generator alone (seeds 101–106). No
  registered fit informed it.

## Added 2026-11-06: how much the other programs contribute

This analysis is descriptive and post hoc, from `results/benchmark_sim/summarize_false_active_contribution.R`
(output `false_active_contribution.csv`). Activity is not stored in the S-14 checkpoints, so every
program other than the matched one is counted, active or not. Standardized |β| is the log hazard
ratio per training SD of the program score. ΔC is the validation C of the full risk score minus
that of the matched program alone.

| Survival | Share | Arm | Mean FA | Risk variance from other programs (median) | Std. \|β\| matched (median) | Largest other std. \|β\| (median) | ΔC |
|---|---|---|---|---|---|---|---|
| strong | 0.05% | two-step | 2.8 | 69% | 0.41 | 0.24 | +0.022 |
| strong | 0.05% | A1 w_S = 1 | 2.5 | 8.7% | 1.84 | 0.33 | +0.012 |
| strong | 0.1% | two-step | 2.4 | 1.4% | 2.80 | 0.27 | +0.008 |
| strong | 0.1% | A1 w_S = 1 | 1.6 | 0.5% | 3.76 | 0.23 | +0.005 |
| strong | 0.2% | two-step | 1.8 | 0.4% | 5.02 | 0.24 | +0.003 |
| strong | 0.2% | A1 w_S = 1 | 1.0 | 0.1% | 4.79 | 0.16 | +0.001 |
| realistic | 0.05% | two-step | 1.6 | 57% | 0.29 | 0.16 | +0.004 |
| realistic | 0.05% | A1 w_S = 1 | 1.9 | 39% | 0.49 | 0.19 | +0.013 |
| realistic | 0.1% | A1 w_S = 1 | 1.0 | 1.3% | 1.05 | 0.14 | +0.004 |
| realistic | 0.2% | A1 w_S = 1 | 0.5 | 0.6% | 1.08 | 0.08 | 0.000 |

All arms and cells are in the CSV. Two cases differ:

- **Where the program is resolved (0.1% and 0.2%).** The other programs carry at most 1.3% of the
  risk-score variance, and removing them changes C by at most 0.008. Their counts fail the
  registered criterion, but they barely affect the risk score.
- **Where the program is only partly resolved (0.05%).** The other programs carry 9% to 69% of the
  risk-score variance and add to C. They are not noise. They hold pieces of the planted
  direction that the matched program missed, and the split is much larger for the two-step fit
  (69%) than for the supervised fit (8.7%) with strong survival.

The registered criterion counts both cases the same way. For prediction, the first case is
negligible and the second is part of the signal.
