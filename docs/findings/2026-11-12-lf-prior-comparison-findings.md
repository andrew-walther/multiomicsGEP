# Task A findings: point-exponential vs point-Laplace priors on L and F

Date: 2026-11-12. Author: Andrew Walther.

**Result.** Point-Laplace stays. Point-exponential L and F, fitted to nonnegative X/sd input
with the Cox projection centered, fails two of the four registered criteria. It averages
2.6 to 4.3 false-active programs per fit in every simulation cell (limit 1). On the real
training data its mean within-cohort CV concordance is 0.573 against 0.658 for point-Laplace,
a paired difference of −0.085 with SE 0.031. None of its five fold fits converged in 300
sweeps. On real data, in every fold, one point-exponential program's scores correlate 0.998
with the TCGA indicator, and every program's correlation is at least 0.56. Its programs
encode cohort membership, and none is survival-active.

The task F rank-transform arm (descriptive) also tracks the cohort label (largest correlation
1.000 in every fold, every program at least 0.72),
converges in no fold, and has within-cohort CV concordance 0.615 (difference from
point-Laplace −0.043, SE 0.016).

Registration: the 2026-11-12 "A" entry in the
[attempts log](../plans/2026-10-16-round2-attempts-log.md), task A of the
[week plan](../plans/2026-11-09-week-plan.md), pooled per DECISIONS.md 2026-11-12.
Runner: `results/benchmark_sim/run_lf_prior_comparison.R`. Aggregate outputs:
`results/benchmark_sim/outputs/lf_prior_2026-11-09/`.

## Design

Both arms use `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)` with Normal β, K = 7, no
pruning, at most 300 sweeps and tolerance 10⁻⁵, which is the working-model call.

| Arm | Input | L and F prior | Cox projection |
|---|---|---|---|
| PL (working model) | z-score within cohort | point-Laplace | as is (input already centered) |
| PE | X/sd within cohort, nonnegative | point-exponential | centered with one pooled mean (`center_projection = TRUE`) |
| RANK (task F, descriptive) | within-patient ranks over the D4 genes, each gene divided by its within-cohort SD of the ranks, nonnegative | point-exponential | centered with one pooled mean |

**Simulation.** Four cells of 10 seeds each: the S-9 generator at p = 2,000 with planted
variance shares 0.25%, 1% and 5%, and the S-14 near-threshold cell (share 0.05%, strong
survival). n = 400 training and 120 validation patients, one cohort, no strata.

**Real data.** The 273 TCGA and CPTAC training patients, the 2,218 D4 genes and the S-11
five folds (balanced by cohort and event, seed 2026), pooled without strata. The primary
metric is within-cohort concordance: only TCGA–TCGA and CPTAC–CPTAC pairs are compared,
pooled by pair count. This keeps a score that separates the two cohorts from being
rewarded. Pooled and per-cohort concordance are reported beside it.

All 98 planned fits returned (80 simulation, 18 real), with no errors and no fit warnings.
Real fits took 52 to 81 seconds and simulation fits about 8 to 15 seconds on one core each.

## Simulation

| Cell | Arm | Recovered | Mean false-active | Converged | Median best \|r\| | Median validation C | Median calibration slope | Median programs with nonzero loadings |
|---|---|---|---|---|---|---|---|---|
| 0.05% (S-14) | PE | 8/10 | 2.6 | 10/10 | 0.762 | 0.839 | 1.26 | 7 |
| 0.05% (S-14) | PL | 5/10 | 3.0 | 10/10 | 0.711 | 0.832 | 1.58 | 7 |
| 0.25% | PE | 10/10 | 3.0 | 10/10 | 0.923 | 0.940 | 1.39 | 6 |
| 0.25% | PL | 10/10 | 0.9 | 10/10 | 0.934 | 0.947 | 1.45 | 7 |
| 1% | PE | 10/10 | 4.0 | 10/10 | 0.973 | 0.971 | 1.82 | 6 |
| 1% | PL | 10/10 | 0.1 | 10/10 | 0.973 | 0.974 | 2.16 | 7 |
| 5% | PE | 10/10 | 4.3 | 10/10 | 0.992 | 0.982 | 2.07 | 5.5 |
| 5% | PL | 10/10 | 0.0 | 10/10 | 0.991 | 0.987 | 3.23 | 6 |

Point-exponential recovers the planted program as often as point-Laplace, and more often
near the threshold (8 against 5 of 10). It converges in every fit and its calibration
slopes are closer to 1. But it marks several other programs survival-active in every
fit: at shares of 1% and above, every point-exponential fit has at least two false-active
programs, while point-Laplace has at most one in any fit. Near the threshold both arms
exceed the false-activity limit (PE 2.6, PL 3.0).

**PL reuse check.** The 40 point-Laplace refits reproduce the saved S-9 and S-14 fits:
the same recovery verdict in 40/40 and validation C equal to three decimals in 40/40.
The false-active count differs in 5 S-14 seeds, by one program each. All five are seeds
where no program matched the planted one (best |r| below 0.7). The S-14 scoring removed
the best-correlated program from the false-active count even when it did not match, and
the S-9 scoring used here removes it only when it matches. The fits are the same, so the
saved fits could have been reused. They were refitted only so that both arms use the same
scoring code.

## Real training data

| Arm | Mean within-cohort C | Difference from PL (SE) | Mean pooled C | Difference from PL (SE) | Folds converged | Survival-active programs per fold | Max \|cor\| of a program with the cohort label | Max \|cor\| with mean expression |
|---|---|---|---|---|---|---|---|---|
| PL | 0.658 | | 0.663 | | 3/5 | 2 to 4 | 0.03 to 0.05 | 0.33 to 0.42 |
| PE | 0.573 | −0.085 (0.031) | 0.502 | −0.161 (0.036) | 0/5 | 0 | 0.998 | 0.88 to 0.91 |
| RANK | 0.615 | −0.043 (0.016) | 0.620 | −0.043 (0.018) | 0/5 | 2 | 1.000 | 0.84 to 0.87 |

Fold-wise values are in `real_fits.csv` and per-program values in `real_programs.csv`.
All seven columns keep nonzero loadings in every arm and fold. The point-Laplace refits
reproduce the saved S-11 fold fits exactly: the same pooled concordance, convergence flag,
number of survival-active programs and programs with nonzero loadings in all five folds.
The two point-Laplace folds that reach 300 sweeps (folds 3 and 4) also did so in S-11.

**What the point-exponential programs encode.** X/sd keeps each cohort's mean level. In
the prepared input the median gene mean is 11.0 in TCGA and 15.5 in CPTAC, against a
within-cohort SD of 1, so the mean level dominates the data. In every fold one
point-exponential program's training scores correlate 0.998 in absolute value with the
TCGA indicator, and every program's correlation is between 0.56 and 0.998. The program
most correlated with each patient's mean expression reaches 0.88 to 0.91.
With pooled fitting and one projection center, the programs describe the difference
between the cohorts. None is survival-active in any fold. The rank transform removes
patient-level means, but each gene's average rank still differs between platforms, so its
programs separate the cohorts in the same way.

**Overlap with the current programs.** Top-100 gene overlap with the current S-11 fit's
programs 2 and 3: the best point-exponential program shares 8 genes with program 2 and 42
with program 3. The best rank-transform program shares 30 genes with program 2 and 37 with
program 3. The point-Laplace full-data fit is the S-11 fit (100/100 for every program).

## Criteria

| Criterion | Result | Pass |
|---|---|---|
| 1. Recovery in every cell ≥ PL's minus 1 of 10 | 8 vs 5, 10 vs 10, 10 vs 10, 10 vs 10 | yes |
| 2. Mean false-active programs ≤ 1 in every cell | 2.6, 3.0, 4.0, 4.3 | **no**, every cell |
| 3. ≥ 9 of 10 converged in every cell | 10/10 in every cell | yes |
| 4. Within-cohort CV C not lower than PL's by more than one SE, all folds converge | −0.085 (SE 0.031); 0/5 converged | **no** |

Point-Laplace stays as the L and F prior. The failure is in false activity in every
simulation cell and in real-data concordance and convergence.

## Limits

- The real-data comparison is pooled by decision (DECISIONS.md 2026-11-12). Point-
  exponential on X/sd could behave differently with cohort strata, which center the
  projection within cohort, or with cohort-specific mean offsets in the reconstruction.
  Neither was tested.
- The simulation has one cohort, so it does not show the cohort-membership problem seen on
  real data. Its false-activity failure is separate from that problem.
- The near-threshold cell fails the false-activity limit for both priors (PL 3.0 per fit).
  That is a property of the working model there, already seen in S-14.
- At the 300-sweep cap, the point-exponential and rank-transform fits are not converged.
  Their reported concordance is from the state at the cap.
- No external cohort was used.

## Independent review

[Reviewed 2026-11-12](../reviews/2026-11-12-task-a-independent-review.md) (`3ff618fc`). The
review recomputed all four criteria from the fit rows, rescored all saved simulation states
and reran two fits from scratch (one PE simulation fit, one PE real fold), which matched
exactly. It agrees with the decision. It adds two points that apply here: keeping
point-Laplace is a comparison result, not evidence that point-Laplace passes an absolute
acceptance gate (its S-14 false activity is 3.0 and 3 of 5 real folds converge), and this is
the registered K = 7 comparison, not a new evaluation at the adopted K = 5. The full-data
point-Laplace fit is now included in `pl_s11_reuse_check.csv`.

## Text for the shared status files

Recorded by this session in DECISIONS.md, PROJECT_STATUS.qmd and the week-plan checklist.
