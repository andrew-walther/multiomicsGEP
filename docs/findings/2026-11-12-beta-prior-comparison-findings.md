# Task B findings: β prior

Date: 2026-11-12. Author: Andrew Walther.

**Result.** Normal β stays, applying the registered rule literally. Neither sparse prior
meets A's false-activity limit in the S-14 near-threshold cell: shared point-Laplace
averages 3.3 false-active programs per fit there, shared point-normal 3.5, and the limit
is 1. Shared point-Laplace also loses on training CV, with within-cohort C 0.649 against
0.658 (paired difference −0.008, SE 0.002). Shared point-normal would qualify if the S-14
cell were left out. Its CV difference is −0.002 (SE 0.004). But Normal itself averages 3.0
false-active programs in that cell, so in this cell the criterion is one that the working
model also fails. In practice the four priors give nearly the same fits. The C values
differ in the third decimal, and the shared point-normal prior keeps most programs
(inclusion probabilities near 1 on real data).

Registration: the 2026-11-12 "B" entry in the
[attempts log](../plans/2026-10-16-round2-attempts-log.md), task B of the
[week plan](../plans/2026-11-09-week-plan.md). Runner:
`results/benchmark_sim/run_lf_prior_comparison.R` (modes `beta_smoke`, `beta_grid`,
`beta_summary`). Aggregate outputs: `beta_*.csv` in
`results/benchmark_sim/outputs/lf_prior_2026-11-09/`.

## Design

All arms use A's result: point-Laplace L and F on per-cohort z-scores,
`canonicalize = TRUE`, `w_S = 1`, K = 7, no pruning, pooled training, no strata. Only the
β prior differs.

| Arm | β prior |
|---|---|
| Normal (working model) | Normal, one prior per program (A's PL fits, reused) |
| Shared point-Laplace | one point-Laplace prior shared across the 7 programs |
| Shared point-normal | one point-normal prior shared across the 7 programs |
| Shared Normal (reference) | one Normal prior shared across the 7 programs |

The sparse priors are shared because a per-program prior is fitted to one
pseudo-observation, and a two-parameter point-mass prior cannot be estimated from one
value. The October 15 comparison also used the shared point-Laplace prior. The shared
Normal arm shows how much of any difference comes from sharing rather than sparsity.

Simulation cells and real folds are A's. 138 new fits (120 simulation, 18 real) returned
with no errors and no fit warnings.

## Simulation

| Cell | Prior | Recovered | Mean false-active | Converged | Median validation C | Median calibration slope |
|---|---|---|---|---|---|---|
| 0.05% (S-14) | Normal | 5/10 | 3.0 | 10/10 | 0.832 | 1.579 |
| | shared point-Laplace | 5/10 | 3.3 | 10/10 | 0.834 | 1.553 |
| | shared point-normal | 4/10 | 3.5 | 10/10 | 0.830 | 1.548 |
| | shared Normal | 4/10 | 4.1 | 10/10 | 0.830 | 1.557 |
| 0.25% | Normal | 10/10 | 0.9 | 10/10 | 0.947 | 1.449 |
| | shared point-Laplace | 10/10 | 0.3 | 10/10 | 0.947 | 1.451 |
| | shared point-normal | 10/10 | 0.3 | 10/10 | 0.947 | 1.454 |
| | shared Normal | 10/10 | 1.5 | 10/10 | 0.947 | 1.438 |
| 1% | Normal | 10/10 | 0.1 | 10/10 | 0.974 | 2.156 |
| | shared point-Laplace | 10/10 | 0.0 | 10/10 | 0.974 | 2.156 |
| | shared point-normal | 10/10 | 0.0 | 10/10 | 0.974 | 2.156 |
| | shared Normal | 10/10 | 0.4 | 10/10 | 0.974 | 2.182 |
| 5% | Normal | 10/10 | 0.0 | 10/10 | 0.987 | 3.234 |
| | shared point-Laplace | 10/10 | 0.0 | 10/10 | 0.987 | 3.235 |
| | shared point-normal | 10/10 | 0.0 | 10/10 | 0.987 | 3.235 |
| | shared Normal | 10/10 | 0.0 | 10/10 | 0.987 | 2.965 |

The sparse priors reduce false activity a little where the planted program is clear
(0.3 against 0.9 at the 0.25% share), and sharing a Normal prior increases it. Near the
threshold all four priors have three to four false-active programs per fit. Calibration
slopes differ between priors by at most 0.03 except in the 5% cell. Both sparse priors
are "closer to 1" than Normal in the S-14 and 1% cells, but at 1% the difference is in the
fourth decimal (2.1556 and 2.1558 against 2.1561). No fit set every β to zero.

**Point-normal inclusion probabilities, simulation.** Mean number of programs with
posterior inclusion probability above 0.5: 4.2 (S-14), 1.3 (0.25%), 1.0 (1%) and 1.0 (5%).
Where one planted program is clear, the inclusion probability picks it out alone.

## Real training data

| Prior | Mean within-cohort C | Difference from Normal (SE) | Mean pooled C | Folds converged | Survival-active programs per fold (mean) |
|---|---|---|---|---|---|
| Normal | 0.658 | | 0.663 | 3/5 | 2.6 |
| shared point-Laplace | 0.649 | −0.008 (0.002) | 0.657 | 3/5 | 2.4 |
| shared point-normal | 0.656 | −0.002 (0.004) | 0.661 | 3/5 | 2.4 |
| shared Normal | 0.651 | −0.007 (0.004) | 0.656 | 3/5 | 2.6 |

All arms reach the 300-sweep cap in the same two folds (3 and 4) as the working model.
In the full-data fit, the shared point-normal prior gives inclusion probabilities from
0.94 to 1.00 for all seven programs. Three programs are survival-active by the 1.96-SD
rule: programs 2 (β −0.37), 3 (0.22) and 5 (0.22). Across folds the inclusion
probabilities are 1.00 for every program in folds 1 to 3, and 0.37 to 1.00 in folds 4
and 5. With n = 273 the shared prior's point mass at zero is small, so the prior acts
close to a Normal prior. Inclusion probabilities above 0.5 would mark nearly every program
as included. They do not reproduce the three-program survival-active set.

## Rule

| Part | Shared point-Laplace | Shared point-normal |
|---|---|---|
| (a) A's recovery, false-activity and convergence criteria in every cell | fails S-14 (false-active 3.3) | fails S-14 (false-active 3.5; recovery 4 vs 5 is within the margin) |
| (a) without the S-14 cell | passes | passes |
| (b) calibration slope closer to 1 in ≥ 2 cells | 2 cells (S-14, 1%) | 2 cells (S-14, 1%) |
| (c) within-cohort CV C not lower by more than 1 SE | fails (−0.008, SE 0.002) | passes (−0.002, SE 0.004) |
| Qualifies | no | no |
| Would qualify without S-14 | no | yes |

**Is the S-14 criterion binding?** Yes, and the baseline fails it too. Normal averages
3.0 false-active programs in the S-14 cell, so a sparse prior is held to a limit there
that the current model does not meet. Under the literal rule, Normal stays. Had the S-14
cell been dropped, point-normal would have replaced Normal on the strength of a
calibration criterion decided by differences of 0.03 and 0.0003 and a CV difference of
−0.002. That would not be a substantive improvement. Andrew's view (2026-11-12) is that
small false-active coefficients are acceptable. Either way, the practical difference
between the priors is negligible.

## Consequences

- The working model keeps per-program Normal β. Task C's methylation fits, which use
  Normal β, need no rerun on account of B.
- Task F and E proceed with point-Laplace L and F and Normal β.
- The ROADMAP item on replacing the 1.96-SD activity rule with point-normal inclusion
  probabilities is not supported on real data. The inclusion probabilities there are
  near 1 for almost every program.

## Limits

- The sparse priors were tested only in shared mode. Per-program sparse priors are not
  estimable from one pseudo-observation.
- No warm-up or supervised start was used, matching the working-model call. The October
  15 shared point-Laplace configuration used a supervised start and a 10-sweep warm-up.
- The calibration criterion is decided by very small differences.
- No external cohort was used.

## Text for the shared status files

Recorded by this session in DECISIONS.md, PROJECT_STATUS.qmd, the attempts log and the
week-plan checklist.
