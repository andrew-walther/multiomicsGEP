# Task B findings: β prior

Date: 2026-11-12. Author: Andrew Walther.

**Result.** Normal β stays, applying the registered rule literally. Neither sparse prior
meets A's false-activity limit in the S-14 near-threshold cell: shared point-Laplace
averages 3.3 false-active programs per fit there, shared point-normal 3.5, and the limit
is 1. Shared point-Laplace also loses on training CV, with within-cohort C 0.649 against
0.658 (paired difference −0.008, SE 0.002). Removing the S-14 cell from the rule does not
change the verdict. Without it, the sparse priors improve calibration in only one of the
three remaining cells, so part (b) fails. Shared point-normal qualifies only if S-14 is
exempted from the false-activity part while its calibration result is kept. Its CV
difference is −0.002 (SE 0.004). Normal itself averages 3.0 false-active programs in the
S-14 cell, so in that cell the criterion is one that the working model also fails. In practice the four priors give nearly the same fits. The C values
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
| (b) without the S-14 cell | 1 cell (1%), fails | 1 cell (1%), fails |
| (c) within-cohort CV C not lower by more than 1 SE | fails (−0.008, SE 0.002) | passes (−0.002, SE 0.004) |
| Qualifies | no | no |
| Qualifies with S-14 removed from (a) and (b) | no | no |
| Qualifies with S-14 exempted from (a) only | no | yes |

**Is the S-14 criterion binding?** Yes, and the baseline fails it too. Normal averages
3.0 false-active programs in the S-14 cell, so a sparse prior is held to a limit there
that the current model does not meet. Under the literal rule, Normal stays. Dropping the
S-14 cell entirely does not change that, because part (b) then fails: the only remaining
cell where the sparse priors are closer to 1 is the 1% cell, by 0.0003. Point-normal would
replace Normal only under a mixed sensitivity that exempts S-14 from false activity but
keeps its calibration difference of 0.03. Even then, the CV difference is −0.002.
That would not be a substantive improvement. Andrew's view (2026-11-12) is that
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

## Independent review

[Reviewed 2026-11-12](../reviews/2026-11-12-task-b-independent-review.md) (`b4dc665e`). Two
fresh fits matched exactly, and the inclusion probabilities were verified. The review agrees
that Normal β stays. It found two errors, both corrected: the "without S-14" sensitivity had
kept S-14's calibration result (with S-14 removed from both parts, neither sparse prior
qualifies), and the real-data all-zero check had counted survival-active programs instead of
testing the coefficients (no fold has every coefficient zero). `beta_rule.csv` and
`beta_real_paired.csv` were regenerated. Future point-normal checkpoints should keep the
final sweep's prior and pseudo-observations so the inclusion probabilities can be audited
without refitting.

## Text for the shared status files

Recorded by this session in DECISIONS.md, PROJECT_STATUS.qmd, the attempts log and the
week-plan checklist.

## Rerun at K = 5

This real-data rerun uses the adopted working rank from training CV. The K = 7
results above remain the record for the original comparison. Only rule (c) is
reevaluated at K = 5. Simulation criteria (a) and (b) were not rerun, so this
comparison cannot establish that a sparse prior qualifies at K = 5 under the
complete registered rule.

The command is
`OMP_NUM_THREADS=1 Rscript results/benchmark_sim/run_training_cv_redo.R brerun 3 5`.
It fits shared point-Laplace, point-normal and Normal β on `D4_Z`, with folds
0 to 5, using three PSOCK workers and one numerical-library thread per worker.
Fold 0 is the full fit. The per-program Normal reference reuses
`data/training_cv_redo_2026-11-12/rank_D4_Z_K05_f*.rds` without refitting.
All fits retain point-Laplace L/F, canonicalization, unit survival weight,
pooled TCGA and CPTAC without strata, no pruning, and config iteration,
tolerance and α settings. No external outcomes enter the comparison.

### Paired real-data results

All 18 new fits converged, with no fit warnings or errors. All five Normal
reference folds also converged. Every arm had three survival-active programs
in the full fit, and no CV fold set every β to zero.

| β prior at K = 5 | Mean within-cohort C | Difference (paired SE) | Mean pooled C | Pooled difference (paired SE) | Folds converged | Mean survival-active per fold | Rule (c) |
|---|---|---|---|---|---|---|---|
| Per-program Normal, reused | 0.654620 | 0 | 0.660322 | 0 | 5/5 | 2.6 | reference |
| Shared point-Laplace | 0.648599 | −0.006021 (0.004924) | 0.658705 | −0.001617 (0.001812) | 5/5 | 2.0 | fails |
| Shared point-normal | 0.647783 | −0.006837 (0.006578) | 0.657493 | −0.002830 (0.005189) | 5/5 | 3.0 | fails |
| Shared Normal | 0.647783 | −0.006837 (0.006578) | 0.657493 | −0.002830 (0.005189) | 5/5 | 3.0 | fails |

Both sparse priors fail part (c) at K = 5. Point-normal misses its one-SE
boundary by only 0.000259 in C, so its pass/fail verdict is sensitive to small
changes in fold estimates. At K = 7 it passed this CV criterion, although it
failed the complete rule because of the S-14 simulation cell. This rerun
provides no basis for replacing per-program Normal β. The broader registered
recovery, false-activity and convergence gate still applies to scientific
adoption.

Full-fit point-normal inclusion probabilities, in program order, are
0.999999700, 1.000000000, 0.999999920, 0.999999978 and 0.999999994.
Shared point-normal and shared Normal have practically identical concordance
here. These PIPs do not provide a useful sparse program selection rule on
these data.

Aggregate results are in
`results/benchmark_sim/outputs/training_cv_redo_2026-11-12/b_rerun_paired.csv`.
Private fits and scores remain in ignored `data/`. The initial driver wrote
all results but exited with a trailing parser error after its source file
was edited while R was reading it. A subsequent command reused all 18
checkpoints, without refitting, to verify a clean exit.

### Implementation and verification

- `f_fit()` accepts optional `prior_beta` and `beta_prior_mode` arguments. Their
  defaults preserve per-program Normal β. It takes a patient-by-gene matrix,
  outcomes, rank, centering flag and parsed config, and returns the existing
  fitter's result. No fitter or modular update code changes.
- `f_task()` dispatches the three new arm names to shared priors. Existing
  scoring and survival-activity definitions stay in use. It saves each fit's
  metrics, coefficient moments, loadings and held-out risk scores under its own
  `b_rerun_*_D4_Z_K05_f*` key in ignored `data/`. For point-normal it also saves
  `lfp_pip()` values from the final prior trace. The activity threshold and
  posterior inclusion probability are distinct summaries.
- `f_b_summary()` reads six checkpoints per arm and six saved references. It
  matches the five held-out folds, computes each fold's C difference from
  Normal, and reports its mean and sample SD divided by √5. Rule (c) requires
  the within-cohort mean difference to be at least minus its paired SE.
  Full fits supply inclusion probabilities and full-fit activity, never CV C.
  Missing checkpoints, failed fits or nonfinite concordances stop aggregation.
  The current implementation also stops for undefined full-fit point-normal
  PIPs. A flat coefficient likelihood can produce such a PIP in future runs.
- The `brerun` command accepts rank from the command line, or reads
  `benchmark$k_working_yfb` when omitted. It limits workers to three and writes
  `b_rerun_paired.csv`. The original `summary` mode excludes these new arms
  to preserve the existing screen, scaling and rank tables.
- `tests/test_training_cv_brerun.R` uses a captured fitter call and synthetic
  checkpoints to check prior dispatch, paired estimates, the one-SE rule,
  zero-coefficient counts, PIPs and failure on a missing fold. It uses base R
  assertions and the runner's existing yaml dependency. No patient data enter
  the tests.

Independent code review found no blocking issue for this rerun. Its optional
flat-likelihood PIP concern is recorded above. All five PIPs in the current
full fit are defined.

The full suite, `Rscript tests/run_tests.R`, passed **816/816** with zero
failures. The separate K = 5 runner regression script also passed.

## Text for the shared status files (K = 5)

The real-data β-prior comparison was repeated at the adopted K = 5 using
five saved training folds. The 18 new shared-prior fits all converged without
fit warnings. Per-program Normal had within-cohort C 0.654620. Shared
point-Laplace had 0.648599, a paired difference of −0.006021 with SE 0.004924.
Shared point-normal had 0.647783, a difference of −0.006837 with SE 0.006578.
Both fail the registered CV criterion, although point-normal is close to
the boundary. Simulation criteria were not rerun at K = 5. Normal β stays.
The K = 7 results remain labeled and preserved. The methylation/expression
K = 5 rerun is paused pending confirmation that the separate K = 7 study
has finished. Its findings still says “in progress.” No shared status file
was edited by this rerun.
