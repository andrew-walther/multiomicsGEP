# Task B independent review: β prior

Date: 2026-11-12. Author: Andrew Walther.

## Verdict

Retaining per-program Normal β follows the registered four-cell rule. Both sparse
priors fail the S-14 false-activity limit. Shared point-Laplace also fails the paired
training-CV concordance criterion. This is a decision under the registered rule, not
evidence that Normal controls false activity near the detection threshold.

The reported sensitivity claim that point-normal would qualify without S-14 needs
correction. The runner removes S-14 from criterion (a) but still counts its calibration
improvement in criterion (b). Removing that cell from both leaves only one calibration
improvement among the three remaining cells. Neither sparse prior then qualifies under
the requirement of at least two improved cells. The main decision is unchanged.

## Scope and independence

The review used the Task B design and reporting rules in the
[week plan](../plans/2026-11-09-week-plan.md), the prospective B registration in the
[attempts log](../plans/2026-10-16-round2-attempts-log.md), the dated Normal-β correction
in `DECISIONS.md`, the fitter and runner, and the aggregate fit-level CSVs. An independent
statistical reviewer separately recomputed the rule from the raw CSVs without reading
the findings or result entries. Both calculations reached the verdict above.

An initial search inadvertently displayed brief result summaries in the week-plan
checklist and attempts log before the primary review calculations were complete. Thus
the primary review was not fully blinded. The full findings note and decision result
entry were read only after the rule had been independently recomputed. This limitation
does not apply to the separate statistical reviewer.

Work stayed on `main`, without a worktree. The audited runner, `fit_ssbmf_yfb.R` and
`update_beta.R` have no diff from Task B result commit `b0c73b1c`. Other sessions had
unrelated working-tree changes, which were left untouched. This review writes only this
note and its index entry. Rerun checkpoints remain under ignored
`data/task_b_review_2026-11-12/`.

The scope is the registered K = 7 Task B comparison. A separate session added a K = 5
real-data follow-up to the findings note during this review. That follow-up was not
audited here and does not replace the K = 7 evidence checked below.

## Arms, implementation and deviations

`lfp_fit()` in `results/benchmark_sim/run_lf_prior_comparison.R`, lines 82–94, uses
identical point-Laplace L/F priors, z-scored matrices, K = 7, 300-sweep cap, tolerance
10⁻⁵, `canonicalize = TRUE`, `w_S = 1`, no pruning and pooled training without strata.
It changes only the β prior and `beta_prior_mode`. The Normal baseline uses one prior
per program. The sparse arms and descriptive Normal reference share a prior across
programs. Simulations use the same generated matrices and seed rule. Real fits use the
same prepared patients, genes and folds.

This is a disclosed refinement of the week plan's statement that only `prior_beta`
changes. The prospective registration explicitly names the sharing modes and reference
arm. Sharing is reasonable because estimating a point mass and slab scale from one
pseudo-observation gives little information about either parameter. It also matches the
earlier shared point-Laplace comparison. It does mean that a contrast with the working
Normal arm combines prior-family and sharing effects. The shared Normal reference helps
describe those effects but does not turn the experiment into a full factorial comparison.

In shared mode, the fitter holds `g_beta` fixed during each sweep, sends it to
`update_beta_k()` as `g_fixed`, and estimates the next prior from the finite-noise
pseudo-observations after the sweep. See `code/fit_ssbmf_yfb.R`, lines 1103–1118 and
1240–1257, and `code/update_beta.R`, lines 189–192. Exactly flat likelihoods are excluded
from prior estimation and retain their prior posterior. No warm-up or supervised start
was added. No external outcomes enter these fits or this review.

The Normal rows in `beta_sim_fits.csv` and `beta_real_fits.csv` exactly reproduce A's
`PL` rows in `sim_fits.csv` and `real_fits.csv` on every shared column. All 160 expected
simulation rows and 24 real rows are present, with unique design keys and no recorded
errors. All simulation slopes are defined. Missing fits did not affect the result.

## Rule recomputation

Counts and medians below were recalculated from `beta_sim_fits.csv`. The cell order is
S-14 0.05%, S-9 0.25%, S-9 1%, S-9 5%.

| Arm | Recovered out of 10 | Mean false-active | Converged out of 10 |
|---|---|---|---|
| Normal, per program | 5, 10, 10, 10 | 3.0, 0.9, 0.1, 0.0 | 10, 10, 10, 10 |
| Shared point-Laplace | 5, 10, 10, 10 | 3.3, 0.3, 0.0, 0.0 | 10, 10, 10, 10 |
| Shared point-normal | 4, 10, 10, 10 | 3.5, 0.3, 0.0, 0.0 | 10, 10, 10, 10 |
| Shared Normal, descriptive | 4, 10, 10, 10 | 4.1, 1.5, 0.4, 0.0 | 10, 10, 10, 10 |

Both sparse arms meet the recovery margin and convergence requirement in every cell.
Both exceed the mean-false-active limit of 1 in S-14, so part (a) fails. Normal also
fails that limit in S-14. All three S-9 cells pass part (a) for both sparse priors.

For part (b), define improvement as
`|median slope Normal − 1| − |median slope sparse − 1|`. Positive values pass the
cell-specific comparison.

| Cell | Normal median slope | Point-Laplace improvement | Point-normal improvement |
|---|---|---|---|
| S-14 0.05% | 1.578794238 | 0.026049456 | 0.030526039 |
| S-9 0.25% | 1.448775994 | −0.001822938 | −0.005342863 |
| S-9 1% | 2.156120982 | 0.000521926 | 0.000319240 |
| S-9 5% | 3.234021853 | −0.000562609 | −0.001005404 |

Each sparse arm passes in two of four cells, S-14 and S-9 1%, so part (b) passes.
The improvement at 1% is small and both slopes remain far above 1. If S-14 is removed
from this calculation, each arm improves only one of three cells and part (b) fails.

Part (c) was recalculated by pairing each sparse arm with Normal by fold in
`beta_real_fits.csv`. The SE is the sample SD of the five fold differences divided by √5.
Concordance uses within-cohort pairs with the model's fixed hazard orientation.

| Arm | Mean within-cohort C | Difference from Normal | Paired SE | Part (c) |
|---|---|---|---|---|
| Normal, per program | 0.657841724 | | | reference |
| Shared point-Laplace | 0.649391280 | −0.008450445 | 0.002211785 | fail |
| Shared point-normal | 0.656137083 | −0.001704641 | 0.003610039 | pass |
| Shared Normal, descriptive | 0.651052137 | −0.006789588 | 0.004463982 | descriptive |

All arms converge in three of five real CV folds. Folds 3 and 4 reach 300 sweeps.
All full-data fits converge. The B registration requires simulation convergence but
does not explicitly repeat A's all-CV-folds-converged condition. Either interpretation
retains Normal here. The CV values are finite iteration-budget results, and should not
be presented as a comparison of fully converged fold fits. No sparse real fit has all
coefficients zero.

| Rule interpretation | Point-Laplace | Point-normal |
|---|---|---|
| Registered four-cell rule | does not qualify | does not qualify |
| Exempt S-14 only from part (a), retain it in calibration | does not qualify | qualifies |
| Remove S-14 from parts (a) and (b), retain two-improved-cell threshold | does not qualify | does not qualify |

The middle row is what the current `qualifies_without_s14` field computes. It is a
different sensitivity question from removing the S-14 cell.

## Inclusion probability and scratch reruns

The formula in `lfp_pip()`, lines 98–117, is the posterior slab probability for
`g = π₀δ₀ + (1 − π₀)N(0, σ²)`. It uses `prior_before`, which is correct because
`prior_after` is estimated only after the coefficients have been updated. The component
means are zero for the priors initialized and fitted in these arms. Flat likelihoods
receive NA by design rather than a data-informed inclusion probability.

The original compact simulation checkpoints save PIPs but omit `history$prior_trace`.
Real checkpoints likewise save the program PIPs without their generating trace. Thus
the original saved trace requested by the prompt was unavailable. The scratch reruns
also saved the full returned fits locally, and the audit read that saved trace back
from disk. No original checkpoint or aggregate output was overwritten.

For the rerun of S-9 1%, seed 1, shared point-normal, the final sweep used
π₀ = 0.830256496376 and slab SD = 6.58423955535. An independent log-density calculation,
using `plogis(log slab density − log spike density)`, gave the following probabilities.

| Program | PIP |
|---|---|
| 1 | 1.000000000000 |
| 2 | 0.009770579992 |
| 3 | 0.003278899550 |
| 4 | 0.001991833015 |
| 5 | NA, flat likelihood |
| 6 | 0.002034186040 |
| 7 | 0.001725626356 |

The maximum absolute difference from `lfp_pip()` was 3.47 × 10⁻¹⁸. A separate
`ebnm(..., g_init = prior_before, fix_g = TRUE)` call reproduced the saved trace's
posterior means within 1.36 × 10⁻¹⁹ and SDs within 1.74 × 10⁻¹⁸. The initial audit
call requested an unsupported EBNM output option and stopped after the two successful
fits. Repeating only the posterior verification with the default output succeeded.
This was an audit-call error, not a fit error.

Both requested fits ran from scratch, sequentially with one numerical worker and
BLAS/OpenMP thread limits set to 1. They used the registered configuration keys,
`lfp_sim_task()` and `lfp_real_task()` with `PRIV` redirected to a new ignored review
directory. The simulation used seed 43 for fitting, as the runner prescribes, and the
real fold used seed 42. Cache files did not exist in the new directory before execution.

| Rerun | Result | Agreement with original |
|---|---|---|
| S-9 1%, seed 1, shared point-normal | recovered, false-active 0, converged in 29 sweeps, validation C 0.9643312, calibration slope 1.912037 | all metrics except elapsed time match within 10⁻¹⁰ |
| Real fold 1, shared point-normal | 218 training and 55 held-out patients, converged in 283 sweeps, 7 programs with nonzero loadings, 2 survival-active, within-cohort C 0.665859564, pooled C 0.668240850 | all metrics except elapsed time match within 10⁻¹⁰ |

The saved EF, EBeta, EBeta2 and Tau arrays had maximum absolute difference zero for
both reruns. Held-out risk scores had maximum absolute difference zero in the real
fold. Neither fit produced warnings. Elapsed fitting times were 27 and 118 seconds.
The environment was R 4.5.2 on arm64 macOS, ebnm 1.1-42, with Accelerate BLAS.

The original point-normal real-data PIPs range from 0.9386884 to 0.9999986 in the
full-data fit. In folds 1–3 all are effectively 1, while fold 4 includes a value
0.3747013 and fold 5 a value 0.5978277. These outputs support the reported limitation:
PIP above 0.5 would include almost every real-data program and does not reproduce the
three-program activity set. No threshold replacement is justified by this experiment.

## Errors and comparison with the findings

1. **[P2] The S-14 sensitivity result retains S-14 in calibration.**
   `results/benchmark_sim/run_lf_prior_comparison.R`, lines 573–584, subsets part (a)
   but uses all four `slope_closer` values in `qualifies_without_s14`. Consequently
   `beta_rule.csv` marks point-normal as qualifying without S-14. The findings note's
   opening paragraph, rule table and interpretation repeat that claim, as do the Task B
   decision and B-result entries. A true cell-removal calculation must also subset the
   calibration comparison. Alternatively, label this result explicitly as exempting
   S-14 from part (a) only. This correction changes the sensitivity conclusion, not the
   registered four-cell decision.

2. **[P3] No active coefficient is not the same as every coefficient being zero.**
   The runner's line 563 sets `folds_all_beta_zero = sum(x$n_active == 0)`, then lines
   579–584 use it for the all-zero failure gate. Activity requires the 1.96-SD criterion
   and nonzero loadings. Nonzero but uncertain coefficients can therefore be counted as
   all-zero. Use the saved `max_abs_beta` and an explicit zero tolerance for the
   registered collapse check, and report inactive programs separately. Every current
   sparse fold has active programs and nonzero coefficients, so this error does not
   change this experiment's result.

The findings accurately report the primary decision, simulation counts, paired CV
differences, shared-mode design, real-fold convergence limitation and nearly saturated
real-data inclusion probabilities. The substantive disagreement is the interpretation
of the S-14 sensitivity. Retaining Normal does not establish that the false-activity
acceptance gate is met. The registered recovery, false-activity and convergence gate
remains relevant to any scientific adoption claim.

## Proposed follow-up changes, not made in this review

- Correct or explicitly relabel the S-14 sensitivity calculation, regenerate its
  aggregate rule output, and correct the corresponding findings, DECISIONS,
  PROJECT_STATUS, week-plan checklist and B-result text. Keep the Normal decision.
- Correct the all-zero gate to test coefficients rather than activity counts.
- Preserve the final sweep's `prior_before`, `x` and `s` in future point-normal
  checkpoints so inclusion probabilities can be audited without refitting.
- Continue to report the two nonconverged real CV folds. Any longer-fit study or
  revised acceptance criterion needs its own prospective design.
