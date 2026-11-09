# Single-modality pruning after the KL fix (attempt KL-1)

Issue #6, item 1. Registration: [`2026-11-06-kl-fix-pruning-rerun-plan.md`](../plans/2026-11-06-kl-fix-pruning-rerun-plan.md),
commit 0aad350. Outputs: `results/benchmark_sim/outputs/kl_fix_pruning_rerun_2026-11-06/`
(`results.csv`, `summary.csv`, one CSV per fit in `cells/`). Runner:
`results/benchmark_sim/run_kl_fix_pruning_rerun.R`, run from a `git archive` of 0aad350 with
4 PSOCK workers.

## Summary

- All 160 fits completed with no errors or warnings, and every returned fit converged under the
  ELBO rule with no intermediate caps.
- Without pruning, the selected canonical fit recovers the planted program in 17/20 and 19/20
  seeds, the same counts as on October 15.
- Ordinary and protected pruning give identical fits in every seed. Each recovers the planted
  program in 11/20 seeds at both noise levels, against 2/20 for protected pruning on October 15.
  Both fail the recovery criterion (≥ 16/20) and the rank criterion (median retained K is 7 and
  6.5, not 3).
- The October 15 conclusion therefore changes in degree, not in direction. Pruning still loses the
  planted program whenever it removes anything: 14 of the 15 fits that had recovered the program
  and then pruned lost it. It loses it less often only because, after the fix, it removes nothing
  in 23 of 40 fits.
- Whether pruning starts is often decided by a numerical defect. Programs with all-zero loadings
  get a deletion gain of about ±10⁻⁵ instead of 0, because the β KL of an empty program is
  computed in a regime where ebnm's point-Laplace solver is inaccurate. A failing test records
  this (`tests/test_empty_program_kl.R`). The fitter is not changed.

## Results

Starting K = 7, 20 seeds per row. False-active is the mean number of credibly survival-active
programs other than the planted one.

| Arm | Noise | Recovered | False-active (ref + 0.5) | ELBO-converged | Median K | Mean validation C |
|---|---:|---:|---:|---:|---:|---:|
| reference (α_F = 0) | 0.25 | 1 | 1.15 | 20 | 7 | 0.643 |
| unpruned | 0.25 | **17** | 1.00 (1.65) | 20 | 7 | 0.764 |
| ordinary pruning | 0.25 | 11 | 0.95 (1.65) | 20 | 7 | n/a |
| protected pruning | 0.25 | 11 | 0.95 (1.65) | 20 | 7 | n/a |
| reference (α_F = 0) | 1 | 1 | 0.40 | 20 | 7 | 0.541 |
| unpruned | 1 | **19** | 0.20 (0.90) | 20 | 7 | n/a |
| ordinary pruning | 1 | 11 | 0.20 (0.90) | 20 | 6.5 | n/a |
| protected pruning | 1 | 11 | 0.20 (0.90) | 20 | 6.5 | n/a |

Mean validation C is n/a where at least one fit has an undefined C (constant risk). Per-fit
values are in `results.csv`.

Registered criteria: the unpruned arm passes recovery, false-active and convergence (rank does not
apply). Both pruning arms pass false-active and convergence and fail recovery and rank. October 15
comparison (Task1-3, K = 7 rows, old code): protected pruning recovered 2/20 at both noise levels
with median K 2 and 1; the reference recovered 1/20 at both levels, as here.

## What pruning did

The first deletion is the program with the largest nullcheck gain, if that gain is ≥ 0. The gains
were read from the unpruned fit of the same seed, which is the first stage of the pruning fits.

| First deletion | Fits (noise 0.25; 1) | Planted program recovered before → after |
|---|---|---|
| none (every gain < 0) | 13; 10 | 21 → 21 (2 not recovered before or after) |
| a nonempty program with gain > 0 | 3; 6 | 7 → 0 (2 not recovered before) |
| an empty program whose residue gain was ≥ 0 | 4; 4 | 8 → 1 |

Two things follow.

1. **Loss during the refit is still the mechanism.** After any deletion, the refit starts from the
   retained means with reset variances, β and noise precisions, as Task1-4 traced. Protection does
   not help, because the deleted program is never the protected planted program. The program is
   lost in the refit, not deleted. Task1-5's full-state refit
   (`prune_preserve_posterior = TRUE`) targeted this and was not rerun here.
2. **Whether pruning starts is decided by rounding.** In the 23 fits with no deletion, the fit
   kept seven columns while only two or three programs had nonzero loadings. Empty programs are not
   removed because their deletion gain is slightly negative.

## The numerical defect

For a program with all-zero scores, `update_beta_k()` floors the β precision at 10⁻¹⁰, so the EBNM
noise SD is 10⁵ and the likelihood for βₖ is flat. Then q(βₖ) equals the shared prior and
KL(q‖g) = 0 up to terms of order 10⁻¹⁰. ebnm 1.1.42's point-Laplace solver is inaccurate there. For
the fitted shared prior (π₀ = 0.5, Laplace scale 0.363, prior E[β²] = 0.132), it returns a posterior
second moment of 65535.87 inside the fit (0 in the standalone test) and a log-likelihood about
3 × 10⁻⁶ off. At s = 1000 its posterior SD is 23.7, against a prior SD of 0.36. `compute_ebnm_kl()`
then returns −KL = +8.2 × 10⁻⁶ for each empty program in seed 20261016 (it cannot be positive),
and its sign and size vary by seed. In 11 of the 40 canonical unpruned fits the empty programs'
residue gain was ≥ 0.

The KL fix itself (31db892) is correct. Before it, an empty program's L-update KL carried a bias of
order n·log(A/2π) with A = 10⁻¹⁰. That made empty programs cheap to delete, so pruning always
started.

`tests/test_empty_program_kl.R` (standalone, `Rscript tests/test_empty_program_kl.R`) asserts that
an empty program's posterior equals the prior and its −KL is 0. It fails 0/2 on current code and is
not added to `tests/run_tests.R`, so the shared suite still passes while a fix is chosen.

## Limits

- Starting K = 7 only. The original panel also used 5, 10 and 15.
- Current code differs from the October 15 code in more than the KL fix (for example
  `sign_correction` was removed). The comparison with Task1-3 is "after the fix, on current code".
- The empty-program defect affects which fits begin pruning. Until it is fixed, the pruning arms'
  recovery and rank counts are not a clean estimate of the corrected rule. The finding that a
  deletion followed by a mean-only refit loses the planted program does not depend on it: it holds
  for deletions of nonempty programs (7 of 7).

## Open decisions for Andrew

1. How to treat an empty program: set its β KL to exactly 0 when the precision is at the floor
   (q = g), drop programs with all-zero loadings before the nullcheck, or make the nullcheck treat
   them as removable. Any of these is a fitter change and needs the failing test to pass first.
2. After that fix, whether to rerun KL-1, and whether to include the full-state refit arm
   (`prune_preserve_posterior = TRUE`). That arm is the one that addresses the refit loss.
3. Until then, the K = 7 no-pruning fallback stands. It recovers 17/20 and 19/20 on current code.
