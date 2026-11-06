# Rerun of the October 15 single-modality pruning cells after the KL fix (KL-1)

Registered 2026-11-06, before any registered fit. Attempt label KL-1 in
[`2026-10-16-round2-attempts-log.md`](2026-10-16-round2-attempts-log.md).

## Why

Issue #6, item 1. Commit 31db892 fixed `compute_ebnm_kl()`: with a scalar
precision (the L update) the Gaussian log-normalizer entered the KL once instead
of once per patient. That KL feeds `elbo_full` and `neg_kl`, which the canonical
single-modality fit uses to stop and which `prune_yf_nullcheck()` uses to decide
which program to drop. The October 15 conclusion that ELBO pruning removes the
planted program in single modality rests on those quantities:

- protected pruning recovered the planted program in 2 of 20 seeds at every
  starting K and both noise levels (Task1-3, `outputs/k_selection/round2/attempt3/`);
- ordinary pruning removed it on the two datasets where it was tried
  (DECISIONS.md, 2026-10-15, open item 2).

None of these fits was rerun after the fix. The question is whether, with the
fix, the canonical single-modality fit on the October 15 test bed recovers the
planted program with pruning off, with ordinary pruning and with protected pruning.

## Design

Test bed, unchanged: `make_testbed()` in `results/benchmark_sim/run_supervision_diagnostics.R`
(low-variance prognostic scenario, expression block only, n_train = 400,
n_validation = 120, p = 50, planted log hazard ratio 3.3, cohorts z-standardized),
seeds 20261016–20261035, noise 0.25 and 1. Starting K = 7 only. The original
panel also used K_init 5, 10 and 15; K = 7 is the configuration used for every
downstream fit, and the question is whether pruning keeps the planted program,
which K = 7 answers. If a pruning arm passes every criterion at K = 7, the other
starting ranks would be needed before any adoption, and that would be a new
registration.

Arms, with settings copied from the Task1-3 runner:

| Arm | Settings |
|---|---|
| reference | `alpha_F = 0` (false-positive reference, as in Task1-3) |
| unpruned | `canonicalize = TRUE, w_S = 1, init_supervised = TRUE, beta_prior_mode = "shared", prior_beta = "point_laplace", beta_warmup = 10`, no pruning |
| ordinary | unpruned + `prune = TRUE` |
| protected | ordinary + `prune_protect_survival = TRUE` |

All arms: point-Laplace L and F priors, `max_iter = 300`, `tol = 1e-5` (from
`config/globals.yml`), `converge_on = "elbo"`, starting K = 7. 4 arms × 2 noise
levels × 20 seeds = 160 fits. No fitter default or code changes.

Code: `results/benchmark_sim/run_kl_fix_pruning_rerun.R`, run from a `git archive`
snapshot of the registration commit (other sessions are editing the fitter in the
working tree), with 4 PSOCK workers. Outputs:
`results/benchmark_sim/outputs/kl_fix_pruning_rerun_2026-11-06/`.

## Criteria (unchanged from Task1-1)

Per arm and noise level:

1. recovery ≥ 16/20 (best |cor| with the planted loading ≥ 0.7 and the matched
   coefficient has the true sign);
2. mean false-active programs ≤ reference + 0.5;
3. 20/20 returned fits converged under the ELBO rule;
4. median retained K = 3 (pruning arms only).

Errors count as failures. Intermediate pruning-refit caps are reported
separately, as in Task1-3.

Diagnostics recorded per fit, not criteria: the number of programs with nonzero
loadings (‖E[F_k]‖² > 1e-8), the squared loading norm of each program, and the
nullcheck's deletion gain for each program of the returned fit. The one-seed
smoke test (seed 20261016, outside the registered results) showed both pruning
arms keeping K = 7 with three nonzero programs and a largest deletion gain of
about −1e-6 to −1e-5, so these diagnostics are needed to read the rank result.

## How the result is read

- The primary result is recovery in the ordinary and protected arms. If either
  recovers ≥ 16/20 at both noise levels, the October 15 conclusion that pruning
  removes the planted program does not hold after the fix. If both stay near
  2/20, the conclusion stands.
- The unpruned arm checks that the supervised result (17/20 and 19/20 on
  October 15) still holds on current code.
- The comparison is with the Task1-3 K = 7 rows. The fitter has changed in other
  ways since October 15 (for example `sign_correction` was removed). Differences
  are reported as "after the fix, on current code", not attributed to the KL fix
  alone.
- No default, threshold or K is adopted from this run. A passing pruning arm
  would be reported to Andrew with the K_init panel as the next registration.
