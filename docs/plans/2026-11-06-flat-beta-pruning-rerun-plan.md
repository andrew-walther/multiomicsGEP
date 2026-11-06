# Repeat KL-1 after handling flat β likelihoods (KL-2)

Question. Does pruning still lose the planted program once the canonical fitter
uses exact β posterior moments and zero β KL for flat likelihoods?

Repeat [KL-1](2026-11-06-kl-fix-pruning-rerun-plan.md) with its existing runner,
seeds 20261016–20261035, noise 0.25 and 1, and starting K = 7. Keep all four
arms and all thresholds, including recovery ≥ 16/20, false-active count at
most reference + 0.5, convergence in 20/20, and median retained K = 3 for
pruning. Report warnings and errors, first-deletion gains and paired losses.
Use current branch code with the exact flat-likelihood update, so other code
changes since KL-1 also limit attribution. Use synthetic data only.

Run from a fixed source snapshot with at most four PSOCK workers.

```sh
Rscript results/benchmark_sim/run_kl_fix_pruning_rerun.R 4 results/benchmark_sim/outputs/flat_beta_pruning_rerun_2026-11-06
```

Do not adopt pruning based on this repeat. The full-state refit is a separate
comparison and is outside this registration.
