# Reconstruction weights for expression and methylation

Date: 2026-11-12. Author: Andrew Walther.

## Scope and progress

The registered study is [the modality-weights prompt](../prompts/2026-11-12-codex-d-modality-weights-prompt.md), revised in `9c4cac23`.
The aim is recovery of joint and modality-exclusive programs with a large CpG set.
Prediction is secondary. Work proceeds on main, with at most two CPU workers.
Task C inputs are read-only. Patient fits and scores remain under ignored `data/`.

## Derivation and independent review

The weighted objective multiplies each entire Gaussian reconstruction log-likelihood
by its modality weight. All priors and KL terms retain weight one. L and the
reconstruction part of F use effective precision w_m τ_mj. Survival updates and
Z = Σ_m Y_m F_m are unchanged. The positive weight cancels from the featurewise
τ point estimate, n / Σ_i E[R²_mij]. The objective includes the weighted log τ term.
Common-program rescaling preserves the objective when priors transform with the
posterior. Data-block rescaling changes the relative survival projections.

The updated derivation verification passed 109/109 checks. A fresh independent
statistical review ran the same checks and found no blocking mathematical errors.
A stale sentence was corrected following review. The tests include brute-force
coordinate coefficients, numerical τ maximization, canonicalization and a negative
control for block-scale dependence of Z.

For the revised zero-survival scaling check, corresponding starts and uncanonicalized
reconstruction coordinates are required. Projection normalization and fresh SVD
starts themselves depend on block scale. Positive rescaling before training-fold
z-scoring is tested separately for the supervised fit.

## Fitter implementation and verification

`code/modality_weights_yfb.R` validates named positive powers and computes both
training-only rules. The leading eigenvalue uses the patient Gram matrix, avoiding
a CpG-by-CpG allocation. Constant training features fail explicitly.
`fit_ssbmf_yfb()` sends w_m τ_mj to reconstruction updates while retaining the
unweighted τ estimate and projection. The objective weights both Gaussian terms.
NULL and all-one weights bypass extra arithmetic.

Andrew approved `F_prior_by_group = TRUE` as a separate option on November 12.
FALSE pools the F prior while retaining feature-group labels for the weights.
The all-one pooled path reproduces the ungrouped posterior and history bit for bit.
Weighted pruning, posterior restarts and burn-in/warm starts are rejected because
this study does not implement or evaluate them. Defaults remain unchanged.

The focused run passed 25/25 tests, including the prior feature-group checks.
The tests failed before the new API was implemented. The weighted objective was
checked using numerical integration of loading KL under the returned prior, not
only its cached KL. A separately recomputed objective also matched. Full suite:
865 passed, zero failed. The existing bounded-β profile-likelihood fallback warning
was retained. The revised zero-survival check uses a controlled reconstruction
coordinate sweep with corresponding transformed moments, separately from the
supervised raw-unit preprocessing test.

## Registered simulation results

The 20-fit smoke check and all 200 scientific fits completed. All scientific fits
converged, with zero fitting errors and zero recorded warnings. K = 5, 120 training
and 120 validation patients, seeds 1–10, point-Laplace L/F and per-program Normal β
were fixed. The runner uses one worker and one numerical thread.

`run_modality_weight_simulation.R` reuses the existing training-only scaling,
common signed-SVD starts and survival regeneration with seed + 12000. Block B's
new generator option zeros methylation program 1 and expression program 2 after
loading generation. Program 3 remains joint. Default generator output is unchanged.
The gate applies to primary program 1. False activity excludes unique loading
matches at absolute correlation ≥ 0.7 to either prognostic truth, regardless of
coefficient sign. Wrong signs still fail survival recovery.

### Block A: balance

Each row has 10 seeds. Empty blocks here count final zero posterior second-moment
blocks, summed over all seeds and programs. Expression-only has 50 possible blocks,
and joint arms have 100. Repeated point-mass resets are distinct and retained in
the precision CSV, along with their defined denominator.

| Expression:CpGs | CpG nonzero | Arm | Recovery | Mean false activity | Converged | Final empty blocks | Gate |
|---|---:|---|---:|---:|---:|---:|---|
| 50:50 | 3% | expression_only | 9/10 | 0.2 | 10/10 | 28 | TRUE |
| 50:250 | 3% | expression_only | 5/10 | 0.3 | 10/10 | 33 | FALSE |
| 50:250 | 15% | expression_only | 7/10 | 0.4 | 10/10 | 31 | FALSE |
| 50:50 | 3% | grouped_count | 10/10 | 0.4 | 10/10 | 29 | TRUE |
| 50:250 | 3% | grouped_count | 0/10 | 1.2 | 10/10 | 82 | FALSE |
| 50:250 | 15% | grouped_count | 0/10 | 1.5 | 10/10 | 79 | FALSE |
| 50:50 | 3% | grouped_eigenvalue | 9/10 | 0.0 | 10/10 | 26 | TRUE |
| 50:250 | 3% | grouped_eigenvalue | 0/10 | 1.1 | 10/10 | 85 | FALSE |
| 50:250 | 15% | grouped_eigenvalue | 0/10 | 1.2 | 10/10 | 82 | FALSE |
| 50:50 | 3% | grouped_unweighted | 10/10 | 0.4 | 10/10 | 29 | TRUE |
| 50:250 | 3% | grouped_unweighted | 5/10 | 0.9 | 10/10 | 18 | FALSE |
| 50:250 | 15% | grouped_unweighted | 4/10 | 1.6 | 10/10 | 19 | FALSE |

Both weighted rules passed the balanced control, but neither recovered program 1
in any seed of either 50:250 cell. Neither rule is effective under the registered
criterion. No rule was selected. The weighted fits had 79–85 final empty blocks
out of 100 in the imbalanced cells, versus 18–19 for unweighted grouped fits.
This is loss of molecular programs under the stronger relative prior, even though
the algorithm converged. It is not a numerical convergence failure.

### Block B: joint and modality-exclusive programs

A full pattern success requires all three loading correlations to reach absolute
0.7 and each matched program to have the registered modality pattern. Survival
activity and sign of programs 1 and 2 are reported separately. The absent modality
must be empty or have precision share below 5%. Both shares must reach 5% for the
joint program. This definition does not require activity of nonprognostic program 3.

| Expression:CpGs | Arm | Full pattern successes | Converged | Final empty blocks |
|---|---|---:|---:|---:|
| 50:50 | grouped_count | 0/10 | 10/10 | 32 |
| 50:250 | grouped_count | 0/10 | 10/10 | 84 |
| 50:50 | grouped_eigenvalue | 0/10 | 10/10 | 39 |
| 50:250 | grouped_eigenvalue | 0/10 | 10/10 | 86 |
| 50:50 | grouped_unweighted | 0/10 | 10/10 | 32 |
| 50:250 | grouped_unweighted | 0/10 | 10/10 | 33 |
| 50:50 | pooled_count | 0/10 | 10/10 | 18 |
| 50:250 | pooled_count | 0/10 | 10/10 | 82 |

Grouped and pooled feature-count arms each had zero full pattern successes in both
cells. This study provides no full-pattern evidence that separate priors serve
the intended purpose. Per-truth loading recovery, modality pattern, activity and
sign are retained in `patterns.csv` and the supplementary summaries. Per-modality,
per-program precision shares and emptiness are retained in `precision.csv`.
The final empty-block definition differs from the historical runner's repeated
reset count, so those two totals should not be compared directly.

### Verification and reproduction

The new generator tests failed on the baseline interface and passed 3/3 on the
implementation. Pattern and rule-selection tests passed 2/2. The full suite passed
870 tests with zero failures. The existing bounded-β fallback diagnostic warning
was retained. Aggregate CSVs and source hashes are in
`results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weights/`.
The smoke output uses a separate `_weights_smoke/` directory and cannot select a rule.

```sh
Rscript results/benchmark_sim/run_modality_weight_simulation.R --smoke
Rscript results/benchmark_sim/run_modality_weight_simulation.R
Rscript results/benchmark_sim/summarize_modality_weight_simulation.R
```

The fit runners refuse to overwrite existing study directories. The summarizer
recomputes aggregate exhibits from the fixed scientific CSVs without refitting.
The original source hashes identify the code used for fitting. Supplemental
aggregation verified all 200 metric keys, 240 pattern rows and 1,850 precision rows.
An independent raw-CSV calculation reproduced all supplemental counts and shares.
At 50:250, neither weighted grouped arm recovered the methylation-exclusive or
joint loading truth in any seed. At 50:50, grouped count recovered expression-only,
methylation-only and joint loadings in 5, 2 and 7 seeds, respectively. The
corresponding pooled counts were 3, 2 and 5. Those marginal recoveries did not
coincide with all required modality patterns in a single full success.

## TCGA fits: not run

Part 4 was not run because neither weight rule met the registered effectiveness
criterion. No new TCGA 10k or 100k fits, promoter-pair analyses, CV estimates or
100k methylation-specific-factor recovery claims were made. The existing unweighted
10k joint fit remains descriptive and not validated by simulation.

## Text for the shared status files

Pending the registered simulation, TCGA decision and external scoring.
