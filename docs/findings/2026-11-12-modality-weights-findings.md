# Reconstruction weights for expression and methylation

Date: 2026-11-12. Author: Andrew Walther.

## Scope

The registered study is [the modality-weights prompt](../prompts/2026-11-12-codex-d-modality-weights-prompt.md), revised in `9c4cac23`.
The aim is recovery of joint and modality-exclusive programs with a large CpG set.
Prediction is secondary. The study used main, with at most two CPU workers.
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

## Fixed-fit ICGC scoring

The loader starts with 79 matched donors. Twelve have unusable survival, leaving
67. The existing primary-PDAC and primary-tumour filter retains **50 distinct
donors with 30 events**. Each donor contributes one row. No model was refitted,
tuned, reoriented or selected on these outcomes.

The four available full-data K = 5 fits were projected once. Expression uses the
fixed Task C PL full fit, and methylation uses its fixed 10k and 100k PL fits.
The joint model is the existing unweighted grouped 10k fit from the preceding
Task D comparison. No weighted joint fits exist because the gate did not allow
Part 4. The full-data checkpoints are all converged.

Each modality is z-scored within these 50 external donors before aligning to the
training feature order. Missing features receive zero after scaling. Both fits
using expression lack 27 D4 genes in ICGC. No CpGs are missing, and none of the
required available features are constant in this subset. ICGC imputation uses the
existing reader over all available CpGs in the frozen design, before the primary
restriction, matching the loader's imputation convention. This imputation does
not use outcomes. This scaling deliberately follows S-12 rather than Task C's
previous external transform based on training means and SDs.

The intervals are donor percentile-bootstrap 95% intervals with 2,000 draws and
seed 1. All concordances use the fixed Cox direction through
`frozen_reverse_cindex()`. Bootstrap helpers receive explicit FALSE flip flags.

| Fixed fit | ICGC C | 95% interval | Missing genes | Missing CpGs |
|---|---:|---|---:|---:|
| expression_only | 0.660 | 0.552 to 0.764 | 27 | 0 |
| methylation_10k | 0.636 | 0.518 to 0.742 | 0 | 0 |
| methylation_100k | 0.594 | 0.470 to 0.709 | 0 | 0 |
| joint_unweighted_10k | 0.692 | 0.575 to 0.795 | 27 | 0 |

The paired joint-minus-expression difference is **+0.032**
with interval **-0.093 to +0.147**.
The bootstrap resamples the same donor indices for both scores. This is
descriptive and does not establish a predictive gain. The joint fit remains
**not validated by simulation**, and these results change no simulation or
weight-selection decision.

For direction checks, each joint program was matched to the expression-only
program with largest absolute TCGA patient-score correlation. All five matches
exceeded 0.5 in absolute value. The same sign was retained for all five ICGC
score correlations. This describes agreement between molecular projections,
not the sign or significance of an external survival effect.

| Joint program | Expression match | TCGA score correlation | ICGC score correlation | Direction retained |
|---:|---:|---:|---:|---|
| 1 | 1 | 0.665 | 0.382 | TRUE |
| 2 | 1 | 0.715 | 0.713 | TRUE |
| 3 | 4 | 0.716 | 0.507 | TRUE |
| 4 | 3 | -0.599 | -0.831 | TRUE |
| 5 | 1 | 0.722 | 0.761 | TRUE |

Aggregates are in `results/multimodal_real/outputs/modality_weights_icgc_2026-11-12/`.
Input hashes were checked before and after scoring. Fixed fits and prepared Task C
inputs were unchanged. Donor-level projections and risks are saved only in
`data/modality_weights_icgc_2026-11-12/fixed_scores.rds`.

Reproduce the analysis with
`Rscript results/multimodal_real/score_modality_weight_external.R`. It refuses to
rescore an existing output directory. The two focused tests cover ordered alignment,
post-scaling zero fill, constant-feature counts and training-frozen direction
matching. Both passed after failing on the missing implementation.

## Implementation walkthrough

The weight helper takes a numeric training matrix and one modality label per
feature. It returns two named positive numeric weight vectors. The count rule uses
block dimensions. The eigenvalue rule z-scores training rows and uses the largest
eigenvalue of the patient Gram matrix. Its memory requirement is n by n rather
than p by p. It never sees held-out outcomes or observations.

The fitter expands those weights to feature order and multiplies the current τ
only when computing reconstruction coordinates. The returned τ remains a noise
precision. Noise updates retain n / Σ E[R²], while objective evaluation weights
both log τ and residual terms. Priors remain at weight one, which allows an entire
modality-program block to collapse. The grouped and pooled paths differ only in
the scope over which the F prior is fitted. Unsupported restart/pruning controls
fail explicitly instead of receiving partially weighted computations.

The simulation runner uses the same observed matrices, outcomes and signed-SVD
starts within each seed/cell. It fits the four registered arms and writes a row
for every fit, including errors if any occur. It computes recovery from standardized
loading truth, not from validation outcome correlations. Weighted precision uses
posterior second moments. A zero total precision has no defined share and remains
NA. The supplemental summarizer verifies exact seed/arm/program keys before
forming counts, so omitted rows cannot silently improve success rates.

The external alignment helper accepts a donor-by-feature numeric matrix and a
fixed training feature vector. It z-scores the external series, then returns a
matrix in training order with explicit zero columns for absent features. The
scoring runner applies the frozen F, normalization and β once. It passes the
resulting donor risk vectors to the established survival bootstrap helpers.
Matching directions uses TCGA scores only to choose program labels, then reads
the correlation sign on ICGC without changing any scores.

Base R supplies scaling, named alignment, Gram eigenvalues and aggregate checks.
`ebnm` fits the existing empirical-Bayes priors, and `survival` supplies the
established fixed-orientation concordance. No dependencies were added.

## Independent final review

A fresh independent reviewer found no blocking statistical or implementation
errors. The review reproduced all simulation tables, recomputed all weighted
precision and final emptiness rows from the 200 private posterior checkpoints,
verified the planted zero blocks and every Block B matched loading correlation,
and confirmed grouped-prior structure and the explicit pooled path.

For ICGC, an independent comparable-pair calculation reproduced all four C estimates
using 845 comparable pairs. The reviewer reproduced every 2,000-draw interval and
SE, the paired-difference interval and all five direction statements from the
private frozen scores. The final findings and the no-Part-4 decision were checked
against those calculations.

The review identified a weak reconstruction-scale fixture whose posterior means
were nearly zero. Nonzero L and both-modality F norm assertions failed on the old
fixture. An aligned rank-one signal replaced it, after which the focused suite
passed 25/25 and the independent reviewer reproduced that pass. This test covers
a corresponding-coordinate reconstruction sweep, not arbitrary-start or
canonicalized full-fitter equivariance. No scientific fits or scores were changed.

Final verification after the review correction: `Rscript tests/run_tests.R`
passed **872 tests with zero failures**. The existing bounded-β profile-likelihood
fallback diagnostic warning remains expected. The independent reviewer reran the
strengthened focused suite and found no remaining blocking findings.

Pooled fits follow the historical return convention and do not save the numerical
pooled F prior values in their checkpoints. Source review and the all-one bitwise
comparison establish that their computation uses a pooled prior. A future audit
could retain those prior values privately without changing fitting.

The historical summary helper's true-risk-C fields were not computed in this
study and remain NA. Defined held-out concordances are reported with their
denominators, but validation C is not a gate or selection criterion.

## Text for the shared status files

Reconstruction-only modality powers were derived and implemented in the grouped
YFB fitter. Both registered rules failed to recover the primary program in either
50:250 simulation cell, despite convergence of all 200 fits. The weighted fits
emptied most modality-program blocks. Neither rule was selected, so no new TCGA
10k or 100k joint fits ran. All Block B arms had zero full pattern successes,
including grouped and pooled feature-count priors. Separate-prior support was
not established. The multimodal method remains future work.

On 50 primary matched ICGC donors with 30 events, the existing expression, 10k
methylation, 100k methylation and unweighted 10k joint fits gave C = 0.660, 0.636,
0.594 and 0.692, respectively. The paired joint-minus-expression difference was
0.032 with donor-bootstrap interval −0.093 to 0.147. All five TCGA-matched joint
program score directions were retained. These results are descriptive and the
joint fit is not validated by simulation. They choose no weights or models.


## Checkpoint commits

- `e51e03c3`: reviewed weighted derivation and 109 numerical checks.
- `4d2ed95f`: failing fitter acceptance tests.
- `d3a95252`: weighted fitter and approved explicit F-prior pooling option.
- `ebc081ff`: registered simulations, aggregates and no-weight-selected decision.
- `03d79d03`: fixed ICGC scoring, aggregate intervals and reviewed test strengthening.

The findings and review disposition are recorded in the final documentation commit.
All checkpoints were pushed. No protected shared status, week-plan or prelim files
were edited. Proposed status text appears above.
