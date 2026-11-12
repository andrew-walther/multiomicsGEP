# Matched-TCGA multimodal comparison

Author: Andrew Walther

Date: 2026-11-12

**not validated by simulation**

## Scope and registration

Andrew authorized this descriptive comparison in the
[next-steps registration](../prompts/2026-11-12-codex-d-next-prompt.md), overriding
the earlier simulation stop. Nothing in this comparison adopts a prior,
program count or model. The separate informative-scenario simulation does not
change this status.

The runner is `results/multimodal_real/run_multimodal_tcga_comparison.R`.
It reads task C's saved `data/methylation_single_modality_2026-11-12/design.rds`
without modifying that preparation or task C's outputs. The fixed design uses
144 matched TCGA patients with 75 events, 2,218 D4 expression genes and the
10,000 highest-ranked CpGs
from the saved preparation and the same five folds. Feature selection remains
the frozen task C selection. Feature means and SDs are fitted on each fold's
training patients only.

Each fold and the full cohort receive expression-only, methylation-only,
joint grouped-prior and joint pooled-prior fits. K = 5, point-Laplace L/F,
per-program Normal β, w_S = 1, canonicalization, no strata and no pruning are
fixed. The iteration limit and tolerance come from `config/globals.yml`.
The frozen fold sizes are 29, 29, 29, 29 and 28. There are 24 distinct fits. The one-fold smoke fits are reused only if the
saved source manifest and input checksums agree exactly. A capped smoke fit
is retained, rather than rerun to obtain convergence.

## Implementation and checks

`mtcga_prepare_fold()` reads the saved folds and returns training/test indices,
training-only standardized blocks and their stored means and SDs. Fold zero
uses all patients for the full-data fit. It calls the existing `fgsim_scale()`
so validation never refits preprocessing. The new tests use an extreme
held-out feature value to distinguish training-only scaling from whole-cohort
scaling and assert the original fold vector exactly.

`mtcga_fit()` initializes a signed SVD on the fitted matrix and delegates to
the existing working fitter through `fgsim_fit()`. The two joint arms use the
same matrix and deterministic SVD. Grouped output must contain five valid
expression priors and five methylation priors. Single-modality arms fit their
own matrices. Each fit returns an aggregate convergence/C row and a private
checkpoint containing moments, final priors, feature order and scores.
Warnings are captured individually. All patient-level state stays under
ignored `data/`. Fold C uses the saved held-out outcomes and greater-hazard
score orientation, with no outcome-based sign flip. Full-cohort C is apparent
training discrimination, not cross-validation.

`mtcga_paired()` subtracts expression C from each joint arm's C within each of
the five folds. The reported difference is the mean of those five differences,
with SE = SD(differences)/√5. An undefined fold makes both paired summaries
undefined. Missing or undefined folds are not silently excluded.

`mtcga_joint_only()` compares each joint program's expression and methylation
loading blocks with the corresponding single-modality programs. A live joint
program qualifies only if both modality comparisons are defined and neither
has a maximum absolute loading correlation of at least 0.5. Dead single-arm
columns contribute undefined correlations and are excluded from the maximum,
with the number of defined comparisons reported. If all comparisons in a
block are undefined, novelty is unclassified. Dead joint columns do not
qualify. The regression checks include a dead single-modality column and an
entirely undefined comparison block. Top features are ranked by absolute
loading, preserving the signed loading in the aggregate table.

`mtcga_summarize()` writes fit rows, individual warnings, fold means/paired
differences, live-program modality precision shares, grouped empty-block
counts, full-data fitted priors, program matches and top-20 joint-only
features. Precision uses Σ_j τ_mj E[F²_mjk]. Empty-block counts include repeated
refitted point-mass resets and are not counts of distinct programs. The prior
table reports point-mass weight 1 − π and Laplace scale 1/rate. A pure point
mass has no defined slab scale. Pooled priors repeat across modality rows and
are explicitly identified by their arm. Program indices are fit-specific.

Every aggregate CSV includes the exact label “not validated by simulation”.
The runner records package versions, source hashes, the configuration and
before/after checksums of the task C preparation. It verifies input and source
hashes after every fit. One sequential worker is used, with BLAS threads set
to one. The synthetic tests were committed before implementation in
`7a8346ce` and initially failed for the missing runner functions. The tested
runners and smoke aggregates are committed in `c2bde939`. The full suite
passes **856 tests with zero failures**, with its existing bounded-β
profile-likelihood fallback warning. No fitter/update kernel changed.

## Results

All **24 distinct fits** are complete. There are **21 converged fits and
three retained caps**: expression folds 1, 3 and 4, each at 300 iterations.
No fit errors, study warnings or empty-block reset events occurred. The four
fold-1 smoke checkpoints were reused exactly after matching input/source
hashes. Every concordance is defined. The final source checksums match the
recorded source manifest, and task C's preparation has identical before/after
checksums. The aggregate results are committed in `078249b1`.

The [aggregate directory](../../results/multimodal_real/outputs/multimodal_tcga_comparison_2026-11-12/)
contains fit rows, fold means and paired differences, precision shares,
empty-block counts, full priors, program matches, the empty joint-only
feature table, warning headers, package versions, hashes and settings.
Patient-level inputs, fits and scores remain under ignored `data/`.

Five-fold TCGA concordance: **not validated by simulation**.

| Fold | Expression | Methylation | Joint grouped | Joint pooled |
|---|---:|---:|---:|---:|
| 1 | 0.695 | 0.611 | 0.655 | 0.637 |
| 2 | 0.685 | 0.708 | 0.625 | 0.713 |
| 3 | 0.552 | 0.629 | 0.564 | 0.544 |
| 4 | 0.594 | 0.473 | 0.455 | 0.464 |
| 5 | 0.430 | 0.306 | 0.358 | 0.321 |
| Mean | 0.591 | 0.545 | 0.531 | 0.536 |

Paired five-fold differences against expression: **not validated by simulation**.

| Joint arm | Mean C difference | SE of difference | Defined paired folds |
|---|---:|---:|---:|
| Grouped | -0.060 | 0.024 | 5 |
| Pooled | -0.055 | 0.030 | 5 |

The mean differences are computed within the same five saved folds and then
averaged. All five differences contribute to each SE. The grouped joint
mean is 0.060 below expression, with SE 0.024. The pooled joint mean is 0.055
below expression, with SE 0.030. These are descriptive internal-CV results,
not a model-adoption decision or an external-validation interval.

Convergence and iterations: **not validated by simulation**. Entries give the iteration count. An asterisk marks the registered 300-iteration cap without convergence.

| Fit | Expression | Methylation | Joint grouped | Joint pooled |
|---|---:|---:|---:|---:|
| Full cohort | 212 | 112 | 210 | 269 |
| Fold 1 | 300* | 97 | 138 | 156 |
| Fold 2 | 156 | 92 | 135 | 140 |
| Fold 3 | 300* | 109 | 235 | 232 |
| Fold 4 | 300* | 111 | 145 | 141 |
| Fold 5 | 94 | 112 | 141 | 146 |

The three capped expression fits were not rerun to obtain convergence.
These comparisons describe the registered fixed-budget fits, rather than
fully converged optima for every arm. The individual caps remain in the
means and paired differences. Every methylation and joint fit converged.
All 60 grouped modality/program reset-count rows are zero. The full-data
results below use apparent training C and are distinct from fold-held-out C.



Full-cohort fit performance: **not validated by simulation**.

| Arm | Apparent C | Status | Iterations |
|---|---:|---|---:|
| Expression | 0.628 | converged | 212 |
| Methylation | 0.594 | converged | 112 |
| Joint grouped | 0.586 | converged | 210 |
| Joint pooled | 0.568 | converged | 269 |

The following shares partition each program’s L precision between expression and methylation. They are reconstruction-precision contributions, not shares of survival prediction or explained molecular variance. CpGs constitute 10,000/12,218 = 81.85% of the features. Program numbers are fit-specific and are not aligned across the two joint fits. The complete CSV retains both shares for every live program in all six fits of each joint arm.

Full-cohort modality precision shares: **not validated by simulation**.

| Program | Grouped expression | Grouped methylation | Pooled expression | Pooled methylation |
|---:|---:|---:|---:|---:|
| 1 | 0.026 | 0.974 | 0.033 | 0.967 |
| 2 | 0.091 | 0.909 | 0.094 | 0.906 |
| 3 | 0.175 | 0.825 | 0.155 | 0.845 |
| 4 | 0.491 | 0.509 | 0.544 | 0.456 |
| 5 | 0.714 | 0.286 | 0.390 | 0.610 |

Prior entries are reported as point-mass weight / Laplace scale. The pooled prior is shared across modalities. Equal program indices across fits do not imply matched programs.

Full-cohort fitted F priors: **not validated by simulation**.

| Program | Grouped expression | Grouped methylation | Pooled |
|---:|---:|---:|---:|
| 1 | 4.86e-07 / 5.39e-05 | 1.83e-07 / 0.00015 | 2.41e-07 / 0.000133 |
| 2 | 6.57e-07 / 0.000193 | 4.86e-07 / 0.000223 | 5.12e-07 / 0.00023 |
| 3 | 8.88e-07 / 0.00029 | 2.81e-06 / 0.000216 | 2.49e-06 / 0.000225 |
| 4 | 6.59e-07 / 0.00062 | 1.09e-06 / 0.000247 | 0.118 / 0.000358 |
| 5 | 9.18e-07 / 0.000761 | 3.59e-06 / 0.000214 | 1.23e-05 / 0.000352 |

All ten joint programs were live. For each, all five single-modality program correlations were defined in each block. Every joint program matched at least one program in each single-modality fit at |correlation| ≥ 0.5. Consequently, no program met the joint-only rule, and the top-20 feature table contains a header and no data rows. These matches do not establish biological equivalence.

Full-cohort maximum absolute loading correlations: **not validated by simulation**.

| Joint arm | Program | Expression | Methylation | Joint-only |
|---|---:|---:|---:|---|
| Joint grouped | 1 | 0.834 | 0.997 | No |
| Joint grouped | 2 | 0.780 | 0.960 | No |
| Joint grouped | 3 | 0.822 | 0.972 | No |
| Joint grouped | 4 | 0.622 | 0.815 | No |
| Joint grouped | 5 | 0.694 | 0.763 | No |
| Joint pooled | 1 | 0.877 | 0.999 | No |
| Joint pooled | 2 | 0.822 | 0.977 | No |
| Joint pooled | 3 | 0.844 | 0.981 | No |
| Joint pooled | 4 | 0.664 | 0.917 | No |
| Joint pooled | 5 | 0.701 | 0.894 | No |

## Interpretation and limits

The two joint fits do not improve mean held-out TCGA concordance over
expression in this fixed analysis, and neither produces a program that
meets the joint-only loading-correlation rule. The fitted group priors allow
different modality scales, but these descriptive results do not show added
prognostic or interpretive information. Matching programs by loading
correlation is a structural summary, not proof of biological equivalence.

The cohort contains 144 patients and 75 events. The five CV training sets
overlap, and the reported SE follows the registered five-difference formula.
The feature sets are the frozen task C selection. Only feature means and
SDs are re-estimated within folds. This is not a new nested feature-screen
analysis. Three expression fits did not reach the ELBO tolerance at the
iteration cap, limiting a comparison of fully converged estimators.
Grouped and pooled fits also retain different loading-KL accounting, so
their absolute monitored objectives must not be compared as though only the
prior grouping changed. No external scoring, weights, retuning or model
adoption occurred. This comparison remains **not validated by simulation**
even though one smaller synthetic cell passes the expression-only gate.

## Independent review

Read-only code review found no blocking implementation or statistical issue.
A fresh reviewer then independently reconstructed all 24 fitted score
vectors from the frozen raw input using each fold's training-only scaling.
It verified their saved scores and C, the fold means, paired SEs, convergence
and caps, every joint precision share, the empty-block counts, full fitted
priors, all ten live/defined loading matches and the empty top-feature table.
It confirmed every aggregate CSV's validation label, all source and input
hashes, and exact MD5 equality of the four reused smoke checkpoints.
The final code, aggregate and prose review returned **PASS**, with no
blocking code, statistical or reporting defect.

## Text for the shared status files

Task D's matched-TCGA comparison is complete and **not validated by simulation**. It reuses task C's 144 patients, 75 events, five folds,
2,218 D4 genes and 10,000 ranked CpGs read-only, with training-fold z-scores
and K = 5. Mean CV C is 0.591 for expression, 0.545 for methylation, 0.531
for joint grouped priors and 0.536 for joint pooled priors. Joint minus
expression is −0.060 with SE 0.024 for grouped and −0.055 with SE 0.030 for
pooled. All five fold differences contribute. There are 21 converged fits
and three expression caps, in folds 1, 3 and 4 at 300 iterations. The caps
were retained rather than rerun. No errors, study warnings or empty-block
resets occurred. Neither full-data joint fit contains a joint-only program
under the registered loading-correlation rule, so no top-20 joint-only
feature list is generated. Different modality F priors are fitted, but this
analysis does not establish added prognostic or interpretive information.
No model or program is adopted. The full suite passes 856 tests with zero
failures, with its existing bounded-β fallback warning. Independent code,
score and aggregate review found no blocking defect.
