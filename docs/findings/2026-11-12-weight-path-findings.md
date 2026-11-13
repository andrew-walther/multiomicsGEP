# Fixed methylation likelihood-weight path

Registered: 2026-11-12. Completed: 2026-11-13. Author: Andrew Walther.

No reduced methylation weight is usable under the registered rule. All 150 simulation fits converged without errors or warnings. The repeated weights 1 and 0.2 reproduce the earlier aggregates exactly. The four approved TCGA fits are descriptive and are not validated by simulation. Task D stops after this report.

## Scope and fixed rules

The [approved prompt](../prompts/2026-11-12-codex-d-weight-path-prompt.md) fixes methylation weights 1, 0.7, 1/√5, 0.3 and 0.2, with expression weight 1. The three cells are the two Block A 50:250 cells with 15% and 3% methylation nonzeros, and the Block B 50:250 cell with 15%. Each has seeds 1–10. The separate smoke study uses seed 1 in each cell and arm, 15 fits in total. It is excluded from every scientific count below.

The existing generator, survival regeneration, training-fold feature z-scores, common signed-SVD starts, K = 5, point-Laplace L/F priors, separate F priors per modality, per-program Normal β and w_S = 1 are unchanged. No pruning is used. Runs are sequential with one numerical thread. The settings and source hashes were written before fitting. All patient-level and synthetic posterior checkpoints remain under ignored `data/`. Task C input hashes are checked without modifying those inputs.

A reduced weight is usable only if both Block A recovery counts are at least their unweighted counts, with mean false activity ≤1 in both cells. Its matched primary methylation block must have precision share ≥5% in at least 8/10 seeds in both cells. Its Block B expression-only and methylation-only loading-plus-pattern counts must each be at least their unweighted counts. Finally, it must improve either summed Block A primary recovery or summed Block B loading plus pattern. Ties select the largest Block B sum, then the larger weight. These rules were recorded before any fit. Convergence is reported separately and is not an additional clause in this prompt's usable-weight rule.

Recovery retains the earlier full-loading |r| ≥0.7, posterior coefficient activity and correct-sign definition. False activity exempts unique threshold-qualified matches to the two prognostic truths, regardless of coefficient sign. Block B loading-plus-pattern recovery does not require coefficient activity. A joint pattern requires both precision shares ≥5%. An exclusive pattern requires the absent block empty or below 5%.

## Reproducibility and completeness

Both the smoke and full studies passed. For every repeated weight-1 and weight-0.2 fit, the new CSV aggregates exactly match the earlier run's metric outcomes, precision rows and Block B pattern rows. Metrics include recovery, secondary recovery, false activity, convergence, iterations, final emptiness, reset events, full-pattern success, validation C, warnings and errors. Wall-clock timings are excluded. Comparisons use both runs' CSV-serialized values at tolerance zero, which avoids treating CSV rounding as a model change. The full study contains all 150 distinct registered fit identities, 1,500 precision rows, 900 truth-by-modality matched rows and 150 Block B pattern rows. The summaries check exact table identities, missing rows and duplicates before selection.

## Simulation results

Every row below contains 10 seeds and 10 converged fits. Final empty counts have denominator 50 blocks per modality: 10 seeds × 5 programs. They count returned all-zero second moments, not repeated point-mass reset events. Raw signed correlations and every fitted-program precision are retained in the aggregate files.

### Cell summaries and matched primary methylation precision

The methylation share is the mean for the program independently matched to truth 1 by the largest absolute full-loading correlation, even if that correlation fails recovery. A mean is over defined shares, with its denominator shown. Retained means share ≥5%. Empty or below 5% is shown separately because a tiny nonzero block is not an exactly empty block. Block B truth 1 is expression-only, so methylation retention there is not a usability target.

| Block | CpG nonzero | w | Primary recovery | Mean false | Converged | Mean meth share (n) | Retained | Matched empty/<5% | Unmatched | Empty expr | Empty meth |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 15% | 1 | 4/10 | 1.6 | 10/10 | 0.523840 (10) | 9/10 | 1/10 | 0 | 11/50 | 8/50 |
| A | 15% | 0.7 | 4/10 | 0.7 | 10/10 | 0.207537 (10) | 5/10 | 5/10 | 0 | 33/50 | 42/50 |
| A | 15% | 1/√5 | 1/10 | 1.1 | 10/10 | 0.031047 (9) | 1/10 | 8/10 | 1 | 34/50 | 47/50 |
| A | 15% | 0.3 | 0/10 | 1.3 | 10/10 | 0.000002 (9) | 0/10 | 9/10 | 1 | 32/50 | 48/50 |
| A | 15% | 0.2 | 0/10 | 1.5 | 10/10 | 4.96e-11 (9) | 0/10 | 9/10 | 1 | 31/50 | 48/50 |
| A | 3% | 1 | 5/10 | 0.9 | 10/10 | 0.599620 (10) | 10/10 | 0/10 | 0 | 13/50 | 5/50 |
| A | 3% | 0.7 | 1/10 | 1.0 | 10/10 | 0.091161 (10) | 2/10 | 8/10 | 0 | 36/50 | 46/50 |
| A | 3% | 1/√5 | 0/10 | 1.2 | 10/10 | 0.000005 (10) | 0/10 | 10/10 | 0 | 37/50 | 48/50 |
| A | 3% | 0.3 | 0/10 | 1.2 | 10/10 | 0.000003 (10) | 0/10 | 10/10 | 0 | 35/50 | 48/50 |
| A | 3% | 0.2 | 0/10 | 1.2 | 10/10 | 0.000001 (10) | 0/10 | 10/10 | 0 | 34/50 | 48/50 |
| B | 15% | 1 | 5/10 | 0.7 | 10/10 | 0.363149 (10) | 5/10 | 5/10 | 0 | 21/50 | 12/50 |
| B | 15% | 0.7 | 5/10 | 0.2 | 10/10 | 0.058427 (9) | 1/10 | 8/10 | 1 | 38/50 | 47/50 |
| B | 15% | 1/√5 | 5/10 | 0.2 | 10/10 | 2.98e-11 (9) | 0/10 | 9/10 | 1 | 37/50 | 48/50 |
| B | 15% | 0.3 | 5/10 | 0.1 | 10/10 | 3.56e-13 (10) | 0/10 | 10/10 | 0 | 36/50 | 48/50 |
| B | 15% | 0.2 | 5/10 | 0.2 | 10/10 | 2.64e-11 (10) | 0/10 | 10/10 | 0 | 35/50 | 49/50 |

At w = 0.7, Block A recovery is 4/10 at 15% and 1/10 at 3%, compared with 4/10 and 5/10 unweighted. Mean false activity falls to 0.7 and 1.0, but primary methylation retention is only 5/10 and 2/10. At 1/√5, retention falls to 1/10 and 0/10. At 0.3 and 0.2 it is 0/10 in both cells. Thus no reduced weight meets retention, no-loss and gain requirements together. The path already loses many methylation blocks at 0.7. It does not identify a usable interval between dominance and collapse.

### Block-wise correlations and precision for every matched truth

Each entry reports mean absolute Pearson correlation within a modality, followed by the number of defined seeds. Pearson correlation renormalizes each block separately. It is NA when the truth block is absent or the fitted block is zero or constant. Undefined correlations are not replaced by zero. Mean methylation precision and its defined denominator refer to the same full-loading matched program. If every fitted loading is zero, no full-loading match exists. Those seeds are explicitly unmatched, excluded from empty-count denominators, and counted as not retained. Empty-or-under-5% counts exclude unmatched seeds, whose count is shown separately. All three Block A truths are shown, including the nonprognostic third truth. Block B truths 1, 2 and 3 are expression-only, methylation-only and joint, respectively.

| Block | CpG nonzero | w | Truth | Expr absolute r (n) | Meth absolute r (n) | Meth share (n) | Meth ≥5% | Meth empty/matched | Meth empty/<5% | Unmatched |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 15% | 1 | 1 | 0.727 (10) | 0.565 (9) | 0.523840 (10) | 9/10 | 1/10 | 1/10 | 0 |
| A | 15% | 1 | 2 | 0.561 (10) | 0.311 (10) | 0.543719 (10) | 10/10 | 0/10 | 0/10 | 0 |
| A | 15% | 1 | 3 | 0.559 (9) | 0.342 (9) | 0.702423 (10) | 9/10 | 1/10 | 1/10 | 0 |
| A | 15% | 0.7 | 1 | 0.720 (10) | 0.545 (6) | 0.207537 (10) | 5/10 | 4/10 | 5/10 | 0 |
| A | 15% | 0.7 | 2 | 0.381 (10) | 0.181 (5) | 0.183870 (10) | 4/10 | 5/10 | 6/10 | 0 |
| A | 15% | 0.7 | 3 | 0.415 (10) | 0.213 (7) | 0.263253 (10) | 6/10 | 3/10 | 4/10 | 0 |
| A | 15% | 1/√5 | 1 | 0.720 (9) | 0.564 (2) | 0.031047 (9) | 1/10 | 7/9 | 8/10 | 1 |
| A | 15% | 1/√5 | 2 | 0.454 (9) | 0.154 (1) | 0.000004 (9) | 0/10 | 8/9 | 9/10 | 1 |
| A | 15% | 1/√5 | 3 | 0.346 (9) | 0.080 (2) | 0.031043 (9) | 1/10 | 7/9 | 8/10 | 1 |
| A | 15% | 0.3 | 1 | 0.761 (9) | 0.393 (1) | 0.000002 (9) | 0/10 | 8/9 | 9/10 | 1 |
| A | 15% | 0.3 | 2 | 0.428 (9) | 0.251 (2) | 0.000002 (9) | 0/10 | 7/9 | 9/10 | 1 |
| A | 15% | 0.3 | 3 | 0.349 (9) | NA (0) | 0.000000 (9) | 0/10 | 9/9 | 9/10 | 1 |
| A | 15% | 0.2 | 1 | 0.761 (9) | 0.393 (1) | 4.96e-11 (9) | 0/10 | 8/9 | 9/10 | 1 |
| A | 15% | 0.2 | 2 | 0.428 (9) | 0.112 (2) | 4.97e-11 (9) | 0/10 | 7/9 | 9/10 | 1 |
| A | 15% | 0.2 | 3 | 0.351 (9) | NA (0) | 0.000000 (9) | 0/10 | 9/9 | 9/10 | 1 |
| A | 3% | 1 | 1 | 0.604 (10) | 0.508 (10) | 0.599620 (10) | 10/10 | 0/10 | 0/10 | 0 |
| A | 3% | 1 | 2 | 0.551 (10) | 0.437 (10) | 0.529382 (10) | 9/10 | 0/10 | 1/10 | 0 |
| A | 3% | 1 | 3 | 0.308 (8) | 0.169 (10) | 0.772671 (10) | 9/10 | 0/10 | 1/10 | 0 |
| A | 3% | 0.7 | 1 | 0.609 (10) | 0.399 (2) | 0.091161 (10) | 2/10 | 8/10 | 8/10 | 0 |
| A | 3% | 0.7 | 2 | 0.538 (10) | 0.700 (4) | 0.179944 (10) | 4/10 | 6/10 | 6/10 | 0 |
| A | 3% | 0.7 | 3 | 0.164 (10) | 0.155 (4) | 0.179944 (10) | 4/10 | 6/10 | 6/10 | 0 |
| A | 3% | 1/√5 | 1 | 0.630 (10) | 0.369 (2) | 0.000005 (10) | 0/10 | 8/10 | 10/10 | 0 |
| A | 3% | 1/√5 | 2 | 0.443 (10) | 0.188 (1) | 9.69e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| A | 3% | 1/√5 | 3 | 0.125 (10) | 0.036 (1) | 9.69e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| A | 3% | 0.3 | 1 | 0.631 (10) | 0.319 (2) | 0.000003 (10) | 0/10 | 8/10 | 10/10 | 0 |
| A | 3% | 0.3 | 2 | 0.445 (10) | 0.275 (1) | 4.97e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| A | 3% | 0.3 | 3 | 0.124 (10) | 0.062 (1) | 4.97e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| A | 3% | 0.2 | 1 | 0.632 (10) | 0.301 (2) | 0.000001 (10) | 0/10 | 8/10 | 10/10 | 0 |
| A | 3% | 0.2 | 2 | 0.444 (10) | NA (0) | 0.000000 (10) | 0/10 | 10/10 | 10/10 | 0 |
| A | 3% | 0.2 | 3 | 0.123 (10) | NA (0) | 0.000000 (10) | 0/10 | 10/10 | 10/10 | 0 |
| B | 15% | 1 | 1 | 0.570 (10) | NA (0) | 0.363149 (10) | 5/10 | 5/10 | 5/10 | 0 |
| B | 15% | 1 | 2 | NA (0) | 0.177 (9) | 0.780527 (10) | 9/10 | 1/10 | 1/10 | 0 |
| B | 15% | 1 | 3 | 0.562 (8) | 0.314 (9) | 0.736691 (10) | 9/10 | 1/10 | 1/10 | 0 |
| B | 15% | 0.7 | 1 | 0.602 (9) | NA (0) | 0.058427 (9) | 1/10 | 7/9 | 8/10 | 1 |
| B | 15% | 0.7 | 2 | NA (0) | 0.101 (3) | 0.100614 (9) | 2/10 | 6/9 | 7/10 | 1 |
| B | 15% | 0.7 | 3 | 0.443 (9) | 0.428 (3) | 0.100614 (9) | 2/10 | 6/9 | 7/10 | 1 |
| B | 15% | 1/√5 | 1 | 0.587 (9) | NA (0) | 2.98e-11 (9) | 0/10 | 8/9 | 9/10 | 1 |
| B | 15% | 1/√5 | 2 | NA (0) | 0.243 (1) | 2.98e-11 (9) | 0/10 | 8/9 | 9/10 | 1 |
| B | 15% | 1/√5 | 3 | 0.437 (9) | 0.225 (2) | 4.26e-11 (9) | 0/10 | 7/9 | 9/10 | 1 |
| B | 15% | 0.3 | 1 | 0.536 (10) | NA (0) | 3.56e-13 (10) | 0/10 | 9/10 | 10/10 | 0 |
| B | 15% | 0.3 | 2 | NA (0) | 0.282 (1) | 3.56e-13 (10) | 0/10 | 9/10 | 10/10 | 0 |
| B | 15% | 0.3 | 3 | 0.474 (10) | 0.167 (2) | 4.35e-13 (10) | 0/10 | 8/10 | 10/10 | 0 |
| B | 15% | 0.2 | 1 | 0.533 (10) | NA (0) | 2.64e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| B | 15% | 0.2 | 2 | NA (0) | 0.018 (1) | 2.64e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |
| B | 15% | 0.2 | 3 | 0.465 (10) | 0.347 (1) | 2.64e-11 (10) | 0/10 | 9/10 | 10/10 | 0 |

Small defined methylation correlations at low weights can describe only one or two nonempty fitted blocks. They cannot establish retained methylation signal. The expression-block correlations remain defined for the primary truth even when the methylation block collapses. Full-loading recovery and within-expression recovery are therefore different measurements.

### Block B per-program loading and pattern

Counts are out of 10. Pattern alone is reported even when loading recovery fails, so a correct pattern by itself does not identify the true program. No weight achieves a full three-program success. At w = 0.7, the expression-only loading-plus-pattern count is 5/10, methylation-only is 0/10 and joint is 1/10, equal to the unweighted counts. At smaller weights, the joint count drops to zero.

| w | Truth | Loading | Pattern | Loading + pattern |
|---|---|---|---|---|
| 1 | 1 | 5/10 | 5/10 | 5/10 |
| 1 | 2 | 0/10 | 4/10 | 0/10 |
| 1 | 3 | 1/10 | 7/10 | 1/10 |
| 0.7 | 1 | 5/10 | 8/10 | 5/10 |
| 0.7 | 2 | 0/10 | 0/10 | 0/10 |
| 0.7 | 3 | 1/10 | 2/10 | 1/10 |
| 1/√5 | 1 | 5/10 | 9/10 | 5/10 |
| 1/√5 | 2 | 0/10 | 0/10 | 0/10 |
| 1/√5 | 3 | 0/10 | 0/10 | 0/10 |
| 0.3 | 1 | 5/10 | 10/10 | 5/10 |
| 0.3 | 2 | 0/10 | 0/10 | 0/10 |
| 0.3 | 3 | 0/10 | 0/10 | 0/10 |
| 0.2 | 1 | 5/10 | 10/10 | 5/10 |
| 0.2 | 2 | 0/10 | 0/10 | 0/10 |
| 0.2 | 3 | 0/10 | 0/10 | 0/10 |

### Usable-weight decision

| w | No loss A | Meth stays | No loss B | Gain | B sum | Usable |
|---|---|---|---|---|---|---|
| 0.7 | FALSE | FALSE | TRUE | FALSE | 6 | FALSE |
| 1/√5 | FALSE | FALSE | TRUE | FALSE | 5 | FALSE |
| 0.3 | FALSE | FALSE | TRUE | FALSE | 5 | FALSE |
| 0.2 | FALSE | FALSE | TRUE | FALSE | 5 | FALSE |

**Selected usable weight: none.** TCGA fitting proceeds under the separate descriptive authorization.

## Descriptive TCGA full-data fits

All fits use the same 144 matched patients, 2,218 D4 genes and 10,000 ranked CpGs, z-scored on all 144 patients. The unweighted grouped fit is reused, not refit. Four new fits use 0.7, 1/√5, 0.3 and the exact count ratio 2,218/10,000 = 0.2218. The prompt rounds the last label to 0.22. No cross-validation or external scoring is performed. Every result below is **not validated by simulation**.

### Convergence and modality classes

| w | Reused | Converged | Iterations | Seconds | Live | Empty expr/meth | Expr-only | Meth-only | Joint | New |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | TRUE | TRUE | 210 | 214.0 | 5 | 0/0 | 0 | 1 | 4 | 0 |
| 0.7 | FALSE | TRUE | 181 | 245.4 | 5 | 0/0 | 0 | 1 | 4 | 0 |
| 1/√5 | FALSE | TRUE | 195 | 237.3 | 5 | 0/1 | 1 | 1 | 3 | 0 |
| 0.3 | FALSE | TRUE | 137 | 165.5 | 5 | 0/1 | 2 | 1 | 2 | 0 |
| 0.2218 | FALSE | TRUE | 138 | 157.9 | 5 | 0/2 | 2 | 0 | 3 | 0 |

### Per-program precision and fixed-reference score correlations

The two reference correlations are each the best absolute patient-score correlation against any program of the fixed expression-only or 10k-CpG methylation-only fit. New programs require both values below 0.5. Program indices are local to each fit, so equal column numbers across weights do not establish identity.

| w | Program | Expr share | Meth share | Class | Best expr absolute r | Best meth absolute r | New |
|---|---|---|---|---|---|---|---|
| 1 | 1 | 0.026050 | 0.973950 | methylation-exclusive | 0.665 | 1.000 | FALSE |
| 1 | 2 | 0.091475 | 0.908525 | joint | 0.715 | 0.957 | FALSE |
| 1 | 3 | 0.175166 | 0.824834 | joint | 0.716 | 0.935 | FALSE |
| 1 | 4 | 0.490725 | 0.509275 | joint | 0.599 | 0.729 | FALSE |
| 1 | 5 | 0.713815 | 0.286185 | joint | 0.722 | 0.548 | FALSE |
| 0.7 | 1 | 0.031292 | 0.968708 | methylation-exclusive | 0.655 | 0.999 | FALSE |
| 0.7 | 2 | 0.143511 | 0.856489 | joint | 0.737 | 0.935 | FALSE |
| 0.7 | 3 | 0.252326 | 0.747674 | joint | 0.722 | 0.930 | FALSE |
| 0.7 | 4 | 0.650503 | 0.349497 | joint | 0.637 | 0.679 | FALSE |
| 0.7 | 5 | 0.949393 | 0.050607 | joint | 0.822 | 0.400 | FALSE |
| 1/√5 | 1 | 0.022007 | 0.977993 | methylation-exclusive | 0.622 | 0.996 | FALSE |
| 1/√5 | 2 | 0.278136 | 0.721864 | joint | 0.746 | 0.862 | FALSE |
| 1/√5 | 3 | 0.329057 | 0.670943 | joint | 0.622 | 0.991 | FALSE |
| 1/√5 | 4 | 0.828462 | 0.171538 | joint | 0.663 | 0.646 | FALSE |
| 1/√5 | 5 | 1.000000 | 0.000000 | expression-exclusive | 0.824 | 0.396 | FALSE |
| 0.3 | 1 | 0.031285 | 0.968715 | methylation-exclusive | 0.630 | 0.997 | FALSE |
| 0.3 | 2 | 0.445734 | 0.554266 | joint | 0.738 | 0.873 | FALSE |
| 0.3 | 3 | 0.484641 | 0.515359 | joint | 0.625 | 0.987 | FALSE |
| 0.3 | 4 | 0.999376 | 0.000624 | expression-exclusive | 0.874 | 0.520 | FALSE |
| 0.3 | 5 | 1.000000 | 0.000000 | expression-exclusive | 0.845 | 0.314 | FALSE |
| 0.2218 | 1 | 0.051330 | 0.948670 | joint | 0.641 | 0.999 | FALSE |
| 0.2218 | 2 | 0.626362 | 0.373638 | joint | 0.735 | 0.885 | FALSE |
| 0.2218 | 3 | 0.623006 | 0.376994 | joint | 0.628 | 0.986 | FALSE |
| 0.2218 | 4 | 1.000000 | 0.000000 | expression-exclusive | 0.874 | 0.514 | FALSE |
| 0.2218 | 5 | 1.000000 | 0.000000 | expression-exclusive | 0.870 | 0.345 | FALSE |

The lower weights change the modality composition without producing a program that is weakly correlated with both single-modality references. Program 1 remains methylation-heavy at every weight. At 0.2218 its expression precision share rises just above 5%, so it is classified as joint rather than methylation-exclusive. This boundary change does not establish a new biological program. Its correlation with the methylation-only reference remains 0.999. The expression-exclusive programs at the lower weights have strong expression-only reference correlations. These are full-data descriptions and support no claim of improved prediction or simulation-validated recovery.

### Exclusive and new program top features

The complete top-20 lists for every exclusive or new program are in [tops.csv](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/tops.csv). The table contains feature identifiers, signed loadings, CpG-associated genes, UCSC gene-region labels, island relation and a promoter indicator. The indicator is true when any annotation is TSS1500, TSS200, 5′UTR or first exon. This is an annotation definition, not evidence of regulation. Expression genes have no CpG annotation. For a new joint program, both modality lists are retained.

| w | Program | Modality | Top five features |
|---|---|---|---|
| 1 | 1 | methylation | cg06546806, cg02860732, cg22975913, cg15230781, cg04342955 |
| 0.7 | 1 | methylation | cg14965220, cg13450005, cg02860732, cg22975913, cg15230781 |
| 1/√5 | 1 | methylation | cg26756083, cg19084726, cg14747563, cg09103619, cg14965220 |
| 1/√5 | 5 | expression | BMPR2, ITGA1, HMCN1, RPLP1, RPLP2 |
| 0.3 | 1 | methylation | cg09103619, cg14965220, cg26756083, cg15245095, cg14747563 |
| 0.3 | 4 | expression | LAMB1, AXL, FRMD6, PDGFC, RDX |
| 0.3 | 5 | expression | KIAA1211, TNIK, SPATA13, SH3D19, RASEF |
| 0.2218 | 4 | expression | LAMB1, FRMD6, AXL, PDGFC, RDX |
| 0.2218 | 5 | expression | KIAA1211, TNIK, RASEF, SPATA13, SH3D19 |

## Implementation and verification

`run_modality_weight_simulation.R` retains its original study entry point and adds `--weight-path`. The new path prepares each cell and seed once, uses the same common starts for every fixed weight, and checkpoints each fit before comparing repeated arms. `mwpath_cor()` returns a scale-invariant within-block correlation and preserves undefined cases. `mwpath_reproduce()` compares repeated aggregate outcomes exactly. `mwpath_summarize()` checks complete identities and reports defined denominators. `mwpath_matched_stats()` distinguishes an unmatched all-zero fit from a matched block that is empty, so its empty counts use an explicit matched denominator. `mwpath_select()` applies the four clauses and tie break without TCGA or external outcomes.

`run_modality_weight_path.R` reads the prepared design and fixed reference fits. `mwpath_programs()` uses weighted τ E[F²] precision to classify programs, and patient-score correlations to flag new ones. `mwpath_top()` extracts top absolute loadings and joins the installed 450K annotations. `mwpath_real_run()` fits the four approved weights from common starts, records warnings, hashes read-only inputs and writes aggregate outputs. Inputs are patients-by-features matrices and fitter posterior moments. Outputs are CSV summaries and private posterior checkpoints. Existing `yaml`, `ebnm`, `survival` and the locally installed Illumina annotation package are reused. No dependency is added.

The fitted objective and coordinate updates are unchanged by this task. Weighting still reduces both the modality's influence on L and each feature's reconstruction evidence for its own F loading. Posterior second moments determine precision, so a very small nonzero variance need not be an exactly empty block. Independent truth matching may reuse a fitted program for multiple truths, as in the earlier registered study. Missing correlations remain NA, and each mean reports its defined denominator.

Tests were written before the runners. They initially failed because the required helpers or runner did not exist, then passed after implementation. They test the no-gain case, every loss/retention clause, the larger-weight tie break, rejection of altered repeated outcomes, undefined block correlations, second-moment-based program classes, the two-reference new-program rule, and exclusive CpG annotations. The first full suite passed 878 tests with zero failures. The updated suite passed 879 tests with zero failures after the unmatched-denominator reporting fix. The final full suite after the nuisance-truth reporting addition also passed 879 tests with zero failures. A new test first failed without the helper, then passed with explicit matched, unmatched and correlation denominators. One intervening full-suite run aborted at the existing Bayesian-optimization dependency with an OpenMP shared-memory error. That test file passed separately, and the suite was retried with single-thread settings applied at process startup. No code or scientific settings were changed for that retry. The expected bounded-β profile-likelihood fallback warning is retained. No model-fitting function was changed. Source hashes retain the execution-time versions. Later edits document helpers and clarify unmatched-block summary denominators. Summaries were regenerated from unchanged raw aggregates without refitting. The saved posteriors, raw metric, precision and pattern rows are unchanged. Existing matched rows are unchanged. Analysis-only extraction from saved fits adds Block A nuisance-truth-3 rows, bringing matched reporting to all three truths in every cell. The runner now includes these rows directly in future reproductions. No extra fits were run.

## Independent review

The preliminary code review checked the registered design, matching, usable-weight rule, TCGA scaling and top-feature annotation. It caught live-double versus serialized-CSV comparison before the smoke run. The comparison was corrected to exact new-CSV versus old-CSV equality, without a numerical tolerance. Review suggestions added explicit A-recovery and B-exclusive mutation tests and complete table-key checks.

A fresh independent review recalculated weighted precision and exact emptiness from all 150 saved simulation posteriors, checked all 900 matched full/block correlation rows including the added nuisance truth, and independently reapplied the usable-weight rule. It verified all 60 repeated metric rows, 600 precision rows and 60 Block B pattern rows against the prior study exactly. Both smoke and full key sets, all per-program pattern summaries and every matched/unmatched denominator were checked.

For TCGA, the reviewer independently rebuilt the 144-patient standardized matrix and checked all five arms, including the reused baseline. Weighted precision shares, modality classes, both fixed-reference score correlations, every requested top-20 feature and its annotation, convergence iterations and read-only input hashes agreed. There were no remaining blocking findings. Review prompted the explicit unmatched-denominator reporting and inclusion of Block A's third truth, both resolved from saved fits without additional fitting.

## Rebuild and evidence

Run these from the repository root only when the registered output directories do not already exist. The runners refuse to overwrite results.

```sh
Rscript tests/run_tests.R
Rscript results/benchmark_sim/run_modality_weight_simulation.R --weight-path --smoke
Rscript results/benchmark_sim/run_modality_weight_simulation.R --weight-path
Rscript results/multimodal_real/run_modality_weight_path.R
```

Scientific outputs:

- [Simulation settings](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/settings.txt), [reproducibility](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/reproducibility.txt), [cell summaries](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/summary.csv), [truth-matched summaries](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/summary_matched.csv), [pattern summaries](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/summary_patterns.csv), [selection](../../results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/selection.csv).
- [TCGA settings](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/settings.txt), [convergence and classes](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/summaries.csv), [programs](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/programs.csv), [top features](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/tops.csv), [read-only input hashes](../../results/multimodal_real/outputs/modality_weight_path_2026-11-12/input_hashes.csv).

## Commit ledger

- `048c2662`: tested simulation and descriptive TCGA runners, suite integration and README.
- `6168acad`: separate smoke study and complete 150-fit simulation aggregates.
- `96780380`: four descriptive TCGA fit aggregates, reused baseline and input hashes.
- Final checkpoint: independently reviewed findings, helper documentation, explicit unmatched denominators and complete three-truth reporting.

## Suggested next steps after the pause

The path does not support adopting a fixed reduced reconstruction weight. Weighting only the L update would require a new model or optimization design, its derivation and an approved evaluation. It is a decision for Andrew after the pause. No alternative update, extra grid, additional cohort fit or follow-up study was run here.

## Text for the shared status files

The final fixed methylation likelihood-weight path completed all 150 simulation fits, with all converged and no errors or warnings. The repeated weights 1 and 0.2 exactly reproduced the earlier aggregate outcomes. None of 0.7, 1/√5, 0.3 or 0.2 met the registered usable-weight rule. At 0.7, primary methylation precision remained at least 5% in only 5/10 and 2/10 seeds in the two imbalance cells, below the required 8/10. Lower weights lost more methylation signal. Four separately authorized descriptive 10k-CpG TCGA fits were evaluated along the fixed path. They are not validated by simulation and choose no model. Task D is paused after this report, with no follow-up runs started.
