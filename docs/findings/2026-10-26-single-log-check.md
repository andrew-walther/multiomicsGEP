# Effect of the double log transform on the single-modality benchmark

Date: 2026-10-26. Runner: `results/benchmark_sim/run_single_log_check.R`. Outputs:
`results/benchmark_sim/outputs/single_log_check/` (fits git-ignored). The benchmark snapshot was
not modified.

## The problem

Three RNA-seq cohorts (TCGA-PAAD, Dijk, PACA-AU RNA-seq) were transformed with log₂(x + 1)
twice. The cohort loader `load_data_internal.R` in the data folder applies it once, and
`preprocess_merged_cohorts()` / `preprocess_desurv_cohort()` apply it again when
`cfg$pdac$platform_log_transform` is `true` for those cohorts. CPTAC and the three microarray
cohorts are logged once. The check sets every platform flag to `false`, so each cohort is logged
exactly once, and changes nothing else.

## Design

Both arms use the benchmark (D4) settings: combined mean/variance rank, top 3,000 genes per
cohort, intersection across TCGA and CPTAC, per-cohort z-scores, K = 7, α_F = 0, point-exponential
priors, seed 42, at most 300 iterations, and the same five stratified folds of the 273 training
patients. The double-log arm is refitted with the same code as a control. The risk direction of
every fit is set on its own training patients. Each full-training fit is scored once on the five
external cohorts; the difference from the frozen benchmark risk scores has a 95% percentile
interval from 2,000 patient-level bootstrap resamples within each cohort. Nothing is selected
from the external results.

## Results

The refitted double-log arm reproduces the frozen benchmark: identical 2,064 genes and identical
external C in every cohort (to four decimals; Puleo differs by 3 × 10⁻⁵). The comparison therefore
isolates the transform.

**Genes.** The single-log screen keeps 2,218 genes, 1,985 of them among the frozen 2,064 (96%).
The share of negative entries after z-scoring is 47% (double) and 49% (single).

**External concordance (frozen benchmark vs single log):**

| Cohort | n | Frozen C | Single-log C | Difference (95% CI) | Risk-score correlation |
|---|---|---|---|---|---|
| Dijk | 90 | 0.634 | 0.630 | −0.004 (−0.012, 0.002) | 0.991 |
| Moffitt array | 123 | 0.549 | 0.541 | −0.008 (−0.018, 0.001) | 0.997 |
| PACA-AU array | 63 | 0.648 | 0.645 | −0.004 (−0.019, 0.012) | 0.998 |
| PACA-AU RNA-seq | 52 | 0.657 | 0.665 | +0.008 (−0.004, 0.024) | 0.997 |
| Puleo array | 288 | 0.645 | 0.645 | +0.001 (−0.003, 0.004) | 0.998 |
| Mean | | 0.627 | 0.625 | −0.001 | |

**Training cross-validation (same folds):** mean C 0.608 (SE 0.046) double, 0.642 (SE 0.029)
single; paired difference +0.033. Three of the five single-log fold fits stopped at the
300-iteration cap without meeting the convergence rule (all five double-log fold fits converged),
so the CV difference compares partly non-converged fits and should not be read as a gain. One
double-log fold (fold 5) gives a nearly constant risk score (SD 3 × 10⁻⁸; zero survival-active
programs) yet C = 0.73 from tiny differences; this shows how unstable the per-fold estimates are
at n ≈ 55 per fold.

## Conclusion

Correcting the double log changes the selected genes slightly and leaves the external performance
of the benchmark unchanged (mean C 0.625 vs 0.627; every per-cohort interval includes zero; risk
scores correlate at 0.99 or more). The published benchmark numbers do not need revision on this
account. The correction should still be made for the rebuilt method, because a double log
compresses the upper range of expression for three cohorts and the two analyses should share one
pipeline. The non-convergence of three single-log fold fits is a property of the current fitter
on these inputs and is one more instance of the convergence problem the rebuild addresses.

For Andrew's decision: adopt single log for all RNA-seq cohorts in new work (recommended), and
leave the frozen benchmark and the interim write-up numbers unchanged with a note.
