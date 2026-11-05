# Attempt S-12 findings: external scoring of the fixed survival-weighted YFB fit

Written 2026-11-05. Registration:
[`2026-11-05-supervised-external-scoring-plan.md`](../plans/2026-11-05-supervised-external-scoring-plan.md).
Script: `results/benchmark_sim/score_yfb_supervised_external.R`. Outputs:
`results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/` (`external_cindex.csv`,
`external_paired_diff.csv`). Each cohort was scored once. No model was refitted or tuned.

## Summary

- **External C of the w_S = 1 fit:** 0.584 (Dijk), 0.537 (Moffitt), 0.656 (PACA-AU array),
  0.681 (PACA-AU seq), 0.644 (Puleo). The unweighted mean is 0.620.
- **Against the two-step fit** (same w_S = 0 loadings, then Cox), the supervised fit scores
  equal or higher in every cohort. The cohort-stratified pooled difference is +0.007 (95%
  bootstrap interval −0.002 to 0.017). This matches the training CV difference of +0.011.
- **Against the historical frozen YFB benchmark** (mean 0.627), which used the earlier pipeline
  (double log on two cohorts, 2,064 genes), the new fit is 0.007 lower on average. That
  comparison is not paired, and every cohort difference is well inside its interval.
- **Calibration on real data is not overconfident.** Calibration slopes range from 0.36 to 1.57.
  Puleo, the largest cohort, gives 0.99 (SE 0.20). The supervised slopes are closer to 1 than the
  two-step slopes in every cohort. The overconfidence seen at 2,000 genes in simulation does not
  appear here, because the real coefficients are small and the Normal prior shrinks them.

## Per-cohort results

| Cohort | n (events) | Missing D4 genes | C, w_S = 1 (95% CI) | C, two-step | Paired difference (95% CI) | Calibration slope, w_S = 1 (SE) | Calibration slope, two-step | Historical frozen YFB C |
|---|---|---|---|---|---|---|---|---|
| Dijk | 90 (81) | 30 | 0.584 (0.513–0.658) | 0.567 | +0.018 (0.000 to 0.036) | 0.55 (0.36) | 0.35 | 0.635 |
| Moffitt_GEO_array | 123 (83) | 71 | 0.537 (0.459–0.611) | 0.531 | +0.007 (−0.014 to 0.026) | 0.36 (0.36) | 0.28 | 0.549 |
| PACA_AU_array | 63 (38) | 13 | 0.656 (0.562–0.751) | 0.656 | 0.000 (−0.026 to 0.025) | 0.78 (0.33) | 0.73 | 0.648 |
| PACA_AU_seq | 52 (31) | 27 | 0.681 (0.579–0.780) | 0.671 | +0.010 (−0.019 to 0.038) | 1.57 (0.43) | 1.26 | 0.657 |
| Puleo_array | 288 (181) | 82 | 0.644 (0.599–0.686) | 0.643 | +0.001 (−0.007 to 0.010) | 0.99 (0.20) | 0.80 | 0.645 |
| Pooled, stratified (5) | | | | | +0.007 (−0.002 to 0.017) | | | |
| Pooled without PACA_AU_seq | | | | | +0.006 (−0.003 to 0.016) | | | |

PACA_AU_array and PACA_AU_seq share donors, so the pooled difference is also given without the
smaller of the two. The historical column is from `results/benchmark_snapshot_2026-10-14/` and is
not a paired comparison.

## Interpretation

Externally, the fixed configuration chosen by training CV performs about as well as the
historical benchmark, and slightly better than unsupervised loadings followed by Cox. The gain is
small and its pooled interval includes zero. On the real data, survival barely changes the
loadings (S-11), so the gain comes from the joint β update with its shrinkage prior. That update
also gives better-calibrated risk scores than the unpenalized two-step Cox fit. Dijk and Moffitt
are the weakest cohorts for every model, as in earlier benchmarks.

Making survival move the loadings remains the open method question. This result is the reference
any stronger-supervision variant has to match or beat.
