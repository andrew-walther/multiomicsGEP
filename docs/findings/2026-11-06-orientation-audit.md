# Risk orientation audit

A fitted Cox risk score is evaluated with `frozen_reverse_cindex()` so larger
risk predicts shorter survival. A below-chance concordance remains below
chance. For cached historical fits, the score direction comes from the
recorded convention or the training data.

Lines refer to the pre-audit scripts. Paths beginning with `run_` in the
first group are under `results/benchmark_sim/`.

| Script | Line | Reporting trace | Action |
|---|---:|---|---|
| `run_cohort_lmm_benchmark.R` | 126 | May 22 cohort-LMM report and DECISIONS | Score fitted risk |
| `run_desurv_comparison.R` | 164 | May 27 alignment report, October writeup and DECISIONS | Score fitted risk, mark new cache direction, rescore saved cache |
| `run_external_ci_analysis.R` | 53 | October writeup and October 2 decision | Read the cache marker when present. Keep the documented historical direction for the old D4 cache |
| `run_joint_bo.R` | 107 | DECISIONS | Score fitted risk |
| `run_k5_vs_k7_bootstrap_ci.R` | 84, 104 | August 21 findings and DECISIONS | Freeze each cached fit's direction on training data |
| `run_k_init_multistart_check.R` | 145 | August 27 and September 4 progress book and DECISIONS | Score fitted risk |
| `run_k_init_multistart_dip_check.R` | 162 | August 21 findings, September 4 progress book and DECISIONS | Score fitted risk |
| `run_k_parsimony_curve.R` | 105 | July 15 report and PROJECT_STATUS | Score fitted risk |
| `run_k_parsimony_followup.R` | 113 | July 15 report and DECISIONS | Score fitted risk |
| `run_merged_benchmark.R` | 118 | May 25 report and DECISIONS | Score fitted risk |
| `run_yfb_beta_fix_diagnostic.R` | 118 | PROJECT_STATUS and DECISIONS | Score fitted risk |
| `run_yfb_dualF_diagnostic_D4.R` | 68 | DECISIONS | Score fitted risk |
| `run_yfb_dualF_frozen_grid_D4.R` | 74 | DECISIONS | Score fitted risk |
| `run_yfb_vs_ebmf_k7_matched_ci.R` | 34, 51 | October and November writeups | Freeze both directions from cache provenance and rescore |
| `results/cohort_lmm_sim/run_synthetic.R` | 219, 240, 256, 274, 297 | May 22 cohort-LMM report and DECISIONS | Score fitted risk |
| `results/multi_cohort_sim/sim_scoring.R` | 146 | August 21 findings, progress book and DECISIONS | Score fitted risk in the shared helper used by the five simulation runners. Preserve their constant-score convention |
| `code/pathway_enrichment.R` | 733 | Signature robustness and negative-control tables | Take apparent C from the within-cohort Cox summary. This describes a local association and is distinct from held-out prediction |
| `code/concordance_ci.R` | 80, 189 | Legacy helper API | Keep optional same-sample orientation for backward compatibility. Audited validation callers pass explicit directions |

The executable complement search found no hits in the legacy or dated benchmark
archive directories. Historical tables and the benchmark snapshot remain
available as provenance.

## Rescoring from saved scores

- [DeSurv alignment rescore](../../results/benchmark_sim/outputs/desurv_comparison/2026-11-06-orientation-rescore.csv)
  compares the historical good-prognosis cache with
  `desurv_comparison_results.csv`. All 25 cells agree at the saved precision.
  The pre-fix D4 sign is documented in the October 2 decision. The historical
  runner generated the other arms using the same post-fit survival direction
  convention.
- [Rank-matched rescore](../../results/benchmark_sim/outputs/ebmf_cox_external/2026-11-06-yfb-vs-ebmf-k7-orientation-rescore.csv)
  reproduces `yfb_vs_ebmf_k7_matched_ci.csv`, including pooled ΔC = 0.0105
  (95% interval −0.0125 to 0.0329). D4 is negated for its historical cache,
  and the EBMF Cox predictor already has the hazard direction.
- The other affected outputs contain concordance summaries or fit objects.
  Rescoring requires the original outcome-aligned scores or reconstruction
  with patient data or regenerated synthetic datasets. No models were refit.

All modified R files parsed. `Rscript tests/run_tests.R` reported 777 passing
and zero failing tests.
