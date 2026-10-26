# Comparator inventory (M7 preparation, 2026-10-15)

> **Historical comparator inventory, not current execution status.**
> Round2 subsequently completed the five-arm simulation panel and Andrew
> deferred the real-data DeSurv head-to-head. Use the
> [current handoff](2026-10-16-round2-return-note.md) for results and authorization;
> retain this survey for provenance rather than restarting its listed analyses.

Read-only survey of the existing comparator analyses for the single-modality (YFB) and multimodal
models. Nothing was rerun. The survey was compiled by a delegated search and spot-checked against
the files cited.

## Single modality, real data (train TCGA + CPTAC; 5 external cohorts)

| Comparator | Script | Mean external C | Output |
|---|---|---|---|
| YFB, K = 7 (α_F = 0) | `results/benchmark_sim/run_desurv_comparison.R` (D4) | 0.627 | `outputs/cohort_beta_comparison/cohort_beta_comparison_results.csv` |
| EBMF → Cox, K = 7 | `run_ebmf_cox_external.R --k 7` | 0.623 | `outputs/ebmf_cox_external/ebmf_cox_external_results_k7.csv` |
| EBMF → Cox, K = 20 / 40 | `run_ebmf_cox_external.R` | 0.581 / 0.578 | `..._results.csv`, `..._results_k40.csv` |
| EBMF (K = 40) → lasso-Cox | `run_ebmf_cox_regularized.R` | 0.601 | `ebmf_cox_regularized_results.csv` |
| YFB − EBMF → Cox, matched K = 7 | `run_yfb_vs_ebmf_k7_matched_ci.R` | +0.011 (−0.013, 0.033) | `yfb_vs_ebmf_k7_matched_ci.csv` |

- Unsupervised EBMF alone has no C-index role; it appears only in factor and pathway
  correspondence (`run_ebmf_factor_correspondence.R`, `run_ebmf_pathway_concordance.R`).
- NMF is not used anywhere in active code.

## DeSurv

No head-to-head run exists. `run_desurv_comparison.R` fits only in-house models (LB/YFB) under
DeSurv-aligned preprocessing (lines 254–267). `run_sbmf_desurv_overlap.R` compares gene lists. See
also DECISIONS.md ("not a head-to-head DeSurv benchmark") and ROADMAP.md.

## Simulation

- `results/multi_cohort_sim/`: YFB vs PCA → Cox vs EBMF → Cox (`survival_strength_sweep_results.csv`,
  `k7_signal_sweep_results.csv`, `multicohort_sim_results.csv`).
- `results/benchmark_sim/benchmark_helpers.R`'s generator has no comparator runs.
- The planted low-variance test bed (`run_supervision_diagnostics.R`) has α_F = 0 (supervision off)
  and a survival-only Cox lasso/ridge (`run_supervision_ceiling.R`), but no EBMF → Cox yet.

## Multimodal

- Real (TCGA → ICGC): `results/multimodal_real/run_multimodal_real_fit.R` arms `two_step_ebmf_cox`
  and `expression_only_yfb`. In `real_sweep_summary.csv` the two-step C is identical (0.418) at
  starting ranks 6, 7 and 8. This suggests the stage-1 rank is not varying; check before reuse.
- Simulation: `code/run_multimodal_yfb_simulation.R` (two-step, stacked-risk, single-modality
  arms); `outputs/multimodal_yfb_smoke/stage4_smoke.csv`.
- No DeSurv, NMF or PCA comparator exists for the multimodal model.

## Gaps for M7

1. DeSurv head-to-head on the same training data and external cohorts (single modality).
2. EBMF → Cox and unsupervised EBMF on the planted test bed (does supervision find what two-step
   misses?).
3. Multimodal comparators rerun with the final method; investigate the constant two-step C.
