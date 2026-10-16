# results/

Analysis runners and their saved outputs. Run every script from the repository
root with `Rscript`. Runners read hyperparameters from `config/globals.yml` and
real data from `PDAC_DATA_ROOT`. Patient-level inputs and fitted `.rds` objects
are git-ignored; aggregate tables and figures are tracked.

## Layout

| Folder | Contents |
|--------|----------|
| `benchmark_sim/` | Single-modality benchmark and follow-up analyses on simulated and PDAC data; outputs in `benchmark_sim/outputs/<analysis>/` |
| `multimodal_real/` | Expression + methylation fits on matched TCGA/ICGC data; outputs in `multimodal_real/outputs/` |
| `multi_cohort_sim/` | Simulation of shared and cohort-specific programs across cohorts |
| `cohort_lmm_sim/` | Synthetic check of the cohort-effect model |
| `figures/`, `tables/` | Per-cohort outputs from the May 2026 analyses |
| `legacy/`, `benchmark_sim/archive/` | Retired runners and outputs, kept for reference and not re-run |

## Main entry points

**Benchmark**

| Script | Purpose |
|--------|---------|
| `benchmark_sim/run_YFB_benchmark.R` | YFB model (η = (YF)β): simulation and real-data benchmark, the recommended configuration |
| `benchmark_sim/run_LB_benchmark.R` | LB model (η = Lβ), same design |
| `benchmark_sim/run_k_init_sweep.R` | Starting-rank consensus across ELBO, BIC, log-likelihood and held-out C |
| `benchmark_sim/run_K_cv.R` | Cross-validated K selection |
| `benchmark_sim/run_desurv_comparison.R` | Preprocessing comparison that chose the DeSurv-style gene selection (configuration D4) |
| `benchmark_sim/run_external_ci_analysis.R` | External-cohort concordance with bootstrap confidence intervals |

**Comparators**

| Script | Purpose |
|--------|---------|
| `benchmark_sim/run_ebmf_cox_external.R` | Two-step comparator: unsupervised EBMF, then Cox on the factor scores |
| `benchmark_sim/run_sbmf_desurv_overlap.R` | Gene overlap between fitted programs and DeSurv's published programs |

**Interpretation**

| Script | Purpose |
|--------|---------|
| `benchmark_sim/run_pathway_enrichment.R` | Gene-set enrichment of fitted programs |
| `benchmark_sim/run_subtype_concordance.R` | Agreement of program scores with published PDAC subtypes |
| `benchmark_sim/make_km_figures.R` | Kaplan–Meier figures for external cohorts |

**Multimodal**

| Script | Purpose |
|--------|---------|
| `multimodal_real/run_multimodal_real_fit.R` | Fit on matched TCGA expression and methylation, validated on ICGC |
| `multimodal_real/summarize_multimodal_real.R` | Tables and figures from the saved fits |
| `benchmark_sim/run_multimodal_pruning_comparison.R` | Simulation comparison of multimodal fitting variants |

The remaining scripts in `benchmark_sim/` and `multimodal_real/` are diagnostic
and sensitivity studies. Each one's design and outcome is recorded in
`DECISIONS.md` or the dated reports in
`docs/reports/`; read the record before re-running one. Most refuse to overwrite
existing outputs.
