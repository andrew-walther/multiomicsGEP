# results/

Analysis runners and their saved outputs. Run every script from the repository
root with `Rscript`. Runners read hyperparameters from `config/globals.yml` and
real data from `PDAC_DATA_ROOT`. Raw inputs and real-data fits are generally git-ignored. The tracked results
include aggregate tables, figures, synthetic checkpoints, and seven legacy
individual risk-score CSVs. Confirm redistribution permissions for the legacy
exports before reusing them in a new release.

The [November 11 inventory](../docs/findings/2026-11-11-results-inventory-findings.md)
maps result scripts, saved locations and report references. It records historical
path mismatches and unresolved writers. The results-rebuild design will use this
inventory before any folder move or runner consolidation.

## Layout

| Folder | Contents |
|--------|----------|
| `benchmark_sim/` | Single-modality benchmark and follow-up analyses on simulated and PDAC data; outputs in `benchmark_sim/outputs/<analysis>/` |
| `multimodal_real/` | Expression + methylation fits on matched TCGA/ICGC data; outputs in `multimodal_real/outputs/` |
| `multi_cohort_sim/` | Simulation of shared and cohort-specific programs across cohorts |
| `cohort_lmm_sim/` | Synthetic check of the cohort-effect model |
| `benchmark_snapshot_2026-10-14/` | Read-only copy of the approved benchmark state; see its README and checksum manifest |
| `figures/`, `tables/` | Per-cohort outputs from the May 2026 analyses |
| `legacy/`, `benchmark_sim/archive/` | Retired runners and outputs, kept for reference and not re-run |

## Main entry points

**Benchmark**

| Script | Purpose |
|--------|---------|
| `benchmark_sim/run_YFB_benchmark.R` | YFB model (η = (YF)β): historical simulation and real-data benchmark |
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
`docs/plans/2026-10-16-round2-attempts-log.md`, `DECISIONS.md` or the dated reports in
`docs/findings/`; read the record before re-running one. Most refuse to overwrite
existing outputs.

`.slurm` files submit the larger runners as job arrays on UNC Longleaf; see
`longleaf_setup/README.md`.

## Expression cohort reference

Every cohort is bulk tumour gene expression, measured by RNA-seq or microarray.
Patient data are not distributed with this repository. Scripts read the PDAC
cohorts from the directory in `PDAC_DATA_ROOT`. Sample sizes and censoring are
for the patients used in the analysis.

| Cohort | Expression platform | n | Censored | Role |
|--------|---------------------|---|----------|------|
| TCGA-PAAD | RNA-seq | 144 | 48% | Training |
| CPTAC-PDAC | RNA-seq | 129 | 50% | Training |
| Dijk | RNA-seq | 90 | 10% | External validation |
| Moffitt | Microarray | 123 | 33% | External validation |
| PACA-AU (array) | Microarray | 63 | 40% | External validation |
| PACA-AU (RNA-seq) | RNA-seq | 52 | 40% | External validation |
| Puleo | Microarray | 288 | 37% | External validation |

The multimodal analysis uses matched TCGA and ICGC expression and methylation,
assembled by `code/load_multiomics_data.R`.

These counts describe the historical analysis records, not unique donors.
The two PACA-AU series overlap. Current donor accounting and prior development
uses are in the [external-cohort audit](../docs/findings/2026-11-06-external-cohort-uses.md).

## Saved-result rebuilds

Run from the repository root. Read the linked findings and registrations before
replaying a study. These commands preserve earlier results and do not authorize
new fits, tuning choices or external model selection.

### Historical benchmark and workflow figures

```sh
Rscript results/benchmark_sim/run_YFB_benchmark.R
Rscript results/benchmark_sim/run_LB_benchmark.R
```

The benchmark runners write to `results/benchmark_sim/outputs/` and read
`config/globals.yml`. The frozen benchmark used `alpha_F = 0`, per-platform
z-scores and 2,064 genes. It reported two prognostic programs and mean external
C = 0.627 on 616 records across five series, including repeated PACA-AU donors.
Its methods, KM figure and interpretation are retained in the
[October writeup](../paper/writeup/2026-10-16-writeup/README.md) and
[project status](../PROJECT_STATUS.qmd). Its estimates are historical.

The October writeup README retains the aggregate-only workflow-diagram builder
and rendering commands. Its [accuracy dispositions](../paper/writeup/2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md)
record the applied reporting fixes and unresolved donor-dependence, preprocessing
and acquisition-provenance questions.

### Historical supervision diagnostics

The experimental normalized-projection fitter did not pass its corrected
scientific gate and remains unadopted. See the
[execution report](../docs/reviews/2026-10-30-codex-execution-report.md) and
[problem statement](../docs/findings/2026-10-30-yfb-supervision-problem-statement.md)
before another study. Its registered commands are retained here:

```sh
Rscript results/benchmark_sim/run_yfb_supervision_decision.R 4 10 attempt7
Rscript results/benchmark_sim/summarize_yfb_supervision.R attempt1 attempt2 attempt3 attempt4 attempt5 attempt6 attempt7
Rscript results/benchmark_sim/run_yfb_saved_fit_diagnostics.R attempt7
Rscript results/benchmark_sim/diagnose_yfb_saved_geometry.R attempt7
Rscript results/benchmark_sim/plot_yfb_saved_traces.R attempt7
```

Use only the registered code revision. Existing checkpoints are reused.
The corrected seventh simulation uses the same training-SD-transformed loading
for reconstruction and survival. Earlier attempts remain numerical/design
diagnostics. `simulation_units.csv` records raw and prepared shares and the
validation loading direction is fixed from training. Saved-fit replays stop if
loading means differ from the registered checkpoints. Keep the patient-level
training cache from the precision audit under ignored `data/`.

The separate [targeted screen](../docs/plans/2026-10-30-yfb-targeted-diagnostics-plan.md)
is reproduced with:

```sh
Rscript results/benchmark_sim/run_yfb_targeted_diagnostics.R
Rscript results/benchmark_sim/run_yfb_beta_uncertainty_screen.R
```

Existing checkpoints are retained and the coefficient screen refuses overwrite.
Its [findings](../docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md)
distinguish fixed-noise stabilization from recovery, false activity and convergence.
The exact coefficient reference is a fixed Gaussian working model, not Cox coverage.

### Progress-book exhibits

The [current-state chapter](../docs/progress_book/current-state.qmd) distinguishes
historical, experimental and current model states and maintains the Completed
Results checklist. Its architecture and fitting/prediction diagrams include a
comparison of two-step and joint fitting. Rebuild the SVGs and render the book:

```sh
python3 docs/progress_book/figures/build_diagrams.py
quarto render docs/progress_book
```

The October methods record is included in the October 16 chapter. The November 6
chapter's saved-fit comparison and simulation figures are rebuilt with:

```sh
Rscript docs/progress_book/figures/verify_october_yfb.R
Rscript docs/progress_book/figures/plot_planted_yfb.R
Rscript docs/progress_book/figures/plot_planted_yfb_comparators.R
```

The saved-fit comparison needs local PDAC training data and historical fits.
The simulation builders read saved replicate tables. Patient-level inputs and
outputs stay out of Git.

### External orientation and donor intervals

```sh
Rscript results/benchmark_sim/run_yfb_vs_ebmf_k7_matched_ci.R
Rscript results/benchmark_sim/score_yfb_supervised_external.R --donor-only
```

The first command reads local historical risk caches and writes aggregate
rescoring beside the original table. The second uses the saved S-11 fits and
local PDAC data (`PDAC_DATA_ROOT` or the configured path) to write
`external_paired_diff_donor.csv` beside the record-based output. Keep every
PACA-AU donor's records together. The [orientation audit](../docs/findings/2026-11-06-orientation-audit.md)
explains historical score directions, fixed before held-out scoring. The
[external-cohort audit](../docs/findings/2026-11-06-external-cohort-uses.md)
records donor linkage and historical development uses.

### Flat coefficient likelihood regression

```sh
Rscript tests/test_empty_program_kl.R
```

For shared Normal, point-normal and point-Laplace β priors, exactly flat
coefficient likelihoods return prior moments and zero β KL, and provide no
information to the prior refit. The check is also in the full suite. Its
[pruning rerun registration](../docs/plans/2026-11-06-flat-beta-pruning-rerun-plan.md)
locates the synthetic evaluation.

### Fixed real-data comparators and program analyses

The [comparator findings](../docs/findings/2026-11-06-real-data-comparators-findings.md#reproduce-and-trace-the-results)
retain the two `run_prompt6_*` commands, pinned DeSurv requirements, local
checkpoint paths and aggregate table definitions. The
[program-analysis findings](../docs/findings/2026-11-06-desurv-style-analyses-findings.md)
retain the four `run_prompt7_*` commands for training-defined risk groups,
external HR/enrichment tables, KM/forest plots, EBMF rank comparisons and the
simulation audit. These commands require the documented local inputs and fits.
The [November writeup README](../paper/writeup/2026-11-06-writeup/README.md)
retains rendering instructions for its aggregate-based synthesis and fixed-cutpoint curves.
