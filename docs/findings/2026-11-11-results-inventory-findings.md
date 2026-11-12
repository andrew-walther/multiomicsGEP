# Result scripts and saved-result locations

Date: 2026-11-11. Author: Andrew Walther.

This inventory identifies the scripts and saved evidence that a results rebuild
would need to account for. It covers the nine top-level result folders, including
historical diagnostics and the frozen benchmark, plus exhibit builders outside
`results/`. No analysis was rerun, and no result folder was moved or renamed.
The separate Gemini source-folder deletion is recorded under item 00.

## Inventory files and coverage

| Inventory | Coverage | How to use it |
|---|---|---|
| [Script table](2026-11-11-results-scripts-inventory.csv) | 204 tracked R, Python, shell and SLURM files | One row per script, with purpose, role, inputs, output roots, provisional retention role, document references and source-line evidence. |
| [Folder table](2026-11-11-results-folders-inventory.csv) | 251 result/exhibit locations and the ignored-data boundary | Direct tracked-file counts and bytes, ignored-file metadata counts where inspected, and candidate writers. Includes 12 ignored-only directories under `results/`. |
| [Report-source table](2026-11-11-results-report-consumers.csv) | All 82 tracked `.qmd` and `.Rmd` sources | Literal result/figure references and include/read expressions. Report output formats remain controlled by source YAML and enclosing projects. |

The [November 11 migration record](2026-11-11-early-benchmark-migration-findings.md)
maps the subsequent move of 267 early artifacts to indexed history. This dated
inventory preserves the original paths at its source snapshot.

The source snapshot is `47690aa`, after item 00 was closed. The script table
covers all 155 tracked R/SLURM files under `results/`, including 113 in
`results/benchmark_sim/` and 42 in its sibling folders. It adds 49 executable
or callable scripts outside `results/`, including standalone fitter examples,
reference implementations, figure builders, derivation checks and publication
helpers. Helpers and tests inside result folders are included so they cannot
silently disappear from a later consolidation.

Of 339 tracked R/Python/shell/SLURM files across the repository, the remaining
135 are 85 files under `tests/` and 50 source-only model/helper libraries under
`code/`. They are excluded as standalone scientific-result producers. They remain
dependencies of callers, not deletion candidates. The input evidence records
sourced modules where detected. Source-only functions that callers use to write
files belong in the dependency design in G.

The folder table records directories that directly contain tracked result or
exhibit artifacts and ignored result artifacts under `results/`. Counts are not
recursive, so nested directory counts can be added without counting a file twice.
It does not list every empty directory or every possible future destination.
The ignored `data/` boundary is included without inspecting its contents.
Ignored artifacts were counted by filename and byte-size metadata only. No RDS
object, patient table or private input was opened.

## The nine result folders

The tracked counts below include source and documentation files, not only outputs.
Sizes are approximate MiB of tracked files in this checkout. Ignored checkpoints
are additional and appear separately in the folder table.

| Folder | Tracked files | MiB | What lives here and its role |
|---|---:|---:|---|
| `results/benchmark_sim/` | 5,112 | 300.7 | Historical benchmarks, supervision simulations, current reported-model evidence, comparator work and mechanism diagnostics. This is the largest mixture of model states. |
| `results/benchmark_snapshot_2026-10-14/` | 3 | < 0.01 | Tracked README, ignore rule and manifest. The associated local fits and other snapshot files are ignored. This is immutable copied evidence, not a live runner destination. |
| `results/cohort_lmm_sim/` | 1 | 0.02 | The synthetic cohort-indicator runner. Its saved RDS path is documented even if the output is absent or ignored. |
| `results/feature_selection/` | 8 | 0.03 | Training-only screening experiment, helper/test files and aggregate outputs. |
| `results/figures/` | 110 | 3.2 | Figures from early benchmark and prior/rank experiments. Several legacy runners still write here. |
| `results/legacy/` | 101 | 17.4 | Archived runners, report sources and saved tables/figures. Some runner literals retain their pre-archive destinations. |
| `results/multi_cohort_sim/` | 30 | 0.49 | Cohort-specific/shared-program simulations, recovery/strength studies, template and scoring helpers, plots. |
| `results/multimodal_real/` | 112 | 1.15 | Matched expression/methylation fits, feature-count and collapse diagnostics, aggregate interpretation and comparison outputs. |
| `results/tables/` | 157 | 0.28 | Early benchmark and prior/rank tables. Their producers may be retired scripts now under `results/legacy/`. |

## Main result families and exhibit consumers

“Current reported evidence” means evidence used by the current reports. It does
not mean that the historical runner can be executed on changed fitter code and
silently replace those results. “Historical” means preserve the study with its
recorded model state. “Retired reference” is used for explicitly archived code.
These are inventory roles, not final retention or deletion decisions. G must
settle which analyses become core rebuild targets and which stay frozen.

| Result family | Producer and output location | Known report/evidence consumer |
|---|---|---|
| Working expression fit, feature counts, screening and null checks | `run_yfb_supervised_real.R` and `summarize_yfb_supervised_real.R` → `outputs/yfb_supervised_real_2026-11-04/`, including ignored fit/risk-score checkpoints in `fits/`; local training cache under `data/` | [Real-data findings](2026-11-04-supervised-real-data-findings.md), November 6 chapter, November writeup and prelim builders. |
| Fixed-model external scoring | `score_yfb_supervised_external.R` → `outputs/yfb_supervised_external_2026-11-05/`; `save_yfb_external_patient_scores.R` → ignored `data/yfb_external_patient_scores_ws1.rds` | [External findings](2026-11-05-supervised-external-scoring-findings.md), donor audit and external-KM builder. |
| Survival-weight simulation, joint coefficient option and near-threshold recovery | `run_yfb_survival_weight_p2000.R`, `run_yfb_supervision_strength.R`, `run_yfb_near_threshold.R` and their summarizers → dated output roots | [Weight findings](2026-11-04-survival-weight-p2000-findings.md), [strength findings](2026-11-05-supervision-strength-findings.md), [near-threshold findings](2026-11-05-near-threshold-findings.md), progress-book and prelim simulation figures. |
| Training-selected real-data survival weight | `run_yfb_supervised_real_refit.R` → `outputs/yfb_supervised_refit_2026-11-06/` | [Refit findings](2026-11-06-supervised-real-refit-findings.md), November writeup/prelim. |
| Fixed comparators, risk groups and program rules | `run_prompt6_*`, `run_prompt7_*` → `outputs/prompt6_comparators/`, `outputs/prompt7_desurv_style/`; local DeSurv/mixed-simulation checkpoints under `data/` | [Comparator findings](2026-11-06-real-data-comparators-findings.md), [program findings](2026-11-06-desurv-style-analyses-findings.md), November writeup and prelim tables. |
| Pruning after the KL fixes | `run_kl_fix_pruning_rerun.R` → default KL-1 root, or the registered KL-2 root supplied as a CLI argument | [KL-1](2026-11-06-kl-fix-pruning-rerun-findings.md), [KL-2](2026-11-06-flat-beta-pruning-rerun-findings.md), November reports. No pruning adoption. |
| Earlier normalized-projection supervision and mechanism screens | `run_yfb_supervision_decision.R`, diagnostic replays and targeted runners → October 30 attempt and targeted roots | [Execution report](../reviews/2026-10-30-codex-execution-report.md), October 30 chapter and historical findings. Retain failed/capped states and original denominators. |
| Matched LB/YFB comparison | `run_lb_yfb_matched_comparison.R` → `outputs/lb_yfb_matched_2026-11-06/` | [Comparison findings](2026-11-06-lb-yfb-matched-comparison-findings.md). Exploratory evidence, not this week's future task E. |
| Earlier preprocessing, rank, pathways and cohort comparisons | May–October runners → undated roots such as `desurv_comparison/`, `k_init_sweep/`, `pathway_enrichment/`, `cohort_beta_comparison/` | Dated findings, decisions and historical meeting chapters. These encode model states that differ from the working fit. |
| Multimodal feature expansion and collapse studies | `results/multimodal_real/` runners → `outputs/feature_sensitivity/`, `outputs/round2/`, `outputs/supervised/` and local ignored failure states | October methods record, attempt log and multimodal findings. The future D fitter must not inherit their adoption status by name alone. |
| Report exhibits | `paper/prelim/make_*`, progress-book and meeting figure builders | `paper/prelim/figures/`, `paper/prelim/tables/`, `docs/progress_book/figures/` and dated presentation assets. Some builders also write aggregate tables back under `results/benchmark_sim/outputs/`. |

In this table, `outputs/` is relative to `results/benchmark_sim/`. Full paths
and source-line evidence are in the script table. The report-source table
retains the actual read/include expressions rather than inferring a complete
pipeline from filename similarity.

## Cleanup and alignment opportunities

1. **Several model states share the same parent folder.** `benchmark_sim/` is a
   location, not a model/version label. A rebuild design needs explicit model
   state, options and input provenance for each result family. It must distinguish
   the current reported evidence from historical benchmarks and failed methods.
2. **Saved paths and current writer defaults can disagree.** The KL-2 pruning
   root is a registered CLI override, and the later single-log/KL rerun uses
   `--out`. Earlier `attempt2`, retry and one-off preflight folders do not become
   reproducible just because a current runner has a similar name. The folder
   table records unresolved writers and historical instances explicitly.
3. **Legacy runners still target pre-archive paths.** `results/legacy/` holds
   reports and outputs, while some scripts still write to `results/tables/`
   and `results/figures/`. G should decide whether each remains an executable
   historical replay or becomes a reference-only source. Moving paths first
   would obscure that choice.
4. **One meeting builder has an overwrite risk.**
   `presentation/walther_lab_meeting_08_27_2026/figs/make_factor_figs.R` hardcodes
   the June meeting's assets directory. Its actual destination is recorded,
   rather than inferred from its own directory. This inventory does not fix it.
5. **Rendering, fitting and deployment have different effects.** Quarto sources
   can read aggregates and create report figures without fitting. Publication
   and prelim sync scripts can write or commit outside their source directory.
   The inventory includes these supporting scripts so a future rebuild does not
   mistake a deployment action for a harmless figure target.
6. **Ignored checkpoints need a retention decision.** The snapshot manifest and
   README are tracked while its fits are not. Other ignored fit directories sit
   beside tracked aggregates. The later design can propose an archive with
   checksums and a restoration test, but no retention rule has been adopted here.
7. **Consolidation should follow dependency checks.** Helpers, SLURM launchers,
   report includes and inherited scripts may be consumers even when no report
   names the producer directly. A duplicate-looking filename is insufficient
   evidence for deletion or merging scripts.

## Limits and handoff to G

The CSVs are a static source-and-filesystem inventory, not an executed dependency
graph or proof that all saved results reproduce. Writer roots are resolved where
source literals, arguments or registrations establish them. Braces denote runtime
templates. Source-line expressions preserve the unresolved detail. Candidate
writers in the folder table can be inherited from a parent output root and are
not file-by-file provenance. Empty or unresolved mappings are retained.

Document consumers were found by literal script-name and output-root references
in tracked source documents. That includes references in historical instructions,
not just active imports. Dynamic paths, sourced helpers and include chains require
further inspection before any relocation. No current/frozen designation should be
inferred solely from a directory date or a script's last modification date.

The nine-folder summary is separate from PROJECT_STATUS.qmd: it describes where
results and their producers live, not the scientific verdicts. The next authorized
step after 0b is a design note for G, with Andrew's go-ahead. That design must settle
core versus historical targets, output provenance/versioning, ignored-checkpoint
retention, Longleaf execution and a proposed layout. Nothing in this inventory
moves files or authorizes a rerun.

## Verification

Coverage was checked against `git ls-files`: 204 script rows plus 135 explicit
source-library/test exclusions account for all 339 executable-source candidates.
All 155 result-folder scripts appear exactly once. All 82 tracked Quarto/R
Markdown sources appear in the report inventory. Local links, source paths,
folder counts and representative path mappings were checked. No model code,
configuration, tests, patient inputs or saved outputs were changed. The full R
suite was not run because this task produced documentation only.

Independent review verified the coverage counts and representative source mappings.
Its duplicate-header finding was corrected and all three CSV headers were checked
for uniqueness. Its request to identify ignored fit/risk-score checkpoints beside
the tracked real-data aggregates was incorporated. All 251 folder rows remain.
No unresolved review defect remains.
