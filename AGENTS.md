# AGENTS.md — multiomicsGEP

Use [PROJECT_STATUS.qmd](PROJECT_STATUS.qmd) for the current research state.
This file holds repository instructions.

## Start here (2026-10-30)

The current plan is [`docs/plans/2026-10-30-realignment-plan.md`](docs/plans/2026-10-30-realignment-plan.md).
Start with its review handoff and the
[`execution report`](docs/reviews/2026-10-30-codex-execution-report.md), then the
October 30 progress-book chapter and issue tracker. The [current-state chapter](docs/progress_book/current-state.qmd) summarizes
priorities and saved examples. Update it alongside the roadmap when the
research state changes. The standalone YFB smoke test still reproduces zero coefficients and activities.
Standalone test scripts now return a failure exit code when assertions fail.
Implementation and registered evaluation are separate from scientific adoption.
The supervision work must receive the independent review before merging.
`rebuild/unified-ssbmf` remains an unmerged reference. Remind Andrew of the
registered recovery/false-active/convergence gate when continuing this work.

The opt-in fitter and numerical kernel are separate from historical defaults.
Use `VBeta`, `VL`, and `VF` as their variance carriers; subtracting large second
moments can erase uncertainty. The seventh simulation registration corrects
shared-loading units; earlier attempts remain numerical/design diagnostics.
The approved S-8 mechanism screen adds opt-in fixed-precision and guarded
updates and a separate exact Gaussian coefficient reference; its registration
is `docs/plans/2026-10-30-yfb-targeted-diagnostics-plan.md`. Stalled guarded
updates are failures, not convergence. Results and limits are in
`docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md`: fixed noise
prevents the tested explosions, but no case passes all screening criteria.
These controls do not adopt a new model.
Reliable survival supervision of the molecular loadings is the top priority.
Start with the S-9 YFB test and its registered recovery, false-active and
convergence criteria. Investigate LB with fold-in prediction if YFB stalls.
Multimodal development is deferred except when a unified fitter run with
M = 1 offers a concrete route to resolving supervision. The single-modality
R package can start once the model meets its acceptance criteria.

## Key Instructions

- **Deliverable and code attribution:** Do not name assistant tools or include tool-attribution lines in presentation deliverables or code. Tool names remain appropriate in explicit review handoffs. Andrew requested this on 2026-10-15; it supersedes earlier template author lines.

- **Canonical CAVI loop:** `code/fit_ssbmf_lb.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **Multimodal YFB proposal:** `derivations/SSBMF-YFB-multimodal-model-derivations/multimodal_YFB_derivation.{qmd,pdf}` contains the concise nine-page linear joint-YFB derivation; `multimodal_YFB_implementation.qmd` in the same directory specifies future R code, and `multimodal_YFB_technical_notes.{qmd,pdf}` preserves the extensive version; `verify_derivation.R` checks the algebra without fitting a model. Read it before multimodal work. An isolated implementation exists on the `codex/multimodal-yfb` branch (`code/*multimodal_yfb*.R`: fitter, updates, prediction, simulation, K-selection runner, with tests); it is work in progress, not the production single-modality model, and its defaults are not approved. Current status and review decisions belong in `PROJECT_STATUS.qmd` and `DECISIONS.md`.
- **No status duplication in instructions:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.qmd` instead.
- **Where documents go:** `docs/` type folders — `plans/`, `findings/`, `reviews/`, `prompts/`, `meetings/`, `notes/` — named `YYYY-MM-DD-topic-type.ext` (all hyphens, date first, the day it was written). No archive folder: mark a replaced document superseded in its folder's README. Add a row (date, link, purpose, status) to that README for every new document. See `docs/README.md`.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. The latest verified count and expected diagnostic warnings are recorded in the execution report; tests must all pass before completion.
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 88/88 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git (`**/PDAC_data/` and the local `/data/` folder of TCGA/ICGC copies are git-ignored; never commit them). Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Model entry points:** LB uses `code/fit_ssbmf_lb.R` (η = Lβ) and YFB uses
  `code/fit_ssbmf_yfb.R` (η = (YF)β). The frozen YFB benchmark uses
  `alpha_F = 0`. See the status report for benchmark estimates and the
  current supervision experiments.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.
- **No mannered prose:** Say what you mean; don't substitute metaphor or flourish for direct statement (e.g. "a dial worth turning" instead of "a parameter worth varying," "this point earns its keep" instead of "this point still matters"). Applies to all writing in this repo — reports, the progress book, plan/decision docs, commit messages, comments. When a literal phrase is available, use it.
- **Punctuation in prose:** In the progress book and reports, avoid semicolons and em or en dashes used as punctuation unless they make the sentence clearer. A period, a comma, or two sentences usually read better. Dashes in numeric ranges (Figures 2–4) and compound names (Kullback–Leibler) are fine. When drafting or revising a progress-book chapter, apply the `human-writer` skill.

- **Protected-pruning evaluation:** opt-in changes and fixed runners are on `m5/k-selection`; adoption was blocked. The controlling fallback decision is recorded at the top of `DECISIONS.md`.

- **Planted comparator evaluation:** runner and partial tables are on `m7/sim-comparators`; see the 2026-10-15 decision for the zero/one-factor EBMF blocker and incomplete high-noise denominator.

- **Multimodal comparator audit:** `results/multimodal_real/run_supervised_multimodal_comparison.R` records actual two-step ranks and the fixed four-arm evaluation on `m7/multimodal-comparators`; see the dated progress-book and decision entries for results.

## Stratified YFB scoring

Use fold-specific training centers for centered YFB CV. Keep concordance
pairs within Cox strata and distinguish `strata_id` from cohort-specific
coefficient labels. See the
[stratified CV review handoff](docs/reviews/2026-11-06-stratified-yfb-cv-handoff.md)
for the regression checks and the separate canonical scaling diagnostic.

## Robustness and comparator follow-up

The writeup's real-data workflow figures are generated by
`paper/writeup/2026-10-16-writeup/figures/make_real_data_workflows.R` from audited aggregate
counts and the historical pipeline description. They distinguish primary
ICGC validation from its nested broader-histology sensitivity set and flag
repeated single-modality PACA-AU donors. Rebuild instructions are in the
writeup README; stage shapes distinguish data, processing, models and
evaluation, and dashed notes explain rank selection and patient matching.
Regenerating these figures requires no patient data or fit.

For the writeup accuracy revision, consult
`paper/writeup/2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md` rather than
treating the original review as a live checklist. Appendix C holds retained
diagnostics/provenance; reporting fixes do not authorize donor-aware reanalysis,
new preprocessing fits or external model selection.

The five-prompt branch/results index is
`docs/findings/2026-10-15-m6-m7-summary.md`; read its limitations before using
the tables. E did not adopt protection or a new K. `m6/robustness` contains the
training-only robustness runner; current findings belong in `PROJECT_STATUS.qmd`
and `DECISIONS.md`, not this instruction file.

The merged Round2 diagnostics are historical evidence. Their specification is
`docs/plans/2026-10-16-round2-handoff-plan.md`; status and attempts are in
`docs/plans/2026-10-16-round2-return-note.md` and `docs/plans/2026-10-16-round2-attempts-log.md`.
The benchmark tag and `results/benchmark_snapshot_2026-10-14/` are immutable.

## Quick Reference

Active pointers only — completed one-off analyses and dated reports are indexed in `PROJECT_STATUS.qmd`,
`ROADMAP.md`, and `DECISIONS.md` instead of here.

| What | Where |
|------|-------|
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_ssbmf_lb.R` (factor-wise, canonical) |
| Hold-out prediction | `code/predict_ssbmf_lb.R` — `predict_ssbmf_lb()` (projection), `predict_ssbmf_lb_foldin()` (LB fold-in under the stored priors) |
| Train/test splitting | `code/train_test_split.R` — `stratified_split()` |
| Feature selection | `code/feature_selection.R` — `cox_feature_selection()` |
| K selection | `code/select_K.R` — `auto_prune_K()`, `select_K_cv()` |
| Full ELBO computation | `code/compute_elbo.R` — `compute_ebnm_kl()`, `compute_survival_elbo()`, `compute_normal_kl()` |
| Joint log-likelihood / BIC for K_init sweeps | `code/compute_bic.R` — `compute_joint_ll_bic()`; source after `fit_ssbmf_yfb.R` |
| Cross-validated held-out log-likelihood | `code/compute_cv_loglik.R` — `cv_survival_loglik()`, `bicv_genomics_loglik()`; companion to `compute_bic.R`, reported separately |
| Cohort-specific survival coefficients | `code/update_beta_cohort.R` — `beta_cohort_id` argument on `fit_ssbmf_yfb()`. Distinct from `cohort_id` (genomics offset) and `strata_id` (baseline hazard) |
| Cohort F update | `code/update_F_cohort.R` — `update_F_cohort_all()` (Normal conjugate) |
| Companion doc for fit_ssbmf_lb.R | `docs/lb_walkthroughs/fit_modular.qmd` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **K_init consensus sweep** | `results/benchmark_sim/run_k_init_sweep.R` — ELBO/BIC/log-likelihood/held-out-C consensus; `--reuse-cache`/`--reuse-cv-cache` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Pathway enrichment functions | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **L-update debugging guide** | `docs/lb_walkthroughs/update_L_fix.md` — read before any L-update work |
| Corrected derivations | `derivations/SSBMF-LB-model-derivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| **YFB full derivation + algorithm summary** | `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.{tex,pdf}`; run `verify_derivation.R` in that folder after changing an update. `components/` holds the superseded April pieces |
| LB vs YFB comparison | `derivations/SSBMF-LB-vs-YFB-comparison/SSBMF-LB-vs-YFB-comparison.{tex,pdf}` — precision-share scaling, fold-in prediction, LB code review |
| **Current working plan** | `docs/plans/2026-10-30-realignment-plan.md` (execution handoff and independent review pending); the October 26 parallel-work and October 22 rebuild plans are historical context. The October 15–16 diagnostic work (merged 2026-10-26, tag `archive/round2-supervised-2026-10`) is summarized in `docs/plans/2026-10-16-round2-return-note.md`, with every registered study and outcome in `docs/plans/2026-10-16-round2-attempts-log.md`; `docs/meetings/2026-10-02-meeting-notes-followup.md` preserves the October 2 advisor requirements. |
| Matched TCGA/ICGC expression + methylation loader | `code/load_multiomics_data.R` — `build_multiomics_cohorts()`; real-data checks `tests/test_real_multiomics_loading.R` (local) |
| Multimodal parsimony framework (intercept, ebnm priors, ELBO pruning) and compiled sweep | `code/fit_multimodal_yfb.R` controls `intercept`/`prior_update`/`prune`/`tau_model`; `code/multimodal_yfb_sweep.cpp` |
| Multimodal simulation and real-data runners | `results/benchmark_sim/run_multimodal_pruning_comparison.R`; `results/multimodal_real/run_multimodal_real_fit.R` + `summarize_multimodal_real.R` |
| Prelim proposal chapter source | `paper/prelim/project3-ssbmf.qmd` → synced to `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` by the installed post-commit hook (or `paper/prelim/tools/sync_to_prelim.sh` by hand; never pushes) |
| Prelim proposal chapter plan | `docs/plans/2026-10-01-prelim-proposal-plan.md` — outline, sync to `bios-dissertation`, α_F framing constraint |
| Advisor-facing writeup (interim version approved 2026-10-14; revise after the method work via Prompt 7) | `paper/writeup/2026-10-16-writeup/writeup.qmd` → `writeup.html`/`writeup.pdf`, appendices in `paper/writeup/2026-10-16-writeup/appendices/`; planning page, source inventory, implementation plan and prompts alongside. Each version has its own dated folder (`paper/writeup/YYYY-MM-DD-writeup/`, indexed in `paper/writeup/README.md`); do not edit a circulated version, start the next one. Read `PROJECT_STATUS.qmd` for the review state. |

**Opt-in diagnostic controls from the October 15–16 work** (merged 2026-10-26; all defaults
unchanged): protected and posterior-preserving pruning (`prune_preserve_posterior`), profiled Cox
intercept (`profile_cox_intercept`, `code/profile_cox_intercept.R`), log-space Cox risk-set sums
(`control$cox_arithmetic = "logspace"`), fixed-posterior loading prior (`loading_prior_update`,
`code/loading_prior_fixed_q.R`), feature-sensitivity helpers (`code/feature_sensitivity.R`,
`code/cpg_feature_sensitivity.R`) and the multimodal program-agreement and loading-balance helpers.
None is adopted. Their runners live in `results/benchmark_sim/` and `results/multimodal_real/`; read
the study's entry in `docs/plans/2026-10-16-round2-attempts-log.md` before running one, and do not
repeat a completed study without a new registered entry.

The tracker includes standalone architecture, update-cycle and prediction
SVGs under `docs/progress_book/figures/`, including the matrix comparison of
two-step and joint fitting. Rebuild them with
`python3 docs/progress_book/figures/build_diagrams.py`. The [code map](code/README.md)
separates LB, historical YFB, experimental YFB, multimodal and shared scripts.

Maintain the tracker’s Completed Results checklist alongside saved evidence.
Count a checked diagnostic as a completed result while keeping its acceptance
verdict explicit. Update the displayed checklist count when adding or
changing items.

The October methods development record is included in the October 16 chapter
from `docs/progress_book/chapters/_2026-10-15-methods-record.qmd`. Keep the
October 15 public URL as an alias when regenerating the book.

Record tracker decisions from the project chat in the chapter source.

The November 6 notebook includes an aggregate comparison of saved canonical
YFB and rank-7 EBMF fits. Rebuild its figures with
`Rscript docs/progress_book/figures/verify_october_yfb.R` from the repository root.
The script uses local PDAC training data and historical preprocessing for
alignment with the saved fits. Keep patient-level inputs and outputs out of Git.

Regenerate the November 6 planted-factor simulation figure with
`Rscript docs/progress_book/figures/plot_planted_yfb.R`. The plot reports
all replicates for recovery and only defined values for concordance.
Use `Rscript docs/progress_book/figures/plot_planted_yfb_comparators.R`
to rebuild the loading-recovery and matched-comparator figures from the
saved synthetic replicates.

Evaluate held-out Cox scores with `frozen_reverse_cindex()`. Historical caches
can store good-prognosis scores. Use the documented cache direction or recover
it on training data before evaluating held-out outcomes. See the
[orientation audit](docs/findings/2026-11-06-orientation-audit.md).

For the S-12 pooled external interval, keep all records from each PACA-AU
donor together. Use ICGC donor metadata to link assays, and write only
aggregate output. The external-cohort audit records the eligible scope and
historical uses.

For the flat β likelihood regression, run
`Rscript tests/test_empty_program_kl.R` (also in the full suite). The shared
prior update treats exactly flat coefficient likelihoods as uninformative.
Use the [KL-2 registration](docs/plans/2026-11-06-flat-beta-pruning-rerun-plan.md)
to locate the repeat of the pruning comparison.


For the supervised Chapter 4 figures, use the three rebuild commands in
`README.md` and trace values through
`paper/prelim/notes/chapter4-number-sources-2026-11-06.md`. The external KM
builder consumes fixed scores under the ignored `data/` directory. Keep
patient-level inputs and outputs there, and commit only aggregate exhibits.

For the fixed S-11 real-data comparators, use the two commands in README.md
and the [aggregate source table](docs/findings/2026-11-06-real-data-comparators-findings.md).
Keep the DeSurv fit checkpoints under ignored `data/`.

For the S-11 risk groups, EBMF rank comparison and program-rule audit, use
the four commands in README.md and the
[review handoff](docs/reviews/2026-11-06-comparator-and-program-review-handoff.md).
The S-11 precision audit uses a final-state bound because the saved fit
does not contain loading second moments. Keep the simulation working
checkpoints under ignored `data/`.
