# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

This is the single source of truth for agent instructions in this repo. `AGENTS.md` is a symlink to
this file (so Codex, Antigravity, and other tools that look for `AGENTS.md` read the same content) —
edit `CLAUDE.md`, never `AGENTS.md` directly. If a `GEMINI.md` or similar is ever added, symlink it
here too rather than maintaining a separate copy.

---

## Start here (2026-10-30)

The current plan is [`docs/plans/2026-10-30-realignment-plan.md`](docs/plans/2026-10-30-realignment-plan.md):
first clean up the open review items (Phase C: KL normalizer, single log, orientation in the comparator scripts,
sign and β-variance checks), then make survival supervise F under η = (YF)β, following
[`docs/plans/2026-10-30-Unified-Architecture-Alignment.md`](docs/plans/2026-10-30-Unified-Architecture-Alignment.md)
(Phase S, with a pre-registered Gate S). The execution prompt is
[`docs/prompts/2026-10-30-codex-execution-prompt.md`](docs/prompts/2026-10-30-codex-execution-prompt.md). The 10/30
progress-book chapter and its issue tracker record the progress. `rebuild/unified-ssbmf` is kept unmerged as a reference.
Remind Andrew of the plan and its gate at the start of a session.

## Key Instructions

- **Deliverable and code attribution:** Do not name assistant tools or include tool-attribution lines in presentation deliverables or code. Tool names remain appropriate in explicit review handoffs. Andrew requested this on 2026-10-15; it supersedes earlier template author lines.

- **Canonical CAVI loop:** `code/fit_modular.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **Multimodal YFB proposal:** `derivations/multimodal_YFB/multimodal_YFB_derivation.{qmd,pdf}` contains the concise nine-page linear joint-YFB derivation; `multimodal_YFB_implementation.qmd` in the same directory specifies future R code, and `multimodal_YFB_technical_notes.{qmd,pdf}` preserves the extensive version; `verify_derivation.R` checks the algebra without fitting a model. Read it before multimodal work. An isolated implementation exists on the `codex/multimodal-yfb` branch (`code/*multimodal_yfb*.R`: fitter, updates, prediction, simulation, K-selection runner, with tests); it is work in progress, not the production single-modality model, and its defaults are not approved. Current status and review decisions belong in `PROJECT_STATUS.qmd` and `DECISIONS.md`.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.qmd` instead.
- **Where documents go:** `docs/` type folders — `plans/`, `findings/`, `reviews/`, `prompts/`, `meetings/`, `notes/` — named `YYYY-MM-DD-topic-type.ext` (all hyphens, date first, the day it was written). No archive folder: mark a replaced document superseded in its folder's README. Add a row (date, link, purpose, status) to that README for every new document. See `docs/README.md`.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: all passing (532 as of 2026-10-14 on `fix/multimodal-monitor-objective`).
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 88/88 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git (`**/PDAC_data/` and the local `/data/` folder of TCGA/ICGC copies are git-ignored; never commit them). Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Current model status:** Both LB (`code/fit_modular.R`, η = Lβ) and YFB (`code/fit_cox_on_yf.R`, η = (YF)β) are fully implemented. **Recommended configuration:** YFB × per-platform z-std × DeSurv gene selection × no cohort indicator, K=7 default (K selection is a two-stage ARD framework — K_init from an ELBO/BIC/log-likelihood/CV-C consensus, then `classify_factors()` shrinkage gives K_eff from that one over-specified fit), K_eff=2, mean external C≈0.627 across 5 held-out PDAC cohorts. Optional `cohort_id` (genomics offset), `strata_id` (stratified baseline hazard, off by default — performance-neutral), and `beta_cohort_id` (cohort-specific survival coefficients) parameters exist; see the Quick Reference table below. Full history, benchmark numbers, and the reasoning behind each choice live in `PROJECT_STATUS.qmd`; the dated rationale entries are in `DECISIONS.md` — read those before revisiting any of this, don't rely on this summary alone.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.
- **No mannered prose:** Say what you mean; don't substitute metaphor or flourish for direct statement (e.g. "a dial worth turning" instead of "a parameter worth varying," "this point earns its keep" instead of "this point still matters"). Applies to all writing in this repo — reports, the progress book, plan/decision docs, commit messages, comments. When a literal phrase is available, use it.

- **Protected-pruning evaluation:** opt-in changes and fixed runners are on `m5/k-selection`; adoption was blocked. The controlling fallback decision is recorded at the top of `DECISIONS.md`.

- **Planted comparator evaluation:** runner and partial tables are on `m7/sim-comparators`; see the 2026-10-15 decision for the zero/one-factor EBMF blocker and incomplete high-noise denominator.

- **Multimodal comparator audit:** `results/multimodal_real/run_supervised_multimodal_comparison.R` records actual two-step ranks and the fixed four-arm evaluation on `m7/multimodal-comparators`; see the dated progress-book and decision entries for results.

## Robustness and comparator follow-up

The writeup's real-data workflow figures are generated by
`paper/writeup/figures/make_real_data_workflows.R` from audited aggregate
counts and the historical pipeline description. They distinguish primary
ICGC validation from its nested broader-histology sensitivity set and flag
repeated single-modality PACA-AU donors. Rebuild instructions are in the
writeup README; stage shapes distinguish data, processing, models and
evaluation, and dashed notes explain rank selection and patient matching.
Regenerating these figures requires no patient data or fit.

For the writeup accuracy revision, consult
`paper/writeup/_process/Accuracy_Review_Dispositions_10_16_26.md` rather than
treating the original review as a live checklist. Appendix C holds retained
diagnostics/provenance; reporting fixes do not authorize donor-aware reanalysis,
new preprocessing fits or external model selection.

The five-prompt branch/results index is
`docs/findings/2026-10-15-m6-m7-summary.md`; read its limitations before using
the tables. E did not adopt protection or a new K. `m6/robustness` contains the
training-only robustness runner; current findings belong in `PROJECT_STATUS.qmd`
and `DECISIONS.md`, not this instruction file.

Round 2 work is cumulative on `round2/supervised`. Its controlling specification is
`docs/plans/2026-10-16-round2-handoff-plan.md`; status and attempts are in
`docs/plans/2026-10-16-round2-return-note.md` and `docs/plans/2026-10-16-round2-attempts-log.md`.
The benchmark tag and `results/benchmark_snapshot_2026-10-14/` are immutable.

## Quick Reference

Active pointers only — completed one-off analyses and dated reports are indexed in `PROJECT_STATUS.qmd`,
`ROADMAP.md`, and `DECISIONS.md` instead of here.

| What | Where |
|------|-------|
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_modular.R` (factor-wise, canonical) |
| Hold-out prediction | `code/predict.R` — `predict_supervised_mf()` |
| Train/test splitting | `code/train_test_split.R` — `stratified_split()` |
| Feature selection | `code/feature_selection.R` — `cox_feature_selection()` |
| K selection | `code/select_K.R` — `auto_prune_K()`, `select_K_cv()` |
| Full ELBO computation | `code/compute_elbo.R` — `compute_ebnm_kl()`, `compute_survival_elbo()`, `compute_normal_kl()` |
| Joint log-likelihood / BIC for K_init sweeps | `code/compute_bic.R` — `compute_joint_ll_bic()`; source after `fit_cox_on_yf.R` |
| Cross-validated held-out log-likelihood | `code/compute_cv_loglik.R` — `cv_survival_loglik()`, `bicv_genomics_loglik()`; companion to `compute_bic.R`, reported separately |
| Cohort-specific survival coefficients | `code/update_beta_cohort.R` — `beta_cohort_id` argument on `fit_cox_on_yf()`. Distinct from `cohort_id` (genomics offset) and `strata_id` (baseline hazard) |
| Cohort F update | `code/update_F_cohort.R` — `update_F_cohort_all()` (Normal conjugate) |
| Companion doc for fit_modular.R | `docs/lb_walkthroughs/fit_modular.qmd` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **K_init consensus sweep** | `results/benchmark_sim/run_k_init_sweep.R` — ELBO/BIC/log-likelihood/held-out-C consensus; `--reuse-cache`/`--reuse-cv-cache` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Pathway enrichment functions | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **L-update debugging guide** | `docs/lb_walkthroughs/update_L_fix.md` — read before any L-update work |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| **Current working plan** | `docs/plans/2026-10-26-parallel-work-sessions-plan.md` and `docs/plans/2026-10-22-unified-method-rebuild-plan.md`. The October 15–16 diagnostic work (merged 2026-10-26, tag `archive/round2-supervised-2026-10`) is summarized in `docs/plans/2026-10-16-round2-return-note.md`, with every registered study and outcome in `docs/plans/2026-10-16-round2-attempts-log.md`; `docs/meetings/2026-10-02-meeting-notes-followup.md` preserves the October 2 advisor requirements. |
| Matched TCGA/ICGC expression + methylation loader | `code/load_multiomics_data.R` — `build_multiomics_cohorts()`; real-data checks `tests/test_real_multiomics_loading.R` (local) |
| Multimodal parsimony framework (intercept, ebnm priors, ELBO pruning) and compiled sweep | `code/fit_multimodal_yfb.R` controls `intercept`/`prior_update`/`prune`/`tau_model`; `code/multimodal_yfb_sweep.cpp` |
| Multimodal simulation and real-data runners | `results/benchmark_sim/run_multimodal_pruning_comparison.R`; `results/multimodal_real/run_multimodal_real_fit.R` + `summarize_multimodal_real.R` |
| Prelim proposal chapter source | `paper/prelim/project3-ssbmf.qmd` → synced to `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` by the installed post-commit hook (or `paper/prelim/tools/sync_to_prelim.sh` by hand; never pushes) |
| Prelim proposal chapter plan | `docs/plans/2026-10-01-prelim-proposal-plan.md` — outline, sync to `bios-dissertation`, α_F framing constraint |
| Advisor-facing writeup (interim version approved 2026-10-14; revise after the method work via Prompt 7) | `paper/writeup/writeup.qmd` → `writeup.html`/`writeup.pdf`, appendices in `paper/writeup/appendices/`; planning page, source inventory, implementation plan and prompts alongside. Read `PROJECT_STATUS.qmd` for the review state. |

**Opt-in diagnostic controls from the October 15–16 work** (merged 2026-10-26; all defaults
unchanged): protected and posterior-preserving pruning (`prune_preserve_posterior`), profiled Cox
intercept (`profile_cox_intercept`, `code/profile_cox_intercept.R`), log-space Cox risk-set sums
(`control$cox_arithmetic = "logspace"`), fixed-posterior loading prior (`loading_prior_update`,
`code/loading_prior_fixed_q.R`), feature-sensitivity helpers (`code/feature_sensitivity.R`,
`code/cpg_feature_sensitivity.R`) and the multimodal program-agreement and loading-balance helpers.
None is adopted. Their runners live in `results/benchmark_sim/` and `results/multimodal_real/`; read
the study's entry in `docs/plans/2026-10-16-round2-attempts-log.md` before running one, and do not
repeat a completed study without a new registered entry.
