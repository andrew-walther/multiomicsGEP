# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

This is the single source of truth for agent instructions in this repo. `AGENTS.md` is a symlink to
this file (so Codex, Antigravity, and other tools that look for `AGENTS.md` read the same content) —
edit `CLAUDE.md`, never `AGENTS.md` directly. If a `GEMINI.md` or similar is ever added, symlink it
here too rather than maintaining a separate copy.

---

## Key Instructions

- **Canonical CAVI loop:** `code/fit_modular.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **Multimodal YFB proposal:** `derivations/multimodal_YFB/multimodal_YFB_derivation.{qmd,pdf}` contains the concise nine-page linear joint-YFB derivation; `multimodal_YFB_implementation.qmd` in the same directory specifies future R code, and `multimodal_YFB_technical_notes.{qmd,pdf}` preserves the extensive version; `verify_derivation.R` checks the algebra without fitting a model. Read it before multimodal work. An isolated implementation exists on the `codex/multimodal-yfb` branch (`code/*multimodal_yfb*.R`: fitter, updates, prediction, simulation, K-selection runner, with tests); it is work in progress, not the production single-modality model, and its defaults are not approved. Current status and review decisions belong in `PROJECT_STATUS.qmd` and `DECISIONS.md`.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.qmd` instead.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **License:** MIT (`LICENSE`).
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: all passing (530 as of 2026-10-02 on `codex/multimodal-yfb`).
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 88/88 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git (`**/PDAC_data/` and the local `/data/` folder of TCGA/ICGC copies are git-ignored; never commit them). Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Current model status:** Both LB (`code/fit_modular.R`, η = Lβ) and YFB (`code/fit_cox_on_yf.R`, η = (YF)β) are fully implemented. **Recommended configuration:** YFB × per-platform z-std × DeSurv gene selection × no cohort indicator, K=7 default (K selection is a two-stage ARD framework — K_init from an ELBO/BIC/log-likelihood/CV-C consensus, then `classify_factors()` shrinkage gives K_eff from that one over-specified fit), K_eff=2, mean external C≈0.627 across 5 held-out PDAC cohorts. Optional `cohort_id` (genomics offset), `strata_id` (stratified baseline hazard, off by default — performance-neutral), and `beta_cohort_id` (cohort-specific survival coefficients) parameters exist; see the Quick Reference table below. Full history, benchmark numbers, and the reasoning behind each choice live in `PROJECT_STATUS.qmd`; the dated rationale entries are in `DECISIONS.md` — read those before revisiting any of this, don't rely on this summary alone.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.
- **No mannered prose:** Say what you mean; don't substitute metaphor or flourish for direct statement (e.g. "a dial worth turning" instead of "a parameter worth varying," "this point earns its keep" instead of "this point still matters"). Applies to all writing in this repo — reports, the progress book, plan/decision docs, commit messages, comments. When a literal phrase is available, use it.

## Quick Reference

Active pointers only — completed one-off analyses and dated reports are indexed in `PROJECT_STATUS.qmd`,
`ROADMAP.md`, and `DECISIONS.md` instead of here.

| What | Where |
|------|-------|
| Full project docs & session log | `PROJECT_STATUS.qmd` |
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
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| Global hyperparameter registry | `config/globals.yml` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **K_init consensus sweep** | `results/benchmark_sim/run_k_init_sweep.R` — ELBO/BIC/log-likelihood/held-out-C consensus; `--reuse-cache`/`--reuse-cv-cache` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Pathway enrichment functions | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **L-update debugging guide** | `docs/update_L_fix.md` — read before any L-update work |
| Test suite | `tests/run_tests.R` (530 passing, 2026-10-02) |
| Real-data test suite | `tests/test_real_data_loading.R` (88/88, local-only) |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| Architectural decisions log | `DECISIONS.md` |
| Prioritized next steps | `ROADMAP.md` |
| **Current working plan** | `docs/plans/Working_Plan_10_1_26.md` — start open work here |
| Matched TCGA/ICGC expression + methylation loader | `code/load_multiomics_data.R` — `build_multiomics_cohorts()`; real-data checks `tests/test_real_multiomics_loading.R` (local) |
| Multimodal parsimony framework (intercept, ebnm priors, ELBO pruning) and compiled sweep | `code/fit_multimodal_yfb.R` controls `intercept`/`prior_update`/`prune`/`tau_model`; `code/multimodal_yfb_sweep.cpp` |
| Multimodal simulation and real-data runners | `results/benchmark_sim/run_multimodal_pruning_comparison.R`; `results/multimodal_real/run_multimodal_real_fit.R` + `summarize_multimodal_real.R` |
| Prelim proposal chapter source | `paper/prelim/project3-ssbmf.qmd` → synced to `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` by running `paper/prelim/tools/sync_to_prelim.sh` by hand after committing (no hook; never pushes) |
| Prelim proposal chapter plan | `docs/plans/Prelim_Proposal_Plan_10_1_26.md` — outline, sync to `bios-dissertation`, α_F framing constraint |
