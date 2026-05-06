# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

---

## Key Instructions

- **Canonical CAVI loop:** `code/fit_modular.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.md` instead.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: 171/171 passing.
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 77/77 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git. Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Current model status:** Both LB (`code/fit_modular.R`, η = Lβ) and YFB (`code/fit_cox_on_yf.R`, η = (YF)β) are fully implemented with training concordance sign correction. LB external C-index: 0.51–0.67 across all three training modes (tcga_only, cptac_only, merged). YFB external C-index: 0.55–0.63 on single-cohort modes (normal prior); merged training collapses β→0 for YFB (platform mixing, RNA-seq + proteomics). See `ROADMAP.md` for open items.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.

## Quick Reference

| What | Where |
|------|-------|
| Full project docs & session log | `PROJECT_STATUS.md` |
| Code quick-reference (math ↔ R) | `code/SupervisedMF_Context.md` |
| **Reusable CAVI fitting function** | `code/fit_modular.R` (factor-wise, canonical) |
| Hold-out prediction | `code/predict.R` — `predict_supervised_mf()` |
| Train/test splitting | `code/train_test_split.R` — `stratified_split()` |
| Feature selection | `code/feature_selection.R` — `cox_feature_selection()` |
| K selection | `code/select_K.R` — `auto_prune_K()`, `select_K_cv()` stub |
| Full ELBO computation | `code/compute_elbo.R` — `compute_ebnm_kl()`, `compute_survival_elbo()` |
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| Global hyperparameter registry | `config/globals.yml` |
| **Cluster A benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **Cluster B benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **Benchmark reports (dated)** | `docs/reports/ssbmf_summary_report_MM_DD_YY.{qmd,pdf,html}` — DeSurv record: `_04_29_26` |
| Phase 1 loading heatmaps | `results/benchmark_sim/run_phase1_diagnostics.R` |
| Archived benchmark runners | `results/benchmark_sim/archive/` — 9 retired scripts (see archive/README.md) |
| **β=0 design doc (Cluster A/B)** | `docs/beta_zero_fix_design.md` — five-phase plan |
| **L-update debugging guide** | `docs/update_L_fix.md` — read before any L-update work |
| Exploratory simulation runner (legacy) | `results/legacy/modular_sim_factor/run_factor_modular_simulation.R` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Test suite (core + predict) | `tests/run_tests.R` (171/171) |
| Real-data test suite | `tests/test_real_data_loading.R` (77/77, local-only) |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| Architectural decisions log | `DECISIONS.md` |
| Prioritized next steps | `ROADMAP.md` |
