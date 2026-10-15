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

The five-prompt branch/results index is
`docs/reports/Codex_M6_M7_Summary_10_15_26.md`; read its limitations before using
the tables. E did not adopt protection or a new K. `m6/robustness` contains the
training-only robustness runner; current findings belong in `PROJECT_STATUS.qmd`
and `DECISIONS.md`, not this instruction file.

Round 2 work is cumulative on `round2/supervised`. Its controlling specification is
`docs/plans/Codex_Handoff_Round2_10_16_26.md`; status and attempts are in
`docs/plans/Codex_to_Claude_Round2.md` and `docs/plans/Round2_attempts_log.md`.
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
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **K_init consensus sweep** | `results/benchmark_sim/run_k_init_sweep.R` — ELBO/BIC/log-likelihood/held-out-C consensus; `--reuse-cache`/`--reuse-cv-cache` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Pathway enrichment functions | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **L-update debugging guide** | `docs/update_L_fix.md` — read before any L-update work |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| **Current working plan** | `docs/plans/Meeting_Notes_Followup_10_2_26.md` — start open work here (write-up, feature-count sensitivity, survival weighting); `Working_Plan_10_1_26.md` records the 10/1–10/2 work |
| Matched TCGA/ICGC expression + methylation loader | `code/load_multiomics_data.R` — `build_multiomics_cohorts()`; real-data checks `tests/test_real_multiomics_loading.R` (local) |
| Multimodal parsimony framework (intercept, ebnm priors, ELBO pruning) and compiled sweep | `code/fit_multimodal_yfb.R` controls `intercept`/`prior_update`/`prune`/`tau_model`; `code/multimodal_yfb_sweep.cpp` |
| Multimodal simulation and real-data runners | `results/benchmark_sim/run_multimodal_pruning_comparison.R`; `results/multimodal_real/run_multimodal_real_fit.R` + `summarize_multimodal_real.R` |
| Prelim proposal chapter source | `paper/prelim/project3-ssbmf.qmd` → synced to `bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` by the installed post-commit hook (or `paper/prelim/tools/sync_to_prelim.sh` by hand; never pushes) |
| Prelim proposal chapter plan | `docs/plans/Prelim_Proposal_Plan_10_1_26.md` — outline, sync to `bios-dissertation`, α_F framing constraint |
| Advisor-facing writeup (interim version approved 2026-10-14; revise after the method work via Prompt 7) | `paper/writeup/writeup.qmd` → `writeup.html`/`writeup.pdf`, appendices in `paper/writeup/appendices/`; planning page, source inventory, implementation plan and prompts in `paper/writeup/_process/`. Read `PROJECT_STATUS.qmd` for the review state. |

The opt-in single pruning observer is
`results/benchmark_sim/run_single_pruning_trace.R`; it records process-local
nullcheck states without editing the fitter. Its fixed diagnosis and the rank
panel's Longleaf wrapper are documented in the Round 2 attempt history.
Andrew approved an additional opt-in Task1 remedy after the stage trace:
`prune_preserve_posterior=TRUE` carries the retained full single-modality state
into pruning refits. `run_posterior_k_panel.R` evaluates the unchanged 160-cell
protected single panel against the saved original reference. The completed
experiment improves recovery but fails the unchanged adoption rule (recovery, convergence and rank); see the attempt log. Defaults and the
benchmark remain unchanged.

Task2's fixed TCGA-only diagnostic runner is
`results/multimodal_real/run_round2_collapse_diagnosis.R`. Its opt-in trace and
diagnostic prior controls preserve the default fit. Read the attempt history
before execution or interpretation; training contrasts do not adopt a remedy,
and external cohorts cannot select one.

Task3's pinned DeSurv package is installed in ignored `data/round2_R_library/`;
its preprocessing/projection adapter, driver and Longleaf wrappers are prepared
and tested. Andrew deferred the head-to-head for review. Read the Round2 handoff for author SHAs, installation details and the
raw-loader log-transform concern before changing any preprocessing.

The ten-arm TCGA diagnosis has completed; consult the Round2 handoff for
benchmark-linked results and the open remedy ladder. Andrew deferred Task3's
real-data DeSurv head-to-head for review. Prepared DeSurv commands above must
not be executed without renewed authorization; no real-data DeSurv fit or
external scoring has run. Task4 omits DeSurv as explicitly approved.

Task4's comparator runner now preserves all method rows and full denominators
for EBMF rank0/rank1 and isolated failures. The Round2 attempt log controls
its five-arm rerun; DeSurv is omitted by Andrew's explicit decision.

Task2 R1's fixed-prior evaluation is `run_round2_fixed_prior.R` in the benchmark
folder. Its committed attempt history governs80 simulations and20 TCGA CV fits;
all criteria and defaults remain unchanged, with no external data for selection.
