# CLAUDE.md — multiomicsGEP

For full project context, see **[`PROJECT_STATUS.qmd`](PROJECT_STATUS.qmd)** (renders to `PROJECT_STATUS.pdf`).

---

## Key Instructions

- **Canonical CAVI loop:** `code/fit_modular.R`. `Supervised_Bayesian_MF_V2.R` is reference-only — do not extend. V1 (`code/legacy/`) — do not modify. (`full_sim/` and `modular_sim_block/` are legacy/deprecated.)
- **Formal benchmark pipeline:** `results/benchmark_sim/` — alpha CV, external validation, DeSurv comparison. Exploratory/development fits lived in `results/modular_sim_factor/` (now archived to `results/legacy/modular_sim_factor/`).
- **Modular updates:** `code/update_beta.R`, `code/update_L.R`, `code/update_F.R`, `code/update_tau.R`.
- **Global constants:** `config/globals.yml` — all hyperparameters (lambda, alpha grid, K thresholds, DGP params). Never hardcode values defined here.
- **No `CLAUDE.md` duplication:** Do not maintain a second copy of project status here — update `PROJECT_STATUS.qmd` instead.
- **Living documents:** Update `DECISIONS.md` when making any architectural choice (algorithm variant, hyperparameter decision, design tradeoff). Update `ROADMAP.md` when completing a milestone or identifying a new priority.
- **Commit style:** Detailed messages explaining what changed and why; no "Co-Authored-By" lines; no "Session N:" prefixes.
- **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: 412/412 passing.
- **Real-data tests:** `Rscript tests/test_real_data_loading.R` — 88/88 passing (auto-skips if `PDAC_DATA_ROOT` not set).
- **Real data:** Not in git. Stored locally at `~/Library/CloudStorage/OneDrive-.../UNC Dissertation (Liu)/PDAC_data`. For Longleaf: `export PDAC_DATA_ROOT=/proj/rashidlab/data/PDAC`.
- **Current model status:** Both LB (`code/fit_modular.R`, η = Lβ) and YFB (`code/fit_cox_on_yf.R`, η = (YF)β) are fully implemented with training concordance sign correction. Both support a `cohort_id` parameter (corner-point encoding) to absorb platform offsets, and an optional `strata_id` parameter (stratified Cox partial likelihood — study-specific baseline hazard, Breslow risk sets within study, no parametric baseline); `strata_id` is performance-neutral on D4 (mean external C 0.6267→0.6263) so it is off by default (`DECISIONS.md` 2026-07-15). **Recommended configuration (DeSurv-aligned, D4):** YFB × per-platform z-std × DeSurv gene selection (combined_rank, top-3000 per cohort before normalization) × no cohort indicator, K_init=7 (re-confirmed 2026-07-12 under corrected code — a genuine, non-noise result, not a stale artifact; see `DECISIONS.md`). **As of 2026-08-27, K selection is a two-stage ARD framework, not cross-validation**: K_init is chosen from a consensus of ELBO, BIC, an ELBO-style joint log-likelihood, and cross-validated external C-index (these need not agree — see `DECISIONS.md` 2026-08-27 for the full K_init=2..20 sweep and where they disagree), and ARD shrinkage (`classify_factors()`) then determines K_eff from that single over-specified fit — no separate per-K validation loop for K_eff itself. **mean external C=0.627 across 5 held-out PDAC cohorts, K_eff=2** (post-Phase-1: `DECISIONS.md` 2026-07-12). DeSurv gene selection yields 2064 genes vs 2000. Per-platform z-standardization remains essential: 10 of 12 non-per-platform configurations collapse to β=0 in the 18-config extended benchmark. Any preprocessing that does not normalize per-platform before merging (joint z-std, log-only, or joint QN without rank transform) collapses β→0 for both LB and YFB. **Phase 1 honest conclusion (superseding an earlier, incorrect "C=0.642/K_eff=4" figure from the same session — see `DECISIONS.md` 2026-07-12):** Phase 1a's objective normalization provides **no performance benefit for YFB** — YFB's L is pure-genomics and β is pure-survival, so they never share a coordinate that needs rebalancing; an earlier implementation detail (boosting β's own precision by p, `boost_beta=TRUE`) was an unjustified side effect that inflated K_eff from 2→4 and the C-index from 0.636→0.642 with no real basis. The corrected default (`boost_beta=FALSE`) leaves YFB's fit mechanically identical to pre-Phase-1 (β values essentially unchanged); the small remaining C-index movement (0.636→0.627) is attributable entirely to Phase 1c's train/test preprocessing fix, a genuine bug fix kept regardless of its small effect on this metric. K=7's gap vs. DeSurv's K=3 is attributable to confirmed methodology differences (DeSurv jointly tunes k/α/λ via Bayesian optimization with a fixed elastic-net penalty; we tune only K via CV) rather than evidence the data itself needs more factors — a fresh K-CV confirmed K=7 is a genuine, non-noise result (`DECISIONS.md` 2026-07-12); K_eff=2 already tracks DeSurv's own ~1 survival-active factor reasonably well. **K=7 is not free to shrink in a single-seed comparison (Phase 3, `DECISIONS.md` 2026-07-13):** refitting YFB D4 directly at K∈{2,3,4,5} and re-running external validation (not just counting active factors from the K=7 fit) shows real performance loss at every smaller K tested here, on every one of the 5 held-out cohorts (best smaller-K mean C=0.596 at K=5, vs. K=7's 0.627, outside the 1-SE margin); the curve is non-monotonic and K_eff does not track external performance (K=5 has K_eff=3 but still underperforms K=7's K_eff=2). Caveat: K=2/K=4 converge suspiciously fast to near-zero β, consistent with the same CAVI factor-collapse failure mode documented for Phase 2 — those two specific numbers need a multistart re-check before being treated as a hard ceiling; K=5 (which converges normally and still falls short) makes the qualitative conclusion likely robust regardless. **K-parsimony follow-up, conclusive (Steps 1-4, `DECISIONS.md` 2026-07-13):** Phase 3's caveat was confirmed real, not resolved as hoped — K=2/K=4's underperformance was specifically the CAVI factor-collapse artifact, and warm-starting K=4 (or K=5) from the converged K=7 fit's top-K PVE-ranked columns reaches K=7-level external performance (mean C=0.6270 vs. K=7's 0.6267, both K_eff=2) — best-ELBO multistart and joint (K, α) Bayesian optimization (which instead found a marginally-better but *larger* K=8, 0.6282) independently confirm the same K_eff=2 ceiling and the same statistical tie across K∈{4,5,7,8}; K=2/K=3 remain a genuine floor under every strategy tried. **K=7 is kept as the recommended default** (not changed in `config/globals.yml`) because it's reachable via a single dependency-free fresh-SVD fit, while K=4/K=5's equivalent performance requires a two-step warm-start procedure for no measurable predictive gain — a reproducibility-simplicity judgment call, not a claim that K=7 out-performs K=4/K=5.
- **Documentation audience:** Write ROADMAP.md, DECISIONS.md, PROJECT_STATUS.qmd, and README.md for biostatistician collaborators reading the project cold — not as implementation logs. Avoid internal session terminology (e.g. "Phase A/B/C", "Cluster A/B", "Session N") in prose descriptions; use those labels only in commit messages or as lookup keys. Describe methods and findings in terms a statistical reader would recognize: model variant, prior, training set, metric, result.
- **No mannered prose:** Say what you mean; don't substitute metaphor or flourish for direct statement (e.g. "a dial worth turning" instead of "a parameter worth varying," "this point earns its keep" instead of "this point still matters"). Applies to all writing in this repo — reports, the progress book, plan/decision docs, commit messages, comments. When a literal phrase is available, use it.

## Quick Reference

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
| Cohort F update | `code/update_F_cohort.R` — `update_F_cohort_all()` (Normal conjugate) |
| Companion doc for fit_modular.R | `docs/fit_modular.qmd` |
| Global hyperparameter registry | `config/globals.yml` |
| **LB benchmark runner** | `results/benchmark_sim/run_LB_benchmark.R` |
| **YFB benchmark runner** | `results/benchmark_sim/run_YFB_benchmark.R` |
| **Progress notebook (meeting-facing, one chapter per advisor meeting)** | `docs/progress_book/` — Quarto book, `quarto render` to build; add a new `chapters/YYYY-MM-DD.qmd` per meeting |
| **Benchmark reports (dated)** | `docs/reports/ssbmf_summary_report_MM_DD_YY.{qmd,pdf,html}` — DeSurv record: `_04_29_26` |
| Phase 1 loading heatmaps | `results/benchmark_sim/run_phase1_diagnostics.R` |
| Archived benchmark runners | `results/benchmark_sim/archive/` — 9 retired scripts (see archive/README.md) |
| **β=0 design doc (LB/YFB)** | `docs/beta_zero_fix_design.md` — five-phase plan |
| **L-update debugging guide** | `docs/update_L_fix.md` — read before any L-update work |
| Exploratory simulation runner (legacy) | `results/legacy/modular_sim_factor/run_factor_modular_simulation.R` |
| Alpha CV selection | `code/select_alpha_cv.R` |
| DeSurv preprocessing | `code/preprocess_desurv.R` |
| Cohort extension benchmark | `results/benchmark_sim/run_cohort_lmm_benchmark.R` — 4-way comparison (LB/YFB × base/cohort) |
| Cohort extension report | `docs/reports/cohort_lmm_benchmark_report_5_22_26.{qmd,pdf,html}` |
| **Merged K-CV (all preprocessing)** | `results/benchmark_sim/run_merged_kcv.R` — CV K selection for all 5 preprocessing configs; K floor K≥3; fills globals.yml |
| **Merged 18-config benchmark** | `results/benchmark_sim/run_merged_benchmark.R` — 18-model fit + external validation (5 preprocessing × 2 models × ±cohort) |
| **Merged benchmark report** | `docs/reports/merged_benchmark_report_5_25_26.{qmd,pdf,html}` — 18-config results, recommended configuration |
| **DeSurv comparison runner** | `results/benchmark_sim/run_desurv_comparison.R` |
| **DeSurv comparison results** | `results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv` |
| **DeSurv comparison report** | `docs/reports/desurv_alignment_report_05_27_26.{qmd,pdf,html}` |
| **Pathway enrichment functions** | `code/pathway_enrichment.R` — `load_d4_weights()`, `run_fgsea_program()`/`run_ora_program()`, `build_pdac_genesets()`, subtype/cohort/DeSurv-overlap concordance functions |
| **Pathway enrichment runners** | `results/benchmark_sim/run_pathway_enrichment.R` (fgsea/ORA, all collections), `run_subtype_concordance.R` (PurIST), `run_external_cohort_robustness.R` (5 cohorts), `run_sbmf_desurv_overlap.R` |
| **Pathway enrichment report** | `docs/reports/pathway_enrichment_report_07_15_26.{qmd,pdf,html}` — DECISIONS.md 2026-07-15 |
| Test suite (core + predict) | `tests/run_tests.R` (412/412) |
| Real-data test suite | `tests/test_real_data_loading.R` (88/88, local-only) |
| Corrected derivations | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf` |
| Architectural decisions log | `DECISIONS.md` |
| Prioritized next steps | `ROADMAP.md` |
