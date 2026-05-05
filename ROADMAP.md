# multiomicsGEP — Research Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the multiomicsGEP project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-05-05:** Core model complete (modular CAVI, 171/171 tests passing).
> Benchmark runners (`run_LB_benchmark.R`, `run_YFB_benchmark.R`) now support `--train-mode
> merged|tcga_only|cptac_only`. All 6 full benchmark runs completed. **Key finding:** merged
> training has never produced C-index > 0.50 (β→0 collapse). The archived 0.60–0.65 results
> were from tcga_only training (v1 preprocessing, K=10, alpha CV) — not merged. Both LB and
> YFB collapse β→0 on merged data. YFB collapses on single-cohort too. The β→0 collapse on
> merged TCGA+CPTAC is the core unresolved problem. **Immediate priorities below.**

---

## How to Use This Document

Each item follows this format:

```
- [ ] **Title** `[Priority: High/Medium/Low]` `[Effort: Small/Medium/Large]`
  Description of what this involves and why it matters.
  *Notes: dependencies, caveats, or relevant files*
```

**Priority:** How much this would strengthen the work (High = strengthens core contribution; Medium = meaningful addition; Low = nice-to-have)

**Effort:** Rough implementation cost (Small = hours; Medium = days; Large = weeks or HPC run)

Move completed items to the [Completed](#-completed) section at the bottom.

---

## 🔥 Immediate Priorities

- [ ] **Re-run benchmarks with top_n=2000 to restore single-cohort baseline** `[Priority: High]` `[Effort: Small]`
  `top_n_genes` was inadvertently raised to 5000. The archived 0.60–0.65 results used the DeSurv
  spec of 2000 genes. Reverted in this session; need a clean run to confirm tcga_only/cptac_only
  recover to C ≈ 0.60 with normal prior. Merged training expected to remain ~0.50 regardless.
  *Files: `results/benchmark_sim/run_LB_benchmark.R`, `results/benchmark_sim/run_YFB_benchmark.R`*

- [ ] **Restore alpha CV to benchmark runners** `[Priority: High]` `[Effort: Small]`
  Archived strong baseline used 5-fold CV to select alpha (0.5 for single-cohort, 0.3 for merged).
  Current runners use fixed `alpha=0.5`. Re-adding `select_alpha_cv()` is low effort and ensures
  merged training isn't penalized by a suboptimal fixed alpha.
  *Files: `results/benchmark_sim/run_LB_benchmark.R`, `code/select_alpha_cv.R`*

- [ ] **Diagnose and fix β→0 collapse on merged TCGA+CPTAC** `[Priority: High]` `[Effort: Large]`
  Core unresolved problem. Both LB and YFB collapse β→0 on merged data regardless of prior, K,
  or preprocessing. Cluster A resolved training-side β=0 but external generalization regressed
  (4/5 cohorts). Cluster B (YFB) also collapses on all train modes including single-cohort.
  Candidate next steps:
  (1) **Instrument YFB A_surv/A_gen ratio** — add the `[iter1, k=?] ratio` printout from Cluster A
  to `fit_cox_on_yf.R` and compare to LB on merged data. Determines if YFB has the same structural
  imbalance, which would mean the reformulation didn't fix the root cause.
  (2) **Cox warm-start for YFB** (`cox_warmstart=TRUE`) — initializes EBeta from a Cox fit on YF
  scores; may prevent early β collapse if CAVI starts near a non-zero solution.
  (3) **Understand single-cohort YFB collapse** — LB gets K_eff=2–4 on tcga_only but YFB collapses
  there too. The problem is in the YF reformulation itself, not just merged-data batch effects.
  The YF matrix scale (ZF = Y·F, sd ~ 10–100×) shrinks effective beta to ~0.003–0.008, making
  even the normal prior place most mass near zero.
  *Files: `code/fit_cox_on_yf.R`, `docs/beta_zero_fix_design.md`, `docs/update_L_fix.md`*

- [ ] **Prior comparison follow-up** `[Priority: Medium]` `[Effort: Small]`
  Once β→0 collapse is fixed, compare point_normal vs normal vs point_laplace on the recovered
  baseline. N_burnin=5 is worth re-testing with normal prior (Cluster A showed modest benefit).
  *Files: `results/benchmark_sim/run_LB_benchmark.R`, `results/benchmark_sim/run_YFB_benchmark.R`*

---

## 🧬 Model Specification

- [ ] **Incorporate F'Y term in survival model** `[Priority: Medium]` `[Effort: Large]`
  Replace $L_i$ in the Cox linear predictor with $(F^\top y_i)$ — the projection of patient $i$'s
  observed expression onto the factor space. This decouples factor learning from survival supervision:
  L is learned purely from genomics, and survival scoring is done post-hoc via the factor projection.
  Requires a new q(L) CAVI derivation (the survival term drops out of the L update) and a post-hoc
  scoring step. May improve identifiability and reduce the sign-flip issue observed in synthetic
  validation.
  *Notes: High effort — requires new derivation and update to `code/update_L.R`. See `derivations/qL/qL_update_derivation.pdf` for current L update derivation.*

- [x] **EBMF warm-start initialization for SSBMF** `[Priority: High]` `[Effort: Medium]` *(Complete — 2026-04-29)*
  Initialize SSBMF L and F from an EBMF solution and run CAVI from that starting point.
  Two experiments run: (1) β-only with EL fixed → β non-zero at iter 1, 6/20 factors active,
  confirms β update is functional. (2) Full CAVI warm-start → β collapses to zero in 23 iters,
  L/F updates wash out EBMF structure. Root cause: `update_L_k()` A_surv dominated by A_gen.
  See DECISIONS.md 2026-04-29 warm-start entry.
  *Notes: `fit_modular.R` extended with `EL_init`/`EF_init` params. Driver: `run_ebmf_warmstart.R`.*

- [x] **Debug `update_L_k()`: A_surv / A_gen imbalance** `[Priority: High]` `[Effort: Medium]` *(Complete — Cluster A, 2026-04-29)*
  Cluster A (`docs/beta_zero_fix_design.md` §4) implemented on branch `fix-L-update-beta-cycle`:
  instrumentation confirmed A_surv/A_gen ~ 0–2e-4 at iter 1 (structural imbalance, ~5000× gap);
  inner-loop reorder β → L → F (unconditional); β-only burn-in (`N_burnin`) and progressive α
  schedule (`alpha_schedule`) added as opt-in parameters; `normalize_AB` rescale of A_surv up to
  match A_gen (Fix 4, reformulated from design doc §4.8 — original formula over-shrunk L). On the
  merged TCGA+CPTAC v2 set: 2/20 factors active (|β|>0.05), ELBO monotone, max|EL| = 1.83e3.
  Training-side β=0 failure resolved.
  *Notes: External-cohort C-index is mixed (1/5 cohorts improved vs. baseline; 4/5 regressed).
  The recovered β favors training-Cox-aligned directions that don't transport. See new entry
  "Cluster B / Cox-on-YF reformulation" below — this is the design doc's structural alternative.*

- [x] **Cluster B — Cox-on-YF reformulation** *(Complete — 2026-05-04, merged to `main`)*
  See completed section below.

- [x] **Add λ scaling parameter to balance genomics vs. survival objectives** `[Priority: High]` `[Effort: Medium]` *(Implemented and evaluated — fixed at λ=1.0)*
  λ is implemented as an exposed parameter in `update_L_k()`, `update_L_all()`, and
  `fit_supervised_mf_modular()` (default 1.0), and registered in `config/globals.yml`.
  A controlled sandbox (n=250, p=1000, K=5) tested λ∈{1, p/n=5, 2p/n=10}: hold-out C-index
  was flat at ≈0.805 across all conditions; β RMSE was *worse* at λ=p/n (+0.25) and λ=2p/n (+0.43).
  The dominant β scale error is L–β scale indeterminacy, not gradient imbalance. λ=1 retained.
  *Notes: To experiment, change `cavi.lambda` in `config/globals.yml`. See `DECISIONS.md` entry 2026-04-24 for full analysis. Sandbox: `results/benchmark_sim/sandbox_lambda_test.R`.*

---

## 🌍 Validation & Generalization

- [x] **Cross-cohort validation: train on TCGA-only or CPTAC-only, evaluate on 5 external cohorts** `[Priority: High]` `[Effort: Medium]` *(Complete — see `results/benchmark_sim/outputs/real_data/`)*
  Merge CPTAC and TCGA_PAAD (batch-corrected) as a training set and hold out Dijk,
  Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, and Puleo_array as independent test cohorts.
  This tests generalisation of learned gene expression programs across platforms (RNA-seq,
  microarray, proteomics) and patient populations. More training data than the current 80/20
  within-cohort split. Requires R4 (inter-dataset normalisation) as a prerequisite.
  *Notes: `code/train_test_split.R` and `code/predict.R` are already in place for the projection step. Batch correction (ComBat or similar) not yet implemented — add to R4. See `results/modular_sim_factor/PDAC/` for current within-cohort results.*

- [x] **Inter-dataset normalisation before cohort merging** `[Priority: High]` `[Effort: Medium]` *(Complete — v2 preprocessing, 2026-04-29)*
  Before merging cohorts (required for R3), apply quantile normalisation or z-score
  standardisation across datasets beyond the current per-cohort column-centring. Addresses
  platform-specific mean/variance shifts between RNA-seq, microarray, and proteomics. Essential
  for any analysis that combines data across assay types.
  *Notes: Implemented as `preprocess_merged_cohorts()` in `code/preprocess_desurv.R`. Pipeline: intersect raw gene universes → log₂(x+1) [RNA-seq only] → `preprocessCore::normalize.quantiles()` across all merged samples → top-2000 by merged-matrix variance → per-subject rank transform. Gated by `preprocessing_version = "v2"` in `run_real_data_benchmark()`. See DECISIONS.md 2026-04-29 entry.*

---

## 📐 Model Selection

- [ ] **K selection via cross-validation** `[Priority: High]` `[Effort: Medium]`
  ARD pruning (PVE > 1% OR |β| > 0.001) is fragile: K=10 collapses betas to 0 on merged PDAC
  while K=20 does not — the threshold-based criterion gives no principled way to choose K_max.
  Replace with CV over K ∈ {2, 4, 6, 8, 10, 15, 20}: fit each K, evaluate on a held-out 20%
  split, pick the K that maximizes the criterion.

  **Preferred criterion: held-out genomic reconstruction MSE** (‖Y_test − L_test·EF'‖²).
  This selects K for factorization quality alone, then EBNM beta shrinkage handles which of
  those K factors are survival-relevant. Cleaner separation of concerns than optimizing C-index
  over K (which conflates factorization quality with prior choice and is noisy on small cohorts).
  Held-out C-index over K is a secondary metric to report but not the primary selection criterion.

  The `select_K_cv()` stub exists in `code/select_K.R` — implement it.
  *Notes: 5-fold CV over 7 K values = 35 fits; feasible locally or trivially parallelizable on
  Longleaf. Do not implement until prior comparison benchmarks (LB/YFB at K=20) are complete.*

- [ ] **Move `load_pdac_raw` and `plot_cohort_loading_heatmap` to `code/`** `[Priority: Low]` `[Effort: Small]`
  Both functions currently live in `results/benchmark_sim/run_ssbmf_benchmark.R` (the data-loading
  hub). They are general enough to belong in `code/` for reuse across runner scripts without
  sourcing the full hub. Requires updating all `source()` calls in runner scripts.
  *Notes: Do in a dedicated commit; coordinate with any active Longleaf paths.*

---

## ⚙️ Infrastructure

- [ ] **Investigate {targets} workflow for pipeline reproducibility** `[Priority: Low]` `[Effort: Medium]`
  Adopt the R `targets` package to manage computational dependencies across the pipeline (data
  loading → fitting → tables → figures → reports). Tracks which outputs are stale and re-runs only
  what is needed. Particularly valuable as the number of cohorts, conditions, and reports grows —
  currently all cohorts are re-run together even when only one has changed. Also facilitates
  Longleaf HPC submission of individual targets as SLURM jobs.
  *Notes: Low priority — investigate feasibility and overhead before committing to a refactor. The `{targets}` package works best when the pipeline DAG is well-defined; this is approximately true for the current runner structure. See `results/modular_sim_factor/run_factor_modular_simulation.R` for the current monolithic runner.*

- [x] **`results/` directory cleanup** `[Priority: High]` `[Effort: Small]` *(complete 2026-05-05)*
  The `results/` directory now mixes code scripts, CSV tables, PNG/PDF figures, and
  QMD/PDF/HTML reports within the same subdirectories, making it hard to navigate and
  prone to path errors. Proposed split (within `results/benchmark_sim/`):
  ```
  results/benchmark_sim/
    scripts/        ← MOVE run_ssbmf_benchmark.R, compute_ph_diagnostics.R
    reports/        ← MOVE ssbmf_summary_report.qmd/.pdf/.html
    outputs/        ← benchmark CSVs + figures (already structured by mode/prior)
  results/modular_sim_factor/
    scripts/        ← MOVE run_factor_modular_simulation.R, run_prior_k_comparison.R
    reports/        ← MOVE .qmd/.pdf/.html out of PDAC/ and synthetic/ subdirs
  ```
  Update all `source()` calls, `quarto render` paths, and README references when moving.
  *Notes: Dedicate a single refactor commit. Don't split across sessions. Update CLAUDE.md
  quick-reference table after moving. Coordinate with any active Longleaf HPC paths.*

- [x] **Repository reorganisation** `[Priority: Low]` `[Effort: Small]` *(results/ portion complete 2026-05-05; code/ and derivations/ portions deferred — low priority)*
  The current directory layout has accumulated structural debt across three simulation generations. 
  Need to clean up the results/ directory as multiple stages of tables/figures/reports are floating 
  around without clear structure Proposed clean structure (no file deletion — move and rename only):

  ```
  code/                     ← algorithm only (update_*.R, fit_modular.R, etc.)
    legacy/                 ← V1, V2 monolithic (already there); ADD Supervised_Bayesian_MF_V2.R
    demos/                  ← MOVE top-level demos/ here (executable examples per module)
  docs/                     ← companion .qmd/.pdf/.html (already clean)
    context/                ← MOVE code/SupervisedMF_Context.md, docs/PDAC_data_audit.qmd here
  derivations/              ← mathematical derivations (already clean per-update)
    archive/                ← MOVE derivations/EBMF/, derivations/SurvivalMF/ early sketches
                               MOVE dated .tex drafts from MF_UpdateDerivations/ (keep only REVISED + Companion)
  results/
    reports/                ← MOVE .qmd/.pdf/.html report files out of modular_sim_factor/
    figures/                ← figures by cohort (already structured)
    tables/                 ← tables by cohort (already structured)
    legacy/                 ← MOVE results/full_sim/, results/modular_sim_block/ here
  tests/                    ← test suite (already clean)
  presentation/             ← slide decks (already clean)
  longleaf_setup/           ← HPC scripts (already clean)
  paper/                    ← manuscript (already clean)
  config/                   ← globals.yml (added in Session 15)
  ```

  Additional housekeeping: add `.obsidian/`, `.claire/`, `.claude/worktrees/` to `.gitignore`
  (IDE/tool artifacts should not be tracked). Remove `.Rhistory` files scattered in subdirectories.
  *Notes: Propose only — do not move files until a dedicated refactor commit. Flag as pending in `DECISIONS.md`. Coordinate with any active branches before moving paths that appear in runner scripts or `.qmd` files.*

---

## ✅ Completed

- [x] **Benchmark consolidation (2026-05-04)** — Merged `cox-on-yf-reformulation` to `main`.
  Created `run_LB_benchmark.R` (Cluster A, 4-section: synthetic + PDAC train + external + CSV)
  and `run_YFB_benchmark.R` (Cluster B, same structure). Both runners test point_normal vs normal
  side-by-side. Archived 7 one-off diagnostic scripts to `results/benchmark_sim/archive/`.
  Updated `config/globals.yml` with `benchmark` section (K=10, N_burnin=0, normalize_AB=FALSE,
  cox_warmstart=FALSE, alpha=0.5, lambda=1.0) and lowered `beta_threshold` from 0.05 to 0.001.
  Updated `code/fit_cox_on_yf.R` defaults: prior_beta="normal", N_burnin=0, cox_warmstart=FALSE,
  normalize_AB=FALSE.

- [x] **Cluster B — Cox-on-YF reformulation** *(Completed 2026-05-04)* — `code/fit_cox_on_yf.R`,
  `code/predict_cox_on_yf.R`, `code/update_L_surv_YFB.R`, `code/update_F_surv_YFB.R`,
  derivations in `derivations/cox_on_YF/`, smoke test 3/3 PASS. Synthetic: C-index 0.605 vs
  PCA 0.471 (3 active factors, alpha_F=0). Real-data β=0 collapse on PDAC with point_normal prior
  (same pattern as Cluster A). Normal prior gives non-zero EBeta — full benchmark comparison pending.
- [x] **Core modular CAVI implementation** — `code/fit_modular.R` with four independently-tested update modules. 171/171 tests passing. *(Completed March 2026)*
- [x] **Hold-out prediction pipeline** — `code/predict.R` (`predict_supervised_mf()`), `code/train_test_split.R` (`stratified_split()`). 80/20 stratified hold-out. *(Completed March 2026)*
- [x] **Prior family comparison (PN vs PL × K=5 vs K_eff)** — Four-condition benchmark across 7 PDAC cohorts. Point-laplace preferred at fixed K; K selection dominates prior choice. *(Completed April 2026)*
- [x] **Full ELBO tracking** — Both proxy and full ELBO (genomics + survival + KL) tracked per iteration. `code/compute_elbo.R`. *(Completed April 2026)*
- [x] **DeSurv benchmark — synthetic + PDAC cross-cohort** — DeSurv-aligned preprocessing, alpha CV (1-SE rule), SVD pseudoinverse fix, prior sensitivity (point_normal vs point_laplace). Synthetic: supervised 0.79 > PCA 0.76. PDAC external median C-index 0.60 across 5 cohorts. *(Completed April 2026)*
- [x] **Prior sensitivity report** — point_normal vs point_laplace compared on synthetic and all PDAC training modes. point_normal recommended as default. *(Completed April 2026)*
- [x] **24-page DeSurv benchmark report** — `results/benchmark_sim/ssbmf_summary_report.pdf`. Includes ARD justification over ELBO K-grid-search, alpha gradient notation, per-figure takeaways, side-by-side prior comparison (Fig 19), multi-modal failure documented. *(Completed April 2026)*
- [x] **C-index honest reporting fix** — `get_cindex_comparison()` uses model's own `EL %*% EBeta`; corrected `concordance()` direction convention. *(Completed April 2026)*
