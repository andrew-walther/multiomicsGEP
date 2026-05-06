# multiomicsGEP — Research Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the multiomicsGEP project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-05-05.** Core model complete (modular CAVI, 171/171 tests passing).
> Two model variants are fully implemented and benchmarked on PDAC cross-cohort validation
> (5 independent external cohorts across RNA-seq, microarray, and proteomics platforms).
>
> **LB model** (linear predictor η = Lβ; `code/fit_modular.R`): All three training configurations
> (TCGA-only, CPTAC-only, merged TCGA+CPTAC) converge to 2–3 active factors. External C-index
> ranges 0.55–0.67 depending on cohort and training set. Alpha mixing parameter selected by
> 5-fold CV with 1-SE rule; all configurations select α=0.50.
>
> **YFB model** (linear predictor η = (YF)β, Cox-on-YF; `code/fit_cox_on_yf.R`): Single-cohort
> training with a normal prior gives competitive external C-index (TCGA-only: 0.55–0.63). Merged
> training (TCGA_PAAD + CPTAC) collapses β→0 due to platform mixing — the factor matrix F still
> captures survival-relevant structure, yielding C=0.57–0.64 with near-zero β. Point-normal prior
> collapses all betas regardless of training set; normal prior is required for YFB on real PDAC data.
>
> Both models use SVD initialization with a post-convergence sign check: if the training
> concordance of η = Lβ (or η = (YF)β) is below 0.5, β is globally negated. This corrects
> sign ambiguity from the SVD initialization, which uses only the positive part of the singular
> vectors and can produce an inversely-oriented linear predictor.

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

- [x] **Re-run benchmarks with top_n=2000 to restore single-cohort baseline** *(complete 2026-05-05)*
  `top_n_genes` was inadvertently raised to 5000. Reverted to 2000. K=20 on single-cohort
  data (tcga_only n=144) still gave K_eff=1 and C=0.37–0.50 — K overfitting, not top_n.
  Fixed by adding `k_pdac_single=10` and re-running. See DECISIONS.md 2026-05-05.

- [x] **Restore alpha CV to benchmark runners** *(complete 2026-05-05)*
  `select_alpha_cv()` added to `run_LB_benchmark.R`. The sign correction must not run inside
  CV folds — doing so causes fold-to-fold sign inconsistency that inflates apparent concordance
  and causes the 1-SE rule to select α=1.0 (degenerate pure-survival mode). Fix: set
  `sign_correction=FALSE` in CV fold fits. All three LB training modes select α=0.50 via 1-SE
  rule. YFB runner uses fixed α=0.50 (YFB-compatible alpha CV not yet implemented).

- [x] **Fix K overfitting and re-run benchmarks** *(complete 2026-05-05)*
  `k_pdac_single=10` added; all 6 modes re-run. **K tuning did not recover archived baseline.**
  LB tcga_only K=10: K_eff=2, external C=0.34–0.43 (worse than K=20 which gave 0.47–0.50).
  More active factors → more anti-concordant external predictions. The archived 0.63–0.65 was
  a lucky PCA direction alignment, not a stable property. Root cause is A_surv/A_gen structural
  imbalance — must be addressed directly. See DECISIONS.md 2026-05-05.

- [ ] **Fix YFB β→0 collapse on merged TCGA+CPTAC** `[Priority: Medium]` `[Effort: Medium]`
  The LB model converges to 3 active factors on merged training (external C=0.51–0.67), but
  the YFB model collapses β→0 (machine-epsilon betas ~10⁻¹³) on the same data. The root cause
  is platform mixing: merged training combines RNA-seq (TCGA_PAAD) and proteomics (CPTAC), and
  the ratio of survival precision to genomics precision in the β CAVI update is ~10⁻³ at
  iteration 1 — the genomics term swamps the survival signal in the factor scores. The factor
  matrix F still captures survival-relevant column structure, yielding external C=0.57–0.64
  despite near-zero β. Cox warm-start (initializing β from a Cox fit on YF scores) was tested
  and did not rescue collapse on merged data.
  Candidate next steps:
  (1) **Rescale survival and genomics precision contributions** to comparable magnitude before
  the β update — analogous to `normalize_AB` in the LB model. Evaluated on LB but regressed
  4/5 external cohorts; may behave differently for YFB.
  (2) **Per-platform normalization** — normalize TCGA and CPTAC separately before merging,
  rather than relying on quantile normalization across mixed platforms.
  (3) **Accept single-cohort training** — LB and YFB single-cohort results are already
  competitive; merged training may not add value given the platform heterogeneity.
  *Files: `code/fit_cox_on_yf.R`, `code/preprocess_desurv.R`*

- [ ] **Prior comparison follow-up** `[Priority: Medium]` `[Effort: Small]`
  Current benchmarks test point_normal vs normal for both models. Key finding: for the YFB model
  on real PDAC data, the normal prior is strictly necessary — the spike-and-slab (point_normal)
  collapses all β to zero because the survival signal is too weak to exceed the spike threshold.
  For the LB model, both priors give identical results at the current signal level. A fuller
  comparison should include point_laplace and evaluate whether a β-only warm-up (fitting β while
  holding L fixed) helps the spike-and-slab prior find the active set.
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

- [x] **Investigate and resolve β=0 collapse in LB model** *(Complete — 2026-04-29)*
  Instrumentation confirmed a ~5000× scale imbalance between survival and genomics precision
  contributions to the L update at iteration 1 (A_surv/A_gen ~ 0–2e-4). Three fixes were
  implemented and evaluated: (1) Gauss-Seidel update ordering (β before L in the inner loop);
  (2) β-only warm-up iterations with L held fixed; (3) rescaling survival and genomics
  precision contributions to comparable magnitude. On merged TCGA+CPTAC: 2/20 factors active,
  ELBO monotone. Training-side β=0 resolved. External C-index was mixed (1/5 cohorts improved),
  motivating the Cox-on-YF reformulation below.
  *Files: `code/fit_modular.R`, `code/update_L.R`*

- [x] **Cox-on-YF (YFB) reformulation** *(Complete — 2026-05-04, merged to `main`)*
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

- [x] **Training concordance sign correction + alpha CV fix (2026-05-05)** — Both LB and YFB
  were producing anti-concordant predictions (C < 0.5) due to sign ambiguity in SVD initialization:
  the positive-part truncation of the SVD can produce factor loadings whose inner product with β
  is inversely related to survival risk. Fix: after CAVI convergence, compute training concordance;
  if C_train < 0.5, negate β globally. Implemented as `sign_correction=TRUE` parameter in both
  `fit_supervised_mf_modular()` (LB) and `fit_cox_on_yf()` (YFB).
  A related bug in the alpha CV routine was fixed in the same pass: applying the sign correction
  inside CV folds caused fold-to-fold inconsistency, making the 1-SE rule select α=1.0 (pure-
  survival, degenerate). Fix: `sign_correction=FALSE` in fold fits.
  YFB synthetic result: C-index 0.12 → 0.91 (K_eff=4). LB external results (all modes, α=0.50):
  TCGA-only K_eff=2, C=0.55–0.66; CPTAC-only K_eff=3, C=0.55–0.67; merged K_eff=3, C=0.51–0.67.
  *Files: `code/fit_modular.R`, `code/select_alpha_cv.R`, `code/fit_cox_on_yf.R`,
  `code/predict_cox_on_yf.R`. Report: `docs/reports/ssbmf_summary_report_05_05_26.pdf`.*

- [x] **Benchmark pipeline consolidation (2026-05-04)** — Replaced scattered diagnostic scripts
  with two canonical benchmark runners: `run_LB_benchmark.R` (LB model, η = Lβ) and
  `run_YFB_benchmark.R` (YFB model, η = (YF)β). Both follow a 4-section structure (synthetic
  validation + PDAC training + external validation + CSV output) and test point_normal vs normal
  prior side-by-side. Archived 7 retired scripts to `results/benchmark_sim/archive/`. Updated
  `config/globals.yml` with benchmark defaults (K=10, α=0.5, λ=1.0, beta_threshold=0.001).

- [x] **Cox-on-YF (YFB) model implementation (2026-05-04)** — Reformulated linear predictor to
  η = (Y·F)·β, replacing η = L·β. Training and prediction both use the observed expression matrix
  projected onto the factor space, eliminating the train/test mismatch where training used
  EBNM-shrunk loadings but prediction used OLS-projected loadings. Implementation: `code/fit_cox_on_yf.R`,
  `code/predict_cox_on_yf.R`, `code/update_L_surv_YFB.R`, `code/update_F_surv_YFB.R`,
  derivations in `derivations/cox_on_YF/`. Synthetic: C-index 0.605 vs PCA 0.471 (3 active factors).
  Full PDAC benchmark comparison in `docs/reports/ssbmf_summary_report_05_05_26.pdf`.
- [x] **Core modular CAVI implementation** — `code/fit_modular.R` with four independently-tested update modules. 171/171 tests passing. *(Completed March 2026)*
- [x] **Hold-out prediction pipeline** — `code/predict.R` (`predict_supervised_mf()`), `code/train_test_split.R` (`stratified_split()`). 80/20 stratified hold-out. *(Completed March 2026)*
- [x] **Prior family comparison (PN vs PL × K=5 vs K_eff)** — Four-condition benchmark across 7 PDAC cohorts. Point-laplace preferred at fixed K; K selection dominates prior choice. *(Completed April 2026)*
- [x] **Full ELBO tracking** — Both proxy and full ELBO (genomics + survival + KL) tracked per iteration. `code/compute_elbo.R`. *(Completed April 2026)*
- [x] **DeSurv benchmark — synthetic + PDAC cross-cohort** — DeSurv-aligned preprocessing, alpha CV (1-SE rule), SVD pseudoinverse fix, prior sensitivity (point_normal vs point_laplace). Synthetic: supervised 0.79 > PCA 0.76. PDAC external median C-index 0.60 across 5 cohorts. *(Completed April 2026)*
- [x] **Prior sensitivity report** — point_normal vs point_laplace compared on synthetic and all PDAC training modes. point_normal recommended as default. *(Completed April 2026)*
- [x] **24-page DeSurv benchmark report** — `results/benchmark_sim/ssbmf_summary_report.pdf`. Includes ARD justification over ELBO K-grid-search, alpha gradient notation, per-figure takeaways, side-by-side prior comparison (Fig 19), multi-modal failure documented. *(Completed April 2026)*
- [x] **C-index honest reporting fix** — `get_cindex_comparison()` uses model's own `EL %*% EBeta`; corrected `concordance()` direction convention. *(Completed April 2026)*
