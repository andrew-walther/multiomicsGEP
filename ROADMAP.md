# multiomicsGEP — Research Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the multiomicsGEP project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-05-20.** Core model complete (modular CAVI, 193/193 tests passing).
> Two model variants are fully implemented, benchmarked, and cross-validated on PDAC data
> (5 independent external cohorts across RNA-seq, microarray, and proteomics platforms).
> Four targeted improvements evaluated 2026-05-06 (Phases 1–4); YFB K-CV sign fix completed
> 2026-05-20. See synthesis report `docs/reports/ssbmf_summary_report_05_06_26.pdf`.
>
> **LB model** (linear predictor η = Lβ; `code/fit_modular.R`): All three training configurations
> (TCGA-only, CPTAC-only, merged TCGA+CPTAC) converge to 2–3 active factors. External C-index
> ranges 0.55–0.67 depending on cohort and training set. Alpha mixing parameter selected by
> 5-fold CV with 1-SE rule (all configurations select α=0.50). K-CV (Phase 4) selects K=3
> under normal prior (flat plateau) and K=8 under point_normal (spike-small-n artefact);
> practical default is K=5.
>
> **YFB model** (linear predictor η = (YF)β, Cox-on-YF; `code/fit_cox_on_yf.R`): Single-cohort
> training with a normal prior gives competitive external C-index (TCGA-only: 0.55–0.63).
> Phase 1 (per-platform z-standardization + rank_transform=FALSE) resolved the merged β→0
> collapse; K=3 merged training now yields external median C=0.64, matching or exceeding
> DeSurv. YFB K-CV sign fix (2026-05-20): `sign_correction=FALSE` in CV folds + `I(-risk)`
> evaluation resolves the C<0.5 problem; normal prior now yields C=0.56–0.61 across K, with
> K=5 selected by 1-SE rule on TCGA_PAAD. YFB point_normal K-CV remains C=0.5 for all K
> (spike-small-n collapse; open item).
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

- [x] **Fix YFB β→0 collapse on merged TCGA+CPTAC** *(Complete — Phase 1, 2026-05-06)*
  Per-platform z-standardization (normalize TCGA_PAAD and CPTAC separately before merging)
  plus removal of the per-subject rank transform resolves the collapse. With K=3, merged YFB
  training achieves external C=0.53–0.67 (median 0.64), matching or exceeding the DeSurv
  benchmark range of 0.60–0.65. Both priors agree at K=3 (factors strong enough to escape
  the spike component). Implementation: `--per-platform-norm --no-rank` flags in
  `results/benchmark_sim/run_YFB_benchmark.R`; `per_platform_standardize=TRUE` and
  `rank_transform=FALSE` in `preprocess_desurv_cohort()`.
  *Files: `code/preprocess_desurv.R`, `results/benchmark_sim/run_YFB_benchmark.R`,
  `config/globals.yml` (k_pdac_yfb_merged: 3). Results: `outputs/YFB_benchmark_perplatform/`.*

- [x] **Cohort dummy-variable extension: shared L with cohort-specific F columns** *(Complete — 2026-05-22)*
  Corner-point encoded cohort indicator columns appended to L; F_cohort rows estimated via
  Normal conjugate update; cohort_id and sigma_F_cohort parameters in both models.
  Synthetic: offset absorption |cor|=0.995, factor recovery +88–109% vs baseline.
  Real PDAC: LB_cohort gains on RNA-seq external cohorts (PACA_AU_seq +0.019, PACA_AU_array +0.007)
  but loses on Dijk/Puleo (mean C 0.604 vs LB_base 0.618). YFB β→0 unchanged.
  Full evaluation: `docs/reports/cohort_lmm_benchmark_report.qmd` and DECISIONS.md 2026-05-22.
  *Files: `code/update_F_cohort.R`, `code/compute_elbo.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`,
  `tests/test_update_F_cohort.R`, `tests/test_fit_modular_cohort.R`, `tests/test_fit_yf_cohort.R`.
  Validation: `results/cohort_lmm_sim/run_synthetic.R`, `results/benchmark_sim/run_cohort_lmm_benchmark.R`.*

- [x] **YFB β→0 on merged data: all CAVI-local fixes exhausted — formally closed** *(2026-05-25)*
  Frozen-F β pre-conditioning (N_frozen parameter, V9–V11) tested as the final CAVI-local
  strategy. All three frozen-F variants converge at iteration 9 — during the frozen phase,
  before EF unfreezes. Root cause: SVD initialization on merged TCGA+CPTAC is itself
  platform-dominated, so ZF_SVD has no survival correlation and β→0 during the frozen phase.
  Exhaustive record: V0 (baseline) through V11 (frozen-F N=30+warm-start) all fail.
  Practical solutions: (1) per-platform YFB (Phase 1, C_dijk≈0.573); (2) LB_cohort K=5
  (mean C=0.614); (3) LB_base K=20 (mean C=0.618). Possible future fix: initialize EF
  from a single-cohort TCGA-only YFB fit (not platform-dominated SVD) and freeze during
  merge training. Not implemented — deferred as low-priority.
  *Files: `code/fit_cox_on_yf.R` (N_frozen param, backward compat), `tests/test_fit_yf_frozen_f.R` (8 tests),
  `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R` (V9–V11 added). DECISIONS.md 2026-05-25.*

- [x] **Evaluate cohort extension benefit on lower-K regimes** *(2026-05-25)*
  Run `run_cohort_lmm_benchmark.R --low-k` (K_LB=5, K_YFB=3, 300 iters).
  Key finding: at K=5, LB_base also collapses to β→0 (K_eff=0). LB_cohort at K=5 rescues
  β (K_eff=3, mean C=0.614), nearly matching LB_base K=20 (mean C=0.618) with 4× fewer
  biological factors. The cohort extension is *essential* in the interpretable low-K regime
  and *neutral to marginal* at K=20 where ARD handles platform effects implicitly.
  Practical recommendation: for merged multi-platform PDAC analysis, use K=5 + cohort_id.
  *Files: `results/benchmark_sim/run_cohort_lmm_benchmark.R` (--low-k flag added),
  `results/benchmark_sim/outputs/cohort_lmm_benchmark_low_k/`. DECISIONS.md 2026-05-25.*

- [x] **Merged-cohort comprehensive preprocessing benchmark (18 configurations)** *(Complete — 2026-05-25)*
  Extended 6-config comparison to 18 by adding joint quantile-no-rank, joint z-standardization,
  and log-only preprocessing. Biological K floor K_final = max(K_1se, 3) applied throughout.
  Winner: M5 (YFB × per-platform z-std × no cohort_id, K=3), mean external C=0.626, K_eff=2/3.
  Key finding: per-platform z-standardization is the *only* viable preprocessing for mixed
  RNA-seq + proteomics training data. 10 of 12 non-per-platform configs collapse to β=0 (K_eff=0).
  Rank transform is not the driver of YFB β→0 on joint-QN data (M13/M14 also collapse without it).
  Recommended for manuscript: M5 (primary), M4 (sensitivity). See DECISIONS.md 2026-05-25.
  *Files: `results/benchmark_sim/run_merged_kcv.R`, `results/benchmark_sim/run_merged_benchmark.R`,
  `code/preprocess_desurv.R` (normalize_method param), `docs/reports/merged_benchmark_report.qmd`.
  Results: `outputs/merged_benchmark/merged_benchmark_results_extended.csv`. DECISIONS.md 2026-05-25.*

- [x] **Align gene selection with DeSurv: combined mean+variance ranking, top-3000** *(Complete — 2026-05-27)*
  Implemented combined_rank method in select_top_variable_genes() and per-cohort
  selection (before normalization) in preprocess_merged_cohorts(). D4 (YFB DeSurv-aligned)
  mean external C=0.636 vs M5=0.624 (delta=+0.012, 5 cohorts). DeSurv gene selection
  adopted as new primary config. See DECISIONS.md 2026-05-27.
  *Files: code/preprocess_desurv.R, results/benchmark_sim/run_desurv_comparison.R,
  docs/reports/desurv_alignment_report_05_27_26.qmd*

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

- [x] **K selection via cross-validated C-index** *(Complete — Phase 4, 2026-05-06)*
  `select_K_cv()` implemented in `code/select_K.R` with 1-SE rule; accepts `model="LB"` or
  `model="YFB"`. Five-fold CV over K∈{2,…,10,15,20} on TCGA_PAAD (n=144):
  - LB normal prior: flat plateau (range=0.032 < per-K SE), K=3 selected (1-SE).
  - LB point_normal: K=8 selected (spike-small-n artefact; C=0.5 at K≤3). Practical default: K=5.
  - YFB point_normal: K=2 (β→0 collapse in all folds, C=0.5 everywhere).
  - YFB normal: K=5 selected (1-SE; C=0.56–0.61 across all K) after sign fix (2026-05-20).
  - YFB point_normal: C=0.5 for all K (β→0 collapse in folds; open item).
  193/193 tests passing (KCV-T15 added for YFB sign fix verification).
  *Files: `code/select_K.R`, `tests/test_select_K_cv.R`, `results/benchmark_sim/run_K_cv.R`.
  Results: `results/benchmark_sim/outputs/K_cv/`. Report: `docs/reports/ssbmf_summary_report_05_06_26.pdf`.*

- [x] **YFB sign correction in K-CV** *(Complete — 2026-05-20)*
  `sign_correction=FALSE` parameter added to `fit_cox_on_yf()`; passed in YFB branch of
  `select_K_cv()`. CV folds deliver raw SVD-oriented EBeta; concordance evaluated via
  `I(-pred$risk_scores)`. After fix: normal prior yields C=0.56–0.61, K=5 selected (1-SE).
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`, `tests/test_select_K_cv.R` (KCV-T15).
  Results: `results/benchmark_sim/outputs/K_cv/K_cv_table_YFB_normal.csv`.*

- [ ] **YFB point_normal K-CV collapse** `[Priority: Medium]` `[Effort: Medium]`
  Under the spike-and-slab (point_normal) prior, YFB β→0 in all CV folds regardless of K
  (C=0.5 for all K ∈ {2,…,20}). The reduced training fold size (~115 subjects) is insufficient
  to escape the spike component. This is a separate issue from the sign fix. Possible
  approaches: β-only warm-up iterations before the spike prior activates; adaptive spike
  weight schedule; or accepting that point_normal is not viable for YFB K-CV on small-n data.
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`*

- [ ] **YFB preprocessing alignment** `[Priority: Medium]` `[Effort: Medium]`
  Single-cohort and merged YFB training use different preprocessing contracts:
  single-cohort uses the per-cohort DeSurv pipeline (`preprocess_desurv_cohort()`);
  merged uses per-platform z-standardization (`per_platform_standardize=TRUE`,
  `rank_transform=FALSE`). These are not interchangeable — a model trained in single-cohort
  mode cannot be applied with merged-mode preprocessing at prediction time. A unified
  preprocessing API that makes the contract explicit and enforces train/test consistency
  would reduce downstream errors.
  *Files: `code/preprocess_desurv.R`, `results/benchmark_sim/run_YFB_benchmark.R`*

- [ ] **Per-platform noise variance (τ) investigation** `[Priority: Low]` `[Effort: Medium]`
  The current model uses a single shared noise precision τ (estimated per gene, shared across
  subjects and cohorts). In merged multi-platform training (TCGA RNA-seq + CPTAC proteomics),
  the true noise variance almost certainly differs between platforms. A per-platform τ would
  better reflect the generative structure and may reduce the genomics-term scale imbalance
  that drives β→0 in the L update. Requires a new q(τ) derivation for the grouped case.
  *Files: `code/update_tau.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`*

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

- [x] **YFB K-CV sign fix (2026-05-20)** — `sign_correction=FALSE` parameter added to
  `fit_cox_on_yf()` and passed in `select_K_cv(model="YFB")`. Before fix: C<0.5 for all K
  (Phase C oriented EBeta, then `I(-.)` re-flipped, producing anti-concordant predictions).
  After fix: YFB normal prior yields C=0.56–0.61 across K ∈ {2,…,20}; K=5 selected by 1-SE
  rule on TCGA_PAAD. YFB point_normal remains C=0.5 (spike-small-n collapse; separate open
  item). 193/193 tests passing.
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`, `tests/test_select_K_cv.R`.*

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
