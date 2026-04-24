# multiomicsGEP — Research Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the multiomicsGEP project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-04-24:** Core model implementation complete (modular CAVI, 171/171 tests
> passing). DeSurv-aligned benchmark (Phases 0–3B) implemented and running. Synthetic validation
> shows supervised C-index (0.79) > PCA (0.76) after DGP fix. Real-data benchmark (TCGA+CPTAC
> training, 5 external cohorts) in progress. See `results/benchmark_sim/` for current outputs.

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

## 🧬 Model Specification

- [ ] **Incorporate F'Y term in survival model** `[Priority: Medium]` `[Effort: Large]`
  Replace $L_i$ in the Cox linear predictor with $(F^\top y_i)$ — the projection of patient $i$'s
  observed expression onto the factor space. This decouples factor learning from survival supervision:
  L is learned purely from genomics, and survival scoring is done post-hoc via the factor projection.
  Requires a new q(L) CAVI derivation (the survival term drops out of the L update) and a post-hoc
  scoring step. May improve identifiability and reduce the sign-flip issue observed in synthetic
  validation.
  *Notes: High effort — requires new derivation and update to `code/update_L.R`. See `derivations/qL/qL_update_derivation.pdf` for current L update derivation.*

- [ ] **Add λ scaling parameter to balance genomics vs. survival objectives** `[Priority: High]` `[Effort: Medium]`
  Add a scalar λ to weight the Cox term relative to genomics: ELBO = E[L_gen] + λ·E[L_Cox] − KL.
  The genomics likelihood sums over p features while Cox sums over n patients; since p >> n in all
  cohorts, the genomics gradient dominates. λ addresses this gradient scale asymmetry without
  normalising the raw likelihoods. Grid-search or CV for λ. Expected to stabilise β estimation and
  reduce the factor sign-flip identifiability issue.
  *Notes: Medium effort — λ enters the L and β update equations as a scalar multiplier on the survival precision terms. See `DECISIONS.md` entry 2026-02-12 for the gradient asymmetry discussion. Implement before R5 (K selection), as λ affects factor stability.*

---

## 🌍 Validation & Generalization

- [x] **Cross-cohort validation: train on CPTAC+TCGA, evaluate on remaining 5 cohorts** `[Priority: High]` `[Effort: Medium]` *(In progress — see `results/benchmark_sim/`)*
  Merge CPTAC and TCGA_PAAD (batch-corrected) as a training set and hold out Dijk,
  Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, and Puleo_array as independent test cohorts.
  This tests generalisation of learned gene expression programs across platforms (RNA-seq,
  microarray, proteomics) and patient populations. More training data than the current 80/20
  within-cohort split. Requires R4 (inter-dataset normalisation) as a prerequisite.
  *Notes: `code/train_test_split.R` and `code/predict.R` are already in place for the projection step. Batch correction (ComBat or similar) not yet implemented — add to R4. See `results/modular_sim_factor/PDAC/` for current within-cohort results.*

- [ ] **Inter-dataset normalisation before cohort merging** `[Priority: High]` `[Effort: Medium]`
  Before merging cohorts (required for R3), apply quantile normalisation or z-score
  standardisation across datasets beyond the current per-cohort column-centring. Addresses
  platform-specific mean/variance shifts between RNA-seq, microarray, and proteomics. Essential
  for any analysis that combines data across assay types.
  *Notes: Current preprocessing is per-cohort column-centring only (see `load_real_data()` in `results/modular_sim_factor/run_factor_modular_simulation.R`). R4 is a prerequisite for R3.*

---

## 📐 Model Selection

- [ ] **Revisit K selection strategy after λ scaling is implemented** `[Priority: Medium]` `[Effort: Medium]`
  The current auto-prune (PVE > 1% or |β| > 0.05) uses a pragmatic screening rule. Three cohorts
  saturate at K_max = 10, suggesting more factors may be present. Revisiting K selection after R2
  (λ) is implemented will give a more stable signal for K, since the survival term will be properly
  scaled. The `select_K_cv()` stub in `code/select_K.R` is the target implementation path.
  *Notes: Depends on R2. The ELBO-based CV stub (`select_K_cv()`) exists at `code/select_K.R` lines 125–127. Three cohorts saturate K_max=10: investigate K_max=15 on Longleaf HPC.*

---

## ⚙️ Infrastructure

- [ ] **Investigate {targets} workflow for pipeline reproducibility** `[Priority: Low]` `[Effort: Medium]`
  Adopt the R `targets` package to manage computational dependencies across the pipeline (data
  loading → fitting → tables → figures → reports). Tracks which outputs are stale and re-runs only
  what is needed. Particularly valuable as the number of cohorts, conditions, and reports grows —
  currently all cohorts are re-run together even when only one has changed. Also facilitates
  Longleaf HPC submission of individual targets as SLURM jobs.
  *Notes: Low priority — investigate feasibility and overhead before committing to a refactor. The `{targets}` package works best when the pipeline DAG is well-defined; this is approximately true for the current runner structure. See `results/modular_sim_factor/run_factor_modular_simulation.R` for the current monolithic runner.*

- [ ] **Repository reorganisation** `[Priority: Low]` `[Effort: Small]`
  The current directory layout has accumulated structural debt across three simulation generations.
  Proposed clean structure (no file deletion — move and rename only):

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

- [x] **Core modular CAVI implementation** — `code/fit_modular.R` with four independently-tested update modules. 139/139 tests passing. *(Completed March 2026)*
- [x] **Hold-out prediction pipeline** — `code/predict.R` (`predict_supervised_mf()`), `code/train_test_split.R` (`stratified_split()`). 80/20 stratified hold-out. *(Completed March 2026)*
- [x] **Prior family comparison (PN vs PL × K=5 vs K_eff)** — Four-condition benchmark across 7 PDAC cohorts. Point-laplace preferred at fixed K; K selection dominates prior choice. *(Completed April 2026)*
- [x] **Full ELBO tracking** — Both proxy and full ELBO (genomics + survival + KL) tracked per iteration. `code/compute_elbo.R`. *(Completed April 2026)*
- [x] **DeSurv benchmark (Phases 0–3B)** — DeSurv-aligned preprocessing (`code/preprocess_desurv.R`), alpha CV (`code/select_alpha_cv.R`), synthetic DGP fix (equal factor amplitudes + 4-signal β), SVD pseudoinverse in `predict.R`. Synthetic: supervised 0.79 > PCA 0.76. TCGA+CPTAC real-data run in progress. *(Completed April 2026)*
- [x] **C-index honest reporting fix** — `get_cindex_comparison()` now uses model's own `EL %*% EBeta` instead of refitting coxph. Corrected `concordance()` convention mismatch: passing `I(-lp)` aligns with Cox direction (higher LP = higher risk). Synthetic C-index: Supervised 0.843 vs. PCA 0.842 with no direction workaround. Factor recovery diagnostic: factors 1, 3, 5 recover well; factor 2 correct sign but shrunk; factor 4 not recovered (b=-0.5 near-zero estimate). *(Completed April 9, 2026)*
