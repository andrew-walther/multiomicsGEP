# Feature-selection comparison: experiment plan

> **For agentic workers:** Use superpowers:executing-plans (or subagent-driven-development) task by task, with test-driven steps. Andrew approves each **Gate** before the next task. Start only after the round2 close-out plan's Tasks 1–3 (`docs/plans/2026-10-16-round2-closeout-and-docs-plan.md`), because Task E4 uses `code/feature_sensitivity.R` ported there, and after Task A1 of `docs/plans/2026-10-16-preprocessing-alignment-and-unified-model-plan.md`, which documents the pipeline this experiment uses. The experiment runs on the current single-modality pipeline as it is; if Task A2/A3 later changes preprocessing, its results are re-run in Task A4.

**Goal:** Test whether a survival-informed gene screen lets the model reach the same prediction and program structure with fewer genes than the current unsupervised mean/variance screen.

**Question (from the 10/16 advisor meeting):** compare (1) unsupervised selection by mean and variance with (2) the same selection followed by a screen for genes associated with survival. If (2) performs as well with a smaller gene set, the model can be fitted on fewer features with the same prediction and program information.

**Architecture:** One branch `exp/feature-selection` off `main`. Three small tested helpers in `code/`, one runner in `results/feature_selection/`, compact committed summaries, fits git-ignored.

**Tech stack:** R; `survival`, `ebnm`, `yaml`; existing `fit_cox_on_yf()`, `predict_cox_on_yf()`, `preprocess_merged_cohorts()`, `create_stratified_folds()`, `classify_factors()`.

**Spec:** this document's Question and Decisions sections.

## Model and data (fixed for every arm)

- Training: TCGA-PAAD + CPTAC-PDAC pooled, n = 273 (`cfg$pdac$training_cohorts`). External: the five cohorts in `cfg$pdac$external_cohorts`.
- Preprocessing as for the recommended configuration (D4): `preprocess_merged_cohorts(train_raw, log_transform_flags = cfg$pdac$platform_log_transform[...], rank_transform = FALSE, per_platform_standardize = TRUE, normalize_method = "none", selection_per_cohort = TRUE, selection_method = "combined_rank", top_n = <per-cohort count>)`. With `top_n = 3000` this gives the 2,064-gene baseline.
- Fit: `set.seed(42L); fit_cox_on_yf(Y, time, status, K = cfg$benchmark$k_merged_yfb_desurv (7), max_iter = cfg$cavi$max_iter (300), alpha = cfg$benchmark$alpha (0.5), prior_beta = "normal", verbose = FALSE)`. α_F stays at its default 0.
- Folds: `create_stratified_folds(status_train, n_folds = 5, seed = 42)`, the same folds for every arm.
- Scoring: `predict_cox_on_yf()` returns η = (YF)β. Its sign is not guaranteed: the frozen D4 fit's η is a good-prognosis score, and the writeup uses risk = −η. So each fit's direction d = ±1 comes from `training_direction()` (Task E1) on the patients that fit was trained on, and held-out risk is d·η, scored with `harrell_c()`. Never `oriented_cindex()`, which picks the direction on the evaluation data.

## Arms

| Arm | Genes | Selection uses survival? |
|-----|-------|--------------------------|
| **U2064** | The D4 baseline, 2,064 genes | No (reference) |
| **U(N)** | Unsupervised mean/variance screen with the per-cohort count lowered until about N genes are shared | No: matched-size control |
| **S(N)** | The 2,064 baseline genes ranked by a cohort-stratified univariate Cox test on **fold-training patients only**; keep the top N | Yes |

The U(N) control is what makes the comparison informative: S(N) must beat an equally small unsupervised set, not only the larger baseline.

## Decisions to confirm at Gate E2

- N grid. Proposed: N ∈ {250, 500, 1000}, using U(N)'s actual shared counts as the S(N) sizes so the pairs match exactly.
- Screen statistic. Proposed: Wald p-value from `coxph(Surv(time, status) ~ y + strata(cohort))`, ranked ascending. Stratifying by cohort keeps TCGA-versus-CPTAC differences from passing as survival association.
- Optional third arm: screen all 17,900 shared genes by survival alone. Proposed: leave out for now.

## Global constraints

- **No leakage:** the survival screen runs inside each fold on that fold's training rows; held-out rows never inform it. The full-training fit for S(N) screens on all 273 patients.
- External cohorts never choose an arm, an N or a setting. They are scored once (Task 5), for all arms, after the CV results are reviewed.
- Per-cohort z-standardization runs once on all training patients before the folds, as in the existing benchmark; it is unsupervised and identical across arms. Report CV C as conditional on that step.
- Fitted objects and per-patient scores are git-ignored (under `results/feature_selection/outputs/fits/`); only the summary CSVs and figures are committed.
- Report non-converged fits, constant-risk folds (C undefined, recorded as NA) and errors in the tables; never drop them.

## Review focus

1. A survival screen given N larger than the number of candidate genes must stop with an error, not return fewer genes or all genes. (`cox_feature_selection()`'s fallback of returning every gene is exactly this failure; do not reuse it.)
2. A gene that is constant within a fold, or whose Cox fit does not converge, must get p = 1 (ranked last) and be counted in the manifest; it must not abort the run or silently disappear.
3. U(N) may not hit N exactly. The manifest must record the actual shared count, and S(N) must use that same count.
4. External cohorts missing some selected genes: score on the shared genes (as the benchmark does) and record the number used. Stop if fewer than 100 are shared.
5. The risk direction must come from training patients. A test should confirm that permuting the held-out patients' outcomes leaves d unchanged.
6. A fold whose fit gives constant risk yields C = NA, kept in the table and excluded from that arm's mean with the count of excluded folds shown.

---

### Task E1: Evaluation and screening helpers

**Files:**
- Modify: `code/concordance_ci.R` (add `harrell_c`, `training_direction`), `code/feature_selection.R` (add `survival_screen_genes`), `code/preprocess_desurv.R` (add `unsupervised_gene_set`)
- Test: `tests/test_feature_selection_comparison.R`, registered in `tests/run_tests.R`

**Interfaces (produced):**
- `harrell_c(risk, time, status)` → numeric scalar C, or `NA_real_` when `sd(risk) == 0`. Higher risk means shorter survival; no direction flipping.
- `training_direction(eta, time, status)` → `1L` if `harrell_c(eta, time, status) >= 0.5`, else `-1L`; `NA_integer_` if C is NA. Called only on training patients.
- `survival_screen_genes(Y, time, status, strata, n_keep)` → integer vector of length `n_keep`: column indices of `Y` in ascending order of the stratified Cox Wald p-value. Attribute `"p_values"` holds all p-values (length `ncol(Y)`); attribute `"n_failed"` counts genes assigned p = 1.
- `unsupervised_gene_set(cohort_raw_list, log_transform_flags, n_shared)` → the `preprocess_merged_cohorts()` result for the smallest per-cohort `top_n` whose shared gene count is ≥ `n_shared`, plus elements `top_n_per_cohort` and `n_shared_actual`.

- [ ] **Step 1:** Write the failing tests:
  - `harrell_c` equals `survival::concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance` on a fixture; returns NA for constant risk; returns a value below 0.5 for an inverted risk (no flipping).
  - `training_direction` returns `-1L` for a good-prognosis score and `1L` for a risk score on the same fixture.
  - `survival_screen_genes` on simulated data where columns 1–5 drive survival returns 1–5 among its top 5; `length(result) == n_keep`; `stop()` when `n_keep > ncol(Y)`; a constant column gets p = 1 and `attr(, "n_failed") == 1`; a column associated with survival only through stratum membership is not ranked in the top 5 when `strata` is given.
  - `unsupervised_gene_set` with `n_shared` = the baseline count (on a small synthetic two-cohort fixture) returns nested gene sets: the set for a smaller `n_shared` is a subset of the set for a larger one.
- [ ] **Step 2:** Run `Rscript tests/run_tests.R`; the new tests fail because the functions don't exist.
- [ ] **Step 3:** Implement the three functions with the signatures above. `unsupervised_gene_set` searches `top_n` by bisection between 1 and the largest cohort's gene count.
- [ ] **Step 4:** Run the full suite; 0 failures.
- [ ] **Step 5:** Commit.

### Task E2: Preflight manifest (no fitting)

**Files:** Create `results/feature_selection/run_feature_selection_comparison.R` with a `manifest` mode; output `results/feature_selection/outputs/manifest.csv`.

- [ ] **Step 1:** In `manifest` mode, build U2064, each U(N), and for each fold plus full training the S(N) gene sets. Write one row per (arm, N, fold): requested N, actual genes, overlap with U2064, overlap between S(N) and U(N), `n_failed`. Load no external cohort.
- [ ] **Step 2:** Run it: `Rscript results/feature_selection/run_feature_selection_comparison.R manifest`. Check that U2064 has 2,064 genes, which reproduces the baseline.
- [ ] **Step 3:** Show Andrew the manifest and the decisions listed above. **Gate E2.**

### Task E3: Cross-validated comparison

**Files:** the same runner, `cv` mode; outputs `results/feature_selection/outputs/cv_results.csv` (committed) and fits under `outputs/fits/` (ignored; add to `.gitignore`).

- [ ] **Step 1:** For each arm × N × fold: build the gene set (S(N) from that fold's training rows), fit on fold-training rows, set d with `training_direction()` on the fold-training η, score fold-test rows as d·η with `harrell_c()`. Record per row: arm, N, n_genes, fold, d, C, converged, n_iter, K_eff = `sum(abs(fit$EBeta) > cfg$k_selection$beta_threshold)`, runtime in seconds.
- [ ] **Step 2:** For each arm × N, also fit on all 273 training patients and save the fit (ignored). Record `classify_factors()` counts per class.
- [ ] **Step 3:** Run with `cv`. Expected scale: 7 arm-size combinations × (5 folds + 1 full) = 42 fits; estimate the runtime from the first fit and report it before continuing if the total exceeds 2 hours locally.
- [ ] **Step 4:** Summarize: mean CV C ± SE per arm × N, paired differences S(N) − U(N) on the same folds, K_eff and convergence counts.

### Task E4: Program agreement

**Files:** the same runner, `agreement` mode; output `results/feature_selection/outputs/agreement.csv`.

**Interfaces (consumed):** `feature_program_agreement(old_fit, old_genes, new_fit, new_genes, old_Y, new_Y)` from `code/feature_sensitivity.R` (ported in close-out Task 3).

- [ ] **Step 1:** For each full-training fit, compare its programs with the U2064 fit: matched-program loading correlation on shared genes, activity correlation on the 273 patients, risk-score correlation. Report constant or unmatched programs as NA, not 0.
- [ ] **Step 2:** Show Andrew the CV and agreement tables. **Gate E4:** confirm the arms are final before any external scoring.

### Task E5: External scoring (once) and findings

**Files:** the same runner, `external` mode; outputs `external_results.csv`, `figures/c_vs_genes.pdf/.png`; a findings entry in `ROADMAP.md` Part 1.

- [ ] **Step 1:** Score every full-training fit on the five external cohorts as d·η, with d from that fit's training patients: per-cohort C with `harrell_c()`, the mean, and the number of genes used per cohort. Use the external preprocessing that the benchmark uses (`preprocess_desurv_cohort(..., top_n = NULL, rank_transform = FALSE, per_platform_standardize = TRUE)`).
- [ ] **Step 2:** Plot CV C and external mean C against the number of genes, U and S arms as separate lines, with U2064 marked.
- [ ] **Step 3:** Write a one-paragraph findings entry in `ROADMAP.md` Part 1 and recommend merge (helpers and the result) or discard (tag `archive/exp-feature-selection`, delete the branch). **Gate E5:** Andrew decides.
