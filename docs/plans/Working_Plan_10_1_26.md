# Working Plan — Multimodal YFB, Real Data, and Open Single-Modality Items (10/1/26)

Brings together the notes from the 8/27 lab meeting, the 9/4 and 9/18 meetings with Yusha, and an inventory of the code and data as of 10/1. Use it as the starting point for open work. Check items off or update them in place as they are done, and record architectural choices in `DECISIONS.md`.

## Step 0 — Review the state of the branches first

We haven't confirmed where the work stopped, so start here before anything else.

1. Compare `main` with `codex/multimodal-yfb`:
   - `git log main..codex/multimodal-yfb` and `git log codex/multimodal-yfb..main`
   - Both branches contain the March `fit_modular` plan commit (`16683ed` on main, `f658acb` on the branch). Check whether it was cherry-picked.
2. List other local and remote branches and whether each is merged. In particular, check whether `fix/2026-09-04-review-findings` was ever merged into main.
3. Go through the uncommitted changes on `codex/multimodal-yfb`: 25 modified files and several untracked ones, including the simulation code, the 9/18 chapter and the simulation outputs.
   - Decide what to commit, what to drop, and what belongs on main.
4. Run `Rscript tests/run_tests.R` on the current working tree and record the pass count.
5. Write down the branch decision: keep developing on `codex/multimodal-yfb`, or merge the finished parts to main first.

### Step 0 results (2026-10-01)
- **Current branch:** `codex/multimodal-yfb`, pushed to `origin`. It is 11 commits ahead of `main`.
  - It contains the multimodal derivation, helpers, updates, fitter and prediction, plus this morning's checkpoint of the in-progress multimodal work (`6a5f0d7`) and the README/CLAUDE.md refresh (`79795ef`).
- **`main`:** its 2 commits that the branch lacks (the March `fit_modular` plan and the MIT license) are the same patches as `f658acb` and `a647b32` on the branch (checked with `git range-diff`). `main` has no content the branch is missing.
- **`fix/2026-09-04-review-findings`:** merged into `main` (commits `0808015` … `268a4f5`). The local branch no longer exists.
- **`codex/multimodal-yfb-docs`:** fully contained in `codex/multimodal-yfb`, so it can be deleted.
- **Old branches:**
  - Already merged into `main`: 16 other local branches, which can be deleted whenever convenient.
  - Not merged: `phase1-platform-norm` (2 commits) and `phase2-init-constraints` (1 commit). Both are from April–May; the second ended with the verdict "Discard". Leave them as they are unless they're needed.
- **Housekeeping (WS0) already done:** `.gitignore` covers `/data/` and `Rplots.pdf`.
- **Tests:** `Rscript tests/run_tests.R`: 509 passed, 0 failed (about 76 s).
- **Decision:** keep developing on `codex/multimodal-yfb`. Merge it to `main` after the 10/2 meeting, once the multimodal real-data prototype is in place.

## Context
- **Multimodal model:**
  - The implementation is written and tested in simulation on `codex/multimodal-yfb`, with much of the work still uncommitted.
  - It has two known problems: K_eff doesn't prune, and F supports only a point-exponential prior.
  - It has not been applied to real data.
- **Single-modality model:** several items from 8/27, 9/4 and 9/18 are still open.
- **Goal:** first multimodal real-data results by the 10/2 advisor meeting, then a steady list of work items leading up to the December prelim.

## Key facts found during exploration (reference)
- **Multimodal code:**
  - Files: `code/fit_multimodal_yfb.R`, `code/multimodal_yfb_updates.R`, `code/multimodal_yfb_helpers.R`, `code/simulate_multimodal_yfb.R`, `code/run_multimodal_yfb_simulation.R`.
  - Priors: point-exponential is the only prior for L and F, implemented as a hand-written closed-form update (`updates.R:76, 145-153`). β has a Normal empirical-Bayes prior.
  - Inputs: the preprocessing **requires nonnegative inputs**. Methylation β-values and log2 expression are nonnegative, so they work as-is. arcsine-transformed or centered data need a point-Laplace or Normal prior for F first.
  - K_eff = number of factors with PVE ≥ 0.01 (`fit_multimodal_yfb.R:462-491`).
  - Canonicalization: each factor is rescaled so sd(Σ_m Y_m F_mk) = 1 (`:162-213`). This step is not yet in the derivation.
- **Simulation results:**
  - The smoke results (`results/benchmark_sim/outputs/multimodal_yfb_smoke/`) are **stale**. All joint fits were unconverged, and the run predates the sweep cap and canonicalization.
  - The K panel was run on one draw only.
- **Data** (`data/`, untracked and about 2.4 GB):
  - **The `data/multiomicsGEP_code/` folder is not gitignored**, so a `git add data` would pick it up.
  - ICGC `icgc/icgc_matched_data.rds`:
    - `expr`: 24,033 genes × 79 donors, **not logged**. DO33168 has two RNA samples averaged into one column.
    - `meth`: 485,512 CpGs × 79, 0.18% NA.
    - `info.expr`: 80 rows, one duplicate donor. Survival is `survival_months` + `censored` (`death` = event). 67 donors have usable survival (40 deaths), 50 of them primary PDAC.
    - The cohort is mixed: 12 IPMN, 4 adenosquamous, 2 acinar, 2 metastatic.
  - ICGC `icgc/icgc_meth_filtered_sex.rds`:
    - 305,352 CpGs.
    - Columns are `submitted_donor_id`, which map in order to DO IDs through `info.expr`.
    - `meth.trans` = asin(2β−1), already computed.
  - TCGA files in `tcga/`:
    - `tcga_meth_filtered_sex.rds`: 305,492 CpGs × 150 tumours.
    - `expr_log_normalized.rds`: 20,501 × 150, log2(TMM-CPM+1).
    - Same 150 samples in the same order across files.
    - **No survival in `tcga/`**. It has to be joined from `data/PDAC_data/original/TCGA_PAAD.survival_data.rds` (4 NA times, 2 ≤ 0, 75 events).
    - `expr_raw_count.rds` has no dimnames.
    - 7 of the 17 TCGA "normals" are actually `-01A` tumour barcodes. This doesn't matter as long as we don't use the normals.
  - Overlap: the ICGC CpG set is a strict subset of TCGA's (305,352 shared; align by name, not position). 17,428 genes are shared between TCGA and ICGC.
  - Existing reference fits: `tcga/tcga_flash_K14.rds` (Yusha's unsupervised flash, K = 14, with L_pm/F_pm) and the K = 14 regression/NNLS objects. **These are a ready-made unsupervised baseline to compare factors against.**
  - The single-modality PDAC cohorts are in `data/PDAC_data/original/`, which matches the `PDAC_DATA_ROOT` layout. Setting `PDAC_DATA_ROOT=data/PDAC_data` should let `tests/test_real_data_loading.R` find them.


## Workstream 0 — Housekeeping (first 30 minutes on 10/1)
1. Add `data/multiomicsGEP_code/` (or all of `data/`) and `Rplots.pdf` to `.gitignore`.
   - Verify: `git status` no longer lists `data/`.
2. Run `Rscript tests/run_tests.R`, then make a checkpoint commit of the in-progress multimodal work on `codex/multimodal-yfb`.
   - Verify: tests pass. The expected count is 509, per the 9/18 chapter.
3. Run `PDAC_DATA_ROOT=data/PDAC_data Rscript tests/test_real_data_loading.R`.
   - Verify: 88/88, confirming the local copy is complete.

## Workstream 1 — Multimodal real-data preprocessing (10/1, priority)

**Which file to use for each piece of data.** From the 10/1 inventory; check each one when writing the loader.

| Cohort | Piece | File (under `data/multiomicsGEP_code/` unless noted) | Notes |
|---|---|---|---|
| TCGA (training) | Expression | `tcga/expr_log_normalized.rds` | 20,501 × 150, already log2(TMM-CPM+1). `expr_unlogged.rds` = `TCGA_PAAD$ex` (checked) |
| TCGA | Methylation | `tcga/tcga_meth_filtered_sex.rds` → `$meth` | 305,492 × 150 β-values. Don't use `meth.normal` or the `.trans` versions |
| TCGA | Survival | `data/PDAC_data/original/TCGA_PAAD.survival_data.rds` | Join on the 16-character barcode. Alternative source: `pc_classifier/data_reformatted_updated.rds` `$TCGA_PAAD$annotation` |
| ICGC (validation) | Expression | `icgc/icgc_matched_data.rds` → `$expr` | 24,033 × 79, **not logged**; apply log2(x+1). DO33168 is already averaged over its two samples |
| ICGC | Methylation | `icgc/icgc_matched_data.rds` → `$meth`, restricted to the TCGA filtered CpGs | Same β-values as `icgc_meth_filtered_sex.rds`, but with DO-ID columns that line up with `$expr` |
| ICGC | Survival + clinical | `icgc/icgc_matched_data.rds` → `$info.expr` | `survival_months`, `censored` (`death` = event). De-duplicate by donor. Cross-check against `data/PDAC_data/original/PACA_AU_seq.survival_data.rds`, which agrees for 69 of 80 times |
| Reference | Unsupervised baseline | `tcga/tcga_flash_K14.rds` | Yusha's flash fit (K = 14) for comparing factors |

Create `code/load_multiomics_data.R`, with roxygen-documented functions, to build the matched objects in the format the existing helpers already expect (`multimodal_yfb_helpers.R` preprocessing contract).

1. **TCGA (training cohort)**
   - Expression: start from `expr_log_normalized.rds`.
   - Methylation: `meth` from `tcga_meth_filtered_sex.rds`, using the β-values on the 0–1 scale.
   - Survival: join from `TCGA_PAAD.survival_data.rds` by 16-character barcode. Drop samples with NA or ≤ 0 time; use the `whitelist` flag (decision point).
2. **ICGC (validation cohort)**
   - Expression: `log2(x+1)` on `expr`.
   - Methylation: `meth` restricted to the TCGA filtered CpG set, matched by name.
   - Survival: de-duplicate `info.expr` by donor, then take `survival_months` and the event `censored == "death"`. Drop NA or 0 times.
   - Histology: optionally restrict to primary PDAC (decision point).
3. **Missing values:** KNN imputation (`impute::impute.knn`, the code Yusha provided) for methylation rows with any NA, done separately in each cohort. `impute` is a new dependency; it's Bioconductor and needs confirming.
4. **Feature screening**
   - Genes: the existing DeSurv gene selection (`code/preprocess_desurv.R`) or top-variance genes on the 17,428 shared genes.
   - CpGs: start with top-variance CpGs (e.g. 5k–20k) for a first fit. All 305k shared CpGs × 150 patients ≈ 370 MB per dense matrix copy, so the full set is a later Longleaf run.
   - Both screens are decided on TCGA only, so there is no leakage into the validation cohort.
5. **Standardization:** decide whether to z-score per feature. The single-modality recommended configuration uses per-platform z-standardization. The multimodal code currently requires nonnegative input, which rules out centering until WS2.1 lands.
6. Verify:
   - Unit tests in `tests/test_load_multiomics_data.R`: dims, matching IDs across modalities, no NAs after imputation, survival time > 0, ICGC CpGs ⊂ TCGA CpGs in the same order, expression is logged (max < ~25).
   - Print a summary table of n, p and number of events for each cohort.

## Workstream 2 — Multimodal model refinements
1. **Prior options for F** (requested 9/18):
   - Add point-Laplace and Normal options next to the existing point-exponential, set per modality through an argument like `prior_F = list(expr = "point_laplace", meth = "point_exponential")`. Keep point-exponential for L.
   - Implementation choice (decision point): call `ebnm::ebnm()` for the new families, or write closed-form updates to match the existing hand-written point-exponential one. `ebnm` is already a dependency and returns the posterior second moments we need. Recommendation: use `ebnm` for the new families.
   - Allow signed input (centered data, the arcsine methylation transform) when F is not point-exponential.
   - Verify: update tests where each prior recovers a known signed or sparse F, and the derivation-algebra checks in `verify_derivation.R` are extended.
2. **Intercept / centering** (9/18 note: "treat 0 as the baseline"):
   - Center features, with or without scaling, once F can be signed.
   - Compare against the uncentered nonnegative fit. Centering should also take away the mean-level structure that inflates K_eff.
3. **K pruning — the key fix for K_eff ≈ K_init.** Try these in order:
   - (a) A flashier-style nullcheck: after convergence, drop factor k if the joint ELBO, including the survival term, doesn't decrease when it's removed. Refit and repeat.
   - (b) Compute PVE on centered data, and set the threshold relative to the largest factor.
   - (c) Report external K selection (1-SE on held-out C) as a fallback.
   - Verify: replicated K panel, K_init = 5…15 × ≥ 10 seeds. The target is K_eff ≈ true K (3) across K_init.
4. **Re-run the smoke study converged**, all 7 scenarios × 10 replicates, replacing the stale CSV.
5. **Modality balance:** on the 300k-CpG scale, test whether methylation dominates. Compare the unweighted model with the a_m = ρ_m/p_m weights from the technical notes. Simulate with p_meth = 10 × p_expr.
6. **Derivation upkeep:**
   - Add canonicalization to `multimodal_YFB_derivation.qmd`.
   - Replace the "K_init reduced to K_eff by EB shrinkage" claim (`:60-62`) with whatever pruning approach we actually use.

## Workstream 3 — Real-data application (10/1 afternoon → 10/2)
1. **First fit for 10/2:**
   - Current code with point-exponential priors on nonnegative inputs (log2 expression, β-value methylation), trained on TCGA with a few K_init values (e.g. 5, 7, 10).
   - Report convergence, K_eff, survival-active factors, and C-index: 5-fold CV within TCGA plus external C on ICGC through `predict` for the multimodal fit.
2. **Baselines on the same data:**
   - (i) Single-modality YFB on expression only (`fit_cox_on_yf.R`).
   - (ii) Two-step model: unsupervised flash/EBMF factors, then Cox.
   - (iii) Yusha's existing `tcga_flash_K14` factors, then Cox.
   - Delta-C table.
3. **Factor comparison:**
   - Correlate the multimodal expression F columns with the single-modality YFB programs (Program 7 basal-like/adverse, Program 3 classical/protective) and with `tcga_flash_K14$F_pm`.
   - Check that methylation loadings run in the opposite direction to expression for matched genes, as expected biologically.
4. **Second pass, after WS2:** compare F priors (point-exponential vs point-Laplace vs Normal) and methylation scale (β vs arcsine) by C-index and interpretability.
5. **Longleaf:** the full 305k-CpG fit. Note memory needs: a few dense copies of 305k × 150, plus per-feature τ and second moments.

## Workstream 4 — Single-modality YFB items still open
1. **Clarify what ZF means in the YFB derivation (9/4, 9/18).** The code is not wrong; the derivation just never says clearly what this quantity is.
   - The derivation writes η = (Y·E[F])β̃ with raw ZF. The code uses unit-L2-normalized columns, `ZF = Y %*% EF_norm` (`fit_cox_on_yf.R:643-647`), and prediction applies the stored norms.
   - Add a definition to `derivations/cox_on_YF/` that states Ẑ = Y F̄ D⁻¹ (the projection score, with columns normalized), and settle the subscript notation Yusha asked about.
   - Write the multimodal raw-vs-normalized difference into both documents.
2. **α / α_F, the Bayesian concern.** This is a large problem and will take time. **Before 10/2: write proposed solutions only.** Start the experiments once the multimodal model has a working prototype on real data.
   - Candidate approaches to write up, each with a mechanism and a test:
     - The multimodal model's factor rescaling.
     - A plain sum with a warm start from the α_F = 0 fit.
     - Turning survival on in the F update gradually over iterations.
     - Changing the order of the updates.
     - Treating α as a likelihood power with a stated justification.
   - Try importing the multimodal canonicalization into single-modality YFB, then test α = 1 (no tempering) and α_F = 1 (plain sum).
   - If that is stable and C is within noise, drop both weights. This answers the 9/18 question "use sum of genomics + survival — does this break it?"
   - Note that YFB has never had a dedicated CV over α; the "validated by CV" comment in `globals.yml` comes from the LB model and should be corrected.
   - **Show that supervision improves F (author's direction, 10/1; added from the bios-dissertation abstract session).**
     - Why: with `alpha_F = 0`, the Cox term never reaches the F (or L) updates and survival informs only β (`DECISIONS.md` ~336–352). The prelim literature review and abstract want to claim that supervision changes which programs are found. The abstract's Ch. 4 sentences are scoped down until this is shown (`bios-dissertation/prelim/abstract/notes/abstract-fact-sheet.md`; `bios-dissertation/ROADMAP.md` Project 3 to-do).
     - Goal: remove `alpha_F` entirely, so the YFB Cox term enters the F update as the derivation implies. At minimum, fit with `alpha_F > 0`.
     - Prior evidence to read first: the corrected 16-configuration (N_frozen × alpha_F) grid found no `alpha_F > 0` setting above the `alpha_F = 0` baseline (mean external C 0.6267), and some below it. The earlier `alpha_F > 0` instability did not reproduce (`DECISIONS.md` ~786–803).
     - Simulation: plant a low-variance prognostic program, then compare recovery of its true loading (and survival-active classification) under `alpha_F = 0` vs `alpha_F > 0` / no `alpha_F`. This is the test that external C alone cannot give.
     - Real data: compare against the `alpha_F = 0` baseline on the 5 held-out PDAC cohorts, using the same pooled-bootstrap ΔC.
     - Write the success criteria down before running. Decide in advance how a null result would be reported.
3. **β prior (8/27 g, 9/4 i):** run Normal vs point-Laplace vs point-normal for β on the current YFB config. The April comparison used the older LB model.
4. **Multiple initializations (9/4 j):**
   - Random starts currently collapse to k_eff ≈ 0, which is why it is SVD-init only.
   - Try jittered SVD starts (SVD + small noise) × 20–100, averaging C for each K. That answers the advisor's initialization-sensitivity question without random starts that are known to collapse.
5. **False-positive rate with `beta_cohort_id` (8/27 c):** one simulation run comparing FP rates with and without it.
6. **PCA/t-SNE/UMAP of L (8/27 h):** colour patients by dominant program and by survival group or subtype, for the prelim.
7. **Formal 1-SE rule in `run_k_init_sweep.R` (8/27 e),** applied in the same way to the multimodal K panel.
8. **The K = 10–15 dip (8/27 d):** low priority. Try warm-starting from the K = 7 solution, or deflation-style initialization.
9. **Doc and code inconsistencies to clean up:**
   - The `DECISIONS.md` K = 11/13 entry contradicts the 9/4 chapter.
   - The stale sign-correction comment in `compute_bic.R:189-209`.
   - `fit_modular.R:834` calls `concordance()` without `reverse = TRUE`.
   - The `auto_prune_K` default `beta_thresh = 0.05` is stale.
10. **CPTAC vs TCGA survival (9/4):** compare the KM curves and the baseline hazard between the two training cohorts. This shows whether a shared baseline hazard and shared β are reasonable, and supports the cohort-specific-coefficient interpretation.

## Workstream 5 — 10/2 progress-book chapter (the end point of this plan)
Pull the work above into a new chapter, `docs/progress_book/chapters/2026-10-02.qmd`. The previous chapter is `2026-09-18.qmd`. Register the new chapter in `docs/progress_book/_quarto.yml`, render with `quarto render`, and check that the GitHub Pages publish workflow picks it up.

This chapter is the end-to-end progress record for both advisors. Planned sections:
1. **Summary.** A short list of what changed since 9/18 and the main results.
2. **Single-modality YFB: supervision in the F update.** This is a short section.
   - Restate the problem: α_F = 0 keeps survival out of the F update, so F is the same as in unsupervised EBMF. α = 0.5 tempers the Cox likelihood for β.
   - Present the proposed solutions (WS4.2) and the planned test design, including the low-variance prognostic program simulation.
   - Add results only if any experiments were run.
   - Include the ZF clarification (WS4.1).
3. **Multimodal model form.** Naim has not seen it yet, so this section should stand on its own:
   - Likelihood for each modality: Y_m = L F_mᵀ + E_m, with per-feature precision τ_mj.
   - Priors: L is shared with a point-exponential prior. Each F_m has its own prior family (point-exponential, point-Laplace or Normal). β has a Normal empirical-Bayes prior.
   - Survival: η = Σ_m Y_m F_m β, a Cox partial likelihood.
   - The factor rescaling (identifiability) step, and how K_eff is determined.
   - Summarize from `derivations/multimodal_YFB/multimodal_YFB_derivation.qmd` and link to the full derivation rather than repeating it.
4. **Multimodal refinements and simulation (WS2).**
   - Prior options for F.
   - Centering.
   - The K-pruning fix, with the replicated K_init panel and K_eff compared to the true K.
   - The re-run of the converged smoke study.
   - Modality balance at a realistic p_meth/p_expr ratio.
5. **Real data (WS1 + WS3).**
   - Data description: TCGA training set, ICGC validation set; n, events, genes and CpGs before and after screening.
   - Preprocessing and transformations: log2, CpG matching, imputation, screening, standardization.
   - Fit results: K_eff and survival-active factors, internal CV C and external ICGC C against the expression-only YFB and two-step baselines (delta-C table), and comparison of factors with the single-modality programs and with `tcga_flash_K14`.
6. **Open questions and next steps,** including the decision points below that are still open.

Also: update `DECISIONS.md` for each decision made, `ROADMAP.md` when milestones are reached, and `PROJECT_STATUS.qmd`. Prelim deadlines are the abstract on 11/2 and the proposal on 11/16. The chapter's method and results sections should be reusable as drafts for the prelim.

## Decision points (to settle together before or while executing)
1. Should ICGC validation be restricted to primary PDAC (n ≈ 50 with survival) or use all usable donors (n ≈ 67)?
2. Methylation scale for the first fit: β-values, which the current nonnegative code supports, or arcsine (asin(2β−1)), which needs WS2.1 first.
3. Feature screening for genes and CpGs: method and size.
4. Whether to add `impute` (Bioconductor) as a dependency.
5. How to implement the F priors: `ebnm` or hand-written closed-form updates.

## Review of the meeting notes (8/27, 9/4, 9/18)
Every item from the notes the user shared on 10/1 is listed here, so that none are dropped. Status as of 10/1:
- **Done:** completed and documented.
- **Planned:** has a workstream in this plan.
- **Considered:** evaluated and deliberately not planned now, with the reason.
Re-check this list before writing the 10/2 chapter, and again before the prelim.

### 8/27 lab meeting
| Item | Status | Where / note |
|---|---|---|
| Under-specify K_init (does the fit merge true factors?) | Done | 9/4 chapter §4; under-specified fits merge factors, over-specified fits do not |
| Are both methods fit to the same simulated dataset? | Done | Yes, both arms use the same draws (`results/multi_cohort_sim/make_delta_fig.R`) |
| Bar chart of ΔC between methods | Done | Same script |
| False-positive rate of survival-active factors; does a cohort-specific effect change it? | Planned | WS4.5. Measured without `beta_cohort_id` only |
| Let the model choose K; use C-index as the model-free metric | Single-modality: done. Multimodal: planned | Two-stage K selection exists for single-modality YFB. WS2.3 covers multimodal pruning |
| C drops and recovers at K = 10–13; what differs at K = 13–15? | Planned (low priority) | WS4.8. Ruled out as a bad starting point; cause still unexplained |
| Amber's approach: max over the full grid, then most parsimonious within 1 SE | Planned | WS4.7 (single-modality sweep) and WS2.3(c) (multimodal) |
| Is the log-likelihood cross-validated? | Done | `code/compute_cv_loglik.R`; the 8/27 numbers were in-sample |
| BIC on the regular log-likelihood; full log-likelihood hard in high dimensions | Done, with a caveat | `compute_bic.R` uses E_q log p plus the Cox partial log-likelihood, with a log(n patients) penalty. State the penalty choice in the chapter or proposal |
| Penalize small β to zero (elastic net, as Amber did) | Considered | The EBNM prior already shrinks β; `DECISIONS.md` gives the argument. Revisit if WS4.3 priors fail |
| Point-normal or Laplace prior on β, data-adaptive | Planned | WS4.3. Point-normal collapsed β to 0 on PDAC before; Laplace not yet tested on YFB |
| Prior on K, or a constraint summed across factors | Considered | Not implemented (IBP-type prior is a literature item). The multimodal K-pruning work (WS2.3) comes first. List under future work in the proposal |
| PCA/t-SNE to explain the factors | Planned | WS4.6 |
| Naim: results and method are a good basis for a paper | n/a | Supports the proposal plan |

### 9/4 meeting with Yusha
| Item | Status | Where / note |
|---|---|---|
| Ridge vs Laplace prior for β (EBNM) | Planned | WS4.3. The current prior is Normal with an empirical-Bayes variance |
| ~100 random initializations, average C when selecting K | Planned | WS4.4, using jittered SVD starts because pure random starts collapse. Apply to the multimodal K panel as well |
| The unsupervised joint model misses low-variance prognostic factors | Planned | WS4.2 simulation design. The multimodal simulator already has a `low_variance` scenario; report it in WS2.4 |
| \|β\| threshold relative to the largest β | Done | Implemented as `rel_thresh`, kept for simulation only |
| Cohort indicator: fixed L column, β set to 0 | Done | `fit_cox_on_yf.R`; β_cohort = 0 by construction |
| Revisit the CPTAC vs TCGA survival comparison | **Planned (new)** | WS4.10 below |
| Cohort-specific coefficients: differences may reflect small samples | Done (interpretation) | State this explicitly when cohort-specific results are presented |
| ZF notation (subscript?) | Planned | WS4.1 |
| Multimodal: methylation needs a different, sparser prior; matched training and validation cohorts; all factors shared; methylation and expression anti-correlated | Planned | WS1, WS2.1, WS3. The anti-correlation check is WS3.3 |
| Yusha to share matched TCGA/ICGC data | Done | In `data/` |
| Cohorts vs modalities: extra columns vs extra rows | Planned | Methods framing in the chapter (WS5 §3) and the proposal |
| Derive the key multimodal equations | Done | `derivations/multimodal_YFB/` |

### 9/18 meeting with Yusha
| Item | Status | Where / note |
|---|---|---|
| Let the user choose the F prior (point-exponential, point-Laplace, Normal) and give scenario-specific recommendations | Planned | WS2.1. Add a simulation comparison by scenario so a recommendation can be made |
| Keep point-exponential for L | Done | Current design |
| Intercept, or center Y so 0 is the baseline | Planned | WS2.2 |
| Methylation on 0–1; check normality; optional asin(2x−1) | Planned | WS1, decision point 2. Add a normality check (histograms and QQ plots before and after the transform) |
| Expression log-normalized (Gaussian errors) | Planned | WS1 |
| ICGC: load the matched object, link by donor ID, average tumours from the same patient, `survival_months`/`censored`, log2(x+1) | Planned | WS1. The matched object already averages DO33168 |
| Cross-check ICGC against PACA_AU_seq | **Planned (new check)** | WS1: survival times agree for 69 of 80 samples. Look at the 11 that don't before choosing a survival source |
| Cross-check TCGA expression against Amber's PDAC_data tables | Done in the inventory | `expr_unlogged` = `TCGA_PAAD$ex` exactly. Record this in the WS1 loader test |
| TCGA: use `meth`, not the normal or transformed versions | Planned | WS1 |
| Filter ICGC CpGs to the TCGA set | Planned | WS1 (ICGC ⊂ TCGA) |
| KNN imputation of missing methylation values | Planned | WS1, decision point 4 |
| Fit TCGA, report C-index, compare gene factors with the earlier method | Planned | WS3 |
| α_F: try the plain sum; how to order updates to avoid collapse | Planned (after the prototype) | WS4.2, optional step 6 |
| Simulation goal: identify the correct number of survival-active factors | Planned | WS2.3/2.4: make survival-active count recovery a primary simulation metric |
| Double-check the model section | Planned | WS5 §3. Reconcile the derivation with the code (rescaling step, K_eff claim, WS2.6) |
| Consolidate the multi-omics data into one OneDrive folder and share it | **Planned (admin)** | Data are local in `data/`. Upload and share once preprocessing is settled |
| Prelim content for this project | Planned | `docs/plans/Prelim_Proposal_Plan_10_1_26.md` |
| Lab-meeting update before the prelim | Done | 10/1 |
| Ask Dr. Lin for the sudden-unexpected-death datasets | Outside this project | Project 2. Status unknown; the bios-dissertation notes say the Project 2 application is waiting on that dataset |

## Order of work (end to end)
The lab-meeting update has been given (10/1). What remains is doing the work. **The priority is a working multimodal prototype on real data.**

1. **Step 0 + WS0:** review the branches, housekeeping, and a test baseline.
2. **WS2, multimodal refinements.** The derivation has advisor approval.
   - Order: WS2.1 (F priors) → WS2.2 (centering) → WS2.3 (K pruning, with the replicated panel) → WS2.4 (converged smoke re-run).
   - WS2.5–2.6 if time allows.
3. **WS1 → WS3, real data.** Preprocess matched TCGA/ICGC, then fit the multimodal model and the baselines.
   - If WS2 is not finished, run WS3.1 on the current point-exponential code with nonnegative inputs, so that some real-data results exist for 10/2.
4. **WS4, short single-modality items.**
   - WS4.1: the ZF clarification.
   - WS4.2: the α/α_F proposed solutions, written up but not run.
   - The α/α_F experiments and WS4.3–4.9 come after the multimodal prototype exists.
5. **WS5:** write, render and publish the 10/2 progress-book chapter.
6. **Optional, only after the first version of the chapter exists and the multimodal work is in good shape:** run the α/α_F experiments (WS4.2) and add the results to the chapter as a follow-up section. The multimodal implementation comes first.
   - Deadline: the α_F question must be answered, or the proposal's framing scoped down, before the prelim proposal is due (11/16). See `docs/plans/Prelim_Proposal_Plan_10_1_26.md`.

If time runs short before 10/2, keep steps 1, 3 (the minimal version) and 5, and report the unfinished parts of step 2 as in progress in the chapter.

## Handoff between Claude Code and Codex
Work may move between Claude Code and Codex when usage limits are reached. Both read this file; `AGENTS.md` is a symlink to `CLAUDE.md`.
- After each step, add a dated line to the **Progress log** below: what was done, which files and commits, and what comes next.
- Commit at the end of each step so the other tool starts from a clean tree.
- Record any decision point that has been settled in `DECISIONS.md` and in the log.

## Progress log
- 2026-10-01: Plan written. Lab-meeting update given.
- 2026-10-01: Step 0 done (results above). WS0 was already done in `6a5f0d7`/`79795ef`. Tests 509/509. Plan files committed. Next: WS2.1 (F prior options) and WS1 (data loader), working on both in parallel.

## Verification (overall)
- `Rscript tests/run_tests.R` passes after any change to the model code, with the new tests added to the count.
- The real-data loader tests pass locally.
- Every reported C-index is computed on held-out data. TCGA uses CV folds; ICGC is a fully external cohort that plays no part in feature screening or preprocessing decisions.
