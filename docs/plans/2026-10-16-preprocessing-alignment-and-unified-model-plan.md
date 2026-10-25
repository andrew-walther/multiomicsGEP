# Preprocessing alignment and one model implementation: plan

> **For agentic workers:** Use superpowers:executing-plans task by task. Andrew approves each **Gate**. Task A1 (documentation only) comes before the feature-selection experiment; Tasks A2–A4 come after it.

**Goal:** Single-modality and multimodal results should come from one documented preprocessing pipeline and one model implementation, so that every result is produced the same way and the two analyses can be compared.

**Why:** The two analyses currently differ in preprocessing as well as in code. Differences in expression scaling, in the TCGA expression source and in the training cohorts may explain some of the multimodal model's convergence problems, and they make expression-only versus joint comparisons hard to interpret.

## Current differences (checked against the code on `main`, 2026-10-16)

| Step | Single-modality (recommended D4) | Multimodal |
|------|----------------------------------|------------|
| Code | `fit_cox_on_yf()` | `fit_multimodal_yfb()`, which already accepts one block (expression only) |
| Training cohorts | TCGA-PAAD + CPTAC-PDAC, n = 273 | TCGA only (matched expression and methylation), n = 144 |
| Validation | Five external expression cohorts | ICGC primary PDAC |
| TCGA expression source | `PDAC_data/original/TCGA_PAAD.rds`, log2(x + 1) applied in code | `tcga/expr_log_normalized.rds`, already log2(TMM-CPM + 1) |
| Gene screen | DeSurv combined mean + variance rank, top 3,000 per cohort, intersected across cohorts → 2,064 | Same rank rule, top 3,000 on TCGA only |
| Expression scaling | Per-cohort z-standardization of each gene | None: log scale, optional per-feature intercept |
| Methylation | — | Beta values, KNN-imputed, top 10,000 CpGs by variance |
| Loading prior | Point-exponential (nonnegative) for L and F | Set per modality: point-exponential, point-Laplace or Normal |

**To examine:** the single-modality fit uses a nonnegative (point-exponential) prior on z-standardized expression, about half of which is negative by construction. In the frozen D4 fit every posterior mean of L and F is ≥ 0, so L Fᵀ can represent only above-average expression. The multimodal code checks that point-exponential blocks are nonnegative; the single-modality code does not. How much this affects the reported fit is assessed in Task A1.

## Global constraints

- Documentation first: no pipeline or default changes until Andrew approves the specification (Gate A1).
- The frozen benchmark (tag `benchmark-interim-writeup-2026-10-14`) stays reproducible. A unified implementation must reproduce the D4 fit's risk scores from its own inputs before any old code is retired.
- One choice per step, applied identically to every modality and cohort, with exceptions written down and justified (for example, CpG beta values need imputation and gene expression does not).

---

### Task A1: Preprocessing and feature-selection specification (documentation only)

**Files:** Create `docs/methods/preprocessing.md`. It becomes the methods source for the pinned book chapter and the writeup.

- [ ] **Step 1:** For each data type and cohort, write the exact steps in order: source file, transformation, missing-value handling, feature screening (rule, count, which patients it sees), scaling and centering, cross-cohort intersection, how validation data are aligned to the training features, and where each step is implemented (function and file). Cover both analyses and the feature-selection experiment's survival screen.
- [ ] **Step 2:** Quantify the prior question above on the frozen D4 fit: the share of negative training entries, how much of the reconstruction error falls on negative entries, and how the fit compares with a signed loading prior (point-Laplace) at the same settings on training data only. Report without changing any default.
- [ ] **Step 3:** End with a short list of proposed alignments, each with its reason and expected effect (for example, apply the same per-gene standardization to expression in both analyses; use one TCGA expression source). **Gate A1.**

### Task A2: Align the multimodal preprocessing (after the feature-selection experiment)

- [ ] **Step 1:** Implement only the alignments approved at Gate A1, in the shared preprocessing functions (`code/preprocess_desurv.R`, `code/load_multiomics_data.R`), with tests that each step matches the specification.
- [ ] **Step 2:** Refit the multimodal model under the aligned preprocessing with otherwise unchanged settings, and compare convergence and C with the current fit. This answers whether preprocessing caused some of the convergence problems.

### Task A3: One implementation for one or more modalities

- [ ] **Step 1:** Design note (one page): a single fitter taking a named list of data blocks, where single-modality is one block. Map every option of `fit_cox_on_yf()` and `fit_multimodal_yfb()` onto it, listing which are kept, merged or dropped. **Gate.**
- [ ] **Step 2:** Golden tests before any code moves: the unified fitter, configured as D4, reproduces the frozen D4 risk scores within 1e-6 on the training data, and reproduces the current multimodal fit's outputs under its settings.
- [ ] **Step 3:** Implement; the old fitters stay until the golden tests pass, then move to `code/legacy/` in a separate commit.

### Task A4: Re-run the headline results

- [ ] **Step 1:** Re-run the single-modality benchmark, the multimodal fit and the feature-selection comparison through the unified pipeline; report any change from the frozen numbers, with its cause.
- [ ] **Step 2:** Update `PROJECT_STATUS.qmd`, the pinned book chapter and `ROADMAP.md` Part 1.
