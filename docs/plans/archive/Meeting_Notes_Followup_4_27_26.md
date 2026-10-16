# Meeting Notes & Action Plan: Bayesian Matrix Factorization — TCGA/CPTAC Integration
**Date:** 2026-04-28  
**Attendees:** Andrew, Rashid, Yusha  
**Project:** multiomicsGEP — Semi-supervised Bayesian NMF with Survival

---

## Problem Summary

Training the semi-supervised Bayesian NMF+Survival model on **merged TCGA + CPTAC** cohorts produces factors penalized to zero, failing to capture prognostic survival signals. Training on each cohort **separately** does yield prognostic factors, strongly implicating **batch effects** as the root cause rather than a fundamental model limitation.

Key context:
- mRNA–protein correlation across genes is only ~0.4–0.6 (exists but not strong)
- DeSurv worked acceptably with rank transforms alone, but this approach performs worse in the multi-cohort setting
- Current preprocessing order is likely the primary culprit (see fix below)
- Current gene intersection yielded only **838 common genes** — Rashid flagged this as unexpectedly low and likely a processing error; expected ~2,000+

---

## Diagnostic Visualizations (Attached Images)

Two heatmaps were shared showing GEP (factor) loading matrices colored by study/cohort:

- **Image 1** (`combined_RNA_seq_gbcd_K25.png`): RNA-seq data, K=25 factors, cohorts: pancurx, mow, yeh, hayashi. Many factors show clear study-specific loading patterns (especially GEP1, GEP5, GEP7, GEP8) — consistent with batch factors dominating the factorization. GEP18–GEP20 show unusual intensity patterns worth investigating.
- **Image 2** (`combined_microarray_snn_binary_K20.png`): Microarray data, K=20, cohorts: moffitt, pukio. GEP1 appears almost entirely study-specific (solid red block for one cohort) — a textbook batch effect factor. GEP2–GEP4 also show study-stratified loading.

**Interpretation:** These heatmaps confirm that multiple factors are capturing study-of-origin (batch) rather than biology. The goal of the diagnostic plots below is to formally verify this and identify which factors have non-zero survival betas (β) vs. which are pure batch factors (β ≈ 0).

---

## Action Items

### 1. Fix Preprocessing Pipeline (PRIORITY 1)

Reorder preprocessing steps as follows — **do not filter or normalize within individual cohorts separately**:

```
Step 1: Intersect common genes across TCGA and CPTAC
         → Check: how many common genes? Expect ~2,000; current 838 suggests a bug
Step 2: Log2(x + 1) transform all values
Step 3: Quantile normalization across all merged samples
         → Matches percentile distributions across cohorts without explicit batch labels
         → Do NOT compute variance or filter before this step
Step 4: Compute variance per gene across ALL merged samples (no batch-aware normalization)
Step 5: Select top 2,000 (or 5,000 if computationally feasible) most variable genes
Step 6: Rank transform (puts all genes on common scale for model input)
Step 7: Fit model
```

**Bug check:** Investigate why common gene intersection currently yields only 838 genes.  
Likely causes: gene ID format mismatch (ENSEMBL vs. symbol), version suffixes (e.g., ENSG00000XXXXX.1), or premature per-cohort filtering removing genes before intersection.

---

### 2. Diagnostic Plot: Factor–Cohort Relationship (L matrix inspection)

Produce plots to visually identify which factors capture batch vs. biology:

#### 2a. Heatmap (already partially done — see attached images)
- Rows = samples, columns = GEPs/factors
- Color = factor loading intensity
- Annotate rows by cohort/study label

#### 2b. Box Plots (new — recommended by Rashid)
For each factor (GEP), create a box plot:
- **X-axis:** Cohort/study label (e.g., TCGA vs. CPTAC)
- **Y-axis:** Factor score (loading value) for that GEP
- **Expectation:** A batch factor will show highly separated distributions between cohorts; a biological factor will show overlapping distributions

```r
# Sketch of intended R code structure
library(ggplot2)
library(tidyr)
library(dplyr)

L_df <- as.data.frame(L_matrix)  # samples x factors
L_df$cohort <- sample_metadata$study  # TCGA or CPTAC label

L_long <- L_df %>%
  pivot_longer(cols = starts_with("GEP"), names_to = "factor", values_to = "score")

ggplot(L_long, aes(x = cohort, y = score, fill = cohort)) +
  geom_boxplot() +
  facet_wrap(~ factor, scales = "free_y") +
  theme_bw() +
  labs(title = "Factor Score Distribution by Cohort",
       x = "Study", y = "Factor Loading")
```

#### 2c. Cross-reference with Beta (β) coefficients
- For each factor, check whether the corresponding survival regression coefficient β is non-zero
- **Expected pattern for batch factors:** high cohort separation in box plot AND β ≈ 0
- **Expected pattern for biological/prognostic factors:** overlapping cohort distributions AND β ≠ 0

---

### 3. Unsupervised EBMF Sanity Check (flashier)

Run a purely unsupervised **EBMF (flashier)** model on the merged genomic data (no survival component) as a diagnostic baseline.

**Core diagnostic question:**  
> *When the semi-supervised NMF+Survival model is fit to merged TCGA+CPTAC data, all survival regression coefficients (β) shrink to zero — meaning no factors are identified as prognostic. Is this because the merged data is dominated by dataset-of-origin differences (batch effects), or because the model implementation itself has a bug?*

**Why EBMF answers this:**  
EBMF is fully unsupervised — it has no survival term and no β to shrink. It will find whatever dominant structure exists in the expression matrix, whether that is biology or batch. By then associating those factors with survival via a downstream Cox model, we can distinguish between two failure modes:

| EBMF result | Interpretation | Next action |
|---|---|---|
| EBMF finds factors associated with survival | Merged data contains survival signal; NMF+Survival model has a bug or the β prior is too aggressive | Debug model; switch β prior to normal |
| EBMF finds **no** survival-associated factors | Merged data is too batch-dominated for any method to recover survival signal | Fix preprocessing first (Action Item 1), then recheck |
| EBMF factors split cleanly by cohort with β ≈ 0 | Confirms batch effects are consuming all variance before survival can be captured | Preprocessing fix is the priority |

**Why EBMF over PCA:** EBMF adaptively learns the appropriate amount of sparsity per factor from the data, making it more sensitive to weak biological signals than PCA, which spreads variance across all components equally. It is a stronger benchmark.

**Purpose:** If EBMF applied to the merged data identifies factors associated with survival (via downstream Cox model), and the NMF+Survival model does NOT, this points to a model implementation bug rather than a data problem.

**Steps:**
```
1. Fit EBMF (flashr) to merged, preprocessed gene expression matrix
2. Extract factor loadings (L matrix from flashr)
3. Associate each factor with survival using a univariate Cox model
4. Compare: which factors are survival-associated?
5. Compare these factors to what the NMF+Survival model recovers
```

```r
# NOTE: use flashier (NOT flashr — flashr is deprecated/unmaintained)
# Current version: flashier 1.0.58
# API changed significantly at v0.2.44 — use parameters as shown below
library(flashier)
library(ebnm)  # for prior family functions (ebnm_point_normal, ebnm_normal, etc.)

# Fit unsupervised EBMF
# - greedy_Kmax: max factors to greedily add (not Kmax as in old flashr API)
# - var_type: 0 = constant, 1 = by row, 2 = by column
# - backfit = TRUE runs backfitting after greedy initialization
flash_fit <- flash(
  data        = gene_expr_matrix,   # samples x genes matrix
  greedy_Kmax = 20,                 # upper bound on K; let data decide
  var_type    = 2,                  # per-gene residual variance (recommended for genomics)
  backfit     = TRUE,
  verbose     = 1
)

# Extract loadings (LDF decomposition: Y ≈ L %*% diag(D) %*% t(F))
ldf        <- ldf(flash_fit, type = "2")  # scale so L and F columns have unit 2-norm
L_flash    <- ldf$L   # samples x K matrix (use this as factor scores)

# Cox model for each factor
library(survival)
cox_results <- lapply(1:ncol(L_flash), function(k) {
  df <- data.frame(
    factor_k = L_flash[, k],
    time     = surv_time,
    event    = surv_event
  )
  summary(coxph(Surv(time, event) ~ factor_k, data = df))
})

# Pull concordance (C-index) and p-value per factor for easy comparison
cox_summary <- data.frame(
  factor  = paste0("GEP", seq_len(ncol(L_flash))),
  concordance = sapply(cox_results, function(x) x$concordance["C"]),
  pval        = sapply(cox_results, function(x) x$coefficients[, "Pr(>|z|)"])
)
print(cox_summary[order(cox_summary$pval), ])
```

**Interpretation guide:**
- EBMF finds survival signal → merged data is clean enough; NMF+Survival model has an implementation issue
- EBMF finds NO survival signal → data itself is too dominated by batch; preprocessing fix (Action Item 1) is needed first
- EBMF >> PCA for this comparison (preferred baseline per team discussion)

---

### 4. Prior Specification Update

**Change:** Use **normal prior** (instead of point-normal/spike-and-slab) for the survival regression coefficients **β**.

**Rationale:** The number of regression parameters (β) is small (one per factor), so sparse-inducing priors are not necessary and may be overly aggressive, shrinking meaningful survival signals to zero.

```r
# In model specification, change Beta prior from:
# beta_prior <- "point_normal"   # spike-and-slab — too sparse for small p
# To:
beta_prior <- "normal"           # soft shrinkage, appropriate for small p
```

---

### 5. Future Direction: Multi-Dataset Extension

Once the current two-cohort model is working correctly:
- Extend to incorporate **multiple modalities** (RNA-seq, microarray, proteomics)
- Goal: learn **common biological structure** across datasets rather than study-specific variation
- Survival supervision provides an anchor for **generalizable, outcome-associated factors**
- Consider explicit multi-study model structure to better separate shared vs. dataset-specific components (Rashid noted this as a longer-term next step after current issues are resolved)

---

## Summary of Root Cause Hypotheses (ranked by priority)

| Priority | Hypothesis | Diagnostic | Fix |
|----------|------------|------------|-----|
| 1 | Preprocessing order wrong (filter before normalize) | Check gene counts; rerun with new order | Reorder pipeline (Action Item 1) |
| 2 | Gene intersection bug (838 vs ~2000 expected) | Print gene ID formats; debug intersection | Fix ID matching |
| 3 | Beta prior too sparse (point-normal) | Check β posteriors — are they all near zero? | Switch to normal prior (Action Item 4) |
| 4 | Model implementation bug | EBMF sanity check (Action Item 3) | Debug NMF+Survival code |
| 5 | Data fundamentally too batch-dominated | EBMF finds no survival signal either | Consider ComBat/RUVseq pre-integration |

---

## Notes for Claude Code Session

When opening this file in a Claude Code session, the immediate tasks are:

1. **Debug gene intersection** — load TCGA and CPTAC gene ID lists, identify format mismatches, recount common genes
2. **Rewrite preprocessing pipeline** in the correct order (see Step 1–7 above)
3. **Generate factor–cohort box plots** from the current L matrix output to confirm batch factor identification
4. **Run EBMF/flashr** on merged data and Cox-associate resulting factors
5. **Update model prior** for β from point-normal to normal

All code should be committed in small incremental git commits with descriptive messages. Update `README.md` and `CLAUDE.md` after each phase.
