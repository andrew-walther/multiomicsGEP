# Preprocessing Comparison & K-Floor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the merged-cohort benchmark to cover three additional preprocessing methods (joint z-standardization, quantile-without-rank, log-only), apply a biological K floor (K ≥ 3) to all K-CV selections, and produce a comprehensive 18-configuration comparison table with a clear recommended configuration.

**Architecture:** Add a `normalize_method` parameter to `preprocess_merged_cohorts()` to support joint z-standardization and no-normalization modes. Rewrite the K-CV and benchmark runners as loops over a preprocessing configuration table instead of hard-coded sequential calls. Apply a K floor of `max(K_1se, 3)` in the K-CV runner after `select_K_cv()` returns. All work continues on the `merged-benchmark` branch.

**Tech Stack:** R, `survival`, `ebnm`, `yaml`, `ggplot2`, `preprocessCore`, Quarto. Branch: `merged-benchmark`.

---

## Background and Design Decisions

### K-floor rationale
The 1-SE rule selects K=2 for YFB per-platform because the merged training folds (n≈55/fold) have high variance, making the 1-SE band very wide. K=3 is within the same band and performs better in CV (mean C=0.634 vs 0.625 for K=2). The biological floor `K_final = max(K_1se, 3)` is documented explicitly in the runner and in comments. Motivation: DeSurv uses K=3–4; the manuscript claims "gene expression programs" (plural, implying ≥3); K=3 has lower SE than K=2 on the YFB per-platform CV table.

The floor is implemented in the runner script (not in `select_K_cv()`) so the underlying function remains reusable. Both K_1se and K_final are logged.

### Preprocessing configs to benchmark
Five total; two already have results from the prior session:

| Label | per_platform_standardize | normalize_method | rank_transform | Status |
|-------|--------------------------|------------------|----------------|--------|
| `joint_quantile_rank` | FALSE | `"quantile"` | TRUE | Done (M1, M2) |
| `perplatform_zstd` | TRUE | `"quantile"` | FALSE | Done (M3–M6); YFB K re-run at K=3 |
| `joint_quantile_norank` | FALSE | `"quantile"` | FALSE | **New** (M7, M8, M13, M14) |
| `joint_zstd` | FALSE | `"z_score"` | FALSE | **New** (M9, M10, M15, M16) |
| `log_only` | FALSE | `"none"` | FALSE | **New** (M11, M12, M17, M18) |

`joint_quantile_norank` is already expressible with existing parameters. `joint_zstd` and `log_only` require the new `normalize_method` parameter (Task 1).

### What re-runs are needed from the prior session
- M5 and M6 (YFB × perplatform_zstd × ±cohort) must be re-run at K=3 (the floor overrides the prior K=2 selection). All other existing M1–M4 K values (6, 6, 3, 3) already satisfy the K≥3 floor and do not change.
- K-CV for `perplatform_zstd × YFB` must be re-run with the K floor logged (K_1se=2 → K_final=3); the resulting globals.yml value changes from 2 to 3.

### Model ID allocation
M1–M6: existing (unchanged except M5/M6 re-run at K=3).
M7–M18: new configurations:

| ID | Model | Preprocessing | cohort_id |
|----|-------|---------------|-----------|
| M7 | LB | joint_quantile_norank | No |
| M8 | LB | joint_quantile_norank | Yes |
| M9 | LB | joint_zstd | No |
| M10 | LB | joint_zstd | Yes |
| M11 | LB | log_only | No |
| M12 | LB | log_only | Yes |
| M13 | YFB | joint_quantile_norank | No |
| M14 | YFB | joint_quantile_norank | Yes |
| M15 | YFB | joint_zstd | No |
| M16 | YFB | joint_zstd | Yes |
| M17 | YFB | log_only | No |
| M18 | YFB | log_only | Yes |

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `code/preprocess_desurv.R` | **Modify** | Add `normalize_method` parameter to `preprocess_merged_cohorts()` |
| `tests/test_preprocess_desurv.R` | **Modify** | Add 4 tests for new `normalize_method` modes |
| `tests/run_tests.R` | No change | Already sources `test_preprocess_desurv.R` |
| `config/globals.yml` | **Modify** | Add 6 new K placeholder keys; update `k_merged_yfb_perplatform` from 2 to 3 |
| `results/benchmark_sim/run_merged_kcv.R` | **Rewrite** | Loop over all 5 preprocessing configs; apply K floor; use `set_key()` (not `replace_null_key()`) |
| `results/benchmark_sim/run_merged_benchmark.R` | **Rewrite** | Loop over all 5 preprocessing configs; re-runs M5/M6 at new K |
| `docs/reports/merged_benchmark_report.qmd` | **Modify** | Extend table/heatmap to 18 models; update conclusions |

---

## Task 1: Add `normalize_method` parameter to `preprocess_merged_cohorts()`

Modifies step 5 (quantile normalization) of the existing pipeline to be conditional on a `normalize_method` argument. `"quantile"` preserves current behavior. `"z_score"` replaces the quantile step with a column-wise center+scale. `"none"` skips normalization entirely. Per-platform z-standardization (step 4) and rank transform (step 7) remain independently controlled by their existing flags.

**Files:**
- Modify: `code/preprocess_desurv.R`
- Modify: `tests/test_preprocess_desurv.R`

- [ ] **Step 1: Write 4 failing tests for `normalize_method` in `tests/test_preprocess_desurv.R`**

Append the following at the end of `tests/test_preprocess_desurv.R` (after the existing tests):

```r
# ---------------------------------------------------------------------------
# Tests for normalize_method parameter in preprocess_merged_cohorts()
# ---------------------------------------------------------------------------

# Minimal synthetic merged input: 2 cohorts, 4 subjects, 6 genes
make_synthetic_raw_list <- function() {
  set.seed(1L)
  gene_names <- paste0("G", 1:6)
  list(
    CohortA = list(
      Y          = matrix(round(2^(matrix(runif(12, 2, 8), 4, 6)), 1),
                          nrow = 4, dimnames = list(NULL, gene_names)),
      gene_names = gene_names,
      n          = 4L,
      time       = c(1.0, 2.0, 3.0, 4.0),
      status     = c(1L, 0L, 1L, 0L)
    ),
    CohortB = list(
      Y          = matrix(round(2^(matrix(runif(12, 1, 9), 4, 6)), 1) * 0.5,
                          nrow = 4, dimnames = list(NULL, gene_names)),
      gene_names = gene_names,
      n          = 4L,
      time       = c(0.5, 1.5, 2.5, 3.5),
      status     = c(0L, 1L, 0L, 1L)
    )
  )
}

run_test("T1.8: normalize_method='quantile' (default) matches prior behavior", {
  raw <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out_q   <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                       normalize_method = "quantile",
                                       rank_transform = FALSE)
  out_def <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                       rank_transform = FALSE)
  # quantile is the default — outputs must be identical
  assert_near(out_q$Y, out_def$Y, tol = 1e-12)
})

run_test("T1.9: normalize_method='z_score' gives column means ~0 and SDs ~1", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out   <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                     normalize_method = "z_score",
                                     rank_transform   = FALSE)
  col_means <- colMeans(out$Y)
  col_sds   <- apply(out$Y, 2, sd)
  # Each gene column should be centered and scaled across all 8 subjects
  assert_near(col_means, rep(0, ncol(out$Y)), tol = 1e-10)
  assert_near(col_sds,   rep(1, ncol(out$Y)), tol = 1e-10)
})

run_test("T1.10: normalize_method='none' preserves log-transformed values (no QN distortion)", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out_none <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                        normalize_method = "none",
                                        rank_transform   = FALSE)
  out_q    <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                        normalize_method = "quantile",
                                        rank_transform   = FALSE)
  # 'none' and 'quantile' must produce DIFFERENT matrices
  # (if they were the same, normalization had no effect — that would be surprising)
  assert_equal(identical(out_none$Y, out_q$Y), FALSE)
  # 'none' matrix gene names and dimensions must still be correct
  assert_equal(ncol(out_none$Y), 6L)
  assert_equal(nrow(out_none$Y), 8L)
  assert_equal(out_none$gene_names, out_q$gene_names)
})

run_test("T1.11: normalize_method invalid argument is caught", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  result <- tryCatch(
    preprocess_merged_cohorts(raw, flags, top_n = 6L,
                              normalize_method = "bad_method"),
    error = function(e) e
  )
  assert_equal(inherits(result, "error"), TRUE)
})
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
Rscript -e "
  source('tests/test_helpers.R')
  source('code/preprocess_desurv.R')
  source('tests/test_preprocess_desurv.R')
" 2>&1 | grep -E "PASS|FAIL|T1\.[89]|T1\.1[01]"
```

Expected: T1.8–T1.11 all FAIL with errors about `unused argument (normalize_method = ...)` or similar.

- [ ] **Step 3: Add `normalize_method` parameter to `preprocess_merged_cohorts()` in `code/preprocess_desurv.R`**

Find the function signature at line 233. Change:

```r
preprocess_merged_cohorts <- function(cohort_raw_list,
                                      log_transform_flags,
                                      top_n                    = 2000,
                                      ties_method              = "average",
                                      rank_transform           = TRUE,
                                      per_platform_standardize = FALSE) {
```

To:

```r
preprocess_merged_cohorts <- function(cohort_raw_list,
                                      log_transform_flags,
                                      top_n                    = 2000,
                                      ties_method              = "average",
                                      rank_transform           = TRUE,
                                      per_platform_standardize = FALSE,
                                      normalize_method         = c("quantile", "z_score", "none")) {
  normalize_method <- match.arg(normalize_method)
```

Then find the quantile normalization block (currently unconditional):

```r
  cat(sprintf("  [v2] Quantile normalising merged matrix (%d x %d) ...\n",
              nrow(Y_merged), ncol(Y_merged)))
  Y_qn <- quantile_normalize_merged(Y_merged)
```

Replace it with:

```r
  if (normalize_method == "quantile") {
    cat(sprintf("  [v2] Quantile normalising merged matrix (%d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- quantile_normalize_merged(Y_merged)
  } else if (normalize_method == "z_score") {
    cat(sprintf("  [v2] Joint z-standardizing merged matrix (%d x %d, colMean=0, colSD=1) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- scale(Y_merged, center = TRUE, scale = TRUE)
    colnames(Y_norm) <- common_genes
  } else {
    # normalize_method == "none": skip normalization; pass log-transformed matrix directly.
    # Gene-level mean differences across platforms are NOT removed here.
    # Downstream column-centering inside fit_modular/fit_cox_on_yf provides
    # partial correction, but platform-scale differences remain.
    cat(sprintf("  [v2] Skipping normalization (log-transform only; %d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- Y_merged
  }
```

Then find all downstream uses of `Y_qn` and replace with `Y_norm`. There are two:

```r
  # Change this:
  selected <- select_top_variable_genes(Y_qn, common_genes, top_n = top_n)
  cat(sprintf("  [v2] Genes retained after top-%d variance filter: %d\n",
              top_n, length(selected$gene_names)))

  # To this:
  selected <- select_top_variable_genes(Y_norm, common_genes, top_n = top_n)
  cat(sprintf("  [v2] Genes retained after top-%d variance filter: %d\n",
              top_n, length(selected$gene_names)))
```

The rank transform block references `selected$Y` directly (no `Y_qn`), so it requires no change.

Also add `normalize_method` to the returned list (after `per_platform_standardize`):

```r
  list(
    Y                        = Y_final,
    gene_names               = selected$gene_names,
    n                        = nrow(Y_final),
    p                        = ncol(Y_final),
    dataset_labels           = dataset_labels,
    n_raw_intersect          = length(common_genes),
    rank_transform           = rank_transform,
    per_platform_standardize = per_platform_standardize,
    normalize_method         = normalize_method
  )
```

- [ ] **Step 4: Run tests — expect T1.8–T1.11 to PASS, all prior tests still PASS**

```bash
Rscript -e "
  source('tests/test_helpers.R')
  source('code/preprocess_desurv.R')
  source('tests/test_preprocess_desurv.R')
" 2>&1 | grep -E "PASS|FAIL|Final"
```

Expected: all tests PASS.

- [ ] **Step 5: Run full test suite to confirm no regressions**

```bash
Rscript tests/run_tests.R 2>&1 | tail -5
```

Expected: `238 passed, 0 failed`.

- [ ] **Step 6: Commit**

```bash
git add code/preprocess_desurv.R tests/test_preprocess_desurv.R
git commit -m "Add normalize_method parameter to preprocess_merged_cohorts() (z_score, none modes)

Adds joint z-standardization and no-normalization preprocessing options
needed for the extended preprocessing comparison benchmark. 'quantile'
(default) preserves existing behavior. Four new tests: T1.8-T1.11."
```

---

## Task 2: Update `config/globals.yml` with new K keys and K-floor fix

Adds six new placeholder K keys for the three new preprocessing options (two models each). Also updates `k_merged_yfb_perplatform` from 2 to 3 (applying the K floor: 1-SE selected K=2 but biological minimum is K=3).

**Files:**
- Modify: `config/globals.yml`

- [ ] **Step 1: Add 6 new K placeholder keys and update YFB perplatform K**

Find the existing K keys block in `config/globals.yml`:

```yaml
  k_merged_lb_joint:        6   # LB, joint quantile+rank preprocessing (M1, M2)
  k_merged_lb_perplatform:  3   # LB, per-platform z-std preprocessing (M3, M4)
  k_merged_yfb_perplatform: 2   # YFB, per-platform z-std preprocessing (M5, M6)
```

Replace with:

```yaml
  k_merged_lb_joint:        6   # LB, joint quantile+rank preprocessing (M1, M2)
  k_merged_lb_perplatform:  3   # LB, per-platform z-std preprocessing (M3, M4)
  k_merged_yfb_perplatform: 3   # YFB, per-platform z-std preprocessing (M5, M6)
                                # NOTE: 1-SE selected K=2; biological floor (K>=3) applied.
                                # K=3 has higher CV mean C (0.634 vs 0.625) and lower SE.
  # New preprocessing options — filled by run_merged_kcv.R (biological floor K>=3 applied)
  k_merged_lb_joint_norank:     null  # LB,  joint quantile, no rank transform (M7, M8)
  k_merged_yfb_joint_norank:    null  # YFB, joint quantile, no rank transform (M13, M14)
  k_merged_lb_zstd:             null  # LB,  joint z-standardization (M9, M10)
  k_merged_yfb_zstd:            null  # YFB, joint z-standardization (M15, M16)
  k_merged_lb_logonly:          null  # LB,  log transform only (M11, M12)
  k_merged_yfb_logonly:         null  # YFB, log transform only (M17, M18)
```

- [ ] **Step 2: Verify globals.yml loads correctly**

```bash
Rscript -e "
  cfg <- yaml::read_yaml('config/globals.yml')
  b   <- cfg\$benchmark
  cat('k_merged_yfb_perplatform:', b\$k_merged_yfb_perplatform, '\n')
  cat('k_merged_lb_zstd:', deparse(b\$k_merged_lb_zstd), '\n')
  cat('k_merged_yfb_logonly:', deparse(b\$k_merged_yfb_logonly), '\n')
"
```

Expected:
```
k_merged_yfb_perplatform: 3
k_merged_lb_zstd: NULL
k_merged_yfb_logonly: NULL
```

- [ ] **Step 3: Commit**

```bash
git add config/globals.yml
git commit -m "Update globals.yml: apply K floor (YFB perplatform K 2->3), add new preprocessing K keys

K=3 biological floor: 1-SE selected K=2 for YFB per-platform but K=3 performs
better in CV (mean C=0.634 vs 0.625) and has lower fold SE (0.022 vs 0.030).
Floor K>=3 motivated by DeSurv K=3-4 and multi-program biological interpretation.
Adds 6 null K keys for new preprocessing options (M7-M18)."
```

---

## Task 3: Rewrite `run_merged_kcv.R` to loop over all preprocessing options

The existing script runs 3 K-CV calls in sequence with hard-coded preprocessing parameters. Rewrite it as a loop over a preprocessing configuration table. Key changes:
1. Config table covers all 5 preprocessing options (2 already have K values; 3 are new).
2. `set_key()` replaces `replace_null_key()` — it overwrites any existing value (not just null), enabling clean re-runs.
3. K floor applied after `select_K_cv()` returns: `K_final = max(K_opt, K_MIN_BIOLOGICAL)`.

**Files:**
- Rewrite: `results/benchmark_sim/run_merged_kcv.R`

- [ ] **Step 1: Replace `results/benchmark_sim/run_merged_kcv.R` with the loop-based version**

```r
# ============================================================
# Script:  results/benchmark_sim/run_merged_kcv.R
# Purpose: Select K via cross-validated C-index (1-SE rule, biological floor)
#          on merged TCGA_PAAD + CPTAC training data for all preprocessing x
#          model configurations used in the merged benchmark.
#
#          Biological K floor: K_final = max(K_1se, K_MIN_BIOLOGICAL) where
#          K_MIN_BIOLOGICAL = 3. Motivation: DeSurv (Young et al. 2025) uses
#          K=3-4; the model identifies biological programs (plural requires K>=3);
#          for YFB per-platform, 1-SE selects K=2 but K=3 has higher CV mean C
#          (0.634 vs 0.625) and lower fold SE (0.022 vs 0.030).
#
#          Preprocessing configurations:
#            joint_quantile_rank    - joint QN + rank transform (M1, M2)
#            perplatform_zstd       - per-platform z-std + joint QN (M3-M6)
#            joint_quantile_norank  - joint QN, no rank transform (M7, M8, M13, M14)
#            joint_zstd             - joint z-standardization (M9, M10, M15, M16)
#            log_only               - log transform only (M11, M12, M17, M18)
#
#          Results written to config/globals.yml under benchmark.k_merged_*.
#          Previously computed values are overwritten (not just nulls).
#
#          NOTE on leakage: preprocessing applied to full training matrix before
#          CV fold split. This is standard for K selection, not final inference.
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_kcv.R [--quick]
#          --quick: K_grid=2:4, n_folds=3, max_iter=30 (~3 min)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
source("code/select_K.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

K_GRID            <- if (QUICK_MODE) 2L:4L  else 2L:10L
N_FOLDS           <- if (QUICK_MODE) 3L      else cfg$cavi$n_cv_folds    # 5
MAX_ITER          <- if (QUICK_MODE) 30L     else cfg$cavi$max_iter      # 300
ALPHA             <- cfg$benchmark$alpha                                  # 0.5
LAMBDA            <- cfg$benchmark$lambda                                 # 1.0
PRIOR_BETA        <- "normal"
K_MIN_BIOLOGICAL  <- 3L   # biological floor: K_final = max(K_1se, K_MIN_BIOLOGICAL)

cat("=== Merged-Cohort K-CV (all preprocessing options) ===\n")
cat(sprintf("    K_grid=%d:%d | n_folds=%d | max_iter=%d | K_floor=%d | QUICK=%s\n\n",
            min(K_GRID), max(K_GRID), N_FOLDS, MAX_ITER, K_MIN_BIOLOGICAL, QUICK_MODE))

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load raw training data (done once; preprocessing varies below)
# --------------------------------------------------------------------------

cat("--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga        <- train_raw$TCGA_PAAD$n
n_cptac       <- train_raw$CPTAC$n
time_train    <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train  <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Preprocessing configuration table
#    Each row: label, globals_key_lb, globals_key_yfb, preprocess_args
# --------------------------------------------------------------------------

PREPROC_CONFIGS <- list(
  list(
    label       = "joint_quantile_rank",
    key_lb      = "k_merged_lb_joint",
    key_yfb     = NULL,  # YFB x joint QN excluded (beta->0 structural, see DECISIONS.md 2026-05-22)
    per_plat    = FALSE,
    norm_method = "quantile",
    rank        = TRUE
  ),
  list(
    label       = "perplatform_zstd",
    key_lb      = "k_merged_lb_perplatform",
    key_yfb     = "k_merged_yfb_perplatform",
    per_plat    = TRUE,
    norm_method = "quantile",
    rank        = FALSE
  ),
  list(
    label       = "joint_quantile_norank",
    key_lb      = "k_merged_lb_joint_norank",
    key_yfb     = "k_merged_yfb_joint_norank",
    per_plat    = FALSE,
    norm_method = "quantile",
    rank        = FALSE
  ),
  list(
    label       = "joint_zstd",
    key_lb      = "k_merged_lb_zstd",
    key_yfb     = "k_merged_yfb_zstd",
    per_plat    = FALSE,
    norm_method = "z_score",
    rank        = FALSE
  ),
  list(
    label       = "log_only",
    key_lb      = "k_merged_lb_logonly",
    key_yfb     = "k_merged_yfb_logonly",
    per_plat    = FALSE,
    norm_method = "none",
    rank        = FALSE
  )
)

# --------------------------------------------------------------------------
# Helper: write any value (null or integer) back to globals.yml key
# --------------------------------------------------------------------------

set_key <- function(lines, key, value) {
  pattern <- paste0("^\\s*", key, ":")
  idx     <- grep(pattern, lines)
  if (length(idx) == 0) stop(sprintf("Key '%s' not found in globals.yml", key))
  # Replace the value token (first non-whitespace after the colon) with the new value,
  # preserving any trailing comment.
  lines[idx] <- sub(
    paste0("(^\\s*", gsub("_", "_", key), ":\\s*)\\S+"),
    paste0("\\1", as.character(value)),
    lines[idx]
  )
  lines
}

# --------------------------------------------------------------------------
# Helper: run K-CV for one (preprocessing, model) combination
# --------------------------------------------------------------------------

run_kcv <- function(Y, time, status, model, label) {
  cat(sprintf("--- K-CV: %s x %s ---\n", label, model))
  set.seed(42L)
  result <- select_K_cv(
    Y, time, status,
    K_grid     = K_GRID,
    n_folds    = N_FOLDS,
    model      = model,
    seed       = 42L,
    verbose    = TRUE,
    max_iter   = MAX_ITER,
    prior_beta = PRIOR_BETA,
    alpha      = ALPHA,
    lambda     = LAMBDA,
    sign_correction = FALSE
  )
  K_1se   <- result$K_opt
  K_final <- max(K_1se, K_MIN_BIOLOGICAL)
  if (K_final > K_1se) {
    cat(sprintf("\n  NOTE: 1-SE selected K=%d; biological floor (K>=%d) applied -> K_final=%d\n",
                K_1se, K_MIN_BIOLOGICAL, K_final))
  } else {
    cat(sprintf("\n  -> K_opt=%d (1-SE rule; no floor applied)\n", K_1se))
  }
  cat("\n  CV table:\n")
  print(result$cv_table[, c("K", "mean_cindex", "se_cindex")], row.names = FALSE)
  out_csv <- file.path(OUT_DIR, sprintf("kcv_%s_%s.csv",
                                        gsub(" ", "_", tolower(label)),
                                        tolower(model)))
  write.csv(result$cv_table, out_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", out_csv))
  list(K_opt = K_final, K_1se = K_1se, cv_table = result$cv_table)
}

# --------------------------------------------------------------------------
# 3. Run K-CV for all configurations, write results back to globals.yml
# --------------------------------------------------------------------------

globals_text <- readLines("config/globals.yml")
k_results    <- list()

for (pcfg in PREPROC_CONFIGS) {

  cat(sprintf("\n========== Preprocessing: %s ==========\n\n", pcfg$label))

  Y_prep <- preprocess_merged_cohorts(
    cohort_raw_list     = train_raw,
    log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
    top_n               = cfg$preprocessing$top_n_genes,
    rank_transform      = pcfg$rank,
    per_platform_standardize = pcfg$per_plat,
    normalize_method    = pcfg$norm_method
  )$Y
  cat(sprintf("  Y: n=%d, p=%d\n\n", nrow(Y_prep), ncol(Y_prep)))

  if (!is.null(pcfg$key_lb)) {
    r_lb <- run_kcv(Y_prep, time_train, status_train, "LB", pcfg$label)
    k_results[[pcfg$key_lb]] <- r_lb$K_opt
    globals_text <- set_key(globals_text, pcfg$key_lb, r_lb$K_opt)
  }
  if (!is.null(pcfg$key_yfb)) {
    r_yfb <- run_kcv(Y_prep, time_train, status_train, "YFB", pcfg$label)
    k_results[[pcfg$key_yfb]] <- r_yfb$K_opt
    globals_text <- set_key(globals_text, pcfg$key_yfb, r_yfb$K_opt)
  }
}

# --------------------------------------------------------------------------
# 4. Write all K values back to globals.yml in one atomic write
# --------------------------------------------------------------------------

writeLines(globals_text, "config/globals.yml")
cat("\n--- K values written to config/globals.yml ---\n")
for (nm in names(k_results)) cat(sprintf("  %-35s = %d\n", nm, k_results[[nm]]))

cat("\n============================================================\n")
cat(" K-CV complete. Run run_merged_benchmark.R to fit all models.\n")
cat("============================================================\n")
```

- [ ] **Step 2: Verify syntax**

```bash
Rscript -e "parse(file='results/benchmark_sim/run_merged_kcv.R'); cat('OK\n')"
```

Expected: `OK`

- [ ] **Step 3: Run in quick mode to verify end-to-end**

```bash
Rscript results/benchmark_sim/run_merged_kcv.R --quick 2>&1 | tail -20
```

Expected: K values printed for each config, globals.yml updated with integers. Verify:

```bash
grep "k_merged" config/globals.yml
```

Expected: All 8 keys have integer values (not null), and `k_merged_yfb_perplatform` is 3 (floor applied) or the quick-mode K, whichever is larger.

- [ ] **Step 4: Reset NEW keys only back to null; keep existing values intact**

The existing M1–M6 K values (6, 3, 3) are already correct. Reset only the 6 new keys so the full run fills them:

```bash
sed -i '' 's/k_merged_lb_joint_norank:     [0-9]*/k_merged_lb_joint_norank:     null/' config/globals.yml
sed -i '' 's/k_merged_yfb_joint_norank:    [0-9]*/k_merged_yfb_joint_norank:    null/' config/globals.yml
sed -i '' 's/k_merged_lb_zstd:             [0-9]*/k_merged_lb_zstd:             null/' config/globals.yml
sed -i '' 's/k_merged_yfb_zstd:            [0-9]*/k_merged_yfb_zstd:            null/' config/globals.yml
sed -i '' 's/k_merged_lb_logonly:          [0-9]*/k_merged_lb_logonly:          null/' config/globals.yml
sed -i '' 's/k_merged_yfb_logonly:         [0-9]*/k_merged_yfb_logonly:         null/' config/globals.yml
grep "k_merged" config/globals.yml
```

Expected: existing 3 keys show integers (6, 3, 3); new 6 keys show `null`.

- [ ] **Step 5: Commit**

```bash
git add results/benchmark_sim/run_merged_kcv.R config/globals.yml
git commit -m "Rewrite run_merged_kcv.R: loop over all 5 preprocessing configs, apply K>=3 biological floor

Replaces hard-coded 3-call sequence with a preprocessing config table.
Key changes:
- set_key() replaces replace_null_key() — overwrites any value, enabling re-runs
- K_final = max(K_1se, 3) biological floor with explicit logging
- Covers all 5 preprocessing modes: joint_quantile_rank, perplatform_zstd,
  joint_quantile_norank, joint_zstd, log_only
- YFB x joint QN excluded (structural beta->0, documented)"
```

- [ ] **Step 6: Run full K-CV in background (~3 hours for new preprocessing options)**

Confirm machine is on AC power. The existing M1–M6 K values will be re-confirmed and the 6 new keys filled:

```bash
caffeinate -s Rscript results/benchmark_sim/run_merged_kcv.R > /tmp/merged_kcv_extended.log 2>&1 &
echo "PID: $!"
```

Monitor: `tail -f /tmp/merged_kcv_extended.log`

Check running: `ps aux | grep "run_merged_kcv" | grep -v grep`

- [ ] **Step 7: Verify and commit after K-CV completes**

```bash
grep "k_merged" config/globals.yml   # all 8 keys should show integers >= 3
git add config/globals.yml results/benchmark_sim/outputs/merged_benchmark/kcv_*.csv
git commit -m "Fill K values for all preprocessing options (biological floor K>=3 applied)"
```

---

## Task 4: Rewrite `run_merged_benchmark.R` to loop over all preprocessing options

Restructure the benchmark script to mirror the config-table approach from the K-CV runner. The same 5 preprocessing configs drive the model fitting. Model IDs M1–M6 are preserved for continuity; M7–M18 are the new configs. M5 and M6 are re-run at K=3 (floor applied).

**Files:**
- Rewrite: `results/benchmark_sim/run_merged_benchmark.R`

- [ ] **Step 1: Replace `results/benchmark_sim/run_merged_benchmark.R`**

```r
# ============================================================
# Script:  results/benchmark_sim/run_merged_benchmark.R
# Purpose: Comprehensive merged-cohort benchmark — all preprocessing options.
#          Fits 18 model configurations at CV-selected K (biological floor K>=3)
#          and evaluates external C-index on 5 held-out PDAC cohorts.
#
#          Model IDs:
#            M1–M6:   existing (joint_quantile_rank, perplatform_zstd) x LB/YFB x ±cohort
#            M7–M18:  new (joint_quantile_norank, joint_zstd, log_only) x LB/YFB x ±cohort
#
#          K values read from globals.yml. Run run_merged_kcv.R first.
#          YFB x joint QN excluded (structural beta->0, DECISIONS.md 2026-05-22).
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-25
# Usage:   Rscript results/benchmark_sim/run_merged_benchmark.R [--quick]
#          --quick: max_iter=30, skip top-gene table (smoke test)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

cfg <- yaml::read_yaml("config/globals.yml")
b   <- cfg$benchmark

# --------------------------------------------------------------------------
# Guard: all K values must be filled
# --------------------------------------------------------------------------

REQUIRED_KEYS <- c("k_merged_lb_joint", "k_merged_lb_perplatform",
                   "k_merged_yfb_perplatform", "k_merged_lb_joint_norank",
                   "k_merged_yfb_joint_norank", "k_merged_lb_zstd",
                   "k_merged_yfb_zstd", "k_merged_lb_logonly", "k_merged_yfb_logonly")
missing_keys <- REQUIRED_KEYS[sapply(REQUIRED_KEYS, function(k) is.null(b[[k]]))]
if (length(missing_keys) > 0) {
  stop(sprintf("K values not set: %s\nRun: Rscript results/benchmark_sim/run_merged_kcv.R",
               paste(missing_keys, collapse = ", ")))
}

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")

ALPHA       <- b$alpha
LAMBDA      <- b$lambda
MAX_ITER    <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA  <- "normal"
SIGMA_COH   <- 1.0
BETA_THRESH <- cfg$k_selection$beta_threshold

OUT_DIR <- "results/benchmark_sim/outputs/merged_benchmark"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# 1. Load training data
# --------------------------------------------------------------------------

TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))

# --------------------------------------------------------------------------
# 2. Model configuration table
#    Each row: model ID, model type, preprocessing params, cohort_id flag, K key
# --------------------------------------------------------------------------

MODEL_CONFIGS <- list(
  # --- existing (joint_quantile_rank) ---
  list(id="M1",  model="LB",  per_plat=FALSE, norm="quantile", rank=TRUE,  cohort=FALSE, k_key="k_merged_lb_joint"),
  list(id="M2",  model="LB",  per_plat=FALSE, norm="quantile", rank=TRUE,  cohort=TRUE,  k_key="k_merged_lb_joint"),
  # --- existing (perplatform_zstd) --- re-run M5/M6 at new K=3
  list(id="M3",  model="LB",  per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_lb_perplatform"),
  list(id="M4",  model="LB",  per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_perplatform"),
  list(id="M5",  model="YFB", per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_perplatform"),
  list(id="M6",  model="YFB", per_plat=TRUE,  norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_perplatform"),
  # --- new (joint_quantile_norank) ---
  list(id="M7",  model="LB",  per_plat=FALSE, norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_lb_joint_norank"),
  list(id="M8",  model="LB",  per_plat=FALSE, norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_joint_norank"),
  list(id="M13", model="YFB", per_plat=FALSE, norm="quantile", rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_joint_norank"),
  list(id="M14", model="YFB", per_plat=FALSE, norm="quantile", rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_joint_norank"),
  # --- new (joint_zstd) ---
  list(id="M9",  model="LB",  per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=FALSE, k_key="k_merged_lb_zstd"),
  list(id="M10", model="LB",  per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_zstd"),
  list(id="M15", model="YFB", per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_zstd"),
  list(id="M16", model="YFB", per_plat=FALSE, norm="z_score",  rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_zstd"),
  # --- new (log_only) ---
  list(id="M11", model="LB",  per_plat=FALSE, norm="none",     rank=FALSE, cohort=FALSE, k_key="k_merged_lb_logonly"),
  list(id="M12", model="LB",  per_plat=FALSE, norm="none",     rank=FALSE, cohort=TRUE,  k_key="k_merged_lb_logonly"),
  list(id="M17", model="YFB", per_plat=FALSE, norm="none",     rank=FALSE, cohort=FALSE, k_key="k_merged_yfb_logonly"),
  list(id="M18", model="YFB", per_plat=FALSE, norm="none",     rank=FALSE, cohort=TRUE,  k_key="k_merged_yfb_logonly")
)

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Preprocess training data — cache by (per_plat, norm, rank) combo
# --------------------------------------------------------------------------

cat("--- Preprocessing training data (all modes) ---\n")
preproc_cache <- list()
gene_set_cache <- list()

for (mcfg in MODEL_CONFIGS) {
  cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
  if (!cache_key %in% names(preproc_cache)) {
    cat(sprintf("  Preprocessing: per_plat=%s, norm=%s, rank=%s ...\n",
                mcfg$per_plat, mcfg$norm, mcfg$rank))
    pp <- preprocess_merged_cohorts(
      train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
      top_n = cfg$preprocessing$top_n_genes,
      rank_transform = mcfg$rank, per_platform_standardize = mcfg$per_plat,
      normalize_method = mcfg$norm
    )
    preproc_cache[[cache_key]]  <- pp$Y
    gene_set_cache[[cache_key]] <- pp$gene_names
    cat(sprintf("    n=%d, p=%d\n", nrow(pp$Y), ncol(pp$Y)))
  }
}

# --------------------------------------------------------------------------
# 4. Fit all configurations
# --------------------------------------------------------------------------

cat("\n=== Fitting all 18 configurations ===\n\n")
fits <- list()

for (mcfg in MODEL_CONFIGS) {
  cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
  Y_train   <- preproc_cache[[cache_key]]
  K         <- b[[mcfg$k_key]]
  cohort_id <- if (mcfg$cohort) cohort_labels else NULL

  label <- sprintf("%s %s per_plat=%s norm=%s rank=%s cohort=%s K=%d",
                   mcfg$id, mcfg$model, mcfg$per_plat, mcfg$norm, mcfg$rank, mcfg$cohort, K)
  cat(sprintf("--- Fitting %s ---\n", label))
  set.seed(42L)

  fit <- suppressMessages(
    if (mcfg$model == "LB")
      fit_supervised_mf_modular(Y_train, time_train, status_train,
                                K = K, max_iter = MAX_ITER, alpha = ALPHA,
                                lambda = LAMBDA, prior_beta = PRIOR_BETA,
                                verbose = TRUE, cohort_id = cohort_id,
                                sigma_F_cohort = SIGMA_COH)
    else
      fit_cox_on_yf(Y_train, time_train, status_train,
                   K = K, max_iter = MAX_ITER, alpha = ALPHA,
                   lambda = LAMBDA, prior_beta = PRIOR_BETA,
                   verbose = TRUE, cohort_id = cohort_id,
                   sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fits[[mcfg$id]] <- fit
}

# --------------------------------------------------------------------------
# 5. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(
    Y = raw_ext$Y, gene_names = raw_ext$gene_names,
    top_n = cfg$preprocessing$top_n_genes,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name = ext_cohort
  )

  for (mcfg in MODEL_CONFIGS) {
    cache_key   <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
    train_genes <- gene_set_cache[[cache_key]]
    common      <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) next

    Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx <- match(common, train_genes)
    fit       <- fits[[mcfg$id]]
    EF_sub    <- fit$EF[train_idx, , drop = FALSE]

    if (mcfg$model == "LB") {
      pred  <- predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
    } else {
      pred  <- predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)
    }
    c_val <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      model      = mcfg$id,
      model_type = mcfg$model,
      preprocess = paste(mcfg$norm, if(mcfg$per_plat) "perplat" else "joint",
                         if(mcfg$rank) "rank" else "norank", sep="_"),
      has_cohort = mcfg$cohort,
      K          = b[[mcfg$k_key]],
      cohort     = ext_cohort,
      c_index    = round(c_val, 4),
      k_eff      = sum(abs(fit$EBeta) > BETA_THRESH),
      beta_max   = round(max(abs(fit$EBeta)), 4),
      n_iters    = fit$history$n_iter,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, results_rows)

# --------------------------------------------------------------------------
# 6. Print summary
# --------------------------------------------------------------------------

cat("\n============================================================\n")
cat(" Merged Benchmark Results - Mean External C-index by Model\n")
cat("============================================================\n")
model_order <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10",
                 "M11","M12","M13","M14","M15","M16","M17","M18")
for (mid in model_order) {
  sub <- results_df[results_df$model == mid, ]
  if (nrow(sub) == 0) next
  cat(sprintf("  %s [%s, %s, cohort=%s, K=%d]: mean C=%.3f | K_eff=%d | beta_max=%.4f\n",
              mid, sub$model_type[1], sub$preprocess[1], sub$has_cohort[1],
              sub$K[1], mean(sub$c_index), sub$k_eff[1], sub$beta_max[1]))
}

# --------------------------------------------------------------------------
# 7. Top-20 gene table
# --------------------------------------------------------------------------

if (!QUICK_MODE) {
  cat("\n--- Factor top-20 genes ---\n")
  top_genes_rows <- list()
  for (mcfg in MODEL_CONFIGS) {
    cache_key <- paste(mcfg$per_plat, mcfg$norm, mcfg$rank, sep = "_")
    genes     <- gene_set_cache[[cache_key]]
    fit       <- fits[[mcfg$id]]
    for (k in seq_len(ncol(fit$EF))) {
      idx <- order(abs(fit$EF[, k]), decreasing = TRUE)[1:min(20, nrow(fit$EF))]
      top_genes_rows[[length(top_genes_rows) + 1]] <- data.frame(
        model = mcfg$id, factor = k, gene = genes[idx],
        loading = round(fit$EF[idx, k], 4), stringsAsFactors = FALSE
      )
    }
  }
  top_genes_df <- do.call(rbind, top_genes_rows)
  write.csv(top_genes_df, file.path(OUT_DIR, "merged_benchmark_top_genes_extended.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved top genes: %s\n",
              file.path(OUT_DIR, "merged_benchmark_top_genes_extended.csv")))
}

# --------------------------------------------------------------------------
# 8. Save results
# --------------------------------------------------------------------------

write.csv(results_df,
          file.path(OUT_DIR, "merged_benchmark_results_extended.csv"),
          row.names = FALSE)

compact_fits <- lapply(fits, function(f) list(
  EBeta = f$EBeta, EBeta2 = f$EBeta2, EF = f$EF,
  EF_cohort = f$EF_cohort, EF2_cohort = f$EF2_cohort,
  EF_norms = if (!is.null(f$EF_norms)) f$EF_norms else NULL,
  history = f$history
))
saveRDS(list(fits = compact_fits, results = results_df,
             gene_sets = gene_set_cache,
             params = list(ALPHA=ALPHA, PRIOR_BETA=PRIOR_BETA, K_floor=3L),
             date = Sys.time()),
        file.path(OUT_DIR, "merged_benchmark_fits_extended.rds"))

cat(sprintf("\nResults saved to: %s\n", OUT_DIR))
cat("============================================================\n")
```

- [ ] **Step 2: Verify syntax**

```bash
Rscript -e "parse(file='results/benchmark_sim/run_merged_benchmark.R'); cat('OK\n')"
```

Expected: `OK`

- [ ] **Step 3: Quick-mode smoke test (requires K values filled from Task 3 Step 3 quick run)**

```bash
grep "k_merged" config/globals.yml  # all 8 keys should have integer values
Rscript results/benchmark_sim/run_merged_benchmark.R --quick 2>&1 | tail -25
```

Expected: mean C-index printed for all 18 model IDs (or fewer if K values not yet filled for new preprocessing options). No errors.

- [ ] **Step 4: Commit**

```bash
git add results/benchmark_sim/run_merged_benchmark.R
git commit -m "Rewrite run_merged_benchmark.R: config-table loop over all 18 model configurations

Covers 5 preprocessing options x 2 models x +-cohort_id.
Caches preprocessed matrices to avoid redundant computation.
Saves extended results CSV and top-gene table with _extended suffix
to preserve prior M1-M6 outputs for reference."
```

- [ ] **Step 5: Run full benchmark in background after K-CV completes (~3–4 hours)**

```bash
caffeinate -s Rscript results/benchmark_sim/run_merged_benchmark.R > /tmp/merged_benchmark_extended.log 2>&1 &
echo "PID: $!"
```

Monitor: `tail -f /tmp/merged_benchmark_extended.log`

Check running: `ps aux | grep "run_merged_benchmark" | grep -v grep`

- [ ] **Step 6: Commit results after benchmark completes**

```bash
git add results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results_extended.csv
git add results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_top_genes_extended.csv
git add results/benchmark_sim/outputs/merged_benchmark/kcv_*.csv
git commit -m "Add extended benchmark results: 18-config comparison across all preprocessing options"
```

---

## Task 5: Update `merged_benchmark_report.qmd` with comprehensive results

Extends the report to cover all 18 configurations. The heatmap, summary table, and conclusions are updated with actual findings.

**Files:**
- Modify: `docs/reports/merged_benchmark_report.qmd`

- [ ] **Step 1: Replace the setup chunk and table sections in `merged_benchmark_report.qmd`**

Replace the existing `{r setup}` chunk:

```r
library(yaml)
library(ggplot2)
library(tidyr)
library(dplyr)

cfg     <- yaml::read_yaml(here::here("config/globals.yml"))
results <- read.csv(here::here(
  "results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results_extended.csv"))
top_genes <- tryCatch(
  read.csv(here::here(
    "results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_top_genes_extended.csv")),
  error = function(e) NULL)

model_labels <- c(
  M1="LB|joint_qn_rank|no_coh",     M2="LB|joint_qn_rank|cohort",
  M3="LB|perplat_zstd|no_coh",      M4="LB|perplat_zstd|cohort",
  M5="YFB|perplat_zstd|no_coh",     M6="YFB|perplat_zstd|cohort",
  M7="LB|joint_qn_norank|no_coh",   M8="LB|joint_qn_norank|cohort",
  M9="LB|joint_zstd|no_coh",        M10="LB|joint_zstd|cohort",
  M11="LB|log_only|no_coh",         M12="LB|log_only|cohort",
  M13="YFB|joint_qn_norank|no_coh", M14="YFB|joint_qn_norank|cohort",
  M15="YFB|joint_zstd|no_coh",      M16="YFB|joint_zstd|cohort",
  M17="YFB|log_only|no_coh",        M18="YFB|log_only|cohort"
)
results$model_label <- model_labels[results$model]
results$model_label <- factor(results$model_label,
                               levels = model_labels[order(names(model_labels))])
```

Replace the Conclusions section placeholder with actual findings after benchmark results are available:

```markdown
## Conclusions

*(Fill in after full benchmark results are available. Address:)*
- *Which configuration achieves the highest mean external C-index?*
- *Does removing the rank transform (joint_quantile_norank vs joint_quantile_rank) improve or degrade performance?*
- *How does joint z-standardization compare to quantile normalization as a normalization method?*
- *What happens with log_only (no normalization beyond log)? Does platform scale confounding dominate?*
- *Is per-platform z-std still the best single preprocessing decision?*
- *Recommended configuration for manuscript analysis, accounting for both C-index and prediction interpretability (YFB preferred for external scoring: ZF_new = Y_new %*% EF is exact, no approximation required).*
```

- [ ] **Step 2: Verify Quarto renders with existing data (before extended results are available)**

```bash
quarto render docs/reports/merged_benchmark_report.qmd --to pdf 2>&1 | tail -5
```

Expected: renders without error (Conclusions section shows placeholder text).

- [ ] **Step 3: After benchmark results land, fill in the Conclusions section and re-render**

After examining `results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results_extended.csv`:
- Identify the top 3 configurations by mean external C-index
- Compare preprocessing methods within each model type (LB and YFB)
- Note which preprocessing options cause β→0 (log_only on merged data may collapse)
- State the recommended configuration with rationale (include YFB prediction advantage)

```bash
quarto render docs/reports/merged_benchmark_report.qmd 2>&1 | tail -5
git add docs/reports/merged_benchmark_report.qmd docs/reports/merged_benchmark_report.pdf
git commit -m "Update merged benchmark report: 18-config results and final recommended configuration"
```

---

## Task 6: Final documentation updates

- [ ] **Step 1: Update `DECISIONS.md` with extended benchmark findings**

Add a new entry at the top of `DECISIONS.md` (after the header block) with:
- The full 18-row C-index table
- Key comparisons: which preprocessing wins, whether rank transform matters, whether joint z-std matches quantile, whether log_only is viable
- Final recommended configuration with rationale (model type, preprocessing, cohort indicator, K)
- The prediction advantage of YFB noted explicitly

- [ ] **Step 2: Update `ROADMAP.md` — mark extended benchmark complete, note any new items**

Change the "Merged-cohort 6-configuration benchmark" entry from Complete to reference the extended 18-config version. Add any follow-up items (e.g., pathway enrichment for top genes if the winning config has interpretable factors).

- [ ] **Step 3: Update `CLAUDE.md` current model status**

Update the "Current model status" bullet to reflect the final recommended configuration from the extended benchmark.

- [ ] **Step 4: Final commit**

```bash
git add DECISIONS.md ROADMAP.md CLAUDE.md
git commit -m "Document extended benchmark findings: preprocessing comparison and final configuration recommendation"
```

---

## Self-Review

**Spec coverage:**
- ✅ New `normalize_method` parameter with `z_score` and `none` modes (Task 1)
- ✅ Tests for all new normalize_method modes (Task 1)
- ✅ K floor K_final = max(K_1se, 3) applied and logged (Task 3)
- ✅ `set_key()` handles both null and existing integer values (Task 3)
- ✅ All 5 preprocessing options in K-CV loop (Task 3)
- ✅ All 18 model configurations in benchmark loop (Task 4)
- ✅ M5/M6 re-run at new K=3 floor (Task 4 model table)
- ✅ YFB × joint QN excluded with documentation (Task 3 config table)
- ✅ Preprocessing caching avoids redundant computation (Task 4)
- ✅ Extended results saved with `_extended` suffix (Task 4)
- ✅ Report updated for 18 configs (Task 5)
- ✅ DECISIONS.md, ROADMAP.md, CLAUDE.md updated (Task 6)

**Placeholder scan:** Task 5 Conclusions section is intentionally a placeholder pending results. All code blocks are complete.

**Type consistency:** `normalize_method` appears as `mcfg$norm` in Task 4's MODEL_CONFIGS and is passed to `preprocess_merged_cohorts(normalize_method = mcfg$norm)` — consistent with the parameter name defined in Task 1. `set_key()` function defined and used within Task 3 only.
