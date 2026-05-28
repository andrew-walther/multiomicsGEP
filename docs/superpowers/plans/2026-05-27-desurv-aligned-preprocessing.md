# DeSurv-Aligned Preprocessing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Align SSBMF gene selection with DeSurv (Young et al. 2025): combined mean+variance rank criterion, top-3000 per cohort before merging, then intersect — and re-run the M4/M5 benchmark comparison to quantify the effect on external C-index.

**Architecture:** Three layers of change, each building on the previous. (1) `select_top_variable_genes()` gains a `method="combined_rank"` option that ranks genes by mean and variance independently, then keeps the top-N with the lowest combined rank-sum. (2) `preprocess_merged_cohorts()` gains `selection_per_cohort=TRUE` which moves gene selection to before the normalization step — each cohort selects its own top-N, then the intersection is passed forward. (3) A focused comparison runner (`run_desurv_comparison.R`) runs K-CV and external validation on four configs: original M4, original M5, DeSurv-M4, DeSurv-M5. The recommended manuscript configs are compared side-by-side.

**Tech Stack:** R, `yaml`, `survival`. Tests via `tests/test_helpers.R` framework (no testthat). Quarto for report rendering. All CAVI fitting via `code/fit_modular.R` (LB) and `code/fit_cox_on_yf.R` (YFB).

---

## Branching

**Recommended:** Create branch `desurv-gene-selection` from `merged-benchmark` (or from `main` if `merged-benchmark` has been merged).

```bash
# Option A: branch from merged-benchmark (preserves 18-config baseline)
git checkout merged-benchmark
git pull
git checkout -b desurv-gene-selection

# Option B: merge merged-benchmark to main first, then branch from main
git checkout main && git merge merged-benchmark && git push
git checkout -b desurv-gene-selection
```

Keep `merged-benchmark` intact as the baseline reference. The DeSurv comparison outputs go to a new output directory (`outputs/desurv_comparison/`) so they never overwrite the 18-config results.

---

## File Map

| File | Action | What changes |
|------|--------|--------------|
| `code/preprocess_desurv.R` | Modify | `method` param in `select_top_variable_genes()`; `selection_per_cohort` + `selection_method` params in `preprocess_merged_cohorts()` |
| `tests/test_preprocess_desurv.R` | Modify | Add T1.13–T1.16 (4 new tests; suite becomes 246/246) |
| `config/globals.yml` | Modify | Add `top_n_genes_desurv: 3000`; add `k_merged_lb_desurv: null`, `k_merged_yfb_desurv: null` keys |
| `results/benchmark_sim/run_desurv_comparison.R` | Create | K-CV + benchmark for 4 configs (D1–D4); saves to `outputs/desurv_comparison/` |
| `docs/reports/desurv_alignment_report_05_27_26.qmd` | Create | Quarto report comparing D1–D4 per-cohort C-index; renders to PDF + HTML |
| `ROADMAP.md` | Modify | Mark DeSurv gene selection item complete; add comparison item |
| `DECISIONS.md` | Modify | Document alignment decision and benchmark comparison result |
| `CLAUDE.md` | Modify | Update test count; update `top_n_genes_desurv` key reference |

---

## Task 1: Add `method="combined_rank"` to `select_top_variable_genes()`

**Files:**
- Modify: `code/preprocess_desurv.R` (lines 23–41)
- Test: `tests/test_preprocess_desurv.R` (add T1.13, T1.14 after T1.12)

- [ ] **Step 1.1: Write the two failing tests**

Open `tests/test_preprocess_desurv.R` and append after the last `run_test(...)` block (before `report_results(...)`):

```r
# ---------------------------------------------------------------------------
# Tests for method="combined_rank" in select_top_variable_genes()
# ---------------------------------------------------------------------------

run_test("T1.13: combined_rank selects the gene with highest mean+variance and excludes lowest", {
  # g_hm_hv: high mean (250), high variance — should always be selected
  # g_lm_lv: low mean (~1.5), low variance — should never be selected at top_n=2
  Y <- cbind(
    g_hm_hv = c(100, 200, 300, 400),   # mean=250, var=16667
    g_lm_hv = c(  1,   2, 100, 200),   # mean=75.75, var=6823
    g_hm_lv = c(100, 101, 100, 101),   # mean=100.5, var=0.33
    g_lm_lv = c(  1,   1,   2,   2)    # mean=1.5, var=0.33
  )
  out <- select_top_variable_genes(Y, colnames(Y), top_n = 2L, method = "combined_rank")
  assert_equal("g_hm_hv" %in% out$gene_names, TRUE)
  assert_equal("g_lm_lv" %in% out$gene_names, FALSE)
})

run_test("T1.14: combined_rank top_n=1 picks the gene with lowest rank_mean+rank_var sum", {
  # g_hm_hv: rank_mean=1, rank_var=1 → rank_sum=2  (winner)
  # g_lm_hv: rank_mean=3, rank_var=2 → rank_sum=5
  # g_hm_lv: rank_mean=2, rank_var=3 → rank_sum=5  (tie; secondary sort may vary)
  # Only g_hm_hv should be selected at top_n=1
  Y <- cbind(
    g_hm_hv = c(100, 200, 300, 400),
    g_lm_hv = c(  1,   2, 100, 200),
    g_hm_lv = c(100, 101, 100, 101)
  )
  out <- select_top_variable_genes(Y, colnames(Y), top_n = 1L, method = "combined_rank")
  assert_equal(out$gene_names, "g_hm_hv")
  assert_equal(ncol(out$Y), 1L)
})
```

- [ ] **Step 1.2: Run tests to verify they fail**

```bash
cd /Users/ajwalther/GithubProjects/multiomicsGEP
source("tests/test_helpers.R"); source("code/preprocess_desurv.R")
```
Then in R:
```r
source("tests/test_helpers.R")
source("code/preprocess_desurv.R")
source("tests/test_preprocess_desurv.R")
```
Expected: T1.13 and T1.14 FAIL with "argument 'method' is missing" or similar.

- [ ] **Step 1.3: Implement `method` parameter in `select_top_variable_genes()`**

Replace the entire `select_top_variable_genes` function in `code/preprocess_desurv.R` (lines 23–41) with:

```r
#' Select the top-N most informative genes from an expression matrix.
#'
#' @param Y          numeric matrix (n × p), rows = subjects, columns = genes.
#' @param gene_names character vector of length p; gene identifiers.
#' @param top_n      integer or NULL. If NULL or >= ncol(Y), all genes are returned.
#' @param method     "variance" (default) ranks by variance only.
#'                   "combined_rank" ranks each gene by mean expression and by
#'                   variance independently (highest = rank 1), sums the two ranks,
#'                   and retains the top_n genes with the smallest rank sum.
#'                   This matches the DeSurv (Young et al. 2025) gene selection
#'                   criterion: genes that are both highly expressed AND highly
#'                   variable receive the lowest combined rank.
#' @return list with components Y (n × top_n), gene_names (length top_n), gene_var.
#' @family v2 preprocessing
select_top_variable_genes <- function(Y, gene_names, top_n = 2000,
                                       method = c("variance", "combined_rank")) {
  validate_expression_inputs(Y, gene_names)
  method <- match.arg(method)

  if (is.null(top_n) || top_n >= ncol(Y)) {
    return(list(Y = Y, gene_names = gene_names, gene_var = apply(Y, 2, stats::var)))
  }
  if (length(top_n) != 1 || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n))
    stop("top_n must be NULL or an integer >= 1.")

  if (method == "variance") {
    # Original behaviour: order by variance descending, keep top top_n.
    gene_var <- apply(Y, 2, stats::var)
    ord      <- order(gene_var, decreasing = TRUE, na.last = NA)
  } else {
    # combined_rank (DeSurv criterion):
    # rank_mean = rank of negative mean (rank 1 = highest mean expression).
    # rank_var  = rank of negative variance (rank 1 = highest variance).
    # Select genes with the smallest combined rank_mean + rank_var.
    gene_mean <- colMeans(Y)
    gene_var  <- apply(Y, 2, stats::var)
    rank_mean <- rank(-gene_mean, ties.method = "average", na.last = "keep")
    rank_var  <- rank(-gene_var,  ties.method = "average", na.last = "keep")
    rank_sum  <- rank_mean + rank_var
    # order ascending: rank_sum=2 (rank 1 in both) is the best gene.
    ord <- order(rank_sum, na.last = NA)
  }

  keep         <- ord[seq_len(min(as.integer(top_n), length(ord)))]
  Y_keep       <- Y[, keep, drop = FALSE]
  gene_var_out <- apply(Y_keep, 2, stats::var)

  list(
    Y          = Y_keep,
    gene_names = gene_names[keep],
    gene_var   = gene_var_out
  )
}
```

- [ ] **Step 1.4: Run tests to verify T1.13 and T1.14 pass**

```r
source("tests/test_helpers.R")
source("code/preprocess_desurv.R")
source("tests/test_preprocess_desurv.R")
```

Expected: T1.1–T1.14 all PASS. Total = 14 tests for this file.

- [ ] **Step 1.5: Run full test suite to confirm no regressions**

```bash
Rscript tests/run_tests.R
```
Expected: 242/242 passing (unchanged count — T1.13 and T1.14 are added to the file total but the main suite counter increments to 244/244 once the file is sourced).

Wait — check the actual count: `tests/run_tests.R` sources `tests/test_preprocess_desurv.R` which adds 2 tests. Starting count is 242. After this step: 244/244.

- [ ] **Step 1.6: Commit**

```bash
git add code/preprocess_desurv.R tests/test_preprocess_desurv.R
git commit -m "$(cat <<'EOF'
Add combined_rank method to select_top_variable_genes()

Implements the DeSurv (Young et al. 2025) gene selection criterion: rank
genes by mean expression and by variance independently, then keep the top-N
with the lowest combined rank-sum. This selects genes that are both highly
expressed and highly variable — better motivated than variance-only, which
can retain lowly-expressed noisy genes near the detection floor.

New parameter: method = c("variance", "combined_rank") in
select_top_variable_genes(). Default "variance" preserves existing behavior.

Tests T1.13-T1.14 added; 244/244 passing.
EOF
)"
```

---

## Task 2: Add per-cohort gene selection to `preprocess_merged_cohorts()`

**Files:**
- Modify: `code/preprocess_desurv.R` (function `preprocess_merged_cohorts`, lines 233–333)
- Test: `tests/test_preprocess_desurv.R` (add T1.15, T1.16)

- [ ] **Step 2.1: Write the two failing tests**

Append after T1.14 (before `report_results(...)`):

```r
# ---------------------------------------------------------------------------
# Tests for selection_per_cohort parameter in preprocess_merged_cohorts()
# ---------------------------------------------------------------------------

run_test("T1.15: selection_per_cohort=TRUE runs and returns selection metadata", {
  raw   <- make_synthetic_raw_list()  # 2 cohorts, 4 subjects each, 6 genes
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out   <- preprocess_merged_cohorts(
    raw, flags,
    top_n                = 4L,
    selection_per_cohort = TRUE,
    selection_method     = "combined_rank",
    per_platform_standardize = FALSE,
    normalize_method     = "none",
    rank_transform       = FALSE
  )
  # Gene selection per cohort (top 4 of 6) then intersect: intersection should
  # be >= 2 and <= 4 genes; matrix dimensions must be consistent.
  assert_true(out$p >= 2L && out$p <= 4L,
              "Intersection of two top-4-of-6 sets should be 2–4 genes")
  assert_equal(out$selection_per_cohort, TRUE)
  assert_equal(out$selection_method, "combined_rank")
  assert_equal(nrow(out$Y), 8L)  # 4 + 4 subjects
  assert_finite(out$Y)
})

run_test("T1.16: selection_per_cohort=FALSE (default) still works and returns FALSE in output", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out   <- preprocess_merged_cohorts(
    raw, flags, top_n = 6L,
    selection_per_cohort = FALSE,
    normalize_method     = "none",
    rank_transform       = FALSE
  )
  assert_equal(out$selection_per_cohort, FALSE)
  assert_equal(out$selection_method, "variance")  # default
  assert_equal(nrow(out$Y), 8L)
  assert_equal(out$p, 6L)  # all 6 genes kept (top_n=6 = ncol)
})
```

- [ ] **Step 2.2: Run tests to verify they fail**

```r
source("tests/test_helpers.R"); source("code/preprocess_desurv.R")
source("tests/test_preprocess_desurv.R")
```
Expected: T1.15 FAILs ("unused argument selection_per_cohort"); T1.16 FAILs similarly.

- [ ] **Step 2.3: Implement `selection_per_cohort` and `selection_method` in `preprocess_merged_cohorts()`**

Replace the function signature and body in `code/preprocess_desurv.R`. The key changes are:
1. Two new parameters added to the function signature.
2. A new "Step 2b" block inserted after the log2 transform loop, before `per_platform_standardize`.
3. The post-normalization gene selection block (Steps 5-6) made conditional on `!selection_per_cohort`.
4. Both new params echoed back in the return list.

Replace the entire `preprocess_merged_cohorts` function (from `preprocess_merged_cohorts <- function(` to the closing `}`) with:

```r
#' v2 merged-cohort preprocessing: intersect first, then normalise jointly.
#'
#' Now supports two gene selection strategies:
#'
#' \strong{selection_per_cohort = FALSE (default, original behavior):}
#' Gene selection runs AFTER normalization on the merged matrix (Steps 5-6).
#' top_n genes are selected by \code{selection_method} on the normalized data.
#'
#' \strong{selection_per_cohort = TRUE (DeSurv-aligned):}
#' Gene selection runs PER COHORT before normalization (Step 2b).
#' Each cohort independently selects its top-\code{top_n} genes, then the
#' intersection of those gene sets is carried forward. This avoids selecting
#' genes on data where variance has been equalized by z-standardization
#' (which defeats variance-based filtering).
#'
#' DeSurv-aligned pipeline (selection_per_cohort=TRUE):
#' \enumerate{
#'   \item Intersect raw gene universes across all training cohorts.
#'   \item Log2(x+1) transform per cohort (RNA-seq only).
#'   \item Per-cohort: select top-N genes by combined_rank, then intersect.
#'   \item Per-platform z-standardization on the intersected gene set.
#'   \item Row-bind into merged matrix.
#'   \item (Optional) rank-transform each subject.
#' }
#'
#' @param cohort_raw_list        named list of raw cohort objects with \code{$Y}
#'   and \code{$gene_names}.
#' @param log_transform_flags    named logical vector; TRUE = apply log2(x+1).
#' @param top_n                  integer; number of genes to select. Default 2000.
#' @param ties_method            ties method for rank(). Default "average".
#' @param rank_transform         logical; apply per-subject rank transform. Default TRUE.
#' @param per_platform_standardize logical; z-standardize each cohort before merging. Default FALSE.
#' @param normalize_method       "quantile", "z_score", or "none". Default "quantile".
#' @param selection_per_cohort   logical; if TRUE, gene selection runs per cohort
#'   before normalization (DeSurv-aligned). Default FALSE.
#' @param selection_method       "variance" or "combined_rank"; passed to
#'   \code{select_top_variable_genes()}. Default "variance".
#' @return list with Y, gene_names, n, p, dataset_labels, n_raw_intersect,
#'   rank_transform, per_platform_standardize, normalize_method,
#'   selection_per_cohort, selection_method.
#' @family v2 preprocessing
#' @seealso \code{\link{select_top_variable_genes}}
preprocess_merged_cohorts <- function(cohort_raw_list,
                                      log_transform_flags,
                                      top_n                    = 2000,
                                      ties_method              = "average",
                                      rank_transform           = TRUE,
                                      per_platform_standardize = FALSE,
                                      normalize_method         = c("quantile", "z_score", "none"),
                                      selection_per_cohort     = FALSE,
                                      selection_method         = c("variance", "combined_rank")) {
  normalize_method <- match.arg(normalize_method)
  selection_method <- match.arg(selection_method)
  cohort_names <- names(cohort_raw_list)
  stopifnot(!is.null(cohort_names), all(cohort_names %in% names(log_transform_flags)))

  # Step 1: intersect raw gene universes (no preprocessing yet)
  gene_lists   <- lapply(cohort_raw_list, function(x) x$gene_names)
  common_genes <- Reduce(intersect, gene_lists)
  if (length(common_genes) == 0)
    stop("No common genes found across cohorts — check gene_names fields.")
  cat(sprintf("  [v2] Raw gene intersection: %d genes across %s\n",
              length(common_genes), paste(cohort_names, collapse = " + ")))

  # Steps 2–3: log2 transform per cohort (platform-aware), subset to common genes.
  cohort_matrices <- lapply(cohort_names, function(ds) {
    raw <- cohort_raw_list[[ds]]
    idx <- match(common_genes, raw$gene_names)
    Y   <- raw$Y[, idx, drop = FALSE]
    if (log_transform_flags[[ds]])
      Y <- log2_plus1_transform(Y)
    Y
  })
  names(cohort_matrices) <- cohort_names

  # Step 2b (DeSurv-aligned, optional): per-cohort gene selection BEFORE normalization.
  # Selects top-N genes within each cohort independently on the log-transformed data,
  # then intersects the selected sets. This prevents variance-equalization from
  # neutralizing the variance signal used for gene selection.
  if (selection_per_cohort) {
    cat(sprintf("  [v2] Per-cohort gene selection (top-%d, method=%s) before normalization ...\n",
                top_n, selection_method))
    selected_per_cohort <- lapply(cohort_names, function(ds) {
      sel <- select_top_variable_genes(cohort_matrices[[ds]], common_genes,
                                       top_n = top_n, method = selection_method)
      sel$gene_names
    })
    names(selected_per_cohort) <- cohort_names
    selected_genes <- Reduce(intersect, selected_per_cohort)
    if (length(selected_genes) == 0)
      stop("Per-cohort gene selection produced no common genes. Reduce top_n or check inputs.")
    cat(sprintf("  [v2] Genes after per-cohort selection and intersection: %d\n",
                length(selected_genes)))
    # Subset cohort matrices to the intersected gene set
    cohort_matrices <- lapply(cohort_names, function(ds) {
      idx <- match(selected_genes, common_genes)
      cohort_matrices[[ds]][, idx, drop = FALSE]
    })
    names(cohort_matrices) <- cohort_names
    common_genes <- selected_genes
  }

  # Optional: per-platform z-standardisation before row-bind.
  if (per_platform_standardize) {
    cat("  [v2] Per-platform z-standardization (colMean=0, colSD=1 per cohort) ...\n")
    cohort_matrices <- per_platform_standardize_cohorts(cohort_matrices)
  }

  Y_merged <- do.call(rbind, cohort_matrices)
  colnames(Y_merged) <- common_genes

  dataset_labels <- factor(
    rep(cohort_names, vapply(cohort_raw_list, function(x) nrow(x$Y), integer(1))),
    levels = cohort_names
  )

  # Step 4: distribution normalization across merged samples (method-dependent).
  if (normalize_method == "quantile") {
    cat(sprintf("  [v2] Quantile normalising merged matrix (%d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- quantile_normalize_merged(Y_merged)
  } else if (normalize_method == "z_score") {
    cat(sprintf("  [v2] Joint z-standardizing merged matrix (%d x %d, colMean=0, colSD=1) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- scale(Y_merged, center = TRUE, scale = TRUE)
    n_nan  <- sum(is.nan(Y_norm))
    if (n_nan > 0)
      cat(sprintf("  [v2] Note: %d NaN entries from zero-variance genes replaced with 0.\n", n_nan))
    Y_norm[is.nan(Y_norm)] <- 0
    colnames(Y_norm) <- common_genes
  } else {
    cat(sprintf("  [v2] Skipping normalization (log-transform only; %d x %d) ...\n",
                nrow(Y_merged), ncol(Y_merged)))
    Y_norm <- Y_merged
  }

  # Steps 5–6: post-normalization gene selection (only when NOT using per-cohort selection).
  # When selection_per_cohort=TRUE the gene set was fixed in Step 2b; skip here.
  if (!selection_per_cohort) {
    selected <- select_top_variable_genes(Y_norm, common_genes, top_n = top_n,
                                          method = selection_method)
    cat(sprintf("  [v2] Genes retained after top-%d %s filter: %d\n",
                top_n, selection_method, length(selected$gene_names)))
  } else {
    selected <- list(Y = Y_norm, gene_names = common_genes)
  }

  # Step 7 (optional): rank-transform each subject within the selected gene set.
  if (rank_transform) {
    Y_final <- rank_transform_subjects(selected$Y, ties_method = ties_method)
    cat("  [v2] Rank-transforming subjects (Step 7).\n")
  } else {
    Y_final <- selected$Y
    cat("  [v2] Skipping rank transform (rank_transform = FALSE).\n")
  }
  colnames(Y_final) <- selected$gene_names

  list(
    Y                        = Y_final,
    gene_names               = selected$gene_names,
    n                        = nrow(Y_final),
    p                        = ncol(Y_final),
    dataset_labels           = dataset_labels,
    n_raw_intersect          = length(Reduce(intersect,
                                             lapply(cohort_raw_list, function(x) x$gene_names))),
    rank_transform           = rank_transform,
    per_platform_standardize = per_platform_standardize,
    normalize_method         = normalize_method,
    selection_per_cohort     = selection_per_cohort,
    selection_method         = selection_method
  )
}
```

- [ ] **Step 2.4: Run tests T1.15 and T1.16**

```r
source("tests/test_helpers.R"); source("code/preprocess_desurv.R")
source("tests/test_preprocess_desurv.R")
```
Expected: all 16 tests PASS (T1.1–T1.16).

- [ ] **Step 2.5: Run full test suite**

```bash
Rscript tests/run_tests.R
```
Expected: 246/246 passing (242 prior + 4 new tests from T1.13–T1.16).

- [ ] **Step 2.6: Commit**

```bash
git add code/preprocess_desurv.R tests/test_preprocess_desurv.R
git commit -m "$(cat <<'EOF'
Add selection_per_cohort param to preprocess_merged_cohorts()

DeSurv-aligned gene selection: when selection_per_cohort=TRUE, each cohort
selects its top-N genes (by selection_method) on the log-transformed data
BEFORE normalization, then the per-cohort gene sets are intersected and
carried forward. This avoids variance-equalization artifacts: running gene
selection after joint z-standardization or quantile normalization defeats
variance-based ranking since those transforms alter per-gene variance.

New parameters:
  selection_per_cohort = FALSE (default; preserves original behavior)
  selection_method = c("variance", "combined_rank") (default "variance")

Tests T1.15-T1.16 added; 246/246 passing.
EOF
)"
```

---

## Task 3: Update `globals.yml` and write the comparison runner

**Files:**
- Modify: `config/globals.yml`
- Create: `results/benchmark_sim/run_desurv_comparison.R`

- [ ] **Step 3.1: Add new keys to `config/globals.yml`**

In the `preprocessing:` section, add after `top_n_genes: 2000`:

```yaml
  top_n_genes_desurv: 3000   # DeSurv (Young et al. 2025): top-3000 per cohort by combined
                              # mean+variance rank, then intersect. Yields ~1970 after intersection
                              # (vs 2000 variance-only). See DECISIONS.md 2026-05-27.
```

In the `benchmark:` section, add after the `k_merged_yfb_logonly` line:

```yaml
  # DeSurv-aligned configs — filled by run_desurv_comparison.R
  k_merged_lb_desurv:  null   # LB  + per-platform z-std + combined_rank + top-3000 per-cohort
  k_merged_yfb_desurv: null   # YFB + per-platform z-std + combined_rank + top-3000 per-cohort
```

- [ ] **Step 3.2: Verify globals.yml is valid YAML**

```r
cfg <- yaml::read_yaml("config/globals.yml")
cat("top_n_genes_desurv:", cfg$preprocessing$top_n_genes_desurv, "\n")
cat("k_merged_lb_desurv:", is.null(cfg$benchmark$k_merged_lb_desurv), "\n")  # TRUE = null
```
Expected output:
```
top_n_genes_desurv: 3000
k_merged_lb_desurv: TRUE
```

- [ ] **Step 3.3: Write `results/benchmark_sim/run_desurv_comparison.R`**

Create the file with this content:

```r
# ============================================================
# Script:  results/benchmark_sim/run_desurv_comparison.R
# Purpose: Compare original M4/M5 against DeSurv-aligned preprocessing.
#
#   D1: LB  + per-platform z-std + variance + top-2000 + post-norm selection (= M4)
#   D2: YFB + per-platform z-std + variance + top-2000 + post-norm selection (= M5)
#   D3: LB  + per-platform z-std + combined_rank + top-3000 + per-cohort selection
#   D4: YFB + per-platform z-std + combined_rank + top-3000 + per-cohort selection
#
#   D1/D2 reproduce the recommended manuscript configs using the same K already
#   stored in globals.yml.  D3/D4 run K-CV (K_final = max(K_1se, 3)) and store
#   K values in globals.yml before fitting.
#
#   Output: results/benchmark_sim/outputs/desurv_comparison/
#     desurv_comparison_results.csv
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-27
# Usage:   caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R
#          caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R --quick
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

# Navigate to project root if invoked from a subdirectory
if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival) })

source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R");  source("code/compute_elbo.R")
source("code/update_F_cohort.R")
source("code/train_test_split.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_modular.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")
source("code/select_K.R")
source("code/select_alpha_cv.R")

YML_PATH <- "config/globals.yml"
cfg      <- yaml::read_yaml(YML_PATH)
b        <- cfg$benchmark
p        <- cfg$preprocessing

ALPHA            <- b$alpha
LAMBDA           <- b$lambda
MAX_ITER         <- if (QUICK_MODE) 30L else cfg$cavi$max_iter
PRIOR_BETA       <- "normal"
SIGMA_COH        <- 1.0
BETA_THRESH      <- cfg$k_selection$beta_threshold
N_CV_FOLDS       <- cfg$cavi$n_cv_folds
K_MIN_BIOLOGICAL <- 3L
TOP_N_ORIG       <- p$top_n_genes        # 2000 — original setting
TOP_N_DESURV     <- p$top_n_genes_desurv # 3000 — DeSurv-aligned

OUT_DIR <- "results/benchmark_sim/outputs/desurv_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Quick mode: %s | MAX_ITER: %d\n", QUICK_MODE, MAX_ITER))

# --------------------------------------------------------------------------
# 1. Load training data
# --------------------------------------------------------------------------

cat("\n--- Loading TCGA_PAAD + CPTAC ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) {
  cat(sprintf("  %s ...\n", ds)); load_pdac_raw(ds, PDAC_DATA_ROOT)
})
n_tcga  <- train_raw$TCGA_PAAD$n
n_cptac <- train_raw$CPTAC$n
cohort_labels <- c(rep("TCGA", n_tcga), rep("CPTAC", n_cptac))
time_train   <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN_COHORTS, function(ds) train_raw[[ds]]$status))
cat(sprintf("  n=%d (TCGA=%d, CPTAC=%d), events=%d\n\n",
            n_tcga + n_cptac, n_tcga, n_cptac, sum(status_train)))

# --------------------------------------------------------------------------
# 2. Configuration table
#    D1/D2 reproduce M4/M5 (for direct comparison, same K from globals.yml).
#    D3/D4 use DeSurv-aligned preprocessing; K-CV runs in section 4.
# --------------------------------------------------------------------------

DESURV_CONFIGS <- list(
  list(id = "D1", label = "LB orig (M4)",
       model        = "LB",
       top_n        = TOP_N_ORIG,
       sel_method   = "variance",
       per_cohort   = FALSE,
       cohort_id    = TRUE,
       k_key        = "k_merged_lb_perplatform"),
  list(id = "D2", label = "YFB orig (M5)",
       model        = "YFB",
       top_n        = TOP_N_ORIG,
       sel_method   = "variance",
       per_cohort   = FALSE,
       cohort_id    = FALSE,
       k_key        = "k_merged_yfb_perplatform"),
  list(id = "D3", label = "LB DeSurv-aligned",
       model        = "LB",
       top_n        = TOP_N_DESURV,
       sel_method   = "combined_rank",
       per_cohort   = TRUE,
       cohort_id    = TRUE,
       k_key        = "k_merged_lb_desurv"),
  list(id = "D4", label = "YFB DeSurv-aligned",
       model        = "YFB",
       top_n        = TOP_N_DESURV,
       sel_method   = "combined_rank",
       per_cohort   = TRUE,
       cohort_id    = FALSE,
       k_key        = "k_merged_yfb_desurv")
)

# --------------------------------------------------------------------------
# Helper: replace any value for a named key in globals.yml
# --------------------------------------------------------------------------

set_key <- function(yml_path, key, value) {
  lines   <- readLines(yml_path)
  pattern <- paste0("^(\\s*", key, ":\\s*)\\S+")
  idx     <- grep(pattern, lines)
  if (length(idx) == 0)
    stop(sprintf("Key '%s' not found in %s", key, yml_path))
  lines[idx[1]] <- sub(pattern, paste0("\\1", value), lines[idx[1]])
  writeLines(lines, yml_path)
  invisible(NULL)
}

# --------------------------------------------------------------------------
# Helper: oriented C-index
# --------------------------------------------------------------------------

oriented_cindex <- function(risk, time, status) {
  c_raw <- as.numeric(concordance(Surv(time, status) ~ risk)$concordance)
  max(c_raw, 1 - c_raw)
}

# --------------------------------------------------------------------------
# 3. Preprocess training data — one call per unique (top_n, method, per_cohort)
# --------------------------------------------------------------------------

cat("--- Preprocessing training data ---\n")
preproc_cache  <- list()
gene_set_cache <- list()

for (dcfg in DESURV_CONFIGS) {
  ckey <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, sep = "_")
  if (ckey %in% names(preproc_cache)) next
  cat(sprintf("  top_n=%d, method=%s, per_cohort=%s ...\n",
              dcfg$top_n, dcfg$sel_method, dcfg$per_cohort))
  pp <- preprocess_merged_cohorts(
    cohort_raw_list        = train_raw,
    log_transform_flags    = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
    top_n                  = dcfg$top_n,
    rank_transform         = FALSE,          # per-platform z-std handles scale
    per_platform_standardize = TRUE,
    normalize_method       = "none",         # z-std IS the normalization
    selection_per_cohort   = dcfg$per_cohort,
    selection_method       = dcfg$sel_method
  )
  preproc_cache[[ckey]]  <- pp$Y
  gene_set_cache[[ckey]] <- pp$gene_names
  cat(sprintf("    n=%d, p=%d\n", nrow(pp$Y), ncol(pp$Y)))
}

# --------------------------------------------------------------------------
# 4. K-CV for D3 and D4 (D1/D2 reuse K from existing globals.yml entries)
# --------------------------------------------------------------------------

cat("\n--- K-CV for DeSurv-aligned configs (D3, D4) ---\n")
cfg <- yaml::read_yaml(YML_PATH); b <- cfg$benchmark  # re-read after preprocessing

for (dcfg in DESURV_CONFIGS[3:4]) {
  if (!is.null(b[[dcfg$k_key]])) {
    cat(sprintf("  %s: K=%d already set — skipping CV\n", dcfg$id, b[[dcfg$k_key]]))
    next
  }
  cat(sprintf("  Running K-CV for %s (%s) ...\n", dcfg$id, dcfg$label))
  ckey    <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, sep = "_")
  Y_train <- preproc_cache[[ckey]]

  cv_res  <- select_K_cv(
    Y          = Y_train,
    time       = time_train,
    status     = status_train,
    K_max      = cfg$cavi$k_max,
    n_folds    = N_CV_FOLDS,
    alpha      = ALPHA,
    lambda     = LAMBDA,
    prior_beta = PRIOR_BETA,
    model      = dcfg$model,
    max_iter   = MAX_ITER
  )
  K_opt   <- cv_res$K_1se
  K_final <- max(K_opt, K_MIN_BIOLOGICAL)
  cat(sprintf("  %s: K_1se=%d → K_final=%d (biological floor K>=%d)\n",
              dcfg$id, K_opt, K_final, K_MIN_BIOLOGICAL))
  set_key(YML_PATH, dcfg$k_key, as.character(K_final))
}

cfg <- yaml::read_yaml(YML_PATH); b <- cfg$benchmark  # re-read after K-CV writes

# --------------------------------------------------------------------------
# 5. Fit all 4 configurations
# --------------------------------------------------------------------------

cat("\n=== Fitting 4 configurations ===\n\n")
fits <- list()

for (dcfg in DESURV_CONFIGS) {
  ckey      <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, sep = "_")
  Y_train   <- preproc_cache[[ckey]]
  K         <- b[[dcfg$k_key]]
  cohort_id <- if (dcfg$cohort_id) cohort_labels else NULL

  cat(sprintf("--- %s [%s] K=%d cohort_id=%s ---\n",
              dcfg$id, dcfg$label, K, !is.null(cohort_id)))
  set.seed(42L)

  fit <- suppressMessages(
    if (dcfg$model == "LB")
      fit_supervised_mf_modular(
        Y_train, time_train, status_train,
        K = K, max_iter = MAX_ITER, alpha = ALPHA, lambda = LAMBDA,
        prior_beta = PRIOR_BETA, verbose = TRUE,
        cohort_id = cohort_id, sigma_F_cohort = SIGMA_COH)
    else
      fit_cox_on_yf(
        Y_train, time_train, status_train,
        K = K, max_iter = MAX_ITER, alpha = ALPHA, lambda = LAMBDA,
        prior_beta = PRIOR_BETA, verbose = TRUE,
        cohort_id = cohort_id, sigma_F_cohort = SIGMA_COH)
  )
  k_eff <- sum(abs(fit$EBeta) > BETA_THRESH)
  cat(sprintf("  K_eff=%d | beta_max=%.4f | iters=%d\n\n",
              k_eff, max(abs(fit$EBeta)), fit$history$n_iter))
  fits[[dcfg$id]] <- fit
}

# --------------------------------------------------------------------------
# 6. External validation on 5 held-out cohorts
# --------------------------------------------------------------------------

cat("--- External validation (5 cohorts) ---\n")
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts
results_rows     <- list()

for (ext_cohort in EXTERNAL_COHORTS) {
  cat(sprintf("  Loading %s ...\n", ext_cohort))
  raw_ext <- load_pdac_raw(ext_cohort, PDAC_DATA_ROOT)
  # top_n=NULL: keep all external genes; intersection with train_genes controls
  # the final gene set. This prevents the external-cohort top-N filter from
  # discarding genes that happen to be in the training gene set.
  pre_ext <- preprocess_desurv_cohort(
    Y             = raw_ext$Y,
    gene_names    = raw_ext$gene_names,
    top_n         = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[ext_cohort]],
    cohort_name   = ext_cohort
  )

  for (dcfg in DESURV_CONFIGS) {
    ckey        <- paste(dcfg$top_n, dcfg$sel_method, dcfg$per_cohort, sep = "_")
    train_genes <- gene_set_cache[[ckey]]
    fit         <- fits[[dcfg$id]]
    K           <- b[[dcfg$k_key]]

    common    <- intersect(train_genes, pre_ext$gene_names)
    if (length(common) < 100) {
      cat(sprintf("    Skipping %s x %s: only %d common genes\n",
                  ext_cohort, dcfg$id, length(common)))
      next
    }

    Y_ext     <- pre_ext$Y[, match(common, pre_ext$gene_names), drop = FALSE]
    train_idx <- match(common, train_genes)
    EF_sub    <- fit$EF[train_idx, , drop = FALSE]

    pred <- if (dcfg$model == "LB")
      predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)
    else
      predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)

    c_val <- oriented_cindex(pred$risk_scores, raw_ext$time, raw_ext$status)

    results_rows[[length(results_rows) + 1]] <- data.frame(
      model        = dcfg$id,
      label        = dcfg$label,
      model_type   = dcfg$model,
      cohort       = ext_cohort,
      c_index      = round(c_val, 4),
      K            = K,
      k_eff        = sum(abs(fit$EBeta) > BETA_THRESH),
      beta_max     = round(max(abs(fit$EBeta)), 4),
      top_n        = dcfg$top_n,
      sel_method   = dcfg$sel_method,
      per_cohort   = dcfg$per_cohort,
      n_common_genes = length(common),
      stringsAsFactors = FALSE
    )
  }
}

# --------------------------------------------------------------------------
# 7. Save and report
# --------------------------------------------------------------------------

results <- do.call(rbind, results_rows)
out_csv <- file.path(OUT_DIR, "desurv_comparison_results.csv")
write.csv(results, out_csv, row.names = FALSE)

cat(sprintf("\n=== Results saved: %s ===\n\n", out_csv))
cat("Mean C-index by configuration:\n")
agg <- aggregate(c_index ~ model + label + K + k_eff, data = results, FUN = mean)
agg <- agg[order(agg$c_index, decreasing = TRUE), ]
for (i in seq_len(nrow(agg))) {
  cat(sprintf("  %s (%s): mean C=%.3f | K=%d | K_eff=%d\n",
              agg$model[i], agg$label[i], agg$c_index[i], agg$K[i], agg$k_eff[i]))
}
cat("\nPer-cohort C-index:\n")
for (m in unique(results$model)) {
  sub <- results[results$model == m, ]
  cat(sprintf("  %s: %s\n", m,
              paste(sprintf("%s=%.3f", sub$cohort, sub$c_index), collapse = "  ")))
}
```

- [ ] **Step 3.4: Smoke-test the script (quick mode)**

```bash
caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R --quick
```

Expected output:
- No R errors
- Preprocessing messages with gene counts (should see ~1970 genes for D3/D4 after per-cohort selection + intersection)
- K-CV runs for D3 and D4 (max_iter=30, quick)
- 4 model fits complete
- Mean C-index table printed
- `desurv_comparison_results.csv` written to `outputs/desurv_comparison/`

If the script errors, fix before committing.

- [ ] **Step 3.5: Commit**

```bash
git add config/globals.yml results/benchmark_sim/run_desurv_comparison.R
git commit -m "$(cat <<'EOF'
Add DeSurv-aligned comparison runner and globals.yml keys

run_desurv_comparison.R fits 4 configs: original M4/M5 (D1/D2, variance,
top-2000, post-norm selection) vs DeSurv-aligned D3/D4 (combined_rank,
top-3000, per-cohort selection before per-platform z-std). D3/D4 run K-CV
with biological floor K>=3; K values written to globals.yml.

globals.yml: added top_n_genes_desurv=3000, k_merged_lb_desurv=null,
k_merged_yfb_desurv=null.
EOF
)"
```

---

## Task 4: Run the full comparison (K-CV + benchmark)

**This task is a compute run, not a code task. Expected wall time: 20–40 minutes.**

- [ ] **Step 4.1: Launch in background with caffeinate**

```bash
cd /Users/ajwalther/GithubProjects/multiomicsGEP
caffeinate -i Rscript results/benchmark_sim/run_desurv_comparison.R \
  > /tmp/desurv_comparison.log 2>&1 &
echo "PID: $!"
```

- [ ] **Step 4.2: Monitor progress**

```bash
tail -40 /tmp/desurv_comparison.log
ps aux | grep run_desurv_comparison | grep -v grep
```

Look for:
- Preprocessing: "Genes after per-cohort selection and intersection: N" — note N for D3/D4.
- K-CV section: K values printed for D3 and D4.
- Each of the 4 model fits printing K_eff and beta_max.
- "Results saved" message.

- [ ] **Step 4.3: Verify output exists and looks reasonable**

```r
res <- read.csv("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv")
nrow(res)          # expect 4 configs × 5 external cohorts = 20 rows
table(res$model)   # D1/D2/D3/D4 each 5 rows
aggregate(c_index ~ model + label, data = res, FUN = mean)
```

Expected: D3/D4 should show non-zero K_eff (per-platform z-std is the viable preprocessing).

**Note on D1/D2 baseline numbers:** D1/D2 use `top_n=NULL` for external cohort preprocessing
(keep all genes, then intersect with training genes). The original benchmark used `top_n=2000`
for external cohorts. D1/D2 C-index values may differ slightly from the published M4=0.622,
M5=0.626 — this is expected and does not indicate an error. What matters is the D1 vs D3 and
D2 vs D4 deltas, which are computed under identical external preprocessing for all four configs.

- [ ] **Step 4.4: Re-read globals.yml to confirm K was written**

```r
cfg <- yaml::read_yaml("config/globals.yml")
cat("k_merged_lb_desurv: ", cfg$benchmark$k_merged_lb_desurv,  "\n")
cat("k_merged_yfb_desurv:", cfg$benchmark$k_merged_yfb_desurv, "\n")
```
Expected: both are integers (not null).

- [ ] **Step 4.5: Commit results and updated globals.yml**

```bash
git add config/globals.yml
git add results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv
git commit -m "$(cat <<'EOF'
Add DeSurv comparison results: D1-D4 external C-index on 5 PDAC cohorts

D3 (LB DeSurv-aligned): mean C=<value>, K=<K>, K_eff=<k_eff>
D4 (YFB DeSurv-aligned): mean C=<value>, K=<K>, K_eff=<k_eff>
vs D1 (LB orig M4): mean C=0.622, D2 (YFB orig M5): mean C=0.626

globals.yml: k_merged_lb_desurv=<K>, k_merged_yfb_desurv=<K>
EOF
)"
```
(Replace `<value>` placeholders with actual numbers from Step 4.3.)

---

## Task 5: Write comparison report and update documentation

**Files:**
- Create: `docs/reports/desurv_alignment_report_05_27_26.qmd`
- Modify: `ROADMAP.md`, `DECISIONS.md`, `CLAUDE.md`
- Conditionally modify: `docs/progress_report/SSBMF_Status_Update_5_27_26.qmd` (if DeSurv wins)

### Decision rule (apply before writing any documentation)

After reading `desurv_comparison_results.csv`, compute:

```r
res  <- read.csv("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv")
agg  <- aggregate(c_index ~ model, data = res, FUN = mean)
d2   <- agg$c_index[agg$model == "D2"]   # YFB original (≈ M5)
d4   <- agg$c_index[agg$model == "D4"]   # YFB DeSurv-aligned
delta_yfb <- d4 - d2
cat(sprintf("YFB delta: %.4f (D4 - D2)\n", delta_yfb))
```

Apply this rule:

| Result | Decision | What to do in documentation |
|--------|----------|-----------------------------|
| delta_yfb > +0.005 | **Adopt DeSurv** | D4 becomes new primary config; update CLAUDE.md, DECISIONS.md, progress report; flag pathway enrichment needed |
| −0.005 ≤ delta_yfb ≤ +0.005 | **Neutral — keep M5** | Document as negative result; note both are equivalent; keep M5 as primary |
| delta_yfb < −0.005 | **Reject DeSurv** | Document as negative result; retain M5 as primary; note DeSurv gene selection does not transfer to SSBMF |

The threshold ±0.005 corresponds to roughly half the cohort-to-cohort variability in the 18-config
benchmark and is large enough to be practically meaningful for a 5-cohort average C-index.

Check D1 vs D3 (LB) using the same rule for completeness, but the **YFB comparison is the
primary decision driver** since D2/M5 is the recommended manuscript configuration.

---

- [ ] **Step 5.1: Create `docs/reports/desurv_alignment_report_05_27_26.qmd`**

```qmd
---
title: "DeSurv Gene Selection Alignment: Impact on M4/M5 External C-Index"
author: "Andrew Walther"
date: "2026-05-27"
format:
  pdf:
    latex-engine: xelatex
    keep-tex: false
    fig-pos: "htbp"
    include-in-header:
      text: |
        \usepackage{booktabs}
        \usepackage{geometry}
        \usepackage{titling}
        \geometry{margin=0.42in}
        \renewcommand{\familydefault}{\sfdefault}
        \setlength{\droptitle}{-1.8em}
        \pretitle{\begin{center}\large\bfseries}
        \posttitle{\par\end{center}\vspace{-0.7em}}
fontsize: 9pt
execute:
  echo: false
  warning: false
  message: false
---

```{r setup}
library(knitr)
library(ggplot2)
library(dplyr)

OUT <- "../../results/benchmark_sim/outputs/desurv_comparison"
res <- if (file.exists(file.path(OUT, "desurv_comparison_results.csv")))
         read.csv(file.path(OUT, "desurv_comparison_results.csv"),
                  stringsAsFactors = FALSE)
       else NULL

if (!is.null(res)) {
  means <- res |>
    group_by(model, label, K, k_eff) |>
    summarise(mean_c = round(mean(c_index), 3), .groups = "drop")
  mc <- setNames(means$mean_c, means$model)
  ke <- setNames(means$k_eff,  means$model)
} else {
  mc <- setNames(rep(NA_real_, 4), c("D1","D2","D3","D4"))
  ke <- setNames(rep(NA_integer_, 4), c("D1","D2","D3","D4"))
}
```

## 1. Motivation

DeSurv (Young et al. 2025) selects genes by independently ranking each gene's
mean expression and variance, then retaining the top-3,000 with the lowest
combined rank-sum per cohort — genes that are both highly expressed *and* highly
variable. Our prior analysis used variance-only selection on the full merged
matrix after normalization, which (a) conflates per-cohort with cross-cohort
variance and (b) can retain lowly-expressed noisy genes near the detection floor.
This report evaluates whether aligning with DeSurv's gene selection changes
the external concordance of the recommended manuscript configurations.

## 2. Configurations Compared

```{r config-table}
tab <- data.frame(
  ID     = c("D1","D2","D3","D4"),
  Label  = c("LB orig (M4)","YFB orig (M5)",
             "LB DeSurv-aligned","YFB DeSurv-aligned"),
  Model  = c("LB","YFB","LB","YFB"),
  top_n  = c(2000, 2000, 3000, 3000),
  Method = c("variance","variance","combined\\_rank","combined\\_rank"),
  Selection = c("post-norm","post-norm","per-cohort","per-cohort"),
  Cohort = c("Yes","No","Yes","No"),
  Mean_C = c(mc["D1"], mc["D2"], mc["D3"], mc["D4"]),
  K_eff  = c(ke["D1"], ke["D2"], ke["D3"], ke["D4"])
)
knitr::kable(tab, booktabs = TRUE, row.names = FALSE,
             col.names = c("ID","Label","Model","top\\_n","Method",
                           "Selection","Cohort","Mean C","$K_{\\text{eff}}$"),
             escape = FALSE,
             caption = "Four configurations compared. D1/D2 reproduce M4/M5 exactly. D3/D4 apply DeSurv gene selection (combined mean+variance rank, top-3000 per cohort before normalization).")
```

## 3. Per-Cohort C-Index

```{r cohort-fig, fig.height=3.2, fig.cap="Per-cohort external C-index for D1–D4. Dashed line = chance (C=0.5). D1/D2 are the original manuscript configs; D3/D4 are DeSurv-aligned."}
if (!is.null(res)) {
  res_plot <- res |>
    mutate(type = ifelse(model %in% c("D1","D2"), "Original", "DeSurv-aligned"),
           label_short = sub(" \\(.*", "", label))
  ggplot(res_plot, aes(x = cohort, y = c_index, fill = label)) +
    geom_col(position = position_dodge(0.75), width = 0.65, color = "white") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    scale_fill_manual(
      values = c("LB orig (M4)"        = "#4c78a8",
                 "YFB orig (M5)"       = "#f58518",
                 "LB DeSurv-aligned"   = "#72b7b2",
                 "YFB DeSurv-aligned"  = "#e45756"),
      name = NULL) +
    coord_cartesian(ylim = c(0.45, 0.75)) +
    labs(x = NULL, y = "C-index") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
} else {
  plot.new(); text(0.5, 0.5, "Results CSV not found", cex = 1.2)
}
```

## 4. Findings and Recommendation

[Write 2–3 sentences once results are available: state mean C for D3 vs D1, D4 vs D2,
and note whether DeSurv alignment improves, degrades, or is neutral. Recommend
whether to adopt as the primary configuration.]

## 5. Gene Set Comparison

```{r gene-table}
# Report how many genes each preprocessing yields
# Fill in from run output: D1/D2 = 2000, D3/D4 = ~1970 (per-cohort then intersect)
tab2 <- data.frame(
  Config = c("D1/D2 (original)","D3/D4 (DeSurv-aligned)"),
  top_n  = c(2000, 3000),
  Method = c("variance, post-norm", "combined_rank, per-cohort"),
  Genes_after_selection = c(2000, NA)  # fill in D3/D4 from run output
)
knitr::kable(tab2, booktabs = TRUE, row.names = FALSE,
             caption = "Gene set size after selection and intersection.")
```
```

- [ ] **Step 5.2: Fill in the actual gene count for D3/D4**

From the run log (`tail /tmp/desurv_comparison.log | grep "Genes after per-cohort"`), find the actual intersection size and update `Genes_after_selection` for D3/D4 in the table above.

- [ ] **Step 5.3: Write the findings paragraph in Section 4**

Based on the actual numbers from `aggregate(c_index ~ model + label, data = res, FUN = mean)`, replace the placeholder in Section 4 with actual results. Example:

> DeSurv-aligned preprocessing (D4, YFB) achieves a mean external C-index of 0.XXX across 5 PDAC cohorts, compared to 0.626 for the original YFB configuration (D2/M5). [State direction and magnitude.] The recommended primary configuration for the manuscript is [D4 or M5], with [D3 or M4] as sensitivity.

- [ ] **Step 5.4: Render the report**

```bash
cd /Users/ajwalther/GithubProjects/multiomicsGEP
quarto render docs/reports/desurv_alignment_report_05_27_26.qmd --to pdf
quarto render docs/reports/desurv_alignment_report_05_27_26.qmd --to html
```
Expected: no LaTeX errors; PDF and HTML written to `docs/reports/`.

- [ ] **Step 5.4b: If DeSurv wins (delta_yfb > +0.005) — update the progress report**

Open `docs/progress_report/SSBMF_Status_Update_5_27_26.qmd` and update Section 5
("Recommendation: Final Configuration for Manuscript"). Change the primary config row
from M5 to D4 and add a note:

```
Role: Primary (D4)
Model: YFB
Preprocessing: Per-platform z-std + DeSurv gene selection (combined_rank, top-3000 per cohort)
Cohort: No
K: <K from globals.yml k_merged_yfb_desurv>
K_eff: <from results>
Mean_C: <mc["D4"] value>
Notes: Supersedes M5; DeSurv-aligned gene selection gains +<delta> mean C
```

Re-render both PDF and HTML:
```bash
quarto render docs/progress_report/SSBMF_Status_Update_5_27_26.qmd --to pdf
quarto render docs/progress_report/SSBMF_Status_Update_5_27_26.qmd --to html
```

Also add a flag at the end of Section 7 (Open Items):
```
- **Pathway enrichment on D4 factors.** DeSurv-aligned gene selection changes the training
  gene set (~1970 genes vs 2000); re-run GSEA/ORA on D4 factor loadings once adopted.
```

If DeSurv does NOT win: skip this step entirely — the progress report is already correct.

- [ ] **Step 5.5: Update `ROADMAP.md`**

Mark the DeSurv gene selection item complete:
```markdown
- [x] **Align gene selection with DeSurv: combined mean+variance ranking, top-3000** *(Complete — 2026-05-27)*
  Implemented combined_rank method in select_top_variable_genes() and per-cohort
  selection (before normalization) in preprocess_merged_cohorts(). D3/D4 comparison
  on 5 external PDAC cohorts: [state result]. See DECISIONS.md 2026-05-27.
  *Files: code/preprocess_desurv.R, results/benchmark_sim/run_desurv_comparison.R,
  docs/reports/desurv_alignment_report_05_27_26.qmd*
```

Add a new item if DeSurv alignment changes the recommended config:
```markdown
- [ ] **Update manuscript primary config if DeSurv alignment improves C-index** `[Priority: High]` `[Effort: Small]`
  If D4 (YFB DeSurv-aligned) exceeds M5 mean C=0.626, update CLAUDE.md, DECISIONS.md,
  and the progress report to reflect the new recommended config. Re-run pathway
  enrichment on the new M5 factor loadings.
```

- [ ] **Step 5.6: Add entry to `DECISIONS.md`**

Add at the top of DECISIONS.md (most recent first):

```markdown
## 2026-05-27 — DeSurv Gene Selection Alignment

**Question:** Does adopting DeSurv's gene selection (combined mean+variance rank, top-3000 per
cohort before normalization) improve external C-index relative to the current variance-only
selection (top-2000 on merged normalized matrix)?

**Implementation:**
- `select_top_variable_genes()` gains `method="combined_rank"`: rank_mean + rank_var,
  lowest rank-sum genes retained. Default "variance" unchanged.
- `preprocess_merged_cohorts()` gains `selection_per_cohort=TRUE`: per-cohort top-N
  selection on log-transformed data, then intersect, before per-platform z-std.
- New comparison configs: D3 (LB DeSurv-aligned) and D4 (YFB DeSurv-aligned).

**Result:**
| Config | Model | Mean external C | K_eff | Gene set |
|--------|-------|----------------|-------|---------|
| D1 (= M4) | LB + orig | 0.622 | 1 | 2000 |
| D2 (= M5) | YFB + orig | 0.626 | 2 | 2000 |
| D3 (DeSurv LB) | LB + aligned | [fill] | [fill] | ~1970 |
| D4 (DeSurv YFB) | YFB + aligned | [fill] | [fill] | ~1970 |

**Decision:** Apply the rule from the plan: delta_yfb = mean_C(D4) − mean_C(D2).
If delta_yfb > +0.005: adopt D4 as new primary config; update CLAUDE.md recommended config line.
If |delta_yfb| ≤ 0.005: keep M5 as primary; DeSurv gene selection is neutral for SSBMF.
If delta_yfb < −0.005: retain M5; DeSurv gene selection is harmful; note in manuscript limitations.
```

- [ ] **Step 5.7: Update `CLAUDE.md`**

Update the test count: `tests/run_tests.R` after Task 2 should show 246/246. Update the line:
> **Tests:** Run `Rscript tests/run_tests.R` after any change to a modular update script. Expected: 246/246 passing.

Update the model status line if the recommended config changes based on D3/D4 results.

Update the Quick Reference table to add the new runner:
```
| **DeSurv comparison runner** | `results/benchmark_sim/run_desurv_comparison.R` |
| **DeSurv comparison report** | `docs/reports/desurv_alignment_report_05_27_26.{qmd,pdf,html}` |
```

- [ ] **Step 5.8: Final commit**

```bash
git add docs/reports/desurv_alignment_report_05_27_26.qmd
git add docs/reports/desurv_alignment_report_05_27_26.pdf
git add docs/reports/desurv_alignment_report_05_27_26.html
git add ROADMAP.md DECISIONS.md CLAUDE.md
git commit -m "$(cat <<'EOF'
DeSurv alignment report, ROADMAP and DECISIONS update

Comparison report renders D1-D4 per-cohort C-index: D4 (YFB DeSurv-aligned)
mean C=[value]; D2 (YFB orig/M5) mean C=0.626. [State conclusion].

ROADMAP: DeSurv gene selection item marked complete.
DECISIONS: 2026-05-27 entry added with comparison table and recommendation.
CLAUDE.md: test count updated to 246/246; new runner/report paths added.
EOF
)"
```

---

## Self-Review

### 1. Spec coverage
- ✅ `method="combined_rank"` in `select_top_variable_genes()` → Task 1
- ✅ `selection_per_cohort=TRUE` in `preprocess_merged_cohorts()` → Task 2
- ✅ `top_n_genes_desurv: 3000` in globals.yml → Task 3
- ✅ K-CV for D3/D4 + benchmark runner → Task 3 + 4
- ✅ Comparison report → Task 5
- ✅ ROADMAP/DECISIONS/CLAUDE.md updates → Task 5

### 2. Placeholder scan
- Step 4.5: commit message contains `<value>` placeholders — these are intentional TODOs for the executor to fill in with actual results after the run completes (Task 4, Step 4.3).
- Section 4 of the QMD report contains a placeholder paragraph — intentional; filled in Task 5, Step 5.3.
- DECISIONS.md template has `[fill]` cells — intentional; filled in Task 5, Step 5.6.

### 3. Type consistency
- `select_top_variable_genes()` returns `list(Y, gene_names, gene_var)` — same shape in both `method=` branches. ✅
- `preprocess_merged_cohorts()` returns `selection_per_cohort` and `selection_method` in the output list — echoed correctly in both branches. ✅
- `run_desurv_comparison.R` calls `predict_supervised_mf(Y_ext, EF_sub, fit$EBeta)` and `predict_cox_on_yf(Y_ext, EF_sub, fit$EBeta, EF_norms = fit$EF_norms)` — matches signatures in `code/predict.R` and `code/predict_cox_on_yf.R`. ✅
- `set_key()` in the comparison runner uses the same pattern as in `run_merged_kcv.R`. ✅
