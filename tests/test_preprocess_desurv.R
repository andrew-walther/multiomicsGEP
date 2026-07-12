# =============================================================================
# tests/test_preprocess_desurv.R
# Unit tests for DeSurv-aligned preprocessing helpers.
# =============================================================================

cat("\n========================================\n")
cat("  Tests: preprocess_desurv.R\n")
cat("========================================\n\n")

run_test("T1.1: log2_plus1_transform matches base formula", {
  Y <- matrix(c(0, 1, 3, 7), nrow = 2)
  out <- log2_plus1_transform(Y)
  assert_near(out, log2(Y + 1), tol = 1e-12)
})

run_test("T1.2: select_top_variable_genes keeps highest-variance genes", {
  Y <- cbind(
    g1 = c(1, 1, 1, 1),
    g2 = c(1, 5, 1, 5),
    g3 = c(2, 9, 2, 9),
    g4 = c(4, 4, 4, 4)
  )
  out <- select_top_variable_genes(Y, colnames(Y), top_n = 2)
  assert_equal(out$gene_names, c("g3", "g2"))
  assert_equal(ncol(out$Y), 2L)
})

run_test("T1.3: rank_transform_subjects ranks within each subject", {
  Y <- rbind(
    c(10, 20, 30),
    c(7, 1, 4)
  )
  colnames(Y) <- c("a", "b", "c")
  out <- rank_transform_subjects(Y)
  expected <- rbind(c(1, 2, 3), c(3, 1, 2))
  colnames(expected) <- colnames(Y)
  assert_equal(unname(out), unname(expected))
})

run_test("T1.4: rank_transform_subjects uses average rank for ties", {
  Y <- matrix(c(5, 5, 9), nrow = 1)
  out <- rank_transform_subjects(Y)
  assert_near(as.numeric(out[1, ]), c(1.5, 1.5, 3), tol = 1e-12)
})

run_test("T1.5: preprocess_desurv_cohort applies log-filter-rank pipeline", {
  Y <- rbind(
    c(0, 10, 30, 5),
    c(0, 20, 5, 15),
    c(0, 30, 10, 25)
  )
  gene_names <- c("g1", "g2", "g3", "g4")
  out <- preprocess_desurv_cohort(Y, gene_names, top_n = 3, cohort_name = "toy")

  assert_equal(out$cohort_name, "toy")
  assert_equal(out$n, 3L)
  assert_equal(out$p, 3L)
  assert_equal(length(out$gene_names), 3L)
  assert_finite(out$Y)
  assert_true(all(apply(out$Y, 1, function(x) setequal(as.numeric(x), 1:3))),
              "Each subject should have ranks 1..p after preprocessing")
})

run_test("T1.5b: preprocess_desurv_cohort defaults to rank_transform=TRUE (backward compatible)", {
  Y <- rbind(
    c(0, 10, 30, 5),
    c(0, 20, 5, 15),
    c(0, 30, 10, 25)
  )
  gene_names <- c("g1", "g2", "g3", "g4")
  out <- preprocess_desurv_cohort(Y, gene_names, top_n = 3, cohort_name = "toy")
  assert_true(all(apply(out$Y, 1, function(x) setequal(as.numeric(x), 1:3))),
              "Default behavior must still rank-transform (backward compatible)")
})

run_test("Phase1c-T1: rank_transform=FALSE, per_platform_standardize=TRUE gives per-gene mean~0, sd~1 (not ranks)", {
  # Matches the training pipeline's per_platform_standardize step
  # (code/preprocess_desurv.R's preprocess_merged_cohorts with
  # per_platform_standardize=TRUE, rank_transform=FALSE) -- external cohort
  # preprocessing must use the SAME transform, not the default rank transform.
  set.seed(11)
  Y <- matrix(rnorm(20 * 5, mean = 10, sd = 3), 20, 5)
  gene_names <- paste0("g", 1:5)
  out <- preprocess_desurv_cohort(Y, gene_names, top_n = NULL,
                                  rank_transform = FALSE,
                                  per_platform_standardize = TRUE,
                                  cohort_name = "ext")
  assert_near(colMeans(out$Y), rep(0, 5), tol = 1e-8,
              msg = "per_platform_standardize should give colMean ~ 0")
  assert_near(apply(out$Y, 2, sd), rep(1, 5), tol = 1e-8,
              msg = "per_platform_standardize should give colSD ~ 1")
  # Explicitly NOT rank-transformed: values should not be a permutation of 1:p
  assert_true(!all(apply(out$Y, 1, function(x) setequal(round(x, 6), round(rank(x), 6)))),
              "Should not be rank-transformed when rank_transform=FALSE")
})

run_test("Phase1c-T2: rank_transform=FALSE with per_platform_standardize=FALSE skips both transforms (raw log-scale passthrough)", {
  set.seed(12)
  Y <- matrix(rnorm(10 * 4, mean = 5, sd = 2), 10, 4)
  gene_names <- paste0("g", 1:4)
  out <- preprocess_desurv_cohort(Y, gene_names, top_n = NULL,
                                  log_transform = FALSE,
                                  rank_transform = FALSE,
                                  per_platform_standardize = FALSE,
                                  cohort_name = "raw")
  assert_near(out$Y, Y, tol = 1e-10,
              msg = "With both transforms off and no log transform, Y should pass through unchanged")
})

run_test("T1.6: intersect_preprocessed_cohorts subsets to common genes in reference order", {
  c1 <- list(Y = matrix(1:6, nrow = 2), gene_names = c("g1", "g2", "g3"), p = 3L)
  c2 <- list(Y = matrix(7:12, nrow = 2), gene_names = c("g3", "g2", "g4"), p = 3L)
  out <- intersect_preprocessed_cohorts(list(c1, c2))

  assert_equal(out[[1]]$gene_names, c("g2", "g3"))
  assert_equal(out[[2]]$gene_names, c("g2", "g3"))
  assert_equal(ncol(out[[1]]$Y), 2L)
  assert_equal(ncol(out[[2]]$Y), 2L)
})

run_test("T1.7: merge_preprocessed_cohorts row-binds aligned cohorts", {
  c1 <- list(Y = matrix(1:4, nrow = 2), gene_names = c("g1", "g2"))
  c2 <- list(Y = matrix(5:10, nrow = 3), gene_names = c("g1", "g2"))
  out <- merge_preprocessed_cohorts(list(TCGA = c1, CPTAC = c2))

  assert_equal(out$n, 5L)
  assert_equal(out$p, 2L)
  assert_equal(levels(out$dataset_labels), c("TCGA", "CPTAC"))
  assert_equal(as.integer(table(out$dataset_labels)), c(2L, 3L))
})

run_test("T1.8: merge_preprocessed_cohorts errors on mismatched gene ordering", {
  c1 <- list(Y = matrix(1:4, nrow = 2), gene_names = c("g1", "g2"))
  c2 <- list(Y = matrix(5:8, nrow = 2), gene_names = c("g2", "g1"))
  err <- tryCatch(merge_preprocessed_cohorts(list(c1, c2)), error = function(e) e)
  assert_true(inherits(err, "error"), "Expected merge_preprocessed_cohorts to error")
})

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

run_test("T1.9: normalize_method='quantile' (default) matches prior behavior", {
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

run_test("T1.10: normalize_method='z_score' gives column means ~0 and SDs ~1", {
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

run_test("T1.11: normalize_method='none' preserves log-transformed values (no QN distortion)", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  out_none <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                        normalize_method = "none",
                                        rank_transform   = FALSE)
  out_q    <- preprocess_merged_cohorts(raw, flags, top_n = 6L,
                                        normalize_method = "quantile",
                                        rank_transform   = FALSE)
  # 'none' and 'quantile' must produce DIFFERENT matrices
  assert_equal(identical(out_none$Y, out_q$Y), FALSE)
  # 'none' matrix gene names and dimensions must still be correct
  assert_equal(ncol(out_none$Y), 6L)
  assert_equal(nrow(out_none$Y), 8L)
  assert_equal(out_none$gene_names, out_q$gene_names)
})

run_test("T1.12: normalize_method invalid argument is caught", {
  raw   <- make_synthetic_raw_list()
  flags <- c(CohortA = TRUE, CohortB = FALSE)
  result <- tryCatch(
    preprocess_merged_cohorts(raw, flags, top_n = 6L,
                              normalize_method = "bad_method"),
    error = function(e) e
  )
  assert_equal(inherits(result, "error"), TRUE)
})

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
  Y <- cbind(
    g_hm_hv = c(100, 200, 300, 400),
    g_lm_hv = c(  1,   2, 100, 200),
    g_hm_lv = c(100, 101, 100, 101)
  )
  out <- select_top_variable_genes(Y, colnames(Y), top_n = 1L, method = "combined_rank")
  assert_equal(out$gene_names, "g_hm_hv")
  assert_equal(ncol(out$Y), 1L)
})

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
  # Fixture uses set.seed(1L): CohortA top-4 = {G1,G2,G4,G5}, CohortB = {G2,G3,G5,G6}.
  # Intersection = {G2, G5} → p=2 (deterministic).
  assert_equal(out$p, 2L)
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

report_results("preprocess_desurv.R")
