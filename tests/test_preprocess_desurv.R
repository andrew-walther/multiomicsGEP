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

report_results("preprocess_desurv.R")
