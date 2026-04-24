# ============================================================
# tests/test_select_alpha_cv.R
# Tests for select_alpha_cv() and its fold-based alpha selection behavior.
#
# Usage: sourced by tests/run_tests.R
# ============================================================

library(survival)

source("code/select_alpha_cv.R")
source("code/train_test_split.R")
source("code/predict.R")

cat("\n========================================\n")
cat("  Tests: select_alpha_cv.R\n")
cat("========================================\n\n")

# Small synthetic dataset for CV tests
.alpha_cv_test_data <- local({
  set.seed(321)
  n <- 60; p <- 40; K_true <- 2
  L_true <- matrix(rexp(n * K_true, rate = 1), n, K_true)
  F_true <- matrix(rexp(p * K_true, rate = 1), p, K_true)
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = 0.4), n, p)
  lp <- 0.8 * L_true[, 1] - 0.5 * L_true[, 2]
  time <- rexp(n, rate = exp(scale(lp)))
  status <- rep(c(1, 0), length.out = n)
  list(Y = Y, time = time, status = status)
})

cat("=== T1: Output Shape and Basic Behavior ===\n")

run_test("T1.1: select_alpha_cv returns alpha_opt, cv_table, fold_results", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 42,
                         max_iter = 8)
  assert_true(!is.null(res$alpha_opt), msg = "alpha_opt missing")
  assert_true(!is.null(res$cv_table), msg = "cv_table missing")
  assert_true(!is.null(res$fold_results), msg = "fold_results missing")
})

run_test("T1.2: cv_table has one row per alpha and expected columns", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 42,
                         max_iter = 8)
  assert_true(nrow(res$cv_table) == 3, msg = "Expected one cv_table row per alpha")
  assert_true(all(c("alpha", "mean_cindex", "se_cindex", "n_folds") %in% names(res$cv_table)),
              msg = "cv_table missing expected columns")
})

run_test("T1.3: fold_results has one row per alpha-fold combination", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 42,
                         max_iter = 8)
  assert_true(nrow(res$fold_results) == 9,
              msg = sprintf("Expected 9 alpha-fold rows, got %d", nrow(res$fold_results)))
})

run_test("T1.4: selected alpha is always one of the supplied grid values", {
  d <- .alpha_cv_test_data
  grid <- c(0.1, 0.5, 0.9)
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = grid,
                         n_folds = 3, K_max = 2, seed = 42,
                         max_iter = 8)
  assert_true(res$alpha_opt %in% grid, msg = "alpha_opt should be chosen from alpha_grid")
})

cat("\n=== T2: Reproducibility and Selection Rule ===\n")

run_test("T2.1: Same seed gives identical CV summaries", {
  d <- .alpha_cv_test_data
  res1 <- select_alpha_cv(d$Y, d$time, d$status,
                          alpha_grid = c(0.1, 0.5, 0.9),
                          n_folds = 3, K_max = 2, seed = 11,
                          max_iter = 8)
  res2 <- select_alpha_cv(d$Y, d$time, d$status,
                          alpha_grid = c(0.1, 0.5, 0.9),
                          n_folds = 3, K_max = 2, seed = 11,
                          max_iter = 8)
  assert_true(identical(res1$cv_table, res2$cv_table),
              msg = "CV tables should match for the same seed")
  assert_true(identical(res1$fold_results, res2$fold_results),
              msg = "Fold-level results should match for the same seed")
})

run_test("T2.2: use_1se=FALSE reports selection_rule='max'", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 21,
                         use_1se = FALSE, max_iter = 8)
  assert_equal(res$selection_rule, "max", msg = "Expected max selection rule")
})

run_test("T2.3: use_1se=TRUE reports selection_rule='1se'", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 21,
                         use_1se = TRUE, max_iter = 8)
  assert_equal(res$selection_rule, "1se", msg = "Expected 1se selection rule")
})

run_test("T2.4: 1-SE rule chooses an alpha whose mean exceeds the threshold", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 21,
                         use_1se = TRUE, max_iter = 8)
  best_idx <- which.max(res$cv_table$mean_cindex)
  threshold <- res$cv_table$mean_cindex[best_idx] - res$cv_table$se_cindex[best_idx]
  chosen_mean <- res$cv_table$mean_cindex[res$cv_table$alpha == res$alpha_opt]
  assert_true(chosen_mean >= threshold,
              msg = "Chosen alpha should satisfy the 1-SE threshold")
})

cat("\n=== T3: Input Validation ===\n")

run_test("T3.1: alpha outside [0,1] errors", {
  d <- .alpha_cv_test_data
  err <- tryCatch(
    select_alpha_cv(d$Y, d$time, d$status,
                    alpha_grid = c(-0.1, 0.5), n_folds = 3, K_max = 2,
                    max_iter = 4),
    error = function(e) conditionMessage(e)
  )
  assert_true(grepl("\\[0, 1\\]", err), msg = "Should reject alpha outside [0,1]")
})

run_test("T3.2: too many folds for event count errors", {
  d <- .alpha_cv_test_data
  status_bad <- c(rep(1, 2), rep(0, nrow(d$Y) - 2))
  err <- tryCatch(
    select_alpha_cv(d$Y, d$time, status_bad,
                    alpha_grid = c(0.1, 0.5), n_folds = 3, K_max = 2,
                    max_iter = 4),
    error = function(e) conditionMessage(e)
  )
  assert_true(grepl("at least 3 events", err), msg = "Should reject insufficient event counts")
})

run_test("T3.3: C-index values are finite and bounded", {
  d <- .alpha_cv_test_data
  res <- select_alpha_cv(d$Y, d$time, d$status,
                         alpha_grid = c(0.1, 0.5, 0.9),
                         n_folds = 3, K_max = 2, seed = 42,
                         max_iter = 8)
  assert_finite(res$fold_results$cindex)
  assert_true(all(res$fold_results$cindex >= 0 & res$fold_results$cindex <= 1),
              msg = "C-index should lie in [0, 1]")
})

report_results("select_alpha_cv.R")
