# ============================================================
# tests/test_predict.R
# Tests for predict_supervised_mf() and stratified_split()
#
# Usage: sourced by tests/run_tests.R (or standalone)
# ============================================================

library(survival)

# Source modules under test
source("code/predict.R")
source("code/train_test_split.R")

cat("\n========================================\n")
cat("  Tests: predict.R + train_test_split.R\n")
cat("========================================\n\n")

# ==============================================================================
# T1: predict_supervised_mf — Output Dimensions ----
# ==============================================================================

cat("=== T1: Output Dimensions ===\n")

run_test("T1.1: L_test has correct dimensions (n_test x K)", {
  n_test <- 20; p <- 100; K <- 5
  Y_test <- matrix(rnorm(n_test * p), n_test, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- rnorm(K)
  pred   <- predict_supervised_mf(Y_test, EF, EBeta)
  # Use == not assert_equal: nrow() returns integer, n_test is double
  assert_true(nrow(pred$L_test) == n_test, msg = "nrow(L_test) != n_test")
  assert_true(ncol(pred$L_test) == K,      msg = "ncol(L_test) != K")
})

run_test("T1.2: risk_scores has correct length (n_test)", {
  n_test <- 15; p <- 50; K <- 3
  Y_test <- matrix(rnorm(n_test * p), n_test, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- rnorm(K)
  pred   <- predict_supervised_mf(Y_test, EF, EBeta)
  assert_length(pred$risk_scores, n_test)
})

run_test("T1.3: Single patient prediction works (n_test = 1)", {
  p <- 50; K <- 3
  Y_test <- matrix(rnorm(p), 1, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- rnorm(K)
  pred   <- predict_supervised_mf(Y_test, EF, EBeta)
  assert_true(nrow(pred$L_test) == 1, msg = "nrow(L_test) != 1 for single patient")
  assert_length(pred$risk_scores, 1)
})

# ==============================================================================
# T2: predict_supervised_mf — Consistency Checks ----
# ==============================================================================

cat("\n=== T2: Consistency Checks ===\n")

run_test("T2.1: Training data projection approximately recovers in-sample risk scores", {
  # If we project Y_train through EF, the resulting L should be close to
  # the original EL (not exact due to regularisation and EBNM posterior ≠ MLE)
  set.seed(42)
  n <- 50; p <- 200; K <- 3
  EL_true <- matrix(rnorm(n * K), n, K)
  EF_true <- matrix(rnorm(p * K), p, K)
  Y <- EL_true %*% t(EF_true) + matrix(rnorm(n * p, sd = 0.5), n, p)
  EBeta <- c(1.0, -0.5, 0.3)
  pred <- predict_supervised_mf(Y, EF_true, EBeta)
  # Check that L_test is correlated with EL_true
  cor_diag <- diag(cor(EL_true, pred$L_test))
  assert_true(all(abs(cor_diag) > 0.5),
    msg = "Projected loadings should correlate with true loadings")
})

run_test("T2.2: risk_scores = L_test %*% EBeta (exact identity)", {
  n_test <- 10; p <- 30; K <- 4
  Y_test <- matrix(rnorm(n_test * p), n_test, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- c(1, -2, 0.5, 0)
  pred   <- predict_supervised_mf(Y_test, EF, EBeta)
  expected <- as.vector(pred$L_test %*% EBeta)
  assert_near(pred$risk_scores, expected, tol = 1e-10)
})

run_test("T2.3: All-zero EBeta produces all-zero risk scores", {
  n_test <- 10; p <- 30; K <- 3
  Y_test <- matrix(rnorm(n_test * p), n_test, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- rep(0, K)
  pred   <- predict_supervised_mf(Y_test, EF, EBeta)
  assert_near(pred$risk_scores, rep(0, n_test), tol = 1e-10)
})

# ==============================================================================
# T3: predict_supervised_mf — Error Handling ----
# ==============================================================================

cat("\n=== T3: Error Handling ===\n")

run_test("T3.1: Dimension mismatch Y_test cols != EF rows errors", {
  Y_test <- matrix(rnorm(10 * 50), 10, 50)
  EF     <- matrix(rnorm(60 * 3), 60, 3)  # 60 != 50
  EBeta  <- rnorm(3)
  err <- tryCatch(predict_supervised_mf(Y_test, EF, EBeta),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("mismatch", err, ignore.case = TRUE),
    msg = "Should error on p mismatch")
})

run_test("T3.2: EBeta length != K errors", {
  p <- 50; K <- 3
  Y_test <- matrix(rnorm(10 * p), 10, p)
  EF     <- matrix(rnorm(p * K), p, K)
  EBeta  <- rnorm(K + 1)  # length 4, not 3
  err <- tryCatch(predict_supervised_mf(Y_test, EF, EBeta),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("mismatch", err, ignore.case = TRUE),
    msg = "Should error on K mismatch")
})

run_test("T3.3: Non-matrix Y_test errors gracefully", {
  EF <- matrix(rnorm(30), 10, 3)
  EBeta <- rnorm(3)
  err <- tryCatch(predict_supervised_mf(c(1,2,3), EF, EBeta),
                  error = function(e) conditionMessage(e))
  assert_true(is.character(err), msg = "Should error on non-matrix input")
})

# ==============================================================================
# T4: stratified_split — Basic Properties ----
# ==============================================================================

cat("\n=== T4: Stratified Split Properties ===\n")

run_test("T4.1: Train and test indices are disjoint and cover all samples", {
  status <- c(rep(1, 30), rep(0, 70))
  sp <- stratified_split(status, test_frac = 0.2, seed = 42)
  all_idx <- sort(c(sp$train_idx, sp$test_idx))
  n <- length(status)
  # Use == comparisons — integer vs double type from seq_len vs 1:n
  assert_true(length(all_idx) == n, msg = "Combined indices should cover all n samples")
  assert_true(all(all_idx == seq_len(n)), msg = "Indices should span 1..n with no gaps")
  assert_true(length(intersect(sp$train_idx, sp$test_idx)) == 0,
              msg = "Train and test sets must be disjoint")
})

run_test("T4.2: Test set size is approximately test_frac of total", {
  status <- c(rep(1, 50), rep(0, 50))
  sp <- stratified_split(status, test_frac = 0.2, seed = 42)
  # Allow ±5% tolerance
  assert_true(sp$n_test >= 15 && sp$n_test <= 25,
    msg = sprintf("Expected ~20 test samples, got %d", sp$n_test))
})

run_test("T4.3: Event rate preserved in both sets (within 15%)", {
  status <- c(rep(1, 40), rep(0, 60))
  sp <- stratified_split(status, test_frac = 0.2, seed = 42)
  overall_rate <- mean(status == 1)
  assert_true(abs(sp$event_rate_train - overall_rate) < 0.15,
    msg = "Training event rate should be close to overall rate")
  assert_true(abs(sp$event_rate_test - overall_rate) < 0.15,
    msg = "Test event rate should be close to overall rate")
})

run_test("T4.4: At least 1 event and 1 censored in test set", {
  status <- c(rep(1, 10), rep(0, 90))
  sp <- stratified_split(status, test_frac = 0.2, seed = 42)
  assert_true(sum(status[sp$test_idx] == 1) >= 1,
    msg = "Test set must contain at least 1 event")
  assert_true(sum(status[sp$test_idx] == 0) >= 1,
    msg = "Test set must contain at least 1 censored")
})

run_test("T4.5: Reproducible with same seed", {
  status <- c(rep(1, 50), rep(0, 50))
  sp1 <- stratified_split(status, test_frac = 0.2, seed = 99)
  sp2 <- stratified_split(status, test_frac = 0.2, seed = 99)
  assert_true(identical(sp1$train_idx, sp2$train_idx),
    msg = "Same seed should give same split")
})

run_test("T4.6: Different seed gives different split", {
  status <- c(rep(1, 50), rep(0, 50))
  sp1 <- stratified_split(status, test_frac = 0.2, seed = 1)
  sp2 <- stratified_split(status, test_frac = 0.2, seed = 2)
  assert_true(!identical(sp1$test_idx, sp2$test_idx),
    msg = "Different seeds should give different splits")
})

# ==============================================================================
# T5: stratified_split — Edge Cases ----
# ==============================================================================

cat("\n=== T5: Edge Cases ===\n")

run_test("T5.1: Small dataset (n=10) works", {
  status <- c(1, 1, 1, 0, 0, 0, 0, 1, 0, 1)
  sp <- stratified_split(status, test_frac = 0.2, seed = 42)
  assert_true(sp$n_test >= 2 && sp$n_test <= 4,
    msg = sprintf("Expected 2-4 test samples from n=10, got %d", sp$n_test))
})

run_test("T5.2: Very small dataset (n=4) still works", {
  status <- c(1, 1, 0, 0)
  sp <- stratified_split(status, test_frac = 0.3, seed = 42)
  assert_true(sp$n_test >= 2,
    msg = "Should have at least 2 test samples (1 event + 1 censored)")
  assert_true(sp$n_train >= 2,
    msg = "Should have at least 2 training samples")
})

run_test("T5.3: Too few events (n_event=1) errors", {
  status <- c(1, 0, 0, 0, 0)
  err <- tryCatch(stratified_split(status, test_frac = 0.2, seed = 42),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("at least 2 events", err),
    msg = "Should error with fewer than 2 events")
})

run_test("T5.4: n < 4 errors", {
  status <- c(1, 0, 0)
  err <- tryCatch(stratified_split(status, test_frac = 0.2, seed = 42),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("at least 4", err),
    msg = "Should error with n < 4")
})

report_results("predict.R + train_test_split.R")
