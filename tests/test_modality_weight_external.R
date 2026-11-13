# ============================================================
# Script: test_modality_weight_external.R
# Purpose: Check frozen external feature preparation and program direction matching.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R and existing test helpers
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/multimodal_real/score_modality_weight_external.R")
run_test("Weighted external: series z-scores precede ordered alignment and zero filling", {
  Y <- cbind(a = 1:3, b = rep(7, 3), unused = c(5, 8, 2))
  x <- mwext_align(Y, c("b", "missing", "a"))
  assert_equal(colnames(x$Y), c("b", "missing", "a"))
  assert_near(x$Y[, 1:2], matrix(0, 3, 2), 0)
  assert_near(x$Y[, 3], c(-1, 0, 1), 1e-12)
  assert_equal(x$missing, 1L); assert_equal(x$constant, 1L)
  x2 <- mwext_align(Y * 4, c("b", "missing", "a"))
  assert_near(x2$Y, x$Y, 1e-12)
})
run_test("Weighted external: program matches use TCGA and direction never flips", {
  E <- cbind(c(-2, -1, 1, 2), c(1, -1, -1, 1))
  J <- cbind(-E[, 2], E[, 1])
  X <- E; XJ <- cbind(-E[, 2], -E[, 1])
  a <- mwext_direction(J, E, XJ, X)
  assert_equal(a$expression_program, c(2L, 1L))
  assert_true(a$kept_direction[1]); assert_false(a$kept_direction[2])
  assert_near(a$tcga_correlation, c(-1, 1), 1e-12)
  assert_near(a$icgc_correlation, c(-1, -1), 1e-12)
})
if (sys.nframe() == 0L) {
  result <- report_results("Weighted external")
  if (result$failed > 0L) quit(status = 1L)
}
