# ============================================================
# Script: test_modality_exclusive_generator.R
# Purpose: Verify opt-in modality-exclusive simulated programs.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R; survival; code/simulate_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")

run_test("Exclusive generator: default and FALSE preserve all output", {
  args <- list(n_train = 20L, n_validation = 12L, p_expression = 10L,
               p_methylation = 15L, seed = 71L)
  original <- do.call(simulate_multimodal_yfb_data, args)
  disabled <- do.call(simulate_multimodal_yfb_data, c(args, list(modality_exclusive_programs = FALSE)))
  assert_true(identical(original, disabled))
})

run_test("Exclusive generator: planted blocks and cohort dimensions match the design", {
  args <- list(n_train = 20L, n_validation = 12L, p_expression = 10L,
               p_methylation = 15L, seed = 71L)
  original <- do.call(simulate_multimodal_yfb_data, args)
  d <- do.call(simulate_multimodal_yfb_data, c(args, list(modality_exclusive_programs = TRUE)))
  e <- d$truth$F$expression
  m <- d$truth$F$methylation
  assert_true(all(m[, 1] == 0) && all(e[, 2] == 0))
  assert_true(any(e[, 1] != 0) && any(m[, 2] != 0))
  assert_true(any(e[, 3] != 0) && any(m[, 3] != 0))
  assert_true(identical(e[, c(1, 3)], original$truth$F$expression[, c(1, 3)]))
  assert_true(identical(m[, 2:3], original$truth$F$methylation[, 2:3]))
  assert_true(identical(d$truth$beta, c(0.70, 0.35, 0)))
  assert_true(identical(dim(d$training$Y$expression), c(20L, 10L)))
  assert_true(identical(dim(d$training$Y$methylation), c(20L, 15L)))
  assert_true(identical(dim(d$validation$Y$expression), c(12L, 10L)))
  assert_true(identical(dim(d$validation$Y$methylation), c(12L, 15L)))
})

run_test("Exclusive generator: invalid option and incompatible scenario fail", {
  for (value in list(NULL, NA, 1, "TRUE", c(TRUE, FALSE))) {
    error <- tryCatch({simulate_multimodal_yfb_data(modality_exclusive_programs = value); FALSE},
                      error = function(e) TRUE)
    assert_true(error)
  }
  error <- tryCatch({simulate_multimodal_yfb_data("expression_only",
    modality_exclusive_programs = TRUE); FALSE}, error = function(e) TRUE)
  assert_true(error)
})

if (sys.nframe() == 0L) {
  result <- report_results("test_modality_exclusive_generator.R")
  if (result$failed > 0L) quit(status = 1L)
}
