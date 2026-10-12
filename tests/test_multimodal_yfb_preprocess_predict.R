# ============================================================
# Script: test_multimodal_yfb_preprocess_predict.R
# Purpose: Test data-contract and frozen feature-alignment helpers for multimodal YFB.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("preprocess_multimodal_yfb_training")) {
  source("code/preprocess_multimodal_yfb.R")
}

cat("\n========================================\n")
cat("  Tests: preprocess_multimodal_yfb.R\n")
cat("========================================\n\n")

multimodal_yfb_fixture <- function() {
  list(
    expression = structure(matrix(c(1, 2, 3, 4), nrow = 2),
                           dimnames = list(c("s2", "s1"), c("g1", "g2"))),
    methylation = structure(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2),
                             dimnames = list(c("s1", "s2"), c("cg1", "cg2")))
  )
}

# T1: Training contract ----

run_test("MMYFB-Preprocess-T1: matched IDs are reordered to expression order", {
  prepared <- preprocess_multimodal_yfb_training(
    Y = multimodal_yfb_fixture(),
    time = c(s1 = 2, s2 = 1), event = c(s1 = 0, s2 = 1)
  )
  assert_equal(rownames(prepared$Y$expression), c("s2", "s1"))
  assert_equal(rownames(prepared$Y$methylation), c("s2", "s1"))
  assert_near(prepared$time, c(1, 2), tol = 1e-12)
  assert_equal(prepared$training_spec$feature_names$expression, c("g1", "g2"))
})

run_test("MMYFB-Preprocess-T2: unmatched subject IDs error instead of intersecting", {
  Y <- multimodal_yfb_fixture()
  rownames(Y$methylation) <- c("s1", "other")
  error <- tryCatch({
    preprocess_multimodal_yfb_training(Y, c(s1 = 1, s2 = 2), c(s1 = 1, s2 = 0))
    NULL
  }, error = conditionMessage)
  assert_true(grepl("subject", error, ignore.case = TRUE))
})

run_test("MMYFB-Preprocess-T3: duplicate feature names and negative data error", {
  duplicate_Y <- multimodal_yfb_fixture()
  colnames(duplicate_Y$expression) <- c("g1", "g1")
  duplicate_error <- tryCatch({
    preprocess_multimodal_yfb_training(duplicate_Y, c(s1 = 1, s2 = 2), c(s1 = 1, s2 = 0))
    NULL
  }, error = conditionMessage)
  assert_true(grepl("unique", duplicate_error, ignore.case = TRUE))

  negative_Y <- multimodal_yfb_fixture()
  negative_Y$expression[1, 1] <- -1
  negative_error <- tryCatch({
    preprocess_multimodal_yfb_training(negative_Y, c(s1 = 1, s2 = 2), c(s1 = 1, s2 = 0))
    NULL
  }, error = conditionMessage)
  assert_true(grepl("nonnegative", negative_error, ignore.case = TRUE))
})

# T2: Frozen prediction alignment ----

run_test("MMYFB-Preprocess-T4: validation features are reordered and extras reported", {
  prepared <- preprocess_multimodal_yfb_training(
    Y = multimodal_yfb_fixture(),
    time = c(s1 = 2, s2 = 1), event = c(s1 = 0, s2 = 1)
  )
  validation <- list(
    expression = structure(matrix(c(30, 10, 20), nrow = 1),
                           dimnames = list("v1", c("g2", "g1", "extra"))),
    methylation = structure(matrix(c(40, 50, 60), nrow = 1),
                             dimnames = list("v1", c("cg2", "extra", "cg1")))
  )
  aligned <- align_multimodal_yfb_prediction(validation, prepared$training_spec)
  assert_near(aligned$Y$expression, matrix(c(10, 30), nrow = 1), tol = 1e-12)
  assert_near(aligned$Y$methylation, matrix(c(60, 40), nrow = 1), tol = 1e-12)
  assert_equal(aligned$ignored_features$expression, "extra")
  assert_equal(aligned$ignored_features$methylation, "extra")
})

run_test("MMYFB-Preprocess-T5: missing training feature and nonfinite validation data error", {
  spec <- preprocess_multimodal_yfb_training(
    Y = multimodal_yfb_fixture(), time = c(s1 = 2, s2 = 1), event = c(s1 = 0, s2 = 1)
  )$training_spec
  missing <- list(expression = matrix(1, nrow = 1, dimnames = list("v1", "g1")),
                  methylation = matrix(c(1, 2), nrow = 1,
                                       dimnames = list("v1", c("cg1", "cg2"))))
  missing_error <- tryCatch({ align_multimodal_yfb_prediction(missing, spec); NULL },
                            error = conditionMessage)
  assert_true(grepl("missing", missing_error, ignore.case = TRUE))

  nonfinite <- list(expression = matrix(c(1, Inf), nrow = 1,
                                        dimnames = list("v1", c("g1", "g2"))),
                    methylation = matrix(c(1, 2), nrow = 1,
                                         dimnames = list("v1", c("cg1", "cg2"))))
  nonfinite_error <- tryCatch({ align_multimodal_yfb_prediction(nonfinite, spec); NULL },
                              error = conditionMessage)
  assert_true(grepl("finite", nonfinite_error, ignore.case = TRUE))
})

if (sys.nframe() == 0L) report_results("test_multimodal_yfb_preprocess_predict.R")
