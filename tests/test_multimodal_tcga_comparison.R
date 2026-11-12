# ============================================================
# Script: test_multimodal_tcga_comparison.R
# Purpose: Verify the matched-TCGA comparison's frozen analysis rules.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: test helpers and matched-TCGA comparison runner
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
if (file.exists("results/multimodal_real/run_multimodal_tcga_comparison.R"))
  source("results/multimodal_real/run_multimodal_tcga_comparison.R")

run_test("Matched TCGA: saved folds and training-only scaling are preserved", {
  design <- list(expression = cbind(gene = c(1, 2, 3, 100, 200)),
    methylation = cbind(cpg = c(4, 6, 8, 100, 200)), ranked = "cpg",
    time = 1:5, event = c(1, 0, 1, 0, 1), fold = c(1L, 2L, 3L, 4L, 5L))
  result <- mtcga_prepare_fold(design, 4L, n_cpg = 1L)
  assert_equal(result$tr, c(1L, 2L, 3L, 5L))
  assert_equal(result$te, 4L)
  assert_equal(result$fold, design$fold)
  assert_near(result$training$expression[, 1],
    (c(1, 2, 3, 200) - mean(c(1, 2, 3, 200))) / sd(c(1, 2, 3, 200)))
  assert_near(result$test_data$methylation[1, 1],
    (100 - mean(c(4, 6, 8, 200))) / sd(c(4, 6, 8, 200)))
})

run_test("Matched TCGA: grouped priors must contain valid separate fitted priors", {
  prior <- list(family = "point_laplace", pi = 0.5, rate = 2, point_mass = FALSE)
  fit <- list(prior_F = list(expression = rep(list(prior), 5),
                             methylation = rep(list(prior), 5)))
  assert_true(isTRUE(mtcga_assert_priors(fit, 5L)))
  fit$prior_F$methylation <- fit$prior_F$methylation[-1]
  assert_true(inherits(try(mtcga_assert_priors(fit, 5L), silent = TRUE), "try-error"))
})

run_test("Matched TCGA: paired SE uses all five fold differences", {
  expression <- c(.5, .6, .5, .6, .5)
  joint <- c(.6, .6, .7, .5, .8)
  differences <- joint - expression
  result <- mtcga_paired(joint, expression)
  assert_near(unname(result), c(mean(differences), sd(differences) / sqrt(5)))
  joint[2] <- NA_real_
  assert_true(all(is.na(mtcga_paired(joint, expression))))
})

run_test("Matched TCGA: joint-only requires no match in either modality", {
  v <- c(-2, -1, 0, 1, 2)
  w <- c(2, -1, -2, -1, 2)
  expression <- matrix(v, ncol = 1)
  methylation <- cbind(w, rep(0, 5))
  joint <- cbind(c(v, w), c(w, v), c(w, w), rep(0, 10))
  result <- mtcga_joint_only(joint, expression, methylation,
    c(rep("expression", 5), rep("methylation", 5)))
  assert_equal(result$joint_only, c(FALSE, TRUE, FALSE, FALSE))
  undefined <- mtcga_joint_only(joint, matrix(0, 5, 1), methylation,
    c(rep("expression", 5), rep("methylation", 5)))
  assert_true(all(undefined$expression_defined == 0))
  assert_true(all(!undefined$correlations_defined))
  assert_true(all(!undefined$joint_only))
})

if (sys.nframe() == 0L) {
  results <- report_results("Matched-TCGA comparison")
  if (results$failed > 0L) quit(status = 1L)
}
