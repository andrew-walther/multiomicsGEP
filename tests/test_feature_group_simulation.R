# ============================================================
# Script: test_feature_group_simulation.R
# Purpose: Check training-only scaling, precision shares and recovery scoring.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: test helpers and feature-group simulation runner
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
if (file.exists("results/benchmark_sim/run_feature_group_simulation.R"))
  source("results/benchmark_sim/run_feature_group_simulation.R")

run_test("Feature-group simulation: validation scaling uses training means and SDs", {
  train <- list(rna = cbind(1:3, 2 * (1:3)), meth = matrix(c(4, 6, 8), 3))
  valid <- list(rna = matrix(c(4, 8), 1), meth = matrix(10, 1))
  out <- fgsim_scale(train, valid)
  assert_near(out$training$rna[, 1], c(-1, 0, 1), 0)
  assert_near(out$validation$rna, matrix(c(2, 2), 1), 0)
  assert_near(out$validation$meth, matrix(2, 1), 0)
  assert_near(out$sd$rna, c(1, 2), 0)
})

run_test("Feature-group simulation: modality precision shares match hand calculation", {
  fit <- list(Tau = c(2, 4, 3), EF2 = matrix(1:6, 3, 2))
  out <- fgsim_precision(fit, c("rna", "rna", "meth"))
  assert_near(out$precision[out$program == 1], c(10, 9), 0)
  assert_near(out$share[out$program == 1], c(10, 9) / 19, 1e-14)
  assert_near(out$precision[out$program == 2], c(28, 18), 0)
  fit$EF2[, 2] <- 0
  dead <- fgsim_precision(fit, c("rna", "rna", "meth"))
  assert_true(all(is.na(dead$share[dead$program == 2])))
})

run_test("Feature-group simulation: recovery requires loading match and survival activity", {
  truth <- list(rna = diag(3), meth = diag(3))
  fit <- list(EF = cbind(rbind(diag(3), diag(3)), 1:6, rep(0, 6)),
              EBeta = c(0, 0, 2, 2, 0), EBeta2 = c(0.01, 0.01, 4.01, 4.01, 0.01))
  result <- fgsim_recovery(fit, truth, list(rna = rep(1, 3), meth = rep(1, 3)), 1.96)
  assert_true(result$recovered)
  assert_equal(result$matched_program, 3L)
  assert_equal(result$false_active, 1L)
  fit$EBeta[3] <- 0; fit$EBeta2[3] <- 0.01
  assert_false(fgsim_recovery(fit, truth,
    list(rna = rep(1, 3), meth = rep(1, 3)), 1.96)$recovered)
  fit$EBeta[3] <- 2; fit$EBeta2[3] <- 4.01
  fit$EF[, 3] <- rep(0, 6)
  assert_false(fgsim_recovery(fit, truth,
    list(rna = rep(1, 3), meth = rep(1, 3)), 1.96)$recovered)
})

if (sys.nframe() == 0L) {
  results <- report_results("Feature-group simulation")
  if (results$failed > 0L) quit(status = 1L)
}
