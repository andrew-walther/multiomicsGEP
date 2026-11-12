# ============================================================
# Script: test_feature_group_informative.R
# Purpose: Verify two-direction survival and scoring for the informative Task D study.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: existing simulation runner and test helpers
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/benchmark_sim/run_feature_group_simulation.R")

run_test("Informative simulation: two unit directions generate matched training risk SD", {
  data <- simulate_multimodal_yfb_data("both_informative", n_train = 30L,
    n_validation = 15L, p_expression = 12L, p_methylation = 20L, seed = 71L)
  x <- fgsim_prepare(data, 71L, programs = 1:2)
  raw <- do.call(rbind, Map(function(F, s) sweep(F[, 1:2, drop = FALSE], 1, s, "/"),
                           data$truth$F, x$scaled$sd))
  expected_g <- sweep(raw, 2, sqrt(colSums(raw^2)), "/")
  assert_near(x$g, expected_g, 1e-12)
  assert_near(colSums(x$g^2), c(1, 1), 1e-12)
  direction <- drop(expected_g %*% data$truth$beta[1:2])
  assert_near(x$eta, drop(sweep(x$Y, 2, colMeans(x$Y)) %*% direction) * x$a, 1e-12)
  assert_near(x$eta_validation,
    drop(sweep(x$Y_validation, 2, colMeans(x$Y)) %*% direction) * x$a, 1e-12)
  assert_near(sd(x$eta), sd(data$training$eta), 1e-12)
  set.seed(71L + 12000L)
  et <- (-log(runif(30)) / (data$truth$weibull_scale * exp(x$eta)))^(1 / data$truth$weibull_shape)
  cb <- rexp(30)
  cs <- multimodal_yfb_calibrate_censoring(et, cb, data$truth$target_censoring)
  assert_near(x$time, pmin(et, cb * cs), 0)
  assert_equal(x$event, as.integer(et <= cb * cs))
  et <- (-log(runif(15)) / (data$truth$weibull_scale * exp(x$eta_validation)))^(1 / data$truth$weibull_shape)
  cb <- rexp(15)
  cs <- multimodal_yfb_calibrate_censoring(et, cb, data$truth$target_censoring)
  assert_near(x$time_validation, pmin(et, cb * cs), 0)
  assert_equal(x$event_validation, as.integer(et <= cb * cs))
  arms <- fgsim_arm_inputs(x, 3L)
  fgsim_assert_inputs(arms$grouped, arms$pooled)
  fgsim_assert_expression_inputs(arms$grouped, arms$expression_only)
  j <- x$groups == "expression"
  assert_near(arms$expression_only$g,
    sweep(x$g[j, , drop = FALSE], 2, sqrt(colSums(x$g[j, , drop = FALSE]^2)), "/"), 1e-12)
  for (arm in arms) assert_equal(arm$eta_validation, x$eta_validation)
})

run_test("Informative simulation: both true survival programs are scored separately", {
  g <- cbind(c(-2, -1, 0, 1, 2, 0), c(2, -1, -2, -1, 2, 0))
  fit <- list(EF = cbind(g[, 2], -g[, 1], c(0, 1, 0, -1, 0, 1)),
    EBeta = c(2, -2, 2), VBeta = rep(.01, 3), EBeta2 = rep(100, 3))
  out <- fgsim_recovery_two(fit, g, rep(c("expression", "methylation"), each = 3), 1.96)
  assert_true(out$recovered)
  assert_true(out$recovered2)
  assert_equal(out$matched_program, 2L)
  assert_equal(out$matched_program2, 1L)
  assert_equal(out$false_active, 1L)
  fit$EBeta[1] <- -2
  wrong <- fgsim_recovery_two(fit, g,
    rep(c("expression", "methylation"), each = 3), 1.96)
  assert_false(wrong$recovered2)
  assert_equal(wrong$false_active, 1L)
  collision <- fgsim_recovery_two(fit, cbind(g[, 1], g[, 1]),
    rep(c("expression", "methylation"), each = 3), 1.96)
  assert_true(collision$match_collision)
  assert_equal(collision$false_active, 2L)
  weak <- fit
  weak$EF <- matrix(rep(c(0, 1, 0, -1, 0, 1), 3), 6)
  no_match <- fgsim_recovery_two(weak, g,
    rep(c("expression", "methylation"), each = 3), 1.96)
  assert_false(no_match$recovered)
  assert_false(no_match$recovered2)
  assert_equal(no_match$false_active, 3L)
})

run_test("Informative simulation: separation is assessed only where expression passes", {
  summary <- expand.grid(feature_ratio = c(1L, 5L), methylation_active_fraction = .03,
    arm = c("grouped", "pooled", "expression_only"), stringsAsFactors = FALSE)
  summary$passed <- c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)
  result <- fgsim_separation(summary)
  assert_true(result$testable[result$feature_ratio == 1])
  assert_equal(result$conclusion[result$feature_ratio == 1], "separate priors helped")
  assert_false(result$testable[result$feature_ratio == 5])
  assert_equal(result$conclusion[result$feature_ratio == 5], "not testable")
  summary$passed[summary$arm == "expression_only"] <- TRUE
  result <- fgsim_separation(summary)
  assert_equal(result$conclusion[result$feature_ratio == 5], "separate priors hurt")
})

if (sys.nframe() == 0L) {
  results <- report_results("Informative feature-group simulation")
  if (results$failed > 0L) quit(status = 1L)
}
