# ============================================================
# Script: test_modality_weight_simulation.R
# Purpose: Check fixed recovery patterns and registered weight-rule selection.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: test helpers and modality-weight simulation runner
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/benchmark_sim/run_modality_weight_simulation.R")
run_test("Weighted simulation: all three planted modality patterns are necessary", {
  set.seed(81)
  g <- matrix(0, 20, 3)
  g[1:10, 1] <- rnorm(10); g[11:20, 2] <- rnorm(10); g[, 3] <- rnorm(20)
  g <- sweep(g, 2, sqrt(colSums(g^2)), "/")
  fit <- list(EF = g, EF2 = g^2, Tau = rep(1, 20), EBeta = c(.7, .35, 0),
    EBeta2 = c(.7, .35, 0)^2 + .0001)
  groups <- rep(c("expression", "methylation"), each = 10)
  a <- mwsim_patterns(fit, g, groups, 1.96)
  assert_true(all(a$loading_recovered & a$pattern_correct))
  assert_true(all(a$survival_recovered[1:2]))
  fit$EF2[11:20, 1] <- 0.1
  b <- mwsim_patterns(fit, g, groups, 1.96)
  assert_false(b$pattern_correct[1])
  fit$EBeta[2] <- -.35
  b <- mwsim_patterns(fit, g, groups, 1.96)
  assert_false(b$survival_recovered[2])
})
run_test("Weighted simulation: selection requires rescued imbalance and passing control", {
  a <- expand.grid(feature_ratio = c(1L, 5L), methylation_active_fraction = c(.03, .15),
    arm = c("grouped_unweighted", "grouped_count", "grouped_eigenvalue"), stringsAsFactors = FALSE)
  a$passed <- a$arm != "grouped_unweighted"
  a$recovered <- 8L
  b <- expand.grid(feature_ratio = c(1L, 5L), arm = c("grouped_count", "grouped_eigenvalue"),
    stringsAsFactors = FALSE)
  b$full_pattern_successes <- ifelse(b$arm == "grouped_count", 3L, 4L)
  assert_equal(mwsim_select(a, b)$selected_rule, "eigenvalue")
  a$passed[a$arm == "grouped_eigenvalue" & a$feature_ratio == 1] <- FALSE
  assert_equal(mwsim_select(a, b)$selected_rule, "feature_count")
  a$passed <- FALSE
  assert_equal(mwsim_select(a, b)$selected_rule, "none")
  a$passed <- a$arm != "grouped_unweighted"
  b$full_pattern_successes <- 3L
  assert_equal(mwsim_select(a, b)$selected_rule, "feature_count")
})
if (sys.nframe() == 0L) {
  result <- report_results("Weighted simulation")
  if (result$failed > 0L) quit(status = 1L)
}
