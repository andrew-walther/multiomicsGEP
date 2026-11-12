# ============================================================
# Script: test_lb_yfb_study_runner.R
# Purpose: Check runner options, registered design and failure retention.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R and shared test helpers
# ============================================================
if (!exists('run_test')) source('tests/test_helpers.R')
source('results/benchmark_sim/run_lb_yfb_study.R')
run_test('Runner validates options and enumerates the registered design', {
  args <- lb_yfb_runner_options(c('--mode=simulation', '--scaling=zscore',
    '--prior_LF=point_laplace', '--prior_beta=normal', '--out=/tmp/example'))
  assert_equal(args$mode, 'simulation')
  assert_equal(nrow(lb_yfb_runner_design()), 104L)
  assert_true(inherits(tryCatch(lb_yfb_runner_options(c('--mode=smoke', '--scaling=bad',
    '--prior_LF=point_laplace', '--prior_beta=normal', '--out=/tmp/example')),
    error = identity), 'error'))
})
run_test('Runner retains failed simulation cells for both models', {
  cell <- lb_yfb_runner_design()[1, ]
  result <- lb_yfb_runner_capture(cell, stop('synthetic failure'))
  assert_equal(result$metrics$model, c('LB', 'YFB'))
  assert_true(all(!result$metrics$converged))
  assert_true(all(!result$metrics$recovery_success))
  assert_true(all(is.na(result$metrics$false_active)))
  assert_true(all(result$metrics$error == 'synthetic failure'))
  assert_equal(result$error, 'synthetic failure')
})
run_test('Runner transforms molecular loading truth by inverse feature SD', {
  X <- matrix(c(1, 3, 2, 6), 2)
  assert_near(lb_yfb_runner_loading_truth(c(1, 2), X), c(1, 2) / apply(X, 2, sd))
})
if (sys.nframe() == 0L) {
  r <- report_results('test_lb_yfb_study_runner.R')
  if (r$failed > 0) quit(status = 1)
}
