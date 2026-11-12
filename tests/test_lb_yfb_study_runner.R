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
run_test('Runner validates rank and survival weight', {
  base <- c('--mode=smoke', '--scaling=zscore', '--prior_LF=point_laplace',
            '--prior_beta=normal', '--out=/tmp/example')
  assert_equal(lb_yfb_runner_options(base)$w_S, 1)
  for (bad in c('--K=0', '--K=-1', '--K=1.5', '--K=Inf', '--K=bad',
                '--w_S=0', '--w_S=-1', '--w_S=Inf', '--w_S=bad')) {
    assert_true(inherits(tryCatch(lb_yfb_runner_options(c(base, bad)), error = identity), 'error'))
  }
})
run_test('Smoke records and applies explicit rank and survival weight', {
  run_id <- paste0('runner_test_', Sys.getpid())
  out <- file.path(tempdir(), run_id)
  checkpoint <- file.path('data', paste0('lb_yfb_study_', run_id))
  on.exit(unlink(c(out, checkpoint), recursive = TRUE), add = TRUE)
  opts <- lb_yfb_runner_options(c('--mode=smoke', '--scaling=zscore',
    '--prior_LF=point_laplace', '--prior_beta=normal', paste0('--out=', out),
    '--K=3', '--w_S=2'))
  suppressWarnings(lb_yfb_runner_main(opts))
  settings <- readLines(file.path(out, 'settings.txt'))
  assert_equal(settings[grepl('^K=', settings)], 'K=3')
  assert_equal(settings[grepl('^w_S=', settings)], 'w_S=2')
  saved <- readRDS(file.path(checkpoint, 'patient_level.rds'))
  assert_equal(ncol(saved$pair$fits$LB$EL), 3L)
  assert_equal(ncol(saved$pair$fits$YFB$EL), 3L)
  assert_equal(saved$pair$fits$LB$w_S, 2)
  assert_equal(saved$pair$fits$YFB$w_S, 2)
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
  assert_true(all(c('C_within', 'C_pooled', 'pairs_within', 'pairs_pooled',
                    'C_simulation', 'pairs_simulation') %in% names(result$metrics)))
})
run_test('Runner transforms molecular loading truth by inverse feature SD', {
  X <- matrix(c(1, 3, 2, 6), 2)
  assert_near(lb_yfb_runner_loading_truth(c(1, 2), X), c(1, 2) / apply(X, 2, sd))
})
if (sys.nframe() == 0L) {
  r <- report_results('test_lb_yfb_study_runner.R')
  if (r$failed > 0) quit(status = 1)
}
