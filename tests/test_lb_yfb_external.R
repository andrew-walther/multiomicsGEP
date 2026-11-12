# ============================================================
# Script: test_lb_yfb_external.R
# Purpose: Check fixed external scoring and donor-aware pooled summaries.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival and shared test helpers
# ============================================================
if (!exists('run_test')) source('tests/test_helpers.R')
source('results/benchmark_sim/score_lb_yfb_external.R')
source('code/concordance_ci.R')
run_test('External alignment retains training order and fills absent genes', {
  a <- lb_yfb_external_align(matrix(c(1, 2, 3, 4), 2), c('b', 'a'), c('a', 'missing', 'b'))
  assert_equal(a$Y, matrix(c(3, 4, 0, 0, 1, 2), 2))
  assert_equal(a$missing, 1L)
})
run_test('External score uses fixed fit and normalization arguments', {
  old_lb <- if (exists('predict_ssbmf_lb_foldin', .GlobalEnv)) get('predict_ssbmf_lb_foldin', .GlobalEnv) else NULL
  old_yf <- if (exists('predict_ssbmf_yfb', .GlobalEnv)) get('predict_ssbmf_yfb', .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_lb)) rm('predict_ssbmf_lb_foldin', envir = .GlobalEnv) else assign('predict_ssbmf_lb_foldin', old_lb, .GlobalEnv)
    if (is.null(old_yf)) rm('predict_ssbmf_yfb', envir = .GlobalEnv) else assign('predict_ssbmf_yfb', old_yf, .GlobalEnv)
  }, add = TRUE)
  assign('predict_ssbmf_lb_foldin', function(Y, fit) list(risk_scores = drop(Y %*% fit$weights), converged = TRUE), .GlobalEnv)
  sys.source('code/predict_ssbmf_yfb.R', envir = .GlobalEnv)
  Y <- matrix(1:6, 3, 2); fits <- list(LB = list(weights = c(1, -1)),
    YFB = list(EF = matrix(c(1, 2), 2), EBeta = 3, EF_norms = 2, projection_center = c(1, 2)))
  scores <- lb_yfb_external_predict(Y, fits)
  assert_near(scores$LB, drop(Y %*% c(1, -1)))
  assert_near(scores$YFB, drop(sweep(Y, 2, c(1, 2)) %*% c(1, 2)) * 1.5)
})
run_test('Donor bootstrap respects overlap and pooled estimands', {
  t1 <- 1:20; t2 <- 1:30
  keep <- list(PACA_AU_array = list(LB = sin(t1), YFB = cos(t1 / 2), time = t1,
                                   status = rep(1, 20), donor = paste0('d', 1:20)),
               PACA_AU_seq = list(LB = cos(t2 / 3), YFB = sin(t2 / 4), time = t2,
                                 status = rep(1, 30), donor = paste0('d', 11:40)))
  assert_equal(lb_yfb_external_overlap(lapply(keep, `[[`, 'donor'), 10L), 10L)
  assert_true(inherits(tryCatch(lb_yfb_external_overlap(lapply(keep, `[[`, 'donor'), 47L), error = identity), 'error'))
  result <- lb_yfb_external_pooled_ci(keep, B = 100, seed = 1)
  assert_true(result$se[result$metric == 'equal_diff'] > 0)
  assert_true(abs(result$estimate[result$metric == 'equal_diff'] -
                  result$estimate[result$metric == 'pair_diff']) > 1e-4)
  reference <- bootstrap_concordance_diff_ci_stratified(lapply(keep, `[[`, 'YFB'), lapply(keep, `[[`, 'LB'),
    lapply(keep, `[[`, 'time'), lapply(keep, `[[`, 'status'), B = 100, seed = 1,
    donor_list = lapply(keep, `[[`, 'donor'))
  assert_near(result$lower[result$metric == 'equal_diff'], reference$lower)
  assert_near(result$upper[result$metric == 'equal_diff'], reference$upper)
})
if (sys.nframe() == 0L) {
  r <- report_results('test_lb_yfb_external.R')
  if (r$failed > 0) quit(status = 1)
}
