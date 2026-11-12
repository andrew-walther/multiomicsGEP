# ============================================================
# Script: test_lb_yfb_study_summary.R
# Purpose: Check simulation gates and program summaries.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R and project test and enrichment helpers
# ============================================================
if (!exists('run_test')) source('tests/test_helpers.R')
source('results/benchmark_sim/summarize_lb_yfb_study.R')
source('code/pathway_enrichment.R')
run_test('Simulation summary uses registered paired exclusion and denominators', {
  x <- expand.grid(seed = 1:10, model = c('LB', 'YFB'))
  x$design <- 'S9'; x$share <- .01; x$amp <- NA_real_; x$w_S <- 1
  x$error <- ''; x$recovery_success <- TRUE; x$false_active <- 0; x$converged <- TRUE
  x$validation_C <- .6; x$matched_score_recovery <- NA_real_
  x$recovery_success[x$model == 'LB' & x$seed <= 3] <- FALSE
  x$validation_C[x$model == 'LB' & x$seed == 1] <- NA_real_
  path <- tempfile(fileext = '.csv'); on.exit(unlink(path), add = TRUE)
  write.csv(x, path, row.names = FALSE)
  z <- lb_yfb_summary_simulation(path)
  assert_equal(z$gate_pass[z$model == 'LB'], FALSE)
  assert_equal(z$relative_exclusion[z$model == 'LB'], TRUE)
  assert_equal(z$relative_exclusion[z$model == 'YFB'], FALSE)
  assert_equal(z$validation_C_defined[z$model == 'LB'], 9L)
  x$error[x$model == 'YFB' & x$seed == 1] <- 'failed'
  write.csv(x, path, row.names = FALSE)
  z <- lb_yfb_summary_simulation(path)
  assert_equal(z$gate_pass[z$model == 'YFB'], FALSE)
  assert_equal(z$failed[z$model == 'YFB'], 1L)
})
run_test('Simulation CSV round trip preserves blank success messages', {
  x <- expand.grid(seed = 1:10, model = c('LB', 'YFB'))
  x$design <- 'S9'; x$share <- .01; x$amp <- NA_real_; x$w_S <- 1
  x$error <- ''; x$recovery_success <- TRUE; x$false_active <- 0; x$converged <- TRUE
  x$validation_C <- .6; x$matched_score_recovery <- NA_real_
  path <- tempfile(fileext = '.csv'); on.exit(unlink(path), add = TRUE)
  write.csv(x, path, row.names = FALSE)
  z <- lb_yfb_summary_simulation(path)
  assert_true(all(z$failed == 0))
  assert_true(all(z$gate_pass))
})
run_test('Top genes use absolute weights and dead programs have no arbitrary genes', {
  fit <- list(EF = cbind(c(-3, 2, 1), c(0, 0, 0)), EBeta = c(1, 0), VBeta = c(.04, 0))
  top <- lb_yfb_summary_top(fit, c('a', 'b', 'c'), 'LB', 2)
  assert_equal(top$top$gene, c('a', 'b'))
  assert_equal(top$top$weight, c(-3, 2))
  assert_equal(top$counts$n_top, c(2L, 0L))
  p <- lb_yfb_summary_programs(fit, 'LB', 'full')
  assert_equal(p$survival_active, c(TRUE, FALSE))
  assert_near(p$hazard_ratio[1], exp(1))
})
run_test('Sparse live programs exclude zero weights from gene overlaps', {
  fit <- list(EF = cbind(c(0, -2, 0), c(0, 0, 0)))
  top <- lb_yfb_summary_top(fit, c('a', 'b', 'c'), 'LB', 20)
  assert_equal(top$top$gene, 'b')
  assert_equal(top$top$weight, -2)
  assert_equal(top$counts$n_top, c(1L, 0L))
  overlap <- compute_geneset_overlap(top$top$gene, c('a', 'c'), c('a', 'b', 'c'))
  assert_equal(overlap$overlap_n, 0L)
})
if (sys.nframe() == 0L) {
  r <- report_results('test_lb_yfb_study_summary.R')
  if (r$failed > 0) quit(status = 1)
}
