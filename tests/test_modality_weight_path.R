# ============================================================
# Script: test_modality_weight_path.R
# Purpose: Check weight-path matching, usability and exact reproducibility.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: existing test helpers
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/benchmark_sim/run_modality_weight_simulation.R")
run_test("Weight path: block correlation is undefined for absent truth or empty fit", {
  assert_equal(mwpath_cor(1:4, 2*(1:4)), 1)
  assert_true(is.na(mwpath_cor(rep(0,4), 1:4)))
  assert_true(is.na(mwpath_cor(1:4, rep(0,4))))
})
run_test("Weight path: all four usability conditions and larger-weight tie break", {
  a <- expand.grid(fraction=c(.15,.03), weight=c(1,.7,.3))
  a$recovered <- 5; a$mean_false_active <- .5; a$methylation_retained <- 8
  b <- expand.grid(truth_program=1:3, weight=c(1,.7,.3))
  b$loading_and_pattern <- 3
  assert_equal(unique(mwpath_select(a,b)$selected_weight), NA_real_)
  b$loading_and_pattern[b$truth_program==3 & b$weight<1] <- 4
  assert_equal(unique(mwpath_select(a,b)$selected_weight), .7)
  bad <- a; bad$recovered[bad$weight<1] <- 4
  assert_true(all(!mwpath_select(bad,b)$usable))
  badB <- b; badB$loading_and_pattern[badB$truth_program==2 & badB$weight<1] <- 2
  assert_true(all(!mwpath_select(a,badB)$usable))
  a$methylation_retained[a$weight==.7] <- 7
  assert_equal(unique(mwpath_select(a,b)$selected_weight), .3)
  a$mean_false_active[a$weight==.3] <- 1.1
  assert_equal(unique(mwpath_select(a,b)$selected_weight), NA_real_)
})
run_test("Weight path: reproducibility rejects a changed outcome", {
  x <- data.frame(block='A',methylation_active_fraction=.15,seed=1,weight=1,
    recovered=TRUE,recovered2=FALSE,false_active=1,converged=TRUE,n_iter=19,
    empty_blocks=2,point_mass_resets=26,full_pattern_success=NA,
    c_validation=.763432073544433,n_warnings=0,error='')
  old <- read.csv('results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weights/metrics.csv',colClasses=c(error='character'))
  assert_true(mwpath_reproduce(x, old))
  x$recovered <- FALSE
  assert_true(grepl('Reproducibility',tryCatch({mwpath_reproduce(x,old); ''},error=conditionMessage)))
})
if (sys.nframe()==0L) {
  r <- report_results('Weight path'); if(r$failed>0L) quit(status=1L)
}
