# ============================================================
# Script: test_multimodal_comparator_rank.R
# Purpose: Verify reporting of actual retained two-step factor counts
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: tests/test_helpers.R; multimodal comparison runner
# Run: Rscript tests/test_multimodal_comparator_rank.R
# ============================================================
source('tests/test_helpers.R')
source('results/multimodal_real/run_supervised_multimodal_comparison.R')
run_test('Two-step rank reports retained columns, not greedy ceiling', {
  x<-list(n_factors=6L,loadings=matrix(1,10,6),K_init=8L)
  assert_equal(retained_two_step_rank(x),6L)
  x$n_factors<-8L
  err<-tryCatch(retained_two_step_rank(x),error=function(e)e)
  assert_true(inherits(err,'error'))
})
r<-report_results('rank metadata');if(r$failed)quit(status=1)
run_test('ICGC primary scoring uses the same platform standardization as all donors', {
  all<-matrix(c(1,3,10,20,2,5,7,30),4,2,
              dimnames=list(c('d1','d2','d3','d4'),c('g1','g2')))
  z<-standardize_icgc_expression(all,c('d3','d1'))
  assert_near(z$validation_all,scale(all),tol=1e-12)
  assert_true(identical(z$validation_primary,z$validation_all[c('d3','d1'),,drop=FALSE]))
  b<-c(1,2)
  assert_near(as.vector(z$validation_primary%*%b),
              as.vector(z$validation_all[c('d3','d1'),,drop=FALSE]%*%b),tol=1e-12)
  err<-tryCatch(standardize_icgc_expression(all,'missing'),error=function(e)e)
  assert_true(inherits(err,'error'))
})
r<-report_results('rank and platform scoring');if(r$failed)quit(status=1)
