# ============================================================
# Script: test_supervision_comparators.R
# Purpose: Verify Cox-z program selection and signed planted recovery
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: tests/test_helpers.R; comparator scoring definitions
# Run: Rscript tests/test_supervision_comparators.R
# ============================================================
source('tests/test_helpers.R')
source('results/benchmark_sim/run_supervision_comparators.R')
run_test('Comparator program follows Cox z rather than truth correlation; sign matters', {
  set.seed(2);truth<-matrix(rnorm(90),30,3);F<-cbind(truth[,3],truth[,1])
  tb<-list(F_true=truth,beta_true=c(0,0,3.3))
  wrong<-comparator_metrics(F,c(1,1),c(1,5),tb)
  right<-comparator_metrics(F,c(1,1),c(5,1),tb)
  reversed<-comparator_metrics(F,c(-1,1),c(5,1),tb)
  uns<-comparator_metrics(F,NULL,NULL,tb,TRUE)
  assert_true(!wrong$recovered && wrong$program==2L)
  assert_true(right$recovered && right$program==1L)
  assert_true(!reversed$recovered)
  assert_true(uns$recovered && is.na(uns$false_active))
})
r<-report_results('comparators');if(r$failed)quit(status=1)
