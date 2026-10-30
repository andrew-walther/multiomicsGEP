source("tests/test_helpers.R")
source("results/benchmark_sim/run_yfb_supervision_decision.R")
run_test("YFBS calibration: dead factors cannot set the precision-share grid", {
  rec<-matrix(c(0,0,2,4),2,2)
  sur<-matrix(c(1,1,.2,.8),2,2)
  g<-yfb_calibrate_gamma(rec,sur,c(FALSE,TRUE),.2,.2)
  assert_finite(g)
  assert_near(median(.2*g*sur[,2]/(rec[,2]+.2*g*sur[,2])),.2, .01)
})
