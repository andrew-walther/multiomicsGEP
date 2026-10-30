source("tests/test_helpers.R")
source("results/benchmark_sim/run_yfb_supervision_decision.R")
run_test("YFBS calibration: dead factors cannot set the precision-share grid", {
  rec<-matrix(c(0,0,2,4),2,2)
  sur<-matrix(c(1,1,.2,.8),2,2)
  g<-yfb_calibrate_gamma(rec,sur,c(FALSE,TRUE),.2,.2)
  assert_finite(g)
  assert_near(median(.2*g*sur[,2]/(rec[,2]+.2*g*sur[,2])),.2, .01)
})

source("results/benchmark_sim/summarize_yfb_supervision.R")
run_test("Gate S: all three boundaries must hold in two paired cells", {
  x<-data.frame(p=2000,variance=rep(c(.0025,.01,.05),2),
    arm=rep(c("unsupervised","share0.2"),each=3),n=10,
    recovered=c(1,1,1,4,4,1),mean_false=1,converged=9,errors=0)
  assert_true(yfb_gate_cells(x)$gate$passed)
  bad<-x;bad$mean_false[4]<-1.01
  assert_false(yfb_gate_cells(bad)$gate$passed)
  bad<-x;bad$converged[5]<-8
  assert_false(yfb_gate_cells(bad)$gate$passed)
  bad<-x;bad$recovered[4]<-3
  assert_false(yfb_gate_cells(bad)$gate$passed)
  bad<-x;bad$errors[5]<-1
  assert_false(yfb_gate_cells(bad)$gate$passed)
})
