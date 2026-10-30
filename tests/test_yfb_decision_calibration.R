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
  assert_true(yfb_gate_cells(bad)$gate$passed)
})

run_test("Final-state precision share is invariant to reciprocal factor scaling", {
  set.seed(82);Y<-matrix(rexp(40*10),40,10)
  F<-matrix(rexp(30),10,3);L<-matrix(rexp(120),40,3)
  f<-list(EF=list(F),EF2=list(F^2+.01),EL2=L^2+.1,Tau=rep(2,10),
          EBeta=c(.2,-.3,.1),EBeta2=c(.2,-.3,.1)^2+.01,gamma=10,rho=.2)
  t<-rexp(40);e<-rep(c(0,1),20)
  a<-yfb_precision_share_final(f,Y,t,e)
  cc<-c(.1,2,5);f$EF[[1]]<-sweep(F,2,cc,"*")
  f$EF2[[1]]<-sweep(f$EF2[[1]],2,cc^2,"*")
  f$EL2<-sweep(f$EL2,2,cc^2,"/")
  assert_near(a,yfb_precision_share_final(f,Y,t,e),1e-12)
  f$gamma<-0
  assert_near(yfb_precision_share_final(f,Y,t,e),rep(0,3),1e-12)
})
run_test("Failed fits receive the registered conservative gate penalty", {
  x<-data.frame(arm=c("share0.05","share0.05","share0.05","lasso"),
    recovered=TRUE,false_active=c(2,0,4,NA),converged=c(TRUE,FALSE,FALSE,TRUE),
    error=c("","","numerical failure",""),flags=c("","max_iter","numerical_failure",""))
  z<-yfb_apply_failure_penalty(x,7)
  assert_equal(z$recovered,c(TRUE,FALSE,FALSE,TRUE))
  assert_near(z$false_active[1:3],c(2,7,7),1e-12)
  assert_true(is.na(z$false_active[4]))
  assert_equal(z$raw_recovered,rep(TRUE,4))
})
