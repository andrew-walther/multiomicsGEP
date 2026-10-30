# Independent checks for diagnostic-only noise and objective controls. ----
source("code/fit_yfb_supervised.R")
if(file.exists("code/yfb_diagnostic_controls.R"))source("code/yfb_diagnostic_controls.R")
run_test("YFB screen: defaults reproduce the pre-change synthetic fixture", local({
  d<-readRDS("tests/fixtures/yfb-screen-default.rds")
  f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=3,gamma=2,max_iter=8)
  for(nm in c("EL","EL2","EBeta","VBeta","Tau","eta"))assert_near(f[[nm]],d$fit[[nm]],1e-12)
  assert_near(f$EF[[1]],d$fit$EF[[1]],1e-12)
  assert_near(f$history$objective,d$fit$history$objective,1e-12)
}))
run_test("YFB screen: fixed precision stays equal to the supplied vector", local({
  d<-readRDS("tests/fixtures/yfb-screen-default.rds");tau<-seq(.5,2,length.out=ncol(d$Y))
  f<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=3,gamma=2,
    max_iter=8,tau_init=tau,tau_fixed=TRUE)
  assert_near(f$Tau,tau,0)
  err<-try(fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=3,
    tau_fixed=TRUE),silent=TRUE)
  assert_true(inherits(err,"try-error"))
}))
run_test("YFB screen: risk-set curvature bound dominates the exact Hessian", local({
  time<-c(2,3,2,4);event<-c(1,0,1,1);eta<-c(-1,.3,.4,2)
  D<-yfb_cox_curvature_bound(time,event)
  assert_near(D,c(1,1,1,1.5),1e-12)
  H<-matrix(0,4,4)
  for(i in which(event==1)) {
    ii<-which(time>=time[i]);p<-exp(eta[ii]-max(eta[ii]));p<-p/sum(p)
    H[ii,ii]<-H[ii,ii]+diag(p)-outer(p,p)
  }
  assert_true(min(eigen(diag(D)-H,symmetric=TRUE,only.values=TRUE)$values)> -1e-12)
  a<-1.2
  exact<--log(2*cosh(a));lower<--log(2)-.5*a^2
  assert_true(lower<=exact)
  assert_false(yfb_accept_bound(1,NaN));assert_false(yfb_accept_bound(1,.9))
}))
run_test("YFB screen: rejected sweeps restore state and cannot report convergence", local({
  d<-readRDS("tests/fixtures/yfb-screen-default.rds")
  first<-fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=3,gamma=2,
    max_iter=1,guard=TRUE)
  original<-yfb_accept_bound
  assign("yfb_accept_bound",function(...)FALSE,envir=.GlobalEnv)
  f<-tryCatch(fit_yfb_supervised(list(expression=d$Y),d$time,d$event,K=3,gamma=2,
    max_iter=8,guard=TRUE),finally=assign("yfb_accept_bound",original,envir=.GlobalEnv))
  assert_false(f$converged);assert_true("backtrack_stalled" %in% f$flags)
  for(nm in c("EL","EL2","EBeta","VBeta","Tau","eta"))assert_near(f[[nm]],first[[nm]],1e-12)
  assert_near(f$EF[[1]],first$EF[[1]],1e-12)
  assert_equal(f$g_beta,first$g_beta)
  assert_equal(f$g_F,first$g_F)
  assert_near(f$kl$F,first$kl$F,1e-12)
  assert_equal(nrow(f$history),1L)
  assert_true(f$rejected_trials>=1)
}))
run_test("YFB screen: normalized criterion equals its independently calculated value", local({
  Y<-matrix(c(1,3),2,1);L<-matrix(c(1,2),2,1);F<-matrix(.5,1,1)
  VL<-matrix(c(.1,.2),2,1);VF<-matrix(.2,1,1)
  st<-list(L=L,F=F,L2=L^2+VL,F2=F^2+VF,VL=VL,VF=VF,
    b=.4,VBeta=.3,tau=2,klL=.2,klF=.3,klB=.4)
  # Residual energy 4.25; factor uncertainty 1.135. Normalized projection
  # means (-1,1) have variance .8, giving Var(η_i)=.3+.8*(.16+.3)=.668.
  expected<-log(2/(2*pi))-5.385-.5+2*(-log1p(exp(.8))-.334)-.4
  got<-yfb_bound_objective(Y,c(1,2),c(1,0),st,2)
  assert_near(got,expected,1e-12)
  cc<-3;st$L<-st$L/cc;st$L2<-st$L2/cc^2;st$VL<-st$VL/cc^2
  st$F<-st$F*cc;st$F2<-st$F2*cc^2;st$VF<-st$VF*cc^2
  assert_near(yfb_bound_objective(Y,c(1,2),c(1,0),st,2),expected,1e-12)
}))
