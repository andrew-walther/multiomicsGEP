# ============================================================
# Script: test_working_cycle_audit.R
# Purpose: Verify independent posterior-density and anchored-objective observation.
# Created: 2026-10-16
# Dependencies: existing test helpers and fitting modules
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("results/benchmark_sim/working_cycle_audit_helpers.R"))
  source("results/benchmark_sim/working_cycle_audit_helpers.R")

run_test("Cycle-T1: mixed posterior entropy agrees with independent quadrature", {
  g<-list(family="point_laplace",pi=.4,rate=1.3,point_mass=FALSE)
  p<-multimodal_yfb_signed_posterior(4,2,g)
  q<-cycle_posterior(4,2,g,p$mean,p$second,p$log_ml)
  x<-.5;s<-.5
  ml<-exp(p$log_ml)
  spike<-(1-g$pi)*dnorm(x,0,s)/ml
  density<-function(f)g$pi*g$rate/2*exp(-g$rate*abs(f))*dnorm(x-f,0,s)/ml
  integrand<-function(f){z<-density(f);ifelse(z>0,z*log(z),0)}
  direct<-spike*log(spike)+integrate(integrand,-Inf,0)$value+
    integrate(integrand,0,Inf)$value
  assert_near(q$elogq,direct,tol=1e-8,msg="mixed measure entropy")
  assert_near(cycle_kl(q,g),multimodal_yfb_coordinate_kl(p$x,p$s2,
    p$mean,p$second,p$log_ml),tol=1e-12,msg="independent KL identity")
})
run_test("Cycle-T2: support conflicts and no-information distributions stay explicit", {
  g<-list(family="point_laplace",pi=.5,rate=1,point_mass=FALSE)
  q<-cycle_posterior(0,0,g,0,1,NA_real_)
  assert_near(cycle_kl(q,g),0,tol=1e-12)
  dead<-list(family="point_laplace",pi=0,rate=NA_real_,point_mass=TRUE)
  assert_true(identical(cycle_kl(q,dead),Inf),msg="support cannot be silently clipped")
  assert_near(cycle_kl(cycle_posterior(1,0,dead,0,0,0),dead),0,tol=0)
})
run_test("Cycle-T3: scaling transforms the continuous entropy and prior together", {
  g<-list(family="point_laplace",pi=.7,rate=.8,point_mass=FALSE)
  p<-multimodal_yfb_signed_posterior(3,-2,g)
  q<-cycle_posterior(3,-2,g,p$mean,p$second,p$log_ml)
  z<-cycle_scale(q,g,5)
  assert_near(cycle_kl(z$q,z$g),cycle_kl(q,g),tol=1e-12)
  assert_near(z$q$mean,5*q$mean,tol=1e-12)
  assert_near(z$q$elogq,q$elogq-sum(q$prob)*log(5),tol=1e-12)
})
run_test("Cycle-T4: anchored expected quadratic equals a dense covariance calculation", {
  Y<-matrix(c(1,3,2,-1,2,4),3,2);F<-matrix(c(.2,-.4),2,1)
  state<-list(Y=Y,EL=matrix(.1,3,1),EL2=matrix(.2,3,1),EF=F,
    EF2=F^2+.05,EBeta=.6,EBeta2=.5,Tau=c(1,2))
  w<-c(.3,1,.7);u<-c(.2,-.1,-.1);eta0<-c(-.2,.4,.1)
  work<-list(w=w,u=u,eta0=eta0,logPL=-2)
  got<-cycle_likelihood(state,work)
  H<-diag(w)-tcrossprod(w)/sum(w);Z<-Y%*%F
  covariance<-.5*Y%*%diag(rep(.05,2))%*%t(Y)+(.5-.6^2)*tcrossprod(Z)
  d<-as.vector(.6*Z)-eta0
  expected<--2+sum(u*d)-sum(d*(H%*%d))/2-sum(H*t(covariance))/2
  assert_near(got$survival,expected,tol=1e-12,msg="anchor and uncertainty retained")
  R2<-Y^2-2*Y*(state$EL%*%t(F))+state$EL2%*%t(state$EF2)
  rec<-3*sum(log(state$Tau))/2-sum(sweep(R2,2,state$Tau,"*"))/2-
    length(Y)*log(2*pi)/2
  assert_near(got$reconstruction,rec,tol=1e-12)
})
run_test("Cycle-T5: Normal and positive-slab densities have correct KL", {
  g<-list(family="normal",variance=2,point_mass=FALSE)
  p<-multimodal_yfb_signed_posterior(3,-1,g)
  q<-cycle_posterior(3,-1,g,p$mean,p$second,p$log_ml)
  v<-p$second-p$mean^2
  assert_near(cycle_kl(q,g),.5*((v+p$mean^2)/2-1-log(v/2)),tol=1e-12)
  g<-list(family="point_exponential",pi=1,rate=2,point_mass=FALSE)
  p<-multimodal_yfb_point_exponential_posterior(2,1,g)
  q<-cycle_posterior(2,1,g,p$mean,p$second,p$log_ml)
  assert_near(q$prob,1,tol=0)
  assert_near(cycle_kl(q,g),multimodal_yfb_coordinate_kl(p$x,p$s2,
    p$mean,p$second,p$log_ml),tol=1e-12)
})
run_test("Cycle-T6: fixed-prior loading posterior improves the full fixed quadratic", {
  Y<-matrix(c(1,3,2,-1,2,4),3,2);w<-c(.3,1,.7)
  s<-list(Y=Y,EL=matrix(.1,3,1),EL2=matrix(.2,3,1),EF=matrix(c(.2,-.4),2,1),
    EF2=matrix(c(.2,-.4)^2+.05,2,1),EBeta=.6,EBeta2=.5,Tau=c(1,2))
  g<-list(family="point_laplace",pi=.5,rate=1,point_mass=FALSE)
  yc<-profile_cox_center(Y,w);h<-c(.2,-.1,-.1)
  out<-multimodal_yfb_update_F_mk(yc,s$Tau,s$EL[,1],s$EL2[,1],Y,w,h,
    as.vector(yc%*%s$EF),as.vector(yc^2%*%rep(.05,2)),s$EF[,1],s$EF2[,1],
    s$EBeta,s$EBeta2,g)
  # Each coordinate's prior is fixed, and B excludes its other features.
  H<-diag(w)-tcrossprod(w)/sum(w)
  A<-s$Tau[1]*sum(s$EL2)+s$EBeta2*sum(Y[,1]*(H%*%Y[,1]))
  B<-s$Tau[1]*sum(Y[,1]*s$EL)+s$EBeta*sum(Y[,1]*h)-
    s$EBeta2*sum(Y[,1]*(H%*%(Y[,2]*s$EF[2,1])))
  assert_near(out$details$A[1],A,tol=1e-12)
  assert_near(out$details$B[1],B,tol=1e-12)
  old<-multimodal_yfb_signed_posterior(A,B-1,g)
  new<-multimodal_yfb_signed_posterior(A,B,g)
  qo<-cycle_posterior(A,B-1,g,old$mean,old$second,old$log_ml)
  qn<-cycle_posterior(A,B,g,new$mean,new$second,new$log_ml)
  score<-function(q)-A*q$second/2+B*q$mean-cycle_kl(q,g)
  assert_true(score(qn)>score(qo),msg="wrong residual or KL sign loses ascent")
})
run_test("Cycle-T7: observer preserves a complete small profiled fit", {
  set.seed(79);Y<-matrix(rnorm(48),12,4)
  a<-list(Y=Y,time=1:12,status=rep(c(1,0),6),K=2L,max_iter=3L,
    prior_LF="point_laplace",prior_beta="normal",canonicalize=TRUE,
    profile_cox_intercept=TRUE,w_S=1,verbose=FALSE)
  set.seed(3);expected<-do.call(fit_cox_on_yf,a)
  audit<-cycle_observer()
  set.seed(3);got<-audit$run(function()do.call(fit_cox_on_yf,a))
  assert_true(identical(got,expected),msg="all posterior fields/history unchanged")
  d<-do.call(rbind,audit$rows)
  assert_true(any(d$complete) && any(d$event=="refresh"),msg="complete-cycle coverage")
  assert_true(any(d$event=="F_prior") && any(d$event=="canonical"))
  assert_true(all(d$moment_gap[d$complete]<=1e-12),msg="densities match live moments")
  assert_false(any(d$implemented_common_prior),msg="transient priors labelled honestly")
})
