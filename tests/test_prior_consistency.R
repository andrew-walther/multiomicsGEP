# ============================================================
# Script: test_prior_consistency.R
# Purpose: Check posterior absolute moments and KL changes when a molecular prior is refitted.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing test helpers, base R, molecular updates
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/fit_multimodal_yfb.R")
source("results/benchmark_sim/round2_prior_consistency_helpers.R")
run_test("PriorConsistency-T1: absolute posterior moments match independent quadrature", {
 g<-list(family="point_laplace",pi=.4,rate=1.3,point_mass=FALSE)
 for(x in c(-.5,.5)){
  a<-round2_laplace_abs(4,4*x,g);s<-.5
  den<-(1-g$pi)*dnorm(x,0,s)+integrate(function(f)
    g$pi*g$rate/2*exp(-g$rate*abs(f))*dnorm(x-f,0,s),-Inf,Inf)$value
  expected<-integrate(function(f)abs(f)*g$pi*g$rate/2*
    exp(-g$rate*abs(f))*dnorm(x-f,0,s),-Inf,Inf)$value/den
  assert_near(a$mean_abs,expected,tol=1e-8,msg="conditional absolute moment")
  q<-multimodal_yfb_signed_posterior(4,4*x,g)
  assert_near(a$slab_prob,q$slab_prob,tol=1e-12,msg="same posterior mixture")
 }
})
run_test("PriorConsistency-T2: scale changes preserve probability and rescale absolute moments", {
 g<-list(family="point_laplace",pi=.4,rate=1.3,point_mass=FALSE)
 a<-round2_laplace_abs(c(4,3),c(2,-2),g);g$rate<-g$rate/5
 b<-round2_laplace_abs(c(4,3)/25,c(2,-2)/5,g)
 assert_near(a$slab_prob,b$slab_prob,tol=1e-12,msg="unit-invariant slab probability")
 assert_near(5*a$mean_abs,b$mean_abs,tol=1e-12,msg="absolute moment units")
})
run_test("PriorConsistency-T3: point-prior support conflict gives infinite corrected KL", {
 old<-list(family="point_laplace",pi=.5,rate=1,point_mass=FALSE)
 new<-list(family="point_laplace",pi=0,rate=NA_real_,point_mass=TRUE)
 q<-round2_laplace_abs(4.4,0,old)
 before<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,old)
 after<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,new)
 assert_true(is.finite(before)&&identical(after,-Inf),msg="nonzero continuous posterior cannot have finite KL to delta0")
 assert_true(round2_prior_log_expectation(0,0,new)==0,msg="delta posterior and prior are compatible")
})
run_test("PriorConsistency-T4: closed-form variational EM improves the fixed-q prior term", {
 old<-list(family="point_laplace",pi=.5,rate=1,point_mass=FALSE)
 q<-round2_laplace_abs(c(4.4,2,8),c(0,1,-1),old)
 new<-old;new$pi<-mean(q$slab_prob);new$rate<-sum(q$slab_prob)/sum(q$mean_abs)
 a<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,old)
 b<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,new)
 assert_true(b>=a && new$pi>0 && new$pi<1,msg="same family, no artificial floor, improved expected prior log density")
 assert_near(round2_prior_log_expectation(q$slab_prob,q$mean_abs,old),a,tol=0,msg="unchanged prior has zero KL adjustment")
})
run_test("PriorConsistency-T5: a prior with no spike has exactly no posterior spike", {
 g<-list(family="point_laplace",pi=1,rate=1.3,point_mass=FALSE)
 q<-round2_laplace_abs(rep(3,21),seq(-2,2,length.out=21),g)
 assert_true(all(q$slab_prob==1),msg="support exact, independent of rounded mixture weights")
 assert_true(is.finite(round2_prior_log_expectation(q$slab_prob,q$mean_abs,g)),
   msg="roundoff does not create a false support conflict")
})
run_test("PriorConsistency-T6: observer returns the original update and detects the known KL gap", {
 prior<-list(family="point_laplace",pi=.5,rate=1,point_mass=FALSE)
 args<-list(matrix(1,4,1),1,rep(1,4),rep(1.1,4),matrix(0,4,1),
   rep(1,4),rep(0,4),rep(0,4),rep(.1,4),0,.1,0,0,prior,prior_update="ebnm")
 expected<-do.call(multimodal_yfb_update_F_mk,args)
 record<-NULL;iter<-2L;k<-1L
 #' Capture one diagnostic row for this software fixture
 #' @param x Row from the observer.
 #' @return Assigned fixture record.
 #' @examples
 #' # save_row(data.frame(iteration=2))
 save_row<-function(x)record<<-x
 fn<-round2_prior_audit_wrapper(multimodal_yfb_update_F_mk,save_row)
 got<-do.call(fn,args)
 assert_true(identical(got,expected),msg="observational wrapper changes no update")
 assert_true(record$iteration==2L && record$factor==1L &&
   is.finite(record$stored_kl) && is.infinite(record$corrected_kl) &&
   record$corrected_kl>0,msg="known support mismatch and caller context preserved")
})
