# Exact finite-pattern reference for the Gaussian working likelihood. ----
source("code/fit_yfb_supervised.R")
if(file.exists("code/yfb_joint_beta.R"))source("code/yfb_joint_beta.R")
run_test("YFB joint beta: orthogonal spike/slab posterior matches scalar formulas", local({
  X<-diag(c(1,2,3));y<-c(-1,.4,2);pi0<-.25;sigma<-1.2
  q<-yfb_joint_beta(X,y,weight=2,prior_pi=pi0,slab_sd=sigma)
  a<-2*colSums(X^2);b<-2*drop(crossprod(X,y))
  scalar<-yfb_normal_posterior(b/a,a,ashr::normalmix(c(1-pi0,pi0),c(0,0),c(0,sigma)))
  assert_near(q$mean,scalar$mean,1e-10)
  assert_near(q$variance,scalar$variance,1e-10)
  assert_near(sum(q$weights),1,1e-12)
  mf<-yfb_meanfield_beta(X,y,2,pi0,sigma)
  assert_true(mf$converged);assert_near(mf$mean,q$mean,1e-8)
  assert_near(mf$variance,q$variance,1e-8)
}))
run_test("YFB joint beta: correlated Gaussian posterior has the exact covariance", local({
  X<-cbind(c(1,2,-1),c(1,1,-2));y<-c(.2,-1,.7)
  V<-solve(3*crossprod(X)+diag(2));m<-V%*%(3*crossprod(X,y))
  q<-yfb_joint_beta(X,y,3,1,1)
  assert_near(q$mean,drop(m),1e-10);assert_near(q$covariance,V,1e-10)
  mf<-yfb_meanfield_beta(X,y,3,1,1)
  assert_true(mf$converged);assert_near(mf$mean,q$mean,1e-7)
  assert_true(all(q$variance>mf$variance))
}))
