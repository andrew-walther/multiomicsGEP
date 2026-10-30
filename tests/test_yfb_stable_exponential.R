# Independent quadrature checks for the diagnosed exponential posterior. ----
source("code/yfb_exponential_stable.R")
run_test("Stable exponential: weak message retains the quadrature variance", local({
  msg<-readRDS("results/benchmark_sim/outputs/yfb_supervision_2026-10-30/point_exponential_failure_fixture.rds")
  j<-which.max(abs(msg$q$second-msg$z$EF2))
  A<-1/msg$z$s2[j];B<-msg$z$x[j]*A
  q<-yfb_stable_exp(A,B,msg$gp$pi[2],1/msg$gp$scale[2])
  exact<-read.csv("results/benchmark_sim/outputs/yfb_supervision_2026-10-30/numerical_failure.csv")
  assert_near(q$mean,exact$quadrature_mean,1e-9)
  assert_near(q$second,exact$quadrature_second,1e-9)
  assert_near(q$kl,exact$quadrature_KL,1e-10)
}))
run_test("Stable exponential: moderate messages agree with independent integrals", local({
  for(B in c(-10,0,2,10)) {
    A<-3;pi<-.3;rate<-2
    I<-vapply(0:2,function(k)integrate(function(t)
      t^k*rate*exp((B-rate)*t-.5*A*t^2),0,Inf,rel.tol=1e-10)$value,0.)
    norm<-1-pi+pi*I[1];m<-pi*I[2]/norm;m2<-pi*I[3]/norm
    q<-yfb_stable_exp(A,B,pi,rate)
    assert_near(q$mean,m,1e-8);assert_near(q$second,m2,1e-8)
    assert_near(q$kl,B*m-.5*A*m2-log(norm),1e-8)
  }
}))

run_test("Stable exponential tests preserve the mathematical pi constant", {
  assert_near(pi,acos(-1),1e-15)
})

run_test("Stable coupled sweep: correlated genes agree with conditional quadrature", local({
  Y<-cbind(c(-1,0,1),c(-1,1,0));w<-c(.2,.3,.4);hm<-c(.5,-.4,.8)
  tau<-c(1.5,.7);LR<-c(.7,.3);sumL2<-2;beta<-.8;beta2<-.9
  oldF<-c(.2,.4);oldV<-c(.03,.04);means<-oldF
  variances<-KL<-numeric(2);prior_pi<-.6;rate<-1.3
  # Integrate each conditional directly, using the other gene's current
  # mean. This independent formula does not use the C++ projection refresh.
  for(j in 1:2) {
    other<-3-j
    A<-tau[j]*sumL2+beta2*sum(w*Y[,j]^2)
    B<-tau[j]*LR[j]+beta*sum(Y[,j]*hm)-
      beta2*sum(w*Y[,j]*Y[,other])*means[other]
    integrals<-vapply(0:2,function(k)integrate(function(t)
      t^k*rate*exp((B-rate)*t-.5*A*t^2),0,Inf,rel.tol=1e-11)$value,0.)
    norm<-1-prior_pi+prior_pi*integrals[1]
    means[j]<-prior_pi*integrals[2]/norm
    second<-prior_pi*integrals[3]/norm
    variances[j]<-second-means[j]^2
    KL[j]<-B*means[j]-.5*A*second-log(norm)
  }
  yfb_stable_compile()
  q<-yfb_exp_sweep_cpp(Y,tau,sumL2,colSums(Y^2*w),LR,
    drop(crossprod(Y,hm)),w,drop(Y%*%oldF),drop(Y^2%*%oldV),
    oldF,oldF^2+oldV,oldV,beta,beta2,prior_pi,rate)
  assert_near(q$EF,means,1e-9);assert_near(q$VF,variances,1e-9)
  assert_near(q$kl,KL,1e-9)
  assert_near(q$EZ,drop(Y%*%means),1e-12)
  assert_near(q$VZ,drop(Y^2%*%variances),1e-12)
}))
