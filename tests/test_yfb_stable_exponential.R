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
