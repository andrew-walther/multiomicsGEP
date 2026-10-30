# ============================================================
# Script: diagnose_yfb_numerics.R
# Purpose: Independently integrate the saved synthetic point-exponential failure.
# Created: 2026-10-30
# Dependencies: base R; saved synthetic failing pseudo-message
# ============================================================
out<-"results/benchmark_sim/outputs/yfb_supervision_2026-10-30"
x<-readRDS(file.path(out,"point_exponential_failure_fixture.rds"))
j<-which.max(abs(x$q$second-x$z$EF2))
A<-1/x$z$s2[j];B<-x$z$x[j]*A
pi<-x$gp$pi[2];rate<-1/x$gp$scale[2]
# Cancel the common N(x;0,s²) term before integration, avoiding subtraction
# of enormous Gaussian log likelihoods. u=λθ keeps quadrature well scaled.
I<-vapply(0:2,function(power)stats::integrate(function(u)
  (u/rate)^power*exp(-u+B*u/rate-.5*A*(u/rate)^2),0,Inf,
  rel.tol=1e-11)$value,0.)
normalizer<-1-pi+pi*I[1]
m<-pi*I[2]/normalizer;m2<-pi*I[3]/normalizer
kl<-B*m-.5*A*m2-log(normalizer)
result<-data.frame(sweep=x$iter,factor=x$k,x=x$z$x[j],se=sqrt(x$z$s2[j]),
  cpp_mean=x$z$EF[j],cpp_second=x$z$EF2[j],ebnm_mean=x$q$mean[j],ebnm_second=x$q$second[j],
  quadrature_mean=m,quadrature_second=m2,quadrature_variance=m2-m^2,quadrature_KL=kl)
write.csv(result,file.path(out,"numerical_failure.csv"),row.names=FALSE)
stopifnot(abs(result$cpp_second-m2)/m2<1e-3,
          result$quadrature_variance>0,result$ebnm_second<.5*m2)
print(result)
