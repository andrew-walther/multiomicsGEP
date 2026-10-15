# ============================================================
# Script: test_supervised_robustness.R
# Purpose: Verify ELBO-only start selection and exact supervised initialization
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: tests/test_helpers.R, ebnm, survival, robustness runner
# Run: Rscript tests/test_supervised_robustness.R
# ============================================================
source('tests/test_helpers.R')
source('results/benchmark_sim/run_supervised_robustness.R')
run_test('Best start is selected by ELBO, not recovery', {
 x<-data.frame(start=1:3,final_elbo=c(-4,-1,-2),recovered=c(TRUE,FALSE,TRUE))
 assert_equal(best_elbo_start(x)$start,2L)
 assert_false(best_elbo_start(x)$recovered)
 y<-x;y$final_elbo<-c(-1,-3,-2)
 assert_equal(best_elbo_start(y)$start,1L)
 assert_equal(best_elbo_start(x[c(3,1,2),])$start,2L)
 y$final_elbo[1]<-NA_real_
 err<-tryCatch(best_elbo_start(y),error=function(e)e)
 assert_true(inherits(err,'error'))
})
run_test('Custom supervised start matches the original score-direction construction', {
 tryCatch(source('code/fit_cox_on_yf.R'),error=function(e) {
  if(!grepl('!is.null(real_Y)',conditionMessage(e),fixed=TRUE))stop(e)
 })
 set.seed(22);Y<-scale(matrix(rnorm(60*12),60,12));tb<-list(Ytr=Y,ttr=rexp(60),str=rbinom(60,1,.7),seed=22L)
 zero<-make_jitter_start(tb,0L)
 u<-calc_cox_taylor_yf(rep(0,60),tb$ttr,tb$str)$u
 direction<-as.vector(crossprod(Y,u))
 assert_near(zero$EF_init[,7],direction,tol=1e-12)
 R<-Y-zero$EL_init[,-7]%*%t(zero$EF_init[,-7])
 assert_near(zero$EL_init[,7],as.vector(R%*%direction)/sum(direction^2),tol=1e-12)
 one<-make_jitter_start(tb,1L);again<-make_jitter_start(tb,1L)
 assert_true(identical(one,again))
 assert_false(identical(one,make_jitter_start(tb,2L)))
 # Empirical RMS relative to column SD must agree with the prescribed 0.1;
 # 420/84 fixed-seed draws give ample margin against a wrong jitter scale.
 rms_L<-sqrt(mean(sweep(one$EL_init-zero$EL_init,2,apply(zero$EL_init,2,sd),'/')^2))
 rms_F<-sqrt(mean(sweep(one$EF_init-zero$EF_init,2,apply(zero$EF_init,2,sd),'/')^2))
 assert_true(rms_L>.075 && rms_L<.125)
 assert_true(rms_F>.075 && rms_F<.125)
 assert_true(max(abs(one$EF_init-zero$EF_init))>0)
 assert_true(max(abs(one$EL_init-zero$EL_init))>0)
})
r<-report_results('robustness');if(r$failed)quit(status=1)
run_test('Unjittered custom start reproduces the default supervised initialization', {
 set.seed(22);Y<-scale(matrix(rnorm(60*12),60,12));tb<-list(Ytr=Y,ttr=rexp(60),str=rbinom(60,1,.7),seed=22L)
 a<-list(Y,tb$ttr,tb$str,K=7L,max_iter=8L,canonicalize=TRUE,w_S=1,
   prior_LF='point_laplace',init_supervised=TRUE,beta_prior_mode='shared',prior_beta='point_laplace',
   beta_warmup=10L,verbose=FALSE)
 def<-do.call(fit_cox_on_yf,a)
 cus<-do.call(fit_cox_on_yf,c(a,list(init_method='custom'),make_jitter_start(tb,0L)))
 assert_true(identical(def$EF,cus$EF)&&identical(def$EL,cus$EL)&&identical(def$EBeta,cus$EBeta))
})
r<-report_results('robustness with initialization equivalence');if(r$failed)quit(status=1)
