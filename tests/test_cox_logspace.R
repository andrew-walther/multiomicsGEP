# ============================================================
# Script: test_cox_logspace.R
# Purpose: Verify equivalent Cox working arithmetic without global exponential underflow.
# Created: 2026-10-15
# Dependencies: existing test helpers and multimodal fitter
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/fit_multimodal_yfb.R")
if(file.exists("code/cox_working_logspace.R"))source("code/cox_working_logspace.R")
run_test("LogCox-T1: tied-time Breslow hand calculation is exact", {
 z<-multimodal_yfb_cox_working_logspace(c(1,1,2),c(1,1,1),c(0,0,0))
 assert_near(z$w,c(2/3,2/3,5/3),tol=1e-14)
 assert_near(z$u,c(1/3,1/3,-2/3),tol=1e-14)
 assert_near(z$log_partial_likelihood,-2*log(3),tol=1e-14)
})
run_test("LogCox-T2: extreme finite predictors retain valid late risk sets", {
 z<-multimodal_yfb_cox_working_logspace(c(1,2),c(1,1),c(1000,-1000))
 assert_near(z$w,c(1,1),tol=1e-12);assert_near(z$u,c(0,0),tol=1e-12)
 assert_near(z$log_partial_likelihood,0,tol=1e-12)
 z<-multimodal_yfb_cox_working_logspace(c(1,2),c(1,1),c(-1000,1000))
 assert_near(z$w,c(0,2),tol=1e-12);assert_near(z$log_partial_likelihood,-2000,tol=1e-12)
})
run_test("LogCox-T3: moderate inputs match original and retain shift identities", {
 t<-c(4,1,3,2,2);d<-c(1,1,0,1,0);eta<-c(.4,-.3,.7,1,-.4)
 a<-multimodal_yfb_cox_working(t,d,eta)
 b<-multimodal_yfb_cox_working_logspace(t,d,eta)
 for(nm in names(a))assert_near(b[[nm]],a[[nm]],tol=1e-12)
 c<-multimodal_yfb_cox_working_logspace(t,d,eta+1000)
 assert_near(c$w,b$w,tol=1e-12);assert_near(c$u,b$u,tol=1e-12)
 assert_near(c$h,b$h+1000*b$w,tol=1e-9)
 assert_near(c$log_partial_likelihood,b$log_partial_likelihood,tol=1e-12)
})
run_test("LogCox-T4: score equals numerical derivative and censored null is valid", {
 t<-c(4,1,3,2,2);d<-c(1,1,0,1,0);eta<-c(.4,-.3,.7,1,-.4);step<-1e-5
 z<-multimodal_yfb_cox_working_logspace(t,d,eta)
 for(i in seq_along(eta)) {
  plus<-minus<-eta;plus[i]<-plus[i]+step;minus[i]<-minus[i]-step
  grad<-(multimodal_yfb_cox_working_logspace(t,d,plus)$log_partial_likelihood-
    multimodal_yfb_cox_working_logspace(t,d,minus)$log_partial_likelihood)/(2*step)
  assert_near(grad,z$u[i],tol=1e-8)
 }
 z<-multimodal_yfb_cox_working_logspace(t,rep(0,5),eta)
 assert_near(z$w,rep(0,5),tol=0);assert_near(z$log_partial_likelihood,0,tol=0)
})
run_test("LogCox-T5: explicit original option leaves default full fit identical", {
 set.seed(2);Y<-list(expression=matrix(rnorm(60),15,4))
 dimnames(Y$expression)<-list(paste0("s",1:15),paste0("g",1:4))
 t<-setNames(1:15,rownames(Y$expression));d<-setNames(rep(1,15),names(t))
 ctrl<-list(prior_F=list(expression="normal"),intercept=TRUE,max_outer=2L,max_inner=1L)
 a<-fit_multimodal_yfb(Y,t,d,2L,ctrl)
 b<-fit_multimodal_yfb(Y,t,d,2L,c(ctrl,list(cox_arithmetic="original")))
 assert_true(identical(a,b))
 c<-fit_multimodal_yfb(Y,t,d,2L,c(ctrl,list(cox_arithmetic="logspace")))
 assert_near(a$EL,c$EL,tol=1e-9);assert_near(a$EBeta,c$EBeta,tol=1e-9)
 assert_true(identical(c$diagnostics$controls$cox_arithmetic,"logspace"))
})
run_test("LogCox-T6: the objective monitor uses the same opt-in arithmetic", {
 Y<-list(expression=matrix(c(1,-1),2,1));L<-matrix(1,2,1);F<-list(expression=matrix(1,1,1))
 z<-multimodal_yfb_elbo(Y,1:2,c(1,1),L,2*L,F,lapply(F,function(x)2*x),
  1000,1e6,0,cox_arithmetic="logspace")
 assert_true(is.finite(z));assert_near(attr(z,"components")['survival'],0,tol=1e-12)
})
run_test("LogCox-T7: default full fit remains identical to protected base-branch function", {
 source_text<-system2("git",c("show","feature/supervised-F:code/fit_multimodal_yfb.R"),stdout=TRUE)
 expr<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name("<-"))&&
   identical(e[[2]],as.name("fit_multimodal_yfb")),parse(text=source_text))
 assert_equal(length(expr),1L);base_fit<-eval(expr[[1]][[3]],.GlobalEnv)
 set.seed(2);Y<-list(expression=matrix(rnorm(60),15,4))
 dimnames(Y$expression)<-list(paste0("s",1:15),paste0("g",1:4))
 t<-setNames(1:15,rownames(Y$expression));d<-setNames(rep(1,15),names(t))
 ctrl<-list(prior_F=list(expression="normal"),intercept=TRUE,max_outer=2L,max_inner=1L)
 assert_true(identical(base_fit(Y,t,d,2L,ctrl),fit_multimodal_yfb(Y,t,d,2L,ctrl)))
})
