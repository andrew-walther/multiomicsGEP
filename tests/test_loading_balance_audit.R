# ============================================================
# Script: test_loading_balance_audit.R
# Purpose: Check isolated-coordinate likelihood contributions and their scale conventions.
# Created: 2026-10-15
# Dependencies: existing test helpers and base R
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("code/loading_balance_audit.R"))source("code/loading_balance_audit.R")
Y<-matrix(c(1,2,3,4),2,2,byrow=TRUE)
a<-list(Y_gen=Y,Y_surv=Y,Tau=c(2,3),EL=c(2,1),EL2=c(5,2),R=Y,
 EZ=c(1,2.5),EF=c(.5,.25),beta=.2,beta2=.09,w=c(.1,.2),h_minus=c(.3,-.4))
run_test("BalanceAudit-T1: contributions equal the explicit two-patient algebra", {
 z<-do.call(loading_balance_terms,a)
 assert_near(z$A_gen,c(14,21),tol=1e-14);assert_near(z$A_surv,c(.171,.324),tol=1e-14)
 assert_near(z$B_gen,c(10,24),tol=1e-14);assert_near(z$B_surv,c(-.2385,-.317),tol=1e-14)
})
run_test("BalanceAudit-T2: likelihood shares survive factor-scale changes", {
 z<-do.call(loading_balance_terms,a);b<-a;cscale<-3
 b$EL<-b$EL/cscale;b$EL2<-b$EL2/cscale^2;b$EF<-b$EF*cscale
 b$EZ<-b$EZ*cscale;b$beta<-b$beta/cscale;b$beta2<-b$beta2/cscale^2
 q<-do.call(loading_balance_terms,b)
 assert_near(q$A_gen,z$A_gen/cscale^2,tol=1e-14)
 assert_near(q$B_surv,z$B_surv/cscale,tol=1e-14)
 assert_near(q$precision_share,z$precision_share,tol=1e-14)
 assert_near(q$signal_share,z$signal_share,tol=1e-14)
})
run_test("BalanceAudit-T3: point-mass beta gates both terms, mean-zero beta alone does not", {
 b<-a;b$beta<-b$beta2<-0;z<-do.call(loading_balance_terms,b)
 assert_near(z$A_surv,c(0,0),tol=0);assert_near(z$B_surv,c(0,0),tol=0)
 b$beta2<-.09;z<-do.call(loading_balance_terms,b)
 assert_true(all(z$A_surv>0));assert_true(all(z$B_surv!=0))
})
run_test("BalanceAudit-T4: zero-information shares remain undefined", {
 b<-a;b$Tau<-c(0,0);b$Y_surv<-0*Y;z<-do.call(loading_balance_terms,b)
 assert_true(all(is.na(z$precision_share)));assert_true(all(is.na(z$signal_share)))
 b<-a;b$Y_surv<-0*Y;z<-do.call(loading_balance_terms,b)
 assert_near(z$precision_share,c(0,0),tol=0);assert_near(z$A_gen,c(14,21),tol=1e-14)
})
run_test("BalanceAudit-T5: first-coordinate terms reproduce the actual F-update pseudo-data", {
 source("code/fit_multimodal_yfb.R")
 z<-do.call(loading_balance_terms,a)
 u<-multimodal_yfb_update_F_mk(Y,a$Tau,a$EL,a$EL2,a$R,a$w,a$h_minus,
   a$EZ,rep(.1,2),a$EF,c(.3,.1),a$beta,a$beta2,
   list(pi=.5,rate=1),prior_update="em")
 assert_near(u$details$A[1],z$A_gen[1]+z$A_surv[1],tol=1e-13)
 assert_near(u$details$B[1],z$B_gen[1]+z$B_surv[1],tol=1e-13)
})
