# ============================================================
# Script: test_multimodal_moment_audit.R
# Purpose: Verify finite-moment diagnostics without changing fitting behavior.
# Created: 2026-10-15
# Dependencies: existing test helpers and multimodal updates
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("code/multimodal_moment_audit.R"))source("code/multimodal_moment_audit.R")
run_test("MomentAudit-T1: legitimate signed working means are accepted", {
 z<-multimodal_moment_guard_summary(c(0,2),c(-4,3),c(-2,1),c(0,.2))
 assert_equal(z$invalid,c(FALSE,FALSE,FALSE,FALSE))
 assert_near(z$minimum,c(0,-4,-2,0),tol=0)
 assert_equal(z$negative,c(0L,1L,1L,0L))
})
run_test("MomentAudit-T2: nonfinite and invalid variance values identify their field", {
 z<-multimodal_moment_guard_summary(c(0,Inf),c(NaN,3),c(NA,1),c(-1e-6,.2))
 assert_equal(z$invalid,rep(TRUE,4));assert_equal(z$nonfinite,c(1L,1L,1L,0L))
 assert_near(z$minimum,c(0,3,1,-1e-6),tol=0)
 z<-multimodal_moment_guard_summary(-1,1,2,-1e-11)
 assert_equal(z$invalid,c(TRUE,FALSE,FALSE,FALSE))
})
run_test("MomentAudit-T3: an entry observer leaves the exact F-update output unchanged", {
 source("code/fit_multimodal_yfb.R")
 a<-list(Y_m=matrix(c(1,2,3,2,3,1),3,2),Tau_m=c(1,2),
 EL_k=rep(1,3),EL2_k=rep(2,3),R_mk=matrix(1,3,2),w=rep(.1,3),
 h_minus_k=c(-.1,.2,.1),EZ_k=rep(.2,3),VZ_k=rep(.4,3),
 EF_k=rep(.1,2),EF2_k=rep(.2,2),EBeta_k=.1,EBeta2_k=.02,
 prior=list(pi=.5,rate=1),prior_update="em")
 old<-do.call(multimodal_yfb_update_F_mk,a)
 trace("multimodal_yfb_update_F_mk",tracer=quote({
   z_audit<-multimodal_moment_guard_summary(w,h_minus_k,EZ_k,VZ_k)
   stopifnot(!any(z_audit$invalid))
 }),print=FALSE)
 new<-tryCatch(do.call(multimodal_yfb_update_F_mk,a),finally=untrace("multimodal_yfb_update_F_mk"))
 assert_true(identical(old,new))
})
