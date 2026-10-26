# ============================================================
# Script: test_framework_fold_trace.R
# Purpose: Verify unchanged framework settings and detection of coefficient extinction/revival.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing test helpers and diagnostic runner
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("results/multimodal_real/run_framework_fold_trace.R"))
  source("results/multimodal_real/run_framework_fold_trace.R")
run_test("FrameworkTrace-T1: tracing changes only the observational flag", {
  a <- framework_trace_controls();a$trace <- FALSE
  assert_true(identical(a,fixed_prior_controls("multimodal","framework")),
    msg="same seed-independent prior, initialization, pruning and convergence controls")
})
run_test("FrameworkTrace-T2: all-zero boundary distinguishes transient loss from persistent loss", {
  x <- data.frame(stage_outer=rep(1:4,each=2),inner=1,
    factor=rep(1:2,4),beta=c(1,0,0,0,1,0,0,0),
    beta_sd=c(.1,0,0,0,.1,0,0,0),x=c(2,.1,.2,.1,2,.1,.2,.1),A=1)
  a <- framework_trace_boundaries(x)
  assert_true(identical(a$all_zero,c(FALSE,TRUE,FALSE,TRUE)) &&
    identical(a$all_pseudo_subthreshold,c(FALSE,TRUE,FALSE,TRUE)),
    msg="revival retained and per-sweep beta moments checked")
  assert_true(all(a$persistent_zero_start==4L),msg="transient zero is not permanent extinction")
})
run_test("FrameworkTrace-T3: zero means with positive uncertainty are not degenerate", {
  x <- data.frame(stage_outer=rep(1:2,each=2),inner=1,factor=rep(1:2,2),
    beta=0,beta_sd=c(.1,.1,0,0),x=0,A=1)
  a <- framework_trace_boundaries(x)
  assert_true(identical(a$all_zero,c(FALSE,TRUE)) &&
    all(a$persistent_zero_start==2L),msg="requires both means and SDs zero")
})
