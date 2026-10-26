# ============================================================
# Script: test_single_pruning_trace.R
# Purpose: Verify pruning instrumentation preserves results and original bindings.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing test helpers and single-modality fitter
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("results/benchmark_sim/run_single_pruning_trace.R")
run_test("SinglePruneTrace-T1: instrumentation leaves pruning and errors unchanged", {
  if(!exists("fit_cox_on_yf"))suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"),
    error=function(e) {if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)}))
  set.seed(3);y <- matrix(rnorm(50*10),50,10)
  a <- list(K=3L,max_iter=6L,tol=1e-5,canonicalize=TRUE,w_S=1,
    prior_LF="point_laplace",init_supervised=TRUE,beta_prior_mode="shared",
    prior_beta="point_laplace",beta_warmup=10L,prune=TRUE,
    prune_protect_survival=TRUE,verbose=FALSE)
  old <- prune_yf_nullcheck
  ref <- do.call(fit_cox_on_yf,c(list(y,1:50,rep(1,50)),a))
  tr <- trace_single_pruning(y,1:50,rep(1,50),a)
  assert_true(identical(ref,tr$fit),msg="every posterior and pruning decision unchanged")
  assert_true(length(tr$stages)>0L,msg="actual pruning nullchecks captured")
  for(s in tr$stages) {
    gain <- old(s$fit,y,1:50,rep(1,50),protect_survival=TRUE)
    assert_true(identical(gain,s$delta),msg="exact unmodified nullcheck gains recorded")
  }
  assert_true(identical(prune_yf_nullcheck,old),msg="binding restored after success")
  bad <- a;bad$K <- -1L
  e <- tryCatch(trace_single_pruning(y,1:50,rep(1,50),bad),error=function(e)e)
  assert_true(inherits(e,"error") && identical(prune_yf_nullcheck,old),
    msg="error surfaced and original binding restored")
})
run_test("SinglePruneTrace-T2: diagnosis refuses to overwrite an existing attempt", {
  p <- tempfile("pruning-observer-");dir.create(p)
  writeLines("keep",file.path(p,"checkpoint"))
  e <- tryCatch(run_single_pruning_trace(p),error=function(e)e)
  assert_true(inherits(e,"error"),msg="existing attempt refused before fitting")
  assert_true(identical(readLines(file.path(p,"checkpoint")),"keep"),
    msg="previous output remains intact")
  unlink(p,recursive=TRUE)
})
