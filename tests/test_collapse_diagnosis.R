# ============================================================
# Script: test_collapse_diagnosis.R
# Purpose: Check the fixed diagnosis arms and faithful trace export.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing test helpers and diagnosis runner definitions
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/multimodal_real/run_round2_collapse_diagnosis.R")
run_test("Collapse-T1: diagnostic arms change exactly their stated option", {
  for(m in c("single","multimodal")) {
    a <- collapse_arm_controls(m,"selected")
    b <- collapse_arm_controls(m,"fixed")
    b$beta_prior_fix <- NULL
    assert_true(identical(a,b),msg="fixed arm changes only whether prior refits")
    b <- collapse_arm_controls(m,"no_start")
    b$init_supervised <- TRUE
    assert_true(identical(a,b),msg="no-start arm changes only supervised initialization")
    b <- collapse_arm_controls(m,"normal")
    if(m=="single")b$prior_beta <- "point_laplace" else b$beta_prior <- "shared_point_laplace"
    assert_true(identical(a,b),msg="Normal arm changes only the shared prior family")
  }
})
run_test("Collapse-T2: exported trace retains prior collapse and raw precision", {
  g <- ebnm::laplacemix(c(.5,.5),c(0,0),c(0,1)); z <- ebnm::laplacemix(1,0,0)
  sweep <- list(prior_before=g,prior_after=z,beta_mean=c(.2,0),beta_sd=c(.05,0),
    x=c(.3,-.1),A=c(10,20),canonical_scale=c(2,3))
  tr <- list(list(K=2L,iteration=11L,prior_sweeps=list(sweep),
    live_programs=1,elbo_components=c(reconstruction=-100,survival=-20,kl=3),
    reconstruction_change=c(expression=.1,methylation=.2),
    full_reconstruction_change=c(expression=.01,methylation=.02),
    centered_reconstruction_change=c(expression=.001,methylation=.002),eta_change=.03,
    centered_eta_change=.01,eta_mean=7,eta_sd=2))
  d <- collapse_trace_rows(tr)
  assert_true(nrow(d)==2L && all(d$prior_pi_before==.5) && all(d$prior_pi_after==0),
    msg="point mass transition is preserved")
  assert_near(d$full_reconstruction_change,c(.02,.02),tol=0,
    msg="intercept-inclusive molecular change retained")
  assert_near(d$centered_reconstruction_change,c(.002,.002),tol=0,
    msg="centered molecular diagnostic retained separately")
  assert_near(d$A,c(10,20),tol=0,msg="precisions are actual update values")
  assert_near(d$canonical_scale,c(2,3),tol=0,msg="canonical rescaling recorded separately")
  assert_true(all(d$outer==11L & d$K==2L & d$live_programs==1),msg="stage identifiers retained")
})

run_test("Collapse-T3: a failed data load is a checkpointed arm failure", {
  p <- tempfile("collapse-error-");dir.create(p)
  x <- run_collapse_cell(1L,p,training_path=file.path(p,"missing_training.rds"))
  assert_true(nrow(x)==1L && !is.na(x$error),msg="failure is a row, not a discarded arm")
  saved <- read.csv(file.path(p,"single_selected_status.csv"))
  assert_true(identical(saved$error,x$error) && saved$arm=="selected",
    msg="exact error and arm checkpointed")
  unlink(p,recursive=TRUE)
})
run_test("Collapse-T4: an invalid diagnosis cell stops before fitting", {
  p <- tempfile("collapse-invalid-")
  e <- tryCatch(run_round2_collapse_diagnosis("11",p),error=function(e)e)
  assert_true(inherits(e,"error"),msg="invalid array cell rejected")
  assert_true(length(list.files(p))==0L,msg="no fit or result fabricated")
  unlink(p,recursive=TRUE)
})

run_test("Collapse-T5: Single prior trace records the actual refit and changes no fit", {
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"),error=function(e) {
    if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)
  }))
  source("code/simulate_multimodal_yfb.R")
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train=60L,n_validation=20L,p_expression=12L,p_methylation=12L,
    seed=7L,truncate=FALSE,noise_scale=.5,beta=c(0,0,3.3))$training
  a <- collapse_arm_controls("single","selected")
  a$K <- 3L;a$max_iter <- 4L;a$beta_warmup <- 1L;a$tol <- 0
  b <- a;b$trace_supervision <- FALSE
  f <- do.call(fit_cox_on_yf,c(list(scale(d$Y$expression),d$time,d$event),a))
  off <- do.call(fit_cox_on_yf,c(list(scale(d$Y$expression),d$time,d$event),b))
  tr <- f$history$prior_trace
  assert_true(length(tr)==4L,msg="one actual shared-prior record per iteration")
  for(i in seq_along(tr)) {
    z <- tr[[i]]
    assert_near(z$A,1/z$s^2,tol=1e-10,msg="actual coefficient precision and noise agree")
    if(i>1L) {
      expected <- ebnm::ebnm(z$x,z$s,prior_family="point_laplace")$fitted_g
      assert_true(identical(z$prior_after,expected),msg="recorded prior is the actual refit")
      assert_true(identical(z$prior_before,tr[[i-1L]]$prior_after),msg="recorded prior is used by the next sweep")
    }
  }
  f$history$trace <- NULL;f$history$prior_trace <- NULL
  assert_true(identical(f,off),msg="prior diagnostics change no fitted quantity")
})
