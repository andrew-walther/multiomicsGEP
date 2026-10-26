# ============================================================
# Script: test_multimodal_yfb_trace.R
# Purpose: Verify opt-in collapse diagnostics and the specified diagnostic priors.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing test helpers, ebnm, ashr, multimodal YFB modules
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("code/fit_multimodal_yfb.R")
source("code/simulate_multimodal_yfb.R")
.trace_d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
  n_train=60L,n_validation=20L,p_expression=12L,p_methylation=12L,
  seed=7L,truncate=FALSE,noise_scale=.5,beta=c(0,0,3.3))$training
.trace_ctrl <- list(intercept=TRUE,prior_update="ebnm",prune=FALSE,
  prior_F=list(expression="point_laplace",methylation="point_laplace"),
  init_supervised=TRUE,beta_prior="shared_point_laplace",beta_warmup=1L,
  max_outer=4L,max_inner=2L)

run_test("MMTrace-T1: tracing observes every coefficient and leaves the fit identical", {
  a <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,.trace_ctrl)
  b <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,
    c(.trace_ctrl,list(trace=TRUE)))
  tr <- b$diagnostics$trace
  assert_true(length(tr)==b$diagnostics$iterations,msg="one trace per outer iteration")
  for(i in seq_along(tr)) {
    x <- tr[[i]]
    assert_true(length(x$beta_mean)==3L && length(x$beta_sd)==3L,
      msg="all coefficient means and SDs recorded")
    assert_true(length(x$x)==3L && all(is.finite(x$x)) && all(x$A>=0),
      msg="finite last-sweep pseudo-observations and nonnegative precisions")
    assert_near(x$elbo_components["reconstruction"]+x$elbo_components["survival"]-
      x$elbo_components["kl"],b$diagnostics$history[[i]]$elbo,tol=1e-8,
      msg="recorded components reconstruct the actual monitor")
    assert_true(length(x$prior_sweeps)==2L,msg="prior refits observable within each outer iteration")
  }
  assert_near(tail(tr,1)[[1]]$beta_mean,b$EBeta,tol=0,msg="trace is the returned posterior")
  b$diagnostics$trace <- NULL; b$diagnostics$controls$trace <- NULL
  assert_true(identical(a,b),msg="trace changes no fitted quantity or convergence flag")
  off <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,
    c(.trace_ctrl,list(trace=FALSE,beta_prior_fix=FALSE)))
  assert_true(identical(a,off),msg="explicitly disabled controls preserve the default object")
})

run_test("MMTrace-T2: fixed shared Laplace prior never refits beyond warmup", {
  a <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,
    c(.trace_ctrl,list(trace=TRUE,beta_prior_fix=TRUE)))
  g0 <- ebnm::laplacemix(c(.5,.5),c(0,0),c(0,1))
  assert_true(identical(a$g_beta,g0),msg="exact starting prior held after warmup")
  assert_true(all(vapply(a$diagnostics$trace,function(x)
    identical(x$prior_after,g0),logical(1))),msg="prior held in every recorded sweep")
})

run_test("MMTrace-T3: shared Normal update matches its conjugate posterior", {
  z <- c(-2,-1,1,2); v <- rep(.1,4); w <- rep(.8,4); h <- .7*z
  g <- ashr::normalmix(1,0,1)
  a <- multimodal_yfb_update_beta_shared_k(z,v,w,h,g,prior_family="normal")
  A <- sum(w*(z^2+v)); B <- sum(z*h)
  assert_near(a$mean,B/(A+1),tol=1e-10,msg="Normal posterior mean")
  assert_near(a$variance,1/(A+1),tol=1e-10,msg="Normal posterior variance")
  ctrl <- modifyList(.trace_ctrl,list(beta_prior="shared_normal",beta_warmup=2L,
    trace=TRUE,max_outer=3L))
  f <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,ctrl)
  assert_true(identical(f$diagnostics$trace[[2]]$prior_after,g),msg="Normal warmup held")
  assert_true(!identical(f$g_beta,g),msg="Normal prior refit after warmup")
})

run_test("MMTrace-T4: tracing preserves every pruning stage without changing the result", {
  ctrl <- modifyList(.trace_ctrl,list(prune=TRUE,max_outer=4L))
  a <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,5L,ctrl)
  b <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,5L,
    c(ctrl,list(trace=TRUE)))
  tr <- b$diagnostics$trace
  assert_true(nrow(b$diagnostics$pruning)>0,msg="fixture actually exercises a pruning refit")
  assert_true(length(tr)>b$diagnostics$iterations && tr[[1]]$K==5L &&
    tail(tr,1)[[1]]$K==ncol(b$EL),msg="all pruning stages retained in order")
  assert_near(tail(tr,1)[[1]]$live_programs,
    sum(colSums(do.call(rbind,b$EF)^2)>1e-8),tol=0,msg="live means nonzero loadings, including survival-only programs")
  b$diagnostics$trace <- NULL;b$diagnostics$controls$trace <- NULL
  assert_true(identical(a,b),msg="pruned fit identical with trace removed")
})
run_test("MMTrace-T5: intercept and centered reconstruction changes match consecutive fits", {
  ctrl <- modifyList(.trace_ctrl,list(max_outer=1L))
  a <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,ctrl)
  b <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,
    modifyList(ctrl,list(max_outer=2L,trace=TRUE)))
  x <- b$diagnostics$trace[[2]]
  assert_true(length(x$full_reconstruction_change)==2L &&
    length(x$centered_reconstruction_change)==2L,msg="both observed modality changes recorded")
  for(m in names(a$EF)) {
    old <- a$EL%*%t(a$EF[[m]]);new <- b$EL%*%t(b$EF[[m]])
    full <- multimodal_yfb_relative_change(sweep(old,2,a$mu[[m]],"+"),
      sweep(new,2,b$mu[[m]],"+"))
    centered <- multimodal_yfb_relative_change(scale(old,scale=FALSE),scale(new,scale=FALSE))
    assert_near(x$full_reconstruction_change[m],full,tol=1e-12,msg="full molecular prediction includes intercept")
    assert_near(x$centered_reconstruction_change[m],centered,tol=1e-12,msg="centered diagnostic removes intercept allocation")
  }
})
run_test("MMTrace-T6: diagnostics also support fits without an intercept", {
  ctrl <- modifyList(.trace_ctrl,list(intercept=FALSE,trace=TRUE,max_outer=2L))
  f <- fit_multimodal_yfb(.trace_d$Y,.trace_d$time,.trace_d$event,3L,ctrl)
  for(x in f$diagnostics$trace)
    assert_near(x$full_reconstruction_change,x$reconstruction_change,tol=0,
      msg="with zero intercept the two molecular changes coincide")
})
