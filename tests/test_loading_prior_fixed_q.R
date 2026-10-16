# ============================================================
# Script: test_loading_prior_fixed_q.R
# Purpose: Verify fixed-posterior point-Laplace prior maximization and KL.
# Created: 2026-10-16
# Dependencies: base R, existing test helpers and molecular updates
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("code/multimodal_yfb_updates.R")
if (file.exists("code/loading_prior_fixed_q.R")) source("code/loading_prior_fixed_q.R")

run_test("FixedQ-T1: prior maximizer uses unconditional absolute moments", {
  q <- list(prob = c(.2, .7, 1), mean_abs = c(.1, .8, 2))
  g <- loading_prior_laplace_em(q$prob, q$mean_abs)
  assert_near(g$pi, 1.9 / 3, tol = 1e-15)
  assert_near(g$rate, 1.9 / 2.9, tol = 1e-15)
  best <- loading_prior_laplace_log(q$prob, q$mean_abs, g)
  for (pi_value in c(.01, .2, .5, .8, .99)) for (rate in c(.1, .5, 1, 4)) {
    other <- list(family = "point_laplace", pi = pi_value, rate = rate, point_mass = FALSE)
    assert_true(best >= loading_prior_laplace_log(q$prob, q$mean_abs, other))
  }
})
run_test("FixedQ-T2: boundary support is explicit without floors", {
  dead <- loading_prior_laplace_em(c(0, 0), c(0, 0))
  assert_true(dead$point_mass && dead$pi == 0 && is.na(dead$rate))
  pure <- loading_prior_laplace_em(c(1, 1), c(.5, 1))
  assert_near(pure$pi, 1, tol = 0)
  assert_near(loading_prior_laplace_log(c(0, 0), c(0, 0), dead), 0, tol = 0)
  assert_equal(loading_prior_laplace_log(.1, .2, dead), -Inf)
  assert_equal(loading_prior_laplace_log(.9, .2, pure), -Inf)
  for (q in list(list(c(0), c(.1)), list(c(.5), c(0)), list(c(-.1), c(.2)),
                 list(numeric(), numeric()), list(c(.5), c(NA_real_)))) {
    err <- tryCatch(loading_prior_laplace_em(q[[1]], q[[2]]), error = identity)
    assert_true(inherits(err, "error"), msg = "invalid moments fail loudly")
  }
})
run_test("FixedQ-T3: mixed absolute moment agrees with independent quadrature", {
  g <- list(family = "point_laplace", pi = .4, rate = 1.3, point_mass = FALSE)
  for (B in c(-5, 0, 2)) {
    q <- loading_prior_laplace_stats(c(4, 0), c(B, 0), g)
    p <- multimodal_yfb_signed_posterior(4, B, g)
    density <- function(f) g$pi * g$rate / 2 * exp(-g$rate * abs(f)) *
      dnorm(B/4 - f, 0, .5) / exp(p$log_ml)
    abs_direct <- integrate(function(f) abs(f)*density(f), -Inf, 0)$value +
      integrate(function(f) abs(f)*density(f), 0, Inf)$value
    assert_near(q$prob[1], p$slab_prob, tol = 1e-13)
    assert_near(q$mean_abs[1], abs_direct, tol = 1e-8)
    assert_near(q$mean_abs[2], g$pi/g$rate, tol = 0)
  }
})
run_test("FixedQ-T4: corrected KL equals direct mixed-measure quadrature", {
  old <- list(family = "point_laplace", pi = .4, rate = 1.3, point_mass = FALSE)
  p <- multimodal_yfb_signed_posterior(4, 2, old)
  d <- lapply(c(p, list(A = 4, B = 2)), as.numeric)
  out <- loading_prior_fixed_q(d, old)
  new <- out$prior
  mass <- (1-old$pi)*dnorm(.5, 0, .5)/exp(p$log_ml)
  density <- function(f) old$pi*old$rate/2*exp(-old$rate*abs(f))*
    dnorm(.5-f, 0, .5)/exp(p$log_ml)
  integrand <- function(f) {
    z <- density(f)
    ratio <- log(z) - log(new$pi*new$rate/2) + new$rate*abs(f)
    ifelse(z > 0, z*ratio, 0)
  }
  direct <- mass*log(mass/(1-new$pi)) + integrate(integrand, -Inf, 0)$value +
    integrate(integrand, 0, Inf)$value
  assert_near(out$kl, direct, tol = 1e-8)
  assert_true(out$prior_gain >= 0)
  assert_near(out$kl, multimodal_yfb_coordinate_kl(p$x, p$s2, p$mean,
    p$second, p$log_ml)-out$prior_gain, tol = 1e-12)
})
run_test("FixedQ-T5: no information and point masses keep zero KL", {
  for (g in list(list(family="point_laplace", pi=.3, rate=2, point_mass=FALSE),
                 list(family="point_laplace", pi=0, rate=NA_real_, point_mass=TRUE))) {
    p <- multimodal_yfb_signed_posterior(0, 0, g)
    out <- loading_prior_fixed_q(c(p, list(A=0, B=0)), g)
    assert_near(out$kl, 0, tol=1e-14)
    assert_near(out$prior_gain, 0, tol=1e-14)
  }
})
run_test("FixedQ-T6: scale transformation preserves gain and KL", {
  g <- list(family="point_laplace", pi=.7, rate=.8, point_mass=FALSE)
  p <- multimodal_yfb_signed_posterior(3, -2, g)
  out <- loading_prior_fixed_q(c(p, list(A=3, B=-2)), g)
  c <- 5; scaled <- g; scaled$rate <- g$rate/c
  ps <- multimodal_yfb_signed_posterior(3/c^2, -2/c, scaled)
  other <- loading_prior_fixed_q(c(ps, list(A=3/c^2, B=-2/c)), scaled)
  assert_near(out$prior_gain, other$prior_gain, tol=1e-12)
  assert_near(out$kl, other$kl, tol=1e-12)
  assert_near(out$prior$rate/c, other$prior$rate, tol=1e-12)
})
run_test("FixedQ-T7: tiny slab mass is retained and extreme tails stay finite", {
  for (g in list(list(family="point_laplace", pi=1e-20, rate=1, point_mass=FALSE),
                 list(family="point_laplace", pi=1, rate=1, point_mass=FALSE))) {
    q <- loading_prior_laplace_stats(c(.001, 4, 4), c(0, -100, 100), g)
    assert_finite(c(q$prob, q$mean_abs))
    assert_true(all(q$prob > 0 & q$mean_abs > 0))
    assert_near(q$mean_abs[2], q$mean_abs[3], tol=1e-12)
    assert_true(all(q$prob <= 1))
  }
})
run_test("FixedQ-T8: informative missing normalizers cannot be silently skipped", {
  g <- list(family="point_laplace", pi=.4, rate=1, point_mass=FALSE)
  p <- multimodal_yfb_signed_posterior(4, 2, g); p$log_ml <- NA_real_
  err <- tryCatch(loading_prior_fixed_q(c(p, list(A=4, B=2)), g), error=identity)
  assert_true(inherits(err, "error") && grepl("actual finite normalizer", conditionMessage(err)))
  for (field in c("x", "s2", "mean", "second")) {
    d <- c(multimodal_yfb_signed_posterior(4,2,g),list(A=4,B=2)); d[[field]] <- NA_real_
    assert_true(inherits(tryCatch(loading_prior_fixed_q(d,g),error=identity),"error"))
  }
})
run_test("FixedQ-T9: opt-in changes only prior and KL after the same sequential sweep", {
  args <- list(Y_m=cbind(c(1,2), c(3,4)), Tau_m=c(2,3), EL_k=c(1,2),
    EL2_k=c(1.2,4.5), R_mk=cbind(c(5,6),c(7,8)), w=c(.5,1),
    h_minus_k=c(2,1), EZ_k=c(7,10), VZ_k=c(.3,.4), EF_k=c(1,2),
    EF2_k=c(1.1,4.2), EBeta_k=.5, EBeta2_k=.8,
    prior=list(family="point_laplace", pi=.6, rate=1.2, point_mass=FALSE))
  base <- do.call(multimodal_yfb_update_F_mk, args)
  off <- do.call(multimodal_yfb_update_F_mk, c(args,list(loading_prior_update="historical")))
  on <- do.call(multimodal_yfb_update_F_mk, c(args,list(loading_prior_update="fixed_q_em")))
  assert_equal(base, off)
  for (nm in c("EF","EF2","EZ","VZ","details")) assert_equal(base[[nm]], on[[nm]])
  q <- loading_prior_laplace_stats(on$details$A, on$details$B, args$prior)
  assert_equal(on$prior, loading_prior_laplace_em(q$prob,q$mean_abs))
  expected <- loading_prior_fixed_q(on$details,args$prior)
  assert_near(on$kl,expected$kl,tol=1e-12)
})
run_test("FixedQ-T10: fitter controls reject unsupported prior and damping combinations", {
  source("code/fit_multimodal_yfb.R")
  tryCatch(source("code/fit_cox_on_yf.R"), error=function(e) {
    if (!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE)) stop(e)
  })
  set.seed(46); Y <- scale(matrix(rnorm(120),20,6)); dimnames(Y) <- list(paste0("s",1:20),paste0("g",1:6))
  tm <- setNames(1:20,rownames(Y)); ev <- setNames(rep(c(1,0),10),rownames(Y))
  bad <- tryCatch(fit_cox_on_yf(Y,tm,ev,loading_prior_update="fixed_q_em",verbose=FALSE),error=identity)
  assert_true(inherits(bad,"error") && grepl("fixed_q_em",conditionMessage(bad)))
  for (extra in list(list(damping=.5), list(prune=TRUE), list(prior_F=list(expression="normal")))) {
    ctl <- modifyList(list(loading_prior_update="fixed_q_em", max_outer=1L,
      prior_F=list(expression="point_laplace"), prune=FALSE),extra)
    bad <- tryCatch(fit_multimodal_yfb(list(expression=Y),tm,ev,2L,ctl),error=identity)
    assert_true(inherits(bad,"error") && grepl("fixed_q_em",conditionMessage(bad)))
  }
})
run_test("FixedQ-T11: actual Single and multimodal callers forward the opt-in only", {
  set.seed(46); Y <- scale(matrix(rnorm(120),20,6)); dimnames(Y) <- list(paste0("s",1:20),paste0("g",1:6))
  tm <- setNames(1:20,rownames(Y)); ev <- setNames(rep(c(1,0),10),rownames(Y))
  original <- multimodal_yfb_update_F_mk; seen <- character()
  multimodal_yfb_update_F_mk <- function(...) {
    a <- list(...); seen <<- c(seen, a$loading_prior_update); original(...)
  }
  tryCatch({
    fit_cox_on_yf(Y,tm,ev,K=2L,max_iter=2L,canonicalize=TRUE,prior_LF="point_laplace",
      w_S=1,loading_prior_update="fixed_q_em",verbose=FALSE)
    ns <- length(seen); assert_true(ns > 0 && all(seen=="fixed_q_em"))
    fit_multimodal_yfb(list(expression=Y),tm,ev,2L,list(max_outer=2L,max_inner=1L,
      prune=FALSE, prior_F=list(expression="point_laplace"),loading_prior_update="fixed_q_em"))
    assert_true(length(seen)>ns && all(seen=="fixed_q_em"))
  },finally={multimodal_yfb_update_F_mk <- original})
})
run_test("FixedQ-T12: named scalar precision/signal match unnamed signed posterior", {
  g <- list(family="point_laplace",pi=.4,rate=1.3,point_mass=FALSE)
  plain <- multimodal_yfb_signed_posterior(4,2,g)
  named <- multimodal_yfb_signed_posterior(c(g1=4),c(g1=2),g)
  assert_equal(plain,named)
})
