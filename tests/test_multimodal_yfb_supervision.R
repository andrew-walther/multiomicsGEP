# ============================================================
# Script: test_multimodal_yfb_supervision.R
# Purpose: Tests for survival supervision in the multimodal YFB fitter
#          (Workstream C6): shared beta prior update, supervised start,
#          beta-prior warm-up.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: ebnm; code/multimodal_yfb_updates.R, code/fit_multimodal_yfb.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("multimodal_yfb_coordinate_kl")) source("code/multimodal_yfb_updates.R")
if (!exists("fit_multimodal_yfb")) source("code/fit_multimodal_yfb.R")
if (!exists("simulate_multimodal_yfb_data")) source("code/simulate_multimodal_yfb.R")

local({
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train = 80L,
                                    n_validation = 20L, p_expression = 20L,
                                    p_methylation = 20L, seed = 3L, truncate = FALSE,
                                    noise_scale = 0.5, beta = c(0, 0, 3))
  .mms <<- list(Y = d$training$Y, time = d$training$time, event = d$training$event)
})
.mms_ctrl <- list(intercept = TRUE, prior_update = "ebnm",
                  prior_F = list(expression = "point_laplace", methylation = "point_laplace"))

run_test("MMYFB-Sup-T1: shared-prior beta update is the posterior under the fixed prior", {
  set.seed(2); EZ <- rnorm(40); VZ <- runif(40, 0, 0.1); w <- runif(40, 0.2, 1)
  h <- 0.6 * w * EZ + rnorm(40, sd = 0.3)
  g <- ebnm::laplacemix(c(0.4, 0.6), c(0, 0), c(0, 0.8))
  res <- multimodal_yfb_update_beta_shared_k(EZ, VZ, w, h, g)
  A <- sum(w * (EZ^2 + VZ)); B <- sum(EZ * h)
  ref <- ebnm::ebnm(B / A, sqrt(1 / A), prior_family = "point_laplace", g_init = g, fix_g = TRUE)
  assert_near(res$mean, ref$posterior$mean, tol = 1e-10, msg = "posterior mean")
  assert_near(res$second, ref$posterior$sd^2 + ref$posterior$mean^2, tol = 1e-10,
              msg = "posterior second moment")
  assert_true(res$kl >= -1e-10, msg = "KL nonnegative")
  zero <- multimodal_yfb_update_beta_shared_k(rep(0, 5), rep(0, 5), rep(1, 5), rep(0, 5), g)
  assert_true(zero$mean == 0 && zero$kl == 0, msg = "no information gives beta = 0")
})

run_test("MMYFB-Sup-T2: supervised start needs signed loadings and aligns factor K with u0", {
  err <- tryCatch(fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                                     control = list(init_supervised = TRUE, max_outer = 1L)),
                  error = function(e) e)
  assert_true(inherits(err, "error"), msg = "point-exponential loadings rejected")
  fit <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                            control = c(.mms_ctrl, list(init_supervised = TRUE, max_outer = 1L,
                                                        max_inner = 1L)))
  u0 <- multimodal_yfb_cox_working(.mms$time, .mms$event, rep(0, length(.mms$time)))$u
  Yc <- lapply(.mms$Y, function(x) sweep(x, 2, colMeans(x)))
  zK <- Reduce(`+`, Map(function(y, f) y %*% f[, 3], Yc, fit$EF))
  assert_true(abs(cor(zK, u0)) > 0.4,
              msg = sprintf("factor 3 score correlates with u0 (|cor| = %.2f)", abs(cor(zK, u0))))
})

run_test("MMYFB-Sup-T3: shared beta prior is held during the warm-up and refit after it", {
  ctrl <- c(.mms_ctrl, list(beta_prior = "shared_point_laplace", beta_warmup = 3L,
                            max_inner = 1L))
  g0 <- ebnm::laplacemix(c(0.5, 0.5), c(0, 0), c(0, 1))
  held <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                             control = c(ctrl, list(max_outer = 3L)))
  assert_true(identical(held$g_beta, g0), msg = "prior unchanged within the warm-up")
  refit <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                              control = c(ctrl, list(max_outer = 4L)))
  assert_true(!identical(refit$g_beta, g0), msg = "prior refit after the warm-up")
  assert_true(all(is.finite(refit$EBeta)), msg = "finite coefficients")
})

run_test("MMYFB-Sup-T4: protected nullcheck excludes survival-active programs", {
  fit <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 3L,
                            control = c(.mms_ctrl, list(max_outer = 2L)))
  fit$EBeta <- c(2, 0.01, 0); fit$EBeta2 <- fit$EBeta^2 + c(0.01, 1, 0)
  fit$kl$beta[1] <- 1e9  # deletion raises the ELBO in the unprotected state
  raw <- multimodal_yfb_nullcheck(fit, .mms$Y, .mms$time, .mms$event,
                                 intercept = TRUE)
  pro <- multimodal_yfb_nullcheck(fit, .mms$Y, .mms$time, .mms$event,
                                 intercept = TRUE, protect_survival = TRUE)
  assert_true(raw[1] > 0, msg = "active factor removal would increase ELBO")
  assert_true(pro[1] == -Inf, msg = "credible coefficient never a candidate")
  assert_near(pro[3], raw[3], tol = 0, msg = "zero-SD coefficient remains eligible")
  assert_near(pro[2], raw[2], tol = 0, msg = "unprotected candidate unchanged")
})

run_test("MMYFB-Sup-T5: default pruning sequence matches feature/supervised-F", {
  old <- new.env(parent = globalenv())
  src <- read_baseline_source("code/fit_multimodal_yfb.R")
  eval(parse(text = src), old)
  for (pr in c(FALSE, TRUE)) {
    a <- list(.mms$Y, .mms$time, .mms$event, 3L,
              control = c(.mms_ctrl, list(prune = pr, max_outer = 8L)))
    before <- do.call(old$fit_multimodal_yfb, a)
    now <- do.call(fit_multimodal_yfb, a)
    a$control$prune_protect_survival <- FALSE
    off <- do.call(fit_multimodal_yfb, a)
    assert_true(identical(before, now) && identical(now, off),
                msg = "default and explicit FALSE reproduce baseline, including pruning")
  }
})


run_test("MMYFB-Sup-T6: a dead program is selected for removal under protection", {
  f <- fit_multimodal_yfb(.mms$Y, .mms$time, .mms$event, 2L,
    control = c(.mms_ctrl, list(max_outer = 2L)))
  f$EL <- cbind(f$EL, 0); f$EL2 <- cbind(f$EL2, 0)
  f$EF <- lapply(f$EF, function(x) cbind(x, 0))
  f$EF2 <- lapply(f$EF2, function(x) cbind(x, 0))
  f$EBeta <- c(2, 2, 0); f$EBeta2 <- f$EBeta^2 + c(.01, .01, 0)
  f$kl <- list(L=c(f$kl$L, 0), F=lapply(f$kl$F, function(x) c(x, 0)), beta=c(f$kl$beta, 0))
  delta <- multimodal_yfb_nullcheck(f, .mms$Y, .mms$time, .mms$event,
    intercept = TRUE, protect_survival = TRUE)
  assert_near(delta[3], 0, tol = 1e-8, msg = "dead factor costs neither fit nor KL")
  assert_true(max(delta) >= 0 && which.max(delta) == 3L,
    msg = "existing greedy rule removes the dead program")
})
