# ============================================================
# Script: test_feature_groups_yfb.R
# Purpose: Verify feature-group priors without changing single-group YFB fits.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: ebnm, survival, existing YFB fitter and test helpers
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("fit_ssbmf_yfb")) source("code/fit_ssbmf_yfb.R")
if (!exists("predict_ssbmf_yfb")) source("code/predict_ssbmf_yfb.R")

# Fixtures ----

#' Reproducible signed input with a common prognostic score
#' @param stratified Whether to add two strata with distinct feature means.
#' @return List of data, outcomes and fixed fitter arguments.
#' @examples
#' fg_fixture(FALSE)
fg_fixture <- function(stratified = FALSE) {
  set.seed(20261112)
  n <- 48L; p <- 24L
  L <- matrix(rnorm(n * 2), n, 2)
  F <- cbind(c(seq(-2, 2, length.out = 12), rep(0, 12)),
             c(rep(0, 12), seq(-1.5, 1.5, length.out = 12)))
  Y <- scale(L %*% t(F) + matrix(rnorm(n * p, sd = 0.4), n, p))
  event <- rexp(n, exp(0.6 * L[, 1]))
  censor <- rexp(n, 0.4)
  args <- list(K = 2L, max_iter = 5L, tol = 1e-12,
               prior_LF = "point_laplace", prior_beta = "normal",
               canonicalize = TRUE, w_S = 1, verbose = FALSE)
  if (stratified) {
    args$strata_id <- rep(c("a", "b"), each = n / 2)
    args$center_projection <- TRUE
    Y[args$strata_id == "b", ] <- Y[args$strata_id == "b", ] + 2
  }
  list(Y = Y, time = pmin(event, censor), status = as.numeric(event <= censor),
       args = args)
}

#' Fit a fixed feature-group fixture with a reproducible start
#' @param d Output of fg_fixture().
#' @param extra Additional named fitter arguments.
#' @return YFB fit.
#' @examples
#' # fg_fit(fg_fixture(), list(feature_group = rep("one", 24)))
fg_fit <- function(d, extra = list()) {
  set.seed(731)
  do.call(fit_ssbmf_yfb, c(list(Y = d$Y, time = d$time, status = d$status),
                          d$args, extra))
}

#' Observe real F updates, KL evaluations and canonical rescaling
#'
#' Wrappers delegate every computation to the original functions. They only
#' record the calls needed to check the objective independently of its cache.
#' @param d Fixture, normally run for one sweep.
#' @param groups One label per feature.
#' @return List containing fit, F update calls/results, scalar KLs and rescaling.
#' @examples
#' # fg_observe(fg_fixture(), rep(c("rna", "meth"), each = 12))
fg_observe <- function(d, groups) {
  original_F <- multimodal_yfb_update_F_mk
  original_kl <- compute_ebnm_kl
  original_cz <- canonicalize_yf_factors
  calls <- list(); kls <- numeric(); scales <- list()
  assign("multimodal_yfb_update_F_mk", function(...) {
    args <- as.list(match.call(original_F, expand.dots = TRUE))[-1L]
    args <- lapply(args, eval, envir = parent.frame())
    ans <- original_F(...)
    calls[[length(calls) + 1L]] <<- list(args = args, result = ans)
    ans
  }, envir = .GlobalEnv)
  assign("compute_ebnm_kl", function(...) {
    ans <- original_kl(...)
    kls <<- c(kls, ans)
    ans
  }, envir = .GlobalEnv)
  assign("canonicalize_yf_factors", function(...) {
    ans <- original_cz(...)
    scales[[length(scales) + 1L]] <<- ans
    ans
  }, envir = .GlobalEnv)
  on.exit({
    assign("multimodal_yfb_update_F_mk", original_F, envir = .GlobalEnv)
    assign("compute_ebnm_kl", original_kl, envir = .GlobalEnv)
    assign("canonicalize_yf_factors", original_cz, envir = .GlobalEnv)
  })
  fit <- fg_fit(d, list(feature_group = groups))
  list(fit = fit, calls = calls, kls = kls, scales = scales)
}

#' Loading KL under a returned point-Laplace prior
#'
#' q is the Gaussian pseudo-likelihood times the incoming spike-and-Laplace
#' prior. Integrating E_q|F| supplies the prior-density correction when the
#' marginal-likelihood refit changes that prior. This catches a stale KL cache.
#' @param call Captured real block update from fg_observe().
#' @param returned Prior returned by the fitter, in canonical units.
#' @param scale Final canonical scale for this program.
#' @return Positive KL summed over the block, or Inf for incompatible support.
#' @examples
#' # fg_loading_kl(observed$calls[[1]], fit$prior_F$rna[[1]], scale[1])
fg_loading_kl <- function(call, returned, scale) {
  old <- call$args$prior
  v <- call$result$details
  new_pi <- if (isTRUE(returned$point_mass)) 0 else returned$pi
  new_rate <- if (new_pi > 0) returned$rate / scale else NA_real_
  total <- numeric(length(v$x))
  for (j in seq_along(total)) {
    slab <- v$slab_prob[j]
    spike <- 1 - slab
    if ((slab > 0 && new_pi == 0) || (spike > 0 && new_pi == 1)) return(Inf)
    correction <- 0
    if (spike > 0) correction <- spike * log((1 - old$pi) / (1 - new_pi))
    if (slab > 0) {
      integrand <- function(t) {
        lp <- log(old$pi * old$rate / 2) - old$rate * t - v$log_ml[j]
        t * (exp(lp + dnorm(v$x[j], t, sqrt(v$s2[j]), log = TRUE)) +
               exp(lp + dnorm(v$x[j], -t, sqrt(v$s2[j]), log = TRUE)))
      }
      absolute <- integrate(integrand, 0, Inf, rel.tol = 1e-9)$value
      correction <- correction + slab * log(old$pi * old$rate / (new_pi * new_rate)) +
        (new_rate - old$rate) * absolute
    }
    old_kl <- dnorm(v$x[j], 0, sqrt(v$s2[j]), log = TRUE) +
      (v$x[j] * v$mean[j] - v$second[j] / 2) / v$s2[j] - v$log_ml[j]
    total[j] <- old_kl + correction
  }
  sum(total)
}

# Baseline equivalence ----
# Frozen before feature-group implementation, from commit 19f910bd. Never
# regenerate these outputs using a fitter that already implements grouping.
# Only posterior moments and the monitored history are compared, so newly
# returned group metadata does not invalidate an unchanged numerical path.
for (stratified in c(FALSE, TRUE)) local({
  s <- stratified
  label <- if (s) "stratified_centered" else "pooled"
  run_test(paste("Feature groups: unchanged fitter matches frozen", label, "fit"), {
    fit <- fg_fit(fg_fixture(s))
    ref <- readRDS("tests/fixtures/yfb-feature-group-reference.rds")[[label]]
    for (field in c("EL", "EF", "EBeta"))
      assert_true(identical(fit[[field]], ref[[field]]), paste(label, field))
    assert_true(identical(fit$history$elbo_full, ref$history$elbo_full))
  })
  run_test(paste("Feature groups: NULL and one group preserve frozen", label, "fit"), {
    d <- fg_fixture(s)
    ref <- readRDS("tests/fixtures/yfb-feature-group-reference.rds")[[label]]
    for (extra in list(list(feature_group = NULL),
                       list(feature_group = rep("one", ncol(d$Y))))) {
      fit <- fg_fit(d, extra)
      for (field in c("EL", "EF", "EBeta"))
        assert_true(identical(fit[[field]], ref[[field]]), paste(label, field))
      assert_true(identical(fit$history$elbo_full, ref$history$elbo_full),
                  paste(label, "ELBO history"))
    }
  })
})

# Independent prior calculation ----

run_test("Feature groups: independent KL calculation accounts for a changed prior", {
  old <- list(family = "point_laplace", pi = 0.6, rate = 1.3, point_mass = FALSE)
  v <- multimodal_yfb_signed_posterior(2, 3, old)
  returned <- list(family = "point_laplace", pi = 0.7, rate = 2.1, point_mass = FALSE)
  call <- list(args = list(prior = old), result = list(details = v))
  log_q <- function(t) dnorm(v$x, t, sqrt(v$s2), log = TRUE) +
    log(old$pi * old$rate / 2) - old$rate * abs(t) - v$log_ml
  integrand <- function(t) {
    a <- log_q(t)
    b <- log(returned$pi * returned$rate / 2) - returned$rate * abs(t)
    answer <- exp(a) * (a - b)
    answer[!is.finite(a) & exp(a) == 0] <- 0
    answer
  }
  spike <- 1 - v$slab_prob
  expected <- spike * log(spike / (1 - returned$pi)) +
    integrate(integrand, -Inf, 0, rel.tol = 1e-9)$value +
    integrate(integrand, 0, Inf, rel.tol = 1e-9)$value
  assert_near(fg_loading_kl(call, returned, 1), expected, 1e-8)
  # Canonical units cannot change KL: F -> F/c implies rate -> c rate.
  rescaled <- returned; rescaled$rate <- rescaled$rate * 3
  assert_near(fg_loading_kl(call, rescaled, 3), expected, 1e-8)
})

run_test("Feature groups: noncanonical F prior equals ebnm within each group", {
  # Deliberately unequal sparsity. A pooled prior cannot satisfy both answers.
  x <- c(rep(0, 18), -4, 4, seq(-6, 6, length.out = 20))
  groups <- rep(c("sparse", "dense"), each = 20)
  update <- update_F_surv_YFB_k(rep(1, 40), 1, 1, matrix(x, 1),
                               prior_family = "point_laplace", alpha = 0,
                               feature_group = groups)
  priors <- list()
  for (g in unique(groups)) {
    j <- which(groups == g)
    ref <- ebnm::ebnm(x[j], rep(1, length(j)), prior_family = "point_laplace")
    assert_near(update$mean[j], ref$posterior$mean, 1e-10)
    assert_near(update$second[j], ref$posterior$mean^2 + ref$posterior$sd^2, 1e-10)
    assert_equal(update$ebnm_result[[g]]$fitted_g, ref$fitted_g)
    priors[[g]] <- ref$fitted_g
  }
  assert_true(abs(priors$sparse$pi[2] - priors$dense$pi[2]) > 0.2,
              "fixture must separate the two sparsity estimates")
})

run_test("Feature groups: canonical F refits use block observations and priors rescale", {
  d <- fg_fixture(); d$args$max_iter <- 1L
  groups <- rep(c("rna", "meth"), each = 12)
  observed <- fg_observe(d, groups)
  assert_length(observed$calls, 4L, "one real update per group and program")
  fit <- observed$fit
  assert_equal(names(fit$prior_F), unique(groups))
  scale <- tail(observed$scales, 1)[[1]]$scale
  for (k in seq_len(2)) for (m in seq_along(unique(groups))) {
    g <- unique(groups)[m]; j <- which(groups == g)
    call <- observed$calls[[(k - 1L) * 2L + m]]
    # The block sweep must use its actual columns, with the joint EZ carried in.
    assert_near(call$args$Y_m, d$Y[, j, drop = FALSE], 0)
    detail <- call$result$details
    prior <- call$result$prior
    ref <- ebnm::ebnm(detail$x, sqrt(detail$s2), prior_family = "point_laplace",
                      mode = 0,
                      g_init = if (isTRUE(call$args$prior$point_mass)) NULL else
                        ebnm::laplacemix(c(1 - call$args$prior$pi, call$args$prior$pi),
                                       c(0, 0), c(0, 1 / call$args$prior$rate)))
    stored <- fit$prior_F[[g]][[k]]
    if (isTRUE(prior$point_mass)) {
      assert_true(ref$fitted_g$pi[2] <= sqrt(.Machine$double.eps))
      assert_true(isTRUE(stored$point_mass))
    } else {
      assert_near(prior$pi, ref$fitted_g$pi[2], 1e-8)
      assert_near(1 / prior$rate, ref$fitted_g$scale[2], 1e-8)
      assert_near(stored$rate, prior$rate * scale[k], 1e-10)
      assert_near(stored$pi, prior$pi, 0)
    }
  }
  assert_true(any(vapply(observed$calls,
                         function(call) !isTRUE(call$result$prior$point_mass), logical(1))),
              "at least one fitted prior must exercise rescaling")
})

run_test("Feature groups: last objective includes every group KL and returned moments", {
  d <- fg_fixture(); d$args$max_iter <- 1L
  observed <- fg_observe(d, rep(c("rna", "meth"), each = 12))
  fit <- observed$fit
  assert_length(observed$calls, 4L)
  # L and per-program Normal beta each contribute one negative KL per program.
  assert_length(observed$kls, 4L)
  scale <- tail(observed$scales, 1)[[1]]$scale
  F_kl <- vapply(seq_along(observed$calls), function(i) {
    k <- (i - 1L) %/% 2L + 1L
    g <- c("rna", "meth")[(i - 1L) %% 2L + 1L]
    j <- which(rep(c("rna", "meth"), each = 12) == g)
    call <- observed$calls[[i]]
    assert_near(fit$EF[j, k], call$result$details$mean / scale[k], 1e-12)
    assert_near(fit$EF2[j, k], call$result$details$second / scale[k]^2, 1e-12)
    fg_loading_kl(call, fit$prior_F[[g]][[k]], scale[k])
  }, numeric(1))
  assert_finite(F_kl)
  neg_kl <- sum(observed$kls) - sum(F_kl)
  assert_near(sum(fit$neg_kl), neg_kl, 1e-8)
  value <- yf_canonical_elbo(d$Y, d$time, d$status, fit$EL, fit$EL2,
                             fit$EF, fit$EF2, fit$EBeta, fit$EBeta2,
                             neg_kl = neg_kl, w_S = 1)
  assert_near(tail(fit$history$elbo_full, 1), value, 1e-7)
})

run_test("Feature groups: concatenated prediction equals sum of block projections", {
  d <- fg_fixture()
  groups <- rep(c("rna", "meth"), each = 12)
  fit <- fg_fit(d, list(feature_group = groups))
  assert_near(fit$EF_norms, rep(1, ncol(fit$EF)), 0)
  set.seed(882); new_Y <- matrix(rnorm(7 * ncol(d$Y)), 7)
  predicted <- predict_ssbmf_yfb(new_Y, fit$EF, fit$EBeta, fit$EF_norms)
  block_Z <- lapply(unique(groups), function(g) {
    j <- which(groups == g)
    new_Y[, j, drop = FALSE] %*% fit$EF[j, , drop = FALSE]
  })
  assert_near(predicted$ZF_test, Reduce(`+`, block_Z), 1e-12)
  assert_near(predicted$risk_scores,
              as.vector(Reduce(`+`, block_Z) %*% fit$EBeta), 1e-12)
  assert_true(max(abs(predicted$risk_scores)) > 1e-6,
              "prediction fixture must have a nonzero predictor")
})

run_test("Feature groups: malformed labels and cohort offsets fail explicitly", {
  d <- fg_fixture()
  for (bad in list(rep("a", 23), c(rep("a", 23), NA_character_),
                   c(rep("a", 23), ""))) {
    error <- tryCatch(fg_fit(d, list(feature_group = bad)), error = identity)
    assert_true(inherits(error, "error"))
    assert_true(grepl("feature_group", conditionMessage(error), fixed = TRUE))
    assert_false(grepl("unused argument", conditionMessage(error), fixed = TRUE))
  }
  error <- tryCatch(fg_fit(d, list(feature_group = rep(c("a", "b"), each = 12),
                                 cohort_id = rep(c("a", "b"), each = 24))),
                    error = identity)
  assert_true(inherits(error, "error"))
  assert_true(grepl("cohort_id", conditionMessage(error), fixed = TRUE))
  assert_false(grepl("unused argument", conditionMessage(error), fixed = TRUE))
})

if (sys.nframe() == 0L) {
  results <- report_results("Feature-group YFB")
  if (results$failed > 0L) quit(status = 1L)
}
