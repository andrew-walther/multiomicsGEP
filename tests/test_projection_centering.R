# ============================================================
# Script: test_projection_centering.R
# Purpose: Tests for fit_ssbmf_yfb(center_projection = TRUE): nonnegative X/sd input in the
#          reconstruction, column-centered input in the Cox projection only.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: ebnm, survival; code/fit_ssbmf_yfb.R, code/predict_ssbmf_yfb.R,
#               code/simulate_multimodal_yfb.R (planted test bed only)
# ============================================================
# Derivation: Section "Centering Only the Cox Projection" of
# derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex. Each test
# checks one property stated there; verify_derivation.R checks the update formulas.

## Fixtures ----

#' Small dataset with one prognostic program, signed (z-scored) or nonnegative (X/sd)
#'
#' @param seed RNG seed.
#' @param n,p Dimensions.
#' @param nonnegative TRUE returns X/sd with X > 0 (column means about 3 SD); FALSE the z-score.
#' @return list(Y, time, status).
#' @examples
#' pc_fixture(1)
pc_fixture <- function(seed, n = 60, p = 20, nonnegative = FALSE) {
  set.seed(seed)
  L <- matrix(rexp(n * 3), n); F <- matrix(rexp(p * 3) * rbinom(p * 3, 1, .5), p)
  X <- L %*% t(F) + matrix(rnorm(n * p, sd = .5), n)
  X <- sweep(X, 2, pmax(-apply(X, 2, min), 0) + .1, "+")       # positive offsets
  eta <- 1.2 * as.vector(scale(L[, 2]))
  t <- rexp(n, .1 * exp(eta)); cz <- rexp(n, .05)
  Y <- if (nonnegative) sweep(X, 2, apply(X, 2, sd), "/") else scale(X)
  list(Y = Y, time = pmin(t, cz), status = as.numeric(t <= cz))
}

#' Fitter calls whose output must not change when center_projection is off
#'
#' Covers the historical path and the opt-in canonical controls that the
#' registered S-9, S-11, S-13 and S-14 runners use.
#' @return Named list of argument lists for fit_ssbmf_yfb() (data added by the caller).
#' @examples
#' names(pc_reference_calls())
pc_reference_calls <- function() {
  base <- list(K = 3, max_iter = 15, tol = 1e-8, verbose = FALSE)
  canon <- c(base, list(canonicalize = TRUE, prior_LF = "point_laplace", prior_beta = "normal"))
  list(
    historical          = c(base, list(prior_LF = "point_laplace", prior_beta = "normal")),
    canonical_w1        = c(canon, list(w_S = 1, trace_supervision = TRUE)),
    canonical_ws20      = c(canon, list(w_S = 20)),
    block_weight        = c(canon, list(w_S = 5, survival_block_weight = TRUE)),
    joint_beta          = c(canon, list(w_S = 5, survival_block_weight = TRUE,
                                        beta_posterior = "joint")),
    shared_prune        = c(canon, list(beta_prior_mode = "shared", prune = TRUE)),
    profiled            = c(canon, list(profile_cox_intercept = TRUE)),
    output_rule         = c(canon, list(converge_on = "outputs")),
    stratified          = c(canon, list(strata_id = rep(1:2, length.out = 60))),
    exponential_nonneg  = c(base, list(canonicalize = TRUE, prior_LF = "point_exponential",
                                       prior_beta = "point_normal", nonnegative = TRUE))
  )
}

#' Run one reference call on its fixture
#' @param args One element of pc_reference_calls().
#' @return The fit.
#' @examples
#' # pc_run(pc_reference_calls()$canonical_w1)
pc_run <- function(args) {
  d <- pc_fixture(21, nonnegative = isTRUE(args$nonnegative))
  args$nonnegative <- NULL
  set.seed(7)
  do.call(fit_ssbmf_yfb, c(list(d$Y, d$time, d$status), args))
}

## Option off: bit for bit ----

# The reference fits were produced by the fitter before center_projection existed
# (commit 3c990d2), on this machine, with pc_run() over pc_reference_calls() and
# saveRDS(). Regenerate only from that commit.
run_test("Projection centering off: every reference configuration matches the earlier fitter bit for bit", local({
  ref <- readRDS("tests/fixtures/yfb-projection-centering-off-reference.rds")
  calls <- pc_reference_calls()
  assert_equal(sort(names(ref)), sort(names(calls)))
  for (nm in names(calls)) assert_true(identical(pc_run(calls[[nm]]), ref[[nm]]),
    paste(nm, "differs from the earlier fitter. If R, BLAS or ebnm changed rather than the",
          "fitter, regenerate the fixture from commit 3c990d2"))
}))

run_test("Projection centering off: an explicit FALSE equals omitting the argument", local({
  args <- pc_reference_calls()$block_weight
  a <- pc_run(args); b <- pc_run(c(args, list(center_projection = FALSE)))
  assert_true(identical(a, b))
  assert_true(is.null(a$projection_center))
}))

## Argument checks ----

run_test("Projection centering: invalid values and unsupported combinations stop", local({
  d <- pc_fixture(22, nonnegative = TRUE)
  f <- function(...) fit_ssbmf_yfb(d$Y, d$time, d$status, K = 2, max_iter = 2, verbose = FALSE, ...)
  for (bad in list(NA, "yes", c(TRUE, TRUE)))
    assert_true(inherits(tryCatch(f(canonicalize = TRUE, center_projection = bad),
                                  error = identity), "error"))
  e <- tryCatch(f(center_projection = TRUE), error = conditionMessage)
  assert_true(grepl("canonicalize", e), e)
  e <- tryCatch(f(canonicalize = TRUE, prior_LF = "point_laplace", prior_beta = "normal",
                  center_projection = TRUE, profile_cox_intercept = TRUE), error = conditionMessage)
  assert_true(grepl("profile_cox_intercept", e), e)
  e <- tryCatch(f(canonicalize = TRUE, center_projection = TRUE, strata_id = rep(1:2, 3)),
                error = conditionMessage)
  assert_true(grepl("strata_id", e), e)
}))

## Shift invariance of the Cox projection ----

#' Fit with F held at a fixed start, so only q(beta) and q(L) move
#' @param Y Data matrix. @param d Fixture (time, status). @param EF0 Fixed loadings.
#' @param center Passed to center_projection.
#' @return The fit.
#' @examples
#' # pc_frozen_fit(d$Y, d, EF0, TRUE)
pc_frozen_fit <- function(Y, d, EF0, center, strata = NULL) {
  EL0 <- Y %*% EF0 / sum(EF0^2)
  fit_ssbmf_yfb(Y, d$time, d$status, K = ncol(EF0), max_iter = 6, tol = 1e-12,
                prior_LF = "point_exponential", prior_beta = "normal", canonicalize = TRUE,
                EL_init = EL0, EF_init = EF0, N_frozen = 6, center_projection = center,
                strata_id = strata, trace_supervision = TRUE, verbose = FALSE)
}

run_test("Projection centering: a constant added to one gene leaves the Cox projection unchanged", local({
  d <- pc_fixture(23, nonnegative = TRUE)
  set.seed(1); EF0 <- matrix(rexp(ncol(d$Y) * 2), ncol(d$Y), 2)
  Ys <- d$Y; Ys[, 5] <- Ys[, 5] + 4                    # shift gene 5 for every patient
  # With F fixed, q(beta) depends on the data only through the Cox projection
  a <- pc_frozen_fit(d$Y, d, EF0, TRUE); b <- pc_frozen_fit(Ys, d, EF0, TRUE)
  assert_near(a$EBeta, b$EBeta, 1e-10, "centered: beta means")
  assert_near(a$EBeta2, b$EBeta2, 1e-10, "centered: beta second moments")
  assert_near(a$projection_center[-5], b$projection_center[-5], 1e-12)
  assert_near(b$projection_center[5] - a$projection_center[5], 4, 1e-12)
  # Risk scores of new patients carrying the same shift are unchanged
  Yt <- d$Y[1:10, ]; Yts <- Yt; Yts[, 5] <- Yts[, 5] + 4
  ra <- predict_ssbmf_yfb(Yt, a$EF, a$EBeta, a$EF_norms, projection_center = a$projection_center)
  rb <- predict_ssbmf_yfb(Yts, b$EF, b$EBeta, b$EF_norms, projection_center = b$projection_center)
  assert_near(ra$risk_scores, rb$risk_scores, 1e-10)
  # The monitored survival term at a fixed state: unchanged by the shift when centered,
  # changed by it when the uncentered matrix is the projection (negative control)
  surv <- function(Yg, Ys) yf_canonical_elbo(Yg, d$time, d$status, a$EL, a$EL2, a$EF, a$EF2,
                                             a$EBeta, a$EBeta2, 0, 1, survival_Y = Ys) -
    update_tau(Yg, a$EL, a$EL2, a$EF, a$EF2)$elbo_proxy
  ctr <- function(Y) sweep(Y, 2, colMeans(Y))
  assert_near(surv(d$Y, ctr(d$Y)), surv(Ys, ctr(Ys)), 1e-8)
  assert_true(abs(surv(d$Y, d$Y) - surv(Ys, Ys)) > 1e-3, "uncentered term depends on the shift")
  # Uncentered, the common component's surrogate curvature shrinks beta to zero here
  a0 <- pc_frozen_fit(d$Y, d, EF0, FALSE)
  assert_true(max(abs(a$EBeta)) > 0.1 && max(abs(a0$EBeta)) < 1e-6,
              sprintf("centered max|beta| %.3g, uncentered %.3g", max(abs(a$EBeta)),
                      max(abs(a0$EBeta))))
}))

run_test("Projection centering: training risk is Yc E[F] E[beta] and differs from Y E[F] E[beta] by a constant", local({
  d <- pc_fixture(24, nonnegative = TRUE)
  set.seed(7)
  f <- fit_ssbmf_yfb(d$Y, d$time, d$status, K = 3, max_iter = 20, prior_LF = "point_exponential",
                     prior_beta = "normal", canonicalize = TRUE, center_projection = TRUE,
                     verbose = FALSE)
  assert_near(f$projection_center, colMeans(d$Y), 1e-14)
  rc <- predict_ssbmf_yfb(d$Y, f$EF, f$EBeta, f$EF_norms, projection_center = f$projection_center)$risk_scores
  ru <- predict_ssbmf_yfb(d$Y, f$EF, f$EBeta, f$EF_norms)$risk_scores
  assert_near(mean(rc), 0, 1e-10, "centered training risk has mean zero")
  assert_near(ru - rc, rep(sum(colMeans(d$Y) * (f$EF %*% f$EBeta)), nrow(d$Y)), 1e-10)
}))

## Stratified baseline hazard: centering within strata ----

run_test("Projection centering with strata: centered within each stratum, and a per-stratum constant leaves the Cox projection unchanged", local({
  d <- pc_fixture(28, n = 80, nonnegative = TRUE)
  g <- rep(c("a", "b"), length.out = nrow(d$Y))
  Y <- d$Y; Y[g == "b", ] <- Y[g == "b", ] + 1          # strata differ in level
  set.seed(1); EF0 <- matrix(rexp(ncol(Y) * 2), ncol(Y), 2)
  a <- pc_frozen_fit(Y, d, EF0, TRUE, g)
  means <- rbind(a = colMeans(Y[g == "a", ]), b = colMeans(Y[g == "b", ]))
  assert_near(a$projection_center, means, 1e-14)
  assert_equal(rownames(a$projection_center), c("a", "b"))
  # Gene 5 shifted by 4 in stratum b only: same beta posterior when centered within strata
  Ys <- Y; Ys[g == "b", 5] <- Ys[g == "b", 5] + 4
  b <- pc_frozen_fit(Ys, d, EF0, TRUE, g)
  assert_near(a$EBeta, b$EBeta, 1e-10); assert_near(a$EBeta2, b$EBeta2, 1e-10)
  # Scores are standardized on the centered projection from the start (initial scaling)
  assert_near(a$history$trace$sd_YF[1, ], rep(1, 2), 1e-10)
  # Training risk is the within-stratum centered projection; scores move within a stratum
  # only by a constant when a different row of means is used
  Yc <- Y - means[g, ]
  r <- predict_ssbmf_yfb(Y, a$EF, a$EBeta, a$EF_norms, projection_center = a$projection_center,
                         projection_strata = g)$risk_scores
  assert_near(r, as.vector(Yc %*% a$EF %*% a$EBeta), 1e-10)
  r_a <- predict_ssbmf_yfb(Y[g == "b", ], a$EF, a$EBeta, a$EF_norms,
                           projection_center = a$projection_center["a", ])$risk_scores
  assert_near(r_a - r[g == "b"], rep(r_a[1] - r[g == "b"][1], sum(g == "b")), 1e-10)
  e <- tryCatch(predict_ssbmf_yfb(Y, a$EF, a$EBeta, a$EF_norms, projection_center = a$projection_center,
                                  projection_strata = rep("c", nrow(Y))), error = conditionMessage)
  assert_true(grepl("not seen in training", e), e)
  # Negative control: one global center leaves the between-stratum shift in the survival term
  surv <- function(Yd, Ys) yf_canonical_elbo(Yd, d$time, d$status, a$EL, a$EL2, a$EF, a$EF2,
    a$EBeta, a$EBeta2, 0, 1, strata = g, survival_Y = Ys) -
    update_tau(Yd, a$EL, a$EL2, a$EF, a$EF2)$elbo_proxy
  Ysc <- Ys - rbind(a = colMeans(Ys[g == "a", ]), b = colMeans(Ys[g == "b", ]))[g, ]
  assert_near(surv(Y, Yc), surv(Ys, Ysc), 1e-8)
  glob <- function(M) sweep(M, 2, colMeans(M))
  assert_true(abs(surv(Y, glob(Y)) - surv(Ys, glob(Ys))) > 1e-3, "global centering depends on the shift")
}))

run_test("Projection centering with strata: the monitored ELBO and the pruning check use the within-stratum centered projection", local({
  d <- pc_fixture(30, n = 80, nonnegative = TRUE)
  g <- rep(1:2, length.out = nrow(d$Y))
  Y <- d$Y; Y[g == 2, ] <- Y[g == 2, ] + 1               # strata differ in level
  seen <- NULL
  real_nullcheck <- prune_yf_nullcheck
  # Record what the pruning check receives (pruning rarely fires on a small fixture)
  assign("prune_yf_nullcheck", function(fit, Y, time, status, strata = NULL,
                                        protect_survival = FALSE, survival_Y = Y) {
    seen <<- list(Y = Y, survival_Y = survival_Y, strata = strata)
    rep(-1, ncol(fit$EF))
  }, envir = globalenv())
  on.exit(assign("prune_yf_nullcheck", real_nullcheck, envir = globalenv()))
  set.seed(7)
  f <- fit_ssbmf_yfb(Y, d$time, d$status, K = 2, max_iter = 8, tol = 1e-12, strata_id = g,
                     prior_LF = "point_exponential", prior_beta = "point_normal",
                     canonicalize = TRUE, center_projection = TRUE, prune = TRUE, verbose = FALSE)
  Yc <- Y - f$projection_center[g, ]
  assert_near(seen$survival_Y, Yc, 1e-14); assert_near(seen$Y, Y, 0)
  assert_equal(seen$strata, g)
  # Survival part of the monitored objective, recomputed from the returned state
  surv <- function(Ys) yf_canonical_elbo(Y, d$time, d$status, f$EL, f$EL2, f$EF, f$EF2, f$EBeta,
    f$EBeta2, 0, f$w_S, strata = g, survival_Y = Ys) - update_tau(Y, f$EL, f$EL2, f$EF, f$EF2)$elbo_proxy
  monitored <- tail(f$history$elbo_full, 1) - update_tau(Y, f$EL, f$EL2, f$EF, f$EF2)$elbo_proxy -
    sum(f$neg_kl)
  assert_true(max(abs(f$EBeta)) > 0.5, "coefficients are live")
  assert_near(surv(Yc), monitored, 1e-6)
  assert_true(abs(surv(sweep(Y, 2, colMeans(Y))) - monitored) > 1e-3,
              "global centering gives a different survival term")
}))

## Compatibility with the other opt-in controls ----

run_test("Projection centering runs with the survival-block weight, joint beta, a shared prior and pruning", local({
  d <- pc_fixture(25, nonnegative = TRUE)
  common <- list(d$Y, d$time, d$status, K = 3, max_iter = 20, prior_LF = "point_exponential",
                 canonicalize = TRUE, center_projection = TRUE, verbose = FALSE)
  fits <- list(
    block = do.call(fit_ssbmf_yfb, c(common, list(prior_beta = "point_normal", w_S = 5,
                                                  survival_block_weight = TRUE))),
    joint = do.call(fit_ssbmf_yfb, c(common, list(prior_beta = "normal", w_S = 5,
                                                  survival_block_weight = TRUE,
                                                  beta_posterior = "joint"))),
    prune = do.call(fit_ssbmf_yfb, c(common, list(prior_beta = "point_normal",
                                                  beta_prior_mode = "shared", prune = TRUE))))
  Yc <- sweep(d$Y, 2, colMeans(d$Y))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    assert_true(all(is.finite(f$history$elbo_full)), nm)
    assert_near(f$projection_center, colMeans(d$Y), 1e-14, nm)
    assert_true(all(f$EF >= 0) && all(f$EL >= 0), paste(nm, "nonnegative loadings"))
  }
  # The last monitored value, recomputed from the returned state with the centered
  # projection (variance terms included), for the block-weight and joint-beta fits
  for (nm in c("block", "joint")) {
    f <- fits[[nm]]
    Z <- Yc %*% f$EF; VZ <- Yc^2 %*% pmax(f$EF2 - f$EF^2, 0)
    tay <- calc_cox_taylor_yf(as.vector(Z %*% f$EBeta), d$time, d$status)
    var_eta <- if (nm == "joint") rowSums((Z %*% f$EBeta_cov) * Z) + as.vector(VZ %*% f$EBeta2)
               else as.vector((Z^2 + VZ) %*% f$EBeta2 - Z^2 %*% f$EBeta^2)
    elbo <- update_tau(d$Y, f$EL, f$EL2, f$EF, f$EF2)$elbo_proxy +
      f$w_S * (tay$logPL - 0.5 * sum(tay$w * var_eta)) + sum(f$neg_kl)
    last <- tail(f$history$elbo_full, 1)
    assert_near(elbo, last, 1e-6 * abs(last), paste(nm, "monitored ELBO"))
    assert_true(sum(VZ) > 0, paste(nm, "variance term exercised"))
  }
}))

run_test("Projection centering: with K = 1 the joint beta fit equals the mean-field Normal fit", local({
  # K = 1 joint q(beta) is the Normal-family EBNM update (S-13 identity), so both beta
  # updates must see the same centered scores and score variances
  d <- pc_fixture(27, n = 200, p = 20, nonnegative = TRUE)
  f <- function(...) { set.seed(7)
    fit_ssbmf_yfb(d$Y, d$time, d$status, K = 1, max_iter = 10, tol = 1e-12, init_supervised = TRUE,
                  prior_LF = "point_exponential", prior_beta = "normal", canonicalize = TRUE,
                  center_projection = TRUE, verbose = FALSE, ...) }
  a <- f(); b <- f(beta_posterior = "joint")
  assert_true(abs(a$EBeta) > 0.1, "coefficient is live, so the variance terms matter")
  assert_near(a$EBeta, b$EBeta, 1e-8); assert_near(a$EF, b$EF, 1e-9)
  assert_near(a$history$elbo_full, b$history$elbo_full, 1e-4)
}))

run_test("Projection centering: the monitored survival term uses the centered projection", local({
  d <- pc_fixture(26, nonnegative = TRUE)
  set.seed(7)
  f <- fit_ssbmf_yfb(d$Y, d$time, d$status, K = 2, max_iter = 8, tol = 1e-12,
                     prior_LF = "point_exponential", prior_beta = "normal", canonicalize = TRUE,
                     center_projection = TRUE, verbose = FALSE)
  # The pruning objective at the returned state, with the centered survival design, equals
  # the last monitored value up to the end-of-iteration rescaling (which leaves it unchanged)
  Yc <- sweep(d$Y, 2, colMeans(d$Y))
  full <- yf_canonical_elbo(d$Y, d$time, d$status, f$EL, f$EL2, f$EF, f$EF2, f$EBeta,
                            f$EBeta2, sum(f$neg_kl), f$w_S, survival_Y = Yc)
  assert_near(full, tail(f$history$elbo_full, 1), 1e-6 * abs(full))
  unc <- yf_canonical_elbo(d$Y, d$time, d$status, f$EL, f$EL2, f$EF, f$EF2, f$EBeta,
                           f$EBeta2, sum(f$neg_kl), f$w_S)
  assert_true(abs(unc - full) > 1e-6 * abs(full), "uncentered survival term differs")
}))

## Planted test bed: point-exponential programs on nonnegative X/sd input ----

#' October 15 planted test bed (DECISIONS.md 2026-10-15), expression block, two preprocessings
#'
#' simulate_multimodal_yfb_data("low_variance_prognostic_factor"), n = 120, p = 50,
#' coefficient 1.10 on program 3, noise 0.25. The raw training block gets a positive
#' offset per gene (as in the S-7 generator) so it is nonnegative; X/sd divides by
#' the training SD, and the z-score is the input the test bed used on 2026-10-15.
#' @param seed RNG seed.
#' @return list(Y_sd, Y_z, time, status, F3).
#' @examples
#' # pc_testbed(1)
pc_testbed <- function(seed) {
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train = 120L,
                                    n_validation = 120L, beta = c(0, 0, 1.10),
                                    p_expression = 50L, p_methylation = 50L, seed = seed,
                                    truncate = FALSE, noise_scale = 0.25)
  X <- d$training$Y$expression
  X <- sweep(X, 2, pmax(-apply(X, 2, min), 0) + .1, "+")
  list(Y_sd = sweep(X, 2, apply(X, 2, sd), "/"), Y_z = scale(X),
       time = d$training$time, status = d$training$event, F3 = d$truth$F$expression[, 3])
}

run_test("Projection centering: point-exponential programs on X/sd input do not collapse on the planted test bed", local({
  # Registered test-bed settings: K = 7, at most 300 sweeps, tol 1e-5; point-exponential
  # L and F, point-normal beta, w_S = 1.
  fit <- function(Y, tb, center) fit_ssbmf_yfb(Y, tb$time, tb$status, K = 7, max_iter = 300,
    tol = 1e-5, prior_LF = "point_exponential", prior_beta = "point_normal", canonicalize = TRUE,
    w_S = 1, center_projection = center, verbose = FALSE)
  live <- function(f) sum(colSums(f$EF^2) > 1e-8 & colSums(f$EL^2) > 1e-8)
  active <- function(f) { sdb <- sqrt(pmax(f$EBeta2 - f$EBeta^2, 0))
                          sum(abs(f$EBeta) >= 1.96 * sdb & sdb > 0) }
  for (s in 1:3) {
    tb <- pc_testbed(s)
    set.seed(s); f <- fit(tb$Y_sd, tb, TRUE)
    assert_true(live(f) >= 2, sprintf("seed %d: %d live programs", s, live(f)))
    assert_true(isTRUE(f$history$converged), sprintf("seed %d converged", s))
    assert_true(all(f$EF >= 0) && all(f$EL >= 0), sprintf("seed %d nonnegative", s))
    # Survival reaches the coefficients with centering; matched control without it
    assert_true(active(f) >= 1, sprintf("seed %d: no survival-active program", s))
    set.seed(s); f0 <- fit(tb$Y_sd, tb, FALSE)
    assert_true(active(f0) == 0, sprintf("seed %d: uncentered control has %d active", s, active(f0)))
  }
  # Not claimed: loading recovery. On seeds 1-6 the best correlation with the planted
  # loadings ranged from 0.08 to 0.84, and in some seeds the active program is not the
  # planted one (2026-11-06 check). For reference, the same prior on the z-scored input
  # left 1, 1 and 0 live programs on seeds 1-3.
}))
