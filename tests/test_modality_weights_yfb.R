# ============================================================
# Script: test_modality_weights_yfb.R
# Purpose: Verify reconstruction-only weights and separate prior pooling scope.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: existing feature-group tests, ebnm, survival
# ============================================================
if (!exists("fg_fixture")) source("tests/test_feature_groups_yfb.R")

#' Require a specific validation error rather than an unused-argument error
#' @param expr Expression expected to fail.
#' @param pattern Required message text.
#' @return Invisible TRUE after checking the exception.
#' @examples
#' mw_expect_error(stop("bad modality_weight"), "modality_weight")
mw_expect_error <- function(expr, pattern) {
  error <- tryCatch(force(expr), error = identity)
  assert_true(inherits(error, "error"))
  assert_true(grepl(pattern, conditionMessage(error), fixed = TRUE))
  assert_false(grepl("unused argument", conditionMessage(error), fixed = TRUE))
}

mw_groups <- rep(c("expression", "methylation"), each = 12)
mw_weights <- c(expression = 1, methylation = 0.2)

run_test("Modality weights: NULL and ones preserve grouped results bit for bit", {
  d <- fg_fixture()
  a <- fg_fit(d, list(feature_group = mw_groups))
  b <- fg_fit(d, list(feature_group = mw_groups,
                      modality_weight = c(expression = 1, methylation = 1)))
  for (field in c("EL", "EL2", "EF", "EF2", "EBeta", "EBeta2", "Tau", "history"))
    assert_true(identical(a[[field]], b[[field]]), field)
})
run_test("Modality weights: pooled ones preserve the ungrouped fit bit for bit", {
  d <- fg_fixture()
  a <- fg_fit(d)
  b <- fg_fit(d, list(feature_group = mw_groups, F_prior_by_group = FALSE,
                      modality_weight = c(expression = 1, methylation = 1)))
  for (field in c("EL", "EL2", "EF", "EF2", "EBeta", "EBeta2", "Tau", "history"))
    assert_true(identical(a[[field]], b[[field]]), field)
})
run_test("Modality weights: invalid names, values and missing groups fail loudly", {
  d <- fg_fixture()
  mw_expect_error(fg_fit(d, list(modality_weight = mw_weights)), "feature_group")
  for (weight in list(c(1, 0.2), c(expression = 1), c(expression = 1, methylation = 0),
      c(expression = 1, methylation = Inf), c(expression = 1, other = 0.2),
      c(expression = 1, expression = 0.2)))
    mw_expect_error(fg_fit(d, list(feature_group = mw_groups, modality_weight = weight)),
                 "modality_weight")
})
run_test("Modality weights: fitter L precision and signal use effective tau", {
  d <- fg_fixture(); d$args$max_iter <- 1L
  observer <- new.env(parent = environment(fit_ssbmf_yfb))
  original <- update_L_surv_YFB_k
  records <- list()
  observer$update_L_surv_YFB_k <- function(Tau, EF_k, EF2_k, R_k, ...) {
    ans <- original(Tau, EF_k, EF2_k, R_k, ...)
    records[[length(records) + 1L]] <<- list(tau = Tau, F = EF_k,
      F2 = EF2_k, R = R_k, answer = ans)
    ans
  }
  fitter <- fit_ssbmf_yfb; environment(fitter) <- observer
  set.seed(731)
  fit <- do.call(fitter, c(list(Y = d$Y, time = d$time, status = d$status), d$args,
    list(feature_group = mw_groups, modality_weight = mw_weights)))
  # The first sweep starts at tau=1/feature variance, which is 1 on this fixture.
  a <- records[[1L]]
  assert_near(a$tau, unname(mw_weights[mw_groups]), 1e-12)
  assert_near(a$answer$A, sum(unname(mw_weights[mw_groups]) * a$F2), 1e-10)
  assert_near(a$answer$B, drop(a$R %*% (unname(mw_weights[mw_groups]) * a$F)), 1e-10)
})
run_test("Modality weights: objective recomputes from posterior and returned priors", {
  d <- fg_fixture(); d$args$max_iter <- 1L
  observed <- fg_observe(d, mw_groups, list(modality_weight = mw_weights))
  fit <- observed$fit
  scales <- tail(observed$scales, 1)[[1]]$scale
  klf <- sum(vapply(seq_along(observed$calls), function(i) {
    k <- (i - 1L) %/% 2L + 1L; m <- (i - 1L) %% 2L + 1L
    call <- observed$calls[[i]]
    prior <- fit$prior_F[[m]][[k]]
    if (isTRUE(prior$point_mass)) 0 else fg_loading_kl(call, prior, scales[k])
  }, numeric(1)))
  r2 <- compute_expected_residual_sq(d$Y, fit$EL, fit$EL2, fit$EF, fit$EF2)
  gen <- sum(unname(mw_weights[mw_groups]) *
    (nrow(d$Y) / 2 * log(fit$Tau) - fit$Tau / 2 * colSums(r2)))
  Z <- d$Y %*% fit$EF
  VZ <- d$Y^2 %*% pmax(fit$EF2 - fit$EF^2, 0)
  tay <- calc_cox_taylor_yf(drop(Z %*% fit$EBeta), d$time, d$status)
  surv <- compute_survival_elbo(tay$logPL, tay$w, Z, Z^2 + VZ, fit$EBeta, fit$EBeta2)
  assert_near(tail(fit$history$elbo_full, 1), gen + surv + sum(observed$kls) - klf, 1e-7)
  recomputed <- yf_canonical_elbo(d$Y, d$time, d$status, fit$EL, fit$EL2,
    fit$EF, fit$EF2, fit$EBeta, fit$EBeta2, sum(fit$neg_kl),
    modality_weight = mw_weights, feature_group = mw_groups)
  assert_near(recomputed, tail(fit$history$elbo_full, 1), 1e-7)
})
run_test("Modality weights: positive raw unit changes vanish before supervised fit", {
  d <- fg_fixture(); raw <- d$Y
  scaled <- raw; scaled[, mw_groups == "methylation"] <- 4 * raw[, mw_groups == "methylation"]
  d$Y <- scale(raw)
  a <- fg_fit(d, list(feature_group = mw_groups, modality_weight = mw_weights))
  d$Y <- scale(scaled)
  b <- fg_fit(d, list(feature_group = mw_groups, modality_weight = mw_weights))
  for (field in c("EL", "EF", "EBeta", "Tau")) assert_near(a[[field]], b[[field]], 1e-9)
})
run_test("Modality weights: reconstruction coordinates are block-scale equivariant at zero survival", {
  set.seed(52)
  Y <- matrix(rnorm(60), 10, 6)
  F <- seq(-0.8, 1, length.out = 6); F2 <- F^2 + 0.1
  tau <- seq(1, 2, length.out = 6); labels <- rep(c("expression", "methylation"), each = 3)
  weight <- unname(mw_weights[labels]); c <- c(rep(1, 3), rep(2, 3))
  # Corresponding posterior states, without projection canonicalization or SVD.
  l <- update_L_surv_YFB_k(weight * tau, F, F2, Y, prior_family = "normal")
  l_scaled <- update_L_surv_YFB_k(weight * tau / c^2, F * c, F2 * c^2,
    sweep(Y, 2, c, "*"), prior_family = "normal")
  assert_near(l$mean, l_scaled$mean, 1e-10)
  f <- update_F_surv_YFB_k(weight * tau, l$mean, l$second, Y,
    prior_family = "normal", w_surv = 0, feature_group = labels)
  f_scaled <- update_F_surv_YFB_k(weight * tau / c^2, l_scaled$mean, l_scaled$second,
    sweep(Y, 2, c, "*"), prior_family = "normal", w_surv = 0, feature_group = labels)
  assert_near(f_scaled$mean, f$mean * c, 1e-9)
  assert_near(f_scaled$second, f$second * c^2, 1e-9)
  t <- update_tau(Y, matrix(l$mean), matrix(l$second), matrix(f$mean), matrix(f$second))$Tau
  ts <- update_tau(sweep(Y, 2, c, "*"), matrix(l_scaled$mean), matrix(l_scaled$second),
    matrix(f_scaled$mean), matrix(f_scaled$second))$Tau
  assert_near(ts, t / c^2, 1e-9)
  assert_near(weight * ts * f_scaled$second, weight * t * f$second, 1e-9)
})
run_test("Modality weights: rules use training z-scores and the sample Gram eigenvalue", {
  d <- fg_fixture()
  rules <- yf_modality_weight_rules(d$Y, mw_groups)
  assert_near(rules$feature_count, c(1, 1), 0)
  j <- mw_groups == "expression"
  eigen_ref <- eigen(crossprod(scale(d$Y[, j])) / nrow(d$Y), symmetric = TRUE)$values[1]
  eigen_m <- eigen(crossprod(scale(d$Y[, !j])) / nrow(d$Y), symmetric = TRUE)$values[1]
  assert_near(rules$eigenvalue, c(1, eigen_ref / eigen_m), 1e-10)
  # Changing units must not alter training-only rules.
  altered <- d$Y; altered[, !j] <- altered[, !j] * 4
  assert_near(yf_modality_weight_rules(altered, mw_groups)$eigenvalue, rules$eigenvalue, 1e-10)
})
run_test("Modality weights: prediction keeps the unweighted sum of block projections", {
  d <- fg_fixture()
  fit <- fg_fit(d, list(feature_group = mw_groups, modality_weight = mw_weights))
  pred <- predict_ssbmf_yfb(d$Y, fit$EF, fit$EBeta, fit$EF_norms)
  Z <- Reduce(`+`, lapply(unique(mw_groups), function(m) {
    j <- mw_groups == m
    d$Y[, j] %*% fit$EF[j, ]
  }))
  assert_near(pred$risk_scores, drop(Z %*% fit$EBeta), 1e-12)
})
if (sys.nframe() == 0L) {
  result <- report_results("Modality weights")
  if (result$failed > 0L) quit(status = 1L)
}
