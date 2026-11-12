# ============================================================
# Script: test_lb_yfb_study.R
# Purpose: Verify matched inputs, weighted LB updates and fold-in CV.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: ebnm, survival and project fitters
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("code/fit_ssbmf_lb.R")
source("code/fit_ssbmf_yfb.R")
source("code/predict_ssbmf_lb.R")
source("code/predict_ssbmf_yfb.R")
source("code/train_test_split.R")
if (file.exists("code/lb_yfb_study.R")) source("code/lb_yfb_study.R")

run_test("Study rejects mismatched matrices, folds and posterior starts", {
  a <- list(Y = matrix(1:12, 4), fold = 1:4, state = list(EL = matrix(1, 4, 2), EBeta = c(0, 0)))
  lb_yfb_study_assert_identical(a, a)
  for (field in names(a)) {
    b <- a
    if (field == "Y") b$Y[1, 1] <- 2L
    if (field == "fold") b$fold <- rev(b$fold)
    if (field == "state") b$state$EBeta[1] <- 1
    assert_true(inherits(tryCatch(lb_yfb_study_assert_identical(a, b), error = identity), "error"), field)
  }
})
run_test("Study scaling uses training cohort statistics for held-out rows", {
  X <- matrix(c(1, 3, 101, 103, 10, 14, 110, 114), 4, 2)
  co <- c("a", "a", "b", "b")
  # Both cohorts need two training records for an estimable SD.
  X <- rbind(X, c(1000, 1000), c(2000, 2000)); co <- c(co, "a", "b")
  z <- lb_yfb_study_scale(X, co, 1:4, "zscore")
  assert_near(z[5, ], (X[5, ] - colMeans(X[1:2, ])) / apply(X[1:2, ], 2, sd))
  assert_near(z[6, ], (X[6, ] - colMeans(X[3:4, ])) / apply(X[3:4, ], 2, sd))
})
run_test("Study concordance pools within-cohort pairs and freezes risk direction", {
  time <- c(1, 2, 3, 10, 11, 12); status <- rep(1, 6)
  risk <- c(3, 2, 1, 30, 20, 10); co <- rep(c("a", "b"), each = 3)
  c <- lb_yfb_study_concordance(risk, time, status, co)
  assert_near(c$C_within, 1)
  assert_true(c$C_pooled < c$C_within)
  bad <- lb_yfb_study_concordance(-risk, time, status, co)
  assert_near(bad$C_within, 0)
  unequal <- lb_yfb_study_concordance(c(4, 3, 2, 1, 10, 20), 1:6,
                                     rep(1, 6), c(rep("a", 4), rep("b", 2)))
  assert_near(unequal$C_within, 6 / 7, tol = 1e-12)
  assert_true(all(c("C_a", "C_b", "pairs_a", "pairs_b", "pairs_within", "pairs_pooled") %in% names(unequal)))
  assert_near(unequal$C_a, 1)
  assert_near(unequal$C_b, 0)
  assert_equal(unequal$pairs_a, 6)
  assert_equal(unequal$pairs_b, 1)
  assert_equal(unequal$pairs_within, 7)
  assert_equal(unequal$pairs_pooled, 15)
  # Near-equal times must be treated as ties, matching task A's timefix.
  tied <- lb_yfb_study_concordance(c(3, 1, 2), c(1, 1 + 1e-10, 2),
                                  rep(1, 3), rep("a", 3))
  assert_true("C_a" %in% names(tied))
  assert_near(tied$C_a, .5)

  assert_true(is.na(lb_yfb_study_concordance(rep(0, 6), time, status, co)$C_within))
})
set.seed(20261112)
.study_X <- outer(rexp(40), rexp(10)) + matrix(rexp(400, 4), 40, 10)
.study_time <- rexp(40) + seq_len(40) / 1000
.study_status <- rep(c(1L, 0L), 20)
.study_cohort <- rep(c("TCGA", "CPTAC"), each = 20)

run_test("Weighted LB first sweep matches untempered coordinate equations", {
  Y <- .study_X; K <- 2L
  st <- lb_yfb_study_start(Y, K, "point_exponential", FALSE)
  st$EBeta <- c(0.4, -0.3); st$EBeta2 <- st$EBeta^2
  f <- fit_ssbmf_lb(Y, .study_time, .study_status, K = K, max_iter = 1,
    prior_LF = "point_exponential", prior_beta = "normal", w_S = 10,
    EL_init = st$EL, EF_init = st$EF, EBeta_init = st$EBeta, verbose = FALSE)
  L <- st$EL; L2 <- L^2; F <- st$EF; F2 <- F^2; b <- st$EBeta; b2 <- b^2
  tau <- st$Tau; negL <- negF <- negB <- numeric(K)
  eta0 <- drop(L %*% b)
  tay <- calc_cox_taylor(eta0, .study_time, .study_status)
  z <- eta0 + tay$u / tay$w
  for (k in seq_len(K)) {
    R <- compute_R_k(Y, L, F, k); zk <- compute_z_no_k(z, L, b, k)
    # Independent equations with nonzero survival information: dropping w_S
    # or tempering genomics must change this reference posterior.
    A <- sum(tau * F2[, k]) + 10 * tay$w * b2[k]
    B <- drop(R %*% (tau * F[, k])) + 10 * tay$w * zk * b[k]
    x <- B / A
    a <- ebnm::ebnm(x, 1 / sqrt(A), prior_family = "point_exponential")
    L[, k] <- a$posterior$mean; L2[, k] <- a$posterior$sd^2 + L[, k]^2
    negL[k] <- compute_ebnm_kl(a$log_likelihood, A, x, L[, k], L2[, k])
    a <- update_F_k(tau, L[, k], L2[, k], compute_R_k(Y, L, F, k),
                    prior_family = "point_exponential", alpha = 0)
    F[, k] <- a$mean; F2[, k] <- a$second
    negF[k] <- compute_ebnm_kl(a$ebnm_result$log_likelihood, a$A, a$x, a$mean, a$second)
    a <- update_beta_k(tay$w, zk, L[, k], L2[, k], prior_family = "normal", alpha = 1)
    b[k] <- a$mean; b2[k] <- a$second
    negB[k] <- if (!is.null(a$neg_kl)) a$neg_kl else
      compute_ebnm_kl(a$ebnm_result$log_likelihood, a$A, a$x, a$mean, a$second)
  }
  assert_near(f$EL, L, tol = 1e-9); assert_near(f$EF, F, tol = 1e-9)
  assert_near(f$EBeta, b, tol = 1e-9)
  assert_near(f$foldin$alpha, 0, tol = 0)
  tay_end <- calc_cox_taylor(drop(L %*% b), .study_time, .study_status)
  expected <- update_tau(Y, L, L2, F, F2)$elbo_proxy +
    10 * (compute_survival_elbo(tay_end$logPL, tay_end$w, L, L2, b, b2) + sum(negB)) +
    sum(negL) + sum(negF)
  assert_near(tail(f$history$elbo_full, 1), expected, tol = 1e-7)
})
run_test("Weighted LB rejects incompatible and invalid controls", {
  for (args in list(list(w_S = -1), list(w_S = 1, normalize_AB = TRUE),
                    list(w_S = 1, alpha_schedule = list(warmup_iters = 1, ramp_iters = 2)),
                    list(w_S = 1, EBeta_init = NA_real_))) {
    err <- tryCatch(do.call(fit_ssbmf_lb, c(list(Y = .study_X, time = .study_time,
      status = .study_status, K = 2L, max_iter = 1, verbose = FALSE), args)), error = identity)
    assert_true(inherits(err, "error"))
  }
})
run_test("Study pair shares the complete starting state and uses fold-in prediction", {
  z <- lb_yfb_study_pair(.study_X, .study_time, .study_status, .study_cohort,
    K = 2L, prior_LF = "point_exponential", prior_beta = "normal", w_S = 1,
    max_iter = 3, tol = 1e-5, center_projection = TRUE)
  f <- z$fits$LB
  assert_near(z$predictions$LB, predict_ssbmf_lb_foldin(.study_X, f)$risk_scores, tol = 1e-10)
  assert_equal(f$initial_state, z$start)
  assert_equal(z$fits$YFB$initial_state, z$start)
  assert_near(z$predictions$YFB, predict_ssbmf_yfb(.study_X, z$fits$YFB$EF,
    z$fits$YFB$EBeta, z$fits$YFB$EF_norms,
    projection_center = z$fits$YFB$projection_center)$risk_scores, tol = 1e-10)
  assert_equal(nrow(z$program), 4L)
})
run_test("Study rejects a mutated recorded YFB start", {
  local({
    real_fit <- fit_ssbmf_yfb
    assign("fit_ssbmf_yfb", function(...) {
      fit <- real_fit(...)
      fit$initial_state$EBeta[1] <- 1
      fit
    }, envir = .GlobalEnv)
    on.exit(assign("fit_ssbmf_yfb", real_fit, envir = .GlobalEnv))
    err <- tryCatch(lb_yfb_study_pair(.study_X, .study_time, .study_status,
      .study_cohort, 2L, "point_exponential", "normal", 1, 1, 1e-5), error = identity)
    assert_true(inherits(err, "error"))
    assert_true(grepl("mismatch in state", conditionMessage(err)))
  })
})
run_test("Study CV rejects outcomes that would otherwise be silently truncated", {
  err <- tryCatch(lb_yfb_study_cv(.study_X, c(.study_time, 1), .study_status,
    .study_cohort, rep(1:5, 8), "sd", "point_exponential", "normal", 2L, 1, 3, 1e-5),
    error = identity)
  assert_true(inherits(err, "error"))
})
run_test("Study five-fold CV evaluates each patient once with the correct predictors", {
  folds <- create_stratified_folds(.study_status, 5, 42, cohort_id = .study_cohort)$fold_id
  cv <- lb_yfb_study_cv(.study_X, .study_time, .study_status, .study_cohort, folds,
    scaling = "sd", prior_LF = "point_exponential", prior_beta = "normal", K = 2L,
    w_S = 1, max_iter = 3, tol = 1e-5)
  assert_equal(nrow(cv$metrics), 10L)
  for (j in seq_len(5)) {
    te <- which(folds == j)
    z <- cv$pairs[[j]]
    Y <- lb_yfb_study_scale(.study_X, .study_cohort, which(folds != j), "sd")
    ref <- predict_ssbmf_lb_foldin(Y[te, , drop = FALSE], z$fits$LB)
    assert_near(z$heldout$LB, ref$risk_scores, tol = 1e-10)
    assert_near(z$heldout$YFB, predict_ssbmf_yfb(Y[te, , drop = FALSE], z$fits$YFB$EF,
      z$fits$YFB$EBeta, z$fits$YFB$EF_norms,
      projection_center = z$fits$YFB$projection_center)$risk_scores, tol = 1e-10)
  }
  assert_equal(sort(unlist(lapply(cv$pairs, function(z) z$test_indices))), seq_len(40))
})
run_test("Study saves earlier folds and explicit errors when a later fit fails", {
  local({
    real_pair <- lb_yfb_study_pair
    calls <- 0L
    assign("lb_yfb_study_pair", function(...) {
      calls <<- calls + 1L
      if (calls == 3L) stop("injected later-fold failure")
      real_pair(...)
    }, envir = .GlobalEnv)
    on.exit(assign("lb_yfb_study_pair", real_pair, envir = .GlobalEnv))
    out <- tempfile("study-fold-checkpoints-"); dir.create(out)
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    folds <- create_stratified_folds(.study_status, 5, 42, cohort_id = .study_cohort)$fold_id
    cv <- suppressWarnings(lb_yfb_study_cv(.study_X, .study_time, .study_status,
      .study_cohort, folds, "sd", "point_exponential", "normal", 2L, 1, 1, 1e-5,
      checkpoint = out))
    assert_true(cv$failed)
    assert_equal(length(list.files(out)), 5L)
    assert_true(!is.null(readRDS(file.path(out, "fold_1.rds"))$pair$fits$LB))
    failed <- cv$metrics[cv$metrics$fold == 3, ]
    assert_equal(failed$error, rep("injected later-fold failure", 2))
    assert_true(all(!failed$converged) && all(is.na(failed$C_within)))
    assert_equal(cv$summary$defined_pairs[1], 4L)
    assert_true(is.na(cv$summary$YFB_minus_LB[1]))
    assert_equal(readRDS(file.path(out, "fold_3.rds"))$pair$error,
                 "injected later-fold failure")
  })
})
if (sys.nframe() == 0L) {
  r <- report_results("test_lb_yfb_study.R")
  if (r$failed > 0) quit(status = 1)
}
