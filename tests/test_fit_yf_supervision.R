# =============================================================================
# tests/test_fit_yf_supervision.R
# Tests for the survival-supervision diagnostics of fit_cox_on_yf():
#   - trace_supervision (opt-in per-iteration record; must not change the fit)
#   - w_S (survival weight in the F update; w_S = 1 is the equal-weight posterior)
# Workstream C, step C1 (docs/plans/Meeting_Notes_Followup_10_2_26.md).
# =============================================================================

local({
  set.seed(311)
  n <- 50L; p <- 15L
  Y <- scale(matrix(rnorm(n * p), n, p))
  eta <- 0.8 * as.vector(Y[, 1:3] %*% c(1, 1, 1)) / sqrt(3)
  .sv <<- list(n = n, p = p, Y = Y,
               time = rexp(n, 0.1 * exp(eta)), status = rbinom(n, 1, 0.8))
})

strip_trace <- function(fit) { fit$history$trace <- NULL; fit }

# ---------------------------------------------------------------------------
# T1: trace off vs on gives a bit-for-bit identical fit
# ---------------------------------------------------------------------------
run_test("Supervision-T1: trace_supervision does not change the fit (alpha_F = 0.5)", {
  args <- list(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 15L,
               prior_beta = "normal", alpha_F = 0.5, verbose = FALSE)
  off <- do.call(fit_cox_on_yf, args)
  on  <- do.call(fit_cox_on_yf, c(args, trace_supervision = TRUE))
  assert_true(is.null(off$history$trace), msg = "no trace when off")
  assert_true(identical(off, strip_trace(on)), msg = "fit identical apart from the trace")
  assert_true(nrow(on$history$trace$beta_mean) == on$history$n_iter,
              msg = "trace rows truncated to n_iter")
  assert_true(all(on$history$trace$beta_mean[on$history$n_iter, ] == on$EBeta |
                  on$history$trace$beta_mean[on$history$n_iter, ] == -on$EBeta),
              msg = "last traced beta equals the returned beta up to the Phase C flip")
})

# ---------------------------------------------------------------------------
# T2: logged shares match a hand computation of the first F update (K = 1)
# With K = 1 and EBeta starting at 0, the first iteration is: SVD init ->
# Cox Taylor at eta = 0 -> beta update -> L update -> F update.
# ---------------------------------------------------------------------------
run_test("Supervision-T2: logged survival share matches a hand computation", {
  Y <- .sv$Y; a_F <- 0.5
  fit <- fit_cox_on_yf(Y, .sv$time, .sv$status, K = 1L, max_iter = 1L,
                       prior_beta = "normal", alpha_F = a_F, verbose = FALSE,
                       trace_supervision = TRUE)
  s  <- svd(Y, nu = 1, nv = 1)
  EL <- abs(s$u * sqrt(s$d[1])); EF <- abs(s$v * sqrt(s$d[1]))
  Tau <- 1 / pmax(apply(Y, 2, var), 1e-8)
  ZF <- Y %*% (EF / sqrt(sum(EF^2) + 1e-10))
  tay <- calc_cox_taylor_yf(rep(0, .sv$n), .sv$time, .sv$status)
  z <- tay$u / tay$w; w <- tay$w
  rb <- update_beta_k(w, z, ZF[, 1], ZF[, 1]^2, prior_family = "normal", alpha = 0.5)
  rl <- update_L_surv_YFB_k(Tau, EF[, 1], EF[, 1]^2, Y, prior_family = "point_exponential")
  rf <- update_F_surv_YFB_k(Tau, rl$mean, rl$second, Y,
                            EBeta_k = rb$mean, EBeta2_k = rb$second,
                            YtWY_diag = as.vector(t(Y^2) %*% w),
                            YtWz_no_k = as.vector(t(Y) %*% (w * z)),
                            alpha = a_F)
  share <- a_F * rf$A_surv / rf$A
  assert_near(fit$history$trace$shareA_median[1, 1], median(share), tol = 1e-10,
               msg = "median precision share")
  assert_near(fit$history$trace$shareA_p90[1, 1], quantile(share, 0.9, names = FALSE),
               tol = 1e-10, msg = "90th percentile precision share")
  assert_near(fit$history$trace$beta_mean[1, 1], rb$mean, tol = 1e-10, msg = "E[beta]")
  assert_near(fit$history$trace$sd_YF[1, 1], sd(as.vector(Y %*% rf$mean)), tol = 1e-10,
               msg = "sd(Y EF_k)")
  assert_true(all(share >= 0 & share <= 1), msg = "share within [0, 1]")
})

# ---------------------------------------------------------------------------
# T3: w_surv = 1 is alpha = 0.5 with both precision terms doubled
# ---------------------------------------------------------------------------
run_test("Supervision-T3: w_surv = 1 equals alpha = 0.5 with doubled precisions", {
  set.seed(5); p <- 12L; n <- 30L
  Y <- matrix(rnorm(n * p), n, p); w <- runif(n, 0.2, 1); zz <- rnorm(n)
  EL <- abs(rnorm(n)); Tau <- runif(p, 0.5, 2)
  common <- list(Tau = Tau, EL_k = EL, EL2_k = EL^2 + 0.1, R_k = Y,
                 EBeta_k = 0.4, EBeta2_k = 0.25,
                 YtWY_diag = as.vector(t(Y^2) %*% w),
                 YtWz_no_k = as.vector(t(Y) %*% (w * zz)),
                 prior_family = "normal")
  half <- do.call(update_F_surv_YFB_k, c(common, alpha = 0.5))
  eq   <- do.call(update_F_surv_YFB_k, c(common, w_surv = 1))
  assert_near(eq$A, 2 * half$A, tol = 1e-12, msg = "precision doubled")
  assert_near(eq$x, half$x, tol = 1e-12, msg = "same pseudo-observation")
  assert_near(eq$A, half$A_gen + half$A_surv, tol = 1e-12, msg = "A = A_gen + A_surv")
})

# ---------------------------------------------------------------------------
# T4: w_S = 0 reproduces the alpha_F = 0 fit exactly; w_S > 0 changes it
# ---------------------------------------------------------------------------
run_test("Supervision-T4: w_S = 0 reproduces alpha_F = 0; w_S = 1 differs", {
  args <- list(.sv$Y, .sv$time, .sv$status, K = 2L, max_iter = 10L,
               prior_beta = "normal", verbose = FALSE)
  base <- do.call(fit_cox_on_yf, args)
  w0   <- do.call(fit_cox_on_yf, c(args, w_S = 0))
  w1   <- do.call(fit_cox_on_yf, c(args, w_S = 1))
  assert_true(identical(base$EF, w0$EF) && identical(base$EBeta, w0$EBeta),
              msg = "w_S = 0 identical to alpha_F = 0")
  assert_true(max(abs(base$EF - w1$EF)) > 1e-8, msg = "w_S = 1 moves F")
  err <- tryCatch(fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 2L, max_iter = 2L,
                                w_S = -1, verbose = FALSE), error = function(e) e)
  assert_true(inherits(err, "error"), msg = "negative w_S rejected")
})

# ---------------------------------------------------------------------------
# T5: one canonicalization step leaves L F', E[L F' second moments], eta and
# the survival ELBO term unchanged, and sets sd(Y EF_k) = 1
# ---------------------------------------------------------------------------
run_test("Supervision-T5: canonicalization preserves LF', eta and the ELBO terms", {
  set.seed(8); n <- 40L; p <- 12L; K <- 3L
  Y <- scale(matrix(rnorm(n * p), n, p))
  EL <- matrix(rnorm(n * K), n, K); EL2 <- EL^2 + 0.1
  EF <- matrix(rnorm(p * K), p, K) * c(1, 5, 0.2)[col(matrix(0, p, K))]
  EF2 <- EF^2 + 0.05
  EBeta <- c(0.3, -0.1, 0.8); EBeta2 <- EBeta^2 + 0.02
  cz <- canonicalize_yf_factors(Y, EL, EL2, EF, EF2, EBeta, EBeta2, K)
  assert_near(cz$EL %*% t(cz$EF), EL %*% t(EF), tol = 1e-10, msg = "L F'")
  assert_near(cz$EL2 %*% t(cz$EF2), EL2 %*% t(EF2), tol = 1e-10, msg = "E[l^2] E[f^2]")
  assert_near(Y %*% cz$EF %*% cz$EBeta, Y %*% EF %*% EBeta, tol = 1e-10, msg = "eta")
  assert_near(apply(Y %*% cz$EF, 2, sd), rep(1, K), tol = 1e-10, msg = "unit SD projections")
  w <- runif(n, 0.2, 1)
  e0 <- compute_survival_elbo(-50, w, Y %*% EF, (Y %*% EF)^2, EBeta, EBeta2)
  e1 <- compute_survival_elbo(-50, w, Y %*% cz$EF, (Y %*% cz$EF)^2, cz$EBeta, cz$EBeta2)
  assert_near(e1, e0, tol = 1e-8, msg = "survival ELBO term")
})

# ---------------------------------------------------------------------------
# T6: with w_S = 0, a canonical fit recovers the alpha_F = 0 factorization up
# to the scale of each program (survival enters neither L nor F in both)
# ---------------------------------------------------------------------------
run_test("Supervision-T6: canonical fit with w_S = 0 matches alpha_F = 0 up to scale", {
  set.seed(21); n <- 80L; p <- 30L
  Ys <- scale(matrix(rnorm(n * 2), n, 2) %*% matrix(rnorm(2 * p), 2, p) +
              matrix(rnorm(n * p, sd = 0.7), n, p))
  # Tight tolerance: the two monitors stop at different iterations at 1e-5
  args <- list(Ys, rexp(n), rbinom(n, 1, 0.7), K = 2L, max_iter = 600L, tol = 1e-9,
               prior_LF = "point_laplace", prior_beta = "normal", verbose = FALSE)
  leg <- do.call(fit_cox_on_yf, args)
  can <- do.call(fit_cox_on_yf, c(args, canonicalize = TRUE, w_S = 0))
  assert_true(all(colSums(leg$EF^2) > 0), msg = "both programs retained")
  r <- abs(diag(cor(leg$EF, can$EF)))
  assert_true(all(r > 0.9999), msg = sprintf("loading columns proportional (min |cor| %.5f)", min(r)))
  rel <- max(abs(leg$EL %*% t(leg$EF) - can$EL %*% t(can$EF))) / max(abs(leg$EL %*% t(leg$EF)))
  assert_true(rel < 1e-3, msg = sprintf("same reconstruction (rel diff %.2e)", rel))
})

# ---------------------------------------------------------------------------
# T7: canonical fit returns unit-SD projections, EF_norms = 1, and
# predict_cox_on_yf() on the training data returns Y EF beta
# ---------------------------------------------------------------------------
run_test("Supervision-T7: canonical fit is self-consistent for prediction", {
  fit <- fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 30L,
                       prior_LF = "point_laplace", prior_beta = "normal",
                       canonicalize = TRUE, verbose = FALSE)
  live <- colSums(fit$EF^2) > 1e-12
  assert_near(apply(.sv$Y %*% fit$EF[, live, drop = FALSE], 2, sd), rep(1, sum(live)),
              tol = 1e-8, msg = "sd(Y EF_k) = 1")
  assert_true(all(fit$EF_norms == 1) && fit$w_S == 1, msg = "EF_norms = 1, default w_S = 1")
  pr <- predict_cox_on_yf(.sv$Y, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
  assert_near(pr, as.vector(.sv$Y %*% fit$EF %*% fit$EBeta), tol = 1e-10, msg = "prediction")
  err <- tryCatch(fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 2L, max_iter = 2L,
                                canonicalize = TRUE, alpha_F = 0.5, verbose = FALSE),
                  error = function(e) e)
  assert_true(inherits(err, "error"), msg = "alpha_F rejected with canonicalize")
})

# ---------------------------------------------------------------------------
# T8: the survival part of the canonical F update is exact coordinate ascent.
# With negligible genomics precision, a flat prior and fixed beta, repeated
# sweeps with the fitter's inputs (working weights w, h = w * z_no_k)
# converge to the weighted least-squares solution of beta * Y f ~ z_no_k,
# f = (Y'WY)^{-1} Y'W z_no_k / beta. The diagonal (vector) update does not.
# ---------------------------------------------------------------------------
run_test("Supervision-T8: canonical F sweep solves the survival least squares", {
  set.seed(14); n <- 60L; p <- 6L
  Y <- matrix(rnorm(n * p), n, p) %*% chol(0.5 * diag(p) + 0.5)  # correlated genes
  w <- runif(n, 0.3, 1); znk <- rnorm(n); beta <- 0.7
  target <- solve(crossprod(Y, w * Y), crossprod(Y, w * znk)) / beta
  Tau <- rep(1e-10, p); EL <- rep(0, n)
  prior <- list(family = "normal", variance = 1e8, point_mass = FALSE)
  EF <- rep(0, p); EF2 <- rep(0, p)
  for (it in 1:400) {
    upd <- multimodal_yfb_update_F_mk(Y, Tau, EL, EL, Y, w, w * znk,
                                      as.vector(Y %*% EF), rep(0, n), EF, EF2,
                                      beta, beta^2, prior, prior_update = "ebnm")
    EF <- upd$EF; EF2 <- upd$EF2
  }
  assert_near(EF, as.vector(target), tol = 1e-5, msg = "sweep reaches the WLS solution")
  diag_update <- as.vector(crossprod(Y, w * znk)) / (beta * colSums(w * Y^2))
  assert_true(max(abs(diag_update - target)) > 0.05,
              msg = "the diagonal update differs from the exact solution")
})

# ---------------------------------------------------------------------------
# T9: shared beta prior -- fixed prior is never refit; warm-up holds the
# initial prior for exactly beta_warmup iterations, then it is refit
# ---------------------------------------------------------------------------
run_test("Supervision-T9: shared beta prior honors beta_prior_fix and beta_warmup", {
  base <- list(.sv$Y, .sv$time, .sv$status, K = 3L, prior_LF = "point_laplace",
               prior_beta = "normal", canonicalize = TRUE, verbose = FALSE,
               beta_prior_mode = "shared")
  g0 <- ashr::normalmix(1, 0, 2)
  fx <- do.call(fit_cox_on_yf, c(base, max_iter = 15L, beta_prior_init = list(g0),
                                 beta_prior_fix = TRUE))
  assert_true(identical(fx$g_beta, g0), msg = "fixed prior unchanged")
  wu <- do.call(fit_cox_on_yf, c(base, max_iter = 4L, beta_prior_init = list(g0),
                                 beta_warmup = 4L))
  assert_true(identical(wu$g_beta, g0), msg = "prior held during the warm-up")
  rf <- do.call(fit_cox_on_yf, c(base, max_iter = 5L, beta_prior_init = list(g0),
                                 beta_warmup = 4L))
  assert_true(!identical(rf$g_beta, g0), msg = "prior refit after the warm-up")
  err <- tryCatch(fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 2L, max_iter = 2L,
                                beta_prior_mode = "shared", verbose = FALSE),
                  error = function(e) e)
  assert_true(inherits(err, "error"), msg = "shared prior requires canonicalize")
})

# ---------------------------------------------------------------------------
# T10: supervised initialization starts program K along the Cox score
# direction: its first-iteration score correlates with the null martingale
# residuals far more than the unsupervised start does
# ---------------------------------------------------------------------------
run_test("Supervision-T10: init_supervised starts a program along the score direction", {
  args <- list(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 1L, prior_LF = "point_laplace",
               prior_beta = "normal", canonicalize = TRUE, verbose = FALSE)
  sup  <- do.call(fit_cox_on_yf, c(args, init_supervised = TRUE))
  uns  <- do.call(fit_cox_on_yf, args)
  u0   <- calc_cox_taylor_yf(rep(0, .sv$n), .sv$time, .sv$status)$u
  r_s  <- abs(cor(.sv$Y %*% sup$EF[, 3], u0))
  r_u  <- max(abs(suppressWarnings(cor(.sv$Y %*% uns$EF, u0))), na.rm = TRUE)
  assert_true(r_s > 0.5 && r_s > r_u + 0.1,
              msg = sprintf("supervised |cor| %.2f vs unsupervised max %.2f", r_s, r_u))
})

# ---------------------------------------------------------------------------
# T11: canonical fit stores -KL terms (<= 0) and its last monitored ELBO
# equals yf_canonical_elbo() at the returned moments (before Phase C)
# ---------------------------------------------------------------------------
run_test("Supervision-T11: canonical ELBO terms have the right sign and recompute", {
  fit <- fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 40L,
                       prior_LF = "point_laplace", canonicalize = TRUE,
                       sign_correction = FALSE, verbose = FALSE)
  assert_true(all(fit$neg_kl <= 1e-8), msg = "stored -KL terms are nonpositive")
  e <- yf_canonical_elbo(.sv$Y, .sv$time, .sv$status, fit$EL, fit$EL2, fit$EF, fit$EF2,
                         fit$EBeta, fit$EBeta2, sum(fit$neg_kl), fit$w_S)
  assert_near(e, tail(fit$history$elbo_full, 1), tol = 1e-6 * abs(e),
              msg = "monitored ELBO recomputes at the returned moments")
})

# ---------------------------------------------------------------------------
# T12: pruning removes programs that carry nothing and keeps a real one
# ---------------------------------------------------------------------------
run_test("Supervision-T12: prune removes empty programs and keeps the signal", {
  set.seed(21); n <- 80L; p <- 30L
  Ys <- scale(matrix(rnorm(n * 2), n, 2) %*% matrix(rnorm(2 * p), 2, p) +
              matrix(rnorm(n * p, sd = 0.7), n, p))
  fit <- fit_cox_on_yf(Ys, rexp(n), rbinom(n, 1, 0.7), K = 5L, max_iter = 200L,
                       prior_LF = "point_laplace", canonicalize = TRUE, prune = TRUE,
                       verbose = FALSE)
  assert_true(ncol(fit$EF) == 2L, msg = sprintf("2 programs kept (got %d)", ncol(fit$EF)))
  assert_true(nrow(fit$pruning) == 3L && all(fit$pruning$elbo_gain >= 0),
              msg = "three removals, each not lowering the ELBO")
})

# Protection must change the candidate set, including the literal 0 >= 0 boundary.
run_test("Supervision-T13: protected pruning excludes credible coefficients", {
  f <- fit_cox_on_yf(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 8L,
                     prior_LF = "point_laplace", canonicalize = TRUE,
                     verbose = FALSE, sign_correction = FALSE)
  f$EBeta <- c(2, 0.01, 0); f$EBeta2 <- f$EBeta^2 + c(0.01, 1, 0)
  f$neg_kl[1] <- -1e9  # removal is ELBO-favored despite credible survival activity
  raw <- prune_yf_nullcheck(f, .sv$Y, .sv$time, .sv$status)
  protected <- prune_yf_nullcheck(f, .sv$Y, .sv$time, .sv$status,
                                 protect_survival = TRUE)
  assert_true(raw[1] > 0, msg = "active program would be removed without protection")
  assert_true(all(is.finite(raw)), msg = "unprotected gains remain available")
  assert_true(all(protected[c(1, 3)] == -Inf), msg = "credible and zero-variance boundary protected")
  assert_near(protected[2], raw[2], tol = 0, msg = "eligible gain unchanged")
})

run_test("Supervision-T14: output convergence uses both fitted quantities", {
  a <- list(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 12L,
            prior_LF = "point_laplace", canonicalize = TRUE, verbose = FALSE)
  el <- do.call(fit_cox_on_yf, c(a, converge_on = "elbo"))
  ou <- do.call(fit_cox_on_yf, c(a, converge_on = "outputs"))
  assert_true(is.logical(ou$history$outputs_converged), msg = "final output criterion recorded")
  assert_true(identical(el$history$converged, el$history$elbo_converged),
              msg = "default stopping flag remains ELBO-based")
  assert_true(all(ou$history$output_change >= 0), msg = "finite output changes recorded")
  assert_true(!ou$history$converged || ou$history$outputs_converged,
              msg = "outputs stopping requires both quantities below threshold")
})

run_test("Supervision-T15: default fits and pruning match feature/supervised-F bit for bit", {
  old <- new.env(parent = globalenv())
  src <- system2("git", c("show", "feature/supervised-F:code/fit_cox_on_yf.R"), stdout = TRUE)
  tryCatch(eval(parse(text = src), old), error = function(e) {
    if (!grepl("!is.null(real_Y)", conditionMessage(e), fixed = TRUE)) stop(e)
  })
  for (pr in c(FALSE, TRUE)) {
    a <- list(.sv$Y, .sv$time, .sv$status, K = 3L, max_iter = 12L,
              prior_LF = "point_laplace", canonicalize = TRUE, verbose = FALSE, prune = pr)
    before <- do.call(old$fit_cox_on_yf, a)
    now <- do.call(fit_cox_on_yf, a)
    off <- do.call(fit_cox_on_yf, c(a, prune_protect_survival = FALSE))
    assert_true(identical(before, now) && identical(now, off),
                msg = "default object and explicit FALSE reproduce the baseline")
  }
})

run_test("Supervision-T16: output-change diagnostic agrees with two consecutive states", {
  a <- list(.sv$Y, .sv$time, .sv$status, K = 3L, prior_LF = "point_laplace",
            canonicalize = TRUE, verbose = FALSE, sign_correction = FALSE,
            converge_on = "elbo")
  one <- do.call(fit_cox_on_yf, c(a, max_iter = 1L))
  two <- do.call(fit_cox_on_yf, c(a, max_iter = 2L))
  expected <- max(multimodal_yfb_relative_change(one$EL %*% t(one$EF), two$EL %*% t(two$EF)),
    multimodal_yfb_relative_change(as.vector(.sv$Y %*% one$EF %*% one$EBeta),
                                   as.vector(.sv$Y %*% two$EF %*% two$EBeta)))
  assert_near(two$history$output_change[2], expected, tol = 1e-12,
              msg = "both fitted quantities enter the recorded maximum")
})
