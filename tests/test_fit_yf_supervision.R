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
