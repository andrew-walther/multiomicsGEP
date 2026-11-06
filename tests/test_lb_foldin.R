# ============================================================
# Script:  test_lb_foldin.R
# Purpose: Tests for LB fold-in prediction: fit_ssbmf_lb() stores the fitted
#          loading priors, and predict_ssbmf_lb_foldin() computes new-patient
#          loadings as the training L update without the survival term.
# Author: Andrew Walther
# Created: 2026-11-06
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_lb_foldin.R  (standalone)
# Dependencies: ebnm, survival
# ============================================================
#
# Fold-in (Yusha's note, docs/notes/2026-10-30-yusha-MF_NewPatient_FoldIn.pdf):
# for a new patient, q(l*_k) = EBNM(x = B*_k / A*_k, s = 1 / sqrt(A*_k)) under
# the frozen training prior g_{L,k}, where
#   A*_k = (1 - alpha) sum_j tau_j E[f_jk^2]
#   B*_k = (1 - alpha) sum_j tau_j fbar_jk (y*_j - sum_{m != k} l*_m fbar_jm),
# i.e. the training update with the survival term removed. Sweeps over k
# repeat until the loadings stop changing.

if (!exists("run_test")) source("tests/test_helpers.R")
suppressPackageStartupMessages({
  library(ebnm); library(survival)
})
if (!exists("fit_ssbmf_lb")) {
  suppressMessages(tryCatch(source("code/fit_ssbmf_lb.R"),
                            error = function(e) invisible(NULL)))
}
if (!exists("predict_ssbmf_lb_foldin")) source("code/predict_ssbmf_lb.R")

cat("\n--- test_lb_foldin.R ---\n")

## Simulated LB data with a planted survival signal ----
.fi_sim <- function(seed, n = 120, p = 100, K = 3, beta = c(1, -0.8, 0)) {
  set.seed(seed)
  L  <- matrix(rexp(n * K), n, K)
  Fm <- matrix(0, p, K)
  for (k in seq_len(K)) Fm[sample.int(p, 15), k] <- 1
  Y   <- L %*% t(Fm) + matrix(rnorm(n * p, sd = 0.5), n, p)
  eta <- as.vector(L %*% beta)
  ev  <- -log(runif(n)) / (0.05 * exp(eta))
  cen <- rexp(n, 1) * median(ev) * 2
  list(Y = Y, time = pmin(ev, cen), status = as.integer(ev <= cen))
}
.fi_fit <- function(d, K = 3, alpha = 0.5, ...) {
  suppressWarnings(fit_ssbmf_lb(d$Y, d$time, d$status, K = K, alpha = alpha,
                                verbose = FALSE, ...))
}

.d      <- .fi_sim(1)
.fit_05 <- .fi_fit(.d, alpha = 0.5, max_iter = 60, tol = 1e-4)

run_test("FI-T1: fit stores the fitted loading priors, prior family and final alpha", {
  fi <- .fit_05$foldin
  assert_true(is.list(fi), "fit$foldin exists")
  assert_equal(length(fi$g_L), 3L, "one fitted prior per factor")
  assert_true(all(vapply(fi$g_L, function(g) inherits(g, "gammamix"), logical(1))),
              "point_exponential priors are stored as ebnm fitted_g objects")
  assert_equal(fi$prior_LF, "point_exponential", "prior family recorded")
  assert_near(fi$alpha, 0.5, tol = 0, msg = "alpha of the final L update recorded")
})

run_test("FI-T2: at alpha = 0, fold-in of the training data reproduces the fitted loadings", {
  # With alpha = 0 the training L update is genomics only, so the converged
  # EL is a fixed point of the fold-in map (up to the CAVI tolerance).
  fit0 <- .fi_fit(.d, alpha = 0, max_iter = 500, tol = 1e-8)
  fi   <- predict_ssbmf_lb_foldin(.d$Y, fit0, tol = 1e-10)
  rel  <- max(abs(fi$L_test - fit0$EL)) / max(abs(fit0$EL))
  assert_true(rel < 1e-3, sprintf("max relative difference %.2e", rel))
})

run_test("FI-T3: a one-factor fold-in equals the tempered EBNM update computed by hand", {
  # One live factor of the K = 3 fit, as a one-factor fit object (fold-in
  # depends only on EF, EF2, Tau, EBeta and $foldin). LB fits with K = 1
  # collapse from the pmax(SVD) start, so a fitted K = 1 fixture is not used.
  k    <- which.max(colSums(.fit_05$EF2))
  fit1 <- list(EF = .fit_05$EF[, k, drop = FALSE], EF2 = .fit_05$EF2[, k, drop = FALSE],
               Tau = .fit_05$Tau, EBeta = .fit_05$EBeta[k],
               foldin = list(g_L = .fit_05$foldin$g_L[k], prior_LF = "point_exponential",
                             alpha = 0.5))
  assert_true(sum(fit1$EF2[, 1]) > 0, "fixture factor is live")
  Ynew <- .d$Y[1:20, ]
  A <- (1 - 0.5) * sum(fit1$Tau * fit1$EF2[, 1])
  B <- (1 - 0.5) * as.vector(Ynew %*% (fit1$Tau * fit1$EF[, 1]))
  ref <- ebnm(x = B / A, s = rep(1 / sqrt(A), 20), prior_family = "point_exponential",
              g_init = fit1$foldin$g_L[[1]], fix_g = TRUE)
  fi <- predict_ssbmf_lb_foldin(Ynew, fit1)
  assert_near(fi$L_test[, 1], ref$posterior$mean, tol = 1e-10, msg = "posterior mean")
  assert_near(fi$L2_test[, 1], ref$posterior$sd^2 + ref$posterior$mean^2, tol = 1e-10,
              msg = "posterior second moment")
})

run_test("FI-T4: risk scores are L_test %*% EBeta and moments are consistent", {
  fi <- predict_ssbmf_lb_foldin(.d$Y[1:30, ], .fit_05)
  assert_equal(dim(fi$L_test), c(30L, 3L), "n_test x K loadings")
  assert_near(fi$risk_scores, as.vector(fi$L_test %*% .fit_05$EBeta), tol = 1e-12,
              msg = "risk = L beta")
  assert_true(all(fi$L2_test >= fi$L_test^2 - 1e-12), "E[l^2] >= E[l]^2")
  assert_true(isTRUE(fi$converged), "converged at the default tolerance")
})

run_test("FI-T5: fold-in fails loudly on unsupported inputs", {
  err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
  assert_true(err(predict_ssbmf_lb_foldin(.d$Y[, -1], .fit_05)), "gene-count mismatch")
  old <- .fit_05; old$foldin <- NULL
  assert_true(err(predict_ssbmf_lb_foldin(.d$Y, old)), "fit without stored priors")
  coh <- .fit_05; coh$cohort_id <- rep(1:2, length.out = nrow(.d$Y))
  assert_true(err(predict_ssbmf_lb_foldin(.d$Y, coh)), "cohort-offset fits are not supported")
  Yna <- .d$Y; Yna[1, 1] <- NA
  assert_true(err(predict_ssbmf_lb_foldin(Yna, .fit_05)), "missing values")
  a1 <- .fit_05; a1$foldin$alpha <- 1
  assert_true(err(predict_ssbmf_lb_foldin(.d$Y, a1)), "alpha = 1 (no genomics information)")
  nm <- .fit_05; rownames(nm$EF) <- paste0("g", seq_len(nrow(nm$EF)))
  Ynm <- .d$Y; colnames(Ynm) <- rev(rownames(nm$EF))
  assert_true(err(predict_ssbmf_lb_foldin(Ynm, nm)), "gene order mismatch by name")
})

run_test("FI-T6: hitting max_sweeps warns instead of returning silently", {
  warned <- FALSE
  withCallingHandlers(predict_ssbmf_lb_foldin(.d$Y, .fit_05, max_sweeps = 1),
                      warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  assert_true(warned, "non-convergence is reported")
})
