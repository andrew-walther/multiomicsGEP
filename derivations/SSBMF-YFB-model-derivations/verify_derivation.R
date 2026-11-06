# ============================================================
# Script: verify_derivation.R
# Purpose: Check the EBNM inputs (A, B) of SSBMF-YFB-model-derivation-full
#          against (1) a brute-force evaluation of the expected objective and
#          (2) the update functions that fit_ssbmf_yfb() calls.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, ebnm; code/fit_ssbmf_yfb.R and the helpers it sources
# Run: Rscript derivations/SSBMF-YFB-model-derivations/verify_derivation.R
#      (from the repository root). Exits with status 1 if any check fails.
# ============================================================
# Small synthetic fixtures only: no data are read and no files are written.
# Every posterior second moment carries a positive variance, so a formula that
# uses a squared mean in place of a second moment fails the checks.

suppressMessages(source("code/fit_ssbmf_yfb.R"))

# Assertion helper ----

results <- data.frame(check = character(), passed = logical(), discrepancy = numeric())

#' Record a numerical equality (scaled absolute discrepancy below tolerance)
#'
#' @param name Description of the identity.
#' @param actual Value computed one way.
#' @param expected Value computed independently.
#' @param tol Tolerance on max|actual - expected| / max(1, max|expected|).
#' @return Invisibly TRUE or FALSE; appends a row to `results`.
#' @examples
#' check_close("one", 1, 1)
check_close <- function(name, actual, expected, tol = 1e-8) {
  d <- max(abs(actual - expected)) / max(1, max(abs(expected)))
  ok <- length(actual) == length(expected) && all(is.finite(actual)) && d < tol
  results[nrow(results) + 1, ] <<- list(name, ok, d)
  invisible(ok)
}

# Fixture: a posterior state with nonzero variances ----

set.seed(20261106)
n <- 40; p <- 25; K <- 3; w_S <- 2.5
Y   <- scale(matrix(rnorm(n * p), n, p))
EL  <- matrix(rnorm(n * K), n, K);       EL2 <- EL^2 + matrix(runif(n * K, .05, .3), n, K)
EF  <- matrix(rnorm(p * K, sd = .3), p, K); EF2 <- EF^2 + matrix(runif(p * K, .01, .05), p, K)
EB  <- c(.4, -.3, .2);                    EB2 <- EB^2 + c(.02, .03, .05)
Tau <- runif(p, .5, 2)
time <- rexp(n); status <- rbinom(n, 1, .7)

## Survival working quantities at eta_hat = Y EF EB (Eqs. taylor, uW) ----
tay <- calc_cox_taylor_yf(as.vector(Y %*% EF %*% EB), time, status)
W <- tay$w; z <- as.vector(Y %*% EF %*% EB) + tay$u / W

#' Expected objective as a function of the posterior moments
#'
#' Genomics term (tau fixed, constants dropped) plus w_S times the expected
#' Taylor surrogate -1/2 sum_i W_i [(z_i - E eta_i)^2 + Var(eta_i)], with
#' E z_ik^2 = zbar_ik^2 + sum_j y_ij^2 Var(f_jk) (Eqs. gen-term, surv-expect).
#' @return Scalar objective (KL terms excluded: they depend only on q of the
#'   block being updated and are handled by EBNM).
#' @examples
#' expected_objective(EL, EL2, EF, EF2, EB, EB2)
expected_objective <- function(EL, EL2, EF, EF2, EB, EB2) {
  R2 <- (Y - EL %*% t(EF))^2 + EL2 %*% t(EF2) - EL^2 %*% t(EF^2)
  gen <- -0.5 * sum(sweep(R2, 2, Tau, "*"))
  Zb <- Y %*% EF
  EZ2 <- Zb^2 + Y^2 %*% (EF2 - EF^2)
  eta <- as.vector(Zb %*% EB)
  var_eta <- as.vector(EZ2 %*% EB2 - Zb^2 %*% EB^2)
  gen + w_S * (-0.5 * sum(W * ((z - eta)^2 + var_eta)))
}

#' Read off (A, B) by brute force: the objective as a function of one
#' coordinate's posterior moments (m, s2 = E[x^2]) has the form
#' -1/2 A s2 + B m + C. Evaluate it at three moment pairs and solve.
#' @param set_moments function(m, s2) returning the objective.
#' @return c(A = , B = ).
#' @examples
#' brute_AB(function(m, s2) -0.5 * 2 * s2 + 3 * m)
brute_AB <- function(set_moments) {
  f00 <- set_moments(0, 0); f01 <- set_moments(0, 1); f11 <- set_moments(1, 1)
  A <- -2 * (f01 - f00)
  B <- f11 - f01
  # Quadratic-form check: a fourth point must agree with -1/2 A s2 + B m + C
  f_test <- set_moments(0.7, 1.3)
  stopifnot(abs(f_test - (f00 - 0.5 * A * 1.3 + B * 0.7)) < 1e-6 * max(1, abs(f_test)))
  c(A = A, B = B)
}

# 1. q(beta_k): Eq. AB-beta ----

for (k in 1:K) {
  Zb <- Y %*% EF; VZ <- Y^2 %*% (EF2 - EF^2)
  z_mk <- z - as.vector(Zb[, -k, drop = FALSE] %*% EB[-k])
  A_doc <- w_S * sum(W * (Zb[, k]^2 + VZ[, k]))
  B_doc <- w_S * sum(W * z_mk * Zb[, k])
  bf <- brute_AB(function(m, s2) { b <- EB; b2 <- EB2; b[k] <- m; b2[k] <- s2
                                   expected_objective(EL, EL2, EF, EF2, b, b2) })
  check_close(sprintf("beta_%d: document A matches the expected objective", k), A_doc, bf[["A"]])
  check_close(sprintf("beta_%d: document B matches the expected objective", k), B_doc, bf[["B"]])
  code <- update_beta_k(W, compute_z_no_k(z, Zb, EB, k), Zb[, k], Zb[, k]^2 + VZ[, k],
                        prior_family = "normal", alpha = w_S, survival_divisor = 1)
  check_close(sprintf("beta_%d: update_beta_k() A", k), code$A, A_doc)
  check_close(sprintf("beta_%d: update_beta_k() B", k), code$B, B_doc)
}

# 2. q(l_k): Eq. AB-L (genomics only) ----

for (k in 1:K) {
  R_k <- Y - EL[, -k, drop = FALSE] %*% t(EF[, -k, drop = FALSE])
  A_doc <- sum(Tau * EF2[, k]); B_doc <- as.vector(R_k %*% (Tau * EF[, k]))
  for (i in c(1, 17, n)) {
    bf <- brute_AB(function(m, s2) { L <- EL; L2 <- EL2; L[i, k] <- m; L2[i, k] <- s2
                                     expected_objective(L, L2, EF, EF2, EB, EB2) })
    check_close(sprintf("l_%d,%d: document A and B match the expected objective", i, k),
                c(A_doc, B_doc[i]), unname(bf))
  }
  code <- update_L_surv_YFB_k(Tau, EF[, k], EF2[, k], R_k, prior_family = "point_laplace")
  check_close(sprintf("l_%d: update_L_surv_YFB_k() A and B", k), c(code$A, code$B), c(A_doc, B_doc))
}

# 3. q(f_jk): Eqs. AB-F-gen, AB-F-surv, AB-F ----

#' Document's (A, B) for gene j of program k at the current state
#' @return c(A, B).
#' @examples
#' doc_AB_F(EL, EL2, EF, EF2, 1, 1)
doc_AB_F <- function(EL, EL2, EF, EF2, j, k) {
  R_k <- Y - EL[, -k, drop = FALSE] %*% t(EF[, -k, drop = FALSE])
  Zb <- Y %*% EF
  z_mk <- z - as.vector(Zb[, -k, drop = FALSE] %*% EB[-k])
  y <- Y[, j]
  A <- Tau[j] * sum(EL2[, k]) + w_S * EB2[k] * sum(W * y^2)
  B <- Tau[j] * sum(R_k[, j] * EL[, k]) +
       w_S * (EB[k] * sum(W * y * z_mk) - EB2[k] * sum(W * y * (Zb[, k] - y * EF[j, k])))
  c(A, B)
}

for (k in 1:K) for (j in c(1, 9, p)) {
  bf <- brute_AB(function(m, s2) { Fm <- EF; F2 <- EF2; Fm[j, k] <- m; F2[j, k] <- s2
                                   expected_objective(EL, EL2, Fm, F2, EB, EB2) })
  check_close(sprintf("f_%d,%d: document A and B match the expected objective", j, k),
              doc_AB_F(EL, EL2, EF, EF2, j, k), unname(bf))
}

## Sequential sweep: the code's per-gene (A, B) equal the document's at each step ----
for (k in 1:K) {
  R_k <- Y - EL[, -k, drop = FALSE] %*% t(EF[, -k, drop = FALSE])
  Zb <- Y %*% EF; VZ_k <- as.vector(Y^2 %*% (EF2[, k] - EF[, k]^2))
  z_mk <- z - as.vector(Zb[, -k, drop = FALSE] %*% EB[-k])
  A_gen <- Tau * sum(EL2[, k]); B_gen <- Tau * as.vector(crossprod(R_k, EL[, k]))
  prior <- multimodal_yfb_fit_signed_prior(B_gen / A_gen, 1 / A_gen, "point_laplace")
  upd <- multimodal_yfb_update_F_mk(Y, Tau, EL[, k], EL2[, k], R_k, w_S * W, w_S * W * z_mk,
                                    Zb[, k], VZ_k, EF[, k], EF2[, k], EB[k], EB2[k], prior,
                                    prior_update = "ebnm")
  Fs <- EF; F2s <- EF2; AB <- matrix(NA, p, 2)
  for (j in 1:p) {           # replay the sweep with the document's formulas
    AB[j, ] <- doc_AB_F(EL, EL2, Fs, F2s, j, k)
    Fs[j, k] <- upd$details$mean[j]; F2s[j, k] <- upd$details$second[j]
  }
  check_close(sprintf("f_%d: sequential sweep A (all genes)", k), upd$details$A, AB[, 1])
  check_close(sprintf("f_%d: sequential sweep B (all genes)", k), upd$details$B, AB[, 2])
}

## Negative controls: the April 2026 formulas must NOT match the expected objective ----
k <- 1; j <- 9
Zb <- Y %*% EF; z_mk <- z - as.vector(Zb[, -k] %*% EB[-k])
bfF <- brute_AB(function(m, s2) { Fm <- EF; F2 <- EF2; Fm[j, k] <- m; F2[j, k] <- s2
                                  expected_objective(EL, EL2, Fm, F2, EB, EB2) })
R_k <- Y - EL[, -k] %*% t(EF[, -k])
B_april <- Tau[j] * sum(R_k[, j] * EL[, k]) + w_S * EB[k] * sum(W * Y[, j] * z_mk)  # erratum Y1
results[nrow(results) + 1, ] <- list("negative control: April B_surv (Y1) disagrees",
  abs(B_april - bfF[["B"]]) > 1e-3 * max(1, abs(bfF[["B"]])), abs(B_april - bfF[["B"]]))
bfB <- brute_AB(function(m, s2) { b <- EB; b2 <- EB2; b[k] <- m; b2[k] <- s2
                                  expected_objective(EL, EL2, EF, EF2, b, b2) })
A_april <- w_S * sum(W * Zb[, k]^2)                                                # erratum Y5
results[nrow(results) + 1, ] <- list("negative control: April A_beta without Var(z) (Y5) disagrees",
  abs(A_april - bfB[["A"]]) > 1e-3 * max(1, abs(bfB[["A"]])), abs(A_april - bfB[["A"]]))

# 4. tau: Eq. tau ----

R2 <- (Y - EL %*% t(EF))^2 + EL2 %*% t(EF2) - EL^2 %*% t(EF^2)
check_close("tau: update_tau() equals n / sum_i Rbar^2", update_tau(Y, EL, EL2, EF, EF2)$Tau,
            n / colSums(R2))

# 5. Identities used in the text ----

tay0 <- calc_cox_taylor_yf(rnorm(n, sd = .5), rexp(n), rbinom(n, 1, .6))
check_close("sum of Cox working weights equals the number of events (no ties, no floor)",
            sum(tay0$w), sum(tay0$u + tay0$w))   # sum u = 0  <=>  sum w = sum delta

cz <- canonicalize_yf_factors(Y, EL, EL2, EF, EF2, EB, EB2, K)
check_close("score standardization leaves L F' unchanged", cz$EL %*% t(cz$EF), EL %*% t(EF))
check_close("score standardization leaves eta unchanged", Y %*% cz$EF %*% cz$EBeta, Y %*% EF %*% EB)
check_close("score standardization gives sd(Y F_k) = 1", apply(Y %*% cz$EF, 2, sd), rep(1, K))
rho <- function(EL2, EF, EB2) outer(1 / Tau, w_S * EB2 / colSums(EL2)) * as.vector(crossprod(Y^2, W))
check_close("precision share rho is invariant to the rescaling", rho(cz$EL2, cz$EF, cz$EBeta2),
            rho(EL2, EF, EB2))

# Report ----

print(results, row.names = FALSE, digits = 3)
cat(sprintf("\n%d of %d checks passed.\n", sum(results$passed), nrow(results)))
if (!all(results$passed)) quit(status = 1)
