# ============================================================
# Script: verify_derivation.R
# Purpose: Check the multimodal YFB derivation using derivatives,
#          exact finite-support expectations, and numerical integrals.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-09-10
# Dependencies: base R (stats functions ship with R)
# ============================================================
# These are mathematical fixtures, not a simulation study or a model fit.
# No production module is sourced and no data or result files are written.
# Second moments include positive variances so plug-in errors are detectable.

# Assertion and differentiation helpers ----

check_results <- data.frame(check = character(), passed = logical())

#' Record a mathematical equality or stop with its numerical discrepancy
#' @param name Character description of the identity.
#' @param actual Numeric result being checked.
#' @param expected Independently computed reference.
#' @param tolerance Maximum scaled absolute discrepancy.
#' @return Invisibly TRUE; appends one row to check_results.
#' @examples
#' check_close("addition", 1 + 1, 2)
check_close <- function(name, actual, expected, tolerance = 1e-7) {
  if (length(actual) != length(expected) ||
      any(!is.finite(actual)) || any(!is.finite(expected))) {
    stop(name, ": nonfinite values or incompatible lengths.")
  }
  discrepancy <- max(abs(actual - expected)) / max(1, abs(expected))
  if (discrepancy > tolerance) {
    stop(name, ": scaled error ", format(discrepancy, digits = 8),
         " exceeds ", tolerance, ".")
  }
  check_results <<- rbind(check_results,
                         data.frame(check = name, passed = TRUE))
  invisible(TRUE)
}

#' Approximate a scalar function's gradient by central differences
#' @param fn Function accepting a numeric vector and returning one number.
#' @param x Evaluation point.
#' @param step Perturbation size.
#' @return Numeric gradient with length(x) entries.
#' @examples
#' fd_gradient(function(x) sum(x^2), c(1, 2))
fd_gradient <- function(fn, x, step = 1e-5) {
  vapply(seq_along(x), function(j) {
    displacement <- rep(0, length(x))
    displacement[j] <- step
    (fn(x + displacement) - fn(x - displacement)) / (2 * step)
  }, numeric(1))
}

#' Approximate a scalar function's Hessian by central differences
#' @param fn Function accepting a numeric vector and returning one number.
#' @param x Evaluation point.
#' @param step Perturbation size.
#' @return Symmetric numeric Hessian.
#' @examples
#' fd_hessian(function(x) sum(x^2), c(1, 2))
fd_hessian <- function(fn, x, step = 1e-4) {
  out <- matrix(0, length(x), length(x))
  for (j in seq_along(x)) {
    ej <- rep(0, length(x))
    ej[j] <- step
    for (h in seq_along(x)) {
      eh <- rep(0, length(x))
      eh[h] <- step
      out[j, h] <- (fn(x + ej + eh) - fn(x + ej - eh) -
                   fn(x - ej + eh) + fn(x - ej - eh)) / (4 * step^2)
    }
  }
  out
}

# Cox likelihood and independent derivative check ----

#' Evaluate a tiny Breslow partial likelihood and analytic derivatives
#' @param eta Numeric risk predictor.
#' @param time Follow-up times.
#' @param status Event indicators.
#' @return List with log likelihood, score, negative Hessian and work matrix.
#' @examples
#' cox_fixture(c(0.1, -0.2, 0.3), c(1, 2, 3), c(1, 1, 0))
cox_fixture <- function(eta, time, status) {
  score <- as.numeric(status)
  hessian <- work <- matrix(0, length(eta), length(eta))
  value <- sum(status * eta)
  for (event_time in sort(unique(time[status == 1]))) {
    event_count <- sum(status[time == event_time])
    risk <- time >= event_time
    probability <- as.numeric(risk) * exp(eta)
    denominator <- sum(probability)
    probability <- probability / denominator
    value <- value - event_count * log(denominator)
    score <- score - event_count * probability
    hessian <- hessian + event_count *
      (diag(probability) - tcrossprod(probability))
    work <- work + event_count * diag(probability)
  }
  list(value = value, score = score, hessian = hessian, work = work)
}

time <- c(1, 2, 2, 4)
status <- c(1, 1, 1, 0)
eta0 <- c(0.2, -0.3, 0.4, -0.1)
cox <- cox_fixture(eta0, time, status)
#' Fixed-data Cox log partial likelihood for derivative checks
#' @param eta Numeric four-subject predictor.
#' @return Scalar log partial likelihood.
#' @examples
#' cox_value(rep(0, 4))
cox_value <- function(eta) cox_fixture(eta, time, status)$value
check_close("Breslow score with tied events",
            cox$score, fd_gradient(cox_value, eta0))
check_close("Full Cox Hessian with tied events",
            cox$hessian, -fd_hessian(cox_value, eta0), 1e-6)
check_close("Cox Hessian annihilates a constant shift",
            as.vector(cox$hessian %*% rep(1, 4)), rep(0, 4))
check_close("Working curvature is locally at least the Hessian",
            min(eigen(cox$work - cox$hessian,
                      symmetric = TRUE)$values) >= -1e-10, TRUE)
check_close("Working curvature is not the exact Hessian diagonal",
            max(abs(diag(cox$work) - diag(cox$hessian))) > 0.01, TRUE)

# Expected quadratic objective ----

n <- 4L
p <- 3L
K <- 2L
Y <- cbind(c(1, 2, 1.5, 0.5),
           c(1.2, 1.8, 1.7, 0.4),
           c(-0.4, 0.2, 0.8, 1.0))
modality <- c(1L, 2L, 2L)
likelihood_weight <- c(0.8, 0.35)[modality]
projection_weight <- c(1.0, 0.7)[modality]
X <- sweep(Y, 2, projection_weight, "*")
tau <- c(1.2, 0.7, 1.8)
survival_weight <- 0.9
L <- matrix(c(0.3, 0.8, 0.6, 0.2, 0.7, 0.4, 0.1, 0.5), n, K)
F <- matrix(c(0.5, 0.25, -0.15, 0.2, 0.4, 0.3), p, K)
beta <- c(0.6, -0.35)
vL <- matrix(seq(0.02, 0.09, length.out = n * K), n, K)
vF <- matrix(seq(0.01, 0.06, length.out = p * K), p, K)
vBeta <- c(0.07, 0.04)

#' Compute risk mean and covariance under independent factor posteriors
#' @param f Loading posterior means, p by K.
#' @param beta_mean Coefficient posterior means.
#' @return List of mean, covariance, projected means and score covariances.
#' @examples
#' risk_moments(F, beta)
risk_moments <- function(f, beta_mean) {
  projected <- X %*% f
  covariance <- matrix(0, n, n)
  score_covariance <- vector("list", K)
  for (k in seq_len(K)) {
    score_covariance[[k]] <- tcrossprod(sweep(X, 2, vF[, k], "*"), X)
    covariance <- covariance +
      (beta_mean[k]^2 + vBeta[k]) * score_covariance[[k]] +
      vBeta[k] * tcrossprod(projected[, k])
  }
  list(mean = as.vector(projected %*% beta_mean),
       covariance = covariance, projected = projected,
       score_covariance = score_covariance)
}

#' Evaluate expected omics and survival quadratic likelihood terms
#' @param l Subject posterior means.
#' @param f Loading posterior means.
#' @param beta_mean Coefficient posterior means.
#' @param precision Feature precisions.
#' @param W Fixed Cox curvature matrix.
#' @param tstar Fixed Cox linear coefficient u + W eta0.
#' @return Scalar expected likelihood, without coordinate-constant terms.
#' @examples
#' expected_objective(L, F, beta, tau, cox$work,
#'                    cox$score + cox$work %*% eta0)
expected_objective <- function(l, f, beta_mean, precision, W, tstar) {
  residual_second <- (Y - tcrossprod(l, f))^2 +
    tcrossprod(l^2 + vL, f^2 + vF) - tcrossprod(l^2, f^2)
  omics <- sum(likelihood_weight *
                 (n / 2 * log(precision) -
                    precision / 2 * colSums(residual_second)))
  risk <- risk_moments(f, beta_mean)
  survival <- sum(tstar * risk$mean) -
    sum(risk$mean * (W %*% risk$mean)) / 2 -
    sum(W * t(risk$covariance)) / 2
  omics + survival_weight * survival
}

for (curvature_name in c("full Hessian", "working diagonal")) {
  W <- if (curvature_name == "full Hessian") cox$hessian else cox$work
  tstar <- as.vector(cox$score + W %*% eta0)
  k <- 1L
  risk <- risk_moments(F, beta)
  residual <- Y - tcrossprod(L[, -k, drop = FALSE],
                            F[, -k, drop = FALSE])
  eta_other <- as.vector(risk$projected[, -k, drop = FALSE] %*% beta[-k])
  h_k <- as.vector(tstar - W %*% eta_other)
  weighted_tau <- likelihood_weight * tau

  A_L <- sum(weighted_tau * (F[, k]^2 + vF[, k]))
  B_L <- as.vector(residual %*% (weighted_tau * F[, k]))
  #' Vary one score column while holding other fixture moments fixed
  #' @param column Numeric n-vector of posterior means.
  #' @return Expected quadratic likelihood at this column.
  #' @examples
  #' l_fn(L[, k])
  l_fn <- function(column) {
    l_new <- L
    l_new[, k] <- column
    expected_objective(l_new, F, beta, tau, W, tstar)
  }
  check_close(paste("L gradient:", curvature_name),
              fd_gradient(l_fn, L[, k]), B_L - A_L * L[, k])
  check_close(paste("L curvature:", curvature_name),
              -diag(fd_hessian(l_fn, L[, k])), rep(A_L, n), 1e-6)

  score_second <- sum(risk$projected[, k] *
                        (W %*% risk$projected[, k])) +
    sum(W * t(risk$score_covariance[[k]]))
  A_beta <- survival_weight * score_second
  B_beta <- survival_weight * sum(risk$projected[, k] * h_k)
  #' Vary one survival coefficient mean with its variance fixed
  #' @param coefficient Numeric scalar posterior mean.
  #' @return Expected quadratic likelihood at this coefficient.
  #' @examples
  #' beta_fn(beta[k])
  beta_fn <- function(coefficient) {
    beta_new <- beta
    beta_new[k] <- coefficient
    expected_objective(L, F, beta_new, tau, W, tstar)
  }
  check_close(paste("Beta gradient with F uncertainty:", curvature_name),
              fd_gradient(beta_fn, beta[k]), B_beta - A_beta * beta[k])
  check_close(paste("Beta curvature with F uncertainty:", curvature_name),
              -fd_hessian(beta_fn, beta[k]), A_beta, 1e-6)

  A_F <- B_F <- omitted_interaction <- numeric(p)
  for (a in seq_len(p)) {
    z_other_feature <- risk$projected[, k] - X[, a] * F[a, k]
    A_F[a] <- weighted_tau[a] * sum(L[, k]^2 + vL[, k]) +
      survival_weight * (beta[k]^2 + vBeta[k]) *
      sum(X[, a] * (W %*% X[, a]))
    omitted_interaction[a] <- survival_weight * (beta[k]^2 + vBeta[k]) *
      sum(X[, a] * (W %*% z_other_feature))
    B_F[a] <- weighted_tau[a] * sum(residual[, a] * L[, k]) +
      survival_weight * beta[k] * sum(X[, a] * h_k) -
      omitted_interaction[a]
  }
  #' Vary one loading column while retaining positive posterior variances
  #' @param column Numeric p-vector of loading posterior means.
  #' @return Expected quadratic likelihood at this column.
  #' @examples
  #' f_fn(F[, k])
  f_fn <- function(column) {
    f_new <- F
    f_new[, k] <- column
    expected_objective(L, f_new, beta, tau, W, tstar)
  }
  check_close(paste("F gradient with same-factor terms:", curvature_name),
              fd_gradient(f_fn, F[, k]), B_F - A_F * F[, k])
  check_close(paste("F diagonal curvature:", curvature_name),
              -diag(fd_hessian(f_fn, F[, k])), A_F, 1e-6)
  check_close(paste("Dropping feature interaction is detected:", curvature_name),
              max(abs(omitted_interaction)) > 0.01, TRUE)
}

# Exact expectations from independent two-point distributions ----

# The quadratic moment identities do not require Gaussian q. Enumerating
# independent +/- deviations produces the specified means and variances
# exactly, without Monte Carlo tolerances or sampling error.
signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), p * K + K)))
eta_draws <- t(apply(signs, 1, function(row) {
  f_draw <- F + matrix(row[seq_len(p * K)], p, K) * sqrt(vF)
  beta_draw <- beta + row[p * K + seq_len(K)] * sqrt(vBeta)
  as.vector(X %*% f_draw %*% beta_draw)
}))
risk <- risk_moments(F, beta)
enumerated_mean <- colMeans(eta_draws)
centered_eta <- sweep(eta_draws, 2, enumerated_mean, "-")
enumerated_covariance <- crossprod(centered_eta) / nrow(eta_draws)
check_close("Risk mean from exact enumeration", risk$mean, enumerated_mean)
check_close("Risk covariance includes shared-F cross-subject uncertainty",
            risk$covariance, enumerated_covariance)
W <- cox$hessian
tstar <- as.vector(cox$score + W %*% eta0)
check_close("Expected full quadratic from exact enumeration",
            mean(as.vector(eta_draws %*% tstar) -
                   rowSums((eta_draws %*% W) * eta_draws) / 2),
            sum(tstar * risk$mean) -
              sum(risk$mean * (W %*% risk$mean)) / 2 -
              sum(W * t(risk$covariance)) / 2)

i <- 2L
a <- 1L
factor_signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), 2 * K)))
squared_residuals <- apply(factor_signs, 1, function(row) {
  l_draw <- L[i, ] + row[seq_len(K)] * sqrt(vL[i, ])
  f_draw <- F[a, ] + row[K + seq_len(K)] * sqrt(vF[a, ])
  (Y[i, a] - sum(l_draw * f_draw))^2
})
residual_second <- (Y - tcrossprod(L, F))^2 +
  tcrossprod(L^2 + vL, F^2 + vF) - tcrossprod(L^2, F^2)
check_close("Expected squared residual from exact enumeration",
            mean(squared_residuals), residual_second[i, a])
tau_mle <- n / colSums(residual_second)
#' Vary the point-estimated feature precisions with q fixed
#' @param precision Positive numeric p-vector.
#' @return Expected likelihood at these precisions.
#' @examples
#' tau_fn(tau)
tau_fn <- function(precision) {
  expected_objective(L, F, beta, precision, W, tstar)
}
check_close("Weighted Tau score is zero at n / expected SSE",
            fd_gradient(tau_fn, tau_mle), rep(0, p))

# Prior moments by independent numerical integration ----

x <- 0.65
noise_sd <- 0.7
pi_slab <- 0.6
slab_sd <- 1.1
spike_density <- dnorm(x, 0, noise_sd)
normal_slab_density <- dnorm(x, 0, sqrt(noise_sd^2 + slab_sd^2))
responsibility <- pi_slab * normal_slab_density /
  ((1 - pi_slab) * spike_density + pi_slab * normal_slab_density)
posterior_variance <- 1 / (1 / noise_sd^2 + 1 / slab_sd^2)
posterior_mean <- posterior_variance * x / noise_sd^2
normal_moments <- responsibility *
  c(posterior_mean, posterior_variance + posterior_mean^2)
normal_integrals <- vapply(0:2, function(power) {
  integrate(function(theta) theta^power *
              dnorm(x, theta, noise_sd) * dnorm(theta, 0, slab_sd),
            -Inf, Inf, rel.tol = 1e-10)$value
}, numeric(1))
normal_denominator <- (1 - pi_slab) * spike_density +
  pi_slab * normal_integrals[1]
check_close("Point-normal posterior moments by integration",
            normal_moments, pi_slab * normal_integrals[2:3] / normal_denominator)

rate <- 1.3
truncated_mean <- x - rate * noise_sd^2
h <- truncated_mean / noise_sd
mills <- dnorm(h) / pnorm(h)
exponential_slab_density <- rate *
  exp(-rate * x + rate^2 * noise_sd^2 / 2) * pnorm(h)
responsibility_exp <- pi_slab * exponential_slab_density /
  ((1 - pi_slab) * spike_density + pi_slab * exponential_slab_density)
exp_moments <- responsibility_exp *
  c(truncated_mean + noise_sd * mills,
    truncated_mean^2 + noise_sd^2 + truncated_mean * noise_sd * mills)
exp_integrals <- vapply(0:2, function(power) {
  integrate(function(theta) theta^power *
              dnorm(x, theta, noise_sd) * dexp(theta, rate),
            0, Inf, rel.tol = 1e-10)$value
}, numeric(1))
exp_denominator <- (1 - pi_slab) * spike_density +
  pi_slab * exp_integrals[1]
check_close("Point-exponential marginal density by integration",
            exponential_slab_density, exp_integrals[1])
check_close("Point-exponential posterior moments by integration",
            exp_moments, pi_slab * exp_integrals[2:3] / exp_denominator)

responsibilities <- c(0.25, 0.7, 0.9)
first_moments <- c(0.1, 0.5, 0.8)
second_moments <- c(0.08, 0.6, 1.1)
pi_new <- mean(responsibilities)
rate_new <- sum(responsibilities) / sum(first_moments)
variance_new <- sum(second_moments) / sum(responsibilities)
prior_derivatives <- c(
  sum(responsibilities) / pi_new -
    sum(1 - responsibilities) / (1 - pi_new),
  sum(responsibilities) / rate_new - sum(first_moments),
  -sum(responsibilities) / (2 * variance_new) +
    sum(second_moments) / (2 * variance_new^2))
check_close("Prior EM mixing, rate, and variance scores", prior_derivatives,
            rep(0, 3))

# The log-normalizer identity must agree with a direct Gaussian KL.
pseudo_A <- 2.4
pseudo_B <- 0.8
prior_variance <- 1.21
q_variance <- 1 / (pseudo_A + 1 / prior_variance)
q_mean <- q_variance * pseudo_B
pseudo_x <- pseudo_B / pseudo_A
pseudo_variance <- 1 / pseudo_A
log_normalizer <- dnorm(pseudo_x, 0,
                        sqrt(pseudo_variance + prior_variance), log = TRUE)
expected_pseudo <- -log(2 * pi * pseudo_variance) / 2 -
  (pseudo_x^2 - 2 * pseudo_x * q_mean + q_variance + q_mean^2) /
  (2 * pseudo_variance)
direct_negative_kl <- -(log(prior_variance / q_variance) +
                          (q_variance + q_mean^2) / prior_variance - 1) / 2
check_close("EBNM log-normalizer identity equals direct Gaussian KL",
            log_normalizer - expected_pseudo, direct_negative_kl)

# Structural reductions and projection derivatives ----

F_second <- F[, 1]^2 + vF[, 1]
partitioned_precision <- sum(vapply(split(seq_len(p), modality), function(idx) {
  sum(tau[idx] * F_second[idx])
}, numeric(1)))
check_close("Unweighted modality sum equals stacked score precision",
            partitioned_precision, sum(tau * F_second))
check_close("Replicated features preserve inverse-count weighted score sum",
            sum(rep(tau, 3) * rep(F_second, 3)) / (3 * p),
            sum(tau * F_second) / p)

feature_scale <- c(2, 3, 3)
scaled_Y <- sweep(Y, 2, feature_scale, "*")
scaled_F <- sweep(F, 1, feature_scale, "*")
scaled_tau <- tau / feature_scale^2
check_close("Feature rescaling preserves precision-weighted mean residual",
            sweep((scaled_Y - tcrossprod(L, scaled_F))^2, 2, scaled_tau, "*"),
            sweep((Y - tcrossprod(L, F))^2, 2, tau, "*"))
check_close("Raw YFB is not invariant to modality rescaling",
            max(abs(scaled_Y %*% scaled_F - Y %*% F)) > 0.1, TRUE)

column <- F[, 1]
epsilon <- 0.03
norm <- sqrt(sum(column^2) + epsilon)
j <- 2L
first_derivative <- X[, j] / norm -
  as.vector(X %*% column) * column[j] / norm^3
second_derivative <- -2 * X[, j] * column[j] / norm^3 -
  as.vector(X %*% column) / norm^3 +
  3 * as.vector(X %*% column) * column[j]^2 / norm^5
step <- 1e-4
displacement <- rep(0, p)
displacement[j] <- step
#' Evaluate a normalized projection for a single factor
#' @param loading Numeric p-vector.
#' @return Numeric n-vector with fixed epsilon regularization.
#' @examples
#' normalized_projection(F[, 1])
normalized_projection <- function(loading) {
  as.vector(X %*% loading) / sqrt(sum(loading^2) + epsilon)
}
check_close("Normalized-projection first derivative",
            first_derivative,
            (normalized_projection(column + displacement) -
               normalized_projection(column - displacement)) / (2 * step),
            1e-6)
check_close("Normalized-projection second derivative",
            second_derivative,
            (normalized_projection(column + displacement) -
               2 * normalized_projection(column) +
               normalized_projection(column - displacement)) / step^2,
            1e-6)

# Report ----

if (sys.nframe() == 0L) {
  print(check_results, row.names = FALSE)
  cat(sprintf("\n%d/%d derivation checks passed.\n",
              sum(check_results$passed), nrow(check_results)))
}
