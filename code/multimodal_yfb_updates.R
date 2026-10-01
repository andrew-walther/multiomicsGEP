# ============================================================
# Script: multimodal_yfb_updates.R
# Purpose: Exact modular updates for the isolated matched multimodal YFB model.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; code/multimodal_yfb_helpers.R
# ============================================================

if (!exists("multimodal_yfb_validate_moments")) {
  source("code/multimodal_yfb_helpers.R")
}

#' Update a point-exponential posterior from one Gaussian pseudo-observation
#'
#' @param A Positive quadratic precision.
#' @param B Linear quadratic coefficient.
#' @param prior List with `pi`, `rate`, and optional `point_mass` fields.
#' @return Posterior mean, second moment, slab probability, pseudo-observation,
#'   and pseudo-observation variance.
#' @examples
#' multimodal_yfb_point_exponential_posterior(2, 1, list(pi = 0.5, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_point_exponential_posterior <- function(A, B, prior) {
  if (!is.numeric(A) || length(A) != 1L || !is.numeric(B) || length(B) != 1L ||
      !is.finite(A) || !is.finite(B) || A < 0) {
    stop("A must be a finite nonnegative scalar and B must be finite.")
  }
  point_mass <- isTRUE(prior$point_mass)
  if (!point_mass && (!is.numeric(prior$pi) || !is.numeric(prior$rate) ||
                      length(prior$pi) != 1L || length(prior$rate) != 1L ||
                      !is.finite(prior$pi) || !is.finite(prior$rate) ||
                      prior$pi <= 0 || prior$pi > 1 || prior$rate <= 0)) {
    stop("A non-point-mass prior requires 0 < pi <= 1 and rate > 0.")
  }
  if (point_mass) {
    return(list(mean = 0, second = 0, slab_prob = 0, x = NA_real_, s2 = NA_real_))
  }
  if (A == 0) {
    # A collapsed factor can leave a subnormal residual B through floating-
    # point cancellation (for example 10^-238). Algebraically B is then zero:
    # no coordinate has information once its quadratic precision vanishes.
    # Treat only this underflow-scale remainder as the stated zero-information
    # boundary; a material B with A = 0 remains an invalid, unbounded update.
    if (abs(B) <= sqrt(.Machine$double.xmin)) B <- 0
    if (B != 0) stop("A = 0 with B != 0 is not a valid quadratic update.")
    return(list(mean = prior$pi / prior$rate,
                second = 2 * prior$pi / prior$rate^2,
                slab_prob = prior$pi, x = NA_real_, s2 = Inf))
  }

  x <- B / A
  s2 <- 1 / A
  s <- sqrt(s2)
  slab_location <- x - prior$rate * s2
  h <- slab_location / s
  log_spike <- if (prior$pi < 1) log1p(-prior$pi) + dnorm(x, 0, s, log = TRUE) else -Inf
  log_slab <- log(prior$pi) + log(prior$rate) - prior$rate * x +
    0.5 * prior$rate^2 * s2 + pnorm(h, log.p = TRUE)
  log_total <- max(log_spike, log_slab) +
    log(exp(log_spike - max(log_spike, log_slab)) + exp(log_slab - max(log_spike, log_slab)))
  slab_prob <- exp(log_slab - log_total)
  if (h < -8) {
    # For a = -h >> 0, direct truncated-Normal moments subtract nearly equal
    # O(a^2) terms. The inverse-Mills expansion below gives the positive
    # O(a^-1) mean and O(a^-2) second moment without cancellation.
    a <- -h
    inverse_a <- 1 / a
    delta <- inverse_a - 2 * inverse_a^3 + 10 * inverse_a^5 - 74 * inverse_a^7
    slab_mean <- s * delta
    slab_second <- s2 * (1 - a * delta)
  } else {
    mills <- exp(dnorm(h, log = TRUE) - pnorm(h, log.p = TRUE))
    slab_mean <- slab_location + s * mills
    slab_second <- slab_location^2 + s2 + slab_location * s * mills
  }
  mean <- slab_prob * slab_mean
  second <- max(slab_prob * slab_second, mean^2)
  list(mean = mean, second = second,
       slab_prob = slab_prob, x = x, s2 = s2)
}

#' Update one shared score column using both reconstruction blocks
#'
#' @param Y Named list of observed modality matrices.
#' @param R_minus_k Named list of factor-excluded residual matrices.
#' @param EF_k Named list of loading posterior means for one factor.
#' @param EF2_k Named list of loading posterior second moments for one factor.
#' @param Tau Named list of feature precision vectors.
#' @param prior Shared score point-exponential prior.
#' @return Posterior moments and pseudo-observation diagnostics for the score column.
#' @examples
#' multimodal_yfb_update_L_k(list(expression = matrix(1, 1, 1)),
#'   list(expression = matrix(1, 1, 1)), list(expression = 1),
#'   list(expression = 1), list(expression = 1), list(pi = 1, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_update_L_k <- function(Y, R_minus_k, EF_k, EF2_k, Tau, prior) {
  modalities <- names(Y)
  if (!identical(names(R_minus_k), modalities) || !identical(names(EF_k), modalities) ||
      !identical(names(EF2_k), modalities) || !identical(names(Tau), modalities)) {
    stop("All shared-score inputs must have identical modality names.")
  }
  n <- nrow(Y[[1]])
  A <- 0
  B <- numeric(n)
  for (modality in modalities) {
    if (!is.matrix(Y[[modality]]) || !is.matrix(R_minus_k[[modality]]) ||
        !identical(dim(Y[[modality]]), dim(R_minus_k[[modality]])) ||
        nrow(Y[[modality]]) != n || length(EF_k[[modality]]) != ncol(Y[[modality]]) ||
        length(Tau[[modality]]) != ncol(Y[[modality]])) {
      stop("Each modality must supply conformable Y, residual, loading, and Tau values.")
    }
    multimodal_yfb_validate_moments(EF_k[[modality]], EF2_k[[modality]],
                                    paste0("EF for ", modality))
    A <- A + sum(Tau[[modality]] * EF2_k[[modality]])
    B <- B + as.vector(R_minus_k[[modality]] %*% (Tau[[modality]] * EF_k[[modality]]))
  }
  posterior <- lapply(seq_len(n), function(i) {
    multimodal_yfb_point_exponential_posterior(A, B[i], prior)
  })
  list(A = A, B = B, x = vapply(posterior, `[[`, numeric(1), "x"),
       s2 = vapply(posterior, `[[`, numeric(1), "s2"),
       mean = vapply(posterior, `[[`, numeric(1), "mean"),
       second = vapply(posterior, `[[`, numeric(1), "second"),
       slab_prob = vapply(posterior, `[[`, numeric(1), "slab_prob"),
       prior = multimodal_yfb_update_loading_prior(
         vapply(posterior, `[[`, numeric(1), "slab_prob"),
         vapply(posterior, `[[`, numeric(1), "mean")
       ))
}

#' Maximize one modality-factor point-exponential loading prior
#'
#' @param slab_prob Posterior nonzero probabilities across the modality's features.
#' @param mean Posterior means across the same features.
#' @return List with empirical-Bayes `pi`, `rate`, and `point_mass` boundary flag.
#' @examples
#' multimodal_yfb_update_loading_prior(c(1, 0), c(1, 0))
#' @family multimodal_yfb_updates
multimodal_yfb_update_loading_prior <- function(slab_prob, mean) {
  if (!is.numeric(slab_prob) || !is.numeric(mean) || length(slab_prob) != length(mean) ||
      length(mean) == 0L || any(!is.finite(slab_prob)) || any(!is.finite(mean)) ||
      any(slab_prob < 0 | slab_prob > 1) || any(mean < 0)) {
    stop("slab_prob and mean must be finite, nonnegative vectors of equal positive length.")
  }
  total_slab <- sum(slab_prob)
  if (total_slab == 0) return(list(pi = 0, rate = NA_real_, point_mass = TRUE))
  total_mean <- sum(mean)
  # A finite slab probability can coexist with numerically zero posterior
  # means after a large factor-level rate. Its empirical-Bayes MLE is the
  # infinite-rate limit, represented here by the existing point-mass boundary.
  if (total_mean <= sqrt(.Machine$double.xmin)) {
    return(list(pi = 0, rate = NA_real_, point_mass = TRUE))
  }
  list(pi = total_slab / length(mean), rate = total_slab / total_mean, point_mass = FALSE)
}

#' Sequentially update one modality's loading column
#'
#' Implements the exact scalar Gauss-Seidel update: reconstruction and Cox
#' terms are added, and current projection moments are refreshed per feature.
#'
#' @param Y_m Subject-by-feature observed matrix for one modality.
#' @param Tau_m Feature precision vector.
#' @param EL_k,EL2_k Shared score posterior moments for one factor.
#' @param R_mk Factor-excluded reconstruction residual for this modality.
#' @param w,h_minus_k Frozen Cox working curvature and other-factor coefficient.
#' @param EZ_k,VZ_k Current joint projection mean and variance for this factor.
#' @param EF_k,EF2_k Current loading posterior moments for this modality-factor.
#' @param EBeta_k,EBeta2_k Coefficient posterior moments for this factor.
#' @param prior Fixed modality-factor point-exponential prior for this sweep.
#' @return Updated loading/projection moments, fitted prior, and per-feature diagnostics.
#' @examples
#' multimodal_yfb_update_F_mk(matrix(1, 1, 1), 1, 1, 1, matrix(1, 1, 1),
#'   1, 1, 1, 0, 1, 1, 1, list(pi = 1, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_update_F_mk <- function(Y_m, Tau_m, EL_k, EL2_k, R_mk, w,
                                        h_minus_k, EZ_k, VZ_k, EF_k, EF2_k,
                                        EBeta_k, EBeta2_k, prior) {
  if (!is.matrix(Y_m) || !is.matrix(R_mk) || !identical(dim(Y_m), dim(R_mk)) ||
      !is.numeric(Tau_m) || length(Tau_m) != ncol(Y_m) || length(EL_k) != nrow(Y_m) ||
      length(w) != nrow(Y_m) || length(h_minus_k) != nrow(Y_m) ||
      length(EZ_k) != nrow(Y_m) || length(VZ_k) != nrow(Y_m) ||
      length(EF_k) != ncol(Y_m) || length(EF2_k) != ncol(Y_m)) {
    stop("F-update inputs must be conformable with Y_m.")
  }
  multimodal_yfb_validate_moments(EL_k, EL2_k, "EL")
  multimodal_yfb_validate_moments(EF_k, EF2_k, "EF")
  multimodal_yfb_validate_moments(EBeta_k, EBeta2_k, "EBeta")
  if (any(!is.finite(w)) || any(w < 0) || any(!is.finite(h_minus_k)) ||
      any(!is.finite(EZ_k)) || any(!is.finite(VZ_k)) || any(VZ_k < -1e-10)) {
    stop("Cox working and projection moments must be finite with nonnegative variances.")
  }

  EF_new <- EF_k
  EF2_new <- EF2_k
  EZ_new <- EZ_k
  VZ_new <- pmax(VZ_k, 0)
  details <- vector("list", ncol(Y_m))
  sum_EL2 <- sum(EL2_k)
  for (j in seq_len(ncol(Y_m))) {
    y <- Y_m[, j]
    old_mean <- EF_new[j]
    old_variance <- pmax(EF2_new[j] - old_mean^2, 0)
    EZ_without_j <- EZ_new - y * old_mean
    A <- Tau_m[j] * sum_EL2 + EBeta2_k * sum(w * y^2)
    B <- Tau_m[j] * sum(EL_k * R_mk[, j]) +
      sum(y * (EBeta_k * h_minus_k - w * EBeta2_k * EZ_without_j))
    posterior <- multimodal_yfb_point_exponential_posterior(A, B, prior)
    new_variance <- pmax(posterior$second - posterior$mean^2, 0)
    EF_new[j] <- posterior$mean
    EF2_new[j] <- posterior$mean^2 + new_variance
    EZ_new <- EZ_without_j + y * posterior$mean
    VZ_new <- VZ_new + y^2 * (new_variance - old_variance)
    details[[j]] <- c(posterior, list(A = A, B = B, EZ_without_j = EZ_without_j))
  }
  list(EF = EF_new, EF2 = EF2_new, EZ = EZ_new, VZ = pmax(VZ_new, 0),
       prior = multimodal_yfb_update_loading_prior(
         vapply(details, `[[`, numeric(1), "slab_prob"), EF_new
       ), details = details)
}

#' Update one shared Normal-prior survival coefficient
#'
#' @param EZ_k,VZ_k Projection posterior mean and variance for one factor.
#' @param w,h_minus_k Frozen Cox working quantities.
#' @param prior_variance Prior variance used only for the zero-information boundary.
#' @return Coefficient posterior moments, fitted prior variance, and A/B diagnostics.
#' @examples
#' multimodal_yfb_update_beta_k(c(1, 2), c(0, 0), c(1, 1), c(1, 1))
#' @family multimodal_yfb_updates
multimodal_yfb_update_beta_k <- function(EZ_k, VZ_k, w, h_minus_k,
                                          prior_variance = 1) {
  if (!is.numeric(EZ_k) || !is.numeric(VZ_k) || !is.numeric(w) || !is.numeric(h_minus_k) ||
      length(EZ_k) == 0L || length(VZ_k) != length(EZ_k) || length(w) != length(EZ_k) ||
      length(h_minus_k) != length(EZ_k) || any(!is.finite(EZ_k)) ||
      any(!is.finite(VZ_k)) || any(VZ_k < 0) || any(!is.finite(w)) || any(w < 0)) {
    stop("Beta-update inputs must be finite conformable vectors with nonnegative variance and curvature.")
  }
  A <- sum(w * (EZ_k^2 + VZ_k))
  B <- sum(EZ_k * h_minus_k)
  if (A == 0) {
    if (B != 0) stop("A = 0 with B != 0 is not a valid beta update.")
    return(list(A = A, B = B, x = NA_real_, prior_variance = prior_variance,
                mean = 0, second = prior_variance, variance = prior_variance))
  }
  x <- B / A
  observed_variance <- 1 / A
  fitted_prior_variance <- max(0, x^2 - observed_variance)
  if (fitted_prior_variance == 0) {
    return(list(A = A, B = B, x = x, prior_variance = 0,
                mean = 0, second = 0, variance = 0))
  }
  variance <- 1 / (A + 1 / fitted_prior_variance)
  mean <- variance * B
  list(A = A, B = B, x = x, prior_variance = fitted_prior_variance,
       mean = mean, second = variance + mean^2, variance = variance)
}

#' Update modality-specific feature precisions by their exact MLEs
#'
#' @param Y Named list of observed modality matrices.
#' @param EL,EL2 Shared score posterior moments.
#' @param EF,EF2 Named lists of loading posterior moments.
#' @return List containing modality-specific Tau vectors and expected residuals.
#' @examples
#' multimodal_yfb_update_tau(list(expression = matrix(1, 1, 1)), matrix(1, 1, 1),
#'   matrix(1, 1, 1), list(expression = matrix(0, 1, 1)),
#'   list(expression = matrix(0, 1, 1)))
#' @family multimodal_yfb_updates
multimodal_yfb_update_tau <- function(Y, EL, EL2, EF, EF2) {
  if (!is.matrix(EL) || !is.matrix(EL2) || !identical(dim(EL), dim(EL2)) ||
      !identical(names(Y), names(EF)) || !identical(names(Y), names(EF2))) {
    stop("Tau-update inputs must include conformable shared and modality moments.")
  }
  multimodal_yfb_validate_moments(EL, EL2, "EL")
  Tau <- expected_residual2 <- vector("list", length(Y))
  names(Tau) <- names(expected_residual2) <- names(Y)
  for (modality in names(Y)) {
    if (!is.matrix(Y[[modality]]) || !is.matrix(EF[[modality]]) ||
        !is.matrix(EF2[[modality]]) || nrow(Y[[modality]]) != nrow(EL) ||
        ncol(Y[[modality]]) != nrow(EF[[modality]]) || ncol(EF[[modality]]) != ncol(EL)) {
      stop("Each modality must be conformable with the shared score moments.")
    }
    multimodal_yfb_validate_moments(EF[[modality]], EF2[[modality]],
                                    paste0("EF for ", modality))
    mean_signal <- EL %*% t(EF[[modality]])
    variance_signal <- EL2 %*% t(EF2[[modality]]) - EL^2 %*% t(EF[[modality]]^2)
    expected_residual2[[modality]] <- (Y[[modality]] - mean_signal)^2 +
      pmax(variance_signal, 0)
    expected_sse <- colSums(expected_residual2[[modality]])
    if (any(expected_sse <= 0)) {
      stop("A zero expected residual sum has no finite Tau MLE.")
    }
    Tau[[modality]] <- nrow(Y[[modality]]) / expected_sse
  }
  list(Tau = Tau, expected_residual2 = expected_residual2)
}
