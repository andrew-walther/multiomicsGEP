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
    return(list(mean = 0, second = 0, slab_prob = 0, x = NA_real_, s2 = NA_real_,
                log_ml = NA_real_))
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
                slab_prob = prior$pi, x = NA_real_, s2 = Inf, log_ml = NA_real_))
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
  slab <- multimodal_yfb_positive_truncnorm_moments(slab_location, s)
  mean <- slab_prob * slab$mean
  second <- max(slab_prob * slab$second, mean^2)
  # log_total = log p(x) under the prior: the marginal likelihood used for
  # the coordinate KL in the ELBO (see multimodal_yfb_coordinate_kl())
  list(mean = mean, second = second,
       slab_prob = slab_prob, x = x, s2 = s2, log_ml = log_total)
}

# Compiled loading sweep ----

.multimodal_yfb_cpp <- new.env()

#' Whether to use the compiled (C++) loading sweep
#'
#' Compiles `code/multimodal_yfb_sweep.cpp` with Rcpp once per R session. On
#' macOS, if the default build cannot find the C++ standard headers (a known
#' Command Line Tools issue), it retries with the SDK's libc++ include path,
#' set only for this build. If compilation still fails, it warns once and the
#' R loop is used (same results, slower). Force the R loop with
#' `options(multimodal_yfb.engine = "R")`.
#'
#' @return TRUE if the compiled sweep is available and enabled.
#' @examples
#' multimodal_yfb_use_cpp()
#' @family multimodal_yfb_updates
multimodal_yfb_use_cpp <- function() {
  if (identical(getOption("multimodal_yfb.engine"), "R")) return(FALSE)
  if (!is.null(.multimodal_yfb_cpp$ok)) return(.multimodal_yfb_cpp$ok)
  src <- "code/multimodal_yfb_sweep.cpp"
  build <- function() {
    utils::capture.output(Rcpp::sourceCpp(src, env = globalenv()))
    TRUE
  }
  ok <- requireNamespace("Rcpp", quietly = TRUE) && file.exists(src) &&
    tryCatch(suppressWarnings(build()), error = function(e) FALSE)
  if (!ok && Sys.info()[["sysname"]] == "Darwin") {
    sdk <- tryCatch(system2("xcrun", "--show-sdk-path", stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
    if (length(sdk) && nzchar(sdk)) {
      old <- Sys.getenv("PKG_CXXFLAGS")
      Sys.setenv(PKG_CXXFLAGS = paste0("-I", sdk, "/usr/include/c++/v1 -isysroot ", sdk))
      ok <- tryCatch(suppressWarnings(build()), error = function(e) FALSE)
      Sys.setenv(PKG_CXXFLAGS = old)
    }
  }
  if (!ok) warning("Compiled multimodal YFB sweep unavailable; using the (slower) R loop.")
  .multimodal_yfb_cpp$ok <- ok
  ok
}

#' Moments of a Normal truncated to the positive half-line
#'
#' For theta ~ N(m, s^2) restricted to theta >= 0, with h = m / s and
#' inverse Mills ratio r(h) = phi(h) / Phi(h):
#'   E[theta] = m + s r(h),  E[theta^2] = m^2 + s^2 + m s r(h).
#' For h < -8 the direct formulas subtract nearly equal O(h^2) terms, so an
#' inverse-Mills asymptotic expansion is used instead.
#'
#' @param location Untruncated mean m.
#' @param s Untruncated standard deviation (> 0).
#' @return List with `mean` and `second` moments.
#' @examples
#' multimodal_yfb_positive_truncnorm_moments(0, 1)
#' @family multimodal_yfb_updates
multimodal_yfb_positive_truncnorm_moments <- function(location, s) {
  h <- location / s
  if (h < -8) {
    # For a = -h >> 0, direct truncated-Normal moments subtract nearly equal
    # O(a^2) terms. The inverse-Mills expansion below gives the positive
    # O(a^-1) mean and O(a^-2) second moment without cancellation.
    a <- -h
    inverse_a <- 1 / a
    delta <- inverse_a - 2 * inverse_a^3 + 10 * inverse_a^5 - 74 * inverse_a^7
    list(mean = s * delta, second = s^2 * (1 - a * delta))
  } else {
    mills <- exp(dnorm(h, log = TRUE) - pnorm(h, log.p = TRUE))
    list(mean = location + s * mills,
         second = location^2 + s^2 + location * s * mills)
  }
}

#' Update a signed (point-Laplace or Normal) posterior from one pseudo-observation
#'
#' Pseudo-observation x = B / A with variance s^2 = 1 / A, x | theta ~ N(theta, s^2).
#'
#' point_laplace: theta ~ (1 - pi) delta_0 + pi (lambda / 2) exp(-lambda |theta|).
#'   The slab posterior is a mixture of a positive truncated N(x - lambda s^2, s^2)
#'   and a negative truncated N(x + lambda s^2, s^2), with log weights
#'   log(pi lambda / 2) -/+ lambda x + lambda^2 s^2 / 2 + log Phi((+/-x - lambda s^2) / s).
#'   The spike has log weight log(1 - pi) + log N(x; 0, s^2).
#' normal: theta ~ N(0, sigma^2). Posterior variance v = 1 / (A + 1 / sigma^2),
#'   mean v B; sigma^2 = 0 is the point mass at zero.
#'
#' @param A Nonnegative quadratic precision.
#' @param B Linear quadratic coefficient.
#' @param prior List with `family` ("point_laplace" or "normal"), `point_mass`,
#'   and `pi`, `rate` (point_laplace) or `variance` (normal).
#' @return Posterior `mean`, `second`, nonzero probability `slab_prob`,
#'   pseudo-observation `x`, and its variance `s2`.
#' @examples
#' multimodal_yfb_signed_posterior(2, -1, list(family = "point_laplace",
#'   pi = 0.5, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_signed_posterior <- function(A, B, prior) {
  if (!is.numeric(A) || length(A) != 1L || !is.numeric(B) || length(B) != 1L ||
      !is.finite(A) || !is.finite(B) || A < 0) {
    stop("A must be a finite nonnegative scalar and B must be finite.")
  }
  family <- prior$family
  if (!family %in% c("point_laplace", "normal")) {
    stop("Signed prior family must be point_laplace or normal.")
  }
  if (isTRUE(prior$point_mass)) {
    return(list(mean = 0, second = 0, slab_prob = 0,
                x = if (A > 0) B / A else NA_real_, s2 = if (A > 0) 1 / A else Inf,
                log_ml = NA_real_))
  }
  if (A == 0) {
    if (abs(B) <= sqrt(.Machine$double.xmin)) B <- 0
    if (B != 0) stop("A = 0 with B != 0 is not a valid quadratic update.")
    # No information: the posterior is the prior (mean zero by symmetry)
    second <- if (family == "normal") prior$variance else 2 * prior$pi / prior$rate^2
    return(list(mean = 0, second = second, slab_prob = if (family == "normal") 1 else prior$pi,
                x = NA_real_, s2 = Inf, log_ml = NA_real_))
  }
  x <- B / A
  s2 <- 1 / A
  s <- sqrt(s2)
  if (family == "normal") {
    variance <- 1 / (A + 1 / prior$variance)
    mean <- variance * B
    return(list(mean = mean, second = variance + mean^2, slab_prob = 1, x = x, s2 = s2,
                log_ml = dnorm(x, 0, sqrt(s2 + prior$variance), log = TRUE)))
  }
  lambda <- prior$rate
  log_slab_constant <- log(prior$pi) + log(lambda / 2) + 0.5 * lambda^2 * s2
  log_pos <- log_slab_constant - lambda * x + pnorm((x - lambda * s2) / s, log.p = TRUE)
  log_neg <- log_slab_constant + lambda * x + pnorm((-x - lambda * s2) / s, log.p = TRUE)
  log_spike <- if (prior$pi < 1) log1p(-prior$pi) + dnorm(x, 0, s, log = TRUE) else -Inf
  log_w <- c(spike = log_spike, pos = log_pos, neg = log_neg)
  log_ml <- max(log_w) + log(sum(exp(log_w - max(log_w))))
  weights <- exp(log_w - log_ml)
  pos <- multimodal_yfb_positive_truncnorm_moments(x - lambda * s2, s)
  # theta < 0 part: -theta is a positive truncated N(-x - lambda s^2, s^2)
  neg <- multimodal_yfb_positive_truncnorm_moments(-x - lambda * s2, s)
  mean <- weights[["pos"]] * pos$mean - weights[["neg"]] * neg$mean
  second <- max(weights[["pos"]] * pos$second + weights[["neg"]] * neg$second, mean^2)
  list(mean = mean, second = second, slab_prob = 1 - weights[["spike"]], x = x, s2 = s2,
       log_ml = log_ml)
}

#' KL divergence of one coordinate's posterior from its prior
#'
#' For a posterior q formed from pseudo-observation x ~ N(theta, s^2) and
#' prior g, Bayes' rule gives the EBNM identity used by flashier:
#'   KL(q || g) = E_q[log N(x; theta, s^2)] - log p(x),
#'   E_q[log N(x; theta, s^2)] = -log(2 pi s^2)/2 - (x^2 - 2 x E[theta] + E[theta^2]) / (2 s^2),
#' with p(x) the marginal likelihood of x under g. Coordinates with no
#' information (s^2 = Inf, so q = g) or a point-mass prior (q = g = delta_0)
#' have KL = 0.
#'
#' @param x,s2 Pseudo-observations and their variances.
#' @param mean,second Posterior first and second moments.
#' @param log_ml Log marginal likelihoods log p(x).
#' @return Sum of the coordinate KL divergences (nonnegative up to roundoff).
#' @examples
#' multimodal_yfb_coordinate_kl(1, 1, 0.5, 0.75, dnorm(1, 0, sqrt(2), log = TRUE))
#' @family multimodal_yfb_updates
multimodal_yfb_coordinate_kl <- function(x, s2, mean, second, log_ml) {
  ok <- is.finite(x) & is.finite(s2) & s2 > 0 & is.finite(log_ml)
  if (!any(ok)) return(0)
  expected_loglik <- -0.5 * log(2 * pi * s2[ok]) -
    (x[ok]^2 - 2 * x[ok] * mean[ok] + second[ok]) / (2 * s2[ok])
  sum(expected_loglik - log_ml[ok])
}

#' Fit a point-exponential prior by marginal likelihood (ebnm)
#'
#' Unlike the moment-based EM update, the marginal-likelihood fit can return
#' the point mass at zero (pi = 0) when the pseudo-observations carry no
#' signal, which removes the whole factor column -- the mechanism flashier
#' uses to prune factors.
#'
#' @param x,s2 Pseudo-observations and their variances.
#' @param init Optional current prior (warm start for the optimizer).
#' @return Prior list (`pi`, `rate`, `point_mass`, `family`).
#' @examples
#' multimodal_yfb_fit_exponential_prior_ebnm(c(0.1, 2, 3), c(1, 1, 1))
#' @family multimodal_yfb_updates
multimodal_yfb_fit_exponential_prior_ebnm <- function(x, s2, init = NULL) {
  ok <- is.finite(x) & is.finite(s2) & s2 > 0
  dead <- list(family = "point_exponential", pi = 0, rate = NA_real_, point_mass = TRUE)
  if (!any(ok)) return(dead)
  g_init <- if (!is.null(init) && !isTRUE(init$point_mass) && is.finite(init$rate %||% NA))
    ebnm::gammamix(pi = c(1 - init$pi, init$pi), shape = c(1, 1), scale = c(0, 1 / init$rate),
                   shift = c(0, 0))
  else NULL
  g <- ebnm::ebnm_point_exponential(x[ok], sqrt(s2[ok]), mode = 0, g_init = g_init)$fitted_g
  # gammamix: component 1 is the point mass (scale 0); component 2 is
  # Exp(rate = 1 / scale)
  if (g$pi[2] <= sqrt(.Machine$double.eps) || g$scale[2] <= 0) return(dead)
  list(family = "point_exponential", pi = g$pi[2], rate = 1 / g$scale[2], point_mass = FALSE)
}

#' Fit a signed modality-factor loading prior by empirical Bayes (ebnm)
#'
#' Maximizes the marginal likelihood of the sweep's pseudo-observations
#' (x_j, s_j), as flashier does for each factor. The mode is fixed at zero.
#' Pseudo-observations without information (s_j = Inf) are excluded.
#'
#' @param x Pseudo-observations, one per feature.
#' @param s2 Their variances.
#' @param family "point_laplace" or "normal".
#' @param init Optional current prior (warm start for the optimizer).
#' @return Prior list for [multimodal_yfb_signed_posterior()].
#' @examples
#' multimodal_yfb_fit_signed_prior(c(-2, 0.1, 3), c(1, 1, 1), "point_laplace")
#' @family multimodal_yfb_updates
multimodal_yfb_fit_signed_prior <- function(x, s2, family, init = NULL) {
  ok <- is.finite(x) & is.finite(s2) & s2 > 0
  if (!any(ok)) {
    return(list(family = family, pi = 0, rate = NA_real_, variance = 0, point_mass = TRUE))
  }
  if (family == "normal") {
    g_init <- if (!is.null(init) && !isTRUE(init$point_mass) && is.finite(init$variance %||% NA))
      ashr::normalmix(1, 0, sqrt(init$variance)) else NULL
    g <- ebnm::ebnm_normal(x[ok], sqrt(s2[ok]), mode = 0, g_init = g_init)$fitted_g
    variance <- g$sd^2
    return(list(family = family, variance = variance, point_mass = variance <= 0))
  }
  # Warm start from the current prior: the same marginal-likelihood problem,
  # started closer to its optimum
  g_init <- if (!is.null(init) && !isTRUE(init$point_mass) && is.finite(init$rate %||% NA))
    ebnm::laplacemix(pi = c(1 - init$pi, init$pi), mean = c(0, 0), scale = c(0, 1 / init$rate))
  else NULL
  g <- ebnm::ebnm_point_laplace(x[ok], sqrt(s2[ok]), mode = 0, g_init = g_init)$fitted_g
  # laplacemix: component 1 is the point mass (scale 0); component 2 has
  # density exp(-|theta| / scale) / (2 scale), so rate lambda = 1 / scale
  slab_pi <- g$pi[2]
  if (slab_pi <= sqrt(.Machine$double.eps) || g$scale[2] <= 0) {
    return(list(family = family, pi = 0, rate = NA_real_, point_mass = TRUE))
  }
  list(family = family, pi = slab_pi, rate = 1 / g$scale[2], point_mass = FALSE)
}

#' Update one shared score column using both reconstruction blocks
#'
#' @param Y Named list of observed modality matrices.
#' @param R_minus_k Named list of factor-excluded residual matrices.
#' @param EF_k Named list of loading posterior means for one factor.
#' @param EF2_k Named list of loading posterior second moments for one factor.
#' @param Tau Named list of feature precision vectors.
#' @param prior Shared score point-exponential prior.
#' @param prior_update "em" (moment update) or "ebnm" (marginal-likelihood fit,
#'   which can shrink the whole column to the point mass).
#' @return Posterior moments and pseudo-observation diagnostics for the score column.
#' @examples
#' multimodal_yfb_update_L_k(list(expression = matrix(1, 1, 1)),
#'   list(expression = matrix(1, 1, 1)), list(expression = 1),
#'   list(expression = 1), list(expression = 1), list(pi = 1, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_update_L_k <- function(Y, R_minus_k, EF_k, EF2_k, Tau, prior,
                                      prior_update = "em") {
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
  if (multimodal_yfb_use_cpp()) {
    # All n coordinates share A, so the posterior is one vectorized call
    posterior <- mmyfb_posterior_vec_cpp(A, B, 0L, isTRUE(prior$point_mass),
                                         prior$pi %||% NA_real_, prior$rate %||% NA_real_,
                                         NA_real_)
    get <- function(field) posterior[[field]]
  } else {
    posterior <- lapply(seq_len(n), function(i) {
      multimodal_yfb_point_exponential_posterior(A, B[i], prior)
    })
    get <- function(field) vapply(posterior, `[[`, numeric(1), field)
  }
  x <- get("x"); s2 <- get("s2"); mean <- get("mean"); second <- get("second")
  new_prior <- if (prior_update == "ebnm") {
    multimodal_yfb_fit_exponential_prior_ebnm(x, s2, prior)
  } else {
    multimodal_yfb_update_loading_prior(get("slab_prob"), mean)
  }
  list(A = A, B = B, x = x, s2 = s2, mean = mean, second = second,
       slab_prob = get("slab_prob"), prior = new_prior,
       kl = multimodal_yfb_coordinate_kl(x, s2, mean, second, get("log_ml")))
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
#' @param prior_update For point-exponential loadings: "em" (moment update)
#'   or "ebnm" (marginal-likelihood fit). Signed priors always use ebnm.
#' @param prior Fixed modality-factor prior for this sweep: point-exponential
#'   (`family` NULL or "point_exponential"), "point_laplace" or "normal".
#' @return Updated loading/projection moments, fitted prior, and per-feature diagnostics.
#' @examples
#' multimodal_yfb_update_F_mk(matrix(1, 1, 1), 1, 1, 1, matrix(1, 1, 1),
#'   1, 1, 1, 0, 1, 1, 1, list(pi = 1, rate = 1))
#' @family multimodal_yfb_updates
multimodal_yfb_update_F_mk <- function(Y_m, Tau_m, EL_k, EL2_k, R_mk, w,
                                        h_minus_k, EZ_k, VZ_k, EF_k, EF2_k,
                                        EBeta_k, EBeta2_k, prior,
                                        prior_update = "em") {
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
  # prior$family is NULL or "point_exponential" for the original nonnegative
  # loadings; "point_laplace" and "normal" give signed loadings
  signed <- !is.null(prior$family) && prior$family != "point_exponential"
  # Terms that do not change during the feature sweep, computed once
  # (algebraically identical to computing them inside the loop):
  #   sum_i w_i y_ij^2,  sum_i EL_ik R_ij,  sum_i y_ij h_i
  wY <- w * Y_m
  sum_wy2 <- unname(colSums(wY * Y_m))
  sum_LR <- as.vector(crossprod(R_mk, EL_k))
  sum_yh <- as.vector(crossprod(Y_m, h_minus_k))
  family_code <- match(if (is.null(prior$family)) "point_exponential" else prior$family,
                       c("point_exponential", "point_laplace", "normal")) - 1L
  if (multimodal_yfb_use_cpp()) {
    # Same loop in C++ (code/multimodal_yfb_sweep.cpp); see the R branch for the math
    res <- mmyfb_F_sweep_cpp(
      Y_m, Tau_m, sum_EL2, sum_wy2, sum_LR, sum_yh, w, EZ_new, VZ_new, EF_new, EF2_new,
      EBeta_k, EBeta2_k, family_code, isTRUE(prior$point_mass),
      prior$pi %||% NA_real_, prior$rate %||% NA_real_, prior$variance %||% NA_real_)
    EF_new <- res$EF; EF2_new <- res$EF2; EZ_new <- res$EZ; VZ_new <- res$VZ
    details <- res[c("mean", "second", "slab_prob", "x", "s2", "log_ml", "A", "B")]
  } else {
    for (j in seq_len(ncol(Y_m))) {
      y <- Y_m[, j]
      old_mean <- EF_new[j]
      old_variance <- max(EF2_new[j] - old_mean^2, 0)
      # EZ_without_j = EZ - y_j * old_mean, so
      # sum_i w_i y_ij EZ_without_j,i = sum_i w_i y_ij EZ_i - old_mean * sum_wy2[j]
      sum_wyEZ_without_j <- sum(wY[, j] * EZ_new) - old_mean * sum_wy2[j]
      A <- Tau_m[j] * sum_EL2 + EBeta2_k * sum_wy2[j]
      B <- Tau_m[j] * sum_LR[j] + EBeta_k * sum_yh[j] - EBeta2_k * sum_wyEZ_without_j
      posterior <- if (signed) {
        multimodal_yfb_signed_posterior(A, B, prior)
      } else {
        multimodal_yfb_point_exponential_posterior(A, B, prior)
      }
      new_variance <- pmax(posterior$second - posterior$mean^2, 0)
      EF_new[j] <- posterior$mean
      EF2_new[j] <- posterior$mean^2 + new_variance
      EZ_new <- EZ_new + y * (posterior$mean - old_mean)
      VZ_new <- VZ_new + y^2 * (new_variance - old_variance)
      details[[j]] <- c(posterior, list(A = A, B = B))
    }
    fields <- c("mean", "second", "slab_prob", "x", "s2", "log_ml", "A", "B")
    details <- stats::setNames(lapply(fields, function(f) {
      vapply(details, function(d) as.numeric(d[[f]]), numeric(1))
    }), fields)
  }
  get <- function(field) details[[field]]
  new_prior <- if (signed) {
    multimodal_yfb_fit_signed_prior(get("x"), get("s2"), prior$family, prior)
  } else if (prior_update == "ebnm") {
    multimodal_yfb_fit_exponential_prior_ebnm(get("x"), get("s2"), prior)
  } else {
    c(multimodal_yfb_update_loading_prior(get("slab_prob"), EF_new),
      if (!is.null(prior$family)) list(family = prior$family))
  }
  list(EF = EF_new, EF2 = EF2_new, EZ = EZ_new, VZ = pmax(VZ_new, 0),
       prior = new_prior, details = details,
       kl = multimodal_yfb_coordinate_kl(get("x"), get("s2"), get("mean"),
                                         get("second"), get("log_ml")))
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
                mean = 0, second = prior_variance, variance = prior_variance, kl = 0))
  }
  x <- B / A
  observed_variance <- 1 / A
  fitted_prior_variance <- max(0, x^2 - observed_variance)
  if (fitted_prior_variance == 0) {
    return(list(A = A, B = B, x = x, prior_variance = 0,
                mean = 0, second = 0, variance = 0, kl = 0))
  }
  variance <- 1 / (A + 1 / fitted_prior_variance)
  mean <- variance * B
  # KL(N(mean, variance) || N(0, sigma^2)) = [(variance + mean^2)/sigma^2 - 1 - log(variance/sigma^2)] / 2
  kl <- 0.5 * ((variance + mean^2) / fitted_prior_variance - 1 -
                 log(variance / fitted_prior_variance))
  list(A = A, B = B, x = x, prior_variance = fitted_prior_variance,
       mean = mean, second = variance + mean^2, variance = variance, kl = kl)
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
