# Empirical Bayes normal means with exact KL ----

#' Fit an empirical Bayes normal-means problem and return posterior moments and KL
#'
#' Solves x_j ~ N(theta_j, s_j^2), theta_j ~ g, with g in `family` fitted by
#' marginal likelihood (`ebnm::ebnm()`), and returns the posterior moments and
#' KL(q || g) from the EBMF identity
#'   KL = sum_j E_q log N(x_j; theta_j, s_j^2) - log p_g(x).
#' Entries with zero precision (s = Inf) carry no information; they are given
#' a very large s so that their posterior equals the fitted prior.
#'
#' @param x pseudo-observations
#' @param s their standard errors (> 0; `Inf` allowed)
#' @param family one of "point_exponential", "point_laplace", "normal", "point_normal"
#' @param g_init optional prior to start the optimization from
#' @param fix_g if TRUE, keep `g_init` fixed (posterior only)
#' @return list(`mean`, `mean2`, `kl` (summed), `g`, `loglik`)
#' @keywords internal
#' @noRd
.ebnm_fit <- function(x, s, family, g_init = NULL, fix_g = FALSE) {
  s[is.na(s)] <- Inf
  s[s > 1e50] <- Inf                       # precision below 1e-100: no information
  if (!any(is.finite(s))) {
    # A column with no information (e.g. a program switched off in L): the
    # optimum over (q, g) is q = g = point mass at 0, with KL = 0.
    z <- numeric(length(x))
    return(list(mean = z, mean2 = z, kl = 0, g = NULL, loglik = 0))
  }
  s_big <- 1e6 * max(c(1, abs(x[is.finite(s) & is.finite(x)]), s[is.finite(s)]))
  s <- ifelse(is.finite(s), s, s_big)
  x[!is.finite(x)] <- 0
  # The normal-means problem is scale equivariant: solve it with standard
  # errors near 1 (ebnm loses accuracy when s is very large or very small)
  # and map the posterior and prior back. The KL is scale invariant.
  cs <- stats::median(s)
  x <- x / cs; s <- s / cs
  if (!is.null(g_init)) g_init <- .rescale_g(g_init, 1 / cs)
  out <- c("posterior_mean", "posterior_second_moment", "fitted_g", "log_likelihood")
  fit <- tryCatch(
    if (fix_g) ebnm::ebnm(x, s, prior_family = family, g_init = g_init, fix_g = TRUE, output = out)
    else ebnm::ebnm(x, s, prior_family = family, mode = 0, g_init = g_init, output = out),
    error = function(e) {
      if (is.null(g_init) || fix_g) stop(e)
      ebnm::ebnm(x, s, prior_family = family, mode = 0, output = out)
    })
  m1 <- fit$posterior$mean
  m2 <- pmax(fit$posterior$second_moment, m1^2)
  ell <- sum(-0.5 * log(2 * pi * s^2) - ((x - m1)^2 + (m2 - m1^2)) / (2 * s^2))
  ll <- as.numeric(fit$log_likelihood)
  kl <- ell - ll
  # Impossible results mean ebnm lost accuracy: this happens when the fitted
  # prior's scale is many orders of magnitude below the standard errors, i.e.
  # the column carries no usable information and the program is switching
  # off. The marginal log-likelihood cannot exceed sum(-log(2 pi s^2)/2), and
  # the KL cannot be negative beyond round-off. Return the no-information
  # solution (q = g = point mass at 0); the caller's objective check decides
  # whether the step is accepted.
  ll_max <- sum(-0.5 * log(2 * pi * s^2))
  if (!is.finite(kl) || !is.finite(ll) || ll > ll_max + 1e-6 * abs(ll_max) || kl < -1e-6 * length(x) - 1e-8) {
    z <- numeric(length(x))
    return(list(mean = z, mean2 = z, kl = 0, g = NULL, loglik = 0, degenerate = TRUE))
  }
  kl <- max(kl, 0)                          # round-off below 1e-6 per entry
  list(mean = m1 * cs, mean2 = m2 * cs^2, kl = kl, g = .rescale_g(fit$fitted_g, cs),
       loglik = ll - length(x) * log(cs))
}

#' Rescale a fitted ebnm prior: theta -> c * theta
#' @keywords internal
#' @noRd
.rescale_g <- function(g, c) {
  for (f in intersect(c("scale", "sd", "mean", "shift"), names(g))) g[[f]] <- g[[f]] * c
  g
}

