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
.ebnm_fit <- function(x, s, family, g_init = NULL, fix_g = FALSE) {
  s_big <- 1e6 * max(c(1, abs(x[is.finite(s) & is.finite(x)]), s[is.finite(s)]))
  s <- ifelse(is.finite(s), s, s_big)
  x[!is.finite(x)] <- 0
  out <- c("posterior_mean", "posterior_second_moment", "fitted_g", "log_likelihood")
  fit <- tryCatch(
    if (fix_g) ebnm::ebnm(x, s, prior_family = family, g_init = g_init, fix_g = TRUE, output = out)
    else ebnm::ebnm(x, s, prior_family = family, mode = 0, g_init = g_init, output = out),
    error = function(e) {
      if (is.null(g_init)) stop(e)
      ebnm::ebnm(x, s, prior_family = family, mode = 0, output = out)
    })
  m1 <- fit$posterior$mean
  m2 <- pmax(fit$posterior$second_moment, m1^2)
  ell <- sum(-0.5 * log(2 * pi * s^2) - ((x - m1)^2 + (m2 - m1^2)) / (2 * s^2))
  ll <- as.numeric(fit$log_likelihood)
  list(mean = m1, mean2 = m2, kl = ell - ll, g = fit$fitted_g, loglik = ll)
}
