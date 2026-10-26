# ============================================================
# Script: loading_prior_fixed_q.R
# Purpose: Maximize a point-Laplace loading prior with the posterior held fixed.
# Created: 2026-10-16
# Dependencies: base R and multimodal_yfb_updates.R posterior helpers
# ============================================================

#' Validate sufficient statistics for a mixed point-Laplace posterior
#' @param prob Posterior nonzero probabilities, one per feature.
#' @param mean_abs Unconditional E[|F_j|], including spike mass.
#' @return Invisibly TRUE, or an error for inconsistent sufficient statistics.
#' @examples
#' loading_prior_validate_stats(c(.2, 1), c(.1, .5))
#' @family loading_prior_fixed_q
loading_prior_validate_stats <- function(prob, mean_abs) {
  if (!is.numeric(prob) || !is.numeric(mean_abs) || !length(prob) ||
      length(prob) != length(mean_abs) || any(!is.finite(prob)) ||
      any(prob < 0 | prob > 1) || any(!is.finite(mean_abs)) ||
      any(mean_abs < 0) || any(prob == 0 & mean_abs != 0) ||
      any(prob > 0 & mean_abs == 0))
    stop("Invalid point-Laplace posterior statistics: positive mass needs positive magnitude.")
  invisible(TRUE)
}

#' Reconstruct point-Laplace posterior mass and absolute moments
#'
#' Uses each coordinate's actual incoming prior and frozen A/B. Positive and
#' negative truncated-Normal pieces are added for E[|F|], rather than taking
#' |E[F]|. The latter is zero for a symmetric posterior even with slab mass.
#' No information (A=B=0) returns prior moments; a delta-zero prior stays zero.
#' @param A,B Matching finite precision (A >= 0) and linear-signal vectors.
#' @param prior Incoming point-Laplace prior list.
#' @return List of posterior `prob` and unconditional `mean_abs` vectors.
#' @examples
#' loading_prior_laplace_stats(4, 2,
#'   list(family="point_laplace", pi=.5, rate=1, point_mass=FALSE))
#' @family loading_prior_fixed_q
#' @seealso [multimodal_yfb_signed_posterior()]
loading_prior_laplace_stats <- function(A, B, prior) {
  if (!identical(prior$family, "point_laplace") || !is.numeric(A) ||
      !is.numeric(B) || !length(A) || length(A) != length(B) ||
      any(!is.finite(A)) || any(A < 0) || any(!is.finite(B)) ||
      any(A == 0 & B != 0))
    stop("Finite matching A>=0, B and a point-Laplace prior required; A=0 needs B=0.")
  pi <- prior$pi
  if (length(pi) != 1L || !is.finite(pi) || pi < 0 || pi > 1 ||
      (isTRUE(prior$point_mass) && pi != 0)) stop("Invalid prior probability.")
  if (pi == 0) return(list(prob = rep(0, length(A)), mean_abs = rep(0, length(A))))
  rate <- prior$rate
  if (length(rate) != 1L || !is.finite(rate) || rate <= 0)
    stop("Positive finite Laplace rate required.")
  prob <- mean_abs <- numeric(length(A))
  for (j in seq_along(A)) {
    if (A[j] == 0) {
      prob[j] <- pi; mean_abs[j] <- pi/rate
    } else {
      x <- B[j]/A[j]; s2 <- 1/A[j]; s <- sqrt(s2)
      common <- log(pi) + log(rate/2) + rate^2*s2/2
      lw <- c(if (pi < 1) log1p(-pi) + dnorm(x, 0, s, log=TRUE) else -Inf,
              common - rate*x + pnorm((x-rate*s2)/s, log.p=TRUE),
              common + rate*x + pnorm((-x-rate*s2)/s, log.p=TRUE))
      normalizer <- max(lw) + log(sum(exp(lw-max(lw))))
      weights <- exp(lw-normalizer)
      # Sum the slab pieces directly: 1-spike loses tiny positive probabilities.
      prob[j] <- if (pi == 1) 1 else min(1, weights[2]+weights[3])
      pos <- multimodal_yfb_positive_truncnorm_moments(x-rate*s2, s)$mean
      neg <- multimodal_yfb_positive_truncnorm_moments(-x-rate*s2, s)$mean
      mean_abs[j] <- weights[2]*pos + weights[3]*neg
    }
  }
  loading_prior_validate_stats(prob, mean_abs)
  list(prob=prob, mean_abs=mean_abs)
}

#' Maximize E_q log g for one point-Laplace modality-factor prior
#'
#' For p features, S=sum Pr_q(F_j!=0), T=sum E_q|F_j|, the exact fixed-q
#' maximizer is pi=S/p and rate=S/T. Zero S implies delta-zero and an
#' unidentified slab rate; positive S with zero T has no finite rate maximizer.
#' This is one EM step, not marginal maximization of changing pseudo-data.
#' @param prob,mean_abs Posterior sufficient-statistic vectors.
#' @return Point-Laplace prior list; no floors, rate caps or posterior changes.
#' @examples
#' loading_prior_laplace_em(c(.2, .8), c(.1, .4))
#' @family loading_prior_fixed_q
loading_prior_laplace_em <- function(prob, mean_abs) {
  loading_prior_validate_stats(prob, mean_abs)
  S <- sum(prob); T <- sum(mean_abs)
  if (S == 0) return(list(family="point_laplace", pi=0, rate=NA_real_, point_mass=TRUE))
  rate <- S/T
  if (!is.finite(rate) || rate <= 0) stop("No representable finite positive EM rate.")
  list(family="point_laplace", pi=S/length(prob), rate=rate, point_mass=FALSE)
}

#' Evaluate a mixed-measure expected point-Laplace log prior
#' @param prob,mean_abs Posterior sufficient-statistic vectors.
#' @param prior Point-Laplace prior being evaluated, not necessarily incoming.
#' @return Sum E_q log g, with -Inf for incompatible spike/slab support.
#' @examples
#' loading_prior_laplace_log(.5, .2,
#'   list(family="point_laplace", pi=.4, rate=1, point_mass=FALSE))
#' @family loading_prior_fixed_q
loading_prior_laplace_log <- function(prob, mean_abs, prior) {
  loading_prior_validate_stats(prob, mean_abs)
  if (!identical(prior$family, "point_laplace") || length(prior$pi) != 1L ||
      !is.finite(prior$pi) || prior$pi < 0 || prior$pi > 1)
    stop("A valid point-Laplace prior is required.")
  pi <- prior$pi
  if (isTRUE(prior$point_mass) && pi != 0) stop("Point mass requires pi=0.")
  if (pi == 0) return(if (any(prob > 0)) -Inf else 0)
  if (pi == 1 && any(prob < 1)) return(-Inf)
  if (length(prior$rate) != 1L || !is.finite(prior$rate) || prior$rate <= 0)
    stop("Positive finite Laplace rate required.")
  value <- sum(prob)*(log(pi)+log(prior$rate/2)) - prior$rate*sum(mean_abs)
  if (pi < 1) value <- value + sum(1-prob)*log1p(-pi)
  value
}

#' Update the loading prior and KL while leaving the sweep posterior fixed
#'
#' KL(q||g_new)=KL(q||g_old)-E_q log(g_new/g_old). The old normalizer
#' belongs to each coordinate's actual sweep posterior, not refreshed pseudo-data.
#' @param details Sweep diagnostics A, B, x, s2, mean, second and log_ml.
#' @param prior Incoming point-Laplace prior.
#' @return New `prior`, corrected positive `kl`, and fixed-q `prior_gain`.
#' @examples
#' p <- multimodal_yfb_signed_posterior(4, 2,
#'   list(family="point_laplace", pi=.5, rate=1, point_mass=FALSE))
#' loading_prior_fixed_q(c(p, list(A=4, B=2)),
#'   list(family="point_laplace", pi=.5, rate=1, point_mass=FALSE))
#' @family loading_prior_fixed_q
#' @seealso [multimodal_yfb_coordinate_kl()]
loading_prior_fixed_q <- function(details, prior) {
  fields <- c("A", "B", "x", "s2", "mean", "second", "log_ml")
  if (!is.list(details) || !all(fields %in% names(details)) ||
      any(vapply(details[fields], length, integer(1)) != length(details$A)))
    stop("Complete matching posterior diagnostics are required.")
  if (any(!is.finite(details$mean)) || any(!is.finite(details$second)))
    stop("Finite native posterior moments are required.")
  q <- loading_prior_laplace_stats(details$A, details$B, prior)
  new <- loading_prior_laplace_em(q$prob, q$mean_abs)
  gain <- loading_prior_laplace_log(q$prob, q$mean_abs, new) -
    loading_prior_laplace_log(q$prob, q$mean_abs, prior)
  informative <- details$A > 0 & prior$pi > 0
  if (any(informative & (!is.finite(details$x) | !is.finite(details$s2) |
                         details$s2 <= 0)))
    stop("Informative posterior pseudo-data must be finite with positive variance.")
  if (any(informative & !is.finite(details$log_ml)))
    stop("An informative posterior requires its actual finite normalizer.")
  old_kl <- multimodal_yfb_coordinate_kl(details$x, details$s2, details$mean,
                                        details$second, details$log_ml)
  kl <- old_kl-gain
  if (!is.finite(gain) || !is.finite(kl)) stop("Nonfinite fixed-q prior gain or KL.")
  list(prior=new, kl=kl, prior_gain=gain)
}
