# ============================================================
# Script: multimodal_yfb_helpers.R
# Purpose: Mathematical helpers for the isolated matched multimodal YFB model.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================

#' Validate posterior mean and second-moment arrays
#'
#' Ensures that a posterior second moment is finite and no smaller than the
#' squared posterior mean, up to numerical roundoff.
#'
#' @param mean Numeric vector or matrix of posterior means.
#' @param second Numeric vector or matrix of posterior second moments.
#' @param name Character label used in an error message.
#' @param tolerance Allowed negative variance from floating-point roundoff.
#' @return Invisibly TRUE.
#' @examples
#' multimodal_yfb_validate_moments(c(1, 2), c(1.1, 4.1), "EL")
#' @family multimodal_yfb_helpers
multimodal_yfb_validate_moments <- function(mean, second, name,
                                            tolerance = 1e-10) {
  if (!is.numeric(mean) || !is.numeric(second) || !identical(dim(mean), dim(second)) ||
      length(mean) != length(second)) {
    stop(name, " mean and second moments must be numeric objects of the same shape.")
  }
  if (any(!is.finite(mean)) || any(!is.finite(second))) {
    stop(name, " posterior moments must be finite.")
  }
  if (any(second < mean^2 - tolerance)) {
    stop(name, " posterior second moments must be at least squared means.")
  }
  invisible(TRUE)
}

#' Compute Breslow Cox diagonal working quantities
#'
#' Computes the score `u`, diagonal working curvature `w`, and linear
#' surrogate coefficient `h = u + w * eta` at a supplied mean predictor.
#'
#' @param time Positive numeric follow-up times.
#' @param event Numeric 0/1 event indicators.
#' @param eta Finite numeric mean linear predictor.
#' @return List with `u`, `w`, `h`, and Breslow `log_partial_likelihood`.
#' @examples
#' multimodal_yfb_cox_working(c(1, 2), c(1, 0), c(0, 0))
#' @family multimodal_yfb_helpers
multimodal_yfb_cox_working <- function(time, event, eta) {
  n <- length(time)
  if (!is.numeric(time) || !is.numeric(event) || !is.numeric(eta) ||
      length(event) != n || length(eta) != n || n == 0L) {
    stop("time, event, and eta must be non-empty numeric vectors of equal length.")
  }
  if (any(!is.finite(time)) || any(time <= 0) || any(!is.finite(eta)) ||
      any(!is.finite(event)) || any(!(event %in% c(0, 1)))) {
    stop("time must be positive and finite; event must be 0/1; eta must be finite.")
  }

  ord <- order(time)
  time_sorted <- time[ord]
  event_sorted <- event[ord]
  eta_sorted <- eta[ord]
  eta_shift <- max(eta_sorted)
  theta <- exp(eta_sorted - eta_shift)
  risk_sum_by_row <- rev(cumsum(rev(theta)))
  first_row <- match(time_sorted, time_sorted)
  risk_sum <- risk_sum_by_row[first_row]
  is_first <- !duplicated(time_sorted)
  event_count <- ave(event_sorted, time_sorted, FUN = sum)
  increment <- event_count[is_first] / risk_sum[is_first]
  cumulative_hazard <- cumsum(increment)[cumsum(is_first)]
  w_sorted <- theta * cumulative_hazard
  u_sorted <- event_sorted - w_sorted

  u <- w <- numeric(n)
  u[ord] <- u_sorted
  w[ord] <- w_sorted
  list(
    u = u,
    w = w,
    h = u + w * eta,
    log_partial_likelihood = sum(event_sorted * eta_sorted) -
      sum(event_count[is_first] * (log(risk_sum[is_first]) + eta_shift))
  )
}

#' Compute raw multimodal projection moments
#'
#' Computes `E[Z]`, `Var(Z)`, and `E[Z^2]` for the raw joint YFB projection.
#'
#' @param Y Named list of subject-by-feature numeric matrices.
#' @param EF Named list of feature-by-factor posterior mean matrices.
#' @param EF2 Named list of feature-by-factor posterior second-moment matrices.
#' @return List with `EZ`, `VZ`, and `EZ2`, each n by K.
#' @examples
#' Y <- list(expression = matrix(1, 2, 1), methylation = matrix(1, 2, 1))
#' F <- list(expression = matrix(1, 1, 1), methylation = matrix(1, 1, 1))
#' multimodal_yfb_projection_moments(Y, F, F)
#' @family multimodal_yfb_helpers
multimodal_yfb_projection_moments <- function(Y, EF, EF2) {
  modalities <- names(Y)
  if (is.null(modalities) || !identical(names(EF), modalities) ||
      !identical(names(EF2), modalities) || length(modalities) == 0L) {
    stop("Y, EF, and EF2 must be named lists with identical modalities.")
  }
  n <- nrow(Y[[1]])
  K <- ncol(EF[[1]])
  EZ <- VZ <- matrix(0, n, K)
  for (modality in modalities) {
    if (!is.matrix(Y[[modality]]) || !is.numeric(Y[[modality]]) ||
        !is.matrix(EF[[modality]]) || !is.numeric(EF[[modality]]) ||
        nrow(Y[[modality]]) != n || ncol(Y[[modality]]) != nrow(EF[[modality]]) ||
        ncol(EF[[modality]]) != K) {
      stop("Each modality must have conformable numeric Y and EF matrices.")
    }
    multimodal_yfb_validate_moments(EF[[modality]], EF2[[modality]],
                                    paste0("EF for ", modality))
    EZ <- EZ + Y[[modality]] %*% EF[[modality]]
    VZ <- VZ + (Y[[modality]]^2) %*%
      (EF2[[modality]] - EF[[modality]]^2)
  }
  if (any(VZ < -1e-10)) stop("Projection variances must be nonnegative.")
  VZ <- pmax(VZ, 0)
  list(EZ = EZ, VZ = VZ, EZ2 = EZ^2 + VZ)
}

#' Compute raw multimodal predictor moments
#'
#' Calculates the mean and variance of `eta = sum_k beta_k Z_k` under the
#' mean-field posterior.
#'
#' @param EZ n by K projection posterior means.
#' @param VZ n by K projection posterior variances.
#' @param EBeta Length-K coefficient posterior means.
#' @param EBeta2 Length-K coefficient posterior second moments.
#' @return List with numeric vectors `mean` and `variance`.
#' @examples
#' multimodal_yfb_predictor_moments(matrix(1, 2, 1), matrix(0, 2, 1), 1, 1)
#' @family multimodal_yfb_helpers
multimodal_yfb_predictor_moments <- function(EZ, VZ, EBeta, EBeta2) {
  if (!is.matrix(EZ) || !is.matrix(VZ) || !identical(dim(EZ), dim(VZ))) {
    stop("EZ and VZ must be matrices of the same shape.")
  }
  multimodal_yfb_validate_moments(EBeta, EBeta2, "EBeta")
  if (ncol(EZ) != length(EBeta) || any(!is.finite(EZ)) || any(!is.finite(VZ)) ||
      any(VZ < -1e-10)) {
    stop("Projection moments must be finite and conformable with beta moments.")
  }
  VZ <- pmax(VZ, 0)
  beta_variance <- EBeta2 - EBeta^2
  list(
    mean = as.vector(EZ %*% EBeta),
    variance = as.vector(VZ %*% EBeta2 + EZ^2 %*% beta_variance)
  )
}

#' Compute modality-specific residuals excluding one factor
#'
#' @param Y Named list of subject-by-feature matrices.
#' @param EL n by K shared score posterior means.
#' @param EF Named list of feature-by-factor loading posterior means.
#' @param k One-based factor index to omit.
#' @return Named list of subject-by-feature residual matrices.
#' @examples
#' multimodal_yfb_residual_minus_k(list(expression = matrix(1, 1, 1)),
#'                                 matrix(1, 1, 1),
#'                                 list(expression = matrix(1, 1, 1)), 1)
#' @family multimodal_yfb_helpers
multimodal_yfb_residual_minus_k <- function(Y, EL, EF, k) {
  if (!is.matrix(EL) || !is.numeric(EL) || k < 1L || k > ncol(EL) ||
      k != as.integer(k) || !identical(names(Y), names(EF))) {
    stop("EL, EF, Y, and k must define a valid shared-factor residual.")
  }
  l_other <- EL
  l_other[, k] <- 0
  lapply(names(Y), function(modality) {
    if (!is.matrix(Y[[modality]]) || !is.matrix(EF[[modality]]) ||
        nrow(Y[[modality]]) != nrow(EL) || ncol(EF[[modality]]) != ncol(EL) ||
        ncol(Y[[modality]]) != nrow(EF[[modality]])) {
      stop("Each modality must be conformable with EL and EF.")
    }
    Y[[modality]] - l_other %*% t(EF[[modality]])
  }) |> stats::setNames(names(Y))
}
