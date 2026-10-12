# ============================================================
# Script: fit_multimodal_yfb.R
# Purpose: Fit the isolated matched expression-methylation raw joint-YFB model.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; yaml; multimodal YFB helper, preprocessing, and update modules
# ============================================================

source("code/multimodal_yfb_helpers.R")
source("code/preprocess_multimodal_yfb.R")
source("code/multimodal_yfb_updates.R")

#' Read the named defaults for the isolated multimodal YFB fitter
#'
#' @return Named list of fitter controls from `config/globals.yml`.
#' @examples
#' multimodal_yfb_default_control()
#' @family multimodal_yfb
multimodal_yfb_default_control <- function() {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The yaml package is required to read config/globals.yml.")
  }
  defaults <- yaml::read_yaml("config/globals.yml")$multimodal_yfb
  required <- c("max_outer", "max_inner", "tol", "pve_threshold",
                "prognostic_z_threshold")
  if (!is.list(defaults) || !all(required %in% names(defaults))) {
    stop("config/globals.yml must define all multimodal_yfb fitter defaults.")
  }
  defaults
}

#' Combine posterior moments with a damping weight
#'
#' @param old_mean,old_second Current posterior moments.
#' @param new_mean,new_second Proposed posterior moments.
#' @param damping Weight on the proposed moments, in `(0, 1]`.
#' @return Damped posterior `mean` and `second` moments.
#' @examples
#' multimodal_yfb_damp_moments(0, 0, 1, 1, 0.5)
#' @family multimodal_yfb
multimodal_yfb_damp_moments <- function(old_mean, old_second, new_mean,
                                        new_second, damping) {
  if (damping == 1) return(list(mean = new_mean, second = new_second))
  old_variance <- pmax(old_second - old_mean^2, 0)
  new_variance <- pmax(new_second - new_mean^2, 0)
  mean <- (1 - damping) * old_mean + damping * new_mean
  variance <- (1 - damping) * old_variance + damping * new_variance
  list(mean = mean, second = mean^2 + variance)
}

#' Calculate a scale-invariant change in an identifiable fitted quantity
#'
#' The raw joint-YFB model intentionally retains modality-specific molecular
#' scales.  A convergence check based only on an absolute difference would
#' therefore demand different numerical precision for two equivalent analyses
#' expressed in different units.  This quantity is
#' `max(abs(new - old)) / max(1, max(abs(old)), max(abs(new)))`.
#'
#' @param old,new Finite numeric vectors or matrices with identical dimensions.
#' @return Nonnegative relative maximum change.
#' @examples
#' multimodal_yfb_relative_change(c(10, -20), c(10.01, -20.02))
#' @family multimodal_yfb
multimodal_yfb_relative_change <- function(old, new) {
  if (!is.numeric(old) || !is.numeric(new) || !identical(dim(old), dim(new)) ||
      length(old) != length(new) || any(!is.finite(old)) || any(!is.finite(new))) {
    stop("old and new must be finite numeric objects with identical dimensions.")
  }
  max(abs(new - old)) / max(1, abs(old), abs(new))
}

#' Initialize raw-scale survival coefficients with a standardized Cox fit
#'
#' The multimodal model must retain its raw joint projection
#' `Z = sum_m Y_m F_m`.  For initialization only, each nonconstant projection
#' column is divided by its sample standard deviation before fitting Cox PH.
#' This keeps the numerical regression problem well scaled; the fitted
#' coefficients and covariance are then transformed back so that
#' `Z_raw %*% beta_raw` equals the standardized Cox linear predictor up to an
#' irrelevant constant.  It therefore changes neither the model nor the
#' production predictor scale.
#'
#' @param EZ n-by-K matrix of raw joint projection means.
#' @param time,event Training survival vectors with event coded 0/1.
#' @return List with raw-scale coefficient `mean`, second moment `second`, and
#'   the projection column standard deviations used by the warm start.
#' @examples
#' multimodal_yfb_cox_warm_start(cbind(1:6, 6:1), 1:6, rep(1, 6))
#' @family multimodal_yfb
multimodal_yfb_cox_warm_start <- function(EZ, time, event) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The survival package is required for the multimodal Cox warm start.")
  }
  if (!is.matrix(EZ) || !is.numeric(EZ) || nrow(EZ) == 0L || ncol(EZ) == 0L ||
      any(!is.finite(EZ)) || !is.numeric(time) || !is.numeric(event) ||
      length(time) != nrow(EZ) || length(event) != nrow(EZ) ||
      any(!is.finite(time)) || any(!event %in% c(0, 1))) {
    stop("Cox warm-start inputs must be finite conformable projections and survival vectors.")
  }
  scale <- apply(EZ, 2, stats::sd)
  active <- is.finite(scale) & scale > sqrt(.Machine$double.eps)
  if (!any(active)) {
    stop("Cox warm start requires at least one nonconstant projection column.")
  }

  standardized <- sweep(EZ[, active, drop = FALSE], 2, scale[active], "/")
  covariates <- as.data.frame(standardized)
  names(covariates) <- paste0("z", seq_len(ncol(covariates)))
  covariates$time <- time
  covariates$event <- event
  cox_fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ ., data = covariates,
                    x = FALSE, model = FALSE),
    error = function(error) {
      stop("Multimodal Cox warm start failed: ", conditionMessage(error))
    }
  )

  standardized_mean <- stats::coef(cox_fit)
  standardized_variance <- diag(stats::vcov(cox_fit))
  valid <- is.finite(standardized_mean) & is.finite(standardized_variance) &
    standardized_variance >= 0
  standardized_mean[!valid] <- 0
  standardized_variance[!valid] <- 0
  mean <- numeric(ncol(EZ))
  variance <- numeric(ncol(EZ))
  mean[active] <- standardized_mean / scale[active]
  variance[active] <- standardized_variance / scale[active]^2
  list(mean = mean, second = mean^2 + pmax(variance, 0), scale = scale)
}

#' Canonicalize multimodal factor scales without changing fitted products
#'
#' The raw multimodal parameterization is scale-exchangeable: for a positive
#' factor-specific scale `a_k`, replacing `(L_k, F_mk, beta_k)` by
#' `(a_k L_k, F_mk / a_k, a_k beta_k)` leaves both `L F_m^T` and
#' `(sum_m Y_m F_m) beta` unchanged.  Repeated coordinate updates can otherwise
#' move along this numerically unidentifiable direction indefinitely.  This
#' function fixes the representative by setting the sample standard deviation
#' of each nonconstant *raw joint projection* to one and transforms posterior
#' second moments and point-exponential rates under the same change of
#' variables. This controls the Cox coefficient scale without normalizing the
#' model predictor or omitting its survival-feedback terms.
#'
#' @param Y Named observed modality matrices used to form raw projections.
#' @param EL,EL2 Shared score posterior mean and second-moment matrices.
#' @param EF,EF2 Named loading posterior mean and second-moment matrices.
#' @param EBeta,EBeta2 Shared survival coefficient posterior moments.
#' @param prior_L Shared score point-exponential priors by factor.
#' @param prior_F Named loading priors by modality and factor (family point_exponential,
#'   point_laplace or normal).
#' @param factors Integer factor indices to canonicalize; defaults to all factors.
#' @return Canonicalized moments, priors, and the applied positive factor scales.
#' @examples
#' multimodal_yfb_canonicalize_factors(
#'   list(expression = matrix(1:2, 2, 1)), matrix(1, 1, 1), matrix(1, 1, 1),
#'   list(expression = matrix(2, 1, 1)),
#'   list(expression = matrix(4, 1, 1)), 1, 1,
#'   list(list(pi = 1, rate = 1, point_mass = FALSE)),
#'   list(expression = list(list(pi = 1, rate = 1, point_mass = FALSE)))
#' )
#' @family multimodal_yfb
multimodal_yfb_canonicalize_factors <- function(Y, EL, EL2, EF, EF2, EBeta,
                                                EBeta2, prior_L, prior_F,
                                                factors = seq_len(ncol(EL))) {
  if (!is.matrix(EL) || !is.matrix(EL2) || !identical(dim(EL), dim(EL2)) ||
      !identical(names(Y), names(EF)) || !identical(names(EF), names(EF2)) ||
      !identical(names(EF), names(prior_F))) {
    stop("Canonicalization requires conformable named posterior and prior blocks.")
  }
  K <- ncol(EL)
  if (length(EBeta) != K || length(EBeta2) != K || length(prior_L) != K ||
      any(!vapply(EF, function(x) is.matrix(x) && ncol(x) == K, logical(1))) ||
      any(!vapply(EF2, function(x) is.matrix(x) && ncol(x) == K, logical(1))) ||
      any(!vapply(prior_F, function(x) length(x) == K, logical(1)))) {
    stop("Canonicalization inputs must have the same factor count.")
  }
  if (!is.numeric(factors) || length(factors) == 0L || any(!is.finite(factors)) ||
      any(factors != as.integer(factors)) || any(factors < 1L | factors > K)) {
    stop("factors must be nonempty integer indices in the fitted factor range.")
  }
  factors <- unique(as.integer(factors))
  multimodal_yfb_validate_moments(EL, EL2, "EL")
  for (modality in names(EF)) {
    multimodal_yfb_validate_moments(EF[[modality]], EF2[[modality]],
                                    paste0("EF for ", modality))
  }
  multimodal_yfb_validate_moments(EBeta, EBeta2, "EBeta")

  scale <- rep(1, K)
  for (k in factors) {
    projection_k <- Reduce(`+`, lapply(names(Y), function(modality) {
      Y[[modality]] %*% EF[[modality]][, k]
    }))
    scale[k] <- stats::sd(projection_k)
  }
  scale[scale <= sqrt(.Machine$double.eps)] <- 1
  for (k in factors) {
    EL[, k] <- EL[, k] * scale[k]
    EL2[, k] <- EL2[, k] * scale[k]^2
    EBeta[k] <- EBeta[k] * scale[k]
    EBeta2[k] <- EBeta2[k] * scale[k]^2
    if (!isTRUE(prior_L[[k]]$point_mass)) prior_L[[k]]$rate <- prior_L[[k]]$rate / scale[k]
    for (modality in names(EF)) {
      EF[[modality]][, k] <- EF[[modality]][, k] / scale[k]
      EF2[[modality]][, k] <- EF2[[modality]][, k] / scale[k]^2
      if (!isTRUE(prior_F[[modality]][[k]]$point_mass)) {
        # F_k -> F_k / c rescales an exponential or Laplace rate by c and a
        # Normal prior variance by 1 / c^2
        if (identical(prior_F[[modality]][[k]]$family, "normal")) {
          prior_F[[modality]][[k]]$variance <- prior_F[[modality]][[k]]$variance / scale[k]^2
        } else {
          prior_F[[modality]][[k]]$rate <- prior_F[[modality]][[k]]$rate * scale[k]
        }
      }
    }
  }
  list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2, EBeta = EBeta,
       EBeta2 = EBeta2, prior_L = prior_L, prior_F = prior_F, scale = scale)
}

#' Compute expected reconstruction residuals in feature chunks
#'
#' This is algebraically identical to `multimodal_yfb_update_tau()`, but avoids
#' materializing an n-by-p expected-residual matrix for a complete modality.
#'
#' @param Y Named list of observed modality matrices.
#' @param EL,EL2 Shared score posterior moments.
#' @param EF,EF2 Named loading posterior moments.
#' @param chunk_size Number of features processed per chunk.
#' @return Modality-specific precision vectors and expected residual sums.
#' @examples
#' multimodal_yfb_update_tau_chunked(list(expression = matrix(1, 1, 1)),
#'   matrix(1, 1, 1), matrix(1, 1, 1), list(expression = matrix(0, 1, 1)),
#'   list(expression = matrix(0, 1, 1)), chunk_size = 1)
#' @family multimodal_yfb
multimodal_yfb_update_tau_chunked <- function(Y, EL, EL2, EF, EF2,
                                              chunk_size) {
  if (length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer.")
  }
  Tau <- expected_residual2 <- vector("list", length(Y))
  names(Tau) <- names(expected_residual2) <- names(Y)
  for (modality in names(Y)) {
    p <- ncol(Y[[modality]])
    residual <- numeric(p)
    for (first in seq.int(1L, p, by = as.integer(chunk_size))) {
      index <- first:min(p, first + as.integer(chunk_size) - 1L)
      mean_signal <- EL %*% t(EF[[modality]][index, , drop = FALSE])
      variance_signal <- EL2 %*% t(EF2[[modality]][index, , drop = FALSE]) -
        EL^2 %*% t(EF[[modality]][index, , drop = FALSE]^2)
      residual[index] <- colSums((Y[[modality]][, index, drop = FALSE] -
        mean_signal)^2 + pmax(variance_signal, 0))
    }
    expected_residual2[[modality]] <- residual
    Tau[[modality]] <- ifelse(residual > 0, nrow(EL) / residual, Inf)
  }
  list(Tau = Tau, expected_residual2 = expected_residual2)
}

#' Calculate a fixed-working quadratic objective for iteration monitoring
#'
#' @param Y Named observed modality matrices.
#' @param EL,EL2 Shared score posterior moments.
#' @param EF,EF2 Named loading posterior moments.
#' @param Tau Modality-specific feature precisions.
#' @param EZ,VZ Projection posterior moments.
#' @param EBeta,EBeta2 Survival coefficient posterior moments.
#' @param w,h Frozen Cox working quantities.
#' @return Scalar reconstruction plus Cox working surrogate objective.
#' @examples
#' # Used internally by fit_multimodal_yfb().
#' @family multimodal_yfb
multimodal_yfb_fixed_working_objective <- function(Y, EL, EL2, EF, EF2, Tau,
                                                   EZ, VZ, EBeta, EBeta2, w, h) {
  tau_state <- multimodal_yfb_update_tau_chunked(Y, EL, EL2, EF, EF2,
                                                  chunk_size = 1000L)
  reconstruction <- sum(vapply(names(Y), function(modality) {
    sum(log(Tau[[modality]]) - Tau[[modality]] *
      tau_state$expected_residual2[[modality]] / nrow(EL))
  }, numeric(1))) * nrow(EL) / 2
  eta <- EZ %*% EBeta
  eta_variance <- VZ %*% EBeta2 + EZ^2 %*% EBeta2 - (EZ %*% EBeta)^2
  reconstruction + sum(h * eta - 0.5 * w * (eta^2 + pmax(eta_variance, 0)))
}

#' Fit the isolated raw-projection multimodal YFB model
#'
#' The fitter alternates exact modal updates under a fixed Cox working
#' approximation. It deliberately does not alter either established single-
#' modality YFB implementation.
#'
#' @param Y Named matched expression and methylation matrices.
#' @param time,event Named training survival vectors.
#' @param K Over-specified initial factor count.
#' @param control Named list overriding entries in `multimodal_yfb` config. In
#'   addition to configured defaults it accepts `damping` and `tau_chunk_size`.
#'   Factor pruning: `prior_update` ("em" or "ebnm" marginal-likelihood prior
#'   fits for point-exponential L and F), `prune` (logical: after convergence,
#'   remove any factor whose removal does not lower the ELBO, then refit), and
#'   `survival_elbo` ("plugin": Breslow partial log-likelihood at E[eta];
#'   "corrected": minus (1/2) sum_i w_i Var(eta_i), a second-order
#'   approximation to E_q[log PL]).
#' @param init Optional warm-start state (EL, EL2, EF, EF2, EBeta, EBeta2, Tau,
#'   prior_L, prior_F), e.g. a fit with one factor removed. Used by pruning.
#' @return A `multimodal_yfb_fit` object with posterior moments and diagnostics.
#' @examples
#' # fit_multimodal_yfb(Y, time, event, K = 7)
#' @family multimodal_yfb
#' @seealso [predict_multimodal_yfb()]
fit_multimodal_yfb <- function(Y, time, event, K, control = list(), init = NULL) {
  if (!is.list(control) || (length(control) > 0L && is.null(names(control)))) {
    stop("control must be an empty or named list.")
  }
  # Signed loading priors allow signed (e.g. centered) data in that modality
  prior_F_requested <- modifyList(multimodal_yfb_default_control(), control)$prior_F
  nonnegative <- if (is.list(prior_F_requested) && !is.null(names(prior_F_requested))) {
    vapply(prior_F_requested, function(f) identical(f, "point_exponential"), logical(1))
  } else NULL
  if (!is.null(nonnegative) && is.list(Y)) nonnegative <- nonnegative[intersect(names(nonnegative), names(Y))]
  data <- preprocess_multimodal_yfb_training(Y, time, event, nonnegative)
  Y <- data$Y
  Y_raw <- Y
  n <- nrow(Y[[1]])
  if (length(K) != 1L || !is.finite(K) || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer.")
  }
  settings <- modifyList(multimodal_yfb_default_control(), control)
  settings$damping <- settings$damping %||% 1
  settings$prior_update <- settings$prior_update %||% "em"
  settings$prune <- isTRUE(settings$prune)
  settings$survival_elbo <- settings$survival_elbo %||% "plugin"
  settings$tau_model <- settings$tau_model %||% "feature"
  # Intercept: Y_m = 1 mu_m^T + L F_m^T + E_m, with mu_m re-estimated each outer
  # sweep as colMeans(Y_m - L F_m^T). F then describes deviations from the
  # baseline, so it should normally have a signed prior. A constant shift of
  # eta does not change the Cox partial likelihood, so survival is unaffected.
  settings$intercept <- isTRUE(settings$intercept)
  mu <- NULL
  if (settings$intercept) {
    mu <- lapply(Y_raw, colMeans)
    Y <- Map(function(x, m) sweep(x, 2, m), Y_raw, mu)
  }
  if (!settings$tau_model %in% c("feature", "modality", "feature_eb")) {
    stop("control$tau_model must be feature, modality or feature_eb.")
  }
  if (!settings$prior_update %in% c("em", "ebnm") ||
      !settings$survival_elbo %in% c("plugin", "corrected")) {
    stop("control$prior_update must be em or ebnm; control$survival_elbo plugin or corrected.")
  }
  settings$prior_F <- settings$prior_F %||%
    list(expression = "point_exponential", methylation = "point_exponential")
  # Keep the prior entries for the modalities actually supplied
  if (is.list(settings$prior_F)) settings$prior_F <- settings$prior_F[intersect(names(settings$prior_F), names(Y))]
  if (!is.list(settings$prior_F) || !setequal(names(settings$prior_F), names(Y)) ||
      !all(unlist(settings$prior_F) %in% c("point_exponential", "point_laplace", "normal"))) {
    stop("control$prior_F must name a prior (point_exponential, point_laplace or normal) for each modality.")
  }
  settings$tau_chunk_size <- settings$tau_chunk_size %||% 1000L
  numeric_controls <- c("max_outer", "max_inner", "tol", "pve_threshold",
                        "prognostic_z_threshold", "damping", "tau_chunk_size")
  if (any(!vapply(settings[numeric_controls], function(x) {
    length(x) == 1L && is.finite(x)
  }, logical(1))) || settings$max_outer < 1L || settings$max_inner < 1L ||
      settings$tol <= 0 || settings$pve_threshold < 0 ||
      settings$prognostic_z_threshold < 0 || settings$damping <= 0 ||
      settings$damping > 1 || settings$tau_chunk_size < 1L) {
    stop("Multimodal YFB controls must have valid positive scalar values.")
  }

  if (!is.null(init)) {
    # Warm start (e.g. a pruned fit): take every moment and prior from init
    K <- ncol(init$EL)
    EL <- init$EL; EL2 <- init$EL2; EF <- init$EF; EF2 <- init$EF2
    EBeta <- init$EBeta; EBeta2 <- init$EBeta2; Tau <- init$Tau
    prior_L <- init$prior_L; prior_F <- init$prior_F
  } else {
    decomposition <- svd(do.call(cbind, Y), nu = min(n, K), nv = min(n, K))
    K <- min(K, ncol(decomposition$u))
    EL <- abs(sweep(decomposition$u[, seq_len(K), drop = FALSE], 2,
                    decomposition$d[seq_len(K)], "*"))
    EF <- list()
    first <- 1L
    for (modality in names(Y)) {
      p <- ncol(Y[[modality]])
      EF[[modality]] <- abs(decomposition$v[first:(first + p - 1L), seq_len(K),
                                             drop = FALSE])
      first <- first + p
    }
    EL2 <- EL^2 + 1e-6
    EF2 <- lapply(EF, function(x) x^2 + 1e-6)
    Tau <- lapply(Y, function(x) rep(1 / max(stats::var(as.vector(x)), 1e-6), ncol(x)))
    prior_L <- rep(list(list(pi = 0.5, rate = 1, point_mass = FALSE)), K)
    # Loading prior family per modality (control$prior_F, default from globals):
    # point_exponential (nonnegative F), point_laplace or normal (signed F)
    prior_F <- lapply(names(Y), function(modality) {
      family <- settings$prior_F[[modality]]
      init <- switch(family,
        point_exponential = list(family = family, pi = 0.5, rate = 1, point_mass = FALSE),
        point_laplace = list(family = family, pi = 0.5, rate = 1, point_mass = FALSE),
        normal = list(family = family, variance = 1, point_mass = FALSE))
      rep(list(init), K)
    })
    names(prior_F) <- names(Y)

    projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
    beta_warm_start <- multimodal_yfb_cox_warm_start(
      projection$EZ, data$time, data$event
    )
    EBeta <- beta_warm_start$mean
    EBeta2 <- pmax(beta_warm_start$second, EBeta^2)
  }
  # Per-factor KL(q || g) terms of the ELBO, refreshed at each factor update
  kl_L <- kl_beta <- numeric(K)
  kl_F <- lapply(Y, function(x) numeric(K))

  history <- vector("list", as.integer(settings$max_outer))
  converged <- FALSE
  for (outer in seq_len(as.integer(settings$max_outer))) {
    if (settings$intercept) {
      mu <- lapply(names(Y_raw), function(m) colMeans(Y_raw[[m]] - EL %*% t(EF[[m]])))
      names(mu) <- names(Y_raw)
      Y <- Map(function(x, m) sweep(x, 2, m), Y_raw, mu)
    }
    projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
    predictor <- multimodal_yfb_predictor_moments(projection$EZ, projection$VZ,
                                                   EBeta, EBeta2)
    working <- multimodal_yfb_cox_working(data$time, data$event, predictor$mean)
    old_eta <- predictor$mean
    old_reconstruction <- lapply(names(Y), function(modality) {
      EL %*% t(EF[[modality]])
    })
    names(old_reconstruction) <- names(Y)
    old_moments <- c(EL, unlist(EF, use.names = FALSE), EBeta)

    for (inner in seq_len(as.integer(settings$max_inner))) {
      for (k in seq_len(K)) {
        residual <- multimodal_yfb_residual_minus_k(Y, EL, EF, k)
        updated_L <- multimodal_yfb_update_L_k(
          Y, residual, lapply(EF, `[`, , k), lapply(EF2, `[`, , k), Tau,
          prior_L[[k]], prior_update = settings$prior_update
        )
        kl_L[k] <- updated_L$kl
        damped_L <- multimodal_yfb_damp_moments(EL[, k], EL2[, k],
                                                 updated_L$mean, updated_L$second,
                                                 settings$damping)
        EL[, k] <- damped_L$mean
        EL2[, k] <- damped_L$second
        prior_L[[k]] <- updated_L$prior

        for (modality in names(Y)) {
          residual <- multimodal_yfb_residual_minus_k(Y, EL, EF, k)
          projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
          current_eta <- multimodal_yfb_predictor_moments(
            projection$EZ, projection$VZ, EBeta, EBeta2
          )$mean
          h_minus_k <- working$h - working$w *
            (current_eta - EBeta[k] * projection$EZ[, k])
          updated_F <- multimodal_yfb_update_F_mk(
            Y[[modality]], Tau[[modality]], EL[, k], EL2[, k], residual[[modality]],
            working$w, h_minus_k, projection$EZ[, k], projection$VZ[, k],
            EF[[modality]][, k], EF2[[modality]][, k], EBeta[k], EBeta2[k],
            prior_F[[modality]][[k]], prior_update = settings$prior_update
          )
          kl_F[[modality]][k] <- updated_F$kl
          damped_F <- multimodal_yfb_damp_moments(
            EF[[modality]][, k], EF2[[modality]][, k], updated_F$EF,
            updated_F$EF2, settings$damping
          )
          EF[[modality]][, k] <- damped_F$mean
          EF2[[modality]][, k] <- damped_F$second
          prior_F[[modality]][[k]] <- updated_F$prior
        }

        projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
        current_eta <- multimodal_yfb_predictor_moments(
          projection$EZ, projection$VZ, EBeta, EBeta2
        )$mean
        h_minus_k <- working$h - working$w *
          (current_eta - EBeta[k] * projection$EZ[, k])
        updated_beta <- multimodal_yfb_update_beta_k(
          projection$EZ[, k], projection$VZ[, k], working$w, h_minus_k
        )
        kl_beta[k] <- updated_beta$kl
        damped_beta <- multimodal_yfb_damp_moments(
          EBeta[k], EBeta2[k], updated_beta$mean, updated_beta$second,
          settings$damping
        )
        EBeta[k] <- damped_beta$mean
        EBeta2[k] <- damped_beta$second
        canonical <- multimodal_yfb_canonicalize_factors(
          Y, EL, EL2, EF, EF2, EBeta, EBeta2, prior_L, prior_F, factors = k
        )
        EL <- canonical$EL
        EL2 <- canonical$EL2
        EF <- canonical$EF
        EF2 <- canonical$EF2
        EBeta <- canonical$EBeta
        EBeta2 <- canonical$EBeta2
        prior_L <- canonical$prior_L
        prior_F <- canonical$prior_F
      }
    }

    tau_state <- multimodal_yfb_update_tau_chunked(
      Y, EL, EL2, EF, EF2, settings$tau_chunk_size
    )
    Tau <- multimodal_yfb_tau_from_residuals(tau_state$expected_residual2, n,
                                             settings$tau_model)$Tau
    projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
    predictor <- multimodal_yfb_predictor_moments(projection$EZ, projection$VZ,
                                                   EBeta, EBeta2)
    new_moments <- c(EL, unlist(EF, use.names = FALSE), EBeta)
    reconstruction_change <- vapply(names(Y), function(modality) {
      multimodal_yfb_relative_change(
        old_reconstruction[[modality]], EL %*% t(EF[[modality]])
      )
    }, numeric(1))
    history[[outer]] <- list(
      fixed_working_objective = multimodal_yfb_fixed_working_objective(
        Y, EL, EL2, EF, EF2, Tau, projection$EZ, projection$VZ, EBeta,
        EBeta2, working$w, working$h
      ),
      partial_log_likelihood = multimodal_yfb_cox_working(
        data$time, data$event, predictor$mean
      )$log_partial_likelihood,
      max_moment_change = max(abs(new_moments - old_moments)),
      elbo = multimodal_yfb_elbo(Y, data$time, data$event, EL, EL2, EF, EF2,
                                 EBeta, EBeta2,
                                 sum(kl_L) + sum(unlist(kl_F)) + sum(kl_beta),
                                 settings$survival_elbo, settings$tau_model),
      kl_parts = c(L = sum(kl_L), F = sum(unlist(kl_F)), beta = sum(kl_beta)),
      max_relative_reconstruction_change = max(reconstruction_change),
      max_relative_eta_change = multimodal_yfb_relative_change(old_eta, predictor$mean)
    )
    # Posterior moments can move along equivalent latent-factor allocations even
    # after the fitted data and raw predictor stabilize.  Convergence is therefore
    # assessed on the identifiable fitted reconstruction and η = Σ_m Y_m F_m β.
    if (history[[outer]]$max_relative_reconstruction_change < settings$tol &&
        history[[outer]]$max_relative_eta_change < settings$tol) {
      converged <- TRUE
      history <- history[seq_len(outer)]
      break
    }
  }

  pve <- vapply(seq_len(K), function(k) {
    numerator <- sum(vapply(names(Y), function(modality) {
      sum((EL[, k] %o% EF[[modality]][, k])^2)
    }, numeric(1)))
    numerator / sum(vapply(Y, function(x) sum(x^2), numeric(1)))
  }, numeric(1))
  beta_sd <- sqrt(pmax(EBeta2 - EBeta^2, 0))
  reconstruction_active <- pve >= settings$pve_threshold
  prognostic_z <- abs(EBeta) / pmax(beta_sd, 1e-12)
  survival_active <- prognostic_z >= settings$prognostic_z_threshold
  factor_class <- ifelse(reconstruction_active & survival_active, "both",
    ifelse(reconstruction_active, "reconstruction_only",
      ifelse(survival_active, "survival_only", "inactive")))
  diagnostics <- list(
    iterations = length(history), converged = converged, history = history,
    objective = vapply(history, `[[`, numeric(1), "partial_log_likelihood"),
    fixed_working_objective = vapply(history, `[[`, numeric(1),
                                     "fixed_working_objective"),
    factor_pve = pve,
    reconstruction_active = reconstruction_active,
    prognostic_z = prognostic_z,
    prognostic = survival_active,
    survival_active = survival_active,
    factor_class = factor_class,
    K_eff = sum(reconstruction_active),
    K_eff_reconstruction = sum(reconstruction_active),
    K_eff_survival = sum(survival_active),
    K_eff_retained = sum(reconstruction_active | survival_active),
    controls = settings
  )
  kl <- list(L = kl_L, F = kl_F, beta = kl_beta)
  diagnostics$elbo <- multimodal_yfb_elbo(Y, data$time, data$event, EL, EL2, EF, EF2,
                                          EBeta, EBeta2,
                                          sum(kl_L) + sum(unlist(kl_F)) + sum(kl_beta),
                                          settings$survival_elbo, settings$tau_model)
  out <- list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2, EBeta = EBeta,
              EBeta2 = EBeta2, Tau = Tau, prior_L = prior_L, prior_F = prior_F,
              kl = kl, mu = mu, training_spec = data$training_spec,
              diagnostics = diagnostics)
  class(out) <- "multimodal_yfb_fit"

  # Factor pruning (flashier-style nullcheck) ----
  # Remove the factor whose removal raises the ELBO most (if any does not
  # lower it), refit from the remaining factors, and repeat.
  if (settings$prune && K > 1L) {
    delta <- multimodal_yfb_nullcheck(out, Y_raw, data$time, data$event,
                                      settings$survival_elbo, settings$tau_model,
                                      settings$intercept)
    if (max(delta) >= 0) {
      drop <- which.max(delta)
      keep <- setdiff(seq_len(K), drop)
      pruned_init <- list(
        EL = out$EL[, keep, drop = FALSE], EL2 = out$EL2[, keep, drop = FALSE],
        EF = lapply(out$EF, function(x) x[, keep, drop = FALSE]),
        EF2 = lapply(out$EF2, function(x) x[, keep, drop = FALSE]),
        EBeta = out$EBeta[keep], EBeta2 = out$EBeta2[keep], Tau = out$Tau,
        prior_L = out$prior_L[keep],
        prior_F = lapply(out$prior_F, function(x) x[keep]))
      refit <- fit_multimodal_yfb(Y_raw, time, event, K - 1L, control, init = pruned_init)
      refit$diagnostics$pruning <- rbind(
        data.frame(K_before = K, removed_factor = drop, elbo_gain = delta[drop],
                   converged_before = converged),
        out$diagnostics$pruning, refit$diagnostics$pruning)
      return(refit)
    }
  }
  out
}

#' Residual precisions and their ELBO term under three precision models
#'
#' r2_mj = E_q||y_mj - L f_mj||^2 is the expected residual sum of squares of
#' feature j in modality m (n subjects).
#' - "feature": tau_mj = n / r2_mj (maximum likelihood; the original model).
#'   Unbounded: a factor that fits one feature exactly sends tau_mj to
#'   infinity, which keeps spare factors alive.
#' - "modality": one tau_m = n p_m / sum_j r2_mj per modality.
#' - "feature_eb": tau_mj ~ Gamma(a_m, b_m), with (a_m, b_m) fit by maximum
#'   marginal likelihood within each modality. Conjugacy gives
#'   q(tau_mj) = Gamma(a_m + n/2, b_m + r2_mj/2), E[tau_mj] = (a_m + n/2) / (b_m + r2_mj/2),
#'   and the reconstruction ELBO term per feature is the log marginal
#'   lgamma(a + n/2) - lgamma(a) + a log b - (a + n/2) log(b + r2/2) - (n/2) log(2 pi).
#'
#' @param r2 Named list (by modality) of expected residual sums of squares.
#' @param n Number of subjects.
#' @param tau_model "feature", "modality" or "feature_eb".
#' @return List with `Tau` (named list of precision vectors, E[tau] for
#'   feature_eb), `reconstruction` (scalar ELBO term), and `hyper` (Gamma
#'   parameters for feature_eb).
#' @examples
#' multimodal_yfb_tau_from_residuals(list(expression = c(10, 12, 0.5)), 10, "feature_eb")
#' @family multimodal_yfb
multimodal_yfb_tau_from_residuals <- function(r2, n, tau_model = "feature") {
  out <- lapply(names(r2), function(modality) {
    r <- pmax(r2[[modality]], .Machine$double.eps)
    p <- length(r)
    if (tau_model == "feature") {
      tau <- n / r
      return(list(Tau = tau, reconstruction = sum(0.5 * n * log(tau / (2 * pi)) - 0.5 * tau * r),
                  hyper = NULL))
    }
    if (tau_model == "modality") {
      tau <- n * p / sum(r)
      return(list(Tau = rep(tau, p),
                  reconstruction = 0.5 * n * p * log(tau / (2 * pi)) - 0.5 * tau * sum(r),
                  hyper = NULL))
    }
    log_marginal <- function(theta) {
      a <- exp(theta[1]); b <- exp(theta[2])
      sum(lgamma(a + n / 2) - lgamma(a) + a * log(b) - (a + n / 2) * log(b + r / 2) -
            0.5 * n * log(2 * pi))
    }
    # start at the moment estimate: tau_j ~ n / r_j has mean m and variance v
    tau_hat <- n / r
    m <- mean(tau_hat); v <- max(stats::var(tau_hat), 1e-8 * m^2)
    start <- log(c(m^2 / v, m / v))
    opt <- stats::optim(start, function(th) -log_marginal(th), method = "BFGS")
    a <- exp(opt$par[1]); b <- exp(opt$par[2])
    list(Tau = (a + n / 2) / (b + r / 2), reconstruction = -opt$value,
         hyper = c(a = a, b = b))
  })
  names(out) <- names(r2)
  list(Tau = lapply(out, `[[`, "Tau"),
       reconstruction = sum(vapply(out, `[[`, numeric(1), "reconstruction")),
       hyper = lapply(out, `[[`, "hyper"))
}

#' Evidence lower bound of a multimodal YFB fit
#'
#' ELBO = sum_m sum_j [ (n/2) log(tau_mj / 2 pi) - (tau_mj / 2) E||y_mj - L f_mj||^2 ]
#'        + S(eta) - KL,
#' with tau at its maximum-likelihood value given the expected residuals, KL
#' the summed coordinate KL(q || g) for L, F and beta, and S the survival
#' term: the Breslow partial log-likelihood at E[eta] ("plugin"), or that
#' minus (1/2) sum_i w_i Var(eta_i) ("corrected"; a second-order
#' approximation to E_q[log PL], with w the Cox working curvature).
#'
#' @param Y Named observed modality matrices.
#' @param time,event Training survival outcomes.
#' @param EL,EL2,EF,EF2,EBeta,EBeta2 Posterior moments.
#' @param kl_total Summed KL divergence of all coordinates.
#' @param survival_elbo "plugin" or "corrected".
#' @param tau_model Residual precision model; see [multimodal_yfb_tau_from_residuals()].
#' @return Scalar ELBO with a `components` attribute.
#' @examples
#' # multimodal_yfb_elbo(Y, time, event, fit$EL, fit$EL2, fit$EF, fit$EF2,
#' #                     fit$EBeta, fit$EBeta2, 0, "plugin")
#' @family multimodal_yfb
multimodal_yfb_elbo <- function(Y, time, event, EL, EL2, EF, EF2, EBeta, EBeta2,
                                kl_total, survival_elbo = "plugin",
                                tau_model = "feature") {
  n <- nrow(EL)
  tau_state <- multimodal_yfb_update_tau_chunked(Y, EL, EL2, EF, EF2, chunk_size = 1000L)
  reconstruction <- multimodal_yfb_tau_from_residuals(
    tau_state$expected_residual2, n, tau_model)$reconstruction
  projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
  predictor <- multimodal_yfb_predictor_moments(projection$EZ, projection$VZ, EBeta, EBeta2)
  working <- multimodal_yfb_cox_working(time, event, predictor$mean)
  survival <- working$log_partial_likelihood
  if (survival_elbo == "corrected") {
    survival <- survival - 0.5 * sum(working$w * pmax(predictor$variance, 0))
  }
  structure(reconstruction + survival - kl_total,
            components = c(reconstruction = reconstruction, survival = survival,
                           kl = kl_total))
}

#' ELBO change from removing each factor (flashier-style nullcheck)
#'
#' For each factor k, sets its scores, loadings and coefficient to zero (their
#' KL terms then vanish), re-estimates tau, and recomputes the ELBO with the
#' other factors held fixed. Holding the others fixed makes this a
#' conservative test: a positive change means the data prefer the fit without
#' factor k even before the remaining factors adapt.
#'
#' @param fit A `multimodal_yfb_fit` (with `kl`).
#' @param Y,time,event Training data used for the fit.
#' @param survival_elbo "plugin" or "corrected".
#' @param tau_model Residual precision model.
#' @param intercept Whether the fit has a per-feature intercept (then `Y` is raw).
#' @return Numeric vector: ELBO(without k) - ELBO(with all factors), per k.
#' @examples
#' # multimodal_yfb_nullcheck(fit, Y, time, event)
#' @family multimodal_yfb
multimodal_yfb_nullcheck <- function(fit, Y, time, event, survival_elbo = "plugin",
                                     tau_model = "feature", intercept = FALSE) {
  K <- ncol(fit$EL)
  kl_k <- fit$kl$L + Reduce(`+`, fit$kl$F) + fit$kl$beta
  # With an intercept, Y is the raw data and the intercept is re-estimated for
  # each candidate set of factors (removing a nonnegative-score factor also
  # moves its mean into the intercept)
  centered <- function(EL, EF) {
    if (!intercept) return(Y)
    Map(function(x, f) sweep(x, 2, colMeans(x - EL %*% t(f))), Y, EF)
  }
  full <- multimodal_yfb_elbo(centered(fit$EL, fit$EF), time, event, fit$EL, fit$EL2,
                              fit$EF, fit$EF2, fit$EBeta, fit$EBeta2, sum(kl_k),
                              survival_elbo, tau_model)
  vapply(seq_len(K), function(k) {
    keep <- setdiff(seq_len(K), k)
    if (length(keep) == 0L) return(-Inf)
    EL_keep <- fit$EL[, keep, drop = FALSE]
    EF_keep <- lapply(fit$EF, function(x) x[, keep, drop = FALSE])
    reduced <- multimodal_yfb_elbo(
      centered(EL_keep, EF_keep), time, event, EL_keep, fit$EL2[, keep, drop = FALSE],
      lapply(fit$EF, function(x) x[, keep, drop = FALSE]),
      lapply(fit$EF2, function(x) x[, keep, drop = FALSE]),
      fit$EBeta[keep], fit$EBeta2[keep], sum(kl_k[keep]), survival_elbo, tau_model)
    reduced - full
  }, numeric(1))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
