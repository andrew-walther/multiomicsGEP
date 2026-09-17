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
#' @return A `multimodal_yfb_fit` object with posterior moments and diagnostics.
#' @examples
#' # fit_multimodal_yfb(Y, time, event, K = 7)
#' @family multimodal_yfb
#' @seealso [predict_multimodal_yfb()]
fit_multimodal_yfb <- function(Y, time, event, K, control = list()) {
  data <- preprocess_multimodal_yfb_training(Y, time, event)
  Y <- data$Y
  n <- nrow(Y$expression)
  if (length(K) != 1L || !is.finite(K) || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer.")
  }
  if (!is.list(control) || is.null(names(control))) {
    stop("control must be a named list.")
  }
  settings <- modifyList(multimodal_yfb_default_control(), control)
  settings$damping <- settings$damping %||% 1
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
  prior_F <- lapply(Y, function(x) rep(list(list(pi = 0.5, rate = 1,
                                               point_mass = FALSE)), K))

  projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
  initial_working <- multimodal_yfb_cox_working(data$time, data$event, rep(0, n))
  EBeta <- vapply(seq_len(K), function(k) {
    multimodal_yfb_update_beta_k(projection$EZ[, k], projection$VZ[, k],
                                 initial_working$w, initial_working$h)$mean
  }, numeric(1))
  EBeta2 <- pmax(EBeta^2, 1e-12)

  history <- vector("list", as.integer(settings$max_outer))
  converged <- FALSE
  for (outer in seq_len(as.integer(settings$max_outer))) {
    projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
    predictor <- multimodal_yfb_predictor_moments(projection$EZ, projection$VZ,
                                                   EBeta, EBeta2)
    working <- multimodal_yfb_cox_working(data$time, data$event, predictor$mean)
    old_eta <- predictor$mean
    old_moments <- c(EL, unlist(EF, use.names = FALSE), EBeta)

    for (inner in seq_len(as.integer(settings$max_inner))) {
      for (k in seq_len(K)) {
        residual <- multimodal_yfb_residual_minus_k(Y, EL, EF, k)
        updated_L <- multimodal_yfb_update_L_k(
          Y, residual, lapply(EF, `[`, , k), lapply(EF2, `[`, , k), Tau,
          prior_L[[k]]
        )
        damped_L <- multimodal_yfb_damp_moments(EL[, k], EL2[, k],
                                                 updated_L$mean, updated_L$second,
                                                 settings$damping)
        EL[, k] <- damped_L$mean
        EL2[, k] <- damped_L$second

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
            prior_F[[modality]][[k]]
          )
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
        damped_beta <- multimodal_yfb_damp_moments(
          EBeta[k], EBeta2[k], updated_beta$mean, updated_beta$second,
          settings$damping
        )
        EBeta[k] <- damped_beta$mean
        EBeta2[k] <- damped_beta$second
      }
    }

    tau_state <- multimodal_yfb_update_tau_chunked(
      Y, EL, EL2, EF, EF2, settings$tau_chunk_size
    )
    Tau <- tau_state$Tau
    projection <- multimodal_yfb_projection_moments(Y, EF, EF2)
    predictor <- multimodal_yfb_predictor_moments(projection$EZ, projection$VZ,
                                                   EBeta, EBeta2)
    new_moments <- c(EL, unlist(EF, use.names = FALSE), EBeta)
    history[[outer]] <- list(
      fixed_working_objective = multimodal_yfb_fixed_working_objective(
        Y, EL, EL2, EF, EF2, Tau, projection$EZ, projection$VZ, EBeta,
        EBeta2, working$w, working$h
      ),
      partial_log_likelihood = multimodal_yfb_cox_working(
        data$time, data$event, predictor$mean
      )$log_partial_likelihood,
      max_moment_change = max(abs(new_moments - old_moments)),
      max_eta_change = max(abs(predictor$mean - old_eta))
    )
    if (history[[outer]]$max_moment_change < settings$tol &&
        history[[outer]]$max_eta_change < settings$tol) {
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
  diagnostics <- list(
    iterations = length(history), converged = converged, history = history,
    objective = vapply(history, `[[`, numeric(1), "partial_log_likelihood"),
    fixed_working_objective = vapply(history, `[[`, numeric(1),
                                     "fixed_working_objective"),
    factor_pve = pve,
    reconstruction_active = pve >= settings$pve_threshold,
    prognostic = abs(EBeta) / pmax(beta_sd, 1e-12) >= settings$prognostic_z_threshold,
    K_eff = sum(pve >= settings$pve_threshold), controls = settings
  )
  out <- list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2, EBeta = EBeta,
              EBeta2 = EBeta2, Tau = Tau, training_spec = data$training_spec,
              diagnostics = diagnostics)
  class(out) <- "multimodal_yfb_fit"
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
