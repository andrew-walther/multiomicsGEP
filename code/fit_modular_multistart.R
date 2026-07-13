# ============================================================
# Script:  fit_modular_multistart.R
# Purpose: Multi-initialization wrapper for fit_supervised_mf_modular().
#          Runs N_init restarts (restart 1 = deterministic SVD; restarts 2..N
#          = random normal inits with reproducible seeds) and selects the fit
#          with the highest final training ELBO.
#          Selection criterion: ELBO (not held-out C-index) — selecting on
#          C-index would leak the validation cohort into model selection.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Dependencies: code/fit_modular.R (must be sourced before this file)
# ============================================================

#' Multi-Initialization Wrapper for Supervised Matrix Factorization
#'
#' Runs fit_supervised_mf_modular() N_init times with different initializations,
#' selects the best fit by final training ELBO, and returns the full restart
#' landscape for diagnostics.
#'
#' Restart 1 always uses deterministic SVD initialization (guarantees the
#' single-init baseline is a member of the candidate set). Restarts 2..N_init
#' use random normal initialization with seed = init_seed_base + i, giving
#' reproducible restarts when init_seed_base is fixed.
#'
#' Selection criterion: highest final training ELBO. ELBO is the variational
#' objective; selecting on it cleanly separates model fitting from validation
#' and avoids leaking the external cohort into model selection. Held-out
#' C-index can be computed outside this function on the returned best fit.
#'
#' @param Y           n x p genomics data matrix
#' @param time        n-vector of survival / censoring times
#' @param status      n-vector of event indicators (1=event, 0=censored)
#' @param ...         Additional arguments passed to fit_supervised_mf_modular()
#'                    (K, max_iter, tol, prior_LF, prior_beta, alpha,
#'                    N_burnin, sign_correction, verbose, etc.)
#'                    init_method is overridden internally — do not pass it.
#' @param n_init      integer: number of random restarts (default 30).
#'                    Restart 1 is always SVD; restarts 2..n_init are random.
#' @param init_seed_base  integer or NULL: seed offset for random restarts.
#'                    Restart i uses set.seed(init_seed_base + i). If NULL,
#'                    seeds are drawn from the current RNG state (not
#'                    reproducible across sessions unless outer seed is set).
#' @param beta_threshold  scalar: |beta| threshold for counting K_eff in the
#'                    restart diagnostics (default 0.001).
#'
#' @return Named list:
#'   $best         The fit_supervised_mf_modular() output with highest ELBO.
#'   $best_idx     Integer index of the winning restart (1 = SVD init).
#'   $restarts     data.frame with one row per restart:
#'                   init_id, init_method, seed, final_elbo, k_eff,
#'                   beta_max, n_iter, converged, train_cindex
#'
#' @examples
#' \dontrun{
#' ms <- fit_supervised_mf_modular_multistart(Y, time, status,
#'         K = 10, prior_beta = "point_normal", n_init = 30)
#' best_fit  <- ms$best
#' landscape <- ms$restarts   # ELBO and C-index for every restart
#' }
fit_supervised_mf_modular_multistart <- function(Y, time, status,
                                                  ...,
                                                  n_init           = 30,
                                                  init_seed_base   = 42,
                                                  beta_threshold   = 0.001) {
  if (!is.numeric(n_init) || n_init < 1 || n_init != round(n_init))
    stop("n_init must be a positive integer.")

  n_obs <- nrow(Y)

  fits    <- vector("list", n_init)
  seeds   <- c(NA_integer_, if (n_init > 1) as.integer(init_seed_base) + 2:n_init)
  methods <- c("svd",       if (n_init > 1) rep("random", n_init - 1))

  for (i in seq_len(n_init)) {
    if (methods[i] == "svd") {
      fits[[i]] <- fit_supervised_mf_modular(
        Y, time, status, ..., init_method = "svd", verbose = FALSE
      )
    } else {
      set.seed(seeds[i])
      fits[[i]] <- fit_supervised_mf_modular(
        Y, time, status, ..., init_method = "random", verbose = FALSE
      )
    }
  }

  # Select best fit by final training ELBO
  elbos    <- sapply(fits, function(f) {
    elbo_vec <- f$history$elbo_full
    if (length(elbo_vec) == 0 || all(is.na(elbo_vec))) -Inf else tail(elbo_vec, 1)
  })
  best_idx <- which.max(elbos)

  # Build restart diagnostics data.frame
  restarts <- do.call(rbind, lapply(seq_len(n_init), function(i) {
    f      <- fits[[i]]
    ebeta  <- f$EBeta
    el     <- f$EL
    eta    <- el %*% ebeta
    c_train <- tryCatch(
      as.numeric(survival::concordance(survival::Surv(time, status) ~ eta)$concordance),
      error = function(e) NA_real_
    )
    data.frame(
      init_id      = i,
      init_method  = methods[i],
      seed         = seeds[i],
      final_elbo   = round(elbos[i], 2),
      k_eff        = sum(abs(ebeta) > beta_threshold),
      beta_max     = round(max(abs(ebeta)), 4),
      n_iter       = f$history$n_iter,
      converged    = f$history$converged,
      train_cindex = round(c_train, 4),
      stringsAsFactors = FALSE
    )
  }))

  list(
    best     = fits[[best_idx]],
    best_idx = best_idx,
    restarts = restarts
  )
}

# ============================================================
# Cluster B (YFB) variant -- mirrors fit_supervised_mf_modular_multistart()
# above for fit_cox_on_yf(). Kept in this same file (not a dedicated
# fit_cox_on_yf_multistart.R) since the two wrappers share the identical
# restart/ELBO-selection structure and differ only in which fit function is
# called and how train_cindex is computed from the result.
# ============================================================

#' Multi-Initialization Wrapper for fit_cox_on_yf() (Cluster B / YFB)
#'
#' Runs fit_cox_on_yf() N_init times with different initializations, selects
#' the fit with the highest final training ELBO, and returns the full restart
#' landscape for diagnostics -- the YFB counterpart to
#' fit_supervised_mf_modular_multistart() above.
#'
#' Selection criterion: highest final training ELBO (same rationale as the
#' LB version -- avoids leaking the external cohort into model selection).
#'
#' @param Y           n x p genomics data matrix
#' @param time        n-vector of survival / censoring times
#' @param status      n-vector of event indicators (1=event, 0=censored)
#' @param ...         Additional arguments passed to fit_cox_on_yf() (K,
#'                    max_iter, tol, prior_beta, alpha, sign_correction,
#'                    verbose, etc.). init_method is overridden internally --
#'                    do not pass it.
#' @param n_init      integer: number of random restarts (default 30).
#'                    Restart 1 is always SVD; restarts 2..n_init are random.
#' @param init_seed_base  integer or NULL: seed offset for random restarts.
#'                    Restart i uses set.seed(init_seed_base + i).
#' @param beta_threshold  scalar: |beta| threshold for counting K_eff in the
#'                    restart diagnostics (default 0.001).
#'
#' @return Named list:
#'   $best         The fit_cox_on_yf() output with highest ELBO.
#'   $best_idx     Integer index of the winning restart (1 = SVD init).
#'   $restarts     data.frame with one row per restart:
#'                   init_id, init_method, seed, final_elbo, k_eff,
#'                   beta_max, n_iter, converged, train_cindex
#'                 train_cindex uses Cluster B's predictor ZF·beta, where
#'                 ZF = Y·EF (normalized by EF_norms) -- NOT EL·beta.
#'
#' @examples
#' \dontrun{
#' ms <- fit_cox_on_yf_multistart(Y, time, status,
#'         K = 5, prior_beta = "normal", n_init = 15)
#' best_fit  <- ms$best
#' landscape <- ms$restarts
#' }
fit_cox_on_yf_multistart <- function(Y, time, status,
                                      ...,
                                      n_init           = 30,
                                      init_seed_base   = 42,
                                      beta_threshold   = 0.001) {
  if (!is.numeric(n_init) || n_init < 1 || n_init != round(n_init))
    stop("n_init must be a positive integer.")

  fits    <- vector("list", n_init)
  seeds   <- c(NA_integer_, if (n_init > 1) as.integer(init_seed_base) + 2:n_init)
  methods <- c("svd",       if (n_init > 1) rep("random", n_init - 1))

  for (i in seq_len(n_init)) {
    if (methods[i] == "svd") {
      fits[[i]] <- fit_cox_on_yf(
        Y, time, status, ..., init_method = "svd", verbose = FALSE
      )
    } else {
      set.seed(seeds[i])
      fits[[i]] <- fit_cox_on_yf(
        Y, time, status, ..., init_method = "random", verbose = FALSE
      )
    }
  }

  # Select best fit by final training ELBO
  elbos    <- sapply(fits, function(f) {
    elbo_vec <- f$history$elbo_full
    if (length(elbo_vec) == 0 || all(is.na(elbo_vec))) -Inf else tail(elbo_vec, 1)
  })
  best_idx <- which.max(elbos)

  # Build restart diagnostics data.frame. Cluster B's predictor is
  # ZF %*% beta, where ZF = Y %*% EF normalized by EF_norms (see the
  # Phase C sign-correction block in fit_cox_on_yf.R for the same formula).
  restarts <- do.call(rbind, lapply(seq_len(n_init), function(i) {
    f      <- fits[[i]]
    ebeta  <- f$EBeta
    ZF     <- Y %*% sweep(f$EF, 2, f$EF_norms, "/")
    eta    <- as.vector(ZF %*% ebeta)
    c_train <- tryCatch(
      as.numeric(survival::concordance(survival::Surv(time, status) ~ eta)$concordance),
      error = function(e) NA_real_
    )
    data.frame(
      init_id      = i,
      init_method  = methods[i],
      seed         = seeds[i],
      final_elbo   = round(elbos[i], 2),
      k_eff        = sum(abs(ebeta) > beta_threshold),
      beta_max     = round(max(abs(ebeta)), 4),
      n_iter       = f$history$n_iter,
      converged    = f$history$converged,
      train_cindex = round(c_train, 4),
      stringsAsFactors = FALSE
    )
  }))

  list(
    best     = fits[[best_idx]],
    best_idx = best_idx,
    restarts = restarts
  )
}
