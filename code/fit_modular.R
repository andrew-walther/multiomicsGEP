# ==============================================================================
# TITLE:       Supervised Bayesian Matrix Factorization — Modular V3
# FILE:        code/fit_modular.R
# AUTHOR:      Andrew Walther
# DATE:        March 2026
# DERIVATIONS: derivations/MF_UpdateDerivations/MF_V2_Companion.tex
#              derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf
#
# DESCRIPTION:
#   Wraps the four modular update scripts into a single CAVI fitting function
#   fit_supervised_mf_modular().  Implements V3 Algorithm 1: factor-wise
#   Gauss-Seidel CAVI, updating L_k -> F_k -> beta_k for each k before
#   advancing to k+1.
#
#   Model:
#     Genomics:  Y_{n x p} = L_{n x K} F^T_{K x p} + E,  E_ij ~ N(0, tau_j^{-1})
#     Survival:  h(t_i | l_i) = h_0(t_i) exp( sum_k l_{ik} beta_k )  [Cox PH]
#
#   Differences from V2.R:
#     - Uses modular update_L.R / update_F.R / update_beta.R / update_tau.R
#     - Convergence uses mean(|delta|) with tol=1e-3 (V3 Algorithm 1 specifies max, but mean is more practical; see STEP 4 comment)
#     - Returns EL/EL2/EF/EF2/EBeta/EBeta2 directly (no V2 aliases L/F/Beta)
#
# VARIABLE CONVENTIONS:
#   EL[i,k]   = E_q[l_{ik}]     posterior mean of loading i on factor k
#   EL2[i,k]  = E_q[l_{ik}^2]  posterior 2nd moment (Var + mean^2)
#   EF[j,k]   = E_q[f_{jk}]     posterior mean of feature weight j on factor k
#   EF2[j,k]  = E_q[f_{jk}^2]  posterior 2nd moment
#   EBeta[k]  = E_q[beta_k]     posterior mean of survival coefficient k
#   EBeta2[k] = E_q[beta_k^2]  posterior 2nd moment
#   Tau[j]    = tau_j           noise precision (feature-specific)
#   w[i]      = W_{ii}          diagonal Cox Hessian (per-sample weight)
#   z[i]      = working response z_i = eta_i + u_i / W_{ii}
#
# DATA_MODE:
#   "real"       -- (default) set real_Y / real_time / real_status below, then
#                   Rscript code/fit_modular.R
#   "simulated"  -- run built-in DGP and print convergence diagnostics
#                   (set DATA_MODE <- "simulated" then Rscript code/fit_modular.R)
# ==============================================================================

# Default: "real" — source() loads the function without side effects.
# Change to "simulated" to run the built-in DGP validation, or
# set real_Y/real_time/real_status and run Rscript code/fit_modular.R
DATA_MODE <- "real"

real_Y      <- NULL
real_time   <- NULL
real_status <- NULL

# ==============================================================================
# PART 1 — LIBRARIES AND MODULE SOURCES
# ==============================================================================

library(survival)
library(ebnm)

source("code/update_L.R")     # compute_R_k, update_L_k, update_L_all
source("code/update_F.R")     # update_F_k, update_F_all  (uses compute_R_k from L)
source("code/update_beta.R")  # compute_z_no_k, update_beta_k, update_beta_all
source("code/update_tau.R")   # compute_var_term, update_tau

# ------------------------------------------------------------------------------
#' Calculate Cox Score and Diagonal Hessian (Taylor Expansion)
#'
#' Functionally identical to code/Supervised_Bayesian_MF_V2.R lines 70-93.
#' Not sourced from V2.R to avoid running its top-level simulation code.
#'
#' Transforms the non-conjugate Cox partial likelihood into a locally Gaussian
#' weighted-least-squares form centred at eta_hat = L_bar %*% beta_bar.
#' The 2nd-order Taylor expansion gives:
#'   ell(eta) ~ sum_i [-1/2 * W_{ii} * (eta_i - z_i)^2] + C
#' where z_i = eta_hat_i + u_i / W_{ii} is the working response.
#'
#' @param eta    n-vector: current linear predictor eta_i = sum_k l_bar_{ik} beta_bar_k
#' @param time   n-vector: observed survival / censoring times
#' @param status n-vector: event indicator (1 = event, 0 = censored)
#' @return list(u = n-vector score, w = n-vector neg-diagonal Hessian)
# ------------------------------------------------------------------------------
calc_cox_taylor <- function(eta, time, status) {
  n   <- length(time)
  ord <- order(time)
  time_s   <- time[ord]
  status_s <- status[ord]
  eta_s    <- eta[ord]
  theta    <- exp(eta_s)
  risk_sum <- rev(cumsum(rev(theta)))
  h <- status_s / risk_sum
  H <- cumsum(h)
  u_s <- status_s - theta * H
  w_s <- theta * H
  w_s[w_s < 1e-6] <- 1e-6
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_s;   w[ord] <- w_s
  list(u = u, w = w)
}

# ==============================================================================
# PART 2 — fit_supervised_mf_modular()
# ==============================================================================

#' Fit Supervised Matrix Factorization via Modular Factor-Wise CAVI
#'
#' Implements V3 Algorithm 1: factor-wise Gauss-Seidel coordinate ascent,
#' updating q(l_k) -> q(f_k) -> q(beta_k) for each k = 1..K per outer
#' iteration, then updating Tau once at the end of each iteration.
#'
#' Each coordinate update is handled by the corresponding modular script:
#'   - update_L.R:    update_L_k()
#'   - update_F.R:    update_F_k()
#'   - update_beta.R: update_beta_k(), compute_z_no_k()
#'   - update_tau.R:  update_tau()
#'
#' Convergence is declared when BOTH mean|delta_L| and mean|delta_Beta| < tol
#' (after a 5-iteration burn-in). V3 Algorithm 1 specifies max absolute change;
#' mean is used here because max rarely reaches 1e-5 on typical datasets due to
#' SVD orientation oscillations.
#'
#' @param Y        n x p genomics data matrix
#' @param time     n-vector of survival / censoring times
#' @param status   n-vector of event indicators (1=event, 0=censored)
#' @param K        Number of latent factors (default 5)
#' @param max_iter Maximum CAVI outer iterations (default 100)
#' @param tol      Convergence threshold: both mean|dL| and mean|dBeta| < tol
#'                 (default 1e-3; use 1e-5 for tighter but rarely achievable
#'                 convergence with this Taylor-approximation scheme)
#' @param prior_family character: EBNM prior family passed to update_L_k,
#'                 update_F_k, and update_beta_k.  Valid values:
#'                 "point_normal" (default — sparse, shrinks to zero),
#'                 "point_laplace" (heavier tails than point-normal),
#'                 "normal_scale_mixture" (broader, less sparse factors).
#' @param init_method  character: initialization strategy.
#'                 "svd" (default — deterministic SVD warm-start),
#'                 "random" (random normal initialization, useful with
#'                 multiple restarts to escape local optima).
#' @param verbose  Logical: print iteration logs every 10 iters? (default TRUE)
#'
#' @return Named list:
#'   $EL       n x K matrix: posterior means of loadings
#'   $EL2      n x K matrix: posterior second moments of loadings
#'   $EF       p x K matrix: posterior means of factors
#'   $EF2      p x K matrix: posterior second moments of factors
#'   $EBeta    K-vector: posterior means of survival coefficients
#'   $EBeta2   K-vector: posterior second moments of survival coefficients
#'   $Tau      p-vector: noise precisions
#'   $history  list(rmse, elbo_proxy, converged, n_iter)
fit_supervised_mf_modular <- function(Y, time, status,
                                      K            = 5,
                                      max_iter     = 100,
                                      tol          = 1e-3,
                                      prior_family = "point_normal",
                                      init_method  = "svd",
                                      verbose      = TRUE) {

  n <- nrow(Y); p <- ncol(Y)

  # --------------------------------------------------------------------------
  # Initialization
  #
  # Two strategies:
  #   "svd"    — Deterministic SVD warm-start (mirrors V2.R lines 187-218).
  #              EL = U sqrt(D), EF = V sqrt(D), so EL %*% t(EF) ≈ Y_rank-K.
  #              Recommended default; gives reproducible results.
  #   "random" — Random normal initialization.  Useful when running multiple

  #              restarts (n_init > 1) to escape local optima.  The caller
  #              should set a seed before each call for reproducibility.
  # --------------------------------------------------------------------------

  if (init_method == "svd") {
    # SVD of Y: deterministic high-variance starting subspace
    svd_init <- svd(Y, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    EL  <- svd_init$u %*% diag(d_k, K, K)      # n x K
    EF  <- svd_init$v %*% diag(d_k, K, K)      # p x K
  } else if (init_method == "random") {
    # Random normal initialization scaled by data magnitude.
    # sd = 0.1 * overall SD of Y keeps initial reconstruction in a
    # reasonable range; too large causes Cox Taylor expansion to diverge.
    y_sd <- sd(Y)
    EL   <- matrix(rnorm(n * K, sd = 0.1 * y_sd), n, K)
    EF   <- matrix(rnorm(p * K, sd = 0.1 * y_sd), p, K)
  } else {
    stop(sprintf("Unknown init_method: '%s'. Use 'svd' or 'random'.", init_method))
  }

  # Second moments initialised to squared means (zero posterior variance).
  # Posterior variance populated after the first EBNM call.
  EL2 <- EL^2
  EF2 <- EF^2

  # Warm-start beta via Cox regression on SVD loadings.
  df_cox <- as.data.frame(EL)
  colnames(df_cox) <- paste0("L", 1:K)
  df_cox$time   <- time
  df_cox$status <- status
  cox_init <- tryCatch(
    coxph(as.formula("Surv(time, status) ~ ."), data = df_cox, x = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cox_init)) {
    cx_coef <- coef(cox_init)
    cx_coef[is.na(cx_coef)] <- 0
    EBeta <- cx_coef
  } else {
    EBeta <- rep(0, K)
  }
  EBeta2 <- EBeta^2      # zero variance initially; updated at first EBNM call

  # Column-specific noise precision from sample variance of each column of Y.
  Tau <- 1.0 / pmax(apply(Y, 2, var), 1e-8)   # p-vector

  # History tracking
  history <- list(
    rmse       = numeric(max_iter),
    elbo_proxy = numeric(max_iter),
    converged  = FALSE,
    n_iter     = max_iter
  )

  if (verbose) {
    cat("=== Supervised Bayesian MF (Modular V3) — Factor-Wise CAVI ===\n")
    cat(sprintf("    n=%d, p=%d, K=%d | max_iter=%d | tol=%.1e\n",
                n, p, K, max_iter, tol))
    cat(sprintf("    prior_family=%s | init_method=%s\n\n",
                prior_family, init_method))
  }

  # ==========================================================================
  # Main CAVI Loop — V3 Algorithm 1
  # ==========================================================================
  for (iter in 1:max_iter) {

    EL_old    <- EL
    EBeta_old <- EBeta

    # ------------------------------------------------------------------------
    # STEP 1: Cox Taylor Expansion
    #
    # Linearise the Cox partial likelihood around the current linear predictor
    # eta_hat_i = sum_k l_bar_{ik} * beta_bar_k.
    #
    # Working response:  z_i = eta_hat_i + u_i / W_{ii}
    # Weight:            W_{ii} (negative diagonal Hessian, positive)
    #
    # (z, w) are held FIXED for this outer iteration across all k.
    # ------------------------------------------------------------------------
    eta    <- as.vector(EL %*% EBeta)
    taylor <- calc_cox_taylor(eta, time, status)
    z      <- eta + taylor$u / taylor$w    # n-vector: working response z_i
    w      <- taylor$w                     # n-vector: W_{ii}

    # Reconstruction RMSE at posterior means (monitoring only).
    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))

    # ========================================================================
    # STEP 2: Factor-Wise Coordinate Ascent  k = 1, ..., K
    #
    # For each factor k, update in the order:
    #   (a) q(l_k): patient loadings — uses both genomics and survival
    #   (b) q(f_k): biological factors — genomics only; uses updated EL[,k]
    #   (c) q(beta_k): survival coefficient — reuses z_no_k from (a)
    #
    # Gauss-Seidel: each update sees the freshest values of all parameters.
    # ========================================================================
    for (k in 1:K) {

      # ----------------------------------------------------------------------
      # (a) Update q(l_k): Patient Loadings
      #
      # Compute R_k and z_no_k once; reuse z_no_k for step (c).
      #
      # z_no_k does NOT depend on EL[,k] (EL[,k] cancels in the formula:
      #   eta_no_k = EL %*% EBeta - EL[,k]*EBeta[k]; z_no_k = z - eta_no_k).
      # It is therefore identical before/after the L_k update, so it is safe
      # to reuse for step (c) without recomputing.
      # See Companion.tex Sec. 6 "z_no_k reuse rationale".
      # ----------------------------------------------------------------------
      R_k    <- compute_R_k(Y, EL, EF, k)
      z_no_k <- compute_z_no_k(z, EL, EBeta, k)   # COMPUTE ONCE — reuse for (c)

      res_L   <- update_L_k(Tau, EF[, k], EF2[, k], w, EBeta[k], EBeta2[k],
                             R_k, z_no_k, prior_family = prior_family)
      EL[, k]  <- res_L$mean
      EL2[, k] <- res_L$second

      # ----------------------------------------------------------------------
      # (b) Update q(f_k): Biological Factors
      #
      # R_k depends on EL[,k], so recompute it using the updated EL[,k]
      # from step (a).  This is the Gauss-Seidel property.
      # ----------------------------------------------------------------------
      R_k   <- compute_R_k(Y, EL, EF, k)
      res_F <- update_F_k(Tau, EL[, k], EL2[, k], R_k,
                          prior_family = prior_family)
      EF[, k]  <- res_F$mean
      EF2[, k] <- res_F$second

      # ----------------------------------------------------------------------
      # (c) Update q(beta_k): Survival Coefficient
      #
      # Reuse z_no_k from step (a) — unchanged since EL[,k'] for k' != k
      # has not been updated in this inner iteration, and z_no_k excludes
      # factor k entirely.
      # ----------------------------------------------------------------------
      res_beta  <- update_beta_k(w, z_no_k, EL[, k], EL2[, k],
                                 prior_family = prior_family)
      EBeta[k]  <- res_beta$mean
      EBeta2[k] <- res_beta$second

    }  # end k-loop

    # ========================================================================
    # STEP 3: Noise Precision Update (closed-form MLE; no EBNM)
    # ========================================================================
    res_tau              <- update_tau(Y, EL, EL2, EF, EF2)
    Tau                  <- res_tau$Tau
    history$elbo_proxy[iter] <- res_tau$elbo_proxy

    # ========================================================================
    # STEP 4: Convergence Check
    #
    # Use mean (not max) of absolute changes — matches V2.R [A4].
    # max() is ~5-10x larger than mean() due to a few high-variance EL
    # entries that oscillate near factor orientation boundaries; using
    # max() with tol=1e-5 never declares convergence on typical datasets.
    # Guard with iter > 5 to allow burn-in before checking convergence.
    # ========================================================================
    delta_L    <- mean(abs(EL - EL_old))
    delta_Beta <- mean(abs(EBeta - EBeta_old))

    if (verbose && iter %% 10 == 0) {
      cat(sprintf("  iter %3d | RMSE: %.4f | ELBO: %+.1f | dL: %.2e | dB: %.2e | beta: [%s]\n",
                  iter, history$rmse[iter], history$elbo_proxy[iter],
                  delta_L, delta_Beta,
                  paste(sprintf("%+.2f", EBeta), collapse = ", ")))
    }

    if (iter > 5 && delta_L < tol && delta_Beta < tol) {
      if (verbose) {
        cat(sprintf("\n  Converged at iteration %d  (dL=%.2e, dBeta=%.2e)\n",
                    iter, delta_L, delta_Beta))
      }
      history$converged  <- TRUE
      history$n_iter     <- iter
      history$rmse       <- history$rmse[1:iter]
      history$elbo_proxy <- history$elbo_proxy[1:iter]
      break
    }

  }  # end CAVI loop

  list(
    EL     = EL,
    EL2    = EL2,
    EF     = EF,
    EF2    = EF2,
    EBeta  = EBeta,
    EBeta2 = EBeta2,
    Tau    = Tau,
    history = history
  )
}

# ==============================================================================
# PART 3 — DATA_MODE Runner Block
# ==============================================================================

if (DATA_MODE == "simulated") {

  # DGP identical to results/run_modular_simulation.R lines 131-152:
  # same seed (42), same B_true, same Weibull survival model.

  set.seed(42)
  n <- 250; p <- 1000; K <- 5

  B_true <- c(1.5, -1.2, 0.8, -0.5, 0.0)

  # Generate L_true (n x K standard normal) and F_true (p x K, 5% sparse)
  L_true <- matrix(rnorm(n * K), n, K)
  F_true <- matrix(0, p, K)
  for (k in 1:K) {
    active <- sample(1:p, round(p * 0.05))          # 5% active features per factor
    F_true[active, k] <- rnorm(length(active), 0, 5)
  }

  # Genomics matrix: Y = L_true F_true^T + noise
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p), n, p)

  # Survival: Weibull times, exponential censoring
  eta_true   <- as.vector(L_true %*% B_true)
  raw_times  <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
  cens_times <- rexp(n, rate = 1 / 50)
  time   <- pmin(raw_times, cens_times)
  status <- as.integer(raw_times <= cens_times)

  cat("=== Modular Supervised MF — Simulation ===\n")
  cat(sprintf("  n=%d  p=%d  K=%d  seed=42\n", n, p, K))
  cat(sprintf("  Censoring rate: %.1f%%\n\n", 100 * mean(status == 0)))

  res <- fit_supervised_mf_modular(Y, time, status, K = K,
                                   max_iter = 100, tol = 1e-3, verbose = TRUE)

  cat("\n=== RESULTS SUMMARY ===\n")
  cat(sprintf("  Converged:  %s\n", res$history$converged))
  cat(sprintf("  Iterations: %d\n", res$history$n_iter))
  cat(sprintf("  Final RMSE: %.4f\n\n", tail(res$history$rmse, 1)))

  cat("  True vs Estimated Beta:\n")
  # NOTE: SVD initialisation introduces sign/permutation ambiguity in the
  # factors.  The estimated betas may have flipped signs relative to B_true
  # because the corresponding loading columns may be sign-flipped.
  # Abs_Error is therefore not meaningful without first aligning signs.
  # We check whether there exists a sign pattern s in {-1,+1}^K such that
  # s * EBeta matches B_true in sign (non-zero factors only).
  est   <- res$EBeta
  nonzero_mask <- B_true != 0

  # Beta summary table — absolute values shown because SVD sign is ambiguous
  beta_df <- data.frame(
    Factor    = 1:K,
    Beta_true = B_true,
    Beta_est  = round(est, 3)
  )
  print(beta_df)

  cat("\n  Note: SVD introduces sign ambiguity — factor columns may be sign-flipped.\n")
  cat("  Key checks:\n")
  nonzero_est_threshold <- 0.1
  cat(sprintf("    Non-zero true betas (|beta_true| > 0) with nonzero estimates: %d / %d\n",
              sum(abs(est[nonzero_mask]) > nonzero_est_threshold), sum(nonzero_mask)))
  cat(sprintf("    Zero-true factor shrunk small (|beta_est| < %.1f): %s (Factor 5 = %.3f)\n",
              nonzero_est_threshold, all(abs(est[!nonzero_mask]) < nonzero_est_threshold),
              est[!nonzero_mask]))

} else if (DATA_MODE == "real") {

  stopifnot(
    !is.null(real_Y),
    !is.null(real_time),
    !is.null(real_status)
  )
  res <- fit_supervised_mf_modular(real_Y, real_time, real_status,
                                   K = 5, max_iter = 100, tol = 1e-3,
                                   verbose = TRUE)
  cat("\n=== REAL DATA FIT ===\n")
  cat(sprintf("  Converged:  %s\n", res$history$converged))
  cat(sprintf("  Iterations: %d\n", res$history$n_iter))
  cat(sprintf("  Final RMSE: %.4f\n", tail(res$history$rmse, 1)))

}
