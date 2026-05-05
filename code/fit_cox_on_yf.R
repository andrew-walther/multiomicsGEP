# =============================================================================
# code/fit_cox_on_yf.R
# Purpose: Cluster B CAVI fitting loop — Cox-on-YF reformulation (eta = YFB)
#
# MODEL CHANGE VS. fit_modular.R (Cluster A — eta = LB):
#   Cluster A:  eta_i = l_i · beta
#               L is latent (n×K), learned via dual-source EBNM (genomics + Cox)
#               Training uses EBNM-shrunk E[L]; test uses OLS projection of Y_test
#               onto F — these are DIFFERENT quantities (train/test mismatch).
#
#   Cluster B:  eta_i = (y_i · F) · beta_tilde    [this file]
#               ZF = Y · EF is an n×K matrix of OBSERVED projection scores.
#               Both training and test prediction use Y · EF · beta_tilde —
#               the SAME formula, closing the train/test mismatch.
#
# CONSEQUENCE FOR EACH UPDATE:
#   q(L):    Pure-genomics only. L no longer appears in Cox likelihood.
#            Sources: update_L_surv_YFB.R  (NOT update_L.R)
#   q(F):    Dual-source in principle (F appears in both Y≈LF' and Cox via ZF).
#            alpha_F=0 used here for stability — see DECISIONS.md 2026-04-30.
#            Sources: update_F_surv_YFB.R  (NOT update_F.R)
#   q(beta): Covariate is ZF[,k] = (Y · EF)[,k], NOT EL[,k] as in Cluster A.
#            Sources: update_beta.R          (shared — same interface, different input)
#   q(tau):  Unchanged. Sources: update_tau.R (fully shared)
#   Prediction: Y_test · EF_train · beta_tilde   (L_test never needed)
#
# AUTHOR:      Andrew Walther
# DATE:        May 2026
# DERIVATIONS: derivations/cox_on_YF/
# =============================================================================

# Default: "real" — source() loads the function without side effects.
DATA_MODE <- "real"

real_Y      <- NULL
real_time   <- NULL
real_status <- NULL

# ==============================================================================
# PART 1 — LIBRARIES AND MODULE SOURCES
# ==============================================================================

library(survival)
library(ebnm)

# Cluster B update files (dedicated — do NOT source update_L.R or update_F.R)
source("code/update_L_surv_YFB.R")  # update_L_surv_YFB_k, update_L_surv_YFB_all
source("code/update_F_surv_YFB.R")  # update_F_surv_YFB_k, update_F_surv_YFB_all

# Shared update files (same interface, different inputs vs. Cluster A)
source("code/update_beta.R")        # compute_z_no_k, update_beta_k, update_beta_all
source("code/update_tau.R")         # compute_var_term, update_tau
source("code/compute_elbo.R")       # compute_ebnm_kl, compute_survival_elbo

# compute_R_k is defined in update_L.R (Cluster A); re-source just that function
# by sourcing update_L.R here. The Cluster B update files use only compute_R_k,
# not update_L_k or update_L_all, so there is no crossover.
suppressPackageStartupMessages(
  source("code/update_L.R")           # provides compute_R_k (used by both _surv_YFB files)
)

# ------------------------------------------------------------------------------
#' Calculate Cox Score and Diagonal Hessian (Taylor Expansion)
#'
#' Identical to the helper in fit_modular.R. Duplicated here so fit_cox_on_yf.R
#' is a self-contained entry point (no hidden dependency on fit_modular.R).
#'
#' @param eta    n-vector: current linear predictor
#' @param time   n-vector: observed survival/censoring times
#' @param status n-vector: event indicator (1=event, 0=censored)
#' @return list(u, w, logPL)
# ------------------------------------------------------------------------------
calc_cox_taylor_yf <- function(eta, time, status) {
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
  logPL <- sum(status_s * (eta_s - log(pmax(risk_sum, 1e-300))))
  list(u = u, w = w, logPL = logPL)
}

# ==============================================================================
# PART 2 — fit_cox_on_yf()
# ==============================================================================

#' Fit Cox-on-YF (Cluster B) Supervised Matrix Factorization via CAVI
#'
#' Implements the Cluster B reformulation where the survival linear predictor is
#' eta_i = (y_i · EF) · beta_tilde = ZF_i · beta_tilde. The matrix ZF = Y·EF
#' consists of OBSERVED projection scores, eliminating the train/test formula
#' mismatch present in Cluster A (where training uses EBNM-shrunk E[L] but
#' testing requires a separate OLS projection).
#'
#' Factor-wise Gauss-Seidel CAVI, updating beta_k -> L_k -> F_k for each k
#' before advancing to k+1, then updating Tau once per outer iteration.
#'
#' @param Y        n x p genomics data matrix (should be column-centered)
#' @param time     n-vector of survival / censoring times
#' @param status   n-vector of event indicators (1=event, 0=censored)
#' @param K        Number of latent factors (default 5)
#' @param max_iter Maximum CAVI outer iterations (default 100)
#' @param tol      Convergence threshold on relative full-ELBO change (default 1e-5)
#' @param prior_LF   character: EBNM prior family for L and F (default "point_exponential")
#' @param prior_beta character: EBNM prior family for beta (default "point_normal")
#' @param alpha      numeric in [0, 1]: survival mixing weight for the beta update
#'                 (default 0.5). Note: F update always uses alpha_F=0 (see DECISIONS.md).
#' @param lambda     numeric: lambda multiplier for the beta update (default 1.0)
#' @param init_method character: "svd" (default), "random", or "custom"
#' @param EL_init  Optional n x K matrix: custom initial loadings
#' @param EF_init  Optional p x K matrix: custom initial factors
#' @param N_burnin integer: beta-only burn-in iterations before joint CAVI (default 0).
#'                 0 = off (Cluster A baseline). N_burnin=5 or 10 may help if betas
#'                 collapse after the first joint CAVI iteration.
#' @param cox_warmstart logical: initialize EBeta via Cox regression on ZF before CAVI
#'                 (default FALSE). FALSE = matches Cluster A behavior (EBeta starts at 0).
#'                 TRUE may help if normal prior produces unstable initial betas.
#' @param normalize_AB logical: rescale survival vs. genomics contributions in F update
#'                 (default FALSE). No-op when alpha_F=0 (current default; F update is
#'                 pure-genomics). Forward-compatible when alpha_F > 0 is enabled.
#' @param alpha_schedule NULL or list(warmup_iters, ramp_iters): ramp alpha from 0
#'                 up to `alpha` over warmup+ramp iterations.
#' @param verbose  Logical: print iteration logs? (default TRUE)
#'
#' @return Named list:
#'   $EL, $EL2, $EF, $EF2, $EBeta, $EBeta2, $Tau, $history
fit_cox_on_yf <- function(Y, time, status,
                           K            = 5,
                           max_iter     = 100,
                           tol          = 1e-5,
                           prior_LF     = "point_exponential",
                           prior_beta   = "normal",
                           alpha        = 0.5,
                           lambda       = 1.0,
                           init_method  = "svd",
                           EL_init      = NULL,
                           EF_init      = NULL,
                           N_burnin     = 0,
                           cox_warmstart  = FALSE,
                           normalize_AB   = FALSE,
                           alpha_schedule = NULL,
                           verbose      = TRUE) {

  n <- nrow(Y); p <- ncol(Y)

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1].")
  }

  # ---- Progressive α schedule ----
  use_alpha_schedule <- !is.null(alpha_schedule)
  if (use_alpha_schedule) {
    if (!is.list(alpha_schedule) ||
        !all(c("warmup_iters", "ramp_iters") %in% names(alpha_schedule)) ||
        alpha_schedule$warmup_iters < 0 ||
        alpha_schedule$ramp_iters   < 0) {
      stop("alpha_schedule must be NULL or list(warmup_iters>=0, ramp_iters>=0).")
    }
  }
  alpha_at <- function(iter) {
    if (!use_alpha_schedule) return(alpha)
    wu <- alpha_schedule$warmup_iters
    rp <- alpha_schedule$ramp_iters
    if (iter <= wu) return(0)
    if (iter <= wu + rp) return(alpha * (iter - wu) / max(rp, 1))
    alpha
  }

  if (!is.null(EL_init) && !is.null(EF_init)) init_method <- "custom"

  # --------------------------------------------------------------------------
  # Initialization
  # --------------------------------------------------------------------------
  if (init_method == "svd") {
    svd_init <- svd(Y, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    EL  <- pmax(svd_init$u %*% diag(d_k, K, K), 0)
    EF  <- pmax(svd_init$v %*% diag(d_k, K, K), 0)
  } else if (init_method == "random") {
    y_sd <- sd(Y)
    EL   <- matrix(rnorm(n * K, sd = 0.1 * y_sd), n, K)
    EF   <- matrix(rnorm(p * K, sd = 0.1 * y_sd), p, K)
  } else if (init_method == "custom") {
    if (is.null(EL_init) || is.null(EF_init))
      stop("init_method='custom' requires both EL_init (n x K) and EF_init (p x K).")
    EL <- EL_init
    EF <- EF_init
  } else {
    stop(sprintf("Unknown init_method: '%s'. Use 'svd', 'random', or 'custom'.", init_method))
  }

  EL2 <- EL^2
  EF2 <- EF^2

  # Optionally warm-start beta via Cox regression on ZF = Y·EF.
  # cox_warmstart=FALSE (default) matches Cluster A behavior: EBeta starts at 0.
  # cox_warmstart=TRUE calibrates EBeta to the ZF scale before CAVI begins.
  if (cox_warmstart) {
    ZF_init <- Y %*% EF                              # n × K
    df_cox <- as.data.frame(ZF_init)
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
    if (verbose) {
      cat(sprintf("    [init] Cox warm-start EBeta range: [%.3e, %.3e]\n",
                  min(EBeta), max(EBeta)))
    }
  } else {
    EBeta <- rep(0, K)
  }
  EBeta2 <- EBeta^2

  if (verbose && !cox_warmstart) {
    cat(sprintf("    [init] EBeta initialized to 0 (cox_warmstart=FALSE)\n"))
  }

  # ==========================================================================
  # β-only burn-in: N_burnin iterations with EL and EF held fixed.
  # Under Cluster B, beta's signal path uses ZF = Y·EF (observed), so
  # A_beta = sum(w * ZF_k^2) is non-zero from SVD init regardless of EBeta.
  # Burn-in is a belt-and-suspenders measure to enter joint CAVI with non-zero
  # EBeta regardless of warm-start quality.
  # ==========================================================================
  if (N_burnin > 0) {
    for (b in seq_len(N_burnin)) {
      # Cluster B: ZF = Y · EF (observed projection scores)
      ZF_b     <- Y %*% EF
      eta_b    <- as.vector(ZF_b %*% EBeta)
      taylor_b <- calc_cox_taylor_yf(eta_b, time, status)
      z_b      <- eta_b + taylor_b$u / taylor_b$w
      w_b      <- taylor_b$w
      for (k in seq_len(K)) {
        z_no_k_b <- compute_z_no_k(z_b, ZF_b, EBeta, k)
        # Cluster B: ZF[,k] is observed, so its "second moment" = ZF[,k]^2 (no posterior variance)
        res_b    <- update_beta_k(w_b, z_no_k_b, ZF_b[, k], ZF_b[, k]^2,
                                  prior_family = prior_beta, alpha = alpha)
        EBeta[k]  <- res_b$mean
        EBeta2[k] <- res_b$second
      }
    }
    if (verbose) {
      cat(sprintf("    [burnin] After %d β-only iters, EBeta range: [%.3e, %.3e]\n",
                  N_burnin, min(EBeta), max(EBeta)))
    }
  }

  Tau <- 1.0 / pmax(apply(Y, 2, var), 1e-8)
  y_frob2 <- sum(Y^2)

  history <- list(
    rmse       = numeric(max_iter),
    elbo_proxy = numeric(max_iter),
    elbo_full  = numeric(max_iter),
    delta_L    = numeric(max_iter),
    delta_Beta = numeric(max_iter),
    delta_elbo_rel = rep(NA_real_, max_iter),
    factor_pve = matrix(NA_real_, nrow = max_iter, ncol = K),
    converged  = FALSE,
    n_iter     = max_iter
  )

  if (verbose) {
    cat("=== Cox-on-YF (Cluster B) — Factor-Wise CAVI ===\n")
    cat(sprintf("    n=%d, p=%d, K=%d | max_iter=%d | tol=%.1e\n",
                n, p, K, max_iter, tol))
    cat(sprintf("    prior_LF=%s | prior_beta=%s | alpha=%.2f | alpha_F=0\n\n",
                prior_LF, prior_beta, alpha))
  }

  # ==========================================================================
  # Main CAVI Loop
  # ==========================================================================
  for (iter in 1:max_iter) {

    EL_old    <- EL
    EBeta_old <- EBeta

    alpha_iter <- alpha_at(iter)

    kl_L    <- numeric(K)
    kl_F    <- numeric(K)
    kl_beta <- numeric(K)

    # ------------------------------------------------------------------------
    # STEP 1: Cox Taylor Expansion
    #
    # Cluster B: ZF = Y · EF is the n×K matrix of OBSERVED projection scores.
    # Under eta = ZF · beta_tilde, ZF replaces EL as the survival predictor.
    # ZF is observed (Y is fixed data), so A_beta = sum(w * ZF_k^2) is
    # non-zero from SVD initialization — no chicken-and-egg for beta.
    # ------------------------------------------------------------------------
    ZF     <- Y %*% EF                          # n × K: observed projection scores
    eta    <- as.vector(ZF %*% EBeta)           # eta = ZF * beta_tilde
    taylor <- calc_cox_taylor_yf(eta, time, status)
    z      <- eta + taylor$u / taylor$w    # working response z_i
    w      <- taylor$w                     # W_{ii} diagonal Hessian weights
    YtWY_diag <- as.vector(t(Y^2) %*% w)  # p-vector: diag(Y'diag(w)Y)

    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))

    # ========================================================================
    # STEP 2: Factor-Wise Coordinate Ascent  k = 1, ..., K
    #
    # Update order per factor k:
    #   (a) q(beta_k): survival coefficient — covariate is ZF[,k] (observed)
    #   (b) q(l_k):   pure-genomics only (L not in Cox under Cluster B)
    #   (c) q(f_k):   pure-genomics (alpha_F=0); survival terms receive zero weight
    # ========================================================================
    for (k in 1:K) {

      R_k    <- compute_R_k(Y, EL, EF, k)
      # Cluster B: z_no_k computed w.r.t. ZF (not EL as in Cluster A)
      z_no_k <- compute_z_no_k(z, ZF, EBeta, k)

      # ---- (a) Update q(beta_k): Survival Coefficient ----
      # Covariate is ZF[,k] = (Y·EF)[,k] — observed projection (not latent).
      # Pass ZF[,k]^2 as EL2_k (the second moment argument): ZF is observed, so
      # its "posterior second moment" = squared value (no variance term).
      # A_beta = sum(w * ZF_k^2) is non-zero from SVD init regardless of EBeta,
      # breaking the chicken-and-egg that plagued Cluster A's L update.
      res_beta    <- update_beta_k(w, z_no_k, ZF[, k], ZF[, k]^2,
                                   prior_family = prior_beta, alpha = alpha_iter)
      EBeta[k]    <- res_beta$mean
      EBeta2[k]   <- res_beta$second
      kl_beta[k]  <- compute_ebnm_kl(res_beta$ebnm_result$log_likelihood,
                                      res_beta$A, res_beta$x,
                                      res_beta$mean, res_beta$second)

      # ---- (b) Update q(l_k): Patient Loadings — PURE GENOMICS ----
      # L does not appear in the Cox likelihood under eta = ZF * beta_tilde.
      # Uses update_L_surv_YFB_k() (NOT update_L_k from Cluster A).
      res_L   <- update_L_surv_YFB_k(Tau, EF[, k], EF2[, k], R_k,
                                      prior_family = prior_LF)
      EL[, k]  <- res_L$mean
      EL2[, k] <- res_L$second
      kl_L[k]  <- compute_ebnm_kl(res_L$ebnm_result$log_likelihood,
                                   res_L$A, res_L$x, res_L$mean, res_L$second)

      # ---- (c) Update q(f_k): Biological Factors — alpha_F=0 ----
      # R_k depends on EL[,k] (updated above); recompute before F update.
      # update_F_surv_YFB_k() defaults to alpha=0 (pure genomics), preventing
      # the positive-feedback instability from the dual-source F update.
      # YtWz_no_k passed for completeness but receives zero weight at alpha=0.
      R_k        <- compute_R_k(Y, EL, EF, k)
      YtWz_no_k  <- as.vector(t(Y) %*% (w * z_no_k))  # p-vector
      res_F <- update_F_surv_YFB_k(Tau, EL[, k], EL2[, k], R_k,
                                    EBeta_k      = EBeta[k],
                                    EBeta2_k     = EBeta2[k],
                                    YtWY_diag    = YtWY_diag,
                                    YtWz_no_k    = YtWz_no_k,
                                    prior_family = prior_LF,
                                    alpha        = 0)   # alpha_F=0: pure genomics
      EF[, k]  <- res_F$mean
      EF2[, k] <- res_F$second
      kl_F[k]  <- compute_ebnm_kl(res_F$ebnm_result$log_likelihood,
                                   res_F$A, res_F$x, res_F$mean, res_F$second)

    }  # end k-loop

    # ========================================================================
    # STEP 3: Noise Precision Update
    # ========================================================================
    res_tau              <- update_tau(Y, EL, EL2, EF, EF2)
    Tau                  <- res_tau$Tau
    history$elbo_proxy[iter] <- res_tau$elbo_proxy
    factor_pve_iter <- vapply(seq_len(K), function(k) {
      sum(EL[, k]^2) * sum(EF[, k]^2) / y_frob2
    }, numeric(1))
    history$factor_pve[iter, ] <- factor_pve_iter

    # Full ELBO under Cluster B: ZF = Y·EF is the survival predictor.
    # Pass ZF^2 as EL2 (element-wise squared): ZF is observed, so its
    # posterior second moment equals its squared value (no variance term).
    surv_elbo               <- compute_survival_elbo(taylor$logPL, w,
                                                     ZF, ZF^2, EBeta, EBeta2)
    history$elbo_full[iter] <- (1 - alpha) * res_tau$elbo_proxy +
                               alpha * surv_elbo +
                               sum(kl_L) + sum(kl_F) + sum(kl_beta)

    # ========================================================================
    # STEP 4: Convergence Check
    # ========================================================================
    delta_L    <- mean(abs(EL - EL_old))
    delta_Beta <- mean(abs(EBeta - EBeta_old))
    history$delta_L[iter]    <- delta_L
    history$delta_Beta[iter] <- delta_Beta

    if (iter > 1) {
      elbo_old <- history$elbo_full[iter - 1]
      elbo_new <- history$elbo_full[iter]
      if (is.finite(elbo_new) && is.finite(elbo_old) && elbo_old != 0) {
        history$delta_elbo_rel[iter] <- abs((elbo_new - elbo_old) / abs(elbo_old))
      }
    }

    if (verbose && iter %% 10 == 0) {
      delta_elbo_str <- if (is.finite(history$delta_elbo_rel[iter])) {
        sprintf("%.2e", history$delta_elbo_rel[iter])
      } else "NA"
      cat(sprintf("  iter %3d | RMSE: %.4f | ELBO: %+.1f | dELBO: %s | dB: %.2e | beta: [%s]\n",
                  iter, history$rmse[iter], history$elbo_full[iter],
                  delta_elbo_str, delta_Beta,
                  paste(sprintf("%+.2f", EBeta), collapse = ", ")))
    }

    if (iter > 5 &&
        is.finite(history$delta_elbo_rel[iter]) &&
        history$delta_elbo_rel[iter] < tol) {
      if (verbose) {
        cat(sprintf("\n  Converged at iteration %d  (relative ELBO change = %.2e)\n",
                    iter, history$delta_elbo_rel[iter]))
      }
      history$converged  <- TRUE
      history$n_iter     <- iter
      history$rmse       <- history$rmse[1:iter]
      history$elbo_proxy <- history$elbo_proxy[1:iter]
      history$elbo_full  <- history$elbo_full[1:iter]
      history$delta_L    <- history$delta_L[1:iter]
      history$delta_Beta <- history$delta_Beta[1:iter]
      history$delta_elbo_rel <- history$delta_elbo_rel[1:iter]
      history$factor_pve <- history$factor_pve[1:iter, , drop = FALSE]
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
# PART 3 — DATA_MODE Runner Block (for standalone Rscript use)
# ==============================================================================

if (DATA_MODE == "real") {
  stopifnot(!is.null(real_Y), !is.null(real_time), !is.null(real_status))
  res <- fit_cox_on_yf(real_Y, real_time, real_status,
                       K = 5, max_iter = 100, tol = 1e-3, verbose = TRUE)
  cat(sprintf("  Converged: %s | Iterations: %d\n",
              res$history$converged, res$history$n_iter))
}
