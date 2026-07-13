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
source("code/compute_elbo.R")       # compute_ebnm_kl, compute_survival_elbo, compute_normal_kl
source("code/update_F_cohort.R")    # update_F_cohort_col, update_F_cohort_all
source("code/deflation_init.R")     # deflation_svd_init -- init_method="deflation"

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
#' @param norm_convention character: Phase 1a objective normalization convention
#'                 (see DECISIONS.md and fit_modular.R). Under YFB, L and F are
#'                 pure-genomics with no survival competition at the default
#'                 alpha_F=0 (unlike LB's L), so the genomics divisor does NOT
#'                 reach update_L_surv_YFB_k/update_F_surv_YFB_k -- there is no
#'                 imbalance to fix there, and rescaling a pure-genomics precision
#'                 would only add a gratuitous change with no benefit. Only the
#'                 ELBO assembly (both divisors, reporting/monitoring only) is
#'                 affected by norm_convention itself; the beta update is
#'                 governed separately by `boost_beta` (see below).
#'                 "per_p" (default) -- boosts the survival ELBO term (reporting
#'                   only) by a factor of p; genomics is left unchanged.
#'                   Verified empirically to leave YFB's EL/EF exactly unchanged
#'                   (they never receive genomics_divisor) and to avoid the
#'                   L<->F collapse seen in LB under a genomics-shrinking
#'                   convention -- see DECISIONS.md 2026-07-12.
#'                 "np_n" -- literal genomics/(n*p), survival/n; retained for
#'                   empirical comparison (genomics-shrinking direction).
#' @param boost_beta logical (default FALSE). Beta's own coordinate update has
#'                 no genomics term competing with it in its own formula (ZF is
#'                 observed, 100% survival) -- the genomics/survival imbalance
#'                 `norm_convention` targets does not structurally exist there,
#'                 so boosting beta's precision does not correct any real
#'                 imbalance for it; it only reduces EBNM shrinkage, which
#'                 inflates K_eff (more factors cross a fixed beta_threshold)
#'                 without necessarily reflecting a genuine gain in survival
#'                 signal. Default FALSE leaves beta's precision unboosted
#'                 regardless of `norm_convention`. TRUE reproduces the earlier
#'                 (superseded) design that boosted it alongside the ELBO
#'                 monitor. See DECISIONS.md 2026-07-12.
#' @param init_method character: "svd" (default), "random", "deflation"
#'                 (sequential rank-1 SVD deflation, code/deflation_init.R --
#'                 candidate fix for the CAVI factor-collapse failure mode,
#'                 DECISIONS.md 2026-07-13), or "custom"
#' @param EL_init  Optional n x K matrix: custom initial loadings
#' @param EF_init  Optional p x K matrix: custom initial factors
#' @param N_burnin integer: beta-only burn-in iterations before joint CAVI (default 0).
#'                 0 = off (Cluster A baseline). N_burnin=5 or 10 may help if betas
#'                 collapse after the first joint CAVI iteration.
#' @param N_frozen integer: number of CAVI iterations during which EF is held fixed at
#'                 its SVD initialization while beta and L update freely (default 0 = off).
#'                 Breaks the beta=0 ↔ B_beta=0 CAVI fixed-point without triggering the
#'                 EF instability from dual-source F: EF is SVD-initialized (non-degenerate),
#'                 ZF = Y·EF_init is fixed and non-zero, so A_beta = sum(w·ZF_k^2) > 0
#'                 and beta can grow to a non-zero value. After N_frozen iters, EF unfreezes
#'                 and the full joint CAVI proceeds from a non-zero beta starting point.
#'                 Try N_frozen=10–30 for merged multi-platform data where beta->0 occurs.
#' @param cox_warmstart logical: initialize EBeta via Cox regression on ZF before CAVI
#'                 (default FALSE). FALSE = matches Cluster A behavior (EBeta starts at 0).
#'                 TRUE may help if normal prior produces unstable initial betas.
#' @param normalize_AB logical: rescale survival vs. genomics contributions in F update
#'                 (default FALSE). No-op when alpha_F=0 (current default; F update is
#'                 pure-genomics). Forward-compatible when alpha_F > 0 is enabled.
#' @param alpha_schedule NULL or list(warmup_iters, ramp_iters): ramp alpha from 0
#'                 up to `alpha` over warmup+ramp iterations.
#' @param sign_correction logical: apply post-convergence Phase C sign check — flip
#'                 EBeta if training concordance of ZF·EBeta is < 0.5 (default TRUE).
#'                 Set FALSE inside cross-validation folds so fold-to-fold EBeta
#'                 orientation stays consistent; CV evaluates via I(-risk_scores).
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
                           norm_convention = c("per_p", "np_n"),
                           boost_beta   = FALSE,
                           init_method  = "svd",
                           EL_init      = NULL,
                           EF_init      = NULL,
                           N_burnin     = 0,
                           N_frozen     = 0,
                           cox_warmstart    = FALSE,
                           normalize_AB     = FALSE,
                           alpha_F          = 0,
                           alpha_schedule   = NULL,
                           sign_correction  = TRUE,
                           verbose      = TRUE,
                           cohort_id       = NULL,
                           sigma_F_cohort  = 1.0) {

  n <- nrow(Y); p <- ncol(Y)

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1].")
  }
  if (!is.numeric(N_frozen) || length(N_frozen) != 1 || !is.finite(N_frozen) ||
      N_frozen < 0) {
    stop("N_frozen must be a non-negative integer.")
  }
  N_frozen <- as.integer(N_frozen)

  # ---- Phase 1a objective normalization (see DECISIONS.md) -----------------
  # "per_p": boosts survival (divisor < 1 => multiplies) instead of shrinking
  # genomics -- see fit_modular.R and DECISIONS.md 2026-07-12 for rationale.
  norm_convention  <- match.arg(norm_convention)
  genomics_divisor <- if (norm_convention == "np_n") n * p else 1
  survival_divisor <- if (norm_convention == "np_n") n   else 1 / p
  # beta's own coordinate update has no genomics term competing with it in its
  # own formula in either model (100% survival) -- the genomics/survival
  # imbalance this normalization targets does not structurally exist there.
  # boost_beta=FALSE (default) leaves beta's precision at its original scale;
  # =TRUE reproduces the earlier (superseded) design that also boosted it.
  # See DECISIONS.md 2026-07-12.
  beta_divisor <- if (boost_beta) survival_divisor else 1

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
  # Cohort indicator pre-initialisation (cohort-cols-L extension)
  # Identical logic to fit_modular.R; see that file for full commentary.
  # --------------------------------------------------------------------------
  if (!is.null(cohort_id)) {
    cohort_id <- factor(cohort_id)
    if (nlevels(cohort_id) < 2) {
      L_cohort <- matrix(0.0, n, 0)
      C_cols   <- 0L
    } else {
      L_cohort <- model.matrix(~ cohort_id)[, -1, drop = FALSE]
      C_cols   <- ncol(L_cohort)
    }
    n_c_vec <- if (C_cols > 0) colSums(L_cohort) else numeric(0)

    if (C_cols > 0) {
      EF_cohort_init  <- t(t(L_cohort) %*% Y / n_c_vec)
      EF2_cohort_init <- EF_cohort_init^2 + sigma_F_cohort^2
      Y_for_svd       <- Y - L_cohort %*% t(EF_cohort_init)
    } else {
      EF_cohort_init  <- matrix(0.0, p, 0)
      EF2_cohort_init <- matrix(0.0, p, 0)
      Y_for_svd       <- Y
    }
  } else {
    Y_for_svd <- Y
    L_cohort  <- NULL
    C_cols    <- 0L
  }

  # --------------------------------------------------------------------------
  # Initialization
  # --------------------------------------------------------------------------
  if (init_method == "svd") {
    svd_init <- svd(Y_for_svd, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    # Use abs() rather than pmax(..., 0): pmax zeros out entire columns when the
    # SVD vector points all-negative, giving degenerate ZF = Y·EF columns of all
    # zeros. Under point_exponential prior (non-negative), abs() is equally valid
    # and prevents the degenerate column pathology.
    EL  <- abs(svd_init$u %*% diag(d_k, K, K))
    EF  <- abs(svd_init$v %*% diag(d_k, K, K))
  } else if (init_method == "random") {
    y_sd <- sd(Y)
    EL   <- matrix(rnorm(n * K, sd = 0.1 * y_sd), n, K)
    EF   <- matrix(rnorm(p * K, sd = 0.1 * y_sd), p, K)
  } else if (init_method == "deflation") {
    # Sequential rank-1 SVD deflation (code/deflation_init.R): each factor is
    # fit to the residual after removing prior factors, so successive
    # factors cannot start out near-tied in amplitude the way batch SVD
    # columns from close singular values can -- see DECISIONS.md 2026-07-13.
    defl <- deflation_svd_init(Y_for_svd, K)
    EL <- abs(defl$EL)
    EF <- abs(defl$EF)
  } else if (init_method == "custom") {
    if (is.null(EL_init) || is.null(EF_init))
      stop("init_method='custom' requires both EL_init (n x K) and EF_init (p x K).")
    EL <- EL_init
    EF <- EF_init
  } else {
    stop(sprintf("Unknown init_method: '%s'. Use 'svd', 'random', 'deflation', or 'custom'.", init_method))
  }

  EL2 <- EL^2
  EF2 <- EF^2

  # Augment EL/EF matrices with cohort indicator columns (same as fit_modular.R).
  if (!is.null(cohort_id) && C_cols > 0) {
    EL_aug  <- cbind(EL,  L_cohort)
    EL2_aug <- cbind(EL2, L_cohort)
    EF_aug  <- cbind(EF,  EF_cohort_init)
    EF2_aug <- cbind(EF2, EF2_cohort_init)
  } else {
    EL_aug <- EL; EL2_aug <- EL2; EF_aug <- EF; EF2_aug <- EF2
  }

  # Optionally warm-start beta via Cox regression on ZF = Y·EF.
  # cox_warmstart=FALSE (default) matches Cluster A behavior: EBeta starts at 0.
  # cox_warmstart=TRUE calibrates EBeta to the ZF scale before CAVI begins.
  if (cox_warmstart) {
    # Normalize EF before projecting: CAVI uses ZF = Y·EF_norm (unit-L2-norm columns).
    # Without normalization, ‖ZF‖ ~ O(√(p·n)) and Cox returns β ~ 1e-9 — machine
    # epsilon on the CAVI scale — so the warm-start has no effect.
    EF_norms_ws  <- sqrt(colSums(EF^2) + 1e-10)
    EF_norm_ws   <- sweep(EF, 2, EF_norms_ws, "/")
    ZF_init <- Y %*% EF_norm_ws                      # n × K, normalized scale
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
      # Normalize EF columns before projection so A_beta = sum(w * ZF_k^2)
      # is O(n) not O(n*p). Without normalization, ‖ZF_k‖ ≈ O(√(p·n)) drives
      # EBeta to zero via spike-and-slab shrinkage.
      EF_norms_b <- sqrt(colSums(EF^2) + 1e-10)
      EF_norm_b  <- sweep(EF, 2, EF_norms_b, "/")
      ZF_b     <- Y %*% EF_norm_b
      eta_b    <- as.vector(ZF_b %*% EBeta)
      taylor_b <- calc_cox_taylor_yf(eta_b, time, status)
      z_b      <- eta_b + taylor_b$u / taylor_b$w
      w_b      <- taylor_b$w
      for (k in seq_len(K)) {
        z_no_k_b <- compute_z_no_k(z_b, ZF_b, EBeta, k)
        # Cluster B: ZF[,k] is observed, so its "second moment" = ZF[,k]^2 (no posterior variance)
        res_b    <- update_beta_k(w_b, z_no_k_b, ZF_b[, k], ZF_b[, k]^2,
                                  prior_family = prior_beta, alpha = alpha,
                                  survival_divisor = beta_divisor)
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
    cat(sprintf("    prior_LF=%s | prior_beta=%s | alpha=%.2f | alpha_F=%.2f\n\n",
                prior_LF, prior_beta, alpha, alpha_F))
  }

  if (verbose && N_frozen > 0) {
    cat(sprintf("    [frozen-F] EF frozen for first %d iters — beta and L update freely\n",
                N_frozen))
  }

  # ==========================================================================
  # Main CAVI Loop
  # ==========================================================================
  for (iter in 1:max_iter) {

    EL_old    <- EL
    EBeta_old <- EBeta

    # TRUE during the frozen phase: EF is held at SVD init, only beta and L update.
    # This breaks the beta=0 ↔ B_beta=0 CAVI fixed-point: ZF = Y·EF_init is non-zero
    # from SVD, so A_beta = sum(w·ZF_k^2) > 0 and beta can grow freely before F adapts.
    freeze_F <- iter <= N_frozen
    if (verbose && N_frozen > 0 && iter == N_frozen + 1L) {
      cat(sprintf("    [frozen-F] Unfreezing EF at iter %d — beta_max=%.4e\n",
                  iter, max(abs(EBeta))))
    }

    alpha_iter <- alpha_at(iter)

    kl_L    <- numeric(K)
    kl_F    <- numeric(K)
    kl_beta <- numeric(K)

    # ------------------------------------------------------------------------
    # STEP 1: Cox Taylor Expansion
    #
    # Cluster B: ZF = Y · EF_norm where EF_norm has unit-L2-norm columns.
    # Normalizing EF before projection keeps A_beta = sum(w * ZF_k^2) at
    # O(n) scale rather than O(n*p). Without normalization, ‖EF_k‖ ≈ O(√p)
    # makes ‖ZF_k‖ ≈ O(√(p*n)), causing spike-and-slab to drive EBeta → 0.
    # EF_norms stored in model object so prediction applies identical scaling.
    # ------------------------------------------------------------------------
    # Location 1: ZF must use only the K global-factor columns of EF_aug.
    # If augmented, EF_aug has K+C_cols columns; restricting to 1:K keeps ZF
    # n x K so ZF %*% EBeta (EBeta is K-vector) is well-defined.
    # With cohort_id=NULL, EF_aug = EF (K columns), so this is a no-op.
    EF_global <- EF_aug[, 1:K, drop = FALSE]          # p x K
    EF_norms  <- sqrt(colSums(EF_global^2) + 1e-10)   # K-vector
    EF_norm   <- sweep(EF_global, 2, EF_norms, "/")    # p x K, unit-norm columns
    ZF        <- Y %*% EF_norm                         # n x K: normalized projection scores
    eta    <- as.vector(ZF %*% EBeta)           # eta = ZF * beta_tilde
    taylor <- calc_cox_taylor_yf(eta, time, status)
    z      <- eta + taylor$u / taylor$w    # working response z_i
    w      <- taylor$w                     # W_{ii} diagonal Hessian weights
    YtWY_diag <- as.vector(t(Y^2) %*% w)  # p-vector: diag(Y'diag(w)Y)

    history$rmse[iter] <- sqrt(mean((Y - EL_aug %*% t(EF_aug))^2))

    # ========================================================================
    # STEP 2: Factor-Wise Coordinate Ascent  k = 1, ..., K
    #
    # Update order per factor k:
    #   (a) q(beta_k): survival coefficient — covariate is ZF[,k] (observed)
    #   (b) q(l_k):   pure-genomics only (L not in Cox under Cluster B)
    #   (c) q(f_k):   pure-genomics (alpha_F=0); survival terms receive zero weight
    # ========================================================================
    for (k in 1:K) {

      R_k    <- compute_R_k(Y, EL_aug, EF_aug, k)
      # Cluster B: z_no_k computed w.r.t. ZF (not EL as in Cluster A)
      z_no_k <- compute_z_no_k(z, ZF, EBeta, k)

      # ---- (a) Update q(beta_k): Survival Coefficient ----
      # Covariate is ZF[,k] = (Y·EF)[,k] — observed projection (not latent).
      # Pass ZF[,k]^2 as EL2_k (the second moment argument): ZF is observed, so
      # its "posterior second moment" = squared value (no variance term).
      # A_beta = sum(w * ZF_k^2) is non-zero from SVD init regardless of EBeta,
      # breaking the chicken-and-egg that plagued Cluster A's L update.
      res_beta    <- update_beta_k(w, z_no_k, ZF[, k], ZF[, k]^2,
                                   prior_family = prior_beta, alpha = alpha_iter,
                                   survival_divisor = beta_divisor)
      EBeta[k]    <- res_beta$mean
      EBeta2[k]   <- res_beta$second
      kl_beta[k]  <- compute_ebnm_kl(res_beta$ebnm_result$log_likelihood,
                                      res_beta$A, res_beta$x,
                                      res_beta$mean, res_beta$second)

      # ---- (b) Update q(l_k): Patient Loadings — PURE GENOMICS ----
      # L does not appear in the Cox likelihood under eta = ZF * beta_tilde.
      # Uses update_L_surv_YFB_k() (NOT update_L_k from Cluster A).
      res_L   <- update_L_surv_YFB_k(Tau, EF_aug[, k], EF2_aug[, k], R_k,
                                      prior_family = prior_LF)
      EL[, k]      <- res_L$mean
      EL2[, k]     <- res_L$second
      EL_aug[, k]  <- res_L$mean    # keep augmented matrix in sync
      EL2_aug[, k] <- res_L$second
      kl_L[k]  <- compute_ebnm_kl(res_L$ebnm_result$log_likelihood,
                                   res_L$A, res_L$x, res_L$mean, res_L$second)

      # ---- (c) Update q(f_k): Biological Factors ----
      # Skipped during frozen-F phase (iter <= N_frozen): EF held at SVD init.
      # R_k depends on EL[,k] (updated above); recompute before F update.
      # alpha_F=0 (default): pure genomics — prevents positive-feedback instability
      #   documented in DECISIONS.md 2026-04-30.
      # alpha_F>0: dual-source — survival gradient stabilises ZF on merged
      #   multi-platform data where genomics-only EF drifts to platform contrast.
      #   The instability required degenerate zero-column SVD init (pmax bug, now
      #   fixed with abs()). Test stability empirically before raising alpha_F.
      if (!freeze_F) {
        R_k        <- compute_R_k(Y, EL_aug, EF_aug, k)
        YtWz_no_k  <- as.vector(t(Y) %*% (w * z_no_k))  # p-vector
        res_F <- update_F_surv_YFB_k(Tau, EL_aug[, k], EL2_aug[, k], R_k,
                                      EBeta_k      = EBeta[k],
                                      EBeta2_k     = EBeta2[k],
                                      YtWY_diag    = YtWY_diag,
                                      YtWz_no_k    = YtWz_no_k,
                                      prior_family = prior_LF,
                                      alpha        = alpha_F)
        EF[, k]      <- res_F$mean
        EF2[, k]     <- res_F$second
        EF_aug[, k]  <- res_F$mean    # keep augmented matrix in sync
        EF2_aug[, k] <- res_F$second
        kl_F[k]  <- compute_ebnm_kl(res_F$ebnm_result$log_likelihood,
                                     res_F$A, res_F$x, res_F$mean, res_F$second)
      }

    }  # end k-loop

    # Cohort F update — identical block to fit_modular.R.
    # Also frozen during the N_frozen phase to keep EF_aug fully fixed.
    if (!is.null(cohort_id) && C_cols > 0 && !freeze_F) {
      Y_hat_global <- EL %*% t(EF)
      res_cohort   <- update_F_cohort_all(
        L_cohort, Y - Y_hat_global, Tau, sigma_F_cohort
      )
      EF_aug[,  K + seq_len(C_cols)] <- res_cohort$EF_cohort
      EF2_aug[, K + seq_len(C_cols)] <- res_cohort$EF2_cohort
    }

    # ========================================================================
    # STEP 3: Noise Precision Update
    # ========================================================================
    res_tau              <- update_tau(Y, EL_aug, EL2_aug, EF_aug, EF2_aug)
    Tau                  <- res_tau$Tau
    history$elbo_proxy[iter] <- res_tau$elbo_proxy
    factor_pve_iter <- vapply(seq_len(K), function(k) {
      sum(EL[, k]^2) * sum(EF[, k]^2) / y_frob2
    }, numeric(1))
    history$factor_pve[iter, ] <- factor_pve_iter

    # Development guard: ZF must be n x K; EBeta must be K-vector.
    stopifnot(ncol(ZF) == K, length(EBeta) == K)

    # Full ELBO under Cluster B: ZF = Y·EF is the survival predictor.
    # Pass ZF^2 as EL2 (element-wise squared): ZF is observed, so its
    # posterior second moment equals its squared value (no variance term).
    surv_elbo               <- compute_survival_elbo(taylor$logPL, w,
                                                     ZF, ZF^2, EBeta, EBeta2)
    history$elbo_full[iter] <- (1 - alpha) * (res_tau$elbo_proxy / genomics_divisor) +
                               alpha * (surv_elbo / survival_divisor) +
                               sum(kl_L) + sum(kl_F) + sum(kl_beta)
    if (!is.null(cohort_id) && C_cols > 0) {
      history$elbo_full[iter] <- history$elbo_full[iter] +
        compute_normal_kl(
          EF_aug[,  K + seq_len(C_cols), drop = FALSE],
          EF2_aug[, K + seq_len(C_cols), drop = FALSE],
          sigma_F_cohort
        )
    }

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

  # Phase C: Check training concordance; flip EBeta sign if anti-concordant.
  # With YFB (eta = ZF·beta), the Gram matrix EF'EF and the pmax(SVD) initialization
  # (which discards sign information from F_true) can invert the relationship between
  # ZF[:,k] and the true survival direction. A global sign flip of EBeta corrects this:
  # risk_score_new = ZF * (-EBeta) = -(risk_score_original), yielding 1 - C_train.
  # Deviation from Plan option A (PC1 correlation): training concordance check is more
  # direct and correct for any dataset, not just synthetic.
  #
  # sign_correction = FALSE during cross-validation (select_K_cv), because
  # per-fold sign flips produce fold-to-fold EBeta orientation inconsistency that
  # inflates apparent C-index variance and corrupts K selection. CV evaluates
  # held-out concordance via I(-pred$risk_scores), which implicitly assumes the
  # raw SVD-initialized EBeta orientation — consistent with the LB model convention.
  if (sign_correction) {
    # Location 2: restrict to K global-factor columns before computing ZF_final.
    # EF_norms was computed from EF_global (K columns) earlier in this iteration.
    ZF_final  <- Y %*% sweep(EF_aug[, 1:K, drop = FALSE], 2, EF_norms, "/")
    eta_final <- as.vector(ZF_final %*% EBeta)
    c_train   <- as.numeric(
      concordance(Surv(time, status) ~ eta_final)$concordance
    )
    if (c_train < 0.5) {
      if (verbose) {
        cat(sprintf("    [Phase C] Training C=%.4f < 0.5 — flipping EBeta sign\n", c_train))
      }
      EBeta  <- -EBeta
      EBeta2 <- EBeta^2
    }
  }

  stopifnot(length(EF_norms) == K)  # guard: must be K-vector for prediction

  result <- list(
    EL       = EL,
    EL2      = EL2,
    EF       = EF,
    EF2      = EF2,
    EBeta    = EBeta,
    EBeta2   = EBeta2,
    Tau      = Tau,
    EF_norms = EF_norms,   # K-vector: final iteration column norms for prediction
    history  = history,
    norm_convention  = norm_convention,
    genomics_divisor = genomics_divisor,
    survival_divisor = survival_divisor,
    beta_divisor     = beta_divisor
  )
  if (!is.null(cohort_id) && C_cols > 0) {
    result$L_cohort   <- L_cohort
    result$EF_cohort  <- EF_aug[, K + seq_len(C_cols), drop = FALSE]
    result$EF2_cohort <- EF2_aug[, K + seq_len(C_cols), drop = FALSE]
    result$cohort_id  <- cohort_id
  }
  result
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
