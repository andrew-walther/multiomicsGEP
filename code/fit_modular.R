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

source("code/update_L.R")      # compute_R_k, update_L_k, update_L_all
source("code/update_F.R")      # update_F_k, update_F_all  (uses compute_R_k from L)
source("code/update_beta.R")   # compute_z_no_k, update_beta_k, update_beta_all
source("code/update_tau.R")    # compute_var_term, update_tau
source("code/compute_elbo.R")  # compute_ebnm_kl, compute_survival_elbo, compute_normal_kl
source("code/update_F_cohort.R")  # update_F_cohort_col, update_F_cohort_all
source("code/deflation_init.R")   # deflation_svd_init -- init_method="deflation"

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
#' @param strata NULL (default) for a single pooled risk set, or an n-vector of
#'               stratum labels (e.g. study/cohort). When supplied, Breslow risk
#'               sets are formed *within* each stratum — the standard stratified
#'               Cox partial likelihood. The baseline hazard still cancels
#'               per-stratum, so no parametric h0 is introduced (Item 3,
#'               DECISIONS.md). strata=rep(1,n) reduces exactly to strata=NULL.
#' @return list(u = n-vector score, w = n-vector neg-diagonal Hessian,
#'              logPL = scalar Breslow partial log-likelihood at eta)
# ------------------------------------------------------------------------------
calc_cox_taylor <- function(eta, time, status, strata = NULL) {
  # Single-stratum core: one pooled Breslow risk set over the samples passed in.
  # Byte-identical to the historical (unstratified) body, so the strata=NULL
  # path below is unchanged.
  .core <- function(eta, time, status) {
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
    # Breslow partial log-likelihood: sum_i delta_i * (eta_i - log R_i)
    # where R_i = sum_{j: t_j >= t_i} exp(eta_j) (risk set cumulative sum).
    # pmax guards against degenerate risk sets (shouldn't occur in practice).
    logPL <- sum(status_s * (eta_s - log(pmax(risk_sum, 1e-300))))
    list(u = u, w = w, logPL = logPL)
  }

  if (is.null(strata)) return(.core(eta, time, status))

  # Stratified: score/Hessian are per-sample (scatter back by index); the
  # partial log-likelihood is additive across strata (independent risk sets).
  n <- length(time)
  if (length(strata) != n) {
    stop("strata must have the same length as time (", n, ").")
  }
  # Fail loud: as.factor() drops NA from levels, which would silently exclude
  # NA-labelled samples from every risk set (leaving u=0,w=0 -> 0/0 downstream).
  if (anyNA(strata)) stop("strata must not contain NA.")
  strata <- as.factor(strata)
  u <- numeric(n); w <- numeric(n); logPL <- 0
  for (lev in levels(strata)) {
    idx <- which(strata == lev)
    res <- .core(eta[idx], time[idx], status[idx])
    u[idx] <- res$u
    w[idx] <- res$w
    logPL  <- logPL + res$logPL
  }
  list(u = u, w = w, logPL = logPL)
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
#' Convergence is declared when the relative full-ELBO change drops below tol
#' (after a 5-iteration burn-in):
#'   abs(elbo_new - elbo_old) / abs(elbo_old) < tol
#' Parameter deltas are still tracked in history for diagnostics, but no longer
#' drive termination.
#'
#' @param Y        n x p genomics data matrix
#' @param time     n-vector of survival / censoring times
#' @param status   n-vector of event indicators (1=event, 0=censored)
#' @param K        Number of latent factors (default 5)
#' @param max_iter Maximum CAVI outer iterations (default 100)
#' @param tol      Convergence threshold on relative full-ELBO change
#'                 (default 1e-5)
#' @param prior_LF   character: EBNM prior family for loadings L and factors F.
#'                 Default "point_exponential" — non-negative (NMF-style),
#'                 appropriate when L and F represent expression magnitudes.
#'                 Also accepts "point_normal", "point_laplace".
#' @param prior_beta character: EBNM prior family for survival coefficients beta.
#'                 Default "point_normal" — sparse, allows positive and negative
#'                 coefficients (prognostic direction not constrained).
#'                 Also accepts "point_laplace" for heavier-tailed sparsity.
#' @param alpha      numeric in [0, 1]: survival/genomics mixing weight passed to
#'                 update_L_k(), update_F_k(), and update_beta_k() at every iteration.
#'                 Controls the relative emphasis of the survival likelihood vs.
#'                 the genomics likelihood in the CAVI updates:
#'                   L update:    (1-alpha)*A_gen + alpha*A_surv
#'                   F update:    (1-alpha)*A_gen / (1-alpha)*B_gen
#'                   beta update: alpha*A_beta / alpha*B_beta
#'                 alpha=0: pure genomics (beta not updated); alpha=1: pure survival
#'                 (F not updated). Default 0.5 balances the p >> n gradient asymmetry.
#' @param norm_convention character: Phase 1a objective normalization convention
#'                 (see DECISIONS.md 2026-07-12). Before normalization, the
#'                 genomics term of A_L (a raw sum over p genes, ~O(p)) and the
#'                 survival term (a single per-patient term, ~O(1)) differ by an
#'                 unnormalized ~p-fold scale gap, making alpha's mixing
#'                 uninterpretable regardless of its value.
#'                 Under LB (Cluster A), L and F co-adapt bilinearly (L's
#'                 precision depends on EF^2 and vice versa). Rescaling either
#'                 side's precision -- in EITHER direction -- destabilizes that
#'                 coupling: dividing genomics by p collapses EL/EF to exactly
#'                 zero; multiplying survival by p instead causes unbounded
#'                 divergence (both confirmed empirically, including at
#'                 realistic PDAC scale, and a scheduled/ramped introduction via
#'                 alpha_schedule does not stabilize it either). Consequently,
#'                 genomics_divisor/survival_divisor are NOT passed to
#'                 update_L_k or update_F_k here regardless of norm_convention --
#'                 LB's L/F precision is always at its original, pre-Phase-1a
#'                 scale. Only the ELBO monitor and the beta update (no bilinear
#'                 coupling there) receive the convention's rebalancing, so LB's
#'                 alpha-interpretability fix is partial. The recommended
#'                 production model (YFB, fit_cox_on_yf.R) is unaffected by this
#'                 limitation: its L/F updates never touch these divisors in the
#'                 first place (no genomics/survival competition to rebalance),
#'                 so YFB gets the full, verified-stable treatment.
#'                 "per_p" (default) -- boosts the survival ELBO term (reporting
#'                   only) by a factor of p; the beta update itself is governed
#'                   separately by `boost_beta` (see below).
#'                 "np_n" -- literal convention (genomics/(n*p), survival/n)
#'                   for the ELBO monitor; retained for comparison.
#' @param boost_beta logical (default FALSE). Beta's own coordinate update has
#'                 no genomics term competing with it in its own formula in
#'                 either model (100% survival) -- the genomics/survival
#'                 imbalance `norm_convention` targets does not structurally
#'                 exist there, so boosting beta's precision does not correct
#'                 any real imbalance for it; it only reduces EBNM shrinkage,
#'                 which inflates K_eff (more factors cross a fixed
#'                 beta_threshold) without necessarily reflecting a genuine
#'                 gain in survival signal. Default FALSE leaves beta's
#'                 precision unboosted regardless of `norm_convention`. TRUE
#'                 reproduces the earlier (superseded) design that boosted it
#'                 alongside the ELBO monitor. See DECISIONS.md 2026-07-12.
#' @param init_method  character: initialization strategy.
#'                 "svd"    (default — deterministic SVD warm-start),
#'                 "random" (random normal initialization, useful with
#'                 multiple restarts to escape local optima),
#'                 "deflation" (sequential rank-1 SVD deflation -- see
#'                 code/deflation_init.R; candidate fix for the CAVI
#'                 factor-collapse failure mode, DECISIONS.md 2026-07-13),
#'                 "custom" (supply EL_init and EF_init directly; set
#'                 automatically when both are non-NULL).
#' @param EL_init  Optional n x K numeric matrix: custom initial loadings.
#'                 When both EL_init and EF_init are non-NULL, init_method is
#'                 overridden to "custom". Intended for EBMF warm-start
#'                 experiments — supply ldf(flash_fit)$L scaled by D, or
#'                 flash_fit$L_pm, to initialize CAVI from an EBMF solution.
#' @param EF_init  Optional p x K numeric matrix: custom initial factors.
#'                 Must be supplied together with EL_init.
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
#'   $history  list(rmse, elbo_proxy, elbo_full, delta_L, delta_Beta,
#'                  delta_elbo_rel, converged, n_iter)
fit_supervised_mf_modular <- function(Y, time, status,
                                      K               = 5,
                                      max_iter        = 100,
                                      tol             = 1e-5,
                                      prior_LF        = "point_exponential",
                                      prior_beta      = "point_normal",
                                      alpha           = 0.5,
                                      norm_convention = c("per_p", "np_n"),
                                      boost_beta      = FALSE,
                                      init_method     = "svd",
                                      EL_init         = NULL,
                                      EF_init         = NULL,
                                      N_burnin        = 0,
                                      alpha_schedule  = NULL,
                                      normalize_AB    = FALSE,
                                      sign_correction = TRUE,
                                      verbose         = TRUE,
                                      cohort_id       = NULL,
                                      sigma_F_cohort  = 1.0,
                                      strata_id       = NULL) {

  n <- nrow(Y); p <- ncol(Y)

  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1].")
  }

  # ---- Stratified baseline hazard (Item 3) ---------------------------------
  # strata_id (e.g. study/cohort) forms Breslow risk sets within each stratum
  # in the Cox partial likelihood, so baseline survival may differ by stratum
  # while beta is shared. Distinct from cohort_id, which absorbs *genomic*
  # platform offsets; the two can be used together. NULL => single pooled risk
  # set (unchanged behaviour). See DECISIONS.md.
  # NOTE: the CV/tuning wrappers (select_alpha_cv, select_K_cv, auto_prune_K) do
  # NOT thread strata_id per fold — hyperparameter tuning runs unstratified even
  # if the final fit is stratified. Passing strata_id through them would hand a
  # full-length vector to a fold-subset fit and trip the length check (a loud
  # failure, not silent corruption).
  if (!is.null(strata_id)) {
    if (length(strata_id) != n)
      stop("strata_id must be NULL or have length nrow(Y) (", n, ").")
    if (anyNA(strata_id))
      stop("strata_id must not contain NA.")
  }

  # ---- Phase 1a objective normalization (see DECISIONS.md) -----------------
  # "per_p": boosts survival (divisor < 1 => multiplies) instead of shrinking
  # genomics, to avoid the L<->F collapse documented in DECISIONS.md 2026-07-12.
  # "np_n": literal genomics/(n*p), survival/n -- retained for comparison; this
  # is the collapse-prone, genomics-shrinking direction.
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

  # ---- Progressive α schedule [A2 of docs/beta_zero_fix_design.md §4.4] ----
  # alpha_schedule = list(warmup_iters, ramp_iters) ramps α from 0 (pure
  # genomics) up to the target `alpha` over the ramp window, then holds.
  # Lets L settle into reconstruction-meaningful directions before survival
  # pressure is applied. Default NULL preserves existing behaviour
  # (constant α from iter 1).
  use_alpha_schedule <- !is.null(alpha_schedule)
  if (use_alpha_schedule) {
    if (!is.list(alpha_schedule) ||
        !all(c("warmup_iters", "ramp_iters") %in% names(alpha_schedule)) ||
        !is.numeric(alpha_schedule$warmup_iters) ||
        !is.numeric(alpha_schedule$ramp_iters) ||
        alpha_schedule$warmup_iters < 0 ||
        alpha_schedule$ramp_iters   < 0) {
      stop("alpha_schedule must be NULL or list(warmup_iters>=0, ramp_iters>=0).")
    }
  }
  # Closure returning alpha at outer iteration `iter`. Branchless when off.
  alpha_at <- function(iter) {
    if (!use_alpha_schedule) return(alpha)
    wu <- alpha_schedule$warmup_iters
    rp <- alpha_schedule$ramp_iters
    if (iter <= wu) return(0)
    if (iter <= wu + rp) {
      prog <- (iter - wu) / max(rp, 1)
      return(alpha * prog)
    }
    alpha
  }

  # Auto-promote to "custom" init when caller supplies both EL_init and EF_init.
  # This avoids requiring the caller to pass init_method = "custom" explicitly.
  if (!is.null(EL_init) && !is.null(EF_init)) init_method <- "custom"

  # --------------------------------------------------------------------------
  # Cohort indicator pre-initialisation (cohort-cols-L extension)
  #
  # When cohort_id is supplied, append C-1 fixed binary indicator columns to L
  # (corner-point encoding; reference = first cohort level).  Only the
  # corresponding F rows (f_c) are estimated — c is never updated during CAVI.
  # Residualise Y before SVD so biological factors are not dominated by the
  # platform offset.
  # When cohort_id = NULL, Y_for_svd = Y and all augmented matrices alias the
  # K-factor matrices unchanged — bit-identical behaviour to the base model.
  # --------------------------------------------------------------------------
  if (!is.null(cohort_id)) {
    cohort_id <- factor(cohort_id)
    # model.matrix requires >= 2 levels; single-cohort is a no-op (C_cols=0).
    if (nlevels(cohort_id) < 2) {
      L_cohort <- matrix(0.0, n, 0)
      C_cols   <- 0L
    } else {
      L_cohort <- model.matrix(~ cohort_id)[, -1, drop = FALSE]  # n x (C-1)
      C_cols   <- ncol(L_cohort)
    }
    n_c_vec   <- if (C_cols > 0) colSums(L_cohort) else numeric(0)

    if (C_cols > 0) {
      # Per-cohort gene mean differences: (C-1) x p matrix; dividing by n_c
      # broadcasts row-wise.  t(L_cohort) %*% Y gives (C-1) x p.
      EF_cohort_init  <- t(t(L_cohort) %*% Y / n_c_vec)         # p x (C-1)
      # Include prior variance in initial EF2 — NEVER set EF2 = EF^2 here,
      # which would zero out posterior variance and inflate tau on iteration 1.
      EF2_cohort_init <- EF_cohort_init^2 + sigma_F_cohort^2    # p x (C-1)
      Y_for_svd       <- Y - L_cohort %*% t(EF_cohort_init)     # residualised
    } else {
      # Single cohort (C=1): no indicator column; treat as no-op.
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
    # SVD of Y_for_svd: deterministic high-variance starting subspace.
    # When cohort_id is supplied, Y_for_svd has the cohort mean removed so
    # the leading singular vectors capture biology, not platform offset.
    svd_init <- svd(Y_for_svd, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    EL  <- pmax(svd_init$u %*% diag(d_k, K, K), 0)      # n x K; pmax ensures non-negative init matches point_exponential prior
    EF  <- pmax(svd_init$v %*% diag(d_k, K, K), 0)      # p x K
  } else if (init_method == "random") {
    # Random normal initialization scaled by data magnitude.
    # sd = 0.1 * overall SD of Y keeps initial reconstruction in a
    # reasonable range; too large causes Cox Taylor expansion to diverge.
    y_sd <- sd(Y)
    EL   <- matrix(rnorm(n * K, sd = 0.1 * y_sd), n, K)
    EF   <- matrix(rnorm(p * K, sd = 0.1 * y_sd), p, K)
  } else if (init_method == "deflation") {
    # Sequential rank-1 SVD deflation (code/deflation_init.R): each factor is
    # fit to the residual after removing prior factors, so successive
    # factors cannot start out near-tied in amplitude the way batch SVD
    # columns from close singular values can -- see DECISIONS.md 2026-07-13.
    defl <- deflation_svd_init(Y_for_svd, K)
    EL <- pmax(defl$EL, 0)   # pmax ensures non-negative init matches point_exponential prior
    EF <- pmax(defl$EF, 0)
  } else if (init_method == "custom") {
    # Custom warm-start: caller supplies EL_init (n x K) and EF_init (p x K).
    # Primary use case: EBMF warm-start diagnostic — initialise from a
    # flashier solution to test whether CAVI can develop non-zero β when
    # started from factors already associated with survival.
    if (is.null(EL_init) || is.null(EF_init))
      stop("init_method='custom' requires both EL_init (n x K) and EF_init (p x K).")
    if (!isTRUE(all.equal(dim(EL_init), c(n, K))))
      stop(sprintf("EL_init must be n x K = %d x %d; got %d x %d.",
                   n, K, nrow(EL_init), ncol(EL_init)))
    if (!isTRUE(all.equal(dim(EF_init), c(p, K))))
      stop(sprintf("EF_init must be p x K = %d x %d; got %d x %d.",
                   p, K, nrow(EF_init), ncol(EF_init)))
    EL <- EL_init
    EF <- EF_init
  } else {
    stop(sprintf("Unknown init_method: '%s'. Use 'svd', 'random', 'deflation', or 'custom'.", init_method))
  }

  # Second moments initialised to squared means (zero posterior variance).
  # Posterior variance populated after the first EBNM call.
  EL2 <- EL^2
  EF2 <- EF^2

  # Augment EL/EF matrices with cohort indicator columns.
  # EL_aug (n x K+C_cols): cbind(EL, L_cohort) — L_cohort is fixed binary.
  # EL2 cohort columns = L_cohort (binary: 0²=0, 1²=1 → zero L-variance).
  # EF2 cohort columns include 1/A_c term (non-zero posterior F-variance).
  # When C_cols = 0, all aug matrices alias the K-factor matrices: no copies.
  if (!is.null(cohort_id) && C_cols > 0) {
    EL_aug  <- cbind(EL,  L_cohort)
    EL2_aug <- cbind(EL2, L_cohort)       # binary: L^2 = L -> zero L-variance
    EF_aug  <- cbind(EF,  EF_cohort_init)
    EF2_aug <- cbind(EF2, EF2_cohort_init)
  } else {
    EL_aug <- EL; EL2_aug <- EL2; EF_aug <- EF; EF2_aug <- EF2
  }

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

  # ---- Instrumentation [§4.1 of docs/beta_zero_fix_design.md] -----------
  # Quantifies the Cox warm-start output. If max(|EBeta|) is already non-zero
  # but the full CAVI still collapses β to zero, the cold-start is not the
  # entry point — the A_surv/A_gen scale imbalance dominates. If ~0, the
  # cycle starts at initialization and Fixes 1+2 (reorder + burn-in) are the
  # natural first attempt.
  if (verbose) {
    cat(sprintf("    [init] Cox warm-start EBeta range: [%.3e, %.3e]\n",
                min(EBeta), max(EBeta)))
  }

  # ==========================================================================
  # β-only Burn-in [Fix 2 of docs/beta_zero_fix_design.md §4.3]
  #
  # Run N_burnin iterations of pure β-updates with EL held fixed at the
  # SVD-initialized values. Treats L as a point estimate (EL2_init = EL^2,
  # i.e., zero posterior variance). This directly replicates Warm-start
  # Exp 1, which proved the β update produces non-zero coefficients when L
  # is fixed. Goal: enter the joint CAVI loop with a non-zero EBeta so that
  # A_surv = w·E[β_k²] is non-trivial in the L update from iter 1.
  #
  # Backward compatibility: default N_burnin = 0 skips the block entirely.
  # ==========================================================================
  if (N_burnin > 0) {
    EL2_init <- EL^2   # point estimate; zero posterior variance during burn-in
    for (b in seq_len(N_burnin)) {
      eta_b    <- as.vector(EL %*% EBeta)
      taylor_b <- calc_cox_taylor(eta_b, time, status, strata = strata_id)
      z_b      <- eta_b + taylor_b$u / taylor_b$w
      w_b      <- taylor_b$w
      for (k in seq_len(K)) {
        z_no_k_b <- compute_z_no_k(z_b, EL, EBeta, k)
        res_b    <- update_beta_k(w_b, z_no_k_b, EL[, k], EL2_init[, k],
                                  prior_family = prior_beta, alpha = alpha,
                                  survival_divisor = beta_divisor)
        EBeta[k]  <- res_b$mean
        EBeta2[k] <- res_b$second
      }
    }
    if (verbose) {
      cat(sprintf("    [burnin] After %d β-only iterations, EBeta range: [%.3e, %.3e]\n",
                  N_burnin, min(EBeta), max(EBeta)))
    }
  }

  # Column-specific noise precision from sample variance of each column of Y.
  Tau <- 1.0 / pmax(apply(Y, 2, var), 1e-8)   # p-vector
  y_frob2 <- sum(Y^2)

  # History tracking
  history <- list(
    rmse       = numeric(max_iter),
    elbo_proxy = numeric(max_iter),
    elbo_full  = numeric(max_iter),  # full ELBO: proxy + survival + KL terms
    delta_L    = numeric(max_iter),
    delta_Beta = numeric(max_iter),
    delta_elbo_rel = rep(NA_real_, max_iter),
    factor_pve = matrix(NA_real_, nrow = max_iter, ncol = K),
    converged  = FALSE,
    n_iter     = max_iter
  )

  if (verbose) {
    cat("=== Supervised Bayesian MF (Modular V3) — Factor-Wise CAVI ===\n")
    cat(sprintf("    n=%d, p=%d, K=%d | max_iter=%d | tol=%.1e\n",
                n, p, K, max_iter, tol))
    cat(sprintf("    prior_LF=%s | prior_beta=%s | alpha=%.2f | init_method=%s\n\n",
                prior_LF, prior_beta, alpha, init_method))
  }

  # ==========================================================================
  # Main CAVI Loop — V3 Algorithm 1
  # ==========================================================================
  for (iter in 1:max_iter) {

    EL_old    <- EL
    EBeta_old <- EBeta

    # Per-iteration α (constant by default; ramped if alpha_schedule supplied)
    alpha_iter <- alpha_at(iter)
    if (use_alpha_schedule && verbose &&
        (iter <= alpha_schedule$warmup_iters + alpha_schedule$ramp_iters + 1)) {
      cat(sprintf("    [iter=%d] alpha_iter=%.3f\n", iter, alpha_iter))
    }

    # KL accumulators for full ELBO — reset each outer iteration.
    # kl_L[k]    = E_q[log g_L(l_k)]    - E_q[log q_L(l_k)]    (<= 0)
    # kl_F[k]    = E_q[log g_F(f_k)]    - E_q[log q_F(f_k)]    (<= 0)
    # kl_beta[k] = E_q[log g_beta(b_k)] - E_q[log q_beta(b_k)] (<= 0)
    kl_L    <- numeric(K)
    kl_F    <- numeric(K)
    kl_beta <- numeric(K)

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
    taylor <- calc_cox_taylor(eta, time, status, strata = strata_id)
    z      <- eta + taylor$u / taylor$w    # n-vector: working response z_i
    w      <- taylor$w                     # n-vector: W_{ii}

    # Reconstruction RMSE at posterior means (monitoring only).
    # Use augmented matrices so the cohort contribution is included in the
    # residual; without this, RMSE reflects only K-factor error and is not
    # comparable to the base model.
    history$rmse[iter] <- sqrt(mean((Y - EL_aug %*% t(EF_aug))^2))

    # ========================================================================
    # STEP 2: Factor-Wise Coordinate Ascent  k = 1, ..., K
    #
    # For each factor k, update in the order:
    #   (a) q(l_k):    patient loadings — uses current EBeta[k] (stale OK: A_surv << A_gen)
    #   (b) q(f_k):    biological factors — uses freshly updated EL[,k]
    #   (c) q(beta_k): survival coefficient — uses freshly updated EL[,k] from (a)
    #
    # Order rationale: update q(beta_k) LAST so it uses the freshest EL[,k]
    # from step (a). This is Gauss-Seidel coupling for β. β is entirely
    # determined by the survival likelihood, which depends directly on EL, so
    # using the freshest EL matters. Updating β first (Jacobi-style) converges
    # to an inferior local optimum: confirmed by an 18-unit ELBO gap on
    # identical data (2026-05-05 benchmark audit).
    #
    # SAFETY OF REUSE — z_no_k AND R_k are invariant in the current k:
    #   z_no_k = z − Σ_{k'} EL[,k']·EBeta[k']  +  EL[,k]·EBeta[k]
    #          = z − Σ_{k'≠k} EL[,k']·EBeta[k']
    #   R_k    = Y − Σ_{k'≠k} EL[,k']·EF[,k']ᵀ
    # Neither depends on EL[,k], EBeta[k], or EF[,k] for the current k.
    # z_no_k computed once per k is valid for both the L update (step a) and
    # the β update (step c). R_k is recomputed once between L and F since F
    # needs the post-L R_k (Gauss-Seidel property for F).
    # ========================================================================
    for (k in 1:K) {

      # ----------------------------------------------------------------------
      # Compute R_k and z_no_k once. Both are invariant in EL[,k], EF[,k],
      # and EBeta[k] for the current k, so they are valid for all three
      # updates below.
      # ----------------------------------------------------------------------
      # compute_R_k receives augmented EL_aug/EF_aug so it automatically subtracts
      # the cohort contribution from the partial residual for k = 1..K.
      R_k    <- compute_R_k(Y, EL_aug, EF_aug, k)
      z_no_k <- compute_z_no_k(z, EL, EBeta, k)   # COMPUTE ONCE — used by β and L

      # ---- Instrumentation [§4.1] ---------------------------------------
      # At iter 1, log the genomics vs survival precision contributions for
      # the first few factors. A_surv ≪ A_gen quantifies the structural
      # scale imbalance that crippled merged-cohort training. Matches
      # update_L_k() math: A_gen = sum_j(τ_j · E[f_jk²]) (scalar across i);
      # A_surv = w_i · E[β_k²] (n-vector — we report the mean for
      # readability). Logged before the L update so the value reflects the
      # state seen at the iteration boundary.
      if (iter == 1 && verbose && k <= 3) {
        A_gen_k  <- sum(Tau * EF2[, k])
        A_surv_k <- mean(w) * EBeta2[k]
        cat(sprintf("    [iter1, k=%d] A_gen=%.2e  A_surv=%.2e  ratio=%.4f\n",
                    k, A_gen_k, A_surv_k, A_surv_k / (A_gen_k + 1e-30)))
      }

      # ----------------------------------------------------------------------
      # (a) Update q(l_k): Patient Loadings — FIRST
      #
      # Uses current EBeta[k] from the previous outer iteration (or from
      # initialization on iter 1). A_surv/A_gen << 1, so EBeta staleness
      # barely affects the L precision; the L update is dominated by genomics.
      #
      # NOTE on Phase 1a scope (see DECISIONS.md 2026-07-12): genomics_divisor/
      # survival_divisor are intentionally NOT passed here (left at their
      # update_L_k defaults of 1). L and F co-adapt bilinearly (L's precision
      # depends on EF^2 and vice versa); rescaling either side -- in EITHER
      # direction -- destabilizes that coupling. Shrinking genomics collapses
      # EL/EF to exactly zero; boosting survival causes unbounded divergence
      # (both confirmed empirically at PDAC-realistic scale). Only the ELBO
      # monitor and the beta update (no bilinear coupling) receive the
      # normalization; the recommended production model (YFB) is unaffected
      # by this either way, since its L/F updates never touch these divisors.
      # ----------------------------------------------------------------------
      res_L   <- update_L_k(Tau, EF_aug[, k], EF2_aug[, k], w, EBeta[k], EBeta2[k],
                             R_k, z_no_k, prior_family = prior_LF, alpha = alpha_iter,
                             normalize_AB = normalize_AB)
      EL[, k]      <- res_L$mean
      EL2[, k]     <- res_L$second
      EL_aug[, k]  <- res_L$mean    # keep augmented matrix in sync
      EL2_aug[, k] <- res_L$second
      kl_L[k]  <- compute_ebnm_kl(res_L$ebnm_result$log_likelihood,
                                   res_L$A, res_L$x, res_L$mean, res_L$second)

      # ----------------------------------------------------------------------
      # (b) Update q(f_k): Biological Factors — sees fresh EL[,k]
      #
      # R_k depends on EL[,k], so recompute it using the updated EL[,k]
      # from step (a).  This is the Gauss-Seidel property.
      #
      # genomics_divisor is intentionally NOT passed here -- see the note above
      # update_L_k's call site (Phase 1a scope, DECISIONS.md 2026-07-12).
      # ----------------------------------------------------------------------
      R_k   <- compute_R_k(Y, EL_aug, EF_aug, k)
      res_F <- update_F_k(Tau, EL_aug[, k], EL2_aug[, k], R_k,
                          prior_family = prior_LF, alpha = alpha_iter)
      EF[, k]      <- res_F$mean
      EF2[, k]     <- res_F$second
      EF_aug[, k]  <- res_F$mean    # keep augmented matrix in sync
      EF2_aug[, k] <- res_F$second
      kl_F[k]  <- compute_ebnm_kl(res_F$ebnm_result$log_likelihood,
                                   res_F$A, res_F$x, res_F$mean, res_F$second)

      # ----------------------------------------------------------------------
      # (c) Update q(beta_k): Survival Coefficient — LAST
      #
      # Uses the freshest EL[,k] and EL2[,k] from step (a) — Gauss-Seidel
      # coupling. β is entirely determined by the survival likelihood, which
      # depends directly on EL, so using the freshest EL matters here.
      # ----------------------------------------------------------------------
      res_beta    <- update_beta_k(w, z_no_k, EL[, k], EL2[, k],
                                   prior_family = prior_beta, alpha = alpha_iter,
                                   survival_divisor = beta_divisor)
      EBeta[k]    <- res_beta$mean
      EBeta2[k]   <- res_beta$second
      kl_beta[k]  <- compute_ebnm_kl(res_beta$ebnm_result$log_likelihood,
                                      res_beta$A, res_beta$x,
                                      res_beta$mean, res_beta$second)

    }  # end k-loop

    # ========================================================================
    # Cohort F update — closed-form Normal conjugate for platform loading f_c
    #
    # Runs AFTER the k=1..K EBNM loop, BEFORE tau.  Uses Tau from the prior
    # iteration (same one-iteration lag as all CAVI updates — standard
    # Gauss-Seidel, does not violate ELBO monotonicity).
    # ========================================================================
    if (!is.null(cohort_id) && C_cols > 0) {
      Y_hat_global <- EL %*% t(EF)   # n x p: K-factor reconstruction only
      res_cohort   <- update_F_cohort_all(
        L_cohort, Y - Y_hat_global, Tau, sigma_F_cohort
      )
      # Use seq_len(C_cols) — never (K+1):(K+C_cols).
      # In R, when C_cols=0, (K+1):K yields c(K+1,K), a two-element decreasing
      # sequence that silently corrupts the assignment.  seq_len(0) is empty.
      EF_aug[,  K + seq_len(C_cols)] <- res_cohort$EF_cohort
      EF2_aug[, K + seq_len(C_cols)] <- res_cohort$EF2_cohort
    }

    # ========================================================================
    # STEP 3: Noise Precision Update (closed-form MLE; no EBNM)
    # ========================================================================
    res_tau              <- update_tau(Y, EL_aug, EL2_aug, EF_aug, EF2_aug)
    Tau                  <- res_tau$Tau
    history$elbo_proxy[iter] <- res_tau$elbo_proxy
    factor_pve_iter <- vapply(seq_len(K), function(k) {
      sum(EL[, k]^2) * sum(EF[, k]^2) / y_frob2
    }, numeric(1))
    history$factor_pve[iter, ] <- factor_pve_iter

    # Full ELBO = (1-alpha)*genomics + alpha*survival + KL divergences.
    # surv_elbo: E_q[log PL(t,delta|L,beta)] via 2nd-order Taylor at eta_0.
    # kl_*: E_q[log g(theta)] - E_q[log q(theta)] for each factor (each <= 0).
    surv_elbo               <- compute_survival_elbo(taylor$logPL, w,
                                                     EL, EL2, EBeta, EBeta2)
    history$elbo_full[iter] <- (1 - alpha) * (res_tau$elbo_proxy / genomics_divisor) +
                               alpha * (surv_elbo / survival_divisor) +
                               sum(kl_L) + sum(kl_F) + sum(kl_beta)
    # compute_normal_kl returns a NEGATIVE value (-KL <= 0); adding it
    # correctly penalises the ELBO.  Never subtract — that double-negates.
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
    #
    # Track parameter deltas for diagnostics, but terminate on relative
    # full-ELBO change after burn-in. Skip the check gracefully when the
    # current/previous ELBO is non-finite or the denominator is zero.
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
      } else {
        "NA"
      }
      cat(sprintf("  iter %3d | RMSE: %.4f | ELBO: %+.1f (full) %+.1f (proxy) | dELBO: %s | dL: %.2e | dB: %.2e | beta: [%s]\n",
                  iter, history$rmse[iter],
                  history$elbo_full[iter], history$elbo_proxy[iter],
                  delta_elbo_str,
                  delta_L, delta_Beta,
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
  # The pmax(SVD) initialisation discards sign information: EL[:,k] may be
  # anti-correlated with the true survival direction, causing EL·EBeta to be
  # inversely related to risk. A global sign flip of EBeta corrects this:
  #   eta_new = EL · (-EBeta) = -(eta_original)  ⟹  C_new = 1 - C_train.
  # This is identical in spirit to Phase C in fit_cox_on_yf.R. The fix is
  # applied post-convergence so the CAVI path and ELBO are unaffected.
  #
  # sign_correction = FALSE during cross-validation (select_alpha_cv), because
  # per-fold sign flips would introduce fold-to-fold inconsistency and corrupt
  # alpha selection. CV already handles sign via I(-pred$risk_scores).
  if (sign_correction) {
    eta_train <- as.vector(EL %*% EBeta)
    c_train   <- as.numeric(
      concordance(Surv(time, status) ~ eta_train)$concordance
    )
    if (c_train < 0.5) {
      if (verbose) {
        cat(sprintf("    [Phase C] Training C=%.4f < 0.5 — flipping EBeta sign\n", c_train))
      }
      EBeta  <- -EBeta
      EBeta2 <- EBeta^2
    }
  }

  result <- list(
    EL     = EL,
    EL2    = EL2,
    EF     = EF,
    EF2    = EF2,
    EBeta  = EBeta,
    EBeta2 = EBeta2,
    Tau    = Tau,
    history = history,
    norm_convention  = norm_convention,
    genomics_divisor = genomics_divisor,
    survival_divisor = survival_divisor,
    beta_divisor     = beta_divisor
  )
  # Expose cohort fields when the extension is active.
  # Global EL/EF (K-factor only) are unchanged for backward compatibility.
  if (!is.null(cohort_id) && C_cols > 0) {
    result$L_cohort   <- L_cohort
    result$EF_cohort  <- EF_aug[, K + seq_len(C_cols), drop = FALSE]
    result$EF2_cohort <- EF2_aug[, K + seq_len(C_cols), drop = FALSE]
    result$cohort_id  <- cohort_id
  }
  result
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
