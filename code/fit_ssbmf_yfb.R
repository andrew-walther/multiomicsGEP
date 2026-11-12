# =============================================================================
# code/fit_ssbmf_yfb.R
# Purpose: Cluster B CAVI fitting loop — Cox-on-YF reformulation (eta = YFB)
#
# MODEL CHANGE VS. fit_ssbmf_lb.R (Cluster A — eta = LB):
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
# DERIVATIONS: derivations/SSBMF-YFB-model-derivations/
# =============================================================================

# Default: "real" — source() loads the function without side effects.
DATA_MODE <- "real"

real_Y      <- NULL
real_time   <- NULL
real_status <- NULL

# ==============================================================================
# PART 1 — LIBRARIES AND MODULE SOURCES
# ==============================================================================

source("code/profile_cox_intercept.R")
source("code/feature_groups_yfb.R")
library(survival)
library(ebnm)

# Cluster B update files (dedicated — do NOT source update_L.R or update_F.R)
source("code/update_L_surv_YFB.R")  # update_L_surv_YFB_k, update_L_surv_YFB_all
source("code/update_F_surv_YFB.R")  # update_F_surv_YFB_k, update_F_surv_YFB_all

# Shared update files (same interface, different inputs vs. Cluster A)
source("code/update_beta.R")        # compute_z_no_k, update_beta_k, update_beta_all
source("code/update_beta_cohort.R") # update_beta_cohort_k/_all, compute_pooled_beta -- beta_cohort_id
source("code/update_tau.R")         # compute_var_term, update_tau
source("code/compute_elbo.R")       # compute_ebnm_kl, compute_survival_elbo, compute_normal_kl
source("code/update_F_cohort.R")    # update_F_cohort_col, update_F_cohort_all
source("code/deflation_init.R")     # deflation_svd_init -- init_method="deflation"
# Exact sequential loading sweep and prior fits shared with the multimodal
# model; used only by canonicalize = TRUE (multimodal_yfb_update_F_mk etc.)
if (!exists("multimodal_yfb_update_F_mk")) source("code/multimodal_yfb_updates.R")

# compute_R_k is defined in update_L.R (Cluster A); re-source just that function
# by sourcing update_L.R here. The Cluster B update files use only compute_R_k,
# not update_L_k or update_L_all, so there is no crossover.
suppressPackageStartupMessages(
  source("code/update_L.R")           # provides compute_R_k (used by both _surv_YFB files)
)

# ------------------------------------------------------------------------------
#' Calculate Cox Score and Diagonal Hessian (Taylor Expansion)
#'
#' Identical to the helper in fit_ssbmf_lb.R. Duplicated here so fit_ssbmf_yfb.R
#' is a self-contained entry point (no hidden dependency on fit_ssbmf_lb.R).
#'
#' @param eta    n-vector: current linear predictor
#' @param time   n-vector: observed survival/censoring times
#' @param status n-vector: event indicator (1=event, 0=censored)
#' @param strata NULL (default) for a single pooled risk set, or an n-vector of
#'               stratum labels (e.g. study/cohort). When supplied, Breslow risk
#'               sets are formed *within* each stratum — the standard stratified
#'               Cox partial likelihood. The baseline hazard still cancels
#'               per-stratum, so no parametric h0 is introduced (Item 3,
#'               DECISIONS.md). strata=rep(1,n) reduces exactly to strata=NULL.
#' @return list(u, w, logPL)
# ------------------------------------------------------------------------------
calc_cox_taylor_yf <- function(eta, time, status, strata = NULL) {
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

    # Breslow tied-event handling: all rows sharing an event time must use the
    # SAME risk-set denominator and the SAME cumulative-hazard increment. The
    # naive per-row reverse-cumsum below only gives the correct denominator to
    # the FIRST row of each tied-time block (since it sums theta from that row
    # to the end, which for the earliest position in a contiguous tied block
    # already covers the whole block + everyone after); later rows in the same
    # block wrongly exclude the earlier tied rows. Fix: broadcast the
    # first-position value to every row in the block, and accumulate the
    # cumulative hazard once per unique time (not once per row).
    risk_sum_pos <- rev(cumsum(rev(theta)))       # per-row, pre-tie-fix
    first_idx    <- match(time_s, time_s)         # first sorted position sharing this row's time
    risk_sum     <- risk_sum_pos[first_idx]       # shared denominator within a tied-time block

    d_grp <- ave(status_s, time_s, FUN = sum)     # event count at this row's time, broadcast to the block
    is_first    <- !duplicated(time_s)            # one row per unique time, in increasing time order
    h_unique    <- d_grp[is_first] / risk_sum[is_first]
    H_unique    <- cumsum(h_unique)               # cumulative hazard, incremented once per unique time
    grp_id      <- cumsum(is_first)               # 1,2,3,... group index per row
    H           <- H_unique[grp_id]

    u_s <- status_s - theta * H
    w_s <- theta * H
    w_s[w_s < 1e-6] <- 1e-6
    u <- numeric(n); w <- numeric(n)
    u[ord] <- u_s;   w[ord] <- w_s
    logPL <- sum(status_s * eta_s) - sum(d_grp[is_first] * log(pmax(risk_sum[is_first], 1e-300)))
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

# ------------------------------------------------------------------------------
#' Rescale each program so that its training projection has unit SD
#'
#' Scale canonicalization for the Cox-on-YF model (Workstream C, step C2; the
#' single-modality port of `multimodal_yfb_canonicalize_factors()`). The
#' factorization Y ≈ L Fᵀ and the predictor η = (Y F) β are unchanged when
#' program k is rescaled by c_k > 0 as
#'   F_k → F_k / c_k,  L_k → c_k L_k,  β_k → c_k β_k
#' (second moments by c_k⁻² , c_k², c_k²). Choosing c_k = sd(Y F_k) fixes the
#' otherwise free scale of F_k, so the survival term in the F update cannot
#' feed back through the L update (problem 3 in the Workstream C plan). The
#' ebnm priors on L, F and β are re-estimated at every update and are scale
#' families, so they need no transformation here.
#'
#' @param Y      n x p data matrix used for the projection.
#' @param EL,EL2 n x K posterior mean / second moment of L (first K columns).
#' @param EF,EF2 p x K posterior mean / second moment of F (first K columns).
#' @param EBeta,EBeta2 K-vectors: posterior mean / second moment of β.
#' @param K      number of programs to rescale (the first K columns).
#' @param prior_F NULL or a length-K list of loading priors held across sweeps
#'               (multimodal format: family, pi, rate or variance, point_mass);
#'               rescaled with F so the prior describes the same loadings.
#' @param feature_group NULL or grouped feature indices, indicating nested priors.
#' @return list(EL, EL2, EF, EF2, EBeta, EBeta2, prior_F, scale): the rescaled moments
#'   and the K-vector of applied scales c_k (1 for a program with sd 0).
#' @examples
#' Y <- matrix(rnorm(20), 10, 2)
#' canonicalize_yf_factors(Y, matrix(1, 10, 1), matrix(1, 10, 1),
#'                         matrix(1, 2, 1), matrix(1, 2, 1), 0.5, 0.3, 1)$scale
#' @family F_surv_YFB_update
canonicalize_yf_factors <- function(Y, EL, EL2, EF, EF2, EBeta, EBeta2, K,
                                    prior_F = NULL, feature_group = NULL) {
  scale <- vapply(seq_len(K), function(k) stats::sd(as.vector(Y %*% EF[, k])), numeric(1))
  # A program whose projection has collapsed (sd < 1e-6, against 1 for a live
  # canonical program) is dead: rescaling it would inflate numerical residue
  # in F to unit scale and destabilize the next updates, so it is left as is.
  scale[!is.finite(scale) | scale < 1e-6] <- 1
  for (k in seq_len(K)) {
    EF[, k]  <- EF[, k] / scale[k];  EF2[, k] <- EF2[, k] / scale[k]^2
    EL[, k]  <- EL[, k] * scale[k];  EL2[, k] <- EL2[, k] * scale[k]^2
    EBeta[k] <- EBeta[k] * scale[k]; EBeta2[k] <- EBeta2[k] * scale[k]^2
    # A loading prior held across sweeps moves with F_k -> F_k / c: an
    # exponential or Laplace rate is multiplied by c, a Normal variance by 1/c^2
    pk <- prior_F[[k]]
    if (!is.null(feature_group)) {
      for (g in names(pk)) {
        pg <- pk[[g]]
        if (!isTRUE(pg$point_mass)) {
          if (identical(pg$family, "normal")) pg$variance <- pg$variance / scale[k]^2
          else pg$rate <- pg$rate * scale[k]
        }
        pk[[g]] <- pg
      }
      prior_F[[k]] <- pk
      next
    }
    if (!is.null(pk) && !isTRUE(pk$point_mass)) {
      if (identical(pk$family, "normal")) pk$variance <- pk$variance / scale[k]^2
      else pk$rate <- pk$rate * scale[k]
      prior_F[[k]] <- pk
    }
  }
  list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2, EBeta = EBeta, EBeta2 = EBeta2,
       prior_F = prior_F, scale = scale)
}

# ------------------------------------------------------------------------------
#' ELBO of a canonical Cox-on-YF fit at given posterior moments
#'
#' The objective monitored by fit_ssbmf_yfb(canonicalize = TRUE):
#' E_q[log p(Y | L, F, tau)] at the optimal tau, plus w_S times the Cox term
#' (partial log-likelihood re-expanded at E[eta], with Var(eta) from q(beta)
#' and q(F)), plus the stored -KL terms. Used by the pruning nullcheck.
#'
#' @param Y,time,status,strata Data as in fit_ssbmf_yfb().
#' @param EL,EL2,EF,EF2 Posterior moments (n x K, p x K).
#' @param EBeta,EBeta2 K-vectors.
#' @param neg_kl Sum of the stored -KL terms of the programs included.
#' @param w_S Survival weight.
#' @param survival_Y Design of the Cox projection (default Y). A fit with
#'   center_projection = TRUE passes the column-centered training matrix here,
#'   while the genomics term keeps Y.
#' @return Scalar ELBO.
#' @family F_surv_YFB_update
#' @examples
#' # Used internally by prune_yf_nullcheck().
yf_canonical_elbo <- function(Y, time, status, EL, EL2, EF, EF2, EBeta, EBeta2,
                              neg_kl, w_S = 1, strata = NULL, survival_Y = Y) {
  genomics <- update_tau(Y, EL, EL2, EF, EF2)$elbo_proxy
  if (ncol(EF) == 0) return(genomics + neg_kl)
  ZF <- survival_Y %*% EF
  VZ <- (survival_Y^2) %*% pmax(EF2 - EF^2, 0)
  tay <- calc_cox_taylor_yf(as.vector(ZF %*% EBeta), time, status, strata = strata)
  genomics + w_S * compute_survival_elbo(tay$logPL, tay$w, ZF, ZF^2 + VZ, EBeta, EBeta2) + neg_kl
}

#' ELBO change from removing each program (flashier-style nullcheck)
#'
#' @param fit A fit_ssbmf_yfb(canonicalize = TRUE) result (needs $neg_kl).
#' @param Y,time,status,strata Training data.
#' @param survival_Y Design of the Cox projection (default Y); see yf_canonical_elbo().
#' @return K-vector: ELBO(without program k) - ELBO(full). A nonnegative value
#'   means the program does not improve the objective.
#' @family F_surv_YFB_update
#' @examples
#' @param protect_survival Protect coefficients satisfying |Eβ| >= 1.96 SD(β) and SD(β)>0.
#' # delta <- prune_yf_nullcheck(fit, Y, time, status)
prune_yf_nullcheck <- function(fit, Y, time, status, strata = NULL,
                                protect_survival = FALSE, survival_Y = Y) {
  K <- ncol(fit$EF)
  nk <- fit$neg_kl
  full <- yf_canonical_elbo(Y, time, status, fit$EL, fit$EL2, fit$EF, fit$EF2,
                            fit$EBeta, fit$EBeta2, sum(nk), fit$w_S, strata, survival_Y)
  vapply(seq_len(K), function(k) {
    # A dead zero-variance coefficient must remain removable.
    beta_sd <- sqrt(pmax(fit$EBeta2[k] - fit$EBeta[k]^2, 0))
    if (protect_survival && beta_sd > 0 &&
        abs(fit$EBeta[k]) >= 1.96 * beta_sd) return(-Inf)
    keep <- setdiff(seq_len(K), k)
    yf_canonical_elbo(Y, time, status, fit$EL[, keep, drop = FALSE],
                      fit$EL2[, keep, drop = FALSE], fit$EF[, keep, drop = FALSE],
                      fit$EF2[, keep, drop = FALSE], fit$EBeta[keep], fit$EBeta2[keep],
                      sum(nk[keep]), fit$w_S, strata, survival_Y) - full
  }, numeric(1))
}

# ==============================================================================
# PART 2 — fit_ssbmf_yfb()
# ==============================================================================

#' Subset a canonical single-modality posterior without resetting its state
#' @param state Posterior state returned by prune_preserve_posterior=TRUE.
#' @param keep Unique retained program indices, in the desired order.
#' @return State with score/loading moments, beta moments and loading priors
#'   subset together; feature precisions and the shared beta prior unchanged.
#' @examples
#' # yf_subset_posterior_state(fit$posterior_state,c(1L,3L))
#' @family supervised_yfb
#' @seealso [fit_ssbmf_yfb()]
yf_subset_posterior_state <- function(state,keep) {
  K <- ncol(state$EL)
  if(!is.numeric(keep) || !length(keep) || anyNA(keep) || any(!is.finite(keep)) ||
    any(keep!=as.integer(keep)) || any(keep<1L | keep>K) || anyDuplicated(keep))
    stop("keep must contain unique retained program indices.")
  out <- state
  for(nm in c("EL","EL2","EF","EF2"))out[[nm]] <- state[[nm]][,keep,drop=FALSE]
  for(nm in c("EBeta","EBeta2","prior_F"))out[[nm]] <- state[[nm]][keep]
  out
}

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
#' @param Y        n x p genomics data matrix (should be column-centered, unless
#'                 center_projection = TRUE, which takes nonnegative uncentered input)
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
#'                 (see DECISIONS.md and fit_ssbmf_lb.R). Under YFB, L and F are
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
#' @param verbose  Logical: print iteration logs? (default TRUE)
#' @param w_S      NULL (default) or a nonnegative scalar survival weight for the
#'                 F update. NULL keeps the alpha_F-weighted update. When supplied,
#'                 alpha_F is ignored and q(f_k) combines both likelihoods
#'                 unweighted, with the survival term multiplied by w_S:
#'                 A_F = A_gen + w_S A_surv, B_F = B_gen + w_S B_surv.
#'                 w_S = 1 is the plain (equal-weight) Bayesian posterior.
#' @param canonicalize logical (default FALSE). When TRUE (Workstream C, step
#'                 C2), the survival predictor is the raw projection
#'                 ZF = Y EF with every program rescaled to sd(Y EF_k) = 1 on
#'                 the training data after initialization and at the end of
#'                 every iteration (`canonicalize_yf_factors()`), instead of
#'                 Y EF / ||EF_k||. alpha_F is then not used: every survival
#'                 term (F update, beta update, monitored ELBO) carries the
#'                 same weight w_S (NULL means w_S = 1, the plain Bayesian
#'                 posterior), and the monitored ELBO is
#'                 E[log p(Y)] + w_S * survival term + KL terms, the objective
#'                 the updates climb. The F update is then the exact
#'                 sequential (gene-by-gene) coordinate update of
#'                 `multimodal_yfb_update_F_mk()`, which removes the other
#'                 genes' contribution to program k's projection from each
#'                 gene's survival term (the vector update used otherwise omits
#'                 it), under a loading prior held fixed during the sweep and
#'                 refit afterwards; prior_LF must be point_exponential,
#'                 point_laplace or normal. The beta update and the survival
#'                 ELBO include the projection variance Var(Y F_k) from q(F),
#'                 and the survival ELBO is evaluated at the end of each
#'                 iteration. `alpha`, `norm_convention` and
#'                 `boost_beta` are ignored. `EF_norms` is returned as 1s, so
#'                 predict_ssbmf_yfb() computes Y_test EF beta. Not supported
#'                 with alpha_F > 0, N_burnin > 0, cox_warmstart,
#'                 alpha_schedule or beta_cohort_id.
#' @param init_supervised logical (default FALSE; canonicalize = TRUE only).
#'                 Workstream C step C3: program K starts from the Cox score
#'                 direction of the null model, F_K = Y' u0 with u0 the
#'                 martingale residuals at eta = 0, and L_K = R F_K / ||F_K||^2
#'                 with R the data minus programs 1..K-1. Its scores then
#'                 correlate with survival from the first sweep, so beta_K is
#'                 estimated away from 0 before F moves.
#' @param beta_prior_mode "per_program" (default): each beta_k update fits its
#'                 own ebnm prior to its single pseudo-observation.
#'                 "shared" (canonicalize = TRUE only; step C4): one prior g on
#'                 beta_1..beta_K, held fixed during each sweep and refit by
#'                 marginal likelihood from the K pseudo-observations after the
#'                 sweep (the pattern used for the loadings).
#' @param beta_prior_init NULL or an ebnm prior object: the starting shared prior.
#'                 NULL gives N(0, 1) for "normal", 0.5 delta_0 + 0.5 N(0, 1) for
#'                 "point_normal" and 0.5 delta_0 + 0.5 Laplace(scale 1) for
#'                 "point_laplace" (beta is a log hazard ratio per SD of score).
#' @param beta_warmup integer (default 0): iterations during which the shared
#'                 prior is held at beta_prior_init before it is refit.
#' @param beta_prior_fix logical (default FALSE): never refit the shared prior.
#' @param beta_slab_sd NULL or positive: for prior_beta = "point_normal", fix the
#'                 slab SD at this value and estimate only the mixing weight.
#' @param prune logical (default FALSE; canonicalize = TRUE only). After
#'                 convergence, remove the program whose removal raises the
#'                 ELBO most if any removal does not lower it
#'                 (`prune_yf_nullcheck()`), refit from the remaining programs
#'                 (warm start, shared beta prior carried over, no warm-up),
#'                 and repeat. The result has fewer than K columns and
#'                 `$pruning` records each removal.
#' @param posterior_init NULL (default) or a full canonical posterior state,
#'   preserving EL/EL2, EF/EF2, beta moments, Tau and persistent loading/shared
#'   beta priors. Requires no cohort offset and N_frozen=0; means need not be
#'   supplied separately. Initial uncertainty and beta affect the first update.
#' @param loading_prior_update "historical" (default) or "fixed_q_em". Opt-in
#'   fixed-posterior point-Laplace F-prior maximization with KL to the returned
#'   prior. Requires canonical per-program Normal beta, no cohort/strata effects
#'   and no pruning; no loading sweep, Cox approximation or stopping rule changes.
#' @param profile_cox_intercept Opt-in profiled working Cox intercept. Weighted
#'   centering applies only to survival updates; reconstruction, raw prediction
#'   and the original ELBO monitor/rule remain unchanged. Requires canonical
#'   per-program Normal beta, no cohort/strata effects and no pruning.
#' @param prune_preserve_posterior Logical, default FALSE. Carry that full state
#'   into each pruning refit rather than reinitializing from L/F means. Also
#'   returns posterior_state for continuation. This changes only initialization;
#'   the prior families, coordinate updates, objective and pruning rule are fixed.
#' @param prune_protect_survival Logical, default FALSE. Protect programs with
#'   |E[β_k]| >= 1.96 sqrt(E[β_k²] - E[β_k]²) and SD(β_k)>0 from the ELBO nullcheck.
#' @param converge_on Stopping criterion: default "elbo" preserves the existing
#'   rule; "outputs" requires relative changes of LF' and YFβ both < 1e-4.
#'   Explicitly passing either value records both final criteria; omitting the
#'   argument preserves the default fit object bit for bit.
#' @param trace_supervision logical (default FALSE). When TRUE, `history$prior_trace`
#'   also records each non-cohort sweep's beta prior before/after refitting and
#'   its raw coefficient pseudo-observations/precisions. `history$trace`
#'                 records, per iteration (rows) and program k (columns):
#'                 E[beta_k], E[beta_k^2]; survival's share of the F-update
#'                 precision, w_surv A_surv / A_F (median and 90th percentile over
#'                 features), and of the signal, |w_surv B_surv| /
#'                 (|w_gen B_gen| + |w_surv B_surv|); the unweighted precision
#'                 share A_surv / (A_gen + A_surv) (median); sd(Y EF_k), ||EF_k||
#'                 and ||EL_k||. Here (w_gen, w_surv) = (1 - alpha_F, alpha_F), or
#'                 (1, w_S) when w_S is supplied. The fit itself is unchanged.
#'
#' @param center_projection logical (default FALSE; canonicalize = TRUE only).
#'                 TRUE centers Y inside the Cox projection only: the survival
#'                 predictor is (Y - 1 ybar') E[F] beta with ybar the training
#'                 column means, while the reconstruction Y ~ L F' keeps the
#'                 uncentered Y. For nonnegative input (X/sd within cohort) with
#'                 point-exponential L and F. The Cox partial likelihood is
#'                 unchanged by the centering; it removes the common component
#'                 ybar' F beta from the second-order surrogate, which would
#'                 otherwise add sum_i W_i ybar_j^2 to every gene's survival
#'                 precision (derivation, Section "Centering Only the Cox
#'                 Projection"). Every survival term (beta and F updates,
#'                 monitored ELBO, pruning, output convergence) uses the centered
#'                 matrix. The fit returns `projection_center` = ybar; pass it to
#'                 predict_ssbmf_yfb(). With strata_id the centering is within
#'                 strata (each stratum's partial likelihood is invariant to its
#'                 own constant): y_i - ybar_s(i), and `projection_center` is a
#'                 strata x p matrix with the stratum labels as row names. Not
#'                 combined with profile_cox_intercept. The monitored ELBO uses a different surrogate
#'                 than with the option off, so ELBO values are not comparable
#'                 between the two settings.
#' @param feature_group NULL or one label per column of Y. A single group uses
#'                 the unchanged numerical path. Multiple groups fit one F prior
#'                 per group and program while keeping joint Z = sum_m Y_m F_m.
#'                 Returns prior_F[[group]][[program]] and feature_group.
#'                 Cohort offsets are unsupported. Multiple groups also reject
#'                 pruning, posterior restarts, frozen F and fixed_q_em.
#' @param beta_cohort_id NULL (default) or an n-vector of cohort labels for
#'                 COHORT-SPECIFIC SURVIVAL COEFFICIENTS beta_k^{c(i)}
#'                 (code/update_beta_cohort.R). A third, independent grouping
#'                 variable, distinct from `cohort_id` (genomics offsets in
#'                 L/F, beta_cohort=0 by construction there) and `strata_id`
#'                 (cohort-specific baseline hazard, beta still shared). Only
#'                 `beta_cohort_id` lets a factor be prognostic in one
#'                 cohort and not another; all three can be combined (in
#'                 practice they are often the same underlying grouping
#'                 vector, e.g. dataset_labels, but they are different model
#'                 components). The C cohort coefficients per factor are
#'                 partially pooled through one shared EBNM prior (one
#'                 vectorized `ebnm()` call per factor across cohorts, not C
#'                 independent calls). `EBeta`/`EBeta2` become K x C
#'                 matrices when supplied. NULL (default) reproduces
#'                 today's shared-beta fit bit-for-bit
#'                 (tests/test_update_beta_cohort.R). Not supported with
#'                 `alpha_F > 0` (loud error).
#' @return Named list:
#'   $EL, $EL2, $EF, $EF2, $EBeta, $EBeta2, $Tau, $history. When
#'   `beta_cohort_id` is supplied: `$EBeta`/`$EBeta2` are K x C matrices,
#'   plus `$EBeta_pooled` (K-vector, for external prediction -- see
#'   `predict_ssbmf_yfb()`'s `cohort_id_test` argument), `$beta_cohort_id`,
#'   `$beta_cohort_levels`.
fit_ssbmf_yfb <- function(Y, time, status,
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
                           verbose      = TRUE,
                           cohort_id       = NULL,
                           sigma_F_cohort  = 1.0,
                           strata_id       = NULL,
                           beta_cohort_id  = NULL,
                           w_S             = NULL,
                           canonicalize    = FALSE,
                           init_supervised = FALSE,
                           beta_prior_mode = c("per_program", "shared"),
                           beta_prior_init = NULL,
                           beta_warmup     = 0L,
                           beta_prior_fix  = FALSE,
                           beta_slab_sd    = NULL,
                           prune           = FALSE,
                           prune_protect_survival = FALSE,
                           converge_on = c("elbo", "outputs"),
                           trace_supervision = FALSE,
                           posterior_init = NULL,
                           prune_preserve_posterior = FALSE,
                           profile_cox_intercept = FALSE,
                           loading_prior_update = "historical",
                           survival_block_weight = FALSE,
                           beta_posterior = c("meanfield", "joint"),
                           center_projection = FALSE,
                           feature_group = NULL) {

  # ---- S-13 options (derivations/SSBMF-YFB-model-derivations/components/survival_block_weight/) ----------------
  # survival_block_weight = TRUE: w_S multiplies the Cox term AND KL(q(beta) || g_beta), so the
  # beta update is the untempered empirical-Bayes posterior (Cox weight 1) while the F update
  # still sees survival at weight w_S. At w_S = 1 the fit is unchanged.
  # beta_posterior = "joint": q(beta) = N_K(m, Sigma) with per-program Normal prior variances
  # set by type-II ML (update_beta_joint_ard()); F updates use E[beta_k beta_l].
  beta_posterior <- match.arg(beta_posterior)
  if (!is.logical(survival_block_weight) || length(survival_block_weight) != 1L ||
      is.na(survival_block_weight)) stop("survival_block_weight must be TRUE or FALSE.")
  if ((survival_block_weight || beta_posterior == "joint") && !canonicalize)
    stop("survival_block_weight and beta_posterior = 'joint' require canonicalize = TRUE.")
  if (beta_posterior == "joint" && (prior_beta != "normal" || beta_prior_mode[1] != "per_program" ||
      !is.null(beta_cohort_id) || profile_cox_intercept || prune))
    stop("beta_posterior = 'joint' requires per-program Normal beta, no cohort beta, no profiled ",
         "intercept and no pruning.")
  joint_beta <- beta_posterior == "joint"

  if(!is.logical(profile_cox_intercept) || length(profile_cox_intercept)!=1L ||
     is.na(profile_cox_intercept))stop("profile_cox_intercept must be TRUE or FALSE.")
  if(profile_cox_intercept && (!canonicalize || beta_prior_mode[1]!="per_program" ||
     prior_beta!="normal" || !is.null(cohort_id) || !is.null(strata_id) ||
     !is.null(beta_cohort_id) || prune))
    stop("Profiled Cox curvature requires canonical per-program Normal beta, no cohort/strata effects, and no pruning.")
  if (!identical(loading_prior_update, "historical") &&
      !identical(loading_prior_update, "fixed_q_em"))
    stop("loading_prior_update must be historical or fixed_q_em.")
  if (identical(loading_prior_update, "fixed_q_em") &&
      (!canonicalize || prior_LF != "point_laplace" || beta_prior_mode[1] != "per_program" ||
       prior_beta != "normal" || !is.null(cohort_id) || !is.null(strata_id) ||
       !is.null(beta_cohort_id) || prune))
    stop("fixed_q_em requires canonical point-Laplace F, per-program Normal beta, no cohort/strata effects or pruning.")
  if (!is.logical(center_projection) || length(center_projection) != 1L ||
      is.na(center_projection)) stop("center_projection must be TRUE or FALSE.")
  if (center_projection && !canonicalize)
    stop("center_projection = TRUE requires canonicalize = TRUE.")
  if (center_projection && profile_cox_intercept)
    stop("center_projection and profile_cox_intercept both center the survival design; ",
         "use one of them.")
  record_output_convergence <- !missing(converge_on)
  call_args <- mget(names(formals(sys.function())), envir = environment())
  n <- nrow(Y); p <- ncol(Y)
  groups <- yf_feature_groups(feature_group, p)
  F_group_priors <- if (!is.null(groups)) vector("list", K) else NULL
  if (!is.null(feature_group) && !is.null(cohort_id))
    stop("feature_group does not support cohort_id offsets.")
  if (!is.null(groups) && (prune || prune_preserve_posterior || !is.null(posterior_init) ||
                           N_frozen > 0 || loading_prior_update != "historical"))
    stop("Multiple feature_group labels do not support prune, prune_preserve_posterior, ",
         "posterior_init, N_frozen or loading_prior_update other than historical.")
  # Design of the Cox projection: Y itself unless center_projection = TRUE, in
  # which case Y - 1 ybar' with the training column means (fixed for the fit).
  # The reconstruction always uses Y. With the option off this is the same
  # object as Y, so the default path is unchanged. With strata_id the Cox risk
  # sets are formed within strata, and each stratum's partial likelihood is
  # invariant to its own constant, so Y is centered within strata: row i minus
  # its stratum's training means (derivation, Eq. center-strata).
  projection_center <- NULL
  projection_Y <- Y
  if (center_projection && is.null(strata_id)) {
    projection_center <- colMeans(Y)
    projection_Y <- sweep(Y, 2, projection_center)
  } else if (center_projection) {
    if (length(strata_id) != n || anyNA(strata_id))
      stop("strata_id must have length nrow(Y) (", n, ") and no NA.")
    stratum <- factor(strata_id)
    projection_center <- rowsum(Y, stratum) / as.vector(table(stratum))
    rownames(projection_center) <- levels(stratum)
    projection_Y <- Y - projection_center[as.integer(stratum), , drop = FALSE]
  }
  if(!is.logical(prune_preserve_posterior) || length(prune_preserve_posterior)!=1L ||
     is.na(prune_preserve_posterior))stop("prune_preserve_posterior must be TRUE or FALSE.")
  if(prune_preserve_posterior || !is.null(posterior_init)) {
    if(!canonicalize || !is.null(cohort_id) || N_frozen>0)
      stop("Posterior-state starts require canonicalize=TRUE, no cohort offset and N_frozen=0.")
  }
  if(!is.null(posterior_init)) {
    required <- c("EL","EL2","EF","EF2","EBeta","EBeta2","Tau","prior_F","g_beta")
    if(!is.list(posterior_init) || !all(required %in% names(posterior_init)))
      stop("posterior_init must contain every posterior moment, precision and persistent prior.")
    for(nm in c("EL","EL2","EF","EF2")) {
      nr <- if(nm %in% c("EL","EL2"))n else p
      x <- posterior_init[[nm]]
      if(!is.matrix(x) || !is.numeric(x) || !identical(dim(x),c(as.integer(nr),as.integer(K))) ||
        any(!is.finite(x)))stop("posterior_init has invalid dimensions or values for ",nm)
    }
    if(length(posterior_init$EBeta)!=K || length(posterior_init$EBeta2)!=K ||
      any(!is.finite(c(posterior_init$EBeta,posterior_init$EBeta2))) ||
      length(posterior_init$Tau)!=p || any(!is.finite(posterior_init$Tau)) ||
      any(posterior_init$Tau<=0) || length(posterior_init$prior_F)!=K)
      stop("posterior_init beta moments, precisions and priors must match the retained programs.")
    multimodal_yfb_validate_moments(posterior_init$EL,posterior_init$EL2,"posterior_init EL")
    multimodal_yfb_validate_moments(posterior_init$EF,posterior_init$EF2,"posterior_init EF")
    multimodal_yfb_validate_moments(posterior_init$EBeta,posterior_init$EBeta2,"posterior_init beta")
    EL_init <- posterior_init$EL;EF_init <- posterior_init$EF
    if(beta_prior_mode[1]=="shared")beta_prior_init <- posterior_init$g_beta
  }
  converge_on <- match.arg(converge_on)
  if (record_output_convergence && !exists("multimodal_yfb_relative_change", mode = "function"))
    source("code/fit_multimodal_yfb.R")

  if (!is.null(w_S) && (!is.numeric(w_S) || length(w_S) != 1 || !is.finite(w_S) || w_S < 0))
    stop("w_S must be NULL or a finite nonnegative scalar.")
  if (canonicalize) {
    if (alpha_F > 0)
      stop("alpha_F is not used with canonicalize = TRUE; set the survival weight with w_S.")
    if (N_burnin > 0 || cox_warmstart || !is.null(alpha_schedule) || !is.null(beta_cohort_id))
      stop("canonicalize = TRUE does not support N_burnin, cox_warmstart, alpha_schedule ",
           "or beta_cohort_id.")
    if (is.null(w_S)) w_S <- 1
    if (!prior_LF %in% c("point_exponential", "point_laplace", "normal"))
      stop("canonicalize = TRUE supports prior_LF point_exponential, point_laplace or normal.")
    prior_F_state <- if(is.null(posterior_init)) vector("list", K) else posterior_init$prior_F
    # Loading priors are the persistent factor-specific priors. The L prior is
    # re-estimated from its pseudo-observations at every coordinate update.
  }
  beta_prior_mode <- match.arg(beta_prior_mode)
  if (!canonicalize && (init_supervised || beta_prior_mode == "shared" || prune))
    stop("init_supervised, beta_prior_mode = 'shared' and prune require canonicalize = TRUE.")
  if (beta_prior_mode == "shared") {
    if (!prior_beta %in% c("normal", "point_normal", "point_laplace"))
      stop("beta_prior_mode = 'shared' supports prior_beta normal, point_normal or point_laplace.")
    g_beta <- if (!is.null(beta_prior_init)) beta_prior_init else switch(prior_beta,
      normal        = ashr::normalmix(1, 0, 1),
      point_normal  = ashr::normalmix(c(0.5, 0.5), c(0, 0),
                                      c(0, if (is.null(beta_slab_sd)) 1 else beta_slab_sd)),
      point_laplace = ebnm::laplacemix(c(0.5, 0.5), c(0, 0), c(0, 1)))
  }
  # Weights applied to the genomics and survival terms of the F update
  # (used by the trace; the update itself receives alpha_F or w_S directly).
  F_w_gen  <- if (is.null(w_S)) 1 - alpha_F else 1
  F_w_surv <- if (is.null(w_S)) alpha_F else w_S

  # ---- Cohort-specific survival coefficients (Stage 2, DECISIONS.md 2026-09-04) --
  # beta_cohort_id is a THIRD, independent grouping variable -- distinct from
  # cohort_id (genomics offsets appended to L/F) and strata_id (cohort-specific
  # baseline hazard, beta still shared). Only beta_cohort_id lets a factor be
  # prognostic in one cohort and not another. See code/update_beta_cohort.R's
  # file header for the full three-way distinction. NULL (default) preserves
  # today's shared-beta behavior exactly -- verified bit-for-bit in
  # tests/test_update_beta_cohort.R.
  use_beta_cohort <- !is.null(beta_cohort_id)
  if (use_beta_cohort) {
    if (length(beta_cohort_id) != n)
      stop("beta_cohort_id must be NULL or have length nrow(Y) (", n, ").")
    if (anyNA(beta_cohort_id))
      stop("beta_cohort_id must not contain NA.")
    beta_cohort_id  <- factor(beta_cohort_id)
    beta_cohort_idx <- as.integer(beta_cohort_id)   # n-vector in 1..C
    C_beta          <- nlevels(beta_cohort_id)
    beta_cohort_levels <- levels(beta_cohort_id)
    if (C_beta < 2)
      stop("beta_cohort_id must have at least 2 distinct levels (got ", C_beta, ").")
    if (alpha_F > 0 || (!is.null(w_S) && w_S > 0))
      stop("beta_cohort_id is not supported with alpha_F > 0 or w_S > 0: update_F_surv_YFB_k()'s ",
           "EBeta_k/EBeta2_k arguments are scalars, and `EBeta[k]` on a K x C_beta ",
           "cohort-beta matrix would silently index the wrong (flattened) element ",
           "rather than error. alpha_F > 0 is not the default and has no cohort-beta ",
           "decomposition implemented -- failing loudly here rather than silently ",
           "scoring the F update with the wrong beta value.")
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
  # genomics -- see fit_ssbmf_lb.R and DECISIONS.md 2026-07-12 for rationale.
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
  # Identical logic to fit_ssbmf_lb.R; see that file for full commentary.
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
  init_nonneg <- prior_LF %in% c("point_exponential", "gamma", "exponential")
  if (init_method == "svd") {
    svd_init <- svd(Y_for_svd, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    # Use abs() rather than pmax(..., 0): pmax zeros out entire columns when the
    # SVD vector points all-negative, giving degenerate ZF = Y·EF columns of all
    # zeros. Under point_exponential prior (non-negative), abs() is equally valid
    # and prevents the degenerate column pathology.
    EL  <- svd_init$u %*% diag(d_k, K, K)
    EF  <- svd_init$v %*% diag(d_k, K, K)
    # abs() only for a nonnegative prior: with a signed prior (point_laplace,
    # normal, point_normal) it discards the signs that the first L update
    # needs and the programs collapse to zero (2026-10-15, Workstream C).
    if (init_nonneg) { EL <- abs(EL); EF <- abs(EF) }
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
    EL <- if (init_nonneg) abs(defl$EL) else defl$EL
    EF <- if (init_nonneg) abs(defl$EF) else defl$EF
  } else if (init_method == "custom") {
    if (is.null(EL_init) || is.null(EF_init))
      stop("init_method='custom' requires both EL_init (n x K) and EF_init (p x K).")
    EL <- EL_init
    EF <- EF_init
  } else {
    stop(sprintf("Unknown init_method: '%s'. Use 'svd', 'random', 'deflation', or 'custom'.", init_method))
  }

  if (init_supervised && init_method != "custom") {
    # C3: program K starts along the null-model Cox score direction Y'u0
    u0 <- calc_cox_taylor_yf(rep(0, n), time, status, strata = strata_id)$u
    d  <- as.vector(crossprod(Y_for_svd, u0))
    EF[, K] <- d
    R_K <- Y_for_svd - EL[, -K, drop = FALSE] %*% t(EF[, -K, drop = FALSE])
    EL[, K] <- as.vector(R_K %*% d) / sum(d^2)
  }

  EL2 <- if(is.null(posterior_init)) EL^2 else posterior_init$EL2
  EF2 <- if(is.null(posterior_init)) EF^2 else posterior_init$EF2

  if (canonicalize && is.null(posterior_init)) {
    # Start every program at sd(Y EF_k) = 1 so beta starts on its final scale.
    # Standardize the scores the Cox model uses (projection_Y is Y unless
    # center_projection; centering changes the SD only through strata)
    cz <- canonicalize_yf_factors(projection_Y, EL, EL2, EF, EF2, rep(0, K), rep(0, K), K)
    EL <- cz$EL; EL2 <- cz$EL2; EF <- cz$EF; EF2 <- cz$EF2
  }

  # Augment EL/EF matrices with cohort indicator columns (same as fit_ssbmf_lb.R).
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
  if (!is.null(posterior_init)) {
    EBeta <- posterior_init$EBeta
  } else if (cox_warmstart) {
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
  EBeta2 <- if(is.null(posterior_init)) EBeta^2 else posterior_init$EBeta2

  if (verbose && !is.null(posterior_init)) {
    cat(sprintf("    [init] Full posterior warm start; retained EBeta range: [%.3e, %.3e]\n",
                min(EBeta), max(EBeta)))
  } else if (verbose && !cox_warmstart) {
    cat(sprintf("    [init] EBeta initialized to 0 (cox_warmstart=FALSE)\n"))
  }

  # Expand EBeta/EBeta2 from K-vectors to K x C_beta matrices when
  # beta_cohort_id is supplied -- every cohort column starts from the same
  # (possibly warm-started) value. NULL preserves the K-vector path exactly.
  if (use_beta_cohort) {
    EBeta  <- matrix(EBeta,  K, C_beta)
    EBeta2 <- matrix(EBeta2, K, C_beta)
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
      if (use_beta_cohort) {
        eta_b    <- rowSums(ZF_b * t(EBeta)[beta_cohort_idx, , drop = FALSE])
        taylor_b <- calc_cox_taylor_yf(eta_b, time, status, strata = strata_id)
        z_b      <- eta_b + taylor_b$u / taylor_b$w
        w_b      <- taylor_b$w
        res_b    <- update_beta_cohort_all(w_b, z_b, ZF_b, EBeta, beta_cohort_idx, C_beta,
                                            prior_family = prior_beta, alpha = alpha,
                                            survival_divisor = beta_divisor)
        EBeta  <- res_b$EBeta
        EBeta2 <- res_b$EBeta2
      } else {
        eta_b    <- as.vector(ZF_b %*% EBeta)
        taylor_b <- calc_cox_taylor_yf(eta_b, time, status, strata = strata_id)
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
    }
    if (verbose) {
      cat(sprintf("    [burnin] After %d β-only iters, EBeta range: [%.3e, %.3e]\n",
                  N_burnin, min(EBeta), max(EBeta)))
    }
  }

  Tau <- if(is.null(posterior_init)) 1.0 / pmax(apply(Y, 2, var), 1e-8) else posterior_init$Tau
  # Record the actual pre-sweep state for exact matched-start comparisons.
  # Keep default fit objects unchanged when no full posterior start was supplied.
  if (!is.null(posterior_init)) {
    initial_state <- list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2,
      EBeta = EBeta, EBeta2 = EBeta2, Tau = Tau,
      prior_F = prior_F_state, g_beta = if (beta_prior_mode == "shared") g_beta else NULL)
  }
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
  if (trace_supervision) {
    trace_names <- c("beta_mean", "beta_second", "shareA_median", "shareA_p90",
                     "shareB_median", "shareB_p90", "shareA_unweighted_median",
                     "sd_YF", "norm_F", "norm_L")
    if (!use_beta_cohort && beta_prior_mode == "shared") history$prior_trace <- vector("list", max_iter)
    history$trace <- setNames(
      lapply(trace_names, function(nm) matrix(NA_real_, max_iter, K)), trace_names)
  }

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
  # Joint q(beta) state (S-13): covariance, prior variances (warm start), E[beta beta']
  EBeta_cov <- matrix(0, K, K); beta_sigma2 <- NULL; Mbeta <- matrix(0, K, K)
  beta_lik_weight <- if (canonicalize && survival_block_weight) 1 else w_S
  for (iter in 1:max_iter) {
    if (trace_supervision && !use_beta_cohort && beta_prior_mode == "shared") {
      trace_prior_before <- g_beta
      trace_x <- trace_s <- trace_A <- numeric(K)
    }

    if (record_output_convergence) {
      old_reconstruction <- EL %*% t(EF)
      old_eta_outputs <- as.vector(projection_Y %*% EF %*% EBeta)
    }
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
    beta_x  <- numeric(K); beta_s <- numeric(K)

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
    EF_norms  <- if (canonicalize) rep(1, K)             # canonical: ZF = Y EF
                 else sqrt(colSums(EF_global^2) + 1e-10)   # K-vector
    EF_norm   <- sweep(EF_global, 2, EF_norms, "/")    # p x K, unit-norm columns
    ZF        <- projection_Y %*% EF_norm              # n x K: normalized projection scores
    eta    <- if (use_beta_cohort) rowSums(ZF * t(EBeta)[beta_cohort_idx, , drop = FALSE])
              else as.vector(ZF %*% EBeta)           # eta = ZF * beta_tilde (or beta_tilde^{c(i)})
    taylor <- calc_cox_taylor_yf(eta, time, status, strata = strata_id)
    z      <- eta + taylor$u / taylor$w    # working response z_i
    w      <- taylor$w                     # W_{ii} diagonal Hessian weights
    survival_Y <- if(profile_cox_intercept)profile_cox_center(Y,w) else projection_Y
    if(profile_cox_intercept)
      z <- profile_cox_center(eta,w) + taylor$u/w
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
      if (canonicalize) {
        # Current projections (earlier programs already moved this sweep) and
        # their variance under q(F): Var(Y F_k) = Y^2 Var(F_k)
        ZF   <- survival_Y %*% EF_aug[, 1:K, drop = FALSE]
        VZ_k <- as.vector(survival_Y^2 %*% pmax(EF2_aug[, k] - EF_aug[, k]^2, 0))
      }
      # Cluster B: z_no_k computed w.r.t. ZF (not EL as in Cluster A)
      z_no_k <- if (use_beta_cohort) compute_z_no_k_cohort(z, ZF, EBeta, beta_cohort_idx, k)
                else compute_z_no_k(z, ZF, EBeta, k)

      # ---- (a) Update q(beta_k): Survival Coefficient ----
      # Covariate is ZF[,k] = (Y·EF)[,k] — observed projection (not latent).
      # Pass ZF[,k]^2 as EL2_k (the second moment argument): ZF is observed, so
      # its "posterior second moment" = squared value (no variance term).
      # A_beta = sum(w * ZF_k^2) is non-zero from SVD init regardless of EBeta,
      # breaking the chicken-and-egg that plagued Cluster A's L update.
      if (use_beta_cohort) {
        # One vectorized ebnm() call across the C_beta cohort columns of
        # factor k (code/update_beta_cohort.R) -- partial pooling, not C
        # independent update_beta_k() calls.
        res_beta     <- update_beta_cohort_k(w, z_no_k, ZF[, k], ZF[, k]^2,
                                             beta_cohort_idx, C_beta,
                                             prior_family = prior_beta, alpha = alpha_iter,
                                             survival_divisor = beta_divisor)
        EBeta[k, ]   <- res_beta$mean
        EBeta2[k, ]  <- res_beta$second
      } else if (joint_beta) {
        # Joint q(beta) given the current scores of all programs (derivation, arm A2)
        VZ_all <- survival_Y^2 %*% pmax(EF2_aug[, 1:K, drop = FALSE] - EF_aug[, 1:K, drop = FALSE]^2, 0)
        jb <- update_beta_joint_ard(ZF, VZ_all, w, z, weight = beta_lik_weight,
                                    sigma2_init = beta_sigma2)
        beta_sigma2 <- jb$sigma2; EBeta_cov <- jb$Sigma; Mbeta <- jb$M
        EBeta <- jb$mean; EBeta2 <- diag(jb$M)
        res_beta <- NULL
      } else {
        res_beta    <- update_beta_k(w, z_no_k, ZF[, k],
                                     if (canonicalize) ZF[, k]^2 + VZ_k else ZF[, k]^2,
                                     prior_family = prior_beta,
                                     alpha = if (canonicalize) beta_lik_weight else alpha_iter,
                                     survival_divisor = if (canonicalize) 1 else beta_divisor,
                                     g_fixed = if (beta_prior_mode == "shared") g_beta else NULL)
        if (trace_supervision && beta_prior_mode == "shared") {
          trace_x[k] <- res_beta$x
          trace_s[k] <- res_beta$s
          trace_A[k] <- res_beta$A
        }
        EBeta[k]    <- res_beta$mean
        EBeta2[k]   <- res_beta$second
        if (beta_prior_mode == "shared") { beta_x[k] <- res_beta$x; beta_s[k] <- res_beta$s }
      }
      kl_beta[k]  <- if (joint_beta) jb$neg_kl / K   # joint -KL, spread so sum(kl_beta) is exact
                     else if (!is.null(res_beta$neg_kl)) res_beta$neg_kl
                     else compute_ebnm_kl(res_beta$ebnm_result$log_likelihood,
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
      if (!freeze_F && canonicalize) {
        # Exact sequential sweep over genes (multimodal_yfb_update_F_mk). The
        # survival terms are scaled by w_S through the working weights:
        # A_j = tau_j sum E[l^2] + w_S E[beta^2] sum_i w_i y_ij^2
        # B_j = tau_j (R_k' E[l])_j + w_S { E[beta] (Y'(w z_-k))_j
        #       - E[beta^2] sum_i w_i y_ij (Z_k - y_j E[f_j])_i }
        R_k    <- compute_R_k(Y, EL_aug, EF_aug, k)
        A_gen  <- Tau * sum(EL2_aug[, k])
        B_gen  <- Tau * as.vector(crossprod(R_k, EL_aug[, k]))
        if (is.null(groups) && is.null(prior_F_state[[k]])) {
          prior_F_state[[k]] <- if (prior_LF == "point_exponential")
            multimodal_yfb_fit_exponential_prior_ebnm(B_gen / A_gen, 1 / A_gen)
          else multimodal_yfb_fit_signed_prior(B_gen / A_gen, 1 / A_gen, prior_LF)
        }
        VZ_k <- as.vector(survival_Y^2 %*% pmax(EF2_aug[, k] - EF_aug[, k]^2, 0))
        A_surv_k <- EBeta2[k] * as.vector(crossprod(w, survival_Y^2))
        if (!is.null(groups)) {
          h_k <- if (joint_beta)
            w_S * w * (EBeta[k] * z - as.vector(ZF[, -k, drop = FALSE] %*% Mbeta[k, -k]))
          else w_S * w * z_no_k
          upd <- yf_update_F_groups(Y, Tau, EL_aug[, k], EL2_aug[, k], R_k,
            survival_Y, w_S * w, h_k, EF_aug[, k], EF2_aug[, k], EBeta[k], EBeta2[k],
            prior_F_state[[k]], groups, prior_LF, h_beta_weighted = joint_beta)
        } else if (max(A_gen + w_S * A_surv_k) < 1e-8) {
          # Dead program: neither likelihood carries information about F_k.
          # With the prior held fixed, its posterior would equal the prior
          # (zero mean, positive second moment), which revives L_k and beta_k
          # from nothing. Set it to exactly zero instead (equivalent to
          # pruning it while keeping K columns).
          upd <- list(EF = rep(0, p), EF2 = rep(0, p), kl = 0,
                      prior = list(family = prior_LF, pi = 0, rate = NA_real_,
                                   variance = 0, point_mass = TRUE),
                      details = list(A = A_gen + w_S * A_surv_k, B = B_gen))
        } else {
          # Joint q(beta): signal w (E[beta_k] z - sum_{l != k} E[beta_k beta_l] Z_l)
          h_k <- if (joint_beta)
            w_S * w * (EBeta[k] * z - as.vector(ZF[, -k, drop = FALSE] %*% Mbeta[k, -k]))
          else w_S * w * z_no_k
          upd <- multimodal_yfb_update_F_mk(
            survival_Y, Tau, EL_aug[, k], EL2_aug[, k], R_k, w_S * w, h_k,
            as.vector(survival_Y %*% EF_aug[, k]), VZ_k, EF_aug[, k], EF2_aug[, k],
            EBeta[k], EBeta2[k], prior_F_state[[k]], prior_update = "ebnm",
            loading_prior_update = loading_prior_update, h_beta_weighted = joint_beta)
        }
        prior_F_state[[k]] <- upd$prior
        EF[, k] <- upd$EF; EF2[, k] <- upd$EF2
        EF_aug[, k] <- upd$EF; EF2_aug[, k] <- upd$EF2
        kl_F[k] <- -upd$kl   # multimodal helper returns +KL; history stores -KL
        # Unweighted survival terms for the trace
        res_F <- list(A = upd$details$A, A_gen = A_gen,
                      A_surv = A_surv_k,
                      B_gen = B_gen,
                      B_surv = if (w_S > 0) (upd$details$B - B_gen) / w_S else 0 * B_gen)
      } else if (!freeze_F) {
        R_k        <- compute_R_k(Y, EL_aug, EF_aug, k)
        YtWz_no_k  <- as.vector(t(Y) %*% (w * z_no_k))  # p-vector
        # EBeta[k]/EBeta2[k] on a K x C_beta cohort-beta matrix would index the
        # wrong (flattened) element, not error -- harmless ONLY because
        # alpha_F > 0 with beta_cohort_id is rejected at the top of this
        # function (alpha=alpha_F=0 here zero-weights A_surv/B_surv below
        # regardless of what EBeta_k/EBeta2_k contain).
        res_F <- update_F_surv_YFB_k(Tau, EL_aug[, k], EL2_aug[, k], R_k,
                                      EBeta_k      = EBeta[k],
                                      EBeta2_k     = EBeta2[k],
                                      YtWY_diag    = YtWY_diag,
                                      YtWz_no_k    = YtWz_no_k,
                                      prior_family = prior_LF,
                                      alpha        = alpha_F,
                                      w_surv       = w_S,
                                      feature_group = feature_group)
        EF[, k]      <- res_F$mean
        EF2[, k]     <- res_F$second
        EF_aug[, k]  <- res_F$mean    # keep augmented matrix in sync
        EF2_aug[, k] <- res_F$second
        kl_F[k]  <- if (!is.null(groups)) sum(res_F$neg_kl) else
          compute_ebnm_kl(res_F$ebnm_result$log_likelihood,
                                     res_F$A, res_F$x, res_F$mean, res_F$second)
        if (!is.null(groups)) F_group_priors[[k]] <-
          lapply(res_F$ebnm_result, function(fit) fit$fitted_g)
      }
      if (!freeze_F) {
        if (trace_supervision) {
          sA  <- F_w_surv * res_F$A_surv / res_F$A
          bg  <- abs(F_w_gen * res_F$B_gen); bs <- abs(F_w_surv * res_F$B_surv)
          sB  <- ifelse(bg + bs > 0, bs / (bg + bs), 0)
          sA0 <- res_F$A_surv / (res_F$A_gen + res_F$A_surv)
          history$trace$shareA_median[iter, k] <- stats::median(sA, na.rm = TRUE)
          history$trace$shareA_p90[iter, k]    <- stats::quantile(sA, 0.9, names = FALSE, na.rm = TRUE)
          history$trace$shareB_median[iter, k] <- stats::median(sB, na.rm = TRUE)
          history$trace$shareB_p90[iter, k]    <- stats::quantile(sB, 0.9, names = FALSE, na.rm = TRUE)
          history$trace$shareA_unweighted_median[iter, k] <- stats::median(sA0, na.rm = TRUE)
        }
      }
      if (trace_supervision) {
        history$trace$beta_mean[iter, k]   <- EBeta[k]
        history$trace$beta_second[iter, k] <- EBeta2[k]
        history$trace$sd_YF[iter, k]       <- stats::sd(as.vector(projection_Y %*% EF[, k]))
        history$trace$norm_F[iter, k]      <- sqrt(sum(EF[, k]^2))
        history$trace$norm_L[iter, k]      <- sqrt(sum(EL[, k]^2))
      }

    }  # end k-loop

    # Shared beta prior: refit from the sweep's K pseudo-observations after
    # the warm-up (held fixed during the sweep, as for the loading priors)
    # Flat β likelihoods have s = Inf and carry no prior information.
    informative_beta <- is.finite(beta_s)
    if (beta_prior_mode == "shared" && !beta_prior_fix && iter > beta_warmup &&
        any(informative_beta)) {
      fit_g <- if (prior_beta == "point_normal" && !is.null(beta_slab_sd))
        ebnm::ebnm(beta_x[informative_beta], beta_s[informative_beta], prior_family = "point_normal", scale = beta_slab_sd)
      else ebnm::ebnm(beta_x[informative_beta], beta_s[informative_beta], prior_family = prior_beta)
      g_beta <- fit_g$fitted_g
    }

    if (trace_supervision && !use_beta_cohort && beta_prior_mode == "shared") {
      # Raw pseudo-observations precede the end-of-sweep canonical rescaling.
      # Record the exact prior used, and the refit passed to the next sweep.
      history$prior_trace[[iter]] <- list(prior_before = trace_prior_before,
        prior_after = g_beta, x = trace_x, s = trace_s, A = trace_A,
        beta_mean = EBeta, beta_sd = sqrt(pmax(EBeta2 - EBeta^2, 0)))
    }

    # Fail loud on divergence: a non-finite loading or coefficient otherwise
    # propagates silently (ebnm returns zeros and the fit looks converged).
    if (!all(is.finite(EF)) || !all(is.finite(EL)) || !all(is.finite(EBeta)))
      stop(sprintf("fit_ssbmf_yfb() diverged at iteration %d: non-finite loadings or ",
                   iter), "survival coefficients (max |E[beta]| before failure: ",
           signif(max(abs(EBeta_old)), 3), ").")

    # Cohort F update — identical block to fit_ssbmf_lb.R.
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

    # Development guard: ZF must be n x K; EBeta must be K-vector (or K x C_beta).
    if (use_beta_cohort) {
      stopifnot(ncol(ZF) == K, all(dim(EBeta) == c(K, C_beta)))
    } else {
      stopifnot(ncol(ZF) == K, length(EBeta) == K)
    }

    # Full ELBO under Cluster B: ZF = Y·EF is the survival predictor.
    # Pass ZF^2 as EL2 (element-wise squared): ZF is observed, so its
    # posterior second moment equals its squared value (no variance term).
    surv_elbo <- if (canonicalize) {
      # End-of-iteration objective: Cox partial log-likelihood re-expanded at
      # the current predictor, with Var(eta) from q(beta) and q(F)
      ZF_end <- projection_Y %*% EF_aug[, 1:K, drop = FALSE]
      VZ_end <- (projection_Y^2) %*% pmax(EF2_aug[, 1:K, drop = FALSE] - EF_aug[, 1:K, drop = FALSE]^2, 0)
      tay_end <- calc_cox_taylor_yf(as.vector(ZF_end %*% EBeta), time, status,
                                    strata = strata_id)
      if (joint_beta) {
        # Var(eta_i) = Zbar_i' Sigma Zbar_i + sum_k VZ_ik E[beta_k^2]
        var_eta <- rowSums((ZF_end %*% EBeta_cov) * ZF_end) + as.vector(VZ_end %*% EBeta2)
        tay_end$logPL - 0.5 * sum(tay_end$w * var_eta)
      } else compute_survival_elbo(tay_end$logPL, tay_end$w, ZF_end, ZF_end^2 + VZ_end,
                                   EBeta, EBeta2)
    } else compute_survival_elbo(
      taylor$logPL, w, ZF, ZF^2, EBeta, EBeta2,
      cohort_idx = if (use_beta_cohort) beta_cohort_idx else NULL
    )
    history$elbo_full[iter] <- if (canonicalize) {
      # Same weighting as the updates: genomics + w_S * survival + KL terms; under the
      # survival-block weight the beta KL is also multiplied by w_S
      res_tau$elbo_proxy + w_S * surv_elbo + sum(kl_L) + sum(kl_F) +
        (if (survival_block_weight) w_S else 1) * sum(kl_beta)
    } else {
      (1 - alpha) * (res_tau$elbo_proxy / genomics_divisor) +
        alpha * (surv_elbo / survival_divisor) +
        sum(kl_L) + sum(kl_F) + sum(kl_beta)
    }
    if (!is.null(cohort_id) && C_cols > 0) {
      history$elbo_full[iter] <- history$elbo_full[iter] +
        compute_normal_kl(
          EF_aug[,  K + seq_len(C_cols), drop = FALSE],
          EF2_aug[, K + seq_len(C_cols), drop = FALSE],
          sigma_F_cohort
        )
    }

    # Scale canonicalization (after the ELBO, which used this iteration's ZF).
    # Leaves L F', eta and the ELBO unchanged; fixes sd(Y EF_k) = 1 for the
    # next iteration's ZF.
    if (canonicalize) {
      cz <- canonicalize_yf_factors(projection_Y, EL, EL2, EF, EF2, EBeta, EBeta2, K,
                                    prior_F = prior_F_state, feature_group = groups)
      EL <- cz$EL; EL2 <- cz$EL2; EF <- cz$EF; EF2 <- cz$EF2
      EBeta <- cz$EBeta; EBeta2 <- cz$EBeta2; prior_F_state <- cz$prior_F
      if (joint_beta) {
        # beta_k -> c_k beta_k: Sigma_kl -> c_k c_l Sigma_kl, sigma2_k -> c_k^2 sigma2_k
        EBeta_cov <- EBeta_cov * tcrossprod(cz$scale)
        Mbeta <- tcrossprod(EBeta) + EBeta_cov
        if (!is.null(beta_sigma2)) beta_sigma2 <- beta_sigma2 * cz$scale^2
      }
      EL_aug[, 1:K] <- EL; EL2_aug[, 1:K] <- EL2
      EF_aug[, 1:K] <- EF; EF2_aug[, 1:K] <- EF2
    }

    # ========================================================================
    # STEP 4: Convergence Check
    # ========================================================================
    if (record_output_convergence) {
      reconstruction_change <- multimodal_yfb_relative_change(
        old_reconstruction, EL %*% t(EF))
      eta_change <- multimodal_yfb_relative_change(old_eta_outputs,
        as.vector(projection_Y %*% EF %*% EBeta))
      if (iter == 1L) history$output_change <- numeric(max_iter)
      history$output_change[iter] <- max(reconstruction_change, eta_change)
      history$outputs_converged <- history$output_change[iter] < 1e-4
    }
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

    # iter > N_frozen: convergence cannot fire while F is still held fixed --
    # beta/L can plateau during the frozen phase on their own, and without this
    # guard the loop would exit before F ever gets a chance to unfreeze,
    # silently defeating the N_frozen mechanism (DECISIONS.md 2026-08-20).
    elbo_converged <- is.finite(history$delta_elbo_rel[iter]) &&
      history$delta_elbo_rel[iter] < tol
    if (record_output_convergence) history$elbo_converged <- elbo_converged
    stop_now <- if (converge_on == "elbo") elbo_converged else history$outputs_converged
    if (iter > 5 && iter > N_frozen && stop_now) {
      if (verbose) {
        cat(sprintf("\n  Converged at iteration %d  (relative ELBO change = %.2e)\n",
                    iter, history$delta_elbo_rel[iter]))
      }
      history$converged  <- TRUE
      if (record_output_convergence)
        history$output_change <- history$output_change[seq_len(iter)]
      history$n_iter     <- iter
      history$rmse       <- history$rmse[1:iter]
      history$elbo_proxy <- history$elbo_proxy[1:iter]
      history$elbo_full  <- history$elbo_full[1:iter]
      history$delta_L    <- history$delta_L[1:iter]
      history$delta_Beta <- history$delta_Beta[1:iter]
      history$delta_elbo_rel <- history$delta_elbo_rel[1:iter]
      history$factor_pve <- history$factor_pve[1:iter, , drop = FALSE]
      if (trace_supervision) {
        history$trace <- lapply(history$trace, function(m) m[1:iter, , drop = FALSE])
        if (!use_beta_cohort && beta_prior_mode == "shared") history$prior_trace <- history$prior_trace[seq_len(iter)]
      }
      break
    }

  }  # end CAVI loop

  # No post-fit sign flip (sign_correction removed 2026-11-06, DECISIONS.md):
  # CAVI fits beta to the Cox likelihood given ZF, so a converged fit is
  # already oriented, and a flip would hide a fit that is not.

  stopifnot(length(EF_norms) == K)  # guard: must be K-vector for prediction

  if (use_beta_cohort) {
    # Column names are how predict_ssbmf_yfb() matches a test cohort's label
    # against the right column of EBeta.
    colnames(EBeta)  <- beta_cohort_levels
    colnames(EBeta2) <- beta_cohort_levels
  }

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
    beta_divisor     = beta_divisor,
    canonicalize     = canonicalize,
    w_S              = w_S,
    g_beta           = if (beta_prior_mode == "shared") g_beta else NULL,
    neg_kl           = kl_L + kl_F +           # per-program -KL terms (last iteration)
                       (if (survival_block_weight) w_S else 1) * kl_beta
  )
  # S-13 fields only when an S-13 option is on, so default fit objects are unchanged
  if (survival_block_weight || joint_beta) {
    result$survival_block_weight <- survival_block_weight
    result$beta_posterior <- beta_posterior
    result$EBeta_cov <- if (joint_beta) EBeta_cov else diag(pmax(EBeta2 - EBeta^2, 0), length(EBeta))
    result$beta_prior_var <- if (joint_beta) beta_sigma2 else NULL
  }
  if (!is.null(cohort_id) && C_cols > 0) {
    result$L_cohort   <- L_cohort
    result$EF_cohort  <- EF_aug[, K + seq_len(C_cols), drop = FALSE]
    result$EF2_cohort <- EF2_aug[, K + seq_len(C_cols), drop = FALSE]
    result$cohort_id  <- cohort_id
  }
  if (use_beta_cohort) {
    # EBeta_pooled: for external prediction, where no beta_cohort_id is
    # available (code/update_beta_cohort.R's compute_pooled_beta()).
    # Computed from POOLED PATIENT-SUMS at the final converged w/z/ZF (the
    # last main-loop iteration's values, still in scope here). The result is
    # NOT simply rowMeans(EBeta), which is wrong under unequal cohort sizes;
    # rowMeans(EBeta) is only the starting value of the single Gauss-Seidel
    # sweep in compute_pooled_beta(), and it is the state that produced
    # w/z/ZF, so the sweep's inputs are consistent.
    EBeta_pooled <- compute_pooled_beta(
      w, z, ZF, rowMeans(EBeta),
      prior_family = prior_beta, alpha = alpha, survival_divisor = beta_divisor
    )
    result$EBeta_pooled <- EBeta_pooled
    result$beta_cohort_id    <- beta_cohort_id
    result$beta_cohort_levels <- beta_cohort_levels
  }

  if(prune_preserve_posterior || !is.null(posterior_init)) {
    if (!is.null(posterior_init)) result$initial_state <- initial_state
    result$posterior_state <- list(EL=result$EL,EL2=result$EL2,EF=result$EF,EF2=result$EF2,
      EBeta=result$EBeta,EBeta2=result$EBeta2,Tau=result$Tau,
      prior_F=prior_F_state,g_beta=result$g_beta)
  }

  if(profile_cox_intercept)result$profile_cox_intercept <- TRUE
  if (center_projection) result$projection_center <- projection_center
  if (!is.null(groups)) {
    result$feature_group <- as.character(feature_group)
    result$prior_F <- stats::setNames(lapply(names(groups), function(g) {
      if (canonicalize) lapply(prior_F_state, function(pk) pk[[g]]) else
        lapply(seq_len(K), function(k) F_group_priors[[k]][[g]])
    }), names(groups))
  }
  if(identical(loading_prior_update,"fixed_q_em"))
    result$loading_prior_update <- loading_prior_update

  # Pruning (canonical only): drop the program whose removal raises the ELBO
  # most, if any removal does not lower it, and refit from a warm start.
  if (prune && K > 1) {
    delta <- prune_yf_nullcheck(result, Y, time, status, strata = strata_id,
                                protect_survival = prune_protect_survival,
                                survival_Y = projection_Y)
    if (max(delta) >= 0) {
      drop <- which.max(delta); keep <- setdiff(seq_len(K), drop)
      refit_args <- modifyList(call_args, list(
        K = K - 1L, EL_init = EL[, keep, drop = FALSE], EF_init = EF[, keep, drop = FALSE],
        init_method = "custom", init_supervised = FALSE, beta_warmup = 0L,
        beta_prior_init = if (beta_prior_mode == "shared") g_beta else NULL))
      refit_args["posterior_init"] <- list(NULL)
      if(prune_preserve_posterior)
        refit_args$posterior_init <- yf_subset_posterior_state(result$posterior_state,keep)
      if (beta_prior_mode == "shared") refit_args["beta_prior_init"] <- list(g_beta)
      if (!record_output_convergence) refit_args$converge_on <- NULL
      refit <- do.call(fit_ssbmf_yfb, refit_args)
      refit$pruning <- rbind(
        data.frame(K_before = K, removed_program = drop, elbo_gain = delta[drop],
                   converged_before = history$converged),
        refit$pruning)
      return(refit)
    }
  }
  result
}

# ==============================================================================
# PART 3 — DATA_MODE Runner Block (for standalone Rscript use)
# Skipped when the file is sourced (sys.nframe() > 0), so source() only
# defines the functions.
# ==============================================================================

if (sys.nframe() == 0L && DATA_MODE == "real") {
  stopifnot(!is.null(real_Y), !is.null(real_time), !is.null(real_status))
  res <- fit_ssbmf_yfb(real_Y, real_time, real_status,
                       K = 5, max_iter = 100, tol = 1e-3, verbose = TRUE)
  cat(sprintf("  Converged: %s | Iterations: %d\n",
              res$history$converged, res$history$n_iter))
}
