# ============================================================
# Script: feature_groups_yfb.R
# Purpose: Fit separate loading priors while retaining the joint YFB projection.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: ebnm and existing multimodal posterior/KL helpers
# ============================================================

#' Validate feature labels and return blocks in first-appearance order
#' @param feature_group NULL or one nonmissing, nonempty label per feature.
#' @param p Number of matrix columns.
#' @return NULL for zero/one group, otherwise a named list of column indices.
#' @examples
#' yf_feature_groups(c("rna", "meth", "rna"), 3)
#' @family F_surv_YFB_update
yf_feature_groups <- function(feature_group, p) {
  if (is.null(feature_group)) return(NULL)
  if (!is.atomic(feature_group) || !is.null(dim(feature_group)) ||
      length(feature_group) != p || anyNA(feature_group) ||
      any(!nzchar(as.character(feature_group))))
    stop("feature_group must have one nonmissing, nonempty label per feature.")
  labels <- as.character(feature_group)
  levels <- unique(labels)
  if (length(levels) < 2L) return(NULL)
  stats::setNames(lapply(levels, function(g) which(labels == g)), levels)
}

#' Solve independent Gaussian empirical-Bayes problems within feature groups
#' @param x,s Pseudo-observations and their standard deviations.
#' @param groups Named list of feature indices.
#' @param family ebnm prior family.
#' @return Combined posterior mean/SD/second moment, group fits and negative KLs.
#' @examples
#' yf_group_ebnm(c(-2, 0, 3, 0), rep(1, 4), list(a = 1:2, b = 3:4), "point_laplace")
#' @family F_surv_YFB_update
yf_group_ebnm <- function(x, s, groups, family) {
  mean <- sd <- numeric(length(x))
  fits <- vector("list", length(groups)); names(fits) <- names(groups)
  neg_kl <- stats::setNames(numeric(length(groups)), names(groups))
  for (g in names(groups)) {
    j <- groups[[g]]
    fit <- ebnm::ebnm(x[j], s[j], prior_family = family)
    mean[j] <- fit$posterior$mean; sd[j] <- fit$posterior$sd
    second <- mean[j]^2 + sd[j]^2
    expected <- sum(dnorm(x[j], 0, s[j], log = TRUE) +
                      (x[j] * mean[j] - second / 2) / s[j]^2)
    neg_kl[g] <- fit$log_likelihood - expected
    fits[[g]] <- fit
  }
  list(mean = mean, sd = sd, second = mean^2 + sd^2,
       ebnm_result = fits, neg_kl = neg_kl)
}

#' Sum the loading KL under the refitted prior without changing q
#'
#' KL(q||g_new) = KL(q||g_old) + E_q[log(g_old/g_new)]. The sequential
#' sweep holds g_old fixed, so refitting g requires this density correction.
#' @param update Real sequential block update with posterior diagnostics.
#' @param old Prior used during that sweep.
#' @return Positive KL under update$prior. Incompatible support raises an error.
#' @examples
#' # yf_group_loading_kl(block_update, incoming_prior)
#' @family F_surv_YFB_update
#' @seealso [loading_prior_laplace_stats()]
yf_group_loading_kl <- function(update, old) {
  new <- update$prior
  v <- update$details
  if (identical(old$family, "normal")) {
    old_log <- if (isTRUE(old$point_mass)) 0 else
      sum(-0.5 * log(2 * pi * old$variance) - v$second / (2 * old$variance))
    new_log <- if (isTRUE(new$point_mass)) {
      if (any(v$second > 0)) -Inf else 0
    } else sum(-0.5 * log(2 * pi * new$variance) - v$second / (2 * new$variance))
    # A continuous normal and a delta-zero distribution have singular support.
    if (isTRUE(old$point_mass) && !isTRUE(new$point_mass)) new_log <- -Inf
  } else {
    q <- if (identical(old$family, "point_laplace"))
      loading_prior_laplace_stats(v$A, v$B, old) else
        list(prob = v$slab_prob, mean_abs = v$mean)
    # The exponential and Laplace differ by a constant log(2) in the slab
    # density. It cancels between priors of the same family.
    old_laplace <- old; old_laplace$family <- "point_laplace"
    new_laplace <- new; new_laplace$family <- "point_laplace"
    old_log <- loading_prior_laplace_log(q$prob, q$mean_abs, old_laplace)
    new_log <- loading_prior_laplace_log(q$prob, q$mean_abs, new_laplace)
  }
  kl <- update$kl + old_log - new_log
  if (!is.finite(kl))
    stop("Nonfinite feature_group F KL: the refitted prior is incompatible with its sweep posterior.")
  kl
}

#' Update one program across feature blocks using the joint survival projection
#'
#' Every block receives the same current Z_k = Σ_m Y_m F_mk and V^Z_k.
#' Its sequential sweep refreshes these before the next block is visited.
#' @param Y,Tau,EL_k,EL2_k,R_k Reconstruction inputs for one program.
#' @param survival_Y Projection design, possibly centered within strata.
#' @param w,h_minus_k Weighted Cox curvature and working signal.
#' @param EF_k,EF2_k Loading moments across all features.
#' @param EBeta_k,EBeta2_k Coefficient moments.
#' @param priors Named list of persistent priors, one per group.
#' @param groups Named index list from yf_feature_groups().
#' @param family Loading prior family.
#' @param h_beta_weighted Whether h_minus_k already includes beta moments.
#' @return Combined moments, joint EZ/VZ, group priors, summed positive KL,
#'   A/B diagnostics and named counts of refitted point-mass block resets.
#' @examples
#' # Called by fit_ssbmf_yfb() for explicit multiple feature groups.
#' @family F_surv_YFB_update
#' @seealso [multimodal_yfb_update_F_mk()]
yf_update_F_groups <- function(Y, Tau, EL_k, EL2_k, R_k, survival_Y, w,
                                h_minus_k, EF_k, EF2_k, EBeta_k, EBeta2_k,
                                priors, groups, family, h_beta_weighted = FALSE) {
  EZ <- as.vector(survival_Y %*% EF_k)
  VZ <- as.vector(survival_Y^2 %*% pmax(EF2_k - EF_k^2, 0))
  A_gen <- Tau * sum(EL2_k)
  B_gen <- Tau * as.vector(crossprod(R_k, EL_k))
  A <- B <- numeric(ncol(Y)); kl <- 0
  empty_blocks <- stats::setNames(integer(length(groups)), names(groups))
  if (is.null(priors)) priors <- stats::setNames(vector("list", length(groups)), names(groups))
  for (g in names(groups)) {
    j <- groups[[g]]
    if (is.null(priors[[g]])) {
      priors[[g]] <- if (family == "point_exponential")
        multimodal_yfb_fit_exponential_prior_ebnm(B_gen[j] / A_gen[j], 1 / A_gen[j]) else
          multimodal_yfb_fit_signed_prior(B_gen[j] / A_gen[j], 1 / A_gen[j], family)
    }
    total_A <- A_gen[j] + EBeta2_k * as.vector(crossprod(w, survival_Y[, j, drop = FALSE]^2))
    if (max(total_A) < 1e-8) {
      # Same dead-program threshold as the single-group path. Remove this
      # block's contribution from the joint projection before visiting another.
      EZ <- EZ - as.vector(survival_Y[, j, drop = FALSE] %*% EF_k[j])
      VZ <- pmax(VZ - as.vector(survival_Y[, j, drop = FALSE]^2 %*%
                                 pmax(EF2_k[j] - EF_k[j]^2, 0)), 0)
      EF_k[j] <- EF2_k[j] <- 0
      priors[[g]] <- list(family = family, pi = 0, rate = NA_real_,
                          variance = 0, point_mass = TRUE)
      A[j] <- total_A; B[j] <- B_gen[j]
    } else {
      update <- multimodal_yfb_update_F_mk(
        Y_m = survival_Y[, j, drop = FALSE], Tau_m = Tau[j], EL_k = EL_k,
        EL2_k = EL2_k, R_mk = R_k[, j, drop = FALSE], w = w,
        h_minus_k = h_minus_k, EZ_k = EZ, VZ_k = VZ,
        EF_k = EF_k[j], EF2_k = EF2_k[j], EBeta_k = EBeta_k,
        EBeta2_k = EBeta2_k, prior = priors[[g]], prior_update = "ebnm",
        h_beta_weighted = h_beta_weighted)
      if (isTRUE(update$prior$point_mass)) {
        # A refitted delta_0 prior makes this modality-program block empty.
        # Remove the moments JUST returned by the sweep, whose EZ/VZ already
        # include them. Removing the incoming moments would leave stale joint
        # projections for the next modality. q = g = delta_0 gives KL = 0.
        update$EZ <- update$EZ - as.vector(survival_Y[, j, drop = FALSE] %*% update$EF)
        update$VZ <- pmax(update$VZ - as.vector(survival_Y[, j, drop = FALSE]^2 %*%
          pmax(update$EF2 - update$EF^2, 0)), 0)
        update$EF[] <- update$EF2[] <- 0
        empty_blocks[g] <- empty_blocks[g] + 1L
      } else {
        # All other support/numerical failures remain explicit errors.
        kl <- kl + yf_group_loading_kl(update, priors[[g]])
      }
      EF_k[j] <- update$EF; EF2_k[j] <- update$EF2
      EZ <- update$EZ; VZ <- update$VZ
      priors[[g]] <- update$prior
      A[j] <- update$details$A; B[j] <- update$details$B
    }
  }
  list(EF = EF_k, EF2 = EF2_k, prior = priors, kl = kl, EZ = EZ, VZ = VZ,
       empty_blocks = empty_blocks, details = list(A = A, B = B))
}
