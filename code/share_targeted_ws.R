# ============================================================
# Script: share_targeted_ws.R
# Purpose: Calculate directional loading precision shares and target weights.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: code/fit_ssbmf_yfb.R (calc_cox_taylor_yf)
# ============================================================

#' Calculate a survival weight from a directional precision share
#'
#' At the supplied canonical fit state, r_k = w_S E[β_k²]
#' f̂_k'Y'WY f̂_k / (sum_i E[l_ik²] · sum_j τ_j f̂_jk²). The target weight
#' inverts s = r/(1+r), holding all posterior moments and W fixed.
#' The closed forms are approximations, not inputs to the weight. Mean τ appears only in the closed-form approximations. This helper supports pooled, unweighted fits only. It rejects
#' stored projection centering and other identifiable unsupported controls.
#' The caller must supply the same training matrix used in the fit. Some fitter
#' options are not stored in the fit, so the caller must verify this scope.
#' @param fit Converged canonical YFB fit with EL, EF, EBeta, EBeta2, Tau and history.
#' @param Y Numeric training matrix, patients by features, used for fitting.
#' @param time,status Training event times and 0/1 event indicators.
#' @param target_share Desired directional share strictly between zero and one.
#' @param pilot_loading Optional fixed pilot loading vector for tracking across passes.
#'   Omit for the initial w_S = 1 pilot.
#' @param w_current Positive weight at the fit state, defaulting to fit$w_S or 1.
#' @return List containing w_S, selected_program, per-program diagnostics, ratio,
#'   no_signal_guard, capped, pilot_loading, lost_activity, matching_correlation,
#'   and loading_second_moment_bound. Missing EL2 is
#'   replaced by EL², giving an upper bound on the precision ratios and shares.
#' @examples
#' # Source fit_ssbmf_yfb.R first, then call with a converged canonical fit.
#' @seealso calc_cox_taylor_yf
share_targeted_ws <- function(fit, Y, time, status, target_share,
                              w_current = if (is.null(fit$w_S)) 1 else fit$w_S,
                              pilot_loading = NULL) {
  #' Check whether an argument is a finite numeric scalar.
  #' @param x Candidate value.
  #' @return Logical scalar.
  #' @examples
  #' # scalar(1)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is.list(fit) || !isTRUE(fit$history$converged)) stop("fit must be converged")
  if (!isTRUE(fit$canonicalize)) stop("fit must use canonical scaling")
  unsupported <- c("projection_center", "cohort_id", "beta_cohort_id",
                   "modality_weight", "feature_group", "strata_id")
  if (any(!vapply(fit[unsupported], is.null, logical(1))) ||
      isTRUE(fit$profile_cox_intercept)) stop("unsupported fit controls: pooled unweighted input required")
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) || nrow(Y) < 2L || ncol(Y) < 1L)
    stop("Y must be a finite numeric matrix with at least two rows")
  n <- nrow(Y); p <- ncol(Y)
  if (!is.numeric(time) || length(time) != n || any(!is.finite(time)) ||
      any(time < 0) || !is.numeric(status) || length(status) != n ||
      any(!is.finite(status)) || any(!status %in% c(0, 1))) stop("invalid time or status")
  if (!scalar(target_share) || target_share <= 0 || target_share >= 1 ||
      !scalar(w_current) || w_current <= 0) stop("invalid target_share or w_current")
  if (!is.null(fit$w_S) && (!scalar(fit$w_S) || fit$w_S != w_current))
    stop("w_current must match fit$w_S")
  if (is.null(pilot_loading) && w_current != 1) stop("initial pilot must have w_S = 1")
  if (!is.null(pilot_loading) && (!is.numeric(pilot_loading) || !is.null(dim(pilot_loading)) ||
      length(pilot_loading) != p || any(!is.finite(pilot_loading)) ||
      stats::sd(pilot_loading) == 0)) stop("invalid pilot_loading")
  EL <- fit$EL; EF <- fit$EF; beta <- fit$EBeta; beta2 <- fit$EBeta2; tau <- fit$Tau
  if (!is.matrix(EL) || !is.numeric(EL) || !is.matrix(EF) || !is.numeric(EF) ||
      nrow(EL) != n || nrow(EF) != p || ncol(EL) < 1L || ncol(EL) != ncol(EF) ||
      any(!is.finite(EL)) || any(!is.finite(EF))) stop("invalid EL or EF")
  K <- ncol(EF)
  if (!is.numeric(beta) || !is.null(dim(beta)) || length(beta) != K || any(!is.finite(beta)) ||
      !is.numeric(beta2) || length(beta2) != K || any(!is.finite(beta2)) ||
      any(beta2 < beta^2) || !is.numeric(tau) || !length(tau) %in% c(1L, p) ||
      any(!is.finite(tau)) || any(tau <= 0)) stop("invalid beta moments or Tau")
  bound <- is.null(fit$EL2)
  EL2 <- if (bound) EL^2 else fit$EL2
  if (!is.matrix(EL2) || !is.numeric(EL2) || !identical(dim(EL2), dim(EL)) ||
      any(!is.finite(EL2)) || any(EL2 < EL^2)) stop("invalid EL2")
  if (!exists("calc_cox_taylor_yf", mode = "function")) stop("source code/fit_ssbmf_yfb.R first")
  scores <- Y %*% EF
  weights <- calc_cox_taylor_yf(as.vector(scores %*% beta), time, status)$w
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0)) stop("invalid Cox working weights")
  tau <- rep(tau, length.out = p)
  norm_f <- sqrt(colSums(EF^2))
  sum_l2 <- colSums(EL2)
  sd_score <- apply(scores, 2, stats::sd)
  b <- beta * sd_score
  beta_sd <- sqrt(pmax(beta2 - beta^2, 0))
  active <- beta_sd > 0 & abs(beta) >= 1.96 * beta_sd &
    colSums(EF^2) > 1e-8 & sd_score > 0
  pve <- colSums(EL^2) * colSums(EF^2) / sum(Y^2)
  directional <- per_gene <- approx_directional <- approx_gene <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    if (norm_f[k] > 0 && sum_l2[k] > 0) {
      projected <- as.vector(Y %*% (EF[, k] / norm_f[k]))
      directional[k] <- w_current * beta2[k] * sum(weights * projected^2) / (sum_l2[k] * sum(tau * (EF[, k] / norm_f[k])^2))
      gene_ratio <- w_current * beta2[k] * colSums(Y^2 * weights) / (tau * sum_l2[k])
      per_gene[k] <- stats::median(gene_ratio / (1 + gene_ratio))
    }
    if (is.finite(pve[k]) && pve[k] > 0) {
      approx_directional[k] <- w_current * mean(weights) * b[k]^2 / (mean(tau) * p * pve[k])
      approx_gene[k] <- approx_directional[k] / (p * pve[k])
    }
  }
  diagnostics <- data.frame(program = seq_len(K), active = active, standardized_beta = b,
    pve = pve, directional_ratio = directional, directional_share = directional / (1 + directional),
    per_gene_median_share = per_gene, closed_directional_ratio = approx_directional,
    closed_directional_share = approx_directional / (1 + approx_directional),
    closed_gene_ratio = approx_gene, closed_gene_share = approx_gene / (1 + approx_gene))
  if (is.null(pilot_loading)) {
    selected <- if (any(active)) which.max(replace(abs(b), !active, -Inf)) else NA_integer_
    if (is.na(selected)) return(list(w_S = 1, selected_program = selected, diagnostics = diagnostics,
      ratio = NA_real_, no_signal_guard = TRUE, capped = FALSE, pilot_loading = NULL,
      lost_activity = FALSE, matching_correlation = NA_real_, loading_second_moment_bound = bound))
    pilot_loading <- EF[, selected]
    if (stats::sd(pilot_loading) == 0) stop("selected pilot loading has undefined correlation")
    matching_correlation <- 1
  } else {
    correlations <- vapply(seq_len(K), function(k) {
      if (stats::sd(EF[, k]) == 0) return(NA_real_)
      abs(stats::cor(pilot_loading, EF[, k]))
    }, numeric(1))
    if (all(is.na(correlations))) stop("no loading has a defined pilot correlation")
    selected <- which.max(replace(correlations, is.na(correlations), -Inf))
    matching_correlation <- correlations[selected]
  }
  ratio <- directional[selected]
  if (!is.finite(ratio) || ratio <= 0) stop("selected program has undefined or zero directional ratio")
  proposed <- w_current * (target_share / (1 - target_share)) / ratio
  list(w_S = min(proposed, 1e5), selected_program = selected, diagnostics = diagnostics,
       ratio = ratio, no_signal_guard = FALSE, capped = proposed > 1e5,
       pilot_loading = pilot_loading, lost_activity = !active[selected],
       matching_correlation = matching_correlation, loading_second_moment_bound = bound)
}
