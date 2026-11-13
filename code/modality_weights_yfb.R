# ============================================================
# Script: modality_weights_yfb.R
# Purpose: Validate reconstruction weights and compute training-only block rules.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R
# ============================================================

#' Expand named positive modality weights to feature order
#'
#' NULL and all-one vectors return NULL so the historical arithmetic is retained.
#' Labels describe data blocks independently of whether F priors are pooled.
#' @param modality_weight NULL or a named positive finite numeric vector.
#' @param feature_group One modality label per feature, required for weights.
#' @param p Number of features.
#' @return NULL or an unnamed numeric p-vector of reconstruction powers.
#' @examples
#' yf_modality_weights(c(rna = 1, meth = 0.2), c("rna", "meth"), 2)
#' @family F_surv_YFB_update
yf_modality_weights <- function(modality_weight, feature_group, p) {
  if (is.null(modality_weight)) return(NULL)
  if (is.null(feature_group)) stop("modality_weight requires feature_group.")
  yf_feature_groups(feature_group, p)
  labels <- unique(as.character(feature_group))
  if (!is.numeric(modality_weight) || !is.null(dim(modality_weight)) ||
      any(!is.finite(modality_weight)) || any(modality_weight <= 0) ||
      is.null(names(modality_weight)) || anyDuplicated(names(modality_weight)) ||
      length(modality_weight) != length(labels) || !setequal(names(modality_weight), labels))
    stop("modality_weight must be a named positive finite numeric vector keyed by feature_group.")
  if (all(modality_weight == 1)) return(NULL)
  unname(modality_weight[as.character(feature_group)])
}

#' Compute both registered weight rules from training rows only
#'
#' Each block is z-scored using its training means and sample SDs. The leading
#' eigenvalue of X'X/n is computed from XX'/n to avoid a p-by-p CpG matrix.
#' @param Y Numeric training matrix, patients by features. Do not include held-out rows.
#' @param feature_group One modality label per column.
#' @param reference Reference label whose reconstruction power is one.
#' @return Named list of named numeric vectors: feature_count and eigenvalue.
#' @examples
#' # yf_modality_weight_rules(Y_train, feature_group)
#' @family F_surv_YFB_update
#' @seealso [yf_modality_weights()]
yf_modality_weight_rules <- function(Y, feature_group, reference = "expression") {
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) || nrow(Y) < 2L)
    stop("Weight rules require a finite numeric training matrix with at least two rows.")
  yf_feature_groups(feature_group, ncol(Y))
  labels <- unique(as.character(feature_group))
  if (!reference %in% labels) stop("Weight reference is absent from feature_group.")
  count <- leading <- stats::setNames(numeric(length(labels)), labels)
  for (m in labels) {
    X <- Y[, feature_group == m, drop = FALSE]
    count[m] <- ncol(X)
    X <- scale(X)
    if (any(!is.finite(X))) stop("Weight rules cannot z-score a constant training feature.")
    leading[m] <- eigen(tcrossprod(X) / nrow(X), symmetric = TRUE,
                        only.values = TRUE)$values[1L]
  }
  if (any(!is.finite(leading)) || any(leading <= 0)) stop("Invalid leading block eigenvalue.")
  list(feature_count = count[[reference]] / count,
       eigenvalue = leading[[reference]] / leading)
}

#' Evaluate the weighted Gaussian reconstruction contribution at supplied tau
#' @param res_tau Result of update_tau(), including R2_bar and Tau.
#' @param feature_weight NULL or positive weight per feature.
#' @return Scalar log-likelihood without the parameter-independent normal constant.
#' @examples
#' # yf_weighted_genomics(res_tau, feature_weight)
#' @family F_surv_YFB_update
yf_weighted_genomics <- function(res_tau, feature_weight = NULL) {
  if (is.null(feature_weight)) return(res_tau$elbo_proxy)
  n <- nrow(res_tau$R2_bar)
  sum(feature_weight * (n / 2 * log(res_tau$Tau) -
                          res_tau$Tau / 2 * colSums(res_tau$R2_bar)))
}
