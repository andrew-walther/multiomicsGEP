# ============================================================
# Script: predict_cox_on_yf.R
# Purpose: Cluster B prediction — direct Y·EF projection (same formula as training)
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-04
# Dependencies: none (base R only)
# ============================================================

#' Cluster B prediction: project test patients via Y_test · EF_norm · beta_tilde
#'
#' Under the Cox-on-YF reformulation (eta = ZF · beta_tilde where ZF = Y·EF_norm),
#' prediction on new data uses the SAME formula as training. EF columns are
#' normalized to unit L2 norm before projection, matching the normalization
#' applied during training in fit_cox_on_yf(). EF_norms is stored in the
#' fitted model object ($EF_norms) and must be passed here to ensure consistency.
#'
#' @param Y_test   n_test × p numeric matrix: test patients' genomics data.
#'                 Must have the same column ordering as training data (or a
#'                 subset matched by gene index before calling).
#' @param EF       p × K numeric matrix: posterior mean factor weights from
#'                 fit_cox_on_yf()$EF. May be a row-subset for external validation.
#' @param EBeta    K-vector: posterior mean survival coefficients (fit_cox_on_yf()$EBeta).
#' @param EF_norms K-vector of column norms from training (fit_cox_on_yf()$EF_norms).
#'                 Required for correct normalization; if NULL, EF is used as-is
#'                 (not recommended — predictions will be on a different scale).
#'
#' @return Named list:
#'   $ZF_test      n_test × K matrix of normalized projection scores
#'   $risk_scores  n_test-vector of Cox linear predictor values
#'
#' @examples
#' # res <- fit_cox_on_yf(Y_train, time_train, status_train, K=5)
#' # pred <- predict_cox_on_yf(Y_test, res$EF, res$EBeta, res$EF_norms)
#' # concordance(Surv(time_test, status_test) ~ pred$risk_scores)
predict_cox_on_yf <- function(Y_test, EF, EBeta, EF_norms = NULL) {

  if (!is.matrix(Y_test) || !is.numeric(Y_test))
    stop("Y_test must be a numeric matrix.")
  if (!is.matrix(EF) || !is.numeric(EF))
    stop("EF must be a numeric matrix.")
  if (!is.numeric(EBeta))
    stop("EBeta must be a numeric vector.")
  if (ncol(Y_test) != nrow(EF))
    stop(sprintf("Dimension mismatch: Y_test has %d columns but EF has %d rows.",
                 ncol(Y_test), nrow(EF)))
  if (length(EBeta) != ncol(EF))
    stop(sprintf("Dimension mismatch: EF has %d columns but EBeta has length %d.",
                 ncol(EF), length(EBeta)))

  # Apply the same EF column normalization used during training.
  # EF_norms is a K-vector (one norm per factor column); it applies uniformly
  # across all rows, so row-subsetting EF for external validation is safe.
  if (!is.null(EF_norms)) {
    if (length(EF_norms) != ncol(EF))
      stop(sprintf("EF_norms length (%d) must match ncol(EF) (%d).",
                   length(EF_norms), ncol(EF)))
    EF <- sweep(EF, 2, EF_norms, "/")
  }

  ZF_test     <- Y_test %*% EF
  risk_scores <- as.vector(ZF_test %*% EBeta)

  list(ZF_test = ZF_test, risk_scores = risk_scores)
}
