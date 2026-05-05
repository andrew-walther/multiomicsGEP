# ============================================================
# Script: predict_cox_on_yf.R
# Purpose: Cluster B prediction — direct Y·EF projection (same formula as training)
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-04
# Dependencies: none (base R only)
# ============================================================

#' Cluster B prediction: project test patients via Y_test · EF · beta_tilde
#'
#' Under the Cox-on-YF reformulation (eta = ZF · beta_tilde where ZF = Y·EF),
#' prediction on new data uses the SAME formula as training: ZF_test = Y_test·EF.
#' This eliminates the train/test mismatch present in Cluster A (which uses an
#' OLS projection at test time instead of EBNM-shrunk loadings from training).
#'
#' @param Y_test  n_test × p numeric matrix: test patients' genomics data.
#'                Must have the same columns (genes) and ordering as training data.
#' @param EF      p × K numeric matrix: posterior mean factor weights from
#'                a trained Cox-on-YF model (fit_cox_on_yf()$EF).
#' @param EBeta   K-vector: posterior mean survival coefficients (fit_cox_on_yf()$EBeta).
#'
#' @return Named list:
#'   $ZF_test      n_test × K matrix of observed projection scores Y_test · EF
#'   $risk_scores  n_test-vector of Cox linear predictor values (ZF_test · beta_tilde)
#'
#' @examples
#' # res <- fit_cox_on_yf(Y_train, time_train, status_train, K=5)
#' # pred <- predict_cox_on_yf(Y_test, res$EF, res$EBeta)
#' # concordance(Surv(time_test, status_test) ~ pred$risk_scores)
predict_cox_on_yf <- function(Y_test, EF, EBeta) {

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

  # Cluster B: direct observed projection — same formula as training
  ZF_test     <- Y_test %*% EF
  risk_scores <- as.vector(ZF_test %*% EBeta)

  list(ZF_test = ZF_test, risk_scores = risk_scores)
}
