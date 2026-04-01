# ============================================================
# Script: predict.R
# Purpose: Project new patients into a trained SBMF factor space
#          and compute survival risk scores for hold-out evaluation.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-03-31
# Dependencies: survival (for concordance())
# ============================================================

# ============================================================
# predict_supervised_mf() ----
# ============================================================

#' Project new patients into the learned factor space and compute risk scores.
#'
#' Given a trained SBMF model (factor weights EF and survival coefficients
#' EBeta from training), projects new patients' genomics data Y_test into
#' the K-dimensional latent space via pseudo-inverse projection, then
#' computes Cox linear predictor risk scores.
#'
#' **Projection method:**
#'   L_test = Y_test %*% EF %*% (EF'EF + λI)^{-1}
#'
#' This is the least-squares projection: find L that minimises
#' ||Y_test - L %*% t(EF)||² subject to regularisation λ = 1e-8
#' (for numerical stability when EF columns are near-collinear).
#'
#' **Risk score:**
#'   risk_i = sum_k L_test[i,k] * EBeta[k]
#'
#' Higher risk score → higher predicted hazard → worse prognosis.
#'
#' @param Y_test  numeric matrix (n_test × p): new patients' genomics data.
#'                Must have same number of columns (genes) and same gene
#'                ordering as the training data used to produce EF.
#' @param EF      numeric matrix (p × K): posterior mean factor weights
#'                from a trained SBMF model (fit_supervised_mf_modular()$EF).
#' @param EBeta   numeric vector of length K: posterior mean survival
#'                coefficients (fit_supervised_mf_modular()$EBeta).
#' @param lambda  numeric scalar: ridge regularisation for (EF'EF)^{-1}.
#'                Default 1e-8; increase if EF has near-zero singular values.
#'
#' @return Named list:
#'   $L_test       n_test × K matrix of projected patient loadings
#'   $risk_scores  n_test-vector of Cox linear predictor values (L_test %*% EBeta)
#'
#' @examples
#' # After training:
#' # res <- fit_supervised_mf_modular(Y_train, time_train, status_train, K=5)
#' # pred <- predict_supervised_mf(Y_test, res$EF, res$EBeta)
#' # cindex <- concordance(Surv(time_test, status_test) ~ pred$risk_scores)
predict_supervised_mf <- function(Y_test, EF, EBeta, lambda = 1e-8) {

  # --- Input validation ---
  if (!is.matrix(Y_test) || !is.numeric(Y_test))
    stop("Y_test must be a numeric matrix.")
  if (!is.matrix(EF) || !is.numeric(EF))
    stop("EF must be a numeric matrix.")
  if (!is.numeric(EBeta))
    stop("EBeta must be a numeric vector.")

  p_test <- ncol(Y_test)
  p_train <- nrow(EF)
  K <- ncol(EF)

  if (p_test != p_train)
    stop(sprintf("Dimension mismatch: Y_test has %d columns but EF has %d rows (p).",
                 p_test, p_train))
  if (length(EBeta) != K)
    stop(sprintf("Dimension mismatch: EF has %d columns but EBeta has length %d.",
                 K, length(EBeta)))

  # --- Pseudo-inverse projection ---
  # L_test = Y_test (n_test × p)  %*%  EF (p × K)  %*%  (EF'EF + λI)^{-1} (K × K)
  #
  # The ridge term λI prevents singularity when some factors have near-zero

  # weight columns in EF (e.g., factors driven to zero by the point-normal prior).
  FtF_inv <- solve(crossprod(EF) + lambda * diag(K))
  L_test  <- Y_test %*% EF %*% FtF_inv

  # --- Risk scores ---
  # Cox linear predictor: higher → worse prognosis
  risk_scores <- as.vector(L_test %*% EBeta)

  list(L_test = L_test, risk_scores = risk_scores)
}
