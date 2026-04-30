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
#' Given a trained Cluster B SBMF model (factor weights EF and survival
#' coefficients EBeta from training), computes the observed projection scores
#' ZF_test = Y_test · EF and then the Cox linear predictor.
#'
#' **Projection (Cluster B):**
#'   ZF_test = Y_test %*% EF
#'
#' This is the same formula used during training (eta = ZF · beta_tilde).
#' No pseudo-inverse is needed: training and test both use the direct
#' Y · EF projection, eliminating the train/test formula mismatch that
#' motivated the Cox-on-YF reformulation.
#'
#' **Risk score:**
#'   risk_i = sum_k ZF_test[i,k] * EBeta[k]
#'          = (Y_test %*% EF %*% EBeta)[i]
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
#'
#' @return Named list:
#'   $L_test       n_test × K matrix of observed projection scores ZF_test = Y_test %*% EF
#'   $risk_scores  n_test-vector of Cox linear predictor values (ZF_test %*% EBeta)
#'
#' @examples
#' # After training:
#' # res <- fit_supervised_mf_modular(Y_train, time_train, status_train, K=5)
#' # pred <- predict_supervised_mf(Y_test, res$EF, res$EBeta)
#' # cindex <- concordance(Surv(time_test, status_test) ~ pred$risk_scores)
predict_supervised_mf <- function(Y_test, EF, EBeta) {

  # --- Input validation ---
  if (!is.matrix(Y_test) || !is.numeric(Y_test))
    stop("Y_test must be a numeric matrix.")
  if (!is.matrix(EF) || !is.numeric(EF))
    stop("EF must be a numeric matrix.")
  if (!is.numeric(EBeta))
    stop("EBeta must be a numeric vector.")

  p_test  <- ncol(Y_test)
  p_train <- nrow(EF)
  K       <- ncol(EF)

  if (p_test != p_train)
    stop(sprintf("Dimension mismatch: Y_test has %d columns but EF has %d rows (p).",
                 p_test, p_train))
  if (length(EBeta) != K)
    stop(sprintf("Dimension mismatch: EF has %d columns but EBeta has length %d.",
                 K, length(EBeta)))

  # --- Direct projection (Cluster B: same formula as training) ---
  # ZF_test = Y_test · EF  (n_test × K observed projection scores)
  # No pseudo-inverse needed: eta_test = ZF_test · beta_tilde matches training.
  L_test <- Y_test %*% EF

  # --- Risk scores ---
  # Cox linear predictor: higher → worse prognosis
  risk_scores <- as.vector(L_test %*% EBeta)

  list(L_test = L_test, risk_scores = risk_scores)
}
