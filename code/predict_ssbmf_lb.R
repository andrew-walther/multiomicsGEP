# ============================================================
# Script: predict_ssbmf_lb.R
# Purpose: Project new patients into a trained SBMF factor space
#          and compute survival risk scores for hold-out evaluation.
#          predict_ssbmf_lb() uses least-squares projection;
#          predict_ssbmf_lb_foldin() (2026-11-06) uses EBNM fold-in.
# Author: Andrew Walther
# Created: 2026-03-31
# Dependencies: survival (for concordance()), ebnm (fold-in)
# ============================================================

# ============================================================
# predict_ssbmf_lb() ----
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
#'                from a trained SBMF model (fit_ssbmf_lb()$EF).
#' @param EBeta   numeric vector of length K: posterior mean survival
#'                coefficients (fit_ssbmf_lb()$EBeta).
#' @param lambda  numeric scalar: ridge regularisation for (EF'EF)^{-1}.
#'                Default 1e-8; increase if EF has near-zero singular values.
#'
#' @return Named list:
#'   $L_test       n_test × K matrix of projected patient loadings
#'   $risk_scores  n_test-vector of Cox linear predictor values (L_test %*% EBeta)
#'
#' @examples
#' # After training:
#' # res <- fit_ssbmf_lb(Y_train, time_train, status_train, K=5)
#' # pred <- predict_ssbmf_lb(Y_test, res$EF, res$EBeta)
#' # cindex <- concordance(Surv(time_test, status_test) ~ pred$risk_scores, reverse = TRUE)
predict_ssbmf_lb <- function(Y_test, EF, EBeta, lambda = 1e-8) {

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

  # --- Pseudo-inverse projection via SVD ---
  # L_test = Y_test (n_test × p)  %*%  EF (p × K)  %*%  pinv(EF'EF)
  #
  # We use the Moore-Penrose pseudoinverse (SVD-based) rather than solve().
  # ARD can drive some factor columns to near-zero, making EF'EF near-singular
  # even with a fixed ridge term.  SVD thresholds tiny singular values at
  # max(d) * lambda (relative), so collapsed factors contribute zero to L_test
  # without causing numerical errors.  Collapsed factors also have EBeta[k] ≈ 0
  # by the same prior shrinkage, so their contribution to risk_scores is correct.
  FtF  <- crossprod(EF)
  sv   <- svd(FtF)
  # Zero out singular values below lambda * max(d): treats collapsed factors as inert
  d_inv <- ifelse(sv$d > lambda * max(sv$d), 1 / sv$d, 0)
  FtF_pinv <- sv$v %*% diag(d_inv, nrow = K, ncol = K) %*% t(sv$u)

  L_test      <- Y_test %*% EF %*% FtF_pinv

  # --- Risk scores ---
  # Cox linear predictor: higher → worse prognosis
  risk_scores <- as.vector(L_test %*% EBeta)

  list(L_test = L_test, risk_scores = risk_scores)
}

# ============================================================
# predict_ssbmf_lb_foldin() ----
# ============================================================

#' Fold-in prediction of new-patient loadings for the LB model.
#'
#' Computes each new patient's loadings as the training L update with the
#' survival term removed (the patient's outcome is unknown), holding the
#' training factors, noise precisions and fitted loading priors fixed
#' (Yusha's note, docs/notes/2026-10-30-yusha-MF_NewPatient_FoldIn.pdf).
#' For factor k and patient i,
#'
#'   A*_k  = (1 - alpha) sum_j tau_j E[f_jk^2]
#'   B*_ik = (1 - alpha) sum_j tau_j fbar_jk (y*_ij - sum_{m != k} l*_im fbar_jm)
#'   q(l*_ik) = EBNM(x = B*_ik / A*_k, s = 1 / sqrt(A*_k); g = g_{L,k} fixed).
#'
#' The (1 - alpha) factor matches the tempered genomics precision of the
#' training update, so new and training loadings are shrunk on the same
#' scale and EBeta applies to both. It cancels in x but not in s. Sweeps over
#' k repeat (from l* = 0) until the largest change in E[l*] is below tol.
#'
#' Compared with predict_ssbmf_lb()'s least-squares projection, fold-in
#' applies the training prior (nonnegativity and sparsity for
#' point_exponential) and gives posterior second moments.
#'
#' @param Y_test     numeric matrix (n_test x p), same genes, order and
#'                   preprocessing as the training data; no missing values.
#'                   If both Y_test and EF have gene names they must match.
#' @param fit        fit_ssbmf_lb() result with $EF, $EF2, $Tau, $EBeta and
#'                   $foldin (fits made before 2026-11-06 lack $foldin).
#' @param tol        relative convergence tolerance: stop when
#'                   max |change in E[l*]| < tol * max(1, max |E[l*]|).
#' @param max_sweeps maximum number of sweeps over the K factors.
#' @return Named list:
#'   $L_test       n_test x K posterior means E[l*]
#'   $L2_test      n_test x K posterior second moments E[l*^2]
#'   $risk_scores  n_test-vector L_test %*% EBeta (larger = higher hazard)
#'   $n_sweeps     sweeps used
#'   $converged    TRUE if the last sweep met the relative tolerance
#' @examples
#' # fit  <- fit_ssbmf_lb(Y_train, time_train, status_train, K = 5)
#' # pred <- predict_ssbmf_lb_foldin(Y_test, fit)
#' # concordance(Surv(time_test, status_test) ~ pred$risk_scores, reverse = TRUE)
#' @seealso predict_ssbmf_lb, update_L_k
predict_ssbmf_lb_foldin <- function(Y_test, fit, tol = 1e-6, max_sweeps = 200) {

  ## Input checks (fail loud) ----
  if (!is.matrix(Y_test) || !is.numeric(Y_test))
    stop("Y_test must be a numeric matrix.")
  if (is.null(fit$foldin))
    stop("fit has no $foldin field; refit with the current fit_ssbmf_lb().")
  if (!is.null(fit$cohort_id))
    stop("Fold-in does not support cohort-offset fits (cohort_id): a new ",
         "patient's cohort offset is not defined.")
  if (anyNA(Y_test))
    stop("Y_test has missing values; fold-in would treat them as zero signal.")
  EF <- fit$EF; EF2 <- fit$EF2; Tau <- fit$Tau
  K  <- ncol(EF); n_test <- nrow(Y_test)
  if (ncol(Y_test) != nrow(EF))
    stop(sprintf("Dimension mismatch: Y_test has %d columns but EF has %d rows.",
                 ncol(Y_test), nrow(EF)))
  if (!is.null(colnames(Y_test)) && !is.null(rownames(EF)) &&
      !identical(colnames(Y_test), rownames(EF)))
    stop("Gene names of Y_test do not match rownames(fit$EF) in order.")
  g_L <- fit$foldin$g_L; prior_LF <- fit$foldin$prior_LF; alpha <- fit$foldin$alpha
  if (alpha >= 1)
    stop("Fold-in needs alpha < 1: at alpha = 1 the loadings ignore genomics, ",
         "so every new patient would get the prior mean.")

  ## Precision per factor (scalar, the same for every new patient) ----
  # Mirrors update_L_k(): genomics term only, tempered by (1 - alpha),
  # with the same A floor and s clamps.
  A <- pmax((1 - alpha) * colSums(Tau * EF2), 1e-10)
  s <- pmax(1 / sqrt(A), 1e-8)

  L <- matrix(0, n_test, K); L2 <- matrix(0, n_test, K)
  R <- Y_test                       # full residual Y* - L* EF' (L* = 0)
  converged <- FALSE

  ## Coordinate sweeps over factors ----
  for (sweep in seq_len(max_sweeps)) {
    L_old <- L
    for (k in seq_len(K)) {
      R_k <- R + outer(L[, k], EF[, k])                        # R^{-k}
      B   <- (1 - alpha) * as.vector(R_k %*% (Tau * EF[, k]))
      x   <- B / A[k]
      x[!is.finite(x)] <- 0
      res <- ebnm::ebnm(x = x, s = rep(s[k], n_test), prior_family = prior_LF,
                        g_init = g_L[[k]], fix_g = TRUE)
      L[, k]  <- res$posterior$mean
      L2[, k] <- res$posterior$sd^2 + res$posterior$mean^2
      R <- R_k - outer(L[, k], EF[, k])
    }
    if (max(abs(L - L_old)) < tol * max(1, max(abs(L)))) { converged <- TRUE; break }
  }
  if (!converged)
    warning(sprintf("Fold-in did not converge in %d sweeps (last change %.2e > tol %.2e).",
                    max_sweeps, max(abs(L - L_old)), tol))

  list(L_test = L, L2_test = L2, risk_scores = as.vector(L %*% fit$EBeta),
       n_sweeps = sweep, converged = converged)
}
