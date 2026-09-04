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
#' @param cohort_id_test NULL (default), or an n_test-vector of cohort labels
#'                 for scoring under COHORT-SPECIFIC survival coefficients
#'                 (fit_cox_on_yf(..., beta_cohort_id = ...); `EBeta` is then
#'                 the fit's K x C `$EBeta` matrix, not a K-vector). A new
#'                 external cohort has no beta^(c) of its own, so:
#'                   - `NULL` (the default, and the only valid choice when
#'                     `EBeta` is an ordinary K-vector): scores with `EBeta`
#'                     directly (unchanged behavior), OR with the fit's
#'                     `$EBeta_pooled` if that is what is passed as `EBeta` --
#'                     this is what "is the shared part what generalizes?"
#'                     means operationally (code/update_beta_cohort.R).
#'                   - supplied: `EBeta` must be the K x C cohort matrix;
#'                     each test patient is scored with THEIR OWN cohort's
#'                     column, matched by name against
#'                     `colnames(EBeta)`/the fit's `$beta_cohort_levels`. An
#'                     unseen level errors loudly (matching the
#'                     `anyNA(strata)` house style in fit_cox_on_yf.R) rather
#'                     than silently falling back to something else -- call
#'                     again with `cohort_id_test = NULL` and
#'                     `EBeta = fit$EBeta_pooled` for that case instead.
#'                 Appended after `EF_norms` (not inserted earlier) so every
#'                 existing positional caller (`select_K.R`, `select_alpha_cv.R`)
#'                 is unaffected.
#'
#' @return Named list:
#'   $ZF_test      n_test × K matrix of normalized projection scores
#'   $risk_scores  n_test-vector of Cox linear predictor values
#'
#' @examples
#' # res <- fit_cox_on_yf(Y_train, time_train, status_train, K=5)
#' # pred <- predict_cox_on_yf(Y_test, res$EF, res$EBeta, res$EF_norms)
#' # concordance(Surv(time_test, status_test) ~ pred$risk_scores)
#' #
#' # Cohort-specific beta, external cohort with no beta^(c) of its own:
#' # pred <- predict_cox_on_yf(Y_ext, res$EF, res$EBeta_pooled, res$EF_norms)
#' #
#' # Cohort-specific beta, scoring an in-sample/matched cohort by name:
#' # pred <- predict_cox_on_yf(Y_test, res$EF, res$EBeta, res$EF_norms,
#' #                            cohort_id_test = my_cohort_labels)
predict_cox_on_yf <- function(Y_test, EF, EBeta, EF_norms = NULL, cohort_id_test = NULL) {

  if (!is.matrix(Y_test) || !is.numeric(Y_test))
    stop("Y_test must be a numeric matrix.")
  if (!is.matrix(EF) || !is.numeric(EF))
    stop("EF must be a numeric matrix.")
  if (!is.numeric(EBeta))
    stop("EBeta must be numeric.")
  if (ncol(Y_test) != nrow(EF))
    stop(sprintf("Dimension mismatch: Y_test has %d columns but EF has %d rows.",
                 ncol(Y_test), nrow(EF)))

  use_cohort <- !is.null(cohort_id_test)
  if (use_cohort) {
    if (!is.matrix(EBeta))
      stop("cohort_id_test was supplied, so EBeta must be a K x C cohort-beta matrix ",
           "(fit_cox_on_yf(..., beta_cohort_id = ...)$EBeta), not a K-vector.")
    if (nrow(EBeta) != ncol(EF))
      stop(sprintf("Dimension mismatch: EF has %d columns but EBeta has %d rows.",
                   ncol(EF), nrow(EBeta)))
    if (length(cohort_id_test) != nrow(Y_test))
      stop(sprintf("cohort_id_test must have length nrow(Y_test) (%d), got %d.",
                   nrow(Y_test), length(cohort_id_test)))
    if (is.null(colnames(EBeta)))
      stop("EBeta must have column names (cohort levels) to match against cohort_id_test.")
    unseen <- setdiff(as.character(cohort_id_test), colnames(EBeta))
    if (length(unseen) > 0)
      stop("cohort_id_test contains level(s) not seen in training: ",
           paste(unseen, collapse = ", "), ". Score with cohort_id_test = NULL and ",
           "EBeta = fit$EBeta_pooled for an unseen cohort instead.")
  } else {
    if (is.matrix(EBeta) && ncol(EBeta) > 1)
      stop("EBeta is a multi-column matrix but cohort_id_test is NULL: pass ",
           "fit$EBeta_pooled (a K-vector) instead, or supply cohort_id_test.")
    if (length(EBeta) != ncol(EF))
      stop(sprintf("Dimension mismatch: EF has %d columns but EBeta has length %d.",
                   ncol(EF), length(EBeta)))
  }

  # Apply the same EF column normalization used during training.
  # EF_norms is a K-vector (one norm per factor column); it applies uniformly
  # across all rows, so row-subsetting EF for external validation is safe.
  if (!is.null(EF_norms)) {
    if (length(EF_norms) != ncol(EF))
      stop(sprintf("EF_norms length (%d) must match ncol(EF) (%d).",
                   length(EF_norms), ncol(EF)))
    EF <- sweep(EF, 2, EF_norms, "/")
  }

  ZF_test <- Y_test %*% EF

  if (use_cohort) {
    cohort_idx_test <- match(as.character(cohort_id_test), colnames(EBeta))
    risk_scores <- rowSums(ZF_test * t(EBeta)[cohort_idx_test, , drop = FALSE])
  } else {
    risk_scores <- as.vector(ZF_test %*% EBeta)
  }

  list(ZF_test = ZF_test, risk_scores = risk_scores)
}
