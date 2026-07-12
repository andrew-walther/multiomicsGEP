# =============================================================================
# results/multi_cohort_sim/fit_pca_cox.R
#
# PCA + Cox two-step baseline for Phase 2's joint-vs-2-step value-add
# simulation (docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md).
#
# Stage 1: unsupervised PCA on training Y (survival-blind).
# Stage 2: Cox regression on the top-K principal component scores.
# Held-out scoring re-uses the SAME (training-fitted) rotation and center, so
# train and test projections are on the same scale -- no train/test mismatch.
#
# DEPENDENCIES: survival (CRAN)
# =============================================================================

suppressPackageStartupMessages(library(survival))

#' Fit the PCA+Cox two-step baseline on training data
#'
#' @param Y_train      n x p training genomics matrix
#' @param time_train   n-vector of survival/censoring times
#' @param status_train n-vector of event indicators (1=event, 0=censored)
#' @param K            integer: number of principal components to retain
#' @return Named list:
#'   $EF        -- p x K matrix: gene loadings (prcomp rotation)
#'   $center    -- p-vector: training column means (for held-out centering)
#'   $cox_coef  -- K-vector: Cox coefficients on the PC scores (0 if the Cox
#'                 fit is degenerate/rank-deficient, rather than NA)
#' @export
fit_pca_cox <- function(Y_train, time_train, status_train, K) {
  pca <- prcomp(Y_train, rank. = K, center = TRUE, scale. = FALSE)
  EF  <- pca$rotation                       # p x K

  cox_fit <- tryCatch(
    coxph(Surv(time_train, status_train) ~ pca$x),
    error = function(e) NULL
  )
  cox_coef <- if (!is.null(cox_fit)) {
    cc <- as.numeric(coef(cox_fit))
    cc[!is.finite(cc)] <- 0
    cc
  } else {
    rep(0, K)
  }

  list(EF = EF, center = pca$center, cox_coef = cox_coef)
}

#' Score held-out data with a fitted PCA+Cox model
#'
#' Projects Y_test onto the TRAINING rotation/center (not a fresh PCA on the
#' test set), so train and test scores are on the same scale.
#'
#' @param fit     result of fit_pca_cox()
#' @param Y_test  n_test x p genomics matrix
#' @return Named list:
#'   $scores      -- n_test x K matrix: projected PC scores
#'   $risk_scores -- n_test-vector: scores %*% cox_coef
#' @export
predict_pca_cox <- function(fit, Y_test) {
  scores <- sweep(Y_test, 2, fit$center, "-") %*% fit$EF
  list(scores = scores, risk_scores = as.vector(scores %*% fit$cox_coef))
}
