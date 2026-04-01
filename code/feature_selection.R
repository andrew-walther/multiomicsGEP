# ============================================================
# Script: feature_selection.R
# Purpose: Survival-aware gene feature selection for SBMF.
#          Provides univariate Cox p-value filtering as an
#          alternative to variance-based selection.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-01
# Dependencies: survival
# ============================================================

library(survival)

# ============================================================
# cox_feature_selection() ----
# ============================================================

#' Select genes by univariate Cox proportional hazards p-value.
#'
#' For each gene (column of Y), fits a univariate Cox model against
#' the survival outcome and returns indices of genes below the specified
#' p-value threshold.
#'
#' **Important — data leakage prevention:**
#' When used with hold-out evaluation, this function MUST be applied to
#' training data only (Y_train, time_train, status_train).  The selected
#' gene indices are then applied to both Y_train and Y_test.  Fitting on
#' all data would bias the test-set C-index upward by selecting genes
#' already known to be prognostic in the test patients.
#'
#' **Computational note:**
#' Fitting p separate Cox models is O(p × n log n).  For p=5,000 and
#' n~100–300, this takes ~10–60 seconds.  On Longleaf, parallelise with
#' parallel::mclapply().
#'
#' @param Y        numeric matrix (n × p): TRAINING data only
#' @param time     numeric vector (n): survival/censoring times
#' @param status   integer vector (n): event indicators (1 = event, 0 = censored)
#' @param p_thresh numeric: p-value threshold (default 0.05)
#'
#' @return integer vector of column indices (1-based) where the univariate
#'   Cox p-value is below p_thresh.  Returns seq_len(ncol(Y)) (all genes)
#'   if no genes pass the filter.
#'
#' @examples
#' # Fit on training data; apply result to both train and test:
#' # sp <- stratified_split(status_all, test_frac = 0.2, seed = 42)
#' # sel <- cox_feature_selection(Y[sp$train_idx,], time[sp$train_idx],
#' #                               status[sp$train_idx], p_thresh = 0.05)
#' # Y_train_sel <- Y[sp$train_idx, sel]
#' # Y_test_sel  <- Y[sp$test_idx,  sel]
cox_feature_selection <- function(Y, time, status, p_thresh = 0.05) {

  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  if (length(time) != nrow(Y))
    stop("length(time) must equal nrow(Y).")
  if (length(status) != nrow(Y))
    stop("length(status) must equal nrow(Y).")

  p <- ncol(Y)

  # Fit univariate Cox model per gene; capture p-value.
  # tryCatch handles degenerate genes (zero variance, perfect separation).
  pvals <- sapply(seq_len(p), function(j) {
    tryCatch({
      fit    <- coxph(Surv(time, status) ~ Y[, j])
      coefs  <- summary(fit)$coefficients
      # summary()$coefficients has column "Pr(>|z|)" at column 5
      coefs[1, "Pr(>|z|)"]
    }, error = function(e) 1.0)
  })

  selected <- which(pvals < p_thresh)

  if (length(selected) == 0) {
    warning(sprintf(
      "cox_feature_selection: no genes pass p < %.3f; returning all %d genes.",
      p_thresh, p))
    return(seq_len(p))
  }

  selected
}
