# ============================================================
# Script: predict_multimodal_yfb.R
# Purpose: Compute frozen raw-projection risk scores from a multimodal YFB fit.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R; code/preprocess_multimodal_yfb.R
# ============================================================

#' Predict frozen multimodal YFB risk scores
#'
#' The prediction path uses only validation molecular measurements and the
#' fitted loading/coefficient moments; it accepts no validation outcomes.
#'
#' @param fit A `multimodal_yfb_fit` object.
#' @param new_Y Named validation expression and methylation matrices.
#' @return List containing named raw-projection `risk_scores` and
#'   `ignored_features` by modality.
#' @examples
#' # predict_multimodal_yfb(fit, Y_validation)
#' @family multimodal_yfb
#' @seealso [fit_multimodal_yfb()]
predict_multimodal_yfb <- function(fit, new_Y) {
  if (!inherits(fit, "multimodal_yfb_fit")) {
    stop("fit must be a multimodal_yfb_fit object.")
  }
  aligned <- align_multimodal_yfb_prediction(new_Y, fit$training_spec)
  projection <- Reduce(`+`, Map(`%*%`, aligned$Y, fit$EF))
  list(
    risk_scores = stats::setNames(
      as.vector(projection %*% fit$EBeta), rownames(aligned$Y[[1]])
    ),
    ignored_features = aligned$ignored_features
  )
}
