# ============================================================
# Script: multimodal_feature_sensitivity.R
# Purpose: Score feature-count studies under the existing multimodal Cox convention.
# Created: 2026-10-15
# Dependencies: survival, existing multimodal Cox helpers
# ============================================================
if(!exists("multimodal_yfb_cox_working"))source("code/multimodal_yfb_helpers.R")

#' Evaluate frozen risk with undefined concordance for constant predictions
#' Uses the model's Breslow partial likelihood, without outcome-based sign flips.
#' @param time,event Aligned positive follow-up times and binary event indicators.
#' @param risk Finite numeric predictor on the same ordered patients.
#' @return List of C, logPL and risk_sd. C remains NA for zero risk variance.
#' @examples
#' multimodal_feature_metrics(1:5,rep(1,5),5:1)
#' @family multimodal_feature_sensitivity
#' @seealso multimodal_yfb_cox_working
multimodal_feature_metrics <- function(time,event,risk) {
  if(length(risk)<2L || length(time)!=length(risk) || length(event)!=length(risk) ||
     any(!is.finite(c(time,event,risk))))
    stop("Finite aligned survival vectors of length at least two required.")
  working <- multimodal_yfb_cox_working(time,event,risk)
  risk_sd <- stats::sd(risk)
  C <- if(risk_sd==0)NA_real_ else as.numeric(survival::concordance(
    survival::Surv(time,event)~risk,reverse=TRUE)$concordance)
  list(C=C,logPL=working$log_partial_likelihood,risk_sd=risk_sd)
}
