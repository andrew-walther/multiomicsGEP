# ============================================================
# Script: predict_cox_on_yf.R
# Purpose: Compatibility alias. The YFB prediction function was renamed
#          to predict_ssbmf_yfb() in code/predict_ssbmf_yfb.R on 2026-11-04.
#          This file keeps historical runners, frozen results and dated
#          documents that source the old path working. New code should
#          source code/predict_ssbmf_yfb.R and call predict_ssbmf_yfb()
#          directly.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: code/predict_ssbmf_yfb.R
# ============================================================

source("code/predict_ssbmf_yfb.R")

#' Deprecated name for predict_ssbmf_yfb()
#'
#' @param ... Arguments passed unchanged to predict_ssbmf_yfb().
#' @return The predict_ssbmf_yfb() result.
#' @seealso predict_ssbmf_yfb
#' @examples
#' # pred <- predict_cox_on_yf(fit, ...)  # same as predict_ssbmf_yfb()
predict_cox_on_yf <- predict_ssbmf_yfb
