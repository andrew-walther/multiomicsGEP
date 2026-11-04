# ============================================================
# Script: predict.R
# Purpose: Compatibility alias. The LB prediction function was renamed
#          to predict_ssbmf_lb() in code/predict_ssbmf_lb.R on 2026-11-04.
#          This file keeps historical runners, frozen results and dated
#          documents that source the old path working. New code should
#          source code/predict_ssbmf_lb.R and call predict_ssbmf_lb()
#          directly.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: code/predict_ssbmf_lb.R
# ============================================================

source("code/predict_ssbmf_lb.R")

#' Deprecated name for predict_ssbmf_lb()
#'
#' @param ... Arguments passed unchanged to predict_ssbmf_lb().
#' @return The predict_ssbmf_lb() result.
#' @seealso predict_ssbmf_lb
#' @examples
#' # pred <- predict_supervised_mf(Y_test, EF, EBeta)  # same as predict_ssbmf_lb()
predict_supervised_mf <- predict_ssbmf_lb
