# ============================================================
# Script: fit_cox_on_yf.R
# Purpose: Compatibility alias. The YFB fitter (eta = (Y F) beta) was renamed
#          to fit_ssbmf_yfb() in code/fit_ssbmf_yfb.R on 2026-11-04. This file
#          keeps historical runners, frozen results and dated documents that
#          source the old path working. New code should source
#          code/fit_ssbmf_yfb.R and call fit_ssbmf_yfb() directly.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: code/fit_ssbmf_yfb.R
# ============================================================

source("code/fit_ssbmf_yfb.R")

#' Deprecated name for fit_ssbmf_yfb()
#'
#' @param ... Arguments passed unchanged to fit_ssbmf_yfb().
#' @return The fit_ssbmf_yfb() result.
#' @seealso fit_ssbmf_yfb
#' @examples
#' # fit <- fit_cox_on_yf(Y, time, status, K = 5)  # same as fit_ssbmf_yfb()
fit_cox_on_yf <- fit_ssbmf_yfb
