# ============================================================
# Script: fit_modular.R
# Purpose: Compatibility alias. The LB fitter (eta = L beta) was renamed to
#          fit_ssbmf_lb() in code/fit_ssbmf_lb.R on 2026-11-04. This file keeps
#          historical runners, frozen results and dated documents that source
#          the old path working. New code should source code/fit_ssbmf_lb.R
#          and call fit_ssbmf_lb() directly.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: code/fit_ssbmf_lb.R
# ============================================================

source("code/fit_ssbmf_lb.R")

#' Deprecated name for fit_ssbmf_lb()
#'
#' @param ... Arguments passed unchanged to fit_ssbmf_lb().
#' @return The fit_ssbmf_lb() result.
#' @seealso fit_ssbmf_lb
#' @examples
#' # fit <- fit_supervised_mf_modular(Y, time, status, K = 5)  # same as fit_ssbmf_lb()
fit_supervised_mf_modular <- fit_ssbmf_lb
