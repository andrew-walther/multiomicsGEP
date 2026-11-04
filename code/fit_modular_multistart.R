# ============================================================
# Script: fit_modular_multistart.R
# Purpose: Compatibility alias. The multi-start wrappers were renamed on
#          2026-11-04 and moved to code/fit_ssbmf_multistart.R:
#            fit_supervised_mf_modular_multistart() -> fit_ssbmf_lb_multistart()
#            fit_cox_on_yf_multistart()             -> fit_ssbmf_yfb_multistart()
#          This file keeps historical runners working. New code should source
#          code/fit_ssbmf_multistart.R directly.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: code/fit_ssbmf_multistart.R (and the fitter it wraps)
# ============================================================

source("code/fit_ssbmf_multistart.R")

#' Deprecated name for fit_ssbmf_lb_multistart()
#'
#' @param ... Arguments passed unchanged to fit_ssbmf_lb_multistart().
#' @return The fit_ssbmf_lb_multistart() result.
#' @seealso fit_ssbmf_lb_multistart
#' @examples
#' # fit_supervised_mf_modular_multistart(Y, time, status, K = 5, N_init = 3)
fit_supervised_mf_modular_multistart <- fit_ssbmf_lb_multistart

#' Deprecated name for fit_ssbmf_yfb_multistart()
#'
#' @param ... Arguments passed unchanged to fit_ssbmf_yfb_multistart().
#' @return The fit_ssbmf_yfb_multistart() result.
#' @seealso fit_ssbmf_yfb_multistart
#' @examples
#' # fit_cox_on_yf_multistart(Y, time, status, K = 5, N_init = 3)
fit_cox_on_yf_multistart <- fit_ssbmf_yfb_multistart
