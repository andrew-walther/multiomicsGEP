# ============================================================
# Script: profile_cox_intercept.R
# Purpose: Profile the arbitrary constant predictor from a frozen Cox working quadratic.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: base R
# ============================================================

#' Center survival covariates under frozen Cox working weights
#'
#' With D=sum(w)>0, Xc=X-1(w'X/D), so Xc' diag(w) Xc equals
#' X' {diag(w)-ww'/D} X. Applying the same transformation to the observed
#' survival design and its loading moments retains all cross-subject covariance
#' needed for the profiled expected quadratic. Reconstruction is not centered
#' by this function. For D=0 the existing zero-information design is returned.
#' @param x Finite numeric vector or subject-by-feature matrix.
#' @param w Conformable finite nonnegative working curvature vector.
#' @return Centered x, preserving dimensions and names; x unchanged if sum(w)=0.
#' @examples
#' profile_cox_center(matrix(1:12,4,3),c(1,2,0,1))
#' @family multimodal_yfb_helpers
#' @seealso [multimodal_yfb_cox_working()]
profile_cox_center <- function(x,w) {
  if(!is.numeric(x) || any(!is.finite(x)) ||
    !is.numeric(w) || any(!is.finite(w)) || any(w<0) ||
    length(w)!=if(is.matrix(x))nrow(x) else length(x))
    stop("Finite conformable x and nonnegative w are required.")
  D <- sum(w)
  if(D==0)return(x)
  if(is.matrix(x))sweep(x,2,colSums(w*x)/D) else x-sum(w*x)/D
}
