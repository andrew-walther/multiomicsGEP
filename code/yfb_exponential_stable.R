# ============================================================
# Script: yfb_exponential_stable.R
# Purpose: Compile and expose the opt-in stable point-exponential calculations.
# Created: 2026-10-30
# Dependencies: Rcpp; base R
# ============================================================
.yfb_stable_cpp <- new.env()

#' Load the experimental numerical kernel once per R session
#' @return Invisibly TRUE; build failure is an error, never a silent fallback.
#' @examples
#' yfb_stable_compile()
#' @family yfb_supervised
 yfb_stable_compile <- function() {
  if(isTRUE(.yfb_stable_cpp$ready)) return(invisible(TRUE))
  old<-Sys.getenv("PKG_CXXFLAGS");on.exit(Sys.setenv(PKG_CXXFLAGS=old))
  if(Sys.info()[["sysname"]]=="Darwin") {
    sdk<-system2("xcrun","--show-sdk-path",stdout=TRUE)
    Sys.setenv(PKG_CXXFLAGS=paste(old,paste0("-I",sdk,"/usr/include/c++/v1 -isysroot ",sdk)))
  }
  Rcpp::sourceCpp("code/yfb_exponential_stable.cpp",env=.GlobalEnv)
  .yfb_stable_cpp$ready<-TRUE
  invisible(TRUE)
}

#' Evaluate exact point-exponential posterior moments and KL stably
#'
#' Cancel the common N(x;0,1/A) term in the spike/slab odds. With
#' I=∫₀∞ λ exp((B−λ)θ−Aθ²/2)dθ and N=1−π+πI,
#' KL(q||g)=B E[θ]−A E[θ²]/2−log N. This is the same EBNM posterior,
#' evaluated without subtracting huge Gaussian log likelihoods.
#' @param A,B Positive quadratic precisions and linear coefficients.
#' @param pi Slab proportion in [0,1].
#' @param rate Positive exponential rate; ignored for a point mass.
#' @return List with mean, second and coordinate KL vectors.
#' @examples
#' yfb_stable_exp(3,2,.3,2)
#' @family yfb_supervised
 yfb_stable_exp <- function(A,B,pi,rate) {
  A<-rep_len(A,length(B))
  yfb_stable_compile()
  yfb_exp_vector_cpp(A,B,pi,if(pi==0)1 else rate)
}
