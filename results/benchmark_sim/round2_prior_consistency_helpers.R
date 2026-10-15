# ============================================================
# Script: round2_prior_consistency_helpers.R
# Purpose: Audit expected prior log densities without changing fitted posteriors.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: base R and existing truncated-Normal moment helper
# ============================================================

#' Calculate absolute moments of the fixed point-Laplace posterior
#' @param A,B Conformable finite Gaussian-quadratic precision and signal vectors.
#' @param prior Existing point-Laplace prior list.
#' @return Data frame of posterior slab probabilities and E[|F|].
#' @examples
#' # round2_laplace_abs(4,2,list(family="point_laplace",pi=.5,rate=1))
round2_laplace_abs <- function(A,B,prior) {
  if(!identical(prior$family,"point_laplace") || length(A)!=length(B) ||
    !length(A) || any(!is.finite(A)) || any(A<0) || any(!is.finite(B)))
    stop("Finite matching A>=0,B and point-Laplace prior required.")
  pi <- if(isTRUE(prior$point_mass))0 else prior$pi
  if(!is.finite(pi)||pi<0||pi>1)stop("Invalid prior probability.")
  if(pi==0)return(data.frame(slab_prob=rep(0,length(A)),mean_abs=rep(0,length(A))))
  lambda <- prior$rate
  if(!is.finite(lambda)||lambda<=0)stop("Positive finite Laplace rate required.")
  rows <- lapply(seq_along(A),function(j) {
    if(A[j]==0) {
      if(B[j]!=0)stop("Zero precision with nonzero signal.")
      return(c(slab_prob=pi,mean_abs=pi/lambda))
    }
    x<-B[j]/A[j];s2<-1/A[j];s<-sqrt(s2)
    common<-log(pi)+log(lambda/2)+lambda^2*s2/2
    lw<-c(if(pi<1)log1p(-pi)+dnorm(x,0,s,log=TRUE)else -Inf,
      common-lambda*x+pnorm((x-lambda*s2)/s,log.p=TRUE),
      common+lambda*x+pnorm((-x-lambda*s2)/s,log.p=TRUE))
    ml<-max(lw)+log(sum(exp(lw-max(lw))));p<-exp(lw-ml)
    pos<-multimodal_yfb_positive_truncnorm_moments(x-lambda*s2,s)$mean
    neg<-multimodal_yfb_positive_truncnorm_moments(-x-lambda*s2,s)$mean
    c(slab_prob=if(pi==1)1 else min(1,p[2]+p[3]),mean_abs=p[2]*pos+p[3]*neg)
  })
  as.data.frame(do.call(rbind,rows))
}

#' Evaluate the point-Laplace prior contribution under a fixed posterior
#' @param slab_prob Posterior slab probabilities.
#' @param mean_abs Posterior absolute moments, including the spike.
#' @param prior Point-Laplace prior whose log density is being averaged.
#' @return Summed E_q[log g(F)], including -Inf for incompatible support.
#' @examples
#' # round2_prior_log_expectation(.3,.1,list(family="point_laplace",pi=.5,rate=1))
round2_prior_log_expectation <- function(slab_prob,mean_abs,prior) {
  if(!identical(prior$family,"point_laplace") || length(slab_prob)!=length(mean_abs) ||
    any(!is.finite(slab_prob)) || any(slab_prob<0|slab_prob>1) ||
    any(!is.finite(mean_abs)) || any(mean_abs<0))stop("Invalid posterior moments/family.")
  pi<-if(isTRUE(prior$point_mass))0 else prior$pi
  if(!is.finite(pi)||pi<0||pi>1)stop("Invalid prior probability.")
  if(pi==0)return(if(any(slab_prob>0))-Inf else 0)
  if(pi==1 && any(slab_prob<1))return(-Inf)
  lambda<-prior$rate
  if(!is.finite(lambda)||lambda<=0)stop("Positive finite rate required.")
  ans<-sum(slab_prob)*(log(pi)+log(lambda/2))-lambda*sum(mean_abs)
  if(pi<1)ans<-ans+sum(1-slab_prob)*log1p(-pi)
  ans
}

#' Observe a Single fitting call without changing its loading-update result
#' @param original Existing multimodal_yfb_update_F_mk function.
#' @param sink Callback accepting one aggregate diagnostic row per loading block.
#' @return Wrapper returning the unmodified original result. Caller must expose
#'   the existing Single fitter's iter and k variables.
#' @examples
#' # wrapped <- round2_prior_audit_wrapper(multimodal_yfb_update_F_mk,print)
round2_prior_audit_wrapper <- function(original,sink) {
  force(original);force(sink)
  #' Audit the old/new molecular priors for one actual loading update
  #' @param ... Original loading-update arguments.
  #' @return Original update, unchanged.
  #' @examples
  #' # wrapped(...)
  function(...) {
    frame<-parent.frame();iteration<-get("iter",frame,inherits=FALSE)
    factor<-get("k",frame,inherits=FALSE);args<-list(...)
    old<-if("prior"%in%names(args))args[["prior"]] else args[[14]]
    out<-original(...)
    q<-round2_laplace_abs(out$details$A,out$details$B,old)
    before<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,old)
    after<-round2_prior_log_expectation(q$slab_prob,q$mean_abs,out$prior)
    delta<-after-before
    sink(data.frame(iteration=iteration,factor=factor,
      old_pi=if(isTRUE(old$point_mass))0 else old$pi,old_rate=old$rate,
      new_pi=if(isTRUE(out$prior$point_mass))0 else out$prior$pi,new_rate=out$prior$rate,
      continuous_mass=sum(q$slab_prob),mean_abs_sum=sum(q$mean_abs),
      slab_probability_discrepancy=max(abs(q$slab_prob-out$details$slab_prob)),
      expected_log_prior_old=before,expected_log_prior_new=after,
      prior_step_delta=delta,stored_kl=out$kl,corrected_kl=out$kl-delta))
    out
  }
}
