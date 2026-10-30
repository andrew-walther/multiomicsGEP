# ============================================================
# Script: yfb_joint_beta.R
# Purpose: Exact and mean-field spike/slab reference posteriors for fixed Gaussian designs.
# Created: 2026-10-30
# Dependencies: base R; ashr and yfb_normal_posterior for mean-field comparison
# ============================================================

#' Enumerate the point-normal posterior under a Gaussian working likelihood
#'
#' y|β ~ N(Xβ, I/weight). For included coordinates S, covariance is
#' V_S=(weight X_S'X_S + I/σ²)^-1. All 2^K inclusion patterns are summed.
#' @param X Fixed n-by-K design; K must be at most ten for this reference.
#' @param y Length-n working response.
#' @param weight Positive likelihood multiplier.
#' @param prior_pi Common inclusion probability in [0,1].
#' @param slab_sd Positive common Normal slab SD.
#' @return Means, covariance, variances, inclusion probabilities and pattern weights.
#' @examples
#' yfb_joint_beta(diag(2),c(0,1),1,.5,1)
#' @family yfb_supervised
#' @seealso yfb_meanfield_beta
 yfb_joint_beta <- function(X,y,weight,prior_pi,slab_sd) {
  K<-ncol(X)
  stopifnot(K<=10,weight>0,prior_pi>=0,prior_pi<=1,slab_sd>0)
  masks<-as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),K)))
  H<-weight*crossprod(X);b<-weight*drop(crossprod(X,y))
  means<-matrix(0,nrow(masks),K);covs<-vector("list",nrow(masks));logs<-numeric(nrow(masks))
  for(i in seq_len(nrow(masks))) {
    S<-which(masks[i,]);k<-length(S);V<-matrix(0,K,K)
    logs[i]<-dbinom(k,K,prior_pi,log=TRUE)-lchoose(K,k)
    if(k && is.finite(logs[i])) {
      P<-H[S,S,drop=FALSE]+diag(1/slab_sd^2,k)
      ch<-chol(P);vs<-chol2inv(ch);m<-drop(vs%*%b[S])
      V[S,S]<-vs;means[i,S]<-m
      logs[i]<-logs[i]+.5*sum(b[S]*m)-sum(log(diag(ch)))-k*log(slab_sd)
    }
    covs[[i]]<-V
  }
  weights<-exp(logs-max(logs));weights<-weights/sum(weights)
  mean<-drop(crossprod(weights,means));covariance<-matrix(0,K,K)
  for(i in seq_along(weights)) {
    delta<-means[i,]-mean
    covariance<-covariance+weights[i]*(covs[[i]]+outer(delta,delta))
  }
  list(mean=mean,variance=diag(covariance),covariance=covariance,
    inclusion=drop(crossprod(weights,masks)),weights=weights)
}

#' Fit independent point-normal coordinates at the same Gaussian working target
#' @param X,y,weight,prior_pi,slab_sd As in yfb_joint_beta.
#' @param max_iter,tol Coordinate-sweep cap and absolute mean-change tolerance.
#' @return Mean, variance and explicit convergence status; capped fits retained.
#' @examples
#' yfb_meanfield_beta(diag(2),c(0,1),1,.5,1)
#' @family yfb_supervised
#' @seealso yfb_joint_beta, yfb_normal_posterior
 yfb_meanfield_beta <- function(X,y,weight,prior_pi,slab_sd,max_iter=1000,tol=1e-8) {
  H<-weight*crossprod(X);b<-weight*drop(crossprod(X,y));K<-ncol(X)
  mean<-variance<-numeric(K);converged<-FALSE
  g<-ashr::normalmix(c(1-prior_pi,prior_pi),c(0,0),c(0,slab_sd))
  for(iter in seq_len(max_iter)) {
    old<-mean
    for(k in seq_len(K)) {
      B<-b[k]-sum(H[k,]*mean)+H[k,k]*mean[k]
      q<-yfb_normal_posterior(B/H[k,k],H[k,k],g)
      mean[k]<-q$mean;variance[k]<-q$variance
    }
    if(max(abs(mean-old))<tol) {converged<-TRUE;break}
  }
  list(mean=mean,variance=variance,converged=converged,iterations=iter)
}
