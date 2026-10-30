# ============================================================
# Script: yfb_diagnostic_controls.R
# Purpose: Fixed risk-set curvature and normalized moment criterion for guarded diagnostics.
# Created: 2026-10-30
# Dependencies: base R; helpers in fit_yfb_supervised.R
# ============================================================

#' Bound the negative Cox Hessian by a fixed diagonal matrix
#'
#' Each event's log-sum-exp Hessian is bounded by I_R/2. Thus
#' D_i = (1/2) sum_events 1{i in R_event} bounds the total Hessian.
#' @param time,event Follow-up times and binary event indicators.
#' @param strata Optional risk-stratum labels.
#' @return Length-n curvature bounds, floored at 10^-6.
#' @examples
#' yfb_cox_curvature_bound(c(1,2),c(1,0))
#' @family yfb_supervised
#' @seealso yfb_bound_objective
 yfb_cox_curvature_bound <- function(time,event,strata=NULL) {
  if(is.null(strata))strata<-rep(1,length(time))
  D<-numeric(length(time))
  for(ii in split(seq_along(time),strata))
    D[ii]<-.5*findInterval(time[ii],sort(time[ii][event[ii]==1]))
  pmax(D,1e-6)
}

#' Evaluate one fixed normalized moment criterion, including matching KL
#'
#' At given normalization s, Jensen/Taylor gives E[log PL] >=
#' log PL(Eη) - (1/2) sum_i D_i Var(η_i). Here s is recomputed from
#' posterior loading means at each state. This explicit q-normalized criterion
#' is not claimed to be a fixed-generative-model ELBO under that normalization.
#' @param Y Nonnegative n-by-p prepared expression matrix.
#' @param time,event,strata Survival outcomes and optional strata.
#' @param st List L,F,L2,F2,VL,VF,b,VBeta,tau,klL,klF,klB.
#' @param gamma Fixed survival likelihood weight.
#' @return Scalar normalized moment criterion, larger is better.
#' @examples
#' # Evaluated internally for every guarded proposal.
#' @family yfb_supervised
#' @seealso yfb_cox_curvature_bound, yfb_accept_bound
 yfb_bound_objective <- function(Y,time,event,st,gamma,strata=NULL) {
  R<-Y-st$L%*%t(st$F)
  Q<-colSums(R^2)+colSums(st$VL%*%t(st$F^2)+st$L^2%*%t(st$VF)+st$VL%*%t(st$VF))
  value<-.5*nrow(Y)*sum(log(st$tau/(2*pi)))-.5*sum(st$tau*Q)-sum(st$klL+st$klF)
  if(gamma>0) {
    p<-yfb_projection(scale(Y,scale=FALSE),st$F,st$F2,VF=st$VF)
    eta<-drop(p$Z%*%st$b)
    v<-yfb_predictor_variance(p$Z,p$V,st$b,st$VBeta)
    ll<-yfb_cox(time,event,eta,strata)$log_partial_likelihood
    value<-value+gamma*(ll-.5*sum(yfb_cox_curvature_bound(time,event,strata)*v))-st$klB
  }
  value
}

#' Accept only finite proposals that preserve the same normalized criterion
#' @param old,new Previous and proposed scalar criterion values.
#' @param tol Relative numerical acceptance tolerance.
#' @return Logical scalar; nonfinite proposals are rejected.
#' @examples
#' yfb_accept_bound(1,.9)
#' @family yfb_supervised
#' @seealso yfb_bound_objective
 yfb_accept_bound <- function(old,new,tol=1e-8) {
  is.finite(new) && new>=old-tol*max(1,abs(old))
}
