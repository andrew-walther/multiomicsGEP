# ============================================================
# Script: cox_working_logspace.R
# Purpose: Evaluate the same Breslow working quantities using stable log-space risk sums.
# Created: 2026-10-15
# Dependencies: base R
# ============================================================

#' Compute Breslow working quantities without globally shifted risk weights
#' log R(t) is a suffix log-sum-exp of η. With d_t deaths at time t,
#' log H_i = log sum_{t<=t_i} d_t exp(-log R(t)); w_i=exp(η_i+log H_i).
#' This changes arithmetic, not the diagonal working approximation:
#' u=δ-w, h=u+wη. There is no clipping, added weight or changed likelihood.
#' @param time Positive finite follow-up times.
#' @param event Numeric event indicators, zero or one.
#' @param eta Finite Cox predictors, in the same patient order.
#' @return List u, w, h and Breslow log_partial_likelihood, in original order.
#' @examples
#' multimodal_yfb_cox_working_logspace(c(1,2),c(1,1),c(1000,-1000))
#' @family multimodal_yfb_helpers
#' @seealso multimodal_yfb_cox_working
multimodal_yfb_cox_working_logspace <- function(time,event,eta) {
 n<-length(time)
 if(!is.numeric(time) || !is.numeric(event) || !is.numeric(eta) ||
    length(event)!=n || length(eta)!=n || n==0L ||
    any(!is.finite(time)) || any(time<=0) || any(!is.finite(eta)) ||
    any(!is.finite(event)) || any(!event%in%c(0,1)))
   stop("Positive finite time, binary event and finite conformable eta are required.")
 ord<-order(time);t<-time[ord];d<-event[ord];e<-eta[ord]
 log_risk<-numeric(n);log_risk[n]<-e[n]
 if(n>1L)for(i in rev(seq_len(n-1L))) {
   high<-max(e[i],log_risk[i+1L])
   log_risk[i]<-high+log1p(exp(-abs(e[i]-log_risk[i+1L])))
 }
 first<-!duplicated(t);group<-cumsum(first)
 log_risk<-log_risk[match(t,t)]
 deaths<-as.numeric(ave(d,t,FUN=sum))[first]
 increment<-rep(-Inf,length(deaths));positive<-deaths>0
 increment[positive]<-log(deaths[positive])-log_risk[first][positive]
 log_hazard<-increment
 if(length(increment)>1L)for(g in 2:length(increment)) {
   a<-log_hazard[g-1L];b<-increment[g]
   log_hazard[g]<-if(a==-Inf)b else if(b==-Inf)a else
     max(a,b)+log1p(exp(-abs(a-b)))
 }
 w_sorted<-exp(e+log_hazard[group])
 w<-u<-numeric(n);w[ord]<-w_sorted;u[ord]<-d-w_sorted
 # Select event rows before multiplication: censored rows contribute exactly0.
 log_pl<-sum((e-log_risk)[d==1])
 h<-u+w*eta
 if(any(!is.finite(c(w,u,h,log_pl))))
   stop("Cox log-space arithmetic exceeded finite numeric representation.")
 list(u=u,w=w,h=h,log_partial_likelihood=log_pl)
}
