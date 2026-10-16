# ============================================================
# Script: working_cycle_audit_helpers.R
# Purpose: Observe complete Single working cycles without changing updates.
# Created: 2026-10-16
# Dependencies: existing fitting modules and prior-consistency helpers
# ============================================================
source("results/benchmark_sim/round2_prior_consistency_helpers.R")

# Posterior densities and units ----

#' Evaluate E_q log g under the mixed spike/continuous reference measure
#' @param q Posterior density summary, including nonzero mass and absolute moments.
#' @param g Point-Laplace, point-exponential or Normal prior in the same units.
#' @return Summed expected log density; -Inf for incompatible support.
#' @examples
#' # cycle_log_prior(q, list(family="normal",variance=1))
cycle_log_prior <- function(q,g) {
  if(isTRUE(g$point_mass))return(if(any(q$second>0)) -Inf else 0)
  if(g$family=="normal") {
    if(any(q$prob<1))return(-Inf)
    return(sum(-log(2*pi*g$variance)/2-q$second/(2*g$variance)))
  }
  if(g$pi==0)return(if(any(q$second>0)) -Inf else 0)
  if(g$pi==1 && any(q$prob<1))return(-Inf)
  ans<-sum(q$prob)*(log(g$pi)+log(g$rate)-
    if(g$family=="point_laplace")log(2) else 0)-g$rate*sum(q$abs)
  if(g$pi<1)ans<-ans+sum(1-q$prob)*log1p(-g$pi)
  ans
}

#' Retain sufficient density information from an actual coordinate posterior
#' E_q log q = E_q log g + E_q log N(x;theta,s²) - log m(x).
#' @param A,B Actual coordinate precision and signal, including any existing floor.
#' @param g Prior actually used to form this posterior.
#' @param mean,second Actual native returned posterior moments.
#' @param log_ml Actual marginal log normalizer, summed or coordinate-wise.
#' @return Posterior moments, nonzero probabilities, absolute moments and E log q.
#' @examples
#' # cycle_posterior(A,B,g,mean,second,log_ml)
cycle_posterior <- function(A,B,g,mean,second,log_ml) {
  A<-rep_len(A,length(mean));B<-rep_len(B,length(mean))
  if(any(!is.finite(c(A,B,mean,second))) || any(A<0))stop("Invalid density inputs")
  if(isTRUE(g$point_mass)) {
    if(any(second!=0))stop("Point-prior posterior has positive native second moment")
    return(list(mean=mean,second=second,prob=0*mean,abs=0*mean,elogq=0))
  }
  if(g$family=="normal") {
    variance<-second-mean^2
    if(any(variance<=0))stop("Nondegenerate Normal posterior variance must be positive")
    return(list(mean=mean,second=second,prob=rep(1,length(mean)),abs=NA_real_,
      elogq=sum(-log(2*pi*variance)/2-.5)))
  }
  if(g$family=="point_laplace") {
    m<-round2_laplace_abs(A,B,g);prob<-m$slab_prob;absolute<-m$mean_abs
  } else {
    p<-lapply(seq_along(A),function(j)
      multimodal_yfb_point_exponential_posterior(A[j],B[j],g))
    prob<-vapply(p,`[[`,numeric(1),"slab_prob");absolute<-mean
  }
  q<-list(mean=mean,second=second,prob=prob,abs=absolute,elogq=NA_real_)
  elog_g<-cycle_log_prior(q,g)
  if(all(A==0)) {q$elogq<-elog_g;return(q)}
  if(any(A==0))stop("Mixed zero/nonzero information requires coordinate normalizers")
  x<-B/A
  elog_lik<-sum(log(A/(2*pi))/2-(A/2)*(x^2-2*x*mean+second))
  if(any(!is.finite(log_ml)))stop("Missing posterior normalizer")
  q$elogq<-elog_g+elog_lik-sum(log_ml)
  q
}

#' Compute KL against a declared current prior, retaining infinite support conflicts
#' @param q Posterior summary from cycle_posterior().
#' @param g Prior in the same units.
#' @return Scalar KL, possibly Inf; never silently truncated.
#' @examples
#' # cycle_kl(q,g)
cycle_kl <- function(q,g) q$elogq-cycle_log_prior(q,g)

#' Rescale a mixed distribution and its reference prior together
#' @param q Posterior density summary.
#' @param g Reference prior; this does not assert it is carried by the fitter.
#' @param c Positive multiplier on the random variable.
#' @return List q,g in transformed units. Spike entropy has no Jacobian.
#' @examples
#' # cycle_scale(q,g,2)
cycle_scale <- function(q,g,c) {
  stopifnot(length(c)==1L,is.finite(c),c>0)
  q$mean<-q$mean*c;q$second<-q$second*c^2;q$abs<-q$abs*c
  q$elogq<-q$elogq-sum(q$prob)*log(c)
  if(!isTRUE(g$point_mass)) {
    if(g$family=="normal")g$variance<-g$variance*c^2 else g$rate<-g$rate/c
  }
  list(q=q,g=g)
}

#' Convert the actual ebnm prior without adding a sparsity or variance threshold
#' @param g Actual fitted_g object from ebnm.
#' @param family Declared family: normal or point_laplace.
#' @return Prior in the molecular helper convention.
#' @examples
#' # cycle_prior(res$ebnm_result$fitted_g,"point_laplace")
cycle_prior <- function(g,family) {
  scale<-if(family=="normal")g$sd else g$scale
  live<-which(scale>0 & g$pi>0)
  if(!length(live))return(list(family=family,point_mass=TRUE,pi=0))
  if(length(live)!=1L)stop("Expected a single continuous slab")
  if(family=="normal")return(list(family=family,variance=scale[live]^2,point_mass=FALSE))
  list(family=family,pi=sum(g$pi[live]),rate=1/scale[live],point_mass=FALSE)
}

# Independent anchored likelihood ----

#' Evaluate reconstruction and the anchored expected profiled Cox quadratic
#' @param s Full posterior moment state and observed standardized Y.
#' @param work Frozen w,u,eta0,logPL from the actual outer expansion.
#' @return Reconstruction and survival expected log likelihoods.
#' @examples
#' # cycle_likelihood(state,working)
cycle_likelihood <- function(s,work) {
  meanY<-s$EL%*%t(s$EF)
  r2<-s$Y^2-2*s$Y*meanY+s$EL2%*%t(s$EF2)
  rec<-nrow(s$Y)*sum(log(s$Tau))/2-
    sum(sweep(r2,2,s$Tau,"*"))/2-length(s$Y)*log(2*pi)/2
  YC<-profile_cox_center(s$Y,work$w)
  EZ<-YC%*%s$EF;VZ<-YC^2%*%pmax(s$EF2-s$EF^2,0)
  eta<-as.vector(s$Y%*%s$EF%*%s$EBeta)
  variance<-as.vector(VZ%*%s$EBeta2+EZ^2%*%pmax(s$EBeta2-s$EBeta^2,0))
  d<-eta-work$eta0
  surv<-work$logPL+sum(work$u*d)-
    sum(work$w*(profile_cox_center(d,work$w)^2+variance))/2
  list(reconstruction=rec,survival=surv)
}

#' Read an independent copy of the current Single moments from its active call
#' @param frame Existing fitter call environment.
#' @return Moment state, with no modification of the fitter frame.
#' @examples
#' # cycle_snapshot(parent.frame())
cycle_snapshot <- function(frame) {
  mget(c("Y","EL","EL2","EF","EF2","EBeta","EBeta2","Tau"),frame,inherits=FALSE)
}

# Process-local observation ----

#' Create reversible update wrappers and an aggregate transition collector
#' L/beta priors are ephemeral: transformed copies are mathematical references,
#' not falsely presented as implemented carried priors. F priors are carried.
#' @return Environment with rows, q/g density ledger and run(callback).
#' @examples
#' # audit<-cycle_observer(); fit<-audit$run(function() fit_cox_on_yf(...))
cycle_observer <- function() {
  a<-new.env(parent=environment());a$rows<-list();a$work<-NULL
  a$q<-list(L=list(),F=list(),beta=list());a$g<-a$q
  a$iteration<-0L;a$working_calls<-0L
  #' Record a complete reference objective or explicit missing coverage
  #' @param event Transition label.
  #' @param frame Active original fitter frame.
  #' @param s Proposed moment snapshot, leaving the real fitter untouched.
  #' @return Invisible NULL.
  #' @examples
  #' # record("F_prior",frame,state)
  record<-function(event,frame,s) {
    K<-ncol(s$EF)
    complete<-!is.null(a$work) && all(vapply(a$q,function(z)
      length(z)==K && all(vapply(z,function(x)!is.null(x),logical(1))),logical(1)))
    kl<-c(L=NA_real_,F=NA_real_,beta=NA_real_);rec<-surv<-NA_real_
    moment_gap<-NA_real_
    if(!is.null(a$work)) {
      ll<-cycle_likelihood(s,a$work);rec<-ll$reconstruction;surv<-ll$survival
    }
    if(complete)for(b in names(kl))kl[b]<-sum(vapply(seq_len(K),function(k)
      cycle_kl(a$q[[b]][[k]],a$g[[b]][[k]]),numeric(1)))
    if(complete) {
      moment_gap<-max(vapply(seq_len(K),function(k)max(
        abs(a$q$L[[k]]$mean-s$EL[,k]),abs(a$q$L[[k]]$second-s$EL2[,k]),
        abs(a$q$F[[k]]$mean-s$EF[,k]),abs(a$q$F[[k]]$second-s$EF2[,k]),
        abs(a$q$beta[[k]]$mean-s$EBeta[k]),
        abs(a$q$beta[[k]]$second-s$EBeta2[k])),numeric(1)))
      if(moment_gap>1e-12)stop("Density ledger differs from actual posterior moments")
    }
    a$rows[[length(a$rows)+1L]]<-data.frame(iteration=get("iter",frame),
      factor=if(exists("k",frame,inherits=FALSE))get("k",frame) else NA_integer_,
      event=event,complete=complete,reconstruction=rec,survival=surv,
      KL_L=kl[1],KL_F=kl[2],KL_beta=kl[3],reference_objective=rec+surv-sum(kl),
      moment_gap=moment_gap,implemented_common_prior=FALSE,stringsAsFactors=FALSE)
    invisible(NULL)
  }
  #' Wrap one existing function without changing its result
  #' @param original Original function.
  #' @param kind Update family or working/canonical/tau boundary.
  #' @return Wrapper with identical returned object and original arguments.
  #' @examples
  #' # wrapper(update_beta_k,"beta")
  wrapper<-function(original,kind) {
    force(original);force(kind)
    #' Observe the boundary of an original update
    #' @param ... Unmodified original arguments.
    #' @return Unmodified original result.
    #' @examples
    #' # observed(...)
    function(...) {
      frame<-parent.frame();args<-list(...)
      if(!exists("iter",frame,inherits=FALSE))return(original(...))
      s<-cycle_snapshot(frame);out<-original(...)
      if(kind=="working") {
        iteration<-get("iter",frame)
        if(a$iteration!=iteration){a$iteration<-iteration;a$working_calls<-0L}
        a$working_calls<-a$working_calls+1L
        if(a$working_calls==1L) {
          record("pre_refresh",frame,s)
          a$work<-list(w=out$w,u=out$u,eta0=args[[1]],logPL=out$logPL)
          record("refresh",frame,s)
        }
      } else if(kind=="canonical") {
        record("pre_canonical",frame,s)
        for(k in seq_along(out$scale))for(b in names(a$q))
          if(length(a$q[[b]])>=k && !is.null(a$q[[b]][[k]])) {
            z<-cycle_scale(a$q[[b]][[k]],a$g[[b]][[k]],
              if(b=="F")1/out$scale[k] else out$scale[k])
            a$q[[b]][[k]]<-z$q;a$g[[b]][[k]]<-z$g
          }
        for(nm in c("EL","EL2","EF","EF2","EBeta","EBeta2"))s[[nm]]<-out[[nm]]
        record("canonical",frame,s)
      } else if(kind=="tau") {
        record("pre_tau",frame,s);s$Tau<-out$Tau;record("tau",frame,s)
      } else {
        k<-get("k",frame);record(paste0("pre_",kind),frame,s)
        if(kind=="F") {
          g<-if("prior"%in%names(args))args$prior else args[[14]]
          q<-cycle_posterior(out$details$A,out$details$B,g,out$EF,out$EF2,out$details$log_ml)
          a$q$F[[k]]<-q;a$g$F[[k]]<-g
          s$EF[,k]<-out$EF;s$EF2[,k]<-out$EF2
          record("F_posterior",frame,s)
          a$g$F[[k]]<-out$prior;record("F_prior",frame,s)
        } else {
          g<-cycle_prior(out$ebnm_result$fitted_g,if(kind=="L")"point_laplace" else "normal")
          a$g[[kind]][[k]]<-g;record(paste0(kind,"_prior"),frame,s)
          q<-cycle_posterior(out$A,out$B,g,out$mean,out$second,out$ebnm_result$log_likelihood)
          a$q[[kind]][[k]]<-q
          if(kind=="L") {s$EL[,k]<-out$mean;s$EL2[,k]<-out$second} else {
            s$EBeta[k]<-out$mean;s$EBeta2[k]<-out$second
          }
          record(paste0(kind,"_posterior"),frame,s)
        }
      }
      out
    }
  }
  #' Install observers for one callback and restore functions even on failure
  #' @param callback Zero-argument original fit invocation.
  #' @return Original fit object; aggregate rows remain in this environment.
  #' @examples
  #' # audit$run(function() fit_cox_on_yf(...))
  a$run<-function(callback) {
    kinds<-c(update_beta_k="beta",update_L_surv_YFB_k="L",
      multimodal_yfb_update_F_mk="F",update_tau="tau",
      canonicalize_yf_factors="canonical",calc_cox_taylor_yf="working")
    originals<-mget(names(kinds),.GlobalEnv,inherits=FALSE)
    on.exit(list2env(originals,.GlobalEnv),add=TRUE)
    for(nm in names(kinds))assign(nm,wrapper(originals[[nm]],kinds[[nm]]),.GlobalEnv)
    callback()
  }
  a
}
