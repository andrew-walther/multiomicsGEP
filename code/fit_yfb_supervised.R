# ============================================================
# Script: fit_yfb_supervised.R
# Purpose: Opt-in M=1 nonnegative EB factorization with scaled YF Cox supervision.
# Created: 2026-10-30
# Dependencies: ebnm, ashr, yaml, Rcpp; isolated stable exponential kernel
# ============================================================
source("code/compute_elbo.R")
source("code/multimodal_yfb_updates.R")
source("code/select_K.R")
source("code/yfb_exponential_stable.R")
source("code/yfb_diagnostic_controls.R")

# Normal means and Cox quantities ----

#' Evaluate centered Normal or point-normal posterior with separate variance
#' @param x,A Pseudo-observations and positive precisions.
#' @param g Centered normalmix prior with at most one positive slab SD.
#' @return Mean, variance, second moment and coordinate KL vectors.
#' @examples
#' yfb_normal_posterior(1,2,ashr::normalmix(1,0,1))
#' @family yfb_supervised
 yfb_normal_posterior <- function(x,A,g) {
  positive<-which(g$sd>0 & g$pi>0)
  if(!length(positive)) return(list(mean=x*0,variance=x*0,second=x*0,kl=x*0))
  if(length(positive)!=1 || any(g$mean!=0)) stop("Centered single-slab Normal prior required")
  slab<-positive[1];sigma2<-g$sd[slab]^2;prob0<-g$pi[slab]
  v<-1/(A+1/sigma2);m<-v*A*x
  r<-A*sigma2
  logBF<--.5*log1p(r)+.5*m*(A*x)
  prob<-if(prob0==1)rep(1,length(x)) else plogis(qlogis(prob0)+logBF)
  mean<-prob*m;variance<-prob*v+prob*(1-prob)*m^2
  # Disjoint spike/slab supports make the KL a Bernoulli KL plus the
  # slab probability times a Normal KL. No huge Gaussian subtraction.
  bern<-x*0
  ii<-prob>0;bern[ii]<-prob[ii]*(log(prob[ii])-log(prob0))
  ii<-prob<1;bern[ii]<-bern[ii]+(1-prob[ii])*(log1p(-prob[ii])-log1p(-prob0))
  normalKL<-.5*(log1p(r)-r/(1+r)+m^2/sigma2)
  list(mean=mean,variance=variance,second=mean^2+variance,kl=bern+prob*normalKL)
}

#' Solve normal means with stable moments and KL at the same fitted prior
#' @param x Pseudo-observations.
#' @param A Positive precision, scalar or vector.
#' @param family One of point_exponential, point_normal or normal (centered).
#' @param g Optional fixed prior; otherwise estimated by ebnm.
#' @return Posterior mean, separate variance, second moment, fitted prior and KL.
#' @examples
#' yfb_ebnm(c(1,2),2,"point_exponential")
#' @family yfb_supervised
#' @seealso yfb_normal_posterior, yfb_stable_exp
 yfb_ebnm <- function(x,A,family,g=NULL) {
  A<-rep_len(A,length(x))
  if(any(!is.finite(c(x,A))) || any(A<=0)) stop("Invalid normal-means message")
  if(!family %in% c("point_exponential","point_normal","normal")) stop("Unsupported prior")
  cs<-median(1/sqrt(A));xs<-x/cs;As<-A*cs^2;gs<-g
  if(!is.null(gs)) for(field in intersect(c("scale","sd","mean","shift"),names(gs)))
    gs[[field]]<-gs[[field]]/cs
  if(is.null(gs)) gs<-ebnm::ebnm(xs,1/sqrt(As),prior_family=family,mode=0,
                                 output="fitted_g")$fitted_g
  stable<-if(family=="point_exponential") {
    yfb_stable_exp(As,xs*As,gs$pi[2],if(gs$pi[2]==0)1 else 1/gs$scale[2])
  } else yfb_normal_posterior(xs,As,gs)
  if(any(!is.finite(unlist(stable))) || any(stable$variance<0) ||
     sum(stable$kl)< -1e-5*length(x)) stop("Invalid normal-means posterior or KL")
  if(family!="point_exponential" && any(gs$pi[gs$sd>0]>0) && any(stable$variance<=0))
    stop("Posterior variance lost at finite precision")
  gout<-gs
  for(field in intersect(c("scale","sd","mean","shift"),names(gout)))
    gout[[field]]<-gout[[field]]*cs
  mean<-stable$mean*cs;variance<-stable$variance*cs^2
  if(any(!is.finite(c(mean,variance))) ||
     (family!="point_exponential" && any(gs$pi[gs$sd>0]>0) && any(variance<=0)))
    stop("Posterior variance underflow after rescaling")
  list(mean=mean,variance=variance,second=mean^2+variance,g=gout,kl=sum(stable$kl))
}

#' Compute predictor variance without subtracting huge second moments
#' @param Z,V Projection means and variances, n by K.
#' @param beta,vbeta Coefficient means and separate variances, length K.
#' @return Nonnegative n-vector Var(η)=Z² vβ + V (β²+vβ).
#' @examples
#' yfb_predictor_variance(matrix(1,2,1),matrix(.1,2,1),2,.2)
#' @family yfb_supervised
 yfb_predictor_variance <- function(Z,V,beta,vbeta) {
  drop(Z^2%*%vbeta+V%*%(beta^2+vbeta))
}

#' Fit a bounded point-normal prior by marginal likelihood
#'
#' π ∈ [1/K, 0.5], σ ≥ 0.05. The prior is then passed to ebnm with fix_g=TRUE;
#' constraining the marginal fit avoids substituting a different prior after q.
#' @param x,A Normal-means observations and positive precisions.
#' @param K Number of factors, at least two.
#' @return An ashr normalmix prior with a point mass and a Normal slab.
#' @examples
#' yfb_beta_prior(rep(0,7),rep(1,7),7)
#' @family yfb_supervised
 yfb_beta_prior <- function(x,A,K) {
  lower_pi <- 1/K
  objective <- function(par) {
    pi <- par[1]; sigma <- exp(par[2]); se <- 1/sqrt(A)
    a <- log1p(-pi)+dnorm(x,0,se,log=TRUE)
    b <- log(pi)+dnorm(x,0,sqrt(se^2+sigma^2),log=TRUE)
    m <- pmax(a,b)
    -sum(m+log(exp(a-m)+exp(b-m)))
  }
  # K=2 fixes π=.5, so optimize only the slab scale in that case.
  if(K==2) {
    z <- optimize(function(z)objective(c(.5,z)),
                  c(log(.05),log(max(1,abs(x),1/sqrt(A))*100)))
    par <- c(.5,z$minimum)
  } else {
    z <- optim(c(max(lower_pi,.25),log(max(.05,sd(x),na.rm=TRUE))),objective,
      method="L-BFGS-B",lower=c(lower_pi,log(.05)),upper=c(.5,Inf))
    if(z$convergence!=0) {
      # The saved fixture triggers a line-search failure at a constrained
      # boundary. Profile π with a bounded scalar search, then search log σ.
      # σ > max|x| cannot maximize this centered Normal-mixture likelihood:
      # every slab density is decreasing there. Search each local grid basin.
      upper<-max(.05,abs(x))
      profile<-function(logsd) {
        val<-function(pi)objective(c(pi,logsd))
        r<-optimize(val,c(lower_pi,.5),tol=1e-10)
        candidates<-c(lower_pi,r$minimum,.5)
        values<-vapply(candidates,val,0.)
        list(value=min(values),pi=candidates[which.min(values)],logsd=logsd)
      }
      grid<-seq(log(.05),log(upper),length.out=17)
      candidates<-lapply(grid,profile)
      if(upper>.05) for(j in 2:16) {
        vals<-vapply(candidates[seq_len(17)],`[[`,0.,"value")
        if(vals[j]<=vals[j-1] && vals[j]<=vals[j+1]) {
          r<-optimize(function(v)profile(v)$value,grid[c(j-1,j+1)],tol=1e-9)
          candidates[[length(candidates)+1]]<-profile(r$minimum)
        }
      }
      best<-candidates[[which.min(vapply(candidates,`[[`,0.,"value"))]]
      if(!is.finite(best$value)) stop("Bounded β profile likelihood is nonfinite")
      par<-c(best$pi,best$logsd)
      warning("Bounded β optimizer used the profile-likelihood fallback")
    } else par <- z$par
  }
  ashr::normalmix(c(1-par[1],par[1]),c(0,0),c(0,exp(par[2])))
}

#' Evaluate floored log-space Breslow quantities, optionally within strata
#' @param time,event Follow-up and event indicator.
#' @param eta Finite hazard linear predictor.
#' @param strata Optional patient stratum labels.
#' @return List u,w,h,log_partial_likelihood; w ≥ 1e-6.
#' @examples
#' yfb_cox(c(1,2),c(1,0),c(0,0))
#' @family yfb_supervised
 yfb_cox <- function(time,event,eta,strata=NULL) {
  if(is.null(strata)) strata <- rep(1,length(time))
  if(length(strata)!=length(time) || anyNA(strata)) stop("Invalid strata")
  u<-w<-numeric(length(time)); ll<-0
  for(ii in split(seq_along(time),strata)) {
    q<-multimodal_yfb_cox_working_logspace(time[ii],event[ii],eta[ii])
    u[ii]<-q$u;w[ii]<-pmax(q$w,1e-6);ll<-ll+q$log_partial_likelihood
  }
  list(u=u,w=w,h=u+w*eta,log_partial_likelihood=ll)
}

#' Compute scale-invariant normalized projection moments
#' @param Yc Centered subject-by-gene matrix.
#' @param F,F2 Gene-by-factor posterior means and second moments.
#' @param s Optional sweep-fixed scaling constants.
#' @param VF Optional separate loading variance matrix; the fitter always supplies it.
#' @return Z,V (subject-by-factor), and s (factor vector).
#' @examples
#' yfb_projection(matrix(1,3,2),matrix(1,2,2),matrix(1.1,2,2))
#' @family yfb_supervised
 yfb_projection <- function(Yc,F,F2,s=NULL,VF=NULL) {
  # Use the exact norm away from zero. An additive epsilon would violate
  # scale invariance for small, live factors; dead factors use a finite floor.
  if(is.null(s)) {
    norms<-sqrt(colSums(F^2));norms[norms==0]<-1
    s<-1/norms
  }
  if(is.null(VF))VF<-pmax(F2-F^2,0)
  list(Z=Yc%*%sweep(F,2,s,"*"),
       V=Yc^2%*%sweep(VF,2,s^2,"*"),s=s)
}

#' Scale nonnegative features within each cohort without centering
#' @param X Nonnegative patient-by-feature matrix with no missing values.
#' @param cohort Cohort labels, one per patient.
#' @return Scaled matrix; zero-SD columns cause an error, never a silent drop.
#' @examples
#' yfb_scale(matrix(1:12,4),rep(1,4))
#' @family yfb_supervised
 yfb_scale <- function(X,cohort) {
  if(!is.matrix(X) || any(!is.finite(X)) || any(X<0)) stop("Nonnegative finite X required")
  if(length(cohort)!=nrow(X) || anyNA(cohort)) stop("Invalid cohort labels")
  Y<-X
  for(ii in split(seq_len(nrow(X)),cohort)) {
    sdj<-apply(X[ii,,drop=FALSE],2,sd)
    if(any(!is.finite(sdj)) || any(sdj<=0)) stop("Zero or invalid within-cohort feature SD")
    Y[ii,]<-sweep(X[ii,,drop=FALSE],2,sdj,"/")
  }
  Y
}

#' Canonicalize activity/loading moments without changing prediction or priors
#'
#' Live activities have mean E[L²]=1 even when their original RMS is tiny.
#' If an activity is exactly zero but its loading survives through YF,
#' normalize that loading to unit norm: its scale is otherwise unidentified.
#' Both transformations push forward the fitted exponential prior scales.
#' @param st List L,L2,F,F2,gL,gF and optional untouched coefficient fields.
#' @return Same state with canonical moments/priors and orphan factor indices.
#' @examples
#' # yfb_canonicalize(state)
#' @family yfb_supervised
#' @seealso yfb_projection
 yfb_canonicalize <- function(st) {
  if(is.null(st$VL))st$VL<-pmax(st$L2-st$L^2,0)
  if(is.null(st$VF))st$VF<-pmax(st$F2-st$F^2,0)
  rms<-sqrt(colMeans(st$L2));live<-rms>0
  if(any(!live & colSums(abs(st$L))>0)) stop("Nonzero activity with zero second moment")
  cc<-rms;cc[!live]<-1
  st$L<-sweep(st$L,2,cc,"/");st$L2<-sweep(st$L2,2,cc^2,"/")
  st$F<-sweep(st$F,2,cc,"*");st$F2<-sweep(st$F2,2,cc^2,"*")
  st$VL<-sweep(st$VL,2,cc^2,"/");st$VF<-sweep(st$VF,2,cc^2,"*")
  for(k in seq_along(cc)) {
    st$gL[[k]]$scale<-st$gL[[k]]$scale/cc[k]
    st$gF[[k]]$scale<-st$gF[[k]]$scale*cc[k]
  }
  norms<-sqrt(colSums(st$F^2));st$orphan<-which(!live & norms>0)
  for(k in st$orphan) {
    st$F[,k]<-st$F[,k]/norms[k];st$F2[,k]<-st$F2[,k]/norms[k]^2
    st$VF[,k]<-st$VF[,k]/norms[k]^2
    st$gF[[k]]$scale<-st$gF[[k]]$scale/norms[k]
  }
  st$L2<-st$L^2+st$VL;st$F2<-st$F^2+st$VF
  st
}

# Fitting ----

#' Fit opt-in nonnegative factorization with centered normalized YF supervision
#'
#' M=1 implementation of the October 30 specification. γ multiplies Cox
#' precision; ρ damps its F message before the posterior is computed. The
#' returned objective is a diagonal-Taylor approximation, not a guaranteed
#' lower bound. Every decrease is recorded; convergence also requires stable
#' fitted reconstruction and risk. No factor deletion or refitting occurs.
#' @param X Named list containing one nonnegative patient-by-gene matrix.
#' @param time,event Positive follow-up and binary event vector.
#' @param K Initial factor count, at least three for the hybrid start.
#' @param alpha Survival weight in [0,1); used only when gamma is NULL.
#' @param gamma Optional nonnegative fixed Cox multiplier.
#' @param rho Survival-message damping in (0,1].
#' @param cohort Optional within-cohort scaling labels.
#' @param strata Optional Cox risk-set strata.
#' @param max_iter,tol Sweep cap and convergence tolerance; default from globals.
#' @param hybrid Use two Cox-score starting columns (FALSE for its ablation).
#' @param beta_bound Bound the β slab (FALSE for its ablation).
#' @param scaled Whether X is already divided by within-cohort SD.
#' @param tau_init Optional positive length-p starting precision vector.
#' @param tau_fixed Hold tau_init fixed (diagnostic noise-feedback ablation).
#' @param guard Use the fixed-curvature normalized criterion with whole-sweep retries.
#'   This opt-in diagnostic is not a fixed-generative-model ELBO guarantee.
#' @return Posterior moments, priors, projection, risk, per-sweep history,
#'   convergence flag, failure flags and inactive factor indices.
#' @examples
#' # fit_yfb_supervised(list(expression=X),time,event,K=7,gamma=10)
#' @family yfb_supervised
#' @seealso yfb_projection, yfb_beta_prior
fit_yfb_supervised <- function(X,time,event,K=NULL,alpha=NULL,gamma=NULL,rho=.2,
  cohort=NULL,strata=NULL,max_iter=NULL,tol=NULL,hybrid=TRUE,beta_bound=TRUE,scaled=FALSE,
  tau_init=NULL,tau_fixed=FALSE,guard=FALSE) {
  cfg<-yaml::read_yaml("config/globals.yml")
  if(is.null(K)) K<-cfg$benchmark$k_merged_yfb_desurv
  if(is.null(alpha)) alpha<-cfg$benchmark$alpha
  if(is.null(max_iter)) max_iter<-cfg$cavi$max_iter
  if(is.null(tol)) tol<-cfg$cavi$tol
  if(!is.list(X) || length(X)!=1L) stop("Only M=1 is implemented")
  raw<-X[[1]]; n<-nrow(raw);p<-ncol(raw)
  if(is.null(cohort)) cohort<-rep(1,n)
  Y<-if(scaled) raw else yfb_scale(raw,cohort)
  if(any(!is.finite(Y)) || any(Y<0) || K<3 || K>min(dim(Y)) ||
     alpha<0 || alpha>=1 || rho<=0 || rho>1 || max_iter<1)
    stop("Invalid dimensions, nonnegative input or fitting controls")
  yfb_cox(time,event,rep(0,n),strata)
  if(length(time)!=n || length(event)!=n || sum(event)<1) stop("Invalid survival outcomes")
  if(is.null(gamma)) gamma<-alpha/(1-alpha)*p/sum(event)
  if(length(gamma)!=1 || !is.finite(gamma) || gamma<0) stop("Invalid gamma")
  if(tau_fixed && is.null(tau_init))stop("Fixed precision requires tau_init")
  if(!is.null(tau_init) && (length(tau_init)!=p || any(!is.finite(tau_init)) || any(tau_init<=0)))
    stop("Invalid tau_init")
  Yc<-scale(Y,scale=FALSE)
  sv<-svd(Y,nu=K,nv=K)
  F<-abs(sv$v[,seq_len(K),drop=FALSE])
  if(hybrid && gamma>0) {
    score<-drop(crossprod(Yc,yfb_cox(time,event,rep(0,n),strata)$u))
    F[,K-1]<-pmax(score,0);F[,K]<-pmax(-score,0)
  }
  F<-sweep(F,2,pmax(sqrt(colSums(F^2)),1e-10),"/")
  L<-pmax(Y%*%F,1e-6); L2<-L^2+.01;F2<-F^2+.001
  VL<-matrix(.01,n,K);VF<-matrix(.001,p,K)
  tau<-if(is.null(tau_init))rep(1,p) else tau_init; b<-rep(0,K);b2<-rep(.01,K);VBeta<-rep(.01,K)
  klL<-klF<-numeric(K);klB<-0;gL<-gF<-vector("list",K)
  projection<-yfb_projection(Yc,F,F2,VF=VF);Z<-projection$Z;V<-projection$V
  gB<-ashr::normalmix(c(1-1/K,1/K),c(0,0),c(0,1))
  # Brief coefficient warm-up avoids an absorbing β=0 initialization.
  for(warm in if(gamma>0)seq_len(10) else integer()) {
    tay<-yfb_cox(time,event,drop(Z%*%b),strata)
    if(guard) {tay$w<-yfb_cox_curvature_bound(time,event,strata);tay$h<-tay$u+tay$w*drop(Z%*%b)}
    for(k in seq_len(K)) {
      A<-max(sum(tay$w*(Z[,k]^2+V[,k])),1e-10)
      B<-sum(Z[,k]*(tay$h-tay$w*(drop(Z%*%b)-Z[,k]*b[k])))
      q<-yfb_ebnm(B/A,A,"point_normal",gB)
      b[k]<-q$mean;b2[k]<-q$second;VBeta[k]<-q$variance
    }
  }
  hist<-list();converged<-FALSE;flags<-character(); old_obj<-NA_real_
  old_rec<-L%*%t(F);old_eta<-drop(Z%*%b)
  failure<-NULL; last_state<-NULL;trial_trace<-list();rejected_trials<-0L
  tryCatch({
  for(iter in seq_len(max_iter)) {
    last_state<-list(L=L,L2=L2,F=F,F2=F2,tau=tau,b=b,b2=b2,
      gL=gL,gF=gF,gB=gB,Z=Z,V=V,VL=VL,VF=VF,VBeta=VBeta,klL=klL,klF=klF,klB=klB,pr=yfb_projection(Yc,F,F2,VF=VF))
    accepted<-FALSE
    for(trial in seq_len(if(guard)7L else 1L)) {
    if(trial>1L)for(nm in names(last_state))assign(nm,last_state[[nm]])
    rho_step<-rho*if(guard)2^(-(trial-1L)) else 1
    projection<-yfb_projection(Yc,F,F2,VF=VF);Z<-projection$Z;V<-projection$V;s<-projection$s
    eta<-drop(Z%*%b);tay<-yfb_cox(time,event,eta,strata)
    if(guard) {tay$w<-yfb_cox_curvature_bound(time,event,strata);tay$h<-tay$u+tay$w*eta}
    w<-gamma*tay$w;h<-gamma*tay$h
    share<-numeric(K)
    # Fit a bounded common prior to current β messages; retain it through all
    # coordinate updates so q and its KL always refer to exactly the same g.
    if(gamma>0) {
      Ab<-pmax(colSums((Z^2+V)*w),1e-10)
      Bb<-colSums(Z*(h-w*eta))+b*colSums(Z^2*w)
      gB<-if(beta_bound) yfb_beta_prior(Bb/Ab,Ab,K) else
        ebnm::ebnm(Bb/Ab,1/sqrt(Ab),prior_family="point_normal")$fitted_g
    }
    klB<-0
    R<-Y-L%*%t(F)
    for(k in seq_len(K)) {
      Rk<-R+outer(L[,k],F[,k])
      A<-max(sum(tau*F2[,k]),1e-10)
      q<-yfb_ebnm(drop(Rk%*%(tau*F[,k]))/A,A,"point_exponential")
      L[,k]<-q$mean;L2[,k]<-q$second;VL[,k]<-q$variance;gL[[k]]<-q$g;klL[k]<-q$kl
      # Refit g_F on vector pseudo-observations, then hold it fixed during
      # the exact gene-wise sweep. Reconstruction and projection matrices
      # differ: use Rk for the former and s_k Yc for the latter.
      aY<-Yc*s[k];aR<-pmax(tau*sum(L2[,k]),1e-10)
      aS<-rho_step*b2[k]*colSums(aY^2*w)
      hm<-h-w*(drop(Z%*%b)-Z[,k]*b[k])
      br<-tau*drop(crossprod(Rk,L[,k]))
      bs<-rho_step*(b[k]*drop(crossprod(aY,hm))-
        b2[k]*drop(crossprod(aY,w*Z[,k]))+
        b2[k]*colSums(aY^2*w)*F[,k])
      prec<-aR+aS
      priorfit<-yfb_ebnm((br+bs)/prec,prec,"point_exponential")
      gp<-priorfit$g
      prior<-list(family="point_exponential",pi=gp$pi[2],rate=1/gp$scale[2],
                  point_mass=gp$pi[2]==0)
      # Opt-in stable Gauss–Seidel engine: its reconstruction statistics
      # are passed separately from the centered projection matrix.
      yfb_stable_compile()
      z<-yfb_exp_sweep_cpp(aY,tau,sum(L2[,k]),colSums(aY^2*(rho_step*w)),
        drop(crossprod(Rk,L[,k])),drop(crossprod(aY,rho_step*hm)),rho_step*w,
        Z[,k],V[,k],F[,k],F2[,k],VF[,k],b[k],b2[k],prior$pi,
        if(isTRUE(prior$point_mass))1 else prior$rate)
      # One fixed-prior evaluation at the recorded sequential messages
      # checks the same posterior and computes its stable, atomic KL.
      q<-if(isTRUE(prior$point_mass)) priorfit else
        yfb_ebnm(z$x,1/z$s2,"point_exponential",gp)
      if(max(abs(q$mean-z$EF))>1e-5*max(1,abs(z$EF)) ||
         max(abs(q$second-z$EF2))>1e-5*max(1,abs(z$EF2)) ||
         max(abs(q$variance-z$VF))>1e-5*max(1,abs(z$VF)))
        stop("Sequential posterior mismatch")
      F[,k]<-q$mean;F2[,k]<-q$second;VF[,k]<-q$variance;gF[[k]]<-q$g;klF[k]<-q$kl
      pr<-yfb_projection(Yc,F,F2,s,VF=VF);Z<-pr$Z;V<-pr$V
      share[k]<-median(aS/(aR+aS))
      if(gamma>0) {
        A<-max(sum(w*(Z[,k]^2+V[,k])),1e-10)
        B<-sum(Z[,k]*(h-w*(drop(Z%*%b)-Z[,k]*b[k])))
        q<-yfb_ebnm(B/A,A,"point_normal",gB)
        b[k]<-q$mean;b2[k]<-q$second;VBeta[k]<-q$variance;klB<-klB+q$kl
      }
      R<-Y-L%*%t(F)
    }
    # Canonicalize only at the sweep boundary: s stays fixed inside the
    # sweep; the next projection is invariant and β/gβ need no transformation.
    canonical<-yfb_canonicalize(list(L=L,L2=L2,F=F,F2=F2,VL=VL,VF=VF,gL=gL,gF=gF,b=b,b2=b2,VBeta=VBeta,gB=gB))
    L<-canonical$L;L2<-canonical$L2;F<-canonical$F;F2<-canonical$F2
    gL<-canonical$gL;gF<-canonical$gF;VL<-canonical$VL;VF<-canonical$VF
    R<-Y-L%*%t(F)
    Q<-colSums(R^2)+colSums(VL%*%t(F^2)+L^2%*%t(VF)+VL%*%t(VF))
    if(!tau_fixed)tau<-n/pmax(Q,n*1e-8)
    pr<-yfb_projection(Yc,F,F2,VF=VF);Z<-pr$Z;V<-pr$V
    eta<-drop(Z%*%b); now<-yfb_cox(time,event,eta,strata)
    var_eta<-yfb_predictor_variance(Z,V,b,VBeta)
    obj<-.5*n*sum(log(tau/(2*pi)))-.5*sum(tau*Q)-sum(klL+klF)
    if(gamma>0) obj<-obj+gamma*(now$log_partial_likelihood-.5*sum(now$w*var_eta))-klB
    if(guard)obj<-yfb_bound_objective(Y,time,event,list(L=L,F=F,L2=L2,F2=F2,
      VL=VL,VF=VF,b=b,VBeta=VBeta,tau=tau,klL=klL,klF=klF,klB=klB),gamma,strata)
    accepted<-!guard || (is.finite(obj) && (iter==1L || yfb_accept_bound(old_obj,obj)))
    if(guard)trial_trace[[length(trial_trace)+1L]]<-data.frame(iter=iter,trial=trial,
      rho=rho_step,objective=obj,accepted=accepted)
    if(accepted)break
    rejected_trials<-rejected_trials+1L
    }
    if(!accepted) {
      for(nm in names(last_state))assign(nm,last_state[[nm]])
      flags<-c(flags,"backtrack_stalled");break
    }
    decreased<-is.finite(old_obj) && obj<old_obj-1e-8*max(1,abs(old_obj))
    rec<-L%*%t(F)
    change<-max(sqrt(sum((rec-old_rec)^2)/max(1,sum(old_rec^2))),
                sqrt(sum((eta-old_eta)^2)/max(1,sum(old_eta^2))))
    live<-colSums(F^2)>0 & colSums(L2)>0
    hist[[iter]]<-data.frame(iter=iter,objective=obj,objective_decreased=decreased,
      change=change,orphan_count=length(canonical$orphan),share_median=if(any(live))median(share[live]) else NA_real_,
      share_max=if(any(live))max(share[live]) else NA_real_,
      max_beta=max(abs(b)),norm_L=sqrt(sum(L^2)),norm_F=sqrt(sum(F^2)),
      min_tau=min(tau),reconstruction_ratio=sum(R^2)/sum(Y^2),
      expected_reconstruction_ratio=sum(Q)/sum(Y^2),rho_used=rho_step,backtracks=trial-1L)
    for(j in seq_len(K)) {
      hist[[iter]][[paste0("share_",j)]]<-if(live[j])share[j] else NA_real_
      hist[[iter]][[paste0("beta_",j)]]<-b[j]
      hist[[iter]][[paste0("norm_L_",j)]]<-sqrt(sum(L[,j]^2))
      hist[[iter]][[paste0("norm_F_",j)]]<-sqrt(sum(F[,j]^2))
    }
    if(any(!is.finite(c(eta,obj,L,F,b,b2)))) {flags<-c(flags,"nonfinite");break}
    if(all(L==0) || all(F==0)) {flags<-c(flags,"all_zero");break}
    if(iter>cfg$cavi$burn_in && !decreased && is.finite(old_obj) &&
       abs(obj-old_obj)/max(1,abs(old_obj))<tol && change<sqrt(tol)) {
      converged<-TRUE;break
    }
    old_obj<-obj;old_rec<-rec;old_eta<-eta
  }
  },error=function(e) {failure<<-list(message=conditionMessage(e),sweep=iter,
    factor=k,norm_L=sqrt(sum(L^2)),norm_F=sqrt(sum(F^2)),max_beta=max(abs(b)),
    pseudo_range=if(exists("z",inherits=FALSE))range(z$x,na.rm=TRUE) else NULL)})
  if(!is.null(failure)) {
    for(nm in names(last_state)) assign(nm,last_state[[nm]])
    flags<-c(flags,"numerical_failure");converged<-FALSE
  }
  if(!converged && !length(flags)) flags<-"max_iter"
  # Unsupervised factorization followed by an untempered Cox fit; survival
  # never entered its L/F updates. Use identical machinery for fair comparison.
  if(gamma==0) for(warm in seq_len(20)) {
    tay<-yfb_cox(time,event,drop(Z%*%b),strata)
    A<-pmax(colSums((Z^2+V)*tay$w),1e-10)
    B<-colSums(Z*(tay$h-tay$w*drop(Z%*%b)))+b*colSums(Z^2*tay$w)
    gB<-yfb_beta_prior(B/A,A,K)
    for(k in seq_len(K)) {
      aa<-A[k];bb<-sum(Z[,k]*(tay$h-tay$w*(drop(Z%*%b)-Z[,k]*b[k])))
      q<-yfb_ebnm(bb/aa,aa,"point_normal",gB);b[k]<-q$mean;b2[k]<-q$second;VBeta[k]<-q$variance
    }
  }
  history<-if(length(hist))do.call(rbind,hist) else data.frame(iter=0L,
    objective=NA_real_,objective_decreased=FALSE,change=NA_real_,
    share_median=NA_real_,share_max=NA_real_,max_beta=max(abs(b)),
    norm_L=sqrt(sum(L^2)),norm_F=sqrt(sum(F^2)))
  list(EL=L,EL2=L2,EF=setNames(list(F),names(X)),EF2=setNames(list(F2),names(X)),
    EBeta=b,EBeta2=b2,VBeta=VBeta,VL=VL,VF=setNames(list(VF),names(X)),g_beta=gB,g_L=gL,g_F=gF,Tau=tau,s=pr$s,
    Z=Z,VZ=V,eta=drop(Z%*%b),eta_variance=yfb_predictor_variance(Z,V,b,VBeta),history=history,converged=converged,
    flags=flags,failure=failure,orphan=which(colSums(L2)==0 & colSums(F^2)>0),inactive=which(colSums(F^2)==0 | colSums(L^2)==0),gamma=gamma,rho=rho,
    tau_fixed=tau_fixed,guard=guard,rejected_trials=rejected_trials,
    trials=if(length(trial_trace))do.call(rbind,trial_trace) else NULL,
    kl=list(L=klL,F=klF,beta=klB),
    roles=classify_factors(list(EL=L,EF=F,EBeta=b),Y,
      beta_thresh=cfg$k_selection$beta_threshold,pve_thresh=cfg$k_selection$pve_threshold))
}
