# ============================================================
# Script: fit_yfb_supervised.R
# Purpose: Opt-in M=1 nonnegative EB factorization with scaled YF Cox supervision.
# Created: 2026-10-30
# Dependencies: ebnm, ashr, yaml, Rcpp; existing compiled loading sweep
# ============================================================
source("code/compute_elbo.R")
source("code/multimodal_yfb_updates.R")
source("code/select_K.R")

# Normal means and Cox quantities ----

#' Solve normal means and extract the exact KL from the same ebnm result
#' @param x Pseudo-observations.
#' @param A Positive precision, scalar or vector.
#' @param family ebnm prior family.
#' @param g Optional fixed ebnm prior.
#' @return List of posterior mean, second moment, prior, and positive KL.
#' @examples
#' yfb_ebnm(c(1,2), 2, "point_exponential")
#' @family yfb_supervised
#' @seealso compute_ebnm_kl
 yfb_ebnm <- function(x, A, family, g = NULL) {
  A <- rep_len(A,length(x))
  if(any(!is.finite(c(x,A))) || any(A<=0)) stop("Invalid normal-means message")
  # Normal-means equivariance keeps the median pseudo-SE near one. This
  # avoids near-point-mass optimizer errors at precisions above 10^20.
  cs<-median(1/sqrt(A)); xs<-x/cs; As<-A*cs^2
  gs<-g
  if(!is.null(gs)) for(field in intersect(c("scale","sd","mean","shift"),names(gs)))
    gs[[field]]<-gs[[field]]/cs
  args <- list(x=xs,s=1/sqrt(As),prior_family=family,
    output=c("posterior_mean","posterior_second_moment","fitted_g","log_likelihood"))
  if(is.null(gs)) args$mode <- 0 else {args$g_init<-gs;args$fix_g<-TRUE}
  q <- do.call(ebnm::ebnm,args)
  mean <- q$posterior$mean; second <- q$posterior$second_moment
  kl <- -compute_ebnm_kl(q$log_likelihood,As,xs,mean,second)
  if(any(!is.finite(c(mean,second,kl))) || kl < -1e-5*length(x))
    stop("Invalid ebnm posterior or KL")
  gout<-q$fitted_g
  for(field in intersect(c("scale","sd","mean","shift"),names(gout)))
    gout[[field]]<-gout[[field]]*cs
  list(mean=mean*cs,second=second*cs^2,g=gout,kl=kl)
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
    if(z$convergence!=0) stop("Bounded β marginal likelihood did not converge")
    par <- z$par
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
#' @return Z,V (subject-by-factor), and s (factor vector).
#' @examples
#' yfb_projection(matrix(1,3,2),matrix(1,2,2),matrix(1.1,2,2))
#' @family yfb_supervised
 yfb_projection <- function(Yc,F,F2,s=NULL) {
  # Use the exact norm away from zero. An additive epsilon would violate
  # scale invariance for small, live factors; dead factors use a finite floor.
  if(is.null(s)) s<-1/pmax(sqrt(colSums(F^2)),1e-10)
  list(Z=sweep(Yc%*%F,2,s,"*"),
       V=sweep(Yc^2%*%pmax(F2-F^2,0),2,s^2,"*"),s=s)
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
#' @return Posterior moments, priors, projection, risk, per-sweep history,
#'   convergence flag, failure flags and inactive factor indices.
#' @examples
#' # fit_yfb_supervised(list(expression=X),time,event,K=7,gamma=10)
#' @family yfb_supervised
#' @seealso yfb_projection, yfb_beta_prior
fit_yfb_supervised <- function(X,time,event,K=7,alpha=.5,gamma=NULL,rho=.2,
  cohort=NULL,strata=NULL,max_iter=NULL,tol=NULL,hybrid=TRUE,beta_bound=TRUE,scaled=FALSE) {
  cfg<-yaml::read_yaml("config/globals.yml")
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
  Yc<-scale(Y,scale=FALSE)
  sv<-svd(Y,nu=K,nv=K)
  F<-abs(sv$v[,seq_len(K),drop=FALSE])
  if(hybrid && gamma>0) {
    score<-drop(crossprod(Yc,yfb_cox(time,event,rep(0,n),strata)$u))
    F[,K-1]<-pmax(score,0);F[,K]<-pmax(-score,0)
  }
  F<-sweep(F,2,pmax(sqrt(colSums(F^2)),1e-10),"/")
  L<-pmax(Y%*%F,1e-6); L2<-L^2+.01;F2<-F^2+.001
  tau<-rep(1,p); b<-rep(0,K);b2<-rep(.01,K)
  klL<-klF<-numeric(K);klB<-0;gL<-gF<-vector("list",K)
  projection<-yfb_projection(Yc,F,F2);Z<-projection$Z;V<-projection$V
  gB<-ashr::normalmix(c(1-1/K,1/K),c(0,0),c(0,1))
  # Brief coefficient warm-up avoids an absorbing β=0 initialization.
  for(warm in seq_len(10)) {
    tay<-yfb_cox(time,event,drop(Z%*%b),strata)
    for(k in seq_len(K)) {
      A<-max(sum(tay$w*(Z[,k]^2+V[,k])),1e-10)
      B<-sum(Z[,k]*(tay$h-tay$w*(drop(Z%*%b)-Z[,k]*b[k])))
      q<-yfb_ebnm(B/A,A,"point_normal",gB)
      b[k]<-q$mean;b2[k]<-q$second
    }
  }
  hist<-list();converged<-FALSE;flags<-character(); old_obj<-NA_real_
  old_rec<-L%*%t(F);old_eta<-drop(Z%*%b)
  failure<-NULL; last_state<-NULL
  tryCatch({
  for(iter in seq_len(max_iter)) {
    last_state<-list(L=L,L2=L2,F=F,F2=F2,tau=tau,b=b,b2=b2,
      gL=gL,gF=gF,gB=gB,Z=Z,V=V,pr=yfb_projection(Yc,F,F2))
    projection<-yfb_projection(Yc,F,F2);Z<-projection$Z;V<-projection$V;s<-projection$s
    eta<-drop(Z%*%b);tay<-yfb_cox(time,event,eta,strata)
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
      L[,k]<-q$mean;L2[,k]<-q$second;gL[[k]]<-q$g;klL[k]<-q$kl
      # Refit g_F on vector pseudo-observations, then hold it fixed during
      # the exact gene-wise sweep. Reconstruction and projection matrices
      # differ: use Rk for the former and s_k Yc for the latter.
      aY<-Yc*s[k];aR<-pmax(tau*sum(L2[,k]),1e-10)
      aS<-rho*b2[k]*colSums(aY^2*w)
      hm<-h-w*(drop(Z%*%b)-Z[,k]*b[k])
      br<-tau*drop(crossprod(Rk,L[,k]))
      bs<-rho*(b[k]*drop(crossprod(aY,hm))-
        b2[k]*drop(crossprod(aY,w*Z[,k]))+
        b2[k]*colSums(aY^2*w)*F[,k])
      prec<-aR+aS
      priorfit<-yfb_ebnm((br+bs)/prec,prec,"point_exponential")
      gp<-priorfit$g
      prior<-list(family="point_exponential",pi=gp$pi[2],rate=1/gp$scale[2],
                  point_mass=gp$pi[2]==0)
      # Existing compiled Gauss–Seidel engine: its reconstruction statistics
      # are passed separately from the centered projection matrix.
      if(!multimodal_yfb_use_cpp()) stop("Compiled sweep required for this experimental fitter")
      z<-mmyfb_F_sweep_cpp(aY,tau,sum(L2[,k]),colSums(aY^2*(rho*w)),
        drop(crossprod(Rk,L[,k])),drop(crossprod(aY,rho*hm)),rho*w,
        Z[,k],V[,k],F[,k],F2[,k],b[k],b2[k],0L,isTRUE(prior$point_mass),
        prior$pi,prior$rate,NA_real_)
      # One ebnm call at the recorded sequential messages returns the same
      # posterior and gives its marginal likelihood for the atomic KL.
      q<-if(isTRUE(prior$point_mass)) priorfit else
        yfb_ebnm(z$x,1/z$s2,"point_exponential",gp)
      if(max(abs(q$mean-z$EF))>1e-5*max(1,abs(z$EF)) ||
         max(abs(q$second-z$EF2))>1e-5*max(1,abs(z$EF2)))
        stop("Sequential posterior mismatch")
      F[,k]<-q$mean;F2[,k]<-q$second;gF[[k]]<-q$g;klF[k]<-q$kl
      pr<-yfb_projection(Yc,F,F2,s);Z<-pr$Z;V<-pr$V
      share[k]<-median(aS/(aR+aS))
      if(gamma>0) {
        A<-max(sum(w*(Z[,k]^2+V[,k])),1e-10)
        B<-sum(Z[,k]*(h-w*(drop(Z%*%b)-Z[,k]*b[k])))
        q<-yfb_ebnm(B/A,A,"point_normal",gB)
        b[k]<-q$mean;b2[k]<-q$second;klB<-klB+q$kl
      }
      R<-Y-L%*%t(F)
    }
    # Canonicalize only at the sweep boundary: s stays fixed inside the
    # sweep; the next projection is invariant and β/gβ need no transformation.
    cc<-sqrt(colMeans(L2));cc[cc<1e-10]<-1
    L<-sweep(L,2,cc,"/");L2<-sweep(L2,2,cc^2,"/")
    F<-sweep(F,2,cc,"*");F2<-sweep(F2,2,cc^2,"*")
    for(k in seq_len(K)) {
      gL[[k]]$scale<-gL[[k]]$scale/cc[k]
      gF[[k]]$scale<-gF[[k]]$scale*cc[k]
    }
    R<-Y-L%*%t(F)
    Q<-colSums(R^2)+colSums(L2%*%t(F2)-L^2%*%t(F^2))
    tau<-n/pmax(Q,n*1e-8)
    pr<-yfb_projection(Yc,F,F2);Z<-pr$Z;V<-pr$V
    eta<-drop(Z%*%b); now<-yfb_cox(time,event,eta,strata)
    var_eta<-drop((Z^2+V)%*%b2-Z^2%*%(b^2))
    obj<-.5*n*sum(log(tau/(2*pi)))-.5*sum(tau*Q)-sum(klL+klF)
    if(gamma>0) obj<-obj+gamma*(now$log_partial_likelihood-.5*sum(now$w*var_eta))-klB
    decreased<-is.finite(old_obj) && obj<old_obj-1e-8*max(1,abs(old_obj))
    rec<-L%*%t(F)
    change<-max(sqrt(sum((rec-old_rec)^2)/max(1,sum(old_rec^2))),
                sqrt(sum((eta-old_eta)^2)/max(1,sum(old_eta^2))))
    hist[[iter]]<-data.frame(iter=iter,objective=obj,objective_decreased=decreased,
      change=change,share_median=median(share),share_max=max(share),
      max_beta=max(abs(b)),norm_L=sqrt(sum(L^2)),norm_F=sqrt(sum(F^2)))
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
      q<-yfb_ebnm(bb/aa,aa,"point_normal",gB);b[k]<-q$mean;b2[k]<-q$second
    }
  }
  history<-if(length(hist))do.call(rbind,hist) else data.frame(iter=0L,
    objective=NA_real_,objective_decreased=FALSE,change=NA_real_,
    share_median=NA_real_,share_max=NA_real_,max_beta=max(abs(b)),
    norm_L=sqrt(sum(L^2)),norm_F=sqrt(sum(F^2)))
  list(EL=L,EL2=L2,EF=setNames(list(F),names(X)),EF2=setNames(list(F2),names(X)),
    EBeta=b,EBeta2=b2,g_beta=gB,g_L=gL,g_F=gF,Tau=tau,s=pr$s,
    Z=Z,VZ=V,eta=drop(Z%*%b),history=history,converged=converged,
    flags=flags,failure=failure,inactive=which(colSums(F^2)==0 | colSums(L^2)==0),gamma=gamma,rho=rho,
    roles=classify_factors(list(EL=L,EF=F,EBeta=b),Y,
      beta_thresh=cfg$k_selection$beta_threshold,pve_thresh=cfg$k_selection$pve_threshold))
}
