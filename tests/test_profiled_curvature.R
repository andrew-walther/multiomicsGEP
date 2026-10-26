# ============================================================
# Script: test_profiled_curvature.R
# Purpose: Verify the opt-in profiled Cox quadratic and unchanged default fits.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing test helpers and fitters, ebnm
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/fit_multimodal_yfb.R")
tryCatch(source("code/fit_cox_on_yf.R"),error=function(e){
 if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)
})
#' Read only a named baseline fitter definition from the protected branch
#' @param file,name Repository file and its fitter function name.
#' @return Baseline function closing over an isolated environment.
#' @examples
#' # profile_base_fitter("code/fit_cox_on_yf.R","fit_cox_on_yf")
profile_base_fitter <- function(file,name) {
  text <- system2("git",c("show",paste0("feature/supervised-F:",file)),stdout=TRUE)
  if(!is.null(attr(text,"status")))stop("Cannot read protected baseline")
  env <- new.env(parent=.GlobalEnv)
  for(e in parse(text=text))
    if(is.call(e) && identical(e[[1]],as.name("<-")) &&
       identical(e[[2]],as.name(name)))eval(e,env)
  get(name,env,inherits=FALSE)
}
.profile_base_mm <- profile_base_fitter("code/fit_multimodal_yfb.R","fit_multimodal_yfb")
.profile_base_single <- profile_base_fitter("code/fit_cox_on_yf.R","fit_cox_on_yf")
if(file.exists("code/profile_cox_intercept.R"))source("code/profile_cox_intercept.R")

run_test("ProfileCox-T1: weighted centering preserves contrasts and gives dense profiled curvature", {
  y<-matrix(c(1:5,6,3,5,2,8),5,2);w<-c(.2,1,2,0,3)
  a<-profile_cox_center(y,w)
  H<-diag(w)-tcrossprod(w)/sum(w)
  assert_near(colSums(w*a),c(0,0),tol=1e-12,msg="zero weighted means")
  assert_near(crossprod(a,w*a),crossprod(y,H%*%y),tol=1e-12,msg="full quadratic")
  assert_near(profile_cox_center(sweep(y,2,c(100,-30),"+"),w),a,tol=1e-12,msg="invariant feature offsets")
})
run_test("ProfileCox-T2: zero information and invalid inputs are explicit", {
  y<-matrix(1:8,4,2)
  assert_true(identical(profile_cox_center(y,rep(0,4)),y),msg="zero-curvature path unchanged")
  assert_true(inherits(try(profile_cox_center(y,c(1,-1,1,1)),silent=TRUE),"try-error"))
  assert_true(inherits(try(profile_cox_center(y,1:3),silent=TRUE),"try-error"))
})
run_test("ProfileCox-T3: expected profiled quadratic equals full covariance trace", {
  set.seed(71);Y<-matrix(rnorm(24),8,3);w<-runif(8,.3,2)
  H<-diag(w)-tcrossprod(w)/sum(w)
  m<-c(.2,-.3,.4);v<-c(.04,.02,.01);mb<-.6;vb<-.1
  EZ<-Y%*%m;Cz<-Y%*%diag(v)%*%t(Y)
  direct<-as.numeric(crossprod(mb*EZ,H%*%(mb*EZ)))+
    sum(H*t((mb^2+vb)*Cz+vb*tcrossprod(EZ)))
  YC<-profile_cox_center(Y,w);zc<-YC%*%m
  efficient<-sum(w*((mb^2+vb)*(zc^2+YC^2%*%v)))
  assert_near(direct,efficient,tol=1e-12,msg="all cross-patient covariance retained")
})
run_test("ProfileCox-T4: multimodal trace A and B equal dense profiled scalar update", {
  Y<-list(expression=cbind(c(-3,-1,0,1,2,4),c(1,3,2,5,4,7)))
  rownames(Y$expression)<-paste0("s",1:6);colnames(Y$expression)<-c("a","b")
  L<-matrix(c(1,2,1,3,2,4),6,1);F<-list(expression=matrix(c(.4,-.2),2,1))
  ini<-list(EL=L,EL2=L^2+.1,EF=F,EF2=lapply(F,function(x)x^2+.05),
    EBeta=.5,EBeta2=.35,Tau=list(expression=c(1,1)),
    prior_L=list(list(pi=.5,rate=1,point_mass=FALSE)),
    prior_F=list(expression=list(list(family="normal",variance=1,point_mass=FALSE))))
  t<-setNames(c(1,3,2,5,4,6),rownames(Y$expression));e<-setNames(c(1,1,1,0,1,0),rownames(Y$expression))
  yc<-sweep(Y$expression,2,colMeans(Y$expression-L%*%t(F$expression)))
  z0<-as.vector(yc%*%F$expression);eta0<-z0*.5
  work<-multimodal_yfb_cox_working(t,e,eta0)
  H<-diag(work$w)-tcrossprod(work$w)/sum(work$w)
  h<-work$u+as.vector(H%*%eta0)
  original_F<-multimodal_yfb_update_F_mk;seen_raw_R<-FALSE
  multimodal_yfb_update_F_mk<-function(...) {
    args<-list(...)
    assert_near(args[[5]],yc,tol=1e-12,msg="genomic residual is not survival-centered")
    seen_raw_R<<-TRUE;original_F(...)
  }
  a<-tryCatch(fit_multimodal_yfb(Y,t,e,1L,list(intercept=TRUE,prune=FALSE,
    prior_F=list(expression="normal"),max_outer=1L,max_inner=1L,
    profile_cox_intercept=TRUE,trace=TRUE),init=ini),
    finally={multimodal_yfb_update_F_mk<-original_F})
  assert_true(seen_raw_R,msg="actual multimodal genomic update checked")
  z<-as.vector(yc%*%a$EF$expression)
  vf<-as.vector(a$EF2$expression-a$EF$expression^2)
  s<-a$diagnostics$trace[[1]]$prior_sweeps[[1]]
  aa<-s$canonical_scale[1]
  expectedA<-aa^2*(as.numeric(crossprod(z,H%*%z))+sum(vf*diag(crossprod(yc,H%*%yc))))
  expectedB<-aa*sum(z*h)
  assert_near(s$A,expectedA,tol=1e-9,msg="actual precision includes projected loading covariance")
  assert_near(s$x*s$A,expectedB,tol=1e-9,msg="actual linear term uses profiled working response")
  residual2<-colSums(yc^2-2*yc*(a$EL%*%t(a$EF$expression))+
    a$EL2%*%t(a$EF2$expression))
  rec<-nrow(yc)*sum(log(a$Tau$expression))/2-sum(a$Tau$expression*residual2)/2
  expectedQ<-rec+a$EBeta*sum(h*z)-
    a$EBeta2*(as.numeric(crossprod(z,H%*%z))+sum(vf*diag(crossprod(yc,H%*%yc))))/2
  assert_near(a$diagnostics$history[[1]]$fixed_working_objective,expectedQ,tol=1e-9,
    msg="monitor equals independent dense full-covariance working objective")
})
run_test("ProfileCox-T5: default and explicit FALSE multimodal fits remain identical", {
  set.seed(3);Y<-list(expression=matrix(rnorm(80),20,4))
  dimnames(Y$expression)<-list(paste0("s",1:20),paste0("g",1:4))
  t<-setNames(1:20,rownames(Y$expression));e<-setNames(rep(c(1,0),10),rownames(Y$expression))
  ctrl<-list(prior_F=list(expression="normal"),intercept=TRUE,max_outer=2L,max_inner=1L)
  a<-.profile_base_mm(Y,t,e,2L,ctrl)
  b<-fit_multimodal_yfb(Y,t,e,2L,ctrl)
  c<-fit_multimodal_yfb(Y,t,e,2L,c(ctrl,list(profile_cox_intercept=FALSE)))
  assert_true(identical(a,b)&&identical(a,c),msg="base/default/off object equality")
})
run_test("ProfileCox-T6: unsupported shared-prior combination fails loud", {
  set.seed(4);Y<-list(expression=matrix(rnorm(40),10,4))
  dimnames(Y$expression)<-list(paste0("s",1:10),paste0("g",1:4))
  assert_true(inherits(try(fit_multimodal_yfb(Y,setNames(1:10,rownames(Y$expression)),setNames(rep(1,10),rownames(Y$expression)),2L,
    list(prior_F=list(expression="normal"),profile_cox_intercept=TRUE,
      beta_prior="shared_point_laplace")),silent=TRUE),"try-error"),msg="shared-prior scale question not combined")
})
run_test("ProfileCox-T7: actual single updates use centered survival design and raw genomic residual", {
  set.seed(5);Y<-scale(matrix(rnorm(180),30,6));time<-1:30;e<-rep(c(1,0),15)
  original_beta<-update_beta_k;original_F<-multimodal_yfb_update_F_mk
  seen_beta<-seen_F<-FALSE
  update_beta_k<-function(...) {
    a<-list(...);z<-a[[3]];w<-a[[1]]
    assert_near(sum(w*z),0,tol=1e-10,msg="actual beta covariate weighted centered")
    seen_beta<<-TRUE;original_beta(...)
  }
  multimodal_yfb_update_F_mk<-function(...) {
    a<-list(...);design<-a[[1]];resid<-a[[5]];w<-a[[6]]
    assert_near(colSums(w*design),rep(0,ncol(Y)),tol=1e-10,msg="actual loading survival design")
    assert_near(resid,Y,tol=1e-12,msg="K1 genomics uses original observed matrix")
    seen_F<<-TRUE;original_F(...)
  }
  tryCatch({
    fit_cox_on_yf(Y,time,e,K=1L,max_iter=1L,prior_LF="normal",canonicalize=TRUE,
      w_S=1,profile_cox_intercept=TRUE,verbose=FALSE)
    assert_true(seen_beta&&seen_F,msg="both real updates exercised")
  },finally={
    update_beta_k<-original_beta;multimodal_yfb_update_F_mk<-original_F
  })
})
run_test("ProfileCox-T8: single default/off objects match and unsupported modes fail", {
  set.seed(6);Y<-scale(matrix(rnorm(180),30,6));time<-1:30;e<-rep(c(1,0),15)
  args<-list(Y,time,e,K=2L,max_iter=2L,prior_LF="normal",canonicalize=TRUE,w_S=1,verbose=FALSE)
  a<-do.call(.profile_base_single,args);b<-do.call(fit_cox_on_yf,args)
  c<-do.call(fit_cox_on_yf,c(args,list(profile_cox_intercept=FALSE)))
  assert_true(identical(a,b)&&identical(a,c),msg="base/default/off object equality")
  assert_true(inherits(try(fit_cox_on_yf(Y,time,e,profile_cox_intercept=TRUE,verbose=FALSE),silent=TRUE),"try-error"),
    msg="legacy noncanonical path unsupported")
  assert_true(inherits(try(do.call(fit_cox_on_yf,c(args,list(profile_cox_intercept=TRUE,beta_prior_mode="shared"))),silent=TRUE),"try-error"),
    msg="shared-prior combination unsupported")
})
