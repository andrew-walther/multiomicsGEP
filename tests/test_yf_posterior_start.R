# ============================================================
# Script: test_yf_posterior_start.R
# Purpose: Check the approved opt-in preservation of pruning posterior states.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing fitter/test helpers, ebnm
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(!exists("fit_cox_on_yf"))suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"),
  error=function(e) {if(!grepl("!is.null(real_Y)",conditionMessage(e),fixed=TRUE))stop(e)}))
.ps_y <- matrix(seq_len(240)/240,40,6)+matrix(sin(seq_len(240)),40,6)
.ps_a <- list(K=3L,max_iter=6L,canonicalize=TRUE,w_S=1,prior_LF="point_laplace",
  init_supervised=TRUE,beta_prior_mode="shared",prior_beta="point_laplace",
  beta_warmup=10L,verbose=FALSE)
run_test("YFPosterior-T1: state extraction keeps all moments and persistent priors", {
  f <- do.call(fit_cox_on_yf,c(list(.ps_y,1:40,rep(1,40)),.ps_a,
    list(prune_preserve_posterior=TRUE)))
  s <- f$posterior_state;k <- c(3L,1L)
  x <- yf_subset_posterior_state(s,k)
  for(nm in c("EL","EL2","EF","EF2"))
    assert_true(identical(x[[nm]],s[[nm]][,k,drop=FALSE]),msg=paste(nm,"kept in requested order"))
  for(nm in c("EBeta","EBeta2","prior_F"))
    assert_true(identical(x[[nm]],s[[nm]][k]),msg=paste(nm,"subset without reset"))
  assert_true(identical(x$Tau,s$Tau) && identical(x$g_beta,s$g_beta),
    msg="feature precisions and shared prior retained exactly")
  s2 <- f;s2$posterior_state <- NULL
  ref <- do.call(fit_cox_on_yf,c(list(.ps_y,1:40,rep(1,40)),.ps_a))
  assert_true(identical(s2,ref),msg="requesting the state changes no unpruned fitted quantity")
})
run_test("YFPosterior-T2: first resumed update uses retained beta and loading variance", {
  f <- do.call(fit_cox_on_yf,c(list(.ps_y,1:40,rep(1,40)),.ps_a,
    list(prune_preserve_posterior=TRUE)))
  s <- f$posterior_state
  # Deliberately nonzero beta and posterior variance distinguish a real resume
  # from the old mean-only initialization (beta=0, second moments=mean squared).
  s$EBeta <- c(.2,-.3,.4);s$EBeta2 <- s$EBeta^2+c(.01,.02,.03)
  s$EF2 <- s$EF^2+.02
  orig_t <- calc_cox_taylor_yf;orig_b <- update_beta_k
  eta0 <- b0 <- NULL
  tw <- function(...) {a<-list(...);if(is.null(eta0))eta0<<-a[[1]];do.call(orig_t,a)}
  bw <- function(...) {a<-list(...);if(is.null(b0))b0<<-a;do.call(orig_b,a)}
  assign("calc_cox_taylor_yf",tw,envir=.GlobalEnv)
  assign("update_beta_k",bw,envir=.GlobalEnv)
  tryCatch({
    a <- modifyList(.ps_a,list(max_iter=1L,posterior_init=s,prune_preserve_posterior=TRUE))
    do.call(fit_cox_on_yf,c(list(.ps_y,1:40,rep(1,40)),a))
  },finally={assign("calc_cox_taylor_yf",orig_t,envir=.GlobalEnv);assign("update_beta_k",orig_b,envir=.GlobalEnv)})
  assert_near(eta0,as.vector(.ps_y%*%s$EF%*%s$EBeta),tol=1e-12,
    msg="working predictor starts from retained beta, not zero")
  z <- as.vector(.ps_y%*%s$EF[,1])
  v <- as.vector((.ps_y^2)%*%(s$EF2[,1]-s$EF[,1]^2))
  assert_near(b0[[3]],z,tol=1e-12,msg="retained projection mean")
  assert_near(b0[[4]],z^2+v,tol=1e-10,msg="retained projection uncertainty in beta precision")
})
run_test("YFPosterior-T3: disabled option is identical and malformed states fail loudly", {
  a <- c(list(.ps_y,1:40,rep(1,40)),.ps_a,list(prune=TRUE))
  x <- do.call(fit_cox_on_yf,a)
  y <- do.call(fit_cox_on_yf,c(a,list(prune_preserve_posterior=FALSE)))
  assert_true(identical(x,y),msg="explicit FALSE equals historical mean-only pruning")
  e <- tryCatch(do.call(fit_cox_on_yf,c(a,list(posterior_init=list(EL=matrix(1,1,1))))),error=function(e)e)
  assert_true(inherits(e,"error"),msg="incomplete state rejected before fitting")
})

run_test("YFPosterior-T4: recursive pruning passes the exact retained state", {
  orig_f <- fit_cox_on_yf;orig_n <- prune_yf_nullcheck
  base <- do.call(orig_f,c(list(.ps_y,1:40,rep(1,40)),.ps_a,
    list(prune_preserve_posterior=TRUE)))
  st <- base$posterior_state
  # An exact point-mass dead column guarantees a real deletion, independent
  # of whether the small low-rank fixture's other programs merit pruning.
  for(nm in c("EL","EL2","EF","EF2"))st[[nm]] <- cbind(st[[nm]],0)
  st$EBeta <- c(st$EBeta,0);st$EBeta2 <- c(st$EBeta2,0)
  st$prior_F <- c(st$prior_F,list(list(family="point_laplace",pi=0,
    rate=NA_real_,variance=0,point_mass=TRUE)))
  st["g_beta"] <- list(NULL)
  expected <- st;resumed <- 0L;calls <- 0L
  nw <- function(fit,...) {
    d <- orig_n(fit,...)
    if(max(d)>=0)expected <<- yf_subset_posterior_state(fit$posterior_state,
      setdiff(seq_len(ncol(fit$EF)),which.max(d)))
    d
  }
  fw <- function(...) {
    a <- list(...);calls <<- calls+1L
    if(!is.null(a$posterior_init)) {
      assert_true(identical(a$posterior_init,expected),msg="no retained moment or prior reset at refit call")
      if(calls>1L)resumed <<- resumed+1L
    }
    do.call(orig_f,a)
  }
  assign("fit_cox_on_yf",fw,envir=.GlobalEnv);assign("prune_yf_nullcheck",nw,envir=.GlobalEnv)
  tryCatch({
    a <- c(list(.ps_y,1:40,rep(1,40)),modifyList(.ps_a,list(K=4L,
      posterior_init=st,beta_prior_mode="per_program",prior_beta="normal")),
      list(prune=TRUE,prune_protect_survival=TRUE,prune_preserve_posterior=TRUE))
    do.call(fit_cox_on_yf,a)
  },finally={assign("fit_cox_on_yf",orig_f,envir=.GlobalEnv);assign("prune_yf_nullcheck",orig_n,envir=.GlobalEnv)})
  assert_true(resumed>0L,msg="fixture actually exercises posterior-state pruning refits")
})
