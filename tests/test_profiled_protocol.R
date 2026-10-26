# ============================================================
# Script: test_profiled_protocol.R
# Purpose: Verify profiled-curvature design, reference reuse and complete training checks.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing test helpers and fixed-prior runner
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("results/benchmark_sim/run_round2_fixed_prior.R")
run_test("ProfileProtocol-T1: Profiled curvature adds two full-training checks and reuses only references", {
  old <- fixed_prior_design();d <- fixed_prior_design("profiled")
  assert_true(nrow(old)==100L && nrow(d)==102L &&
    sum(d$kind=="simulation")==80L && sum(d$kind=="training_cv")==20L &&
    sum(d$kind=="training_full")==2L,msg="all specified simulation,paired CV and full-training checks")
  id <- which(d$kind=="training_cv" & d$arm=="framework")
  assert_true(length(id)==10L && identical(d[id,],old[id,]),msg="only unchanged framework folds are reusable")
  assert_true(sum(!(d$kind=="training_cv" & d$arm=="framework"))==92L,
    msg="92 new independent fits plus10 frozen CV references")
})
run_test("ProfileProtocol-T2: Candidate changes only the framework working curvature", {
  for(m in c("single","multimodal")) {
    a <- fixed_prior_controls(m,"profiled");a$profile_cox_intercept <- NULL
    b <- collapse_arm_controls(m,"framework")
    if(m=="single")b$trace_supervision <- FALSE else b$trace <- FALSE
    assert_true(identical(a,b),msg=paste(m,"no simultaneous prior fixation,warmup,weight or objective change"))
  }
})
run_test("ProfileProtocol-T3: cached references require identical folds,cache and all source cells", {
  out <- tempfile("r2-floor-reuse-");src <- tempfile("r2-fixed-source-")
  for(p in c(out,src)) {
    dir.create(file.path(p,"cells"),recursive=TRUE)
    writeLines("training_md5 synthetic-software-fixture",file.path(p,"training_source.txt"))
    write.csv(data.frame(patient_id=paste0("s",1:10),fold=rep(1:5,2)),file.path(p,"folds.csv"),row.names=FALSE)
  }
  d <- fixed_prior_design();d$error <- NA_character_;d$C <- .65
  d$C[d$kind=="training_cv" & d$model=="multimodal" & d$arm=="framework" & d$fold==1L] <- NA_real_
  write.csv(d,file.path(src,"results.csv"),row.names=FALSE)
  id <- which(d$kind=="training_cv" & d$arm=="framework")
  for(i in id)write.csv(d[i,],file.path(src,"cells",sprintf("cell_%03d.csv",i)),row.names=FALSE)
  r <- fixed_prior_reuse_cv(out,src)
  assert_true(nrow(r)==10L && length(list.files(file.path(out,"cells")))==10L,
    msg="exactly ten unchanged reference metrics reused")
  a <- read.csv(file.path(out,"cells",sprintf("cell_%03d.csv",id[1])))
  assert_true(a$C==.65 && !is.na(a$reference_md5),msg="unchanged C with source checksum")
  null_id <- which(d$kind=="training_cv" & d$model=="multimodal" & d$arm=="framework" & d$fold==1L)
  null <- read.csv(file.path(out,"cells",sprintf("cell_%03d.csv",null_id)))
  assert_true(is.na(null$C),msg="collapsed cached reference remains undefined")
  e <- tryCatch(fixed_prior_reuse_cv(out,src),error=function(e)e)
  assert_true(inherits(e,"error"),msg="reuse never overwrites checkpoints")
  writeLines("different training cache",file.path(out,"training_source.txt"))
  e <- tryCatch(fixed_prior_reuse_cv(out,src),error=function(e)e)
  assert_true(inherits(e,"error") && grepl("cache or folds differ",conditionMessage(e)),
    msg="reference cannot be reused for a changed input")
  unlink(c(out,src),recursive=TRUE)
})
run_test("ProfileProtocol-T4: nonconstant but noncredible training fits remain collapsed", {
  d <- fixed_prior_design("profiled");d$error <- NA_character_;d$recovered <- TRUE
  d$false_active <- 0;d$converged <- TRUE;d$C <- .7;d$risk_sd <- .1;d$n_active <- 1L
  d$n_active[d$kind=="training_full" & d$model=="multimodal"] <- 0L
  ref <- expand.grid(model=c("single","multimodal"),noise=c(.25,1),seed=1:20);ref$false_active <- 0
  s <- fixed_prior_summary(d,ref,"profiled")
  assert_true(all(s$simulation$recovery_pass) && all(s$training_cv$pass),msg="simulation and paired CV gates retained")
  f <- s$training_full
  assert_true(!f$no_collapse[f$model=="multimodal"] && f$no_collapse[f$model=="single"],
    msg="positive risk variance alone cannot substitute for a credible program")
})
run_test("ProfileProtocol-T5: full-training checks use every synthetic patient and frozen scaling", {
  source("code/predict_cox_on_yf.R");source("code/concordance_ci.R")
  set.seed(71)
  y <- matrix(rnorm(60*12),60,12,dimnames=list(paste0("s",1:60),paste0("g",1:12)))
  tr <- list(Y=list(expression=y,methylation=y/10),
    time=setNames(exp(rnorm(60)),rownames(y)),event=setNames(rep(c(1,0,1),20),rownames(y)))
  metadata <- list(joint=list(training_spec=list(subject_ids=rownames(y),feature_names=lapply(tr$Y,colnames))))
  F <- matrix(rnorm(12*7),12,7);beta <- c(1,rep(0,6))
  had_read <- exists("readRDS",.GlobalEnv,inherits=FALSE)
  old_read <- get("readRDS",.GlobalEnv,inherits=TRUE)
  had_fit <- exists("fit_cox_on_yf",.GlobalEnv,inherits=FALSE)
  old_fit <- if(had_fit)get("fit_cox_on_yf",.GlobalEnv) else NULL
  assign("readRDS",function(file,...) {
    if(identical(file,"data/round2_tcga_training.rds"))return(tr)
    if(grepl("benchmark_snapshot",file,fixed=TRUE))return(metadata)
    stop("Unexpected fixture input")
  },.GlobalEnv)
  assign("fit_cox_on_yf",function(Y,time,status,...) {
    args <- list(...)
    stopifnot(nrow(Y)==60L,identical(Y,scale(y)),isTRUE(args$profile_cox_intercept))
    list(EL=matrix(0,60,7),EF=F,EF_norms=rep(1,7),EBeta=beta,
      EBeta2=beta^2+.01,profile_cox_intercept=TRUE,
      history=list(converged=TRUE,n_iter=1L))
  },.GlobalEnv)
  out <- tempfile("full-training-routing-");dir.create(file.path(out,"cells"),recursive=TRUE)
  a <- tryCatch(fixed_prior_cell(101L,out,candidate="profiled"),finally={
    if(had_read)assign("readRDS",old_read,.GlobalEnv) else rm("readRDS",envir=.GlobalEnv)
    if(had_fit)assign("fit_cox_on_yf",old_fit,.GlobalEnv) else rm("fit_cox_on_yf",envir=.GlobalEnv)
  })
  expected <- frozen_reverse_cindex(as.vector(scale(y)%*%F%*%beta),tr$time,tr$event)
  assert_true(is.na(a$error) && a$kind=="training_full" && a$n_train==60L &&
    a$n_validation==60L && a$n_active==1L && isTRUE(a$profiled_curvature),msg="full training, not a partial fold or empty validation")
  assert_near(a$C,expected,tol=0,msg="training risk uses the fitted projection in the same units")
  unlink(out,recursive=TRUE)
})

run_test("ProfileProtocol-T6: a collapsed reference keeps paired C undefined without masking candidate results", {
  d <- fixed_prior_design("profiled");d$error<-NA_character_;d$recovered<-TRUE
  d$false_active<-0;d$converged<-TRUE;d$C<-.7;d$risk_sd<-.1;d$n_active<-1L
  d$C[d$kind=="training_cv" & d$model=="multimodal" & d$arm=="framework" & d$fold==1]<-NA_real_
  ref<-expand.grid(model=c("single","multimodal"),noise=c(.25,1),seed=1:20);ref$false_active<-0
  a<-fixed_prior_summary(d,ref,"profiled")$training_cv
  z<-a[a$model=="multimodal",]
  assert_true(!z$pass && is.na(z$delta) && is.na(z$paired_SE) &&
    z$candidate_finite_folds==5L && z$reference_finite_folds==4L &&
    z$candidate_mean_if_all5==.7 && is.na(z$reference_mean_if_all5),
    msg="no .5 substitution, deletion or concealed candidate mean")
})
