# ============================================================
# Script: test_fixed_prior_panel.R
# Purpose: Verify fixed-prior evaluation design, training isolation and failure accounting.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: existing test helpers, yaml
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("results/benchmark_sim/run_round2_fixed_prior.R"))
  source("results/benchmark_sim/run_round2_fixed_prior.R")
run_test("FixedPrior-T1: the design contains80 simulation and20 CV fits without duplicates", {
  d <- fixed_prior_design()
  assert_true(nrow(d)==100L && sum(d$kind=="simulation")==80L &&
    sum(d$kind=="training_cv")==20L && !anyDuplicated(d),msg="fixed complete design")
  s <- d[d$kind=="simulation",]
  assert_true(all(table(s$model,s$noise)==20L) && identical(sort(unique(s$seed)),20261015L+1:20),
    msg="twenty original seeds for each model/noise")
})
run_test("FixedPrior-T2: event-stratified folds are deterministic and cover every patient", {
  e <- rep(c(0L,1L),72)
  a <- fixed_prior_folds(e);b <- fixed_prior_folds(e)
  assert_true(identical(a,b) && length(a)==144L && setequal(a,1:5),msg="one fixed fold per patient")
  t <- table(e,a)
  assert_true(all(apply(t,1,function(x)max(x)-min(x))<=1),msg="events balanced within one patient")
})
run_test("FixedPrior-T3: validation values cannot change training-fold scaling", {
  a <- matrix(1:24,8,3);v <- matrix(c(100,200,300,400,500,600),2,3)
  x <- fixed_prior_scale(a,v);y <- fixed_prior_scale(a,v*100)
  expected <- sweep(sweep(v,2,colMeans(a)),2,apply(a,2,sd),"/")
  assert_true(identical(x$training,y$training),msg="held-out values never influence training statistics")
  assert_near(x$validation,expected,tol=0,msg="apply frozen training-fold center and SD")
})
run_test("FixedPrior-T4: R1 changes only fixed shared prior, leaving other controls intact", {
  for(m in c("single","multimodal")) {
    a <- fixed_prior_controls(m,"fixed");b <- collapse_arm_controls(m,"selected")
    a$beta_prior_fix <- NULL
    if(m=="single")b$trace_supervision <- FALSE else b$trace <- FALSE
    assert_true(identical(a,b),msg=paste(m,"same selected options except prior fixation"))
  }
})
run_test("FixedPrior-T5: errors remain in recovery and convergence denominators", {
  d <- fixed_prior_design()
  d$error <- NA_character_;d$recovered <- TRUE;d$false_active <- 0
  d$converged <- TRUE;d$C <- .7
  ix <- which(d$kind=="simulation" & d$model=="single" & d$noise==.25)[1]
  d$error[ix] <- "planned";d$recovered[ix] <- NA;d$false_active[ix] <- NA;d$converged[ix] <- NA
  ref <- expand.grid(model=c("single","multimodal"),noise=c(.25,1),seed=1:20)
  ref$false_active <- 0
  s <- fixed_prior_summary(d,ref);x <- s$simulation[s$simulation$model=="single" & s$simulation$noise==.25,]
  assert_true(x$n==20L && x$errors==1L && x$recovered==19L &&
    !x$recovery_pass && !x$false_positive_pass && !x$convergence_pass,
    msg="a failed fit cannot pass by conditional denominators")
})
run_test("FixedPrior-T6: paired SE uses the same folds and rejects an undefined C", {
  d <- fixed_prior_design();d$error <- NA_character_;d$recovered <- TRUE
  d$false_active <- 0;d$converged <- TRUE;d$C <- .7
  ix <- which(d$kind=="training_cv" & d$model=="single" & d$arm=="fixed")
  delta <- c(-.01,-.02,0,.01,-.005);d$C[ix] <- .7+delta
  ref <- expand.grid(model=c("single","multimodal"),noise=c(.25,1),seed=1:20);ref$false_active <- 0
  s <- fixed_prior_summary(d,ref)$training_cv;x <- s[s$model=="single",]
  assert_near(x$delta,mean(delta),tol=1e-15,msg="fold-paired mean difference")
  assert_near(x$paired_SE,sd(delta)/sqrt(5),tol=1e-15,msg="one paired SE, not pooled errors")
  d$C[ix[1]] <- NA_real_
  x <- fixed_prior_summary(d,ref)$training_cv
  assert_true(!x$pass[x$model=="single"],msg="undefined held-out C cannot be silently removed")
})
run_test("FixedPrior-T7: each failed fit checkpoints the exact message and refuses overwrite", {
  out <- tempfile("fixed-prior-failure-");dir.create(file.path(out,"cells"),recursive=TRUE)
  had <- exists("make_testbed",.GlobalEnv,inherits=FALSE)
  old <- if(had)get("make_testbed",.GlobalEnv) else NULL
  assign("make_testbed",function(...)stop("planned fixed-prior software failure"),.GlobalEnv)
  a <- tryCatch(fixed_prior_cell(1L,out),finally={
    if(had)assign("make_testbed",old,.GlobalEnv) else rm("make_testbed",envir=.GlobalEnv)
  })
  saved <- read.csv(file.path(out,"cells","cell_001.csv"))
  assert_true(identical(a$error,saved$error) && grepl("planned fixed-prior",a$error),msg="exact failed-fit checkpoint")
  e <- tryCatch(fixed_prior_cell(1L,out),error=function(e)e)
  assert_true(inherits(e,"error"),msg="failed checkpoints also refuse overwrite")
  unlink(out,recursive=TRUE)
})
run_test("FixedPrior-T8: an unknown orchestration mode fails before loading real data", {
  e <- tryCatch(run_fixed_prior_panel(mode="invalid"),error=function(e)e)
  assert_true(inherits(e,"error") && grepl("mode must be",conditionMessage(e)),msg="no accidental fit on an invalid mode")
})
run_test("FixedPrior-T9: CSV type inference does not reject an unchanged checkpoint", {
  d <- fixed_prior_design()[1,,drop=FALSE]
  p <- tempfile();write.csv(d,p,row.names=FALSE);r <- read.csv(p);unlink(p)
  assert_true(fixed_prior_checkpoint_matches(r,d),msg="all-NA integer columns may reload as logical")
  r$seed <- r$seed+1L
  assert_true(!fixed_prior_checkpoint_matches(r,d),msg="changed seed must still be rejected")
})
