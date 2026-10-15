# ============================================================
# Script: test_desurv_h2h.R
# Purpose: Verify reported author settings and isolated seed checkpoints.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: DeSurv in the isolated Round2 library, existing test helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
.libPaths(c(normalizePath("data/round2_R_library"),.libPaths()))
if(file.exists("results/benchmark_sim/run_desurv_h2h.R"))
  source("results/benchmark_sim/run_desurv_h2h.R")
run_test("DeSurvH2H-T1: author settings retain full precision and the full seed protocol", {
  a <- desurv_reported_controls(2L)
  assert_near(c(a$alpha,a$lambda,a$nu),
    c(.333935556177224,.349448737941329,.0556512168045738),tol=0,
    msg="reported training parameters copied at full precision")
  assert_true(a$k==3L && a$ninit==1L && a$maxit==4000L && a$imaxit==4000L &&
    a$tol==1e-7 && a$tol_init==1e-5 && a$seed==2L && !a$parallel_init,
    msg="full paper seed protocol, without quick-mode changes")
})
run_test("DeSurvH2H-T2: seed checkpoint matches a direct author fit and refuses overwrite", {
  set.seed(2026)
  X <- matrix(rexp(12*30),12,30,dimnames=list(paste0("g",1:12),paste0("s",1:30)))
  d <- list(ex=X,sampInfo=data.frame(time=rexp(30),event=rep(c(1,1,0),10)))
  out <- tempfile("desurv-seed-");dir.create(out)
  row <- desurv_seed_fit(1L,d,out)
  direct <- do.call(DeSurv::desurv_fit,c(list(X=d$ex,y=d$sampInfo$time,d=d$sampInfo$event),desurv_reported_controls(1L)))
  got <- readRDS(file.path(out,"seed_001.rds"))
  for(n in c("W","H","beta","loss"))assert_true(identical(got[[n]],direct[[n]]),
    msg=paste("same author optimizer output",n))
  assert_true(is.na(row$error) && row$converged==direct$convergence,
    msg="actual convergence flag is preserved")
  e <- tryCatch(desurv_seed_fit(1L,d,out),error=function(e)e)
  assert_true(inherits(e,"error"),msg="completed seed never overwritten")
  unlink(out,recursive=TRUE)
})
run_test("DeSurvH2H-T3: a failed seed remains an explicit checkpoint", {
  out <- tempfile("desurv-bad-seed-");dir.create(out)
  row <- desurv_seed_fit(1L,list(ex=matrix(-1,4,4),sampInfo=data.frame(time=1:4,event=c(1,0,1,0))),out)
  saved <- read.csv(file.path(out,"seed_001_status.csv"))
  assert_true(!is.na(row$error) && identical(row$error,saved$error),
    msg="invalid nonnegative input retains the exact method error")
  unlink(out,recursive=TRUE)
})
run_test("DeSurvH2H-T4: full BO receives raw training values and filters inside its folds", {
  x <- matrix(c(4,9,2,6,8,3,1,7),2,4)
  d <- list(ex=x,sampInfo=data.frame(time=1:4,event=c(1,0,1,0),dataset=c("A","A","B","B")))
  a <- desurv_bo_controls(d,3L)
  assert_true(identical(a$X,x) && isTRUE(a$preprocess) && a$method_trans_train=="rank",
    msg="no full-training feature screen or rank transform before cross-validation")
  assert_true(a$nfolds==5L && a$bo_fixed$n_starts==30L &&
    a$coarse_control$n_init==50L && a$coarse_control$n_iter==100L,
    msg="paper full search, not its quick mode")
})
run_test("DeSurvH2H-T5: an invalid entry mode stops before loading data or fitting", {
  e <- tryCatch(run_desurv_h2h("invalid"),error=function(e)e)
  assert_true(inherits(e,"error") && grepl("mode must be",conditionMessage(e)),
    msg="entry guard fails loudly before any analysis side effect")
})
run_test("DeSurvH2H-T6: full BO postprocessing follows the paper's one-SE rank rule", {
  h <- data.frame(eval_id=1:3,k_grid=c(3L,5L,7L),mean_cindex=c(.68,.70,.69),
    c_se=rep(.03,3),status="ok",alpha_grid=c(.1,.2,.3))
  x <- list(history=h,overall_best=list(params=list(k_grid=5L,alpha_grid=.2)))
  a <- desurv_bo_summary(x)
  assert_true(a$k_selection$k_selected==3L && a$params$k==3L && a$params$alpha==.2,
    msg="paper03 overrides overall-best rank with its one-SE rank, preserving other reported parameters")
  x$history$mean_cindex[1] <- .60
  b <- desurv_bo_summary(x)
  assert_true(b$params$k==5L,msg="a materially inferior smaller rank is not selected")
})
