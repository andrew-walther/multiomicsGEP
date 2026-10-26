# ============================================================
# Script: test_comparator_small_rank.R
# Purpose: Verify fixed EBMF zero/one-factor handling and arm error isolation.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, flashier, ebnm, glmnet, existing test helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("results/benchmark_sim/run_supervision_comparators.R")
source("code/concordance_ci.R")
set.seed(84)
.rank_F <- rnorm(20)
.rank_Y <- matrix(rnorm(100*20),100,20)
.rank_V <- matrix(rnorm(50*20),50,20)
.rank_tb <- list(Ytr=.rank_Y,Yva=.rank_V,
  ttr=exp(-.2*as.vector(.rank_Y%*%.rank_F)+rnorm(100,sd=.5)),str=rep(1,100),
  tva=exp(-.2*as.vector(.rank_V%*%.rank_F)+rnorm(50,sd=.5)),sva=rep(1,50),
  F_true=cbind(rnorm(20),rnorm(20),.rank_F),beta_true=c(0,0,3.3))
run_test("ComparatorRank-T1: zero-factor datasets retain all arms and the denominator", {
  a <- comparator_ebmf_rows(matrix(numeric(),20,0),.rank_tb,1L)
  assert_true(nrow(a)==3L && all(a$K_final==0L) && all(!a$recovered) && all(is.na(a$C)),
    msg="constant risk has undefined C and counts as not recovered in every method")
})
run_test("ComparatorRank-T2: one factor uses ordinary Cox and records lasso inapplicability", {
  a <- comparator_ebmf_rows(matrix(.rank_F,20,1),.rank_tb,1L)
  c <- a[a$method=="ebmf_cox",];l <- a[a$method=="ebmf_lasso_cox",]
  assert_true(is.na(c$error) && c$recovered && is.finite(c$C),msg="scalar Cox is fitted and can recover the true program")
  assert_true(!l$applicable && !l$recovered && is.na(l$error) && is.na(l$C) &&
    grepl("one factor",l$reason),msg="lasso is explicitly not applicable, not silently padded")
})
run_test("ComparatorRank-T3: a Cox error does not erase molecular recovery", {
  tb <- .rank_tb;tb$ttr[] <- NA_real_
  a <- comparator_ebmf_rows(cbind(.rank_F,rnorm(20)),tb,1L)
  u <- a[a$method=="ebmf_alone",];c <- a[a$method=="ebmf_cox",]
  assert_true(u$recovered && is.na(u$error) && !is.na(c$error),
    msg="an invalid survival arm leaves successful EBMF-alone metrics intact")
})
run_test("ComparatorRank-T4: summaries keep failed and inapplicable datasets in n", {
  r <- data.frame(method="test",noise=1,recovered=c(TRUE,rep(FALSE,19)),
    C=c(.7,rep(NA_real_,19)),cor_planted=c(.9,rep(NA_real_,19)),
    false_active=NA_real_,K_final=c(1,rep(NA_real_,19)),converged=NA,
    applicable=c(TRUE,rep(FALSE,19)),error=c(rep(NA_character_,19),"planned"))
  s <- comparator_panel_summary(r)
  assert_true(s$n==20L && s$recovered==1L && s$errors==1L && s$n_C_defined==1L,
    msg="no conditional denominator after failures or structural C undefined")
  assert_near(s$mean_C_defined,.7,tol=0,msg="defined-only C is labelled with its denominator")
})
run_test("ComparatorRank-T5: retained-rank counts distinguish errors from zero factors", {
  r <- data.frame(method="ebmf_alone",noise=1,K_final=c(0,1,2,NA),
    error=c(NA,NA,NA,"planned"))
  x <- comparator_ebmf_rank_counts(r)
  assert_true(x$n==4L && x$n_zero==1L && x$n_one==1L && x$n_two_or_more==1L && x$n_error==1L,
    msg="an EBMF failure is not misclassified as a fitted zero-factor model")
})
run_test("ComparatorRank-T6: a supervised-arm error retains every other method", {
  load_comparator_code()
  oldfit <- fit_cox_on_yf;oldbed <- make_testbed
  assign("cfg",yaml::read_yaml("config/globals.yml"),.GlobalEnv)
  assign("fit_cox_on_yf",function(...)stop("planned supervised-arm failure"),.GlobalEnv)
  assign("make_testbed",function(seed,noise).rank_tb,.GlobalEnv)
  a <- tryCatch(run_comparator_dataset(1L,1),finally={
    assign("fit_cox_on_yf",oldfit,.GlobalEnv);assign("make_testbed",oldbed,.GlobalEnv)
  })
  assert_true(nrow(a)==5L && setequal(a$method,c("ebmf_alone","ebmf_cox","ebmf_lasso_cox","pca_cox","supervised")),
    msg="all method rows survive one failed fitter")
  assert_true(grepl("planned supervised-arm failure",a$error[a$method=="supervised"]) &&
    is.na(a$error[a$method=="ebmf_alone"]),msg="the exact error stays with its own arm")
})
