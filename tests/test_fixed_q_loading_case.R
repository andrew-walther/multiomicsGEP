# ============================================================
# Script: test_fixed_q_loading_case.R
# Purpose: Verify single-change controls, independent density checks and recovery reporting.
# Created: 2026-10-16
# Dependencies: existing test helpers and registered case driver
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/multimodal_yfb_updates.R")
if(file.exists("results/benchmark_sim/run_fixed_q_loading_case.R"))
  source("results/benchmark_sim/run_fixed_q_loading_case.R")
run_test("FixedQCase-T1: candidate has exactly one change from saved profiled controls", {
  a <- fixed_q_loading_case_controls()
  assert_equal(a$loading_prior_update,"fixed_q_em")
  a$loading_prior_update <- NULL
  assert_equal(a,fixed_prior_controls("single","profiled"))
  assert_true(length(formals(fixed_q_loading_case_controls))==0L)
})
run_test("FixedQCase-T2: quadrature agrees with independently integrated mixed density", {
  g <- list(pi=.4,rate=1.3,family="point_laplace",point_mass=FALSE)
  for(B in c(-5,0,2)) {
    q <- fixed_q_integrated_posterior(4,B,g)
    den <- function(f)g$pi*g$rate/2*exp(-g$rate*abs(f))*dnorm(B/4-f,0,.5)
    normalizer <- (1-g$pi)*dnorm(B/4,0,.5)+integrate(den,-Inf,0,rel.tol=1e-12)$value+integrate(den,0,Inf,rel.tol=1e-12)$value
    assert_near(q$log_ml,log(normalizer),tol=1e-10)
    for(r in 1:2) {
      fn <- function(f)f^r*den(f)/normalizer
      direct <- integrate(fn,-Inf,0,rel.tol=1e-12)$value+integrate(fn,0,Inf,rel.tol=1e-12)$value
      assert_near(q[[c("mean","second")[r]]],direct,tol=1e-9)
    }
  }
})
run_test("FixedQCase-T3: independent tail integration resolves the known approximation gap", {
  A <- 1/8.01^2;g <- list(pi=1,rate=1,family="point_laplace",point_mass=FALSE)
  q <- fixed_q_integrated_posterior(A,0,g)
  z <- integrate(function(f)exp(-A*f^2/2-f),0,Inf,rel.tol=1e-12)$value
  direct <- integrate(function(f)f*exp(-A*f^2/2-f),0,Inf,rel.tol=1e-12)$value/z
  assert_near(q$mean_abs,direct,tol=1e-10)
  assert_near(q$mean,0,tol=1e-14)
  approximate <- loading_prior_laplace_stats(A,0,g)
  assert_true(abs(q$mean_abs-approximate$mean_abs)>1e-5,
    msg="the diagnostic must detect the existing kernel approximation")
})
run_test("FixedQCase-T4: recovery preserves sign, credibility and false-program rules", {
  source("code/predict_cox_on_yf.R");source("code/concordance_ci.R")
  set.seed(19);truth <- cbind(c(1,0,0,0),c(0,1,0,0),c(1,-1,0,1))
  f <- list(EF=cbind(truth[,3],c(1,1,-1,-1)),EBeta=c(1,.3),
    EBeta2=c(1,.3)^2+.01,EF_norms=c(1,1),EL=matrix(.1,6,2),
    history=list(n_iter=3L,converged=TRUE))
  tb <- list(F_true=truth,Yva=matrix(rnorm(20),5,4),tva=1:5,sva=rep(1,5))
  a <- fixed_q_case_metrics(f,tb)
  assert_true(a$recovered && a$false_active==1L && a$n_active==2L && a$converged)
  assert_near(a$cor_planted,1,tol=1e-14)
  f$EBeta[1] <- -1
  assert_false(fixed_q_case_metrics(f,tb)$recovered,msg="opposite beta sign fails recovery")
})
