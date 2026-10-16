# ============================================================
# Script: test_multimodal_feature_metrics.R
# Purpose: Check feature-study risk scoring without assigning concordance to collapsed fits.
# Created: 2026-10-15
# Dependencies: existing test helpers, survival
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("code/multimodal_feature_sensitivity.R"))source("code/multimodal_feature_sensitivity.R")
run_test("MMFeature-T1: perfect risk ordering has C1 and shifted risk preserves partial likelihood", {
  z<-multimodal_feature_metrics(1:5,rep(1,5),5:1)
  shifted<-multimodal_feature_metrics(1:5,rep(1,5),5:1+100)
  assert_near(z$C,1,tol=1e-14);assert_near(z$logPL,shifted$logPL,tol=1e-12)
})
run_test("MMFeature-T2: constant risk stays C=NA with the actual null partial likelihood", {
  z<-multimodal_feature_metrics(1:5,rep(1,5),rep(0,5))
  assert_true(is.na(z$C));assert_near(z$logPL,-log(factorial(5)),tol=1e-12)
  assert_near(z$risk_sd,0,tol=0)
})
run_test("MMFeature-T3: nonfinite scores fail before a plausible-looking metric is returned", {
  e<-tryCatch({multimodal_feature_metrics(1:5,rep(1,5),c(1:4,NA));""},error=conditionMessage)
  assert_equal(e,"Finite aligned survival vectors of length at least two required.")
})
