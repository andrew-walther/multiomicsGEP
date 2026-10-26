# ============================================================
# Script: test_cpg_feature_manifest.R
# Purpose: Verify unchanged variance-ranked CpG sets and exact frozen-input checks.
# Created: 2026-10-15
# Dependencies: base R and existing test helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/cpg_feature_sensitivity.R")
run_test("CpGFeature-T1: nested counts preserve variance rank and exclude constants", {
  Y<-cbind(cg1=c(0,.4,.8),cg2=c(.1,.2,.3),cg3=c(.7,.75,.8),cg4=rep(.4,3))
  z<-cpg_count_manifest(Y,c(1L,2L,3L,4L),"cg1")
  assert_equal(z$design$actual,c(1L,2L,3L,NA_integer_))
  assert_equal(z$features[[2]],c("cg1","cg2"));assert_equal(z$features[[3]],c("cg1","cg2","cg3"))
  assert_true(is.null(z$features[[4]]));assert_equal(z$design$eligible,rep(3L,4))
})
run_test("CpGFeature-T2: a different frozen screen fails rather than being silently replaced", {
  Y<-cbind(cg1=c(0,.4,.8),cg2=c(.1,.2,.3))
  e<-tryCatch({cpg_count_manifest(Y,c(1L,2L),"cg2");""},error=conditionMessage)
  assert_equal(e,"Variance ranking does not reproduce the frozen CpG screen.")
})
run_test("CpGFeature-T3: reconstruction is exact by patient order and probe identity", {
  Y<-matrix(seq_len(12)/20,3,4,dimnames=list(paste0("p",1:3),paste0("cg",1:4)))
  frozen<-Y[,c(3,1),drop=FALSE]
  assert_true(validate_cpg_reconstruction(Y[,4:1],frozen))
  bad<-Y;bad[1,3]<-bad[1,3]+1e-12
  e<-tryCatch({validate_cpg_reconstruction(bad,frozen);""},error=conditionMessage)
  assert_equal(e,"Reconstructed frozen CpG values differ; stop before fitting.")
})
run_test("CpGFeature-T4: patient-order and nonfinite mismatches cannot pass", {
  Y<-matrix(seq_len(12)/20,3,4,dimnames=list(paste0("p",1:3),paste0("cg",1:4)))
  frozen<-Y[,1:2]
  for(bad in list(Y[3:1,],replace(Y,1,NA_real_))) {
    e<-tryCatch({validate_cpg_reconstruction(bad,frozen);""},error=conditionMessage)
    assert_equal(e,"Finite uniquely named matrices on the identical patient order required.")
  }
})
