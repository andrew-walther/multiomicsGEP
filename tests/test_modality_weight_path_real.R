# ============================================================
# Script: test_modality_weight_path_real.R
# Purpose: Check descriptive weight-path precision and patient-score matching.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: existing test helpers
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
source("results/multimodal_real/run_modality_weight_path.R")
run_test("TCGA weight path: weighted second moments determine modality classification", {
  fit <- list(EF = matrix(1, 4, 3), EF2 = cbind(c(1,1,0,0),c(0,0,1,1),rep(1,4)),
    Tau = rep(1,4), modality_weight = c(expression=1,methylation=.3))
  x <- mwpath_programs(fit, rep(c("expression","methylation"),each=2),
    cbind(1:6,(1:6)^2,sin(1:6)), cbind(1:6), cbind((1:6)^2))
  assert_equal(x$class, c("expression-exclusive","methylation-exclusive","joint"))
  assert_near(x$methylation_share,c(0,1,.3/1.3),1e-12)
  assert_true(all(x$expression_empty == c(FALSE,TRUE,FALSE)))
})
run_test("TCGA weight path: new programs require low patient-score correlations to both references", {
  Q <- qr.Q(qr(cbind(rep(1,6),1:6,(1:6)^2,sin(1:6))))[,-1,drop=FALSE]
  fit <- list(EF=matrix(1,4,3),EF2=matrix(1,4,3),Tau=rep(1,4))
  x <- mwpath_programs(fit,rep(c("expression","methylation"),each=2),Q,Q[,1,drop=FALSE],Q[,2,drop=FALSE])
  assert_equal(x$new_program,c(FALSE,FALSE,TRUE))
  assert_near(x$expression_abs_r,c(1,0,0),1e-12)
  assert_near(x$methylation_abs_r,c(0,1,0),1e-12)
})
run_test("TCGA weight path: exclusive CpG top lists preserve island and promoter annotations", {
  fit <- list(EF=matrix(c(0,0,-3,2),4,1))
  p <- data.frame(program=1L,live=TRUE,class="methylation-exclusive",new_program=FALSE)
  ann <- data.frame(UCSC_RefGene_Name=c("G1","G2"),UCSC_RefGene_Group=c("Body;TSS200","Body"),
    Relation_to_Island=c("Island","OpenSea"),row.names=c("cg1","cg2"))
  x <- mwpath_top(fit,p,c("A","B","cg1","cg2"),rep(c("expression","methylation"),each=2),ann)
  assert_equal(x$feature,c("cg1","cg2"));assert_equal(x$promoter,c(TRUE,FALSE))
  assert_equal(x$island_relation,c("Island","OpenSea"));assert_near(x$loading,c(-3,2),0)
})
if (sys.nframe()==0L) {
  result <- report_results("TCGA weight path")
  if (result$failed>0L) quit(status=1L)
}
