# ============================================================
# Script: test_feature_sensitivity.R
# Purpose: Verify nested feature manifests and program agreement with dead factors.
# Created: 2026-10-15
# Dependencies: base R, existing preprocessing/prediction and test helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
source("code/preprocess_desurv.R")
if(file.exists("code/feature_sensitivity.R"))source("code/feature_sensitivity.R")
run_test("Feature-T1: shared counts follow the original cohort intersection and reach all genes", {
  g<-paste0("g",1:6);x<-outer(1:5,1:6);colnames(x)<-g
  z<-feature_count_manifest(list(a=x,b=x),g,c(2L,4L,6L,7L),2L)
  assert_equal(z$design$actual_shared,c(2L,4L,6L,NA_integer_))
  assert_true(setequal(z$genes[[1]],c("g5","g6")))
  assert_true(all(z$genes[[1]]%in%z$genes[[2]]) && setequal(z$genes[[3]],g))
  assert_true(is.null(z$genes[[4]]),msg="unavailable target is explicit")
})
run_test("Feature-T2: entry ranks reproduce combined-rank screening including ties", {
  set.seed(31);g<-paste0("g",1:8)
  a<-matrix(runif(48),6,8,dimnames=list(NULL,g));b<-a[,8:1];colnames(b)<-g
  base<-intersect(select_top_variable_genes(a,g,4,"combined_rank")$gene_names,
    select_top_variable_genes(b,g,4,"combined_rank")$gene_names)
  z<-feature_count_manifest(list(a=a,b=b),g,c(length(base),8L),4L)
  for(q in 1:8){
    wanted<-intersect(select_top_variable_genes(a,g,q,"combined_rank")$gene_names,
      select_top_variable_genes(b,g,q,"combined_rank")$gene_names)
    assert_true(setequal(g[z$entry_rank<=q],wanted))
  }
})
run_test("Feature-T3: finite correlations match base R and constant columns remain undefined", {
  a<-cbind(c(1,3,2,5),rep(0,4));b<-cbind(c(2,1,4,3),c(4,2,5,1))
  z<-feature_cor(a,b)
  assert_near(z[1,],stats::cor(a[,1],b),tol=1e-14)
  assert_true(all(is.na(z[2,])))
})
run_test("Feature-T4: matched common-feature loadings survive reordering and scaling", {
  L<-cbind(c(1,2,4,3,5),c(3,1,2,4,5))
  F<-cbind(c(.2,.5,.8),c(.7,.3,.1));g<-c("a","b","c")
  Y<-matrix(c(2,1,5,3,4,4,1,2,5,3,1,5,3,2,4),5,3)
  old<-list(EL=L,EF=F,EBeta=c(.3,-.2),EF_norms=sqrt(colSums(F^2)))
  NF<-sweep(rbind(F,c(0,0))[,2:1],2,c(3,2),"*")
  new<-list(EL=sweep(L[,2:1],2,c(3,2),"/"),EF=NF,
    EBeta=old$EBeta[2:1],EF_norms=sqrt(colSums(NF^2)))
  z<-feature_program_agreement(old,g,new,c(g,"d"),Y,cbind(Y,1:5))
  paired<-z$pairs[!is.na(z$pairs$old_factor),]
  assert_equal(paired$new_factor,c(2L,1L))
  assert_near(paired$activity_correlation,c(1,1),tol=1e-14)
  assert_near(paired$loading_correlation,c(1,1),tol=1e-14)
  assert_near(z$risk_correlation,1,tol=1e-14)
})
run_test("Feature-T5: dead factors stay unmatched with NA correlations", {
  L<-cbind(c(1,2,4,3,5),c(3,1,2,4,5),0)
  F<-cbind(c(.2,.5,.8),c(.7,.3,.1),0);Y<-matrix(1:15,5,3)
  f<-list(EL=L,EF=F,EBeta=c(.3,-.2,0),EF_norms=c(sqrt(colSums(F[,1:2]^2)),1))
  z<-feature_program_agreement(f,c("a","b","c"),f,c("a","b","c"),Y,Y)
  assert_true(all(is.na(z$activity_correlation[3,])))
  assert_true(is.na(z$pairs$new_factor[z$pairs$old_factor%in%3L]))
  matched<-z$pairs$new_factor[!is.na(z$pairs$old_factor)&!is.na(z$pairs$new_factor)]
  assert_true(length(matched)==2L && !anyDuplicated(matched))
  assert_true(any(is.na(z$pairs$old_factor)&z$pairs$new_factor%in%3L))
})
run_test("Feature-T6: common loadings are aligned by identifier, not row position", {
  L<-cbind(1:5,c(1,4,2,5,3));F<-cbind(c(.1,.5,.8),c(.8,.4,.2))
  Y<-matrix(c(2,1,5,3,4,4,1,2,5,3,1,5,3,2,4),5,3)
  old<-list(EL=L,EF=F,EBeta=c(.3,-.2),EF_norms=sqrt(colSums(F^2)))
  new<-old;new$EF<-F[c(3,1,2),]
  z<-feature_program_agreement(old,c("a","b","c"),new,c("c","a","b"),Y,Y[,c(3,1,2)])
  assert_near(z$pairs$loading_correlation,c(1,1),tol=1e-14)
  assert_near(z$risk_correlation,1,tol=1e-14)
})
run_test("Feature-T7: coefficient comparisons use projection-SD units", {
  L<-cbind(1:5,c(1,4,2,5,3));F<-cbind(c(.1,.5,.8),c(.8,.4,.2))
  Y<-matrix(c(2,1,5,3,4,4,1,2,5,3,1,5,3,2,4),5,3)
  f<-list(EL=L,EF=F,EBeta=c(.3,-.2),EF_norms=sqrt(colSums(F^2)))
  z<-feature_program_agreement(f,c("a","b","c"),f,c("a","b","c"),Y,Y)
  expected<-f$EBeta*apply(Y%*%sweep(F,2,f$EF_norms,"/"),2,sd)
  assert_near(z$pairs$old_beta_per_projection_sd,expected,tol=1e-14)
})
run_test("Feature-T8: invalid patient dimensions and gene identities fail loudly", {
  L<-cbind(1:5,c(1,4,2,5,3));F<-cbind(c(.1,.5,.8),c(.8,.4,.2))
  Y<-matrix(1:15,5,3)
  f<-list(EL=L,EF=F,EBeta=c(.3,-.2),EF_norms=sqrt(colSums(F^2)))
  feature_program_agreement(f,c("a","b","c"),f,c("a","b","c"),Y,Y)
  assert_true(inherits(tryCatch(feature_program_agreement(f,c("a","a","c"),
    f,c("a","b","c"),Y,Y),error=identity),"error"))
  assert_true(inherits(tryCatch(feature_program_agreement(f,c("a","b","c"),
    f,c("a","b","c"),Y,Y[1:4,]),error=identity),"error"))
})
run_test("Feature-T9: named patient-order mismatches cannot yield agreement statistics", {
  L<-cbind(1:5,c(1,4,2,5,3));F<-cbind(c(.1,.5,.8),c(.8,.4,.2))
  Y<-matrix(1:15,5,3,dimnames=list(paste0("p",1:5),NULL))
  f<-list(EL=L,EF=F,EBeta=c(.3,-.2),EF_norms=sqrt(colSums(F^2)))
  err<-tryCatch(feature_program_agreement(f,c("a","b","c"),
    f,c("a","b","c"),Y,Y[5:1,]),error=identity)
  assert_true(inherits(err,"error"),msg="same number of rows does not establish same patients/order")
})
if(sys.nframe()==0L){report_results("feature sensitivity");if(.test_env$failed>0L)quit(status=1L)}
run_test("Feature-T10: survival metrics preserve risk direction and actual partial likelihood", {
  z<-feature_survival_metrics(1:5,rep(1,5),5:1)
  assert_near(z$C,1,tol=1e-14)
  assert_near(z$logPL,calc_cox_taylor_yf(5:1,1:5,rep(1,5))$logPL,tol=1e-12)
  assert_near(z$risk_sd,stats::sd(5:1),tol=1e-14)
})
run_test("Feature-T11: constant risk stays undefined for C, with finite null partial likelihood", {
  z<-feature_survival_metrics(1:5,rep(1,5),rep(0,5))
  assert_true(is.na(z$C));assert_near(z$logPL,-log(factorial(5)),tol=1e-12)
})
run_test("Feature-T12: survival metrics reject nonfinite predictions", {
  failed<-tryCatch({feature_survival_metrics(1:5,rep(1,5),c(1:4,NA));FALSE},error=function(e)identical(conditionMessage(e),"Finite aligned survival outcomes and risk required."))
  assert_true(failed)
})
