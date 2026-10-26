# ============================================================
# Script: test_multimodal_program_agreement.R
# Purpose: Verify raw-projection agreement and program-role accounting across modalities.
# Created: 2026-10-15
# Dependencies: base R, existing test helpers and feature correlation helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
if(file.exists("code/multimodal_program_agreement.R"))source("code/multimodal_program_agreement.R")
# One molecular block, two live factors and one exactly dead factor.
Y<-matrix(c(2,1,5,3,4,4,1,2,5,3,1,5,3,2,4),5,3,
  dimnames=list(paste0("p",1:5),c("cg1","cg2","cg3")))
L<-cbind(1:5,c(3,1,2,4,5),0)
F<-cbind(c(.2,.5,.8),c(.7,.3,.1),0)
old<-structure(list(EL=L,EF=list(methylation=F),EBeta=c(.3,-.2,0),
 EBeta2=c(.3^2+.0001,.2^2+.0001,0),training_spec=list(subject_ids=rownames(Y),
 modalities="methylation",feature_names=list(methylation=colnames(Y)))),class="multimodal_yfb_fit")
new<-old;perm<-c(2,1,3);scale<-c(2,3,1)
new$EL<-sweep(L[,perm],2,scale,"/")
new$EF$methylation<-rbind(sweep(F[,perm],2,scale,"*"),0)
new$EBeta<-old$EBeta[perm]/scale;new$EBeta2<-old$EBeta2[perm]/scale^2
new$training_spec$feature_names$methylation<-c(colnames(Y),"cg4")
newY<-cbind(Y,cg4=1:5)
run_test("MMProgram-T1: permutation/rescaling preserves raw risk and matched projection units", {
  z<-multimodal_feature_agreement(old,new,list(methylation=Y),list(methylation=newY))
  a<-z$pairs[z$pairs$old_factor%in%1:2,]
  assert_equal(a$new_factor,c(2L,1L));assert_near(a$activity_correlation,c(1,1),tol=1e-14)
  assert_near(a$loading_correlation,c(1,1),tol=1e-14)
  assert_near(a$old_beta_per_projection_sd,a$new_beta_per_projection_sd,tol=1e-14)
  assert_near(z$risk_correlation,1,tol=1e-14)
})
run_test("MMProgram-T2: dead factors remain undefined and unmatched", {
  z<-multimodal_feature_agreement(old,new,list(methylation=Y),list(methylation=newY))
  assert_true(all(is.na(z$activity_correlation[3,])))
  a<-z$pairs[z$pairs$old_factor%in%3L,]
  assert_true(is.na(a$new_factor));assert_true(is.na(a$loading_correlation))
})
run_test("MMProgram-T3: block norm shares are scale-invariant and zero-SD beta is not credible", {
  a<-multimodal_feature_roles(old,list(methylation=Y),2)
  b<-multimodal_feature_roles(new,list(methylation=newY),2)
  assert_near(sum(a$norm_share),1,tol=1e-14)
  assert_near(a$norm_share,b$norm_share[c(2,1,3)],tol=1e-14)
  assert_equal(a$credible,c(TRUE,TRUE,FALSE))
  assert_equal(a$nonzero_loading_means,c(3L,3L,0L))
})
run_test("MMProgram-T4: training patient-order mismatches cannot yield program agreement", {
  e<-tryCatch({multimodal_feature_agreement(old,new,list(methylation=Y[5:1,]),
    list(methylation=newY));""},error=conditionMessage)
  assert_equal(e,"Program comparison requires the exact ordered training patients and features.")
})
run_test("MMProgram-T5: two-block adapter preserves the raw joint projection and modality units", {
  joint<-old;joint$EF<-list(expression=F/4,methylation=F)
  joint$training_spec$modalities<-c("expression","methylation")
  joint$training_spec$feature_names<-list(expression=paste0("g",1:3),methylation=colnames(Y))
  expression<-Y*2;colnames(expression)<-paste0("g",1:3)
  blocks<-list(expression=expression,methylation=Y)
  a<-multimodal_feature_adapter(joint,blocks)
  assert_near(a$projection,1.5*(Y%*%F),tol=1e-14)
  assert_equal(length(unique(a$genes)),6L)
  z<-multimodal_feature_agreement(joint,joint,blocks,blocks)
  assert_equal(names(z$loading_correlation),c("expression","methylation"))
  assert_near(z$risk_correlation,1,tol=1e-14)
  roles<-multimodal_feature_roles(joint,blocks,2)
  for(m in names(blocks))assert_near(sum(roles$norm_share[roles$modality==m]),1,tol=1e-14)
})
