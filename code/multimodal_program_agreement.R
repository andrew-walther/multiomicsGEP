# ============================================================
# Script: multimodal_program_agreement.R
# Purpose: Compare saved programs on shared features and report molecular block roles without refitting.
# Created: 2026-10-15
# Dependencies: base R and existing single-block program-agreement helpers
# ============================================================
source("code/feature_sensitivity.R")

#' Align a saved multimodal fit and construct its raw-projection adapter
#' Stacking blocks gives YF=sum_m Y_m F_m. EF_norms=1 therefore preserves
#' the raw multimodal predictor exactly when using the existing matcher.
#' @param fit Saved multimodal fit with training_spec, EL, EF and EBeta.
#' @param Y Named raw training matrices on the fit's exact ordered patients.
#' @return Aligned blocks, prefixed feature identities, stacked data/fit and raw projection.
#' @examples
#' # multimodal_feature_adapter(fit,training$Y)
#' @family multimodal_program_agreement
#' @seealso feature_program_agreement
multimodal_feature_adapter <- function(fit,Y) {
  spec <- fit$training_spec
  if(!inherits(fit,"multimodal_yfb_fit") || !is.list(Y) ||
     !all(spec$modalities%in%names(Y)))
    stop("Program comparison requires the exact ordered training patients and features.")
  blocks <- Y[spec$modalities];features <- spec$feature_names
  for(m in spec$modalities) {
    if(!is.matrix(blocks[[m]]) || any(!is.finite(blocks[[m]])) ||
       !identical(rownames(blocks[[m]]),spec$subject_ids) ||
       anyDuplicated(colnames(blocks[[m]])) || !all(features[[m]]%in%colnames(blocks[[m]])))
      stop("Program comparison requires the exact ordered training patients and features.")
    blocks[[m]] <- blocks[[m]][,features[[m]],drop=FALSE]
    if(nrow(fit$EF[[m]])!=length(features[[m]]) || ncol(fit$EF[[m]])!=ncol(fit$EL))
      stop("Loading dimensions do not match the saved training specification.")
  }
  if(nrow(fit$EL)!=length(spec$subject_ids) || length(fit$EBeta)!=ncol(fit$EL))
    stop("Score/coefficient dimensions do not match the saved training specification.")
  genes <- unlist(lapply(spec$modalities,function(m)paste(m,features[[m]],sep="::")),use.names=FALSE)
  stacked <- do.call(cbind,blocks)
  adapted <- list(EL=fit$EL,EF=do.call(rbind,fit$EF[spec$modalities]),
    EBeta=fit$EBeta,EF_norms=rep(1,ncol(fit$EL)))
  list(blocks=blocks,genes=genes,Y=stacked,fit=adapted,
    projection=stacked%*%adapted$EF)
}

#' Match saved programs and compare loadings within each shared modality
#' Activity matching reuses the existing non-reused greedy absolute-correlation
#' convention. Constant factors stay unmatched. Combined cross-modality loading
#' correlations are discarded because their units differ; report each block.
#' @param old_fit,new_fit Saved multimodal fits on identical training patients.
#' @param old_Y,new_Y Raw training matrices on the corresponding feature sets.
#' @return Long matched/unmatched pairs, activity/projection correlation matrices,
#'   loading matrices by modality and raw-risk correlation. No scoring outcomes used.
#' @examples
#' # multimodal_feature_agreement(fit_small,fit_large,Y_small,Y_large)
#' @family multimodal_program_agreement
#' @seealso multimodal_feature_adapter, feature_program_agreement
multimodal_feature_agreement <- function(old_fit,new_fit,old_Y,new_Y) {
  a <- multimodal_feature_adapter(old_fit,old_Y)
  b <- multimodal_feature_adapter(new_fit,new_Y)
  if(!identical(old_fit$training_spec$subject_ids,new_fit$training_spec$subject_ids))
    stop("Program comparison requires the exact ordered training patients and features.")
  z <- feature_program_agreement(a$fit,a$genes,b$fit,b$genes,a$Y,b$Y)
  pairs <- matrices <- list()
  for(m in union(names(old_fit$EF),names(new_fit$EF))) {
    common <- intersect(old_fit$training_spec$feature_names[[m]],new_fit$training_spec$feature_names[[m]])
    r <- matrix(NA_real_,ncol(old_fit$EL),ncol(new_fit$EL))
    if(length(common)>=2L)r <- feature_cor(
      old_fit$EF[[m]][match(common,old_fit$training_spec$feature_names[[m]]),,drop=FALSE],
      new_fit$EF[[m]][match(common,new_fit$training_spec$feature_names[[m]]),,drop=FALSE])
    p <- z$pairs;p$loading_correlation <- NA_real_
    for(i in seq_len(nrow(p)))if(!is.na(p$old_factor[i]) && !is.na(p$new_factor[i]))
      p$loading_correlation[i] <- r[p$old_factor[i],p$new_factor[i]]
    p$modality <- m;p$n_common_features <- length(common)
    pairs[[m]] <- p;matrices[[m]] <- r
  }
  list(pairs=do.call(rbind,pairs),activity_correlation=z$activity_correlation,
    projection_correlation=z$projection_correlation,loading_correlation=matrices,
    risk_correlation=z$risk_correlation)
}

#' Describe each program's molecular norm share and prognostic coefficient
#' The share is ||L_k F_mk'||_F^2 divided by the sum over rank-one programs
#' within that block. It is neither additive PVE nor survival contribution.
#' @param fit Saved multimodal fit with coefficient second moments.
#' @param Y Raw matrices on the exact training patients/features.
#' @param z_threshold Existing configured prognostic-z criterion, supplied explicitly.
#' @return Long factor-by-modality table with shares, nonzero means, coefficient
#'   uncertainty, per-projection-SD coefficients and positive-SD credibility.
#' @examples
#' # multimodal_feature_roles(fit,training$Y,cfg$multimodal_yfb$prognostic_z_threshold)
#' @family multimodal_program_agreement
#' @seealso multimodal_feature_adapter
multimodal_feature_roles <- function(fit,Y,z_threshold) {
  a <- multimodal_feature_adapter(fit,Y)
  if(length(z_threshold)!=1L || !is.finite(z_threshold) || z_threshold<=0)
    stop("Supply the existing positive prognostic-z criterion.")
  sd_beta <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
  credible <- sd_beta>0 & abs(fit$EBeta)>=z_threshold*sd_beta
  units <- fit$EBeta*apply(a$projection,2,stats::sd)
  rows <- list()
  for(m in names(a$blocks)) {
    norm2 <- colSums(fit$EL^2)*colSums(fit$EF[[m]]^2)
    rows[[m]] <- data.frame(factor=seq_along(fit$EBeta),modality=m,
      norm_share=if(sum(norm2)>0)norm2/sum(norm2) else NA_real_,
      nonzero_loading_means=as.integer(colSums(fit$EF[[m]]!=0)),beta=fit$EBeta,
      sd_beta=sd_beta,credible=credible,beta_per_projection_sd=units)
  }
  do.call(rbind,rows)
}
