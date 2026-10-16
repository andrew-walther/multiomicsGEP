# ============================================================
# Script: feature_sensitivity.R
# Purpose: Build nested feature manifests and compare fitted programs without forcing dead matches.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: base R; existing DeSurv preprocessing and YFB prediction
# ============================================================
if(!exists("select_top_variable_genes"))source("code/preprocess_desurv.R")
if(!exists("predict_cox_on_yf"))source("code/predict_cox_on_yf.R")

# Nested feature sets ----

#' Find cohort-screening sizes that reach requested shared-feature counts
#' Uses the existing combined mean/variance ranks on log-transformed training
#' matrices before standardization. A gene enters the intersection at the
#' largest of its cohort rank positions. Ties can cross a target by more than1;
#' the actual count is reported, never silently trimmed or re-ranked.
#' @param cohort_matrices Named list of finite patients x common-genes matrices.
#' @param genes Unique common gene identifiers in matrix-column order.
#' @param targets Increasing requested shared counts; first equals the baseline.
#' @param baseline_top_n Existing per-cohort screening size for the first set.
#' @return List with design data.frame, gene-vector list and named entry ranks.
#'   Unavailable targets have NA counts and NULL gene vectors.
#' @examples
#' # feature_count_manifest(cohort_matrices,genes,c(2064,5000,10000,17900),3000)
#' @family feature_sensitivity
#' @seealso select_top_variable_genes
feature_count_manifest <- function(cohort_matrices,genes,targets,baseline_top_n) {
  p <- length(genes)
  if(!length(cohort_matrices) || is.null(names(cohort_matrices)) || p<1L ||
     anyNA(genes) || anyDuplicated(genes) || any(!nzchar(genes)) ||
     !length(targets) || any(!is.finite(targets)) || any(targets<1) ||
     any(targets!=as.integer(targets)) || any(diff(targets)<0) ||
     length(baseline_top_n)!=1L || !is.finite(baseline_top_n) ||
     baseline_top_n<1L || baseline_top_n>p || baseline_top_n!=as.integer(baseline_top_n))
    stop("Valid named cohorts, unique genes and increasing integer feature counts required.")
  ranks <- matrix(NA_integer_,p,length(cohort_matrices))
  for(c in seq_along(cohort_matrices)) {
    x <- cohort_matrices[[c]]
    if(!is.matrix(x) || nrow(x)<2L || !identical(colnames(x),genes) || any(!is.finite(x)))
      stop("Cohort matrices must be finite with matching common-gene columns.")
    # Same rank definition/tie order as select_top_variable_genes(combined_rank).
    score <- rank(-colMeans(x),ties.method="average")+
      rank(-apply(x,2,stats::var),ties.method="average")
    ranks[,c] <- match(seq_len(p),order(score,na.last=NA))
  }
  entry <- apply(ranks,1,max);names(entry) <- genes
  q <- actual <- rep(NA_integer_,length(targets));sets <- vector("list",length(targets))
  for(i in seq_along(targets)) {
    if(targets[i]>p)next
    q[i] <- if(i==1L)as.integer(baseline_top_n) else as.integer(sort(entry)[targets[i]])
    selected <- vector("list",length(cohort_matrices))
    for(c in seq_along(cohort_matrices))selected[[c]] <-
      select_top_variable_genes(cohort_matrices[[c]],genes,q[i],"combined_rank")$gene_names
    sets[[i]] <- Reduce(intersect,selected);actual[i] <- length(sets[[i]])
    if(!setequal(sets[[i]],genes[entry<=q[i]]))stop("Rank/selector mismatch.")
    if(i==1L && actual[i]!=targets[i])stop("Requested baseline differs from existing cohort screening.")
    if(i>1L && !all(sets[[i-1L]]%in%sets[[i]]))stop("Feature sets are not nested.")
  }
  list(design=data.frame(target_shared=as.integer(targets),screen_n_per_cohort=q,
    actual_shared=actual,maximum_available=p),genes=sets,entry_rank=entry)
}

# Correlations ----

#' Calculate correlations while preserving undefined constant-column cases
#' @param a,b Finite matrices with the same patient/feature rows and at least2 rows.
#' @return Column-by-column Pearson correlation matrix, NA for constant columns.
#' @examples
#' feature_cor(cbind(1:4,0),matrix(c(2,1,4,3),4,1))
#' @family feature_sensitivity
#' @seealso feature_program_agreement
feature_cor <- function(a,b) {
  a <- as.matrix(a);b <- as.matrix(b)
  if(nrow(a)!=nrow(b) || nrow(a)<2L || any(!is.finite(a)) || any(!is.finite(b)))
    stop("Finite matrices on matching rows required for correlation.")
  ac <- sweep(a,2,colMeans(a),"-");bc <- sweep(b,2,colMeans(b),"-")
  den <- outer(sqrt(colSums(ac^2)),sqrt(colSums(bc^2)))
  r <- crossprod(ac,bc)/den;r[den==0] <- NA_real_
  r[] <- pmax(-1,pmin(1,r));r
}

# Program agreement ----

#' Compare programs at two feature counts on the same ordered patients
#' Greedy non-reused absolute activity-correlation matching follows the existing
#' simulation matcher. Undefined entries cannot be selected; dead factors stay
#' unmatched. Loading correlations use only shared gene identities. Beta is
#' reported per SD of its correctly normalized projection, not in raw units.
#' @param old_fit,new_fit Existing YFB fits with EL,EF,EBeta,EF_norms fields.
#' @param old_genes,new_genes Unique gene identifiers in the corresponding EF order.
#' @param old_Y,new_Y Matrices on identical ordered patients and each fit's genes.
#' @return List of matched/unmatched pairs, all activity/loading/projection
#'   correlation matrices and one risk correlation. No model or selection is run.
#' @examples
#' # feature_program_agreement(old_fit,old_genes,new_fit,new_genes,old_Y,new_Y)
#' @family feature_sensitivity
#' @seealso feature_cor, predict_cox_on_yf
feature_program_agreement <- function(old_fit,old_genes,new_fit,new_genes,old_Y,new_Y) {
  if(anyNA(old_genes) || anyNA(new_genes) || anyDuplicated(old_genes) ||
     anyDuplicated(new_genes) || length(old_genes)!=nrow(old_fit$EF) ||
     length(new_genes)!=nrow(new_fit$EF) || nrow(old_Y)!=nrow(new_Y) ||
     nrow(old_Y)!=nrow(old_fit$EL) || nrow(new_Y)!=nrow(new_fit$EL) ||
     ncol(old_Y)!=length(old_genes) || ncol(new_Y)!=length(new_genes) ||
     !identical(rownames(old_Y),rownames(new_Y)))
    stop("Fits, patient order and unique gene identities must be aligned.")
  r <- feature_cor(old_fit$EL,new_fit$EL)
  common <- intersect(old_genes,new_genes)
  rf <- matrix(NA_real_,ncol(old_fit$EF),ncol(new_fit$EF))
  if(length(common)>=2L)rf <- feature_cor(
    old_fit$EF[match(common,old_genes),,drop=FALSE],
    new_fit$EF[match(common,new_genes),,drop=FALSE])
  a <- predict_cox_on_yf(old_Y,old_fit$EF,old_fit$EBeta,old_fit$EF_norms)
  b <- predict_cox_on_yf(new_Y,new_fit$EF,new_fit$EBeta,new_fit$EF_norms)
  rz <- feature_cor(a$ZF_test,b$ZF_test)
  old_units <- old_fit$EBeta*apply(a$ZF_test,2,stats::sd)
  new_units <- new_fit$EBeta*apply(b$ZF_test,2,stats::sd)
  score <- abs(r);score[!is.finite(score)] <- -Inf
  match <- rep(NA_integer_,nrow(r))
  while(any(is.finite(score))) {
    chosen <- which(score==max(score),arr.ind=TRUE)[1,]
    match[chosen[1]] <- chosen[2]
    score[chosen[1],] <- -Inf;score[,chosen[2]] <- -Inf
  }
  new_only <- setdiff(seq_len(ncol(r)),match[!is.na(match)])
  pairs <- data.frame(old_factor=c(seq_len(nrow(r)),rep(NA_integer_,length(new_only))),
    new_factor=c(match,new_only),activity_correlation=NA_real_,loading_correlation=NA_real_,
    projection_correlation=NA_real_,old_beta_per_projection_sd=NA_real_,
    new_beta_per_projection_sd=NA_real_,activity_best_margin=NA_real_)
  for(i in seq_len(nrow(pairs))) {
    o <- pairs$old_factor[i];n <- pairs$new_factor[i]
    if(!is.na(o)) {
      pairs$old_beta_per_projection_sd[i] <- old_units[o]
      ranked <- sort(abs(r[o,is.finite(r[o,])]),decreasing=TRUE)
      if(length(ranked)>=2L)pairs$activity_best_margin[i] <- ranked[1]-ranked[2]
    }
    if(!is.na(n))pairs$new_beta_per_projection_sd[i] <- new_units[n]
    if(!is.na(o) && !is.na(n)) {
      pairs$activity_correlation[i] <- r[o,n];pairs$loading_correlation[i] <- rf[o,n]
      pairs$projection_correlation[i] <- rz[o,n]
    }
  }
  list(pairs=pairs,activity_correlation=r,loading_correlation=rf,projection_correlation=rz,
    risk_correlation=as.numeric(feature_cor(a$risk_scores,b$risk_scores)))
}
