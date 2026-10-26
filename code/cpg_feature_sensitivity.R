# ============================================================
# Script: cpg_feature_sensitivity.R
# Purpose: Preserve original CpG ranking and exact frozen training inputs during expansion.
# Created: 2026-10-15
# Dependencies: base R
# ============================================================

# CpG feature sets ----

#' Build nested CpG counts with the original positive-variance ranking
#' Uses the existing screen_multiomics_features() variance/sort convention,
#' after the original whole-universe imputation. No new mean or SD filter.
#' @param Y Finite patients x CpGs matrix with unique probe column names.
#' @param counts Increasing positive integer requested probe counts.
#' @param frozen_features Frozen variance-ranked probe names at the first count.
#' @return List of count design, ordered feature vectors and full eligible rank.
#'   Unavailable counts are explicit NA/NULL, never silently truncated.
#' @examples
#' # cpg_count_manifest(Y,c(10000L,25000L,50000L),frozen_names)
#' @family feature_sensitivity
#' @seealso screen_multiomics_features, validate_cpg_reconstruction
cpg_count_manifest <- function(Y,counts,frozen_features) {
  if(!is.matrix(Y) || nrow(Y)<2L || any(!is.finite(Y)) ||
     is.null(colnames(Y)) || anyNA(colnames(Y)) || anyDuplicated(colnames(Y)) ||
     any(!nzchar(colnames(Y))) || !length(counts) || any(!is.finite(counts)) ||
     any(counts<1) || any(counts!=as.integer(counts)) || any(diff(counts)<=0) ||
     length(frozen_features)!=counts[1])stop("Finite named CpGs and increasing integer counts required.")
  variance <- apply(Y,2,stats::var)
  ranked <- names(sort(variance[variance>0],decreasing=TRUE))
  if(!identical(head(ranked,counts[1]),frozen_features))
    stop("Variance ranking does not reproduce the frozen CpG screen.")
  features <- lapply(counts,function(n)if(n<=length(ranked))ranked[seq_len(n)] else NULL)
  actual <- vapply(features,function(x)if(is.null(x))NA_integer_ else length(x),integer(1))
  list(design=data.frame(requested=as.integer(counts),actual=actual,eligible=length(ranked)),
    features=features,ranked=ranked)
}

#' Require exact frozen CpG values before expanding a training feature set
#' @param Y Finite full reconstructed patients x CpGs matrix, uniquely named.
#' @param frozen Finite archived subset on identical ordered named patients.
#' @return TRUE invisibly; errors on any value, identity or order mismatch.
#' @examples
#' # validate_cpg_reconstruction(full_Y,cached_training$Y$methylation)
#' @family feature_sensitivity
#' @seealso cpg_count_manifest
validate_cpg_reconstruction <- function(Y,frozen) {
  for(x in list(Y,frozen)) {
    if(!is.matrix(x) || any(!is.finite(x)) || is.null(rownames(x)) ||
       is.null(colnames(x)) || anyNA(rownames(x)) || anyNA(colnames(x)) ||
       anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
       any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
       !identical(rownames(Y),rownames(frozen)))
      stop("Finite uniquely named matrices on the identical patient order required.")
  }
  if(!all(colnames(frozen)%in%colnames(Y)) ||
     !identical(unname(Y[,match(colnames(frozen),colnames(Y)),drop=FALSE]),unname(frozen)))
    stop("Reconstructed frozen CpG values differ; stop before fitting.")
  invisible(TRUE)
}
