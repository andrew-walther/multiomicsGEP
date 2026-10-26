# ============================================================
# Script: desurv_author_adapter.R
# Purpose: Preserve the authors' preprocessing and projection on matched patients.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: yaml; DeSurv 1.0.1, commit cc9d4c2486b3125704ab71eea61f35e040a40080
# Paper pipeline: rashidlab/DeSurv-paper, commit 0a870a78673296a3802aeddfee0514e4a8b0f9ae
# ============================================================
# The original loader already applies the authors' platform-specific log2(x+1).
# Never pass D4's centered matrix here or apply another log transformation.
# The union/zero-fill below follows paper R/load_data.R and preprocess_helpers.R;
# the actual filtering and within-sample ranks are DeSurv::preprocess_data().

# Author data preparation ----

#' Merge author-loaded cohorts without feature screening or transformation
#' @param cohorts Named list returned by the authors' raw data loader.
#' @return List with raw genes-by-patient ex and aligned eligible sampInfo.
#' @examples
#' # desurv_author_merge(list(TCGA=tcga_raw,CPTAC=cptac_raw))
#' @family desurv_author_pipeline
#' @seealso desurv_author_prepare
desurv_author_merge <- function(cohorts) {
  for(x in cohorts) {
    if(!identical(colnames(x$ex),as.character(x$sampInfo$ID)))
      stop("Author expression columns and survival IDs are not aligned")
  }
  meta <- do.call(rbind,lapply(cohorts,function(x)x$sampInfo[x$sampInfo$keep==1,,drop=FALSE]))
  if(length(cohorts)>1L) {
    # Match the authors' outer join, including last assignment for duplicate IDs.
    universe <- Reduce(union,lapply(cohorts,`[[`,"featInfo"))
    blocks <- lapply(cohorts,function(x) {
      keep <- which(x$sampInfo$keep==1)
      z <- matrix(0,length(universe),length(keep),dimnames=list(universe,colnames(x$ex)[keep]))
      z[x$featInfo,] <- x$ex[,keep,drop=FALSE]
      z
    })
    X <- do.call(cbind,blocks)
  } else {
    x <- cohorts[[1L]]
    X <- x$ex[,x$sampInfo$keep==1,drop=FALSE]
  }
  list(ex=X,sampInfo=meta)
}

#' Prepare author-loaded cohorts with DeSurv's unchanged preprocessing
#' @param cohorts Named list of author loader results: ex (genes x patients),
#'   featInfo and sampInfo (ID, time, event, dataset, keep).
#' @param genes NULL for training top3000-per-cohort intersection; frozen
#'   training genes for validation, with missing genes zero-filled before ranking.
#' @return DeSurv preprocessing result, with the exact retained patient IDs/order.
#' @examples
#' # desurv_author_prepare(list(TCGA=tcga_raw,CPTAC=cptac_raw))
#' @family desurv_author_pipeline
#' @seealso DeSurv::preprocess_data
desurv_author_prepare <- function(cohorts,genes=NULL) {
  merged <- desurv_author_merge(cohorts)
  X <- merged$ex;meta <- merged$sampInfo
  if(!is.null(genes)) {
    missing <- setdiff(genes,rownames(X))
    if(length(missing)) {
      message("DeSurv validation: zero-filling ",length(missing),"/",length(genes)," training genes before ranking")
      X <- rbind(X,matrix(0,length(missing),ncol(X),dimnames=list(missing,colnames(X))))
    }
  }
  ngene <- yaml::read_yaml("config/globals.yml")$preprocessing$top_n_genes_desurv
  DeSurv::preprocess_data(X,meta$time,meta$event,meta$dataset,ngene=ngene,
    genes=genes,method_trans_train="rank",verbose=FALSE)
}

# Frozen validation projection ----

#' Project the union of top genes exactly as paper code/05_external_validation.R
#' @param fit DeSurv fit with gene-named p x K W and K coefficients beta.
#' @param X Nonnegative, author-preprocessed p x n validation matrix.
#' @param ntop Reported training-selected top-gene count per program.
#' @return Numeric n-vector of Cox linear predictors; no external reorientation.
#' @examples
#' # desurv_author_risk(fit,validation$ex,270L)
#' @family desurv_author_pipeline
#' @seealso DeSurv::desurv_get_top_genes
desurv_author_risk <- function(fit,X,ntop) {
  genes <- unique(unlist(DeSurv::desurv_get_top_genes(fit$W,ntop)$top_genes))
  if(!length(genes) || any(!genes %in% rownames(X)))
    stop("Frozen DeSurv projection genes are absent")
  as.vector(t(X[genes,,drop=FALSE])%*%fit$W[genes,,drop=FALSE]%*%fit$beta)
}
