# ============================================================

source("code/simulate_multimodal_yfb.R", local = TRUE)
# Script: share_targeted_structure.R
# Purpose: Summarize training program structure and signature enrichment.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, fgsea, code/pathway_enrichment.R,
#   code/simulate_multimodal_yfb.R
# ============================================================

#' Match nonconstant loadings one-to-one, retaining original indices
#' @param ref,estimate Feature-by-program loading matrices.
#' @return Match table, or NULL if either matrix has no variable columns.
#' @examples
#' # s16_structure_matches(diag(2), diag(2))
s16_structure_matches <- function(ref, estimate) {
  i <- which(apply(ref, 2, stats::sd) > 0)
  j <- which(apply(estimate, 2, stats::sd) > 0)
  if (!length(i) || !length(j)) return(NULL)
  m <- multimodal_yfb_match_factors(ref[, i, drop = FALSE], estimate[, j, drop = FALSE])
  m$true_factor <- i[m$true_factor]
  m$estimated_factor <- j[m$estimated_factor]
  m
}

#' Compute program structure and fold-stability rule for one real-data arm
#' @param full Full-training real result.
#' @param folds List of five held-out-fold real results.
#' @param d Prepared training data.
#' @param arm Arm name.
#' @return Program, correlation and fold-match tables.
#' @examples
#' # s16_structure_arm(full, folds, d, "w1")
s16_structure_arm <- function(full, folds, d, arm) {
  fit <- full$arms[[arm]]$fit; ref <- full$arms$twostep$fit
  if (is.null(fit) || is.null(ref)) stop("missing full fit or two-step reference")
  Y <- sweep(d$Y, 2, full$center)
  projection <- if (is.null(fit$EF_norms)) fit$EF else sweep(fit$EF, 2, fit$EF_norms, "/")
  Z <- Y %*% projection
  component <- sweep(Z, 2, fit$EBeta, "*"); eta <- rowSums(component)
  se <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  active <- se > 0 & abs(fit$EBeta) >= 1.96 * se & colSums(fit$EF^2) > 1e-8
  sdcomponent <- apply(component, 2, stats::sd)
  ratio <- if (max(sdcomponent) > 0) sdcomponent / max(sdcomponent) else rep(NA_real_, ncol(Z))
  allocation <- if (stats::var(eta) > 0) as.numeric(stats::cov(component, eta)) / stats::var(eta) else rep(NA_real_, ncol(Z))
  correlations <- abs(stats::cor(fit$EF, ref$EF))
  matches <- s16_structure_matches(fit$EF, ref$EF)
  if (length(folds) != 5L) stop("exactly five fold results are required")
  fold_matches <- do.call(rbind, lapply(seq_along(folds), function(i) {
    # One row per full program makes an unavailable fold explicitly unknown.
    rows <- data.frame(fold = i, true_factor = seq_len(ncol(fit$EF)),
      estimated_factor = NA_integer_, correlation = NA_real_,
      abs_correlation = NA_real_, active_match = NA)
    ff <- folds[[i]]$arms[[arm]]$fit
    if (is.null(ff) || (!is.null(ff$history$converged) && !isTRUE(ff$history$converged))) return(rows)
    m <- s16_structure_matches(fit$EF, ff$EF)
    # An observed fit without a match cannot meet the matching criterion.
    rows$active_match <- FALSE
    if (is.null(m)) return(rows)
    fs <- sqrt(pmax(ff$EBeta2 - ff$EBeta^2, 0))
    fa <- fs > 0 & abs(ff$EBeta) >= 1.96 * fs & colSums(ff$EF^2) > 1e-8
    idx <- m$true_factor
    rows$estimated_factor[idx] <- m$estimated_factor
    rows$correlation[idx] <- m$correlation
    rows$abs_correlation[idx] <- m$abs_correlation
    rows$active_match[idx] <- m$abs_correlation >= .7 & fa[m$estimated_factor]
    rows
  }))
  full_cox <- survival::coxph(survival::Surv(d$time, d$status) ~ Z,
    control = survival::coxph.control(timefix = FALSE))
  if (any(!is.finite(stats::coef(full_cox)))) stop("full program Cox fit not estimable")
  log_full <- as.numeric(stats::logLik(full_cox))
  recon <- fit$EL %*% t(fit$EF); residual <- sum((Y - recon)^2)
  program <- do.call(rbind, lapply(seq_len(ncol(Z)), function(k) {
    keep <- setdiff(seq_len(ncol(Z)), k)
    reduced <- if (length(keep)) {
      Zdrop <- Z[, keep, drop = FALSE]
      survival::coxph(survival::Surv(d$time, d$status) ~ Zdrop,
        control = survival::coxph.control(timefix = FALSE))
    } else survival::coxph(survival::Surv(d$time, d$status) ~ 1,
                            control = survival::coxph.control(timefix = FALSE))
    if (any(!is.finite(stats::coef(reduced)))) stop("reduced program Cox fit not estimable")
    fold_status <- fold_matches$active_match[fold_matches$true_factor == k]
    nfold <- sum(fold_status, na.rm = TRUE)
    unknown <- sum(is.na(fold_status))
    base <- active[k] && !is.na(ratio[k]) && ratio[k] >= .25 &&
      !is.na(allocation[k]) && allocation[k] >= .05
    main <- if (!base) FALSE else if (nfold >= 4) TRUE else if (nfold + unknown < 4) FALSE else NA
    drop_recon <- recon - tcrossprod(fit$EL[, k], fit$EF[, k])
    data.frame(arm = arm, program = k, active = active[k], component_ratio = ratio[k],
      signed_risk_share = allocation[k], active_folds = nfold,
      unknown_folds = unknown, main_program = main,
      delta_logPL = log_full - as.numeric(stats::logLik(reduced)),
      semipartial_R2 = (sum((Y - drop_recon)^2) - residual) / sum(Y^2))
  }))
  list(program = program, correlations = correlations, matches = matches, fold_matches = fold_matches)
}

#' Load existing tumor and fibroblast gene signatures
#' @param pdac_cache Existing cached PDAC gene-set RDS.
#' @param subtype_file Local cmbSubtypes.RData containing human Elyada CAF sets.
#' @return Named list of Classical, basal-like, iCAF and myCAF gene sets.
#' @examples
#' # sets <- s16_structure_genesets("path/pdac_genesets.rds", "path/cmbSubtypes.RData")
s16_structure_genesets <- function(pdac_cache, subtype_file) {
  sets <- readRDS(pdac_cache)
  names_needed <- c("Moffitt_Classical", "Moffitt_BasalLike")
  if (!all(names_needed %in% names(sets))) stop("Moffitt signatures absent from cache")
  e <- new.env(parent = emptyenv())
  load(subtype_file, envir = e)
  idx <- which(e$schemaList == "Elyada_CAF")
  if (length(idx) != 1L) stop("human Elyada_CAF signature unavailable")
  tab <- e$subtypeGeneList[[idx]]
  if (!all(c("geneSymbol", "iCAF", "myCAF") %in% names(tab))) stop("invalid Elyada_CAF schema")
  out <- c(sets[names_needed], list(Elyada_iCAF = tab$geneSymbol[which(tab$iCAF)],
                                     Elyada_myCAF = tab$geneSymbol[which(tab$myCAF)]))
  if (any(vapply(out, length, integer(1)) == 0L)) stop("empty signature")
  out
}

#' Annotate loading programs against fixed training gene signatures
#' @param fit Full-training fitted state containing EF.
#' @param genes Gene symbols in loading-row order, unique and nonempty.
#' @param genesets Named signatures from s16_structure_genesets().
#' @return Ranked enrichment table for each program and signature.
#' @examples
#' # s16_structure_annotation(fit, genes, sets)
#' @seealso run_fgsea_program
s16_structure_annotation <- function(fit, genes, genesets) {
  if (length(genes) != nrow(fit$EF) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes))
    stop("genes must uniquely name every loading row")
  if (!length(genesets) || is.null(names(genesets))) stop("named signatures required")
  overlap <- vapply(genesets, function(x) length(intersect(x, genes)), integer(1))
  if (any(overlap < 10L)) stop("signature has fewer than 10 measured genes: ", paste(names(overlap)[overlap < 10L], collapse = ", "))
  do.call(rbind, lapply(seq_len(ncol(fit$EF)), function(k) {
    weights <- stats::setNames(abs(fit$EF[, k]), genes)
    if (!any(weights > 0)) stop("zero loading program cannot be annotated: ", k)
    result <- run_fgsea_program(weights, genesets, seed = 1, collection = "PDAC signatures")
    result$program <- k
    result
  }))
}
