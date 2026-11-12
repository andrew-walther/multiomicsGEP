# ============================================================
# Script: score_lb_yfb_external.R
# Purpose: Score fixed LB and YFB fits on the S-12 external cohorts.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, yaml and the study prediction and preprocessing scripts
# ============================================================

#' Score fixed fits on one prepared external matrix
#' @param Y Patients by training genes, standardized within the external series.
#' @param fits Named list containing fixed LB and YFB fits.
#' @return Named list of Cox-oriented risk vectors.
#' @examples
#' # lb_yfb_external_predict(Y, fits)
lb_yfb_external_predict <- function(Y, fits) {
  lb <- predict_ssbmf_lb_foldin(Y, fits$LB)
  if (!isTRUE(lb$converged)) stop('LB external fold-in did not converge')
  yf <- predict_ssbmf_yfb(Y, fits$YFB$EF, fits$YFB$EBeta,
    EF_norms = fits$YFB$EF_norms, projection_center = fits$YFB$projection_center)
  list(LB = as.numeric(lb$risk_scores), YFB = as.numeric(yf$risk_scores))
}

#' Align a preprocessed external matrix to training genes
#' @param Y Standardized patients by measured genes. @param measured Gene names.
#' @param genes Training gene names in fitted order.
#' @return List containing aligned matrix and missing gene count.
#' @examples
#' # lb_yfb_external_align(matrix(1, 2, 1), 'a', c('a', 'b'))
lb_yfb_external_align <- function(Y, measured, genes) {
  if (anyDuplicated(measured) || anyDuplicated(genes)) stop('Duplicate gene names')
  hit <- match(genes, measured)
  out <- matrix(0, nrow(Y), length(genes))
  out[, !is.na(hit)] <- Y[, hit[!is.na(hit)], drop = FALSE]
  list(Y = out, missing = sum(is.na(hit)))
}

#' Check audited cross-assay donor overlap
#' @param donor_list Named donor vectors for external series.
#' @param expected Required shared PACA-AU donor count.
#' @return Number of shared donors.
#' @examples
#' # lb_yfb_external_overlap(list(PACA_AU_array = 'a', PACA_AU_seq = 'a'), 1)
lb_yfb_external_overlap <- function(donor_list, expected = 47L) {
  if (!all(c('PACA_AU_array', 'PACA_AU_seq') %in% names(donor_list))) stop('Missing PACA-AU series')
  shared <- length(intersect(donor_list$PACA_AU_array, donor_list$PACA_AU_seq))
  if (shared != expected) stop('PACA-AU donor overlap differs from audit: ', shared)
  shared
}

#' Compute equal-series and within-series pair-weighted concordance
#' @param keep Named series list with LB, YFB, time, status and donor vectors.
#' @param indices Optional list of row indices for a bootstrap sample.
#' @return Named vector of model concordances and paired differences.
#' @examples
#' # lb_yfb_external_pooled_stat(keep)
lb_yfb_external_pooled_stat <- function(keep, indices = lapply(keep, function(x) seq_along(x$time))) {
  scores <- lapply(seq_along(keep), function(i) {
    x <- keep[[i]]; idx <- indices[[i]]
    if (length(idx) < 2 || sum(x$status[idx]) < 2) stop('Bootstrap stratum has fewer than two events')
    lapply(c('LB', 'YFB'), function(model) {
      risk <- x[[model]][idx]
      if (length(unique(risk)) < 2) stop('Constant risk in pooled bootstrap')
      fit <- survival::concordance(survival::Surv(x$time[idx], x$status[idx]) ~ risk, reverse = TRUE)
      pairs <- sum(fit$count[1:3])
      if (!is.finite(fit$concordance) || pairs <= 0) stop('Undefined pooled concordance')
      c(C = unname(fit$concordance), pairs = unname(pairs))
    })
  })
  C <- do.call(rbind, lapply(scores, function(x) vapply(x, function(y) y['C'], numeric(1))))
  pairs <- vapply(scores, function(x) x[[1]]['pairs'], numeric(1))
  equal <- colMeans(C); weighted <- colSums(C * pairs) / sum(pairs)
  c(equal_LB = equal[1], equal_YFB = equal[2], equal_diff = equal[2] - equal[1],
    pair_LB = weighted[1], pair_YFB = weighted[2], pair_diff = weighted[2] - weighted[1])
}

#' Donor bootstrap pooled model concordances and paired differences
#' @param keep Named series list with aligned scores, outcomes and donor IDs.
#' @param B Bootstrap replicates. @param seed Random seed.
#' @return Aggregate estimates, percentile limits and bootstrap standard errors.
#' @examples
#' # lb_yfb_external_pooled_ci(keep, B = 2000)
lb_yfb_external_pooled_ci <- function(keep, B = 2000L, seed = 1L) {
  donors <- lapply(keep, `[[`, 'donor')
  if (any(vapply(seq_along(keep), function(i) length(donors[[i]]) != length(keep[[i]]$time) ||
                 anyNA(donors[[i]]) || any(!nzchar(donors[[i]])), logical(1)))) stop('Invalid donor IDs')
  rows <- lapply(donors, function(x) split(seq_along(x), as.character(x)))
  ids <- unique(unlist(donors, use.names = FALSE))
  membership <- vapply(ids, function(id) paste(which(vapply(rows, function(x) id %in% names(x), logical(1))), collapse = ','), character(1))
  groups <- split(ids, membership)
  estimate <- lb_yfb_external_pooled_stat(keep)
  set.seed(seed)
  boot <- replicate(B, {
    drawn <- unlist(lapply(groups, function(x) x[sample.int(length(x), length(x), replace = TRUE)]), use.names = FALSE)
    idx <- lapply(rows, function(x) unlist(x[drawn], use.names = FALSE))
    lb_yfb_external_pooled_stat(keep, idx)
  })
  data.frame(metric = names(estimate), estimate = unname(estimate),
    lower = apply(boot, 1, quantile, probs = .025), upper = apply(boot, 1, quantile, probs = .975),
    se = apply(boot, 1, sd), row.names = NULL)
}

#' Run external scoring and write aggregate outputs
#' @param checkpoint Fixed patient-level training checkpoint path.
#' @param input Training input path with ordered genes.
#' @param out Fresh aggregate output directory. @param B Bootstrap replicates.
#' @return Invisibly, aggregate tables.
#' @examples
#' # lb_yfb_external_main('data/checkpoint.rds', 'data/input.rds', '/tmp/external')
lb_yfb_external_main <- function(checkpoint, input, out, B = 2000L) {
  if (file.exists(out)) stop('Output path already exists')
  if (grepl('(^|/)data(/|$)', out)) stop('Aggregate output must be outside data/')
  cfg <- yaml::read_yaml('config/globals.yml')
  source('results/benchmark_sim/benchmark_helpers.R', local = TRUE)
  source('code/preprocess_desurv.R', local = TRUE)
  source('code/predict_ssbmf_lb.R', local = TRUE)
  source('code/predict_ssbmf_yfb.R', local = TRUE)
  source('code/concordance_ci.R', local = TRUE)
  # Prediction functions are resolved by the scoring helper in the global environment.
  for (path in c('code/predict_ssbmf_lb.R', 'code/predict_ssbmf_yfb.R')) source(path, local = .GlobalEnv)
  genes <- readRDS(input)$genes
  fits <- readRDS(checkpoint)$pair$fits
  if (is.null(genes) || anyDuplicated(genes) || nrow(fits$YFB$EF) != length(genes) || nrow(fits$LB$EF) != length(genes)) stop('Training genes do not match fits')
  keep <- list(); rows <- list(); diffs <- list()
  for (co in cfg$pdac$external_cohorts) {
    raw <- load_pdac_raw(co, PDAC_DATA_ROOT)
    pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
      log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
      rank_transform = FALSE, per_platform_standardize = TRUE)
    aligned <- lb_yfb_external_align(pre$Y, pre$gene_names, genes)
    risk <- lb_yfb_external_predict(aligned$Y, fits)
    ids <- as.character(raw$sampID)
    if (length(ids) != raw$n || anyNA(ids) || any(!nzchar(ids))) stop('Invalid sample IDs in ', co)
    if (co %in% c('PACA_AU_array', 'PACA_AU_seq')) {
      meta <- readRDS(file.path(PDAC_DATA_ROOT, 'original', paste0(co, '.rds')))$sampInfo
      if (!all(c('icgc_sample_id', 'icgc_donor_id') %in% names(meta)) || anyDuplicated(meta$icgc_sample_id)) stop('Invalid donor metadata')
      donor <- as.character(meta$icgc_donor_id[match(sub('_seq$', '', ids), meta$icgc_sample_id)])
      if (anyNA(donor) || any(!nzchar(donor))) stop('Unmatched donors in ', co)
      donor <- paste0('PACA_AU:', donor)
    } else donor <- paste0(co, ':', ids)
    keep[[co]] <- c(risk, list(time = raw$time, status = raw$status, donor = donor))
    for (model in c('LB', 'YFB')) {
      ci <- bootstrap_concordance_ci(risk[[model]], raw$time, raw$status, B = B, seed = 1, flip = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(cohort = co, model = model, n = raw$n,
        events = sum(raw$status), missing_genes = aligned$missing, C = ci$estimate, lower = ci$lower, upper = ci$upper)
    }
    ci <- bootstrap_concordance_diff_ci(risk$YFB, risk$LB, raw$time, raw$status,
      B = B, seed = 1, flip_a = FALSE, flip_b = FALSE)
    diffs[[co]] <- data.frame(cohort = co, difference = ci$estimate, lower = ci$lower, upper = ci$upper, se = ci$se)
  }
  shared <- lb_yfb_external_overlap(lapply(keep, `[[`, 'donor'))
  pooled <- lb_yfb_external_pooled_ci(keep, B = B)
  dir.create(out, recursive = TRUE)
  tables <- list(external_cindex = do.call(rbind, rows), external_paired_diff = do.call(rbind, diffs), external_pooled = pooled)
  for (name in names(tables)) write.csv(tables[[name]], file.path(out, paste0(name, '.csv')), row.names = FALSE)
  writeLines(c(paste0('checkpoint=', checkpoint), paste0('input=', input), paste0('B=', B),
    paste0('shared_paca_au_donors=', shared), 's12_score_correlation=unavailable_aggregate_scores',
    paste0('git_head=', system2('git', c('rev-parse', 'HEAD'), stdout = TRUE))), file.path(out, 'settings.txt'))
  invisible(tables)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (any(!grepl('^--(checkpoint|input|out)=.+$', args))) stop('Use --checkpoint=, --input= or --out=')
  keys <- sub('^--([^=]+)=.*', '\\1', args)
  if (anyDuplicated(keys)) stop('Duplicate option')
  opts <- as.list(setNames(sub('^--[^=]+=', '', args), keys))
  day <- format(Sys.Date(), '%Y-%m-%d')
  if (is.null(opts$checkpoint)) opts$checkpoint <- file.path('data', paste0('lb_yfb_study_lb_yfb_training_', day), 'patient_level.rds')
  if (is.null(opts$input)) opts$input <- 'data/lf_prior_real_2026-11-12.rds'
  if (is.null(opts$out)) opts$out <- file.path('results/benchmark_sim/outputs', paste0('lb_yfb_external_', day))
  lb_yfb_external_main(opts$checkpoint, opts$input, opts$out)
}
