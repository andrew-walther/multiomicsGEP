# ============================================================
# Script: summarize_lb_yfb_study.R
# Purpose: Summarize registered LB/YFB simulations and training fits.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm and project study and enrichment helpers
# ============================================================

#' Summarize each simulation cell with explicit denominators
#' The registered gate requires recovery >= 8/10, false activity <= 1 and
#' convergence >= 9/10. Relative exclusion requires the other model to pass.
#' @param path Aggregate per-seed CSV, including explicit failure rows.
#' @return Per-cell aggregates and the registered relative exclusion decision.
#' @examples
#' # lb_yfb_summary_simulation('simulation_metrics.csv')
lb_yfb_summary_simulation <- function(path) {
  # Force message storage to character: read.csv otherwise turns an all-blank
  # success column into logical NA, which nzchar() treats as a nonempty string.
  metrics <- read.csv(path, colClasses = c(error = 'character'))
  keys <- c('design', 'share', 'amp', 'w_S', 'model')
  labels <- lapply(metrics[keys], function(x) ifelse(is.na(x), '<NA>', as.character(x)))
  groups <- split(seq_len(nrow(metrics)), do.call(interaction, c(labels, list(drop = TRUE, lex.order = TRUE))))
  rows <- lapply(groups, function(idx) {
    x <- metrics[idx, , drop = FALSE]
    expected <- if (x$design[1] == 'matched') 6L else 10L
    if (anyDuplicated(x$seed)) stop('Duplicate seed within simulation cell')
    #' Average defined descriptive values while reporting the denominator below
    #' @param z Numeric per-seed measurements.
    #' @return Mean of finite values or NA when none is defined.
    #' @examples
    #' # finite_mean(c(1, NA_real_))
    finite_mean <- function(z) if (any(is.finite(z))) mean(z[is.finite(z)]) else NA_real_
    complete <- setequal(x$seed, seq_len(expected)) && !any(nzchar(x$error))
    gate <- if (x$design[1] == 'matched') NA else
      complete && all(is.finite(x$false_active)) && all(!is.na(x$recovery_success)) &&
      sum(x$recovery_success) >= 8 && mean(x$false_active) <= 1 && sum(x$converged, na.rm = TRUE) >= 9
    data.frame(x[1, keys, drop = FALSE], expected_seeds = expected, observed_seeds = nrow(x),
      failed = sum(nzchar(x$error)), recovery_successes = sum(x$recovery_success, na.rm = TRUE),
      recovery_defined = sum(!is.na(x$recovery_success)), false_active_mean = finite_mean(x$false_active),
      false_active_defined = sum(is.finite(x$false_active)), converged = sum(x$converged, na.rm = TRUE),
      validation_C_mean = finite_mean(x$validation_C), validation_C_defined = sum(is.finite(x$validation_C)),
      matched_score_mean = finite_mean(x$matched_score_recovery), matched_score_defined = sum(is.finite(x$matched_score_recovery)),
      gate_pass = gate)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out$relative_exclusion <- NA
  for (i in seq_len(nrow(out))) {
    if (is.na(out$gate_pass[i])) next
    other <- which(out$design == out$design[i] &
      (is.na(out$share) & is.na(out$share[i]) | !is.na(out$share) & !is.na(out$share[i]) & out$share == out$share[i]) &
      (is.na(out$amp) & is.na(out$amp[i]) | !is.na(out$amp) & !is.na(out$amp[i]) & out$amp == out$amp[i]) &
      out$w_S == out$w_S[i] & out$model != out$model[i])
    if (length(other) == 1L && !is.na(out$gate_pass[other]))
      out$relative_exclusion[i] <- !out$gate_pass[i] && out$gate_pass[other]
  }
  out
}

#' Describe program coefficients and loading activity for a fit
#' The plug-in hazard ratio is exp(E[β_k]), not the posterior mean of exp(β_k).
#' Activity requires a positive posterior SD and |E[β_k]| >= 1.96 SD(β_k).
#' @param fit One fitted model. @param model Model label. @param fold Fold label or full.
#' @return Per-program coefficient, uncertainty, hazard ratio and activity.
#' @examples
#' # lb_yfb_summary_programs(fit, 'LB', 'full')
lb_yfb_summary_programs <- function(fit, model, fold) {
  variance <- if (!is.null(fit$VBeta)) fit$VBeta else pmax(fit$EBeta2 - fit$EBeta^2, 0)
  sd <- sqrt(variance)
  data.frame(fold = fold, model = model, program = seq_along(fit$EBeta),
    beta = as.numeric(fit$EBeta), beta_sd = as.numeric(sd), hazard_ratio = exp(as.numeric(fit$EBeta)),
    nonzero_loading = colSums(abs(fit$EF)) > 0,
    survival_active = sd > 0 & abs(fit$EBeta) >= 1.96 * sd)
}

#' Rank nonzero program loadings by absolute magnitude
#' @param fit One fitted model. @param genes Ordered training gene symbols.
#' @param model Model label. @param n Maximum genes per program.
#' @return Signed top genes and counts including exactly dead programs.
#' @examples
#' # lb_yfb_summary_top(fit, genes, 'LB')
lb_yfb_summary_top <- function(fit, genes, model, n = 20L) {
  EF <- fit$EF
  if (nrow(EF) != length(genes) || anyDuplicated(genes)) stop('Invalid ordered genes')
  rownames(EF) <- genes
  live <- which(colSums(abs(EF)) > 0)
  counts <- data.frame(model = model, program = seq_len(ncol(EF)), n_top = integer(ncol(EF)))
  if (!length(live)) return(list(top = data.frame(model = character(), program = integer(), rank = integer(), gene = character(), weight = numeric()), counts = counts))
  labels <- as.list(setNames(as.character(seq_len(ncol(EF))), as.character(seq_len(ncol(EF)))))
  top <- top_n_genes_table(abs(EF), live, labels, n = n)
  top <- top[top$weight > 0, , drop = FALSE]
  counts$n_top <- tabulate(top$program, nbins = ncol(EF))
  top$weight <- EF[cbind(match(top$gene, genes), top$program)]
  list(top = data.frame(model = model, top[, c('program', 'rank', 'gene', 'weight')]), counts = counts)
}

#' Aggregate saved study outputs without writing patient-level information
#' @param date Run date used in the output and checkpoint paths.
#' @return Invisibly, aggregate output tables.
#' @examples
#' # lb_yfb_summary_main('2026-11-12')
lb_yfb_summary_main <- function(date = '2026-11-12') {
  if (!grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}$', date)) stop('Invalid date')
  source('code/pathway_enrichment.R', local = .GlobalEnv)
  source('code/lb_yfb_study.R', local = .GlobalEnv)
  source('code/predict_ssbmf_yfb.R', local = .GlobalEnv)
  simdir <- file.path('results/benchmark_sim/outputs', paste0('lb_yfb_simulation_', date))
  traindir <- file.path('results/benchmark_sim/outputs', paste0('lb_yfb_training_', date))
  checkpoint <- file.path('data', paste0('lb_yfb_study_lb_yfb_training_', date), 'patient_level.rds')
  input <- readRDS('data/lf_prior_real_2026-11-12.rds'); saved <- readRDS(checkpoint)
  simulation <- lb_yfb_summary_simulation(file.path(simdir, 'simulation_metrics.csv'))
  write.csv(simulation, file.path(simdir, 'simulation_summary.csv'), row.names = FALSE)
  programs <- list()
  for (j in seq_along(saved$cv$pairs)) {
    pair <- saved$cv$pairs[[j]]
    if (is.null(pair$fits)) stop('CV fold ', j, ' has no fits: ', pair$error)
    for (model in c('LB', 'YFB')) programs[[length(programs) + 1L]] <- lb_yfb_summary_programs(pair$fits[[model]], model, as.character(j))
  }
  for (model in c('LB', 'YFB')) programs[[length(programs) + 1L]] <- lb_yfb_summary_programs(saved$pair$fits[[model]], model, 'full')
  programs <- do.call(rbind, programs)
  activity <- aggregate(cbind(nonzero_loading, survival_active) ~ fold + model, programs, sum)
  top <- lapply(c('LB', 'YFB'), function(m) lb_yfb_summary_top(saved$pair$fits[[m]], input$genes, m))
  names(top) <- c('LB', 'YFB')
  top_table <- do.call(rbind, lapply(top, `[[`, 'top'))
  top_counts <- do.call(rbind, lapply(top, `[[`, 'counts'))
  desurv <- readRDS('results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds')
  desurv <- desurv[grep('^DeSurv_D[123]_', names(desurv))]
  if (length(desurv) != 3L || any(lengths(desurv) != 270L)) stop('Expected three published DeSurv 270-gene sets')
  sets <- lapply(top, function(x) split(x$top$gene, factor(x$top$program, levels = seq_len(nrow(x$counts)))))
  #' Compare gene sets in the fixed training-gene background
  #' @param a,b Gene-symbol vectors.
  #' @return Named vector of overlap count, Jaccard index and hypergeometric p-value.
  #' @examples
  #' # overlap(c("A", "B"), "A")
  overlap <- function(a, b) unlist(compute_geneset_overlap(a, b, input$genes))
  cross_grid <- expand.grid(LB_program = seq_along(sets$LB), YFB_program = seq_along(sets$YFB))
  across <- do.call(rbind, lapply(seq_len(nrow(cross_grid)), function(k) {
    i <- cross_grid$LB_program[k]; j <- cross_grid$YFB_program[k]
    data.frame(LB_program = i, YFB_program = j, as.list(overlap(sets$LB[[i]], sets$YFB[[j]])))
  }))
  comparator_grid <- expand.grid(model = names(sets), program = seq_along(sets$LB), desurv = names(desurv), stringsAsFactors = FALSE)
  comparator <- do.call(rbind, lapply(seq_len(nrow(comparator_grid)), function(k) {
    m <- comparator_grid$model[k]; i <- comparator_grid$program[k]; d <- comparator_grid$desurv[k]
    data.frame(model = m, program = i, desurv = d, top_n = length(sets[[m]][[i]]),
               as.list(overlap(sets[[m]][[i]], desurv[[d]])))
  }))
  lb <- saved$pair$fits$LB; yf <- saved$pair$fits$YFB
  Y <- lb_yfb_study_scale(input$X, input$cohort, scaling = 'zscore')
  yfb_score <- predict_ssbmf_yfb(Y, yf$EF, yf$EBeta, yf$EF_norms, projection_center = yf$projection_center)$ZF_test
  #' Correlate varying programs, preserving undefined dead-program comparisons
  #' @param a,b Aligned numeric program vectors.
  #' @return Pearson correlation or NA for a constant vector.
  #' @examples
  #' # safe_cor(1:3, 3:1)
  safe_cor <- function(a, b) if (sd(a) == 0 || sd(b) == 0) NA_real_ else cor(a, b)
  cross_cor <- do.call(rbind, lapply(seq_len(ncol(lb$EF)), function(i) do.call(rbind, lapply(seq_len(ncol(yf$EF)), function(j)
    data.frame(LB_program = i, YFB_program = j, cor_LB_L_YFB_YF = safe_cor(lb$EL[, i], yfb_score[, j]),
      cor_EF = safe_cor(lb$EF[, i], yf$EF[, j]))))))
  tables <- list(program_coefficients = programs, program_activity = activity, top_genes = top_table,
    top_gene_counts = top_counts, cross_model_top20_overlap = across, desurv_top20_overlap = comparator,
    cross_model_program_correlations = cross_cor)
  for (name in names(tables)) write.csv(tables[[name]], file.path(traindir, paste0(name, '.csv')), row.names = FALSE)
  invisible(c(list(simulation_summary = simulation), tables))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) > 1L || (length(args) && !grepl('^--date=[0-9]{4}-[0-9]{2}-[0-9]{2}$', args))) stop('Use --date=YYYY-MM-DD')
  lb_yfb_summary_main(if (length(args)) sub('^--date=', '', args) else '2026-11-12')
}
