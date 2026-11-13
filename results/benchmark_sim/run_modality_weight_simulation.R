# ============================================================
# Script: run_modality_weight_simulation.R
# Purpose: Run registered balance and joint/exclusive modality-weight simulations.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, ebnm, survival and existing feature-group simulation helpers
# ============================================================
# Sequential worker, one numerical thread. Synthetic fits stay under ignored
# data/. Aggregate checkpoints retain failed fits and every warning. No tuning.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_feature_group_simulation.R")

#' Compute weighted modality precision, final emptiness and repeated reset counts
#' @param fit Fitter result.
#' @param groups Modality labels in feature order.
#' @return Long program-by-modality data frame. Undefined zero-total shares stay NA.
#' @examples
#' # mwsim_precision(fit, labels)
mwsim_precision <- function(fit, groups) {
  weights <- if (is.null(fit$modality_weight)) rep(1, length(groups)) else
    unname(fit$modality_weight[groups])
  do.call(rbind, lapply(seq_len(ncol(fit$EF)), function(k) {
    precision <- vapply(unique(groups), function(m) {
      j <- groups == m
      sum(weights[j] * fit$Tau[j] * fit$EF2[j, k])
    }, numeric(1))
    data.frame(program = k, modality = names(precision), precision = unname(precision),
      share = if (sum(precision) > 0) unname(precision) / sum(precision) else NA_real_,
      empty = vapply(unique(groups), function(m) all(fit$EF2[groups == m, k] == 0), logical(1)),
      point_mass_resets = if (is.null(fit$empty_blocks)) NA_integer_ else fit$empty_blocks[, k])
  }))
}

#' Score independently matched loading directions and all three modality patterns
#' @param fit Returned posterior moments.
#' @param g Feature-by-three standardized, unit-norm reconstruction truth.
#' @param groups Modality labels in feature order.
#' @param z_threshold Posterior coefficient activity threshold.
#' @return One row per truth, with loading/pattern and separate survival recovery.
#' @examples
#' # mwsim_patterns(fit, truth, groups, 1.96)
mwsim_patterns <- function(fit, g, groups, z_threshold) {
  precision <- mwsim_precision(fit, groups)
  do.call(rbind, lapply(1:3, function(t) {
    r <- vapply(seq_len(ncol(fit$EF)), function(k)
      if (sd(fit$EF[, k]) > 0) cor(g[, t], fit$EF[, k]) else NA_real_, numeric(1))
    k <- if (any(is.finite(r))) which.max(replace(abs(r), !is.finite(r), -Inf)) else NA_integer_
    matched <- !is.na(k)
    p <- precision[precision$program == k & !is.na(k), ]
    shares <- stats::setNames(p$share, p$modality)
    empty <- stats::setNames(p$empty, p$modality)
    pattern <- if (!matched) FALSE else if (t < 3L) {
      absent <- if (t == 1L) "methylation" else "expression"
      isTRUE(empty[[absent]]) || isTRUE(shares[[absent]] < 0.05)
    } else all(is.finite(shares)) && all(shares >= 0.05)
    recovered <- matched && abs(r[k]) >= 0.7
    variance <- if (is.null(fit$VBeta)) pmax(fit$EBeta2 - fit$EBeta^2, 0) else fit$VBeta
    active <- matched && variance[k] > 0 && abs(fit$EBeta[k]) >= z_threshold * sqrt(variance[k])
    data.frame(truth_program = t, matched_program = k,
      loading_correlation = if (matched) r[k] else NA_real_,
      loading_recovered = recovered, pattern_correct = pattern,
      expression_share = if (matched) shares[["expression"]] else NA_real_,
      methylation_share = if (matched) shares[["methylation"]] else NA_real_,
      active = active, correct_sign = matched && fit$EBeta[k] * r[k] > 0,
      survival_recovered = t < 3L && recovered && active && fit$EBeta[k] * r[k] > 0)
  }))
}

#' Apply the preregistered gate and tie breakers without external outcomes
#' @param A Block A cell/arm summary including gate and primary recovery counts.
#' @param B Block B summary including full pattern success counts.
#' @return Single-row data frame with effectiveness and selected rule.
#' @examples
#' # mwsim_select(summary_A, summary_B)
mwsim_select <- function(A, B) {
  effective <- stats::setNames(logical(2), c("feature_count", "eigenvalue"))
  arms <- c(feature_count = "grouped_count", eigenvalue = "grouped_eigenvalue")
  for (rule in names(arms)) {
    x <- A[A$arm == arms[[rule]], ]
    baseline <- A[A$arm == "grouped_unweighted", ]
    key <- paste(x$feature_ratio, x$methylation_active_fraction)
    baseline <- baseline[match(key, paste(baseline$feature_ratio, baseline$methylation_active_fraction)), ]
    control <- x$feature_ratio == 1 & x$methylation_active_fraction == .03
    effective[rule] <- sum(control) == 1 && isTRUE(x$passed[control]) &&
      any(x$feature_ratio == 5 & x$passed & !baseline$passed)
  }
  candidates <- names(effective)[effective]
  selected <- "none"
  if (length(candidates) == 1L) selected <- candidates
  if (length(candidates) == 2L) {
    successes <- vapply(arms, function(arm) sum(B$full_pattern_successes[B$arm == arm]), numeric(1))
    recovery <- vapply(arms, function(arm) sum(A$recovered[A$arm == arm & A$feature_ratio == 5]), numeric(1))
    selected <- names(arms)[order(-successes, -recovery, seq_along(arms))[1L]]
  }
  data.frame(count_effective = effective[[1]], eigenvalue_effective = effective[[2]],
    selected_rule = selected, part4_allowed = selected != "none")
}

#' Fit one registered arm from fixed data and common initial moments
#' @param x Prepared input with Y, outcomes and common starts.
#' @param weights NULL or named reconstruction powers.
#' @param pooled Whether to pool F priors across modalities.
#' @param cfg Global configuration.
#' @return Fitter result.
#' @examples
#' # mwsim_fit(x, rules$feature_count, FALSE, cfg)
mwsim_fit <- function(x, weights, pooled, cfg) {
  fit_ssbmf_yfb(x$Y, x$time, x$event, K = cfg$benchmark$k_working_yfb,
    max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol, alpha = cfg$benchmark$alpha,
    prior_LF = "point_laplace", prior_beta = "normal", canonicalize = TRUE,
    w_S = 1, feature_group = x$groups, modality_weight = weights,
    F_prior_by_group = !pooled, EL_init = x$EL_init, EF_init = x$EF_init, verbose = FALSE)
}

#' Run and checkpoint the complete registered study or a one-seed smoke study
#' @param smoke Use seed 1 only. Its aggregates cannot select a weight rule.
#' @return Invisible A/B summaries and selection, also written to new directories.
#' @examples
#' # mwsim_run(TRUE)
mwsim_run <- function(smoke = FALSE) {
  suffix <- if (smoke) "2026-11-12_weights_smoke" else "2026-11-12_weights"
  out <- paste0("results/benchmark_sim/outputs/feature_group_simulation_", suffix)
  private <- paste0("data/feature_group_simulation_", suffix)
  if (dir.exists(out) || dir.exists(private)) stop("Refusing to overwrite a registered run.")
  dir.create(out, recursive = TRUE); dir.create(private, recursive = TRUE)
  cfg <- yaml::read_yaml("config/globals.yml")
  if (cfg$benchmark$k_working_yfb != 5L) stop("The registered working rank is five.")
  seeds <- if (smoke) 1L else 1:10
  # Freeze every rule before the first scientific fit.
  writeLines(c("registration=docs/prompts/2026-11-12-codex-d-modality-weights-prompt.md",
    paste0("seeds=", paste(seeds, collapse = ",")), "n_train=120; n_validation=120; K=5",
    "A_cells=50:250/15%,50:250/3%,50:50/3%",
    "A_arms=grouped_unweighted,grouped_count,grouped_eigenvalue,expression_only",
    "B_cells=50:50/15%,50:250/15%",
    "B_arms=grouped_count,grouped_eigenvalue,pooled_count,grouped_unweighted",
    "generator=both_informative; existing defaults; exclusive option only for B",
    "B_truth=expression_only1,methylation_only2,joint3; beta=0.70,0.35,0",
    "scaling=per-feature training-only z-score; weights from training rows only",
    "survival=fgsim_prepare programs1:2; seed+12000; preserve original training eta SD",
    "initialization=common joint signed SVD; expression uses corresponding F rows",
    "fit=point_laplace L/F; per-program Normal beta; canonicalize TRUE; wS1; no pruning",
    paste0("max_iter=", cfg$cavi$max_iter, "; tol=", cfg$cavi$tol, "; alpha=", cfg$benchmark$alpha),
    "A_recovery=program1 loading abs-cor>=0.7, posterior activity, correct beta sign",
    "A_false_active=active outside unique abs-cor>=0.7 matches to truths1,2, regardless beta sign",
    "gate=recovery>=8/10; mean false-active<=1; convergence>=9/10; no errors",
    "effective=passes a 50:250 cell unweighted fails AND passes 50:50/3% control",
    "selection=only effective rule; else B grouped full successes, then A 50:250 primary recoveries, then count",
    "B_matching=independent largest abs loading correlation per truth",
    "B_pattern=absent modality empty or precision share<0.05; joint both>=0.05",
    "B_full_success=all three loading recoveries and correct modality patterns; survival separately reported",
    "empty=final all-zero second moments per block; reset events separately, NA when pooled not instrumented",
    "warnings=retain each event; failed fits retained; zero-total precision shares NA",
    "workers=1; numerical_threads=1"), file.path(out, "settings.txt"))
  sources <- c("code/fit_ssbmf_yfb.R", "code/modality_weights_yfb.R", "code/feature_groups_yfb.R",
    "code/simulate_multimodal_yfb.R", "code/multimodal_yfb_updates.R", "config/globals.yml",
    "results/benchmark_sim/run_feature_group_simulation.R",
    "results/benchmark_sim/run_modality_weight_simulation.R")
  write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
    file.path(out, "source_hashes.csv"), row.names = FALSE)
  tables <- list(metrics = NULL, precision = NULL, patterns = NULL, warnings = NULL)
  cells <- rbind(data.frame(block = "A", ratio = c(5L, 5L, 1L), fraction = c(.15, .03, .03)),
                 data.frame(block = "B", ratio = c(1L, 5L), fraction = .15))
  for (cell in seq_len(nrow(cells))) for (seed in seeds) {
    row <- cells[cell, ]; exclusive <- row$block == "B"
    data <- simulate_multimodal_yfb_data("both_informative", p_methylation = 50L * row$ratio,
      seed = seed, methylation_active_fraction = row$fraction,
      modality_exclusive_programs = exclusive)
    x <- fgsim_prepare(data, seed, programs = 1:2)
    arms <- if (exclusive) c("grouped_count", "grouped_eigenvalue", "pooled_count", "grouped_unweighted") else
      c("grouped_unweighted", "grouped_count", "grouped_eigenvalue", "expression_only")
    inputs <- fgsim_arm_inputs(x, cfg$benchmark$k_working_yfb)
    rules <- yf_modality_weight_rules(x$Y, x$groups)
    g3 <- unlist(Map(function(F, s) F[, 3] / s, data$truth$F, x$scaled$sd), use.names = FALSE)
    g3 <- g3 / sqrt(sum(g3^2))
    for (arm in arms) {
      arm_x <- if (arm == "expression_only") inputs$expression_only else inputs$grouped
      weights <- if (grepl("_count$", arm)) rules$feature_count else
        if (arm == "grouped_eigenvalue") rules$eigenvalue else NULL
      base <- data.frame(block = row$block, feature_ratio = row$ratio,
        methylation_active_fraction = row$fraction, seed = seed, arm = arm)
      start <- proc.time()[["elapsed"]]; warnings <- character()
      answer <- tryCatch(withCallingHandlers({
        fit <- mwsim_fit(arm_x, weights, arm == "pooled_count", cfg)
        recovery <- fgsim_recovery_two(fit, arm_x$g, arm_x$groups, cfg$multimodal_yfb$prognostic_z_threshold)
        precision <- mwsim_precision(fit, arm_x$groups)
        precision$weight <- if (is.null(weights)) 1 else unname(weights[precision$modality])
        pattern <- if (exclusive) mwsim_patterns(fit, cbind(x$g, g3), x$groups,
          cfg$multimodal_yfb$prognostic_z_threshold) else NULL
        score <- predict_ssbmf_yfb(arm_x$Y_validation, fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
        saveRDS(list(fit = fit, inputs = arm_x, truth = data$truth), file.path(private,
          sprintf("%s_ratio%d_fraction%s_seed%d_%s.rds", row$block, row$ratio, row$fraction, seed, arm)))
        list(metrics = data.frame(error = "", recovered = recovery$recovered,
          recovered2 = recovery$recovered2, false_active = recovery$false_active,
          converged = isTRUE(fit$history$converged), n_iter = fit$history$n_iter,
          empty_blocks = sum(precision$empty),
          point_mass_resets = if (is.null(fit$empty_blocks)) NA_integer_ else sum(fit$empty_blocks),
          full_pattern_success = if (exclusive) all(pattern$loading_recovered & pattern$pattern_correct) else NA,
          c_validation = frozen_reverse_cindex(score, arm_x$time_validation, arm_x$event_validation)),
          precision = precision, patterns = pattern)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) list(metrics = data.frame(error = conditionMessage(e),
        recovered = FALSE, recovered2 = FALSE, false_active = NA_real_, converged = FALSE,
        n_iter = NA_integer_, empty_blocks = NA_integer_, point_mass_resets = NA_integer_,
        full_pattern_success = if (exclusive) FALSE else NA, c_validation = NA_real_)))
      answer$metrics$elapsed_seconds <- proc.time()[["elapsed"]] - start
      answer$metrics$n_warnings <- length(warnings)
      answer$warnings <- if (length(warnings)) data.frame(warning_number = seq_along(warnings), warning = warnings) else NULL
      for (type in names(tables)) if (!is.null(answer[[type]])) {
        rows <- cbind(base[rep(1, nrow(answer[[type]])), ], answer[[type]])
        tables[[type]] <- rbind(tables[[type]], rows)
        write.csv(tables[[type]], file.path(out, paste0(type, ".csv")), row.names = FALSE)
      }
      cat(sprintf("%s ratio%d fraction%.2f seed%d %s recovered=%s converged=%s error=%s\n",
        row$block, row$ratio, row$fraction, seed, arm, answer$metrics$recovered,
        answer$metrics$converged, answer$metrics$error)); flush.console()
    }
  }
  A <- fgsim_gate(tables$metrics[tables$metrics$block == "A", ])
  bmetrics <- tables$metrics[tables$metrics$block == "B", ]
  B <- do.call(rbind, lapply(split(bmetrics, interaction(bmetrics$feature_ratio, bmetrics$arm)), function(x)
    data.frame(feature_ratio = x$feature_ratio[1], arm = x$arm[1], seeds = nrow(x),
      failed = sum(nzchar(x$error)), full_pattern_successes = sum(x$full_pattern_success),
      converged = sum(x$converged), empty_blocks = sum(x$empty_blocks), mean_false_active = mean(x$false_active))))
  selection <- if (smoke) data.frame(count_effective = FALSE, eigenvalue_effective = FALSE,
    selected_rule = "smoke_only", part4_allowed = FALSE) else mwsim_select(A, B)
  write.csv(A, file.path(out, "summary_A.csv"), row.names = FALSE)
  write.csv(B, file.path(out, "summary_B.csv"), row.names = FALSE)
  write.csv(selection, file.path(out, "selection.csv"), row.names = FALSE)
  print(A); print(B); print(selection)
  invisible(list(A = A, B = B, selection = selection))
}
#' Correlate nonconstant fitted and truth loading blocks
#' @param truth,fitted Equal-length loading vectors, normalized implicitly by cor().
#' @return Signed Pearson correlation, NA for an absent/constant block.
#' @examples
#' mwpath_cor(1:4, 2*(1:4))
mwpath_cor <- function(truth, fitted) {
  if (sd(truth) == 0 || sd(fitted) == 0) return(NA_real_)
  cor(truth, fitted)
}

#' Apply all four registered usability conditions and fixed tie break
#' @param A Two-cell summaries with fraction, weight and primary retention counts.
#' @param B Per-truth loading-plus-pattern summaries.
#' @return One row per candidate, with component checks and selected_weight.
#' @examples
#' # mwpath_select(A, B)
mwpath_select <- function(A, B) {
  candidates <- sort(unique(A$weight[A$weight < 1]), decreasing = TRUE)
  base <- A[A$weight == 1, ]; baseB <- B[B$weight == 1, ]
  rows <- lapply(candidates, function(w) {
    a <- A[A$weight == w, ]; a <- a[match(base$fraction, a$fraction), ]
    b <- B[B$weight == w, ]; b <- b[match(baseB$truth_program, b$truth_program), ]
    data.frame(weight = w,
      no_loss_A = nrow(a)==2 && all(a$recovered >= base$recovered & a$mean_false_active <= 1),
      methylation_stays = all(a$methylation_retained >= 8),
      no_loss_B = all(b$loading_and_pattern[baseB$truth_program < 3] >=
        baseB$loading_and_pattern[baseB$truth_program < 3]),
      gain = sum(a$recovered)>sum(base$recovered) ||
        sum(b$loading_and_pattern)>sum(baseB$loading_and_pattern),
      B_total = sum(b$loading_and_pattern))
  })
  result <- do.call(rbind, rows)
  result$usable <- apply(result[c('no_loss_A','methylation_stays','no_loss_B','gain')], 1, all)
  eligible <- result[result$usable, ]
  selected <- if(nrow(eligible)) eligible$weight[order(-eligible$B_total,-eligible$weight)[1]] else NA_real_
  result$selected_weight <- selected
  result
}

#' Require exactly equal repeated-arm outcomes, excluding wall-clock timing
#' @param metrics New path metrics, complete or one fit.
#' @param old Earlier registered metrics.
#' @return TRUE or an error identifying the mismatching repeated fit.
#' @examples
#' # mwpath_reproduce(metrics, old_metrics)
mwpath_reproduce <- function(metrics, old) {
  fields <- c('recovered','recovered2','false_active','converged','n_iter','empty_blocks',
    'point_mass_resets','full_pattern_success','c_validation','n_warnings','error')
  repeats <- metrics[metrics$weight %in% c(1,.2), ]
  for(i in seq_len(nrow(repeats))) {
    x <- repeats[i, ]; arm <- if(x$weight==1) 'grouped_unweighted' else 'grouped_count'
    y <- old[old$block==x$block & old$feature_ratio==5 &
      old$methylation_active_fraction==x$methylation_active_fraction & old$seed==x$seed & old$arm==arm, ]
    if(nrow(y)!=1 || !isTRUE(all.equal(x[fields], y[fields], tolerance=0, check.attributes=FALSE)))
      stop('Reproducibility check failed: ',x$block,' fraction ',x$methylation_active_fraction,
           ' seed ',x$seed,' weight ',x$weight)
  }
  TRUE
}

#' Run the approved fixed weight path, checking repeated arms before interpretation
#' @param smoke Run one seed per cell and weight in separate directories.
#' @return Invisible aggregate tables, also checkpointed after every fit.
#' @examples
#' # mwpath_run(TRUE)
mwpath_run <- function(smoke=FALSE) {
  suffix <- paste0('2026-11-12_weight_path',if(smoke) '_smoke' else '')
  out <- paste0('results/benchmark_sim/outputs/feature_group_simulation_',suffix)
  private <- paste0('data/feature_group_simulation_',suffix)
  if(dir.exists(out)||dir.exists(private)) stop('Refusing to overwrite a path run.')
  dir.create(out,recursive=TRUE); dir.create(private,recursive=TRUE)
  cfg <- yaml::read_yaml('config/globals.yml')
  if(cfg$benchmark$k_working_yfb!=5) stop('Registered K must be five.')
  weights <- c(1,.7,1/sqrt(5),.3,.2); seeds <- if(smoke) 1L else 1:10
  writeLines(c('registration=docs/prompts/2026-11-12-codex-d-weight-path-prompt.md',
    paste0('weights=',paste(format(weights,digits=17),collapse=',')),
    paste0('seeds=',paste(seeds,collapse=',')),
    'cells=A ratio5 fraction.15,A ratio5 fraction.03,B ratio5 fraction.15',
    'all generator, survival regeneration, training scaling and signed SVD starts unchanged from weights run',
    'grouped point-Laplace L/F, Normal beta, K5,wS1,no pruning; workers1,numerical_threads1',
    paste0('max_iter=',cfg$cavi$max_iter,';tol=',cfg$cavi$tol,';alpha=',cfg$benchmark$alpha),
    'matching=independent largest absolute full loading Pearson correlation per truth',
    'block_cor=signed Pearson within modality; NA for zero/constant fitted or truth block; report defined denominator',
    'precision=weighted Tau*EF2; retention>=.05; zero-total share NA counts as not retained',
    'usable=w<1; BOTH A recovery>=w1 and mean false<=1; BOTH A primary retention>=8/10;',
    'B exclusive1,2 loading+pattern>=w1; strict gain summed A recovery OR summed B loading+pattern',
    'selection=largest B loading+pattern sum, then larger weight; convergence reported separately',
    'reproducibility=exact per-fit prior metric outcomes, precision and pattern tables for w1 and w.2; elapsed excluded'),
    file.path(out,'settings.txt'))
  sources <- c('results/benchmark_sim/run_modality_weight_simulation.R',
    'results/benchmark_sim/run_feature_group_simulation.R','code/fit_ssbmf_yfb.R',
    'code/modality_weights_yfb.R','code/simulate_multimodal_yfb.R','config/globals.yml')
  write.csv(data.frame(path=sources,md5=unname(tools::md5sum(sources))),file.path(out,'source_hashes.csv'),row.names=FALSE)
  previous <- 'results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weights'
  old <- read.csv(file.path(previous,'metrics.csv'),colClasses=c(error='character'))
  oldp <- read.csv(file.path(previous,'precision.csv')); oldb <- read.csv(file.path(previous,'patterns.csv'))
  tables <- list(metrics=NULL,precision=NULL,matched=NULL,patterns=NULL,warnings=NULL)
  cells <- data.frame(block=c('A','A','B'),fraction=c(.15,.03,.15))
  for(cell in seq_len(nrow(cells))) for(seed in seeds) {
    row <- cells[cell, ]; exclusive <- row$block=='B'
    data <- simulate_multimodal_yfb_data('both_informative',p_methylation=250L,seed=seed,
      methylation_active_fraction=row$fraction,modality_exclusive_programs=exclusive)
    x <- fgsim_prepare(data,seed,programs=1:2)
    inputs <- fgsim_arm_inputs(x,5L)$grouped
    g3 <- unlist(Map(function(F,s) F[,3]/s,data$truth$F,x$scaled$sd),use.names=FALSE)
    g3 <- g3/sqrt(sum(g3^2)); truth <- if(exclusive) cbind(x$g,g3) else x$g
    for(w in weights) {
      base <- data.frame(block=row$block,feature_ratio=5,methylation_active_fraction=row$fraction,seed=seed,weight=w)
      start <- proc.time()[['elapsed']]; warnings <- character()
      answer <- tryCatch(withCallingHandlers({
        fit <- mwsim_fit(inputs,c(expression=1,methylation=w),FALSE,cfg)
        rec <- fgsim_recovery_two(fit,inputs$g,inputs$groups,cfg$multimodal_yfb$prognostic_z_threshold)
        precision <- mwsim_precision(fit,inputs$groups)
        matched <- do.call(rbind,lapply(seq_len(ncol(truth)),function(t) {
          r <- vapply(seq_len(ncol(fit$EF)),function(k) mwpath_cor(truth[,t],fit$EF[,k]),numeric(1))
          k <- if(any(is.finite(r))) which.max(replace(abs(r),!is.finite(r),-Inf)) else NA_integer_
          do.call(rbind,lapply(unique(inputs$groups),function(m) {
            j <- inputs$groups==m
            p <- precision[precision$program==k & precision$modality==m & !is.na(k), ]
            data.frame(truth_program=t,matched_program=k,modality=m,
              full_correlation=if(is.na(k)) NA_real_ else r[k],
              block_correlation=if(is.na(k)) NA_real_ else mwpath_cor(truth[j,t],fit$EF[j,k]),
              share=if(nrow(p)) p$share else NA_real_,empty=if(nrow(p)) p$empty else NA,
              retained=if(nrow(p)) isTRUE(p$share>=.05) else FALSE)
          }))
        }))
        pattern <- if(exclusive) mwsim_patterns(fit,truth,inputs$groups,cfg$multimodal_yfb$prognostic_z_threshold) else NULL
        score <- predict_ssbmf_yfb(inputs$Y_validation,fit$EF,fit$EBeta,fit$EF_norms)$risk_scores
        metrics <- data.frame(error='',recovered=rec$recovered,recovered2=rec$recovered2,
          false_active=rec$false_active,converged=isTRUE(fit$history$converged),n_iter=fit$history$n_iter,
          empty_blocks=sum(precision$empty),point_mass_resets=sum(fit$empty_blocks),
          full_pattern_success=if(exclusive) all(pattern$loading_recovered & pattern$pattern_correct) else NA,
          c_validation=frozen_reverse_cindex(score,inputs$time_validation,inputs$event_validation))
        saveRDS(list(fit=fit,inputs=inputs,truth=data$truth),file.path(private,
          sprintf('%s_fraction%s_seed%d_weight%.10f.rds',row$block,row$fraction,seed,w)))
        list(metrics=metrics,precision=precision,matched=matched,patterns=pattern)
      },warning=function(e) {warnings <<- c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),
      error=function(e) list(metrics=data.frame(error=conditionMessage(e),recovered=FALSE,recovered2=FALSE,
        false_active=NA_real_,converged=FALSE,n_iter=NA_integer_,empty_blocks=NA_integer_,
        point_mass_resets=NA_integer_,full_pattern_success=if(exclusive) FALSE else NA,c_validation=NA_real_)))
      answer$metrics$elapsed_seconds <- proc.time()[['elapsed']]-start
      answer$metrics$n_warnings <- length(warnings)
      if(length(warnings)) answer$warnings <- data.frame(warning_number=seq_along(warnings),warning=warnings)
      for(type in names(tables)) if(!is.null(answer[[type]])) {
        rows <- cbind(base[rep(1,nrow(answer[[type]])), ],answer[[type]])
        tables[[type]] <- rbind(tables[[type]],rows)
        write.csv(tables[[type]],file.path(out,paste0(type,'.csv')),row.names=FALSE)
      }
      # Compare repeats immediately; never interpret changed reproductions.
      # Compare serialized aggregates with serialized aggregates at identical CSV precision.
      current <- tail(read.csv(file.path(out,'metrics.csv'),colClasses=c(error='character')),1L)
      mwpath_reproduce(current,old)
      if(w %in% c(1,.2) && !nzchar(answer$metrics$error)) {
        arm <- if(w==1) 'grouped_unweighted' else 'grouped_count'
        for(type in c('precision',if(exclusive) 'patterns')) {
          prior <- if(type=='precision') oldp else oldb
          prior <- prior[prior$block==row$block & prior$feature_ratio==5 &
            prior$methylation_active_fraction==row$fraction & prior$seed==seed & prior$arm==arm, ]
          current_rows <- tail(read.csv(file.path(out,paste0(type,'.csv'))),nrow(answer[[type]]))
          fields <- names(answer[[type]])
          # Earlier precision includes an explicit weight column absent here.
          if(!isTRUE(all.equal(current_rows[fields],prior[fields],tolerance=0,check.attributes=FALSE)))
            stop('Reproducibility ',type,' check failed.')
        }
      }
      cat(sprintf('%s %.2f seed%d w%.10f converged=%s error=%s\n',row$block,row$fraction,seed,w,
        answer$metrics$converged,answer$metrics$error));flush.console()
    }
  }
  writeLines('PASS: all repeated metric, precision and pattern rows match exactly',file.path(out,'reproducibility.txt'))
  mwpath_summarize(out,smoke)
  invisible(tables)
}

#' Summarize path cells and truth-specific block recovery without refitting
#' @param out Completed path aggregate directory.
#' @param smoke Whether selection is prohibited for a smoke run.
#' @return Invisible list of summaries.
#' @examples
#' # mwpath_summarize(out)
mwpath_summarize <- function(out,smoke=FALSE) {
  read <- function(name) read.csv(file.path(out,paste0(name,'.csv')),colClasses=if(name=='metrics') c(error='character') else NA)
  m <- read('metrics'); p <- read('precision'); t <- read('matched'); b <- read('patterns')
  if(nrow(m)!=(if(smoke) 15L else 150L) || anyDuplicated(m[c('block','methylation_active_fraction','seed','weight')]))
    stop('Incomplete or duplicate path metrics.')
  group <- function(x,cols) split(x,interaction(x[cols],drop=TRUE,lex.order=TRUE))
  expected <- expand.grid(cell=1:3,seed=if(smoke) 1L else 1:10,
    weight=c(1,.7,1/sqrt(5),.3,.2))
  expected$block <- c('A','A','B')[expected$cell]
  expected$methylation_active_fraction <- c(.15,.03,.15)[expected$cell]
  #' Check exact table identities after common CSV precision rounding
  #' @param actual,wanted Data frames whose keys must agree exactly.
  #' @param cols Key column names.
  #' @return NULL, or an error for missing, extra or duplicate rows.
  #' @examples
  #' # check_rows(m, expected, c('block','seed','weight'))
  check_rows <- function(actual,wanted,cols) {
    key <- function(x) do.call(paste,lapply(x[cols],function(z)
      if(is.numeric(z)) format(z,digits=14,trim=TRUE) else z))
    a <- key(actual); w <- key(wanted)
    if(anyDuplicated(a) || length(a)!=length(w) || !setequal(a,w))
      stop('Aggregate rows do not match the registered path design.')
  }
  fitkeys <- c('block','methylation_active_fraction','seed','weight')
  check_rows(m,expected,fitkeys)
  successful <- m[!nzchar(m$error),fitkeys]
  wantp <- merge(successful,expand.grid(program=1:5,modality=c('expression','methylation')))
  check_rows(p,wantp,c(fitkeys,'program','modality'))
  wantt <- do.call(rbind,lapply(seq_len(nrow(successful)),function(i)
    merge(successful[i,,drop=FALSE],expand.grid(truth_program=if(successful$block[i]=='B') 1:3 else 1:2,
      modality=c('expression','methylation')))))
  check_rows(t,wantt,c(fitkeys,'truth_program','modality'))
  wantb <- merge(successful[successful$block=='B',],data.frame(truth_program=1:3))
  check_rows(b,wantb,c(fitkeys,'truth_program'))
  keys <- c('block','methylation_active_fraction','weight')
  summary <- do.call(rbind,lapply(group(m,keys),function(x) {
    primary <- t[t$block==x$block[1] & t$methylation_active_fraction==x$methylation_active_fraction[1] &
      t$weight==x$weight[1] & t$truth_program==1 & t$modality=='methylation', ]
    data.frame(x[1,keys],seeds=nrow(x),failed=sum(nzchar(x$error)),recovered=sum(x$recovered),
      mean_false_active=mean(x$false_active),converged=sum(x$converged),
      methylation_retained=sum(primary$retained),n_warnings=sum(x$n_warnings))
  }))
  matched <- do.call(rbind,lapply(group(t,c(keys,'truth_program','modality')),function(x) {
    defined <- is.finite(x$block_correlation); shares <- is.finite(x$share)
    data.frame(x[1,c(keys,'truth_program','modality')],seeds=nrow(x),
      mean_signed_correlation=if(any(defined)) mean(x$block_correlation[defined]) else NA_real_,
      mean_abs_correlation=if(any(defined)) mean(abs(x$block_correlation[defined])) else NA_real_,
      correlation_defined=sum(defined),mean_share=if(any(shares)) mean(x$share[shares]) else NA_real_,
      share_defined=sum(shares),retained=sum(x$retained),empty=sum(x$empty),
      empty_or_under5=sum(x$empty | (is.finite(x$share) & x$share<.05),na.rm=TRUE))
  }))
  empty <- do.call(rbind,lapply(group(p,c(keys,'modality')),function(x)
    data.frame(x[1,c(keys,'modality')],empty=sum(x$empty),denominator=nrow(x))))
  patterns <- do.call(rbind,lapply(group(b,c(keys,'truth_program')),function(x)
    data.frame(x[1,c(keys,'truth_program')],seeds=nrow(x),loading=sum(x$loading_recovered),
      pattern=sum(x$pattern_correct),loading_and_pattern=sum(x$loading_recovered & x$pattern_correct))))
  a <- summary[summary$block=='A', ]; a$fraction <- a$methylation_active_fraction
  selection <- if(smoke || any(summary$failed>0)) data.frame(selected_weight=NA_real_,status='smoke_or_failed') else mwpath_select(a,patterns)
  result <- list(summary=summary,summary_matched=matched,summary_empty=empty,summary_patterns=patterns,selection=selection)
  for(name in names(result)) write.csv(result[[name]],file.path(out,paste0(name,'.csv')),row.names=FALSE)
  print(summary);print(selection);invisible(result)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if('--weight-path' %in% args) mwpath_run('--smoke' %in% args) else mwsim_run('--smoke' %in% args)
}
