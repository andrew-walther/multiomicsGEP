# ============================================================
# Script: run_k_selection_panel.R
# Purpose: Apply the fixed protected-pruning adoption rule on calibrated beds
# Author: Andrew Walther
# Created: 2026-10-15; Round 2 design 2026-10-16
# Dependencies: survival, ebnm, yaml, parallel; existing benchmark fitters
# Run: Rscript results/benchmark_sim/run_k_selection_panel.R [cores]
# ============================================================
# Attempt Task1-3: identical science; handle absent pruning histories in reporting.
# Fixed BEFORE execution and committed in Round2_attempts_log.md.
# Only scientific change: protection now requires SD(beta)>0 as Andrew decided.
# Seeds20261015+1:20; noise {.25,1}; K_init{5,7,10,15}; n_train400,
# n_validation120,p50 per modality,beta3=3.3, unchanged calibrated beds.
# Single: selected + protected pruning; alpha_F=0 signed-L/F reference.
# MM: selected with/without protection; 10/1 per-factor-beta framework reference.
# All defaults unchanged. Convergence: single ELBO, MM its existing outputs;
# single final output criterion also recorded, not adopted. The fitter-returned
# convergence flag is the panel criterion; caps in intermediate pruning fits
# are separately reported as pruning_caps (not silently hidden).
# Adoption at EVERY K/noise: recovery>=16/20, mean false<=reference+.5,
# every fit converged under its own rule, median final K=3. Missing/error rows
# cannot pass. No tuning. Frozen risk orientation. No real/external data loaded.
# 800 independent cells, socket workers only (one serial fit per worker).
# Each completed cell, including errors/warnings, is written before aggregation.
# Refuse overwrite of an existing attempt. Never edit this or its fitters while running.

#' Load existing fitter and diagnostic definitions without executing their grids
#' @return Nothing; functions and constants loaded in the global environment.
#' @examples
#' load_panel_code()
load_panel_code <- function() {
  # Evaluate definitions only; the original runner's CLI/grid is not sourced.
  ex <- parse("results/benchmark_sim/run_supervision_diagnostics.R")
  wanted <- c("load_code", "make_testbed", "score_fit")
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                    as.character(e[[2]]) %in% wanted) eval(e, .GlobalEnv)
  load_code()
  if (!exists("fit_cox_on_yf")) stop("Single-modality fitter was not loaded")
  source("code/multimodal_yfb_helpers.R", local = .GlobalEnv)
  source("code/preprocess_multimodal_yfb.R", local = .GlobalEnv)
  source("code/fit_multimodal_yfb.R", local = .GlobalEnv)
  source("code/predict_multimodal_yfb.R", local = .GlobalEnv)
  assign("N_TRAIN", 400L, .GlobalEnv); assign("BETA3", 3.3, .GlobalEnv)
}

#' Fit and score a fixed panel cell, retaining a failure as a row
#' @param i Row number in the fixed grid.
#' @return One data-frame row of metrics or an explicit error.
#' @examples
#' # panel_cell(1L)
panel_cell <- function(i) {
  x <- grid[i, ]; t0 <- proc.time()[["elapsed"]]
  id <- x
  fit_warnings <- character()
  answer <- withCallingHandlers(tryCatch({
    if (x$model == "single") {
      tb <- make_testbed(x$seed, x$noise)
      a <- list(K = x$K, max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                prior_LF = "point_laplace", verbose = FALSE, trace_supervision = TRUE, converge_on = "elbo")
      a <- modifyList(a, if (x$arm == "protected")
        list(canonicalize = TRUE, w_S = 1, init_supervised = TRUE,
             beta_prior_mode = "shared", prior_beta = "point_laplace", beta_warmup = 10L,
             prune = TRUE, prune_protect_survival = TRUE) else list(alpha_F = 0))
      fit <- do.call(fit_cox_on_yf, c(list(tb$Ytr, tb$ttr, tb$str), a))
      sc <- score_fit(fit, tb)
      metrics <- data.frame(recovered = sc$recovered, false_active = sc$false_z,
        K_final = ncol(fit$EF), converged = fit$history$converged,
        elbo_converged = fit$history$elbo_converged,
        outputs_converged = fit$history$outputs_converged,
        c_validation = sc$c_validation, cor_planted = sc$cor_planted,
        n_active = sc$n_active_z, iterations = sc$n_iter,
        pruning_caps = sum(!as.logical(fit$pruning$converged_before)),
        runaway_F = sc$runaway_F, runaway_L = sc$runaway_L,
        max_abs_beta = sc$max_abs_beta)
    } else {
      d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train = 400L,
        n_validation = 120L, p_expression = 50L, p_methylation = 50L,
        seed = x$seed, truncate = FALSE, noise_scale = x$noise, beta = c(0, 0, 3.3))
      ctrl <- list(intercept = TRUE, prior_update = "ebnm", prune = TRUE,
        prior_F = list(expression = "point_laplace", methylation = "point_laplace"),
        survival_elbo = "plugin", init_supervised = TRUE,
        beta_prior = "shared_point_laplace", beta_warmup = 10L,
        prune_protect_survival = x$arm == "protected")
      if(x$arm=="reference") {
        ctrl$init_supervised <- NULL; ctrl$beta_prior <- NULL; ctrl$beta_warmup <- NULL
      }
      fit <- fit_multimodal_yfb(d$training$Y, d$training$time, d$training$event, x$K, ctrl)
      r <- suppressWarnings(cor(do.call(rbind, d$truth$F), do.call(rbind, fit$EF)))
      r[!is.finite(r)] <- 0; k <- which.max(abs(r[3, ])); rho <- r[3, k]
      active <- fit$diagnostics$survival_active
      h <- tail(fit$diagnostics$history, 1)[[1]]
      risk <- predict_multimodal_yfb(fit, d$validation$Y)$risk_scores
      metrics <- data.frame(recovered = abs(rho) >= .7 && sign(rho * fit$EBeta[k]) == 1,
        false_active = sum(active[setdiff(seq_along(active), if (abs(rho) >= .7) k else integer())]),
        K_final = ncol(fit$EL), converged = fit$diagnostics$converged,
        elbo_converged = NA,
        outputs_converged = max(h$max_relative_reconstruction_change, h$max_relative_eta_change) < 1e-4,
        c_validation = multimodal_yfb_survival_metrics(risk, d$validation$time, d$validation$event)$c_index,
        cor_planted = abs(rho), n_active = sum(active), iterations = fit$diagnostics$iterations,
        pruning_caps = sum(!as.logical(fit$diagnostics$pruning$converged_before)),
        runaway_F = NA_real_, runaway_L = NA_real_, max_abs_beta = max(abs(fit$EBeta)))
    }
    cbind(id, error = NA_character_, runtime = proc.time()[["elapsed"]] - t0, metrics)
  }, error = function(e) cbind(id, error = conditionMessage(e),
                              runtime = proc.time()[["elapsed"]] - t0)),
    warning=function(w) {fit_warnings <<- c(fit_warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  answer$warnings <- paste(unique(fit_warnings),collapse=" | ")
  write.csv(answer,file.path(out_dir,"cells",sprintf("cell_%04d.csv",i)),row.names=FALSE)
  answer
}

#' Construct every prescribed seed, rank, noise and arm combination
#' @return Data frame of 800 unique cells, 320 single and 480 multimodal.
#' @examples
#' nrow(k_panel_grid())
k_panel_grid <- function() {
  single <- expand.grid(seed=20261015L+1:20,noise=c(.25,1),K=c(5L,7L,10L,15L),
    model="single",arm=c("reference","protected"),stringsAsFactors=FALSE)
  multi <- expand.grid(seed=20261015L+1:20,noise=c(.25,1),K=c(5L,7L,10L,15L),
    model="multimodal",arm=c("reference","selected","protected"),stringsAsFactors=FALSE)
  rbind(single,multi)
}

#' Summarize the fixed adoption criteria without dropping failed cells
#' @param res Per-cell metric/error rows, including model, arm, K and noise.
#' @return Group table, reference comparisons and separate criterion flags.
#' @examples
#' # k_panel_summary(results)
k_panel_summary <- function(res) {
  groups <- split(res,list(res$model,res$arm,res$K,res$noise),drop=TRUE)
  summ <- do.call(rbind,lapply(groups,function(d) {
    ok<-d[is.na(d$error),,drop=FALSE]
    data.frame(model=d$model[1],arm=d$arm[1],K=d$K[1],noise=d$noise[1],n=nrow(d),
      errors=sum(!is.na(d$error)),complete=nrow(d)==20L && all(is.na(d$error)),
      recovered=sum(ok$recovered),false_active=mean(ok$false_active),
      median_K=median(ok$K_final),converged=sum(ok$converged),
      elbo_converged=if(all(is.na(ok$elbo_converged)))NA_integer_ else sum(ok$elbo_converged),
      outputs_converged=sum(ok$outputs_converged),C=mean(ok$c_validation))
  }))
  ref<-summ[summ$arm=="reference",c("model","K","noise","false_active","complete")]
  names(ref)[4:5]<-c("reference_false_active","reference_complete")
  summ<-merge(summ,ref,by=c("model","K","noise"),all.x=TRUE)
  summ$convergence_rule<-ifelse(summ$model=="single","elbo","outputs")
  summ$recovery_pass<-summ$complete & summ$recovered>=16L
  summ$false_positive_pass<-summ$complete & summ$reference_complete &
    summ$false_active<=summ$reference_false_active+.5
  summ$convergence_pass<-summ$complete & summ$converged==20L
  summ$rank_pass<-summ$complete & summ$median_K==3L
  summ
}

#' Run the pre-registered panel once, with independent error isolation
#' @param cores Positive socket-worker count; one serial cell per worker.
#' @param output_dir Fresh attempt folder; existing cell folders are never overwritten.
#' @return Invisibly the summary; writes cells, results, summary and adoption CSVs.
#' @examples
#' # run_k_selection_panel(6L)
run_k_selection_panel <- function(cores=6L,
  output_dir="results/benchmark_sim/outputs/k_selection/round2/attempt3") {
  out_dir <- output_dir
  if(dir.exists(file.path(out_dir,"cells")))stop("Attempt output already exists; do not overwrite")
  dir.create(file.path(out_dir,"cells"),recursive=TRUE)
  cfg <- yaml::read_yaml("config/globals.yml");grid <- k_panel_grid()
  write.csv(grid,file.path(out_dir,"design.csv"),row.names=FALSE)
  cl<-parallel::makeCluster(cores)
  parallel::clusterExport(cl,c("load_panel_code","panel_cell","grid","cfg","out_dir"),
    envir=environment())
  invisible(parallel::clusterEvalQ(cl,load_panel_code()))
  rows<-tryCatch(parallel::parLapplyLB(cl,seq_len(nrow(grid)),panel_cell),
    finally=parallel::stopCluster(cl))
  cols<-unique(unlist(lapply(rows,names)))
  res<-do.call(rbind,lapply(rows,function(r) {r[setdiff(cols,names(r))]<-NA;r[cols]}))
  write.csv(res,file.path(out_dir,"results.csv"),row.names=FALSE)
  summ<-k_panel_summary(res)
  write.csv(summ,file.path(out_dir,"summary.csv"),row.names=FALSE)
  flags<-c("recovery_pass","false_positive_pass","convergence_pass","rank_pass")
  adoption<-do.call(rbind,lapply(c("single","multimodal"),function(model) {
    x<-summ[summ$model==model & summ$arm=="protected",]
    pass<-nrow(x)==8L && all(!is.na(as.matrix(x[flags])) & as.matrix(x[flags]))
    failed<-vapply(flags,function(f)any(is.na(x[[f]]) | !x[[f]]),logical(1))
    data.frame(model=model,adopted=pass,reason=if(pass)"All simulation criteria pass; real-data criteria pending"
      else paste(names(failed)[failed],collapse=";"))
  }))
  write.csv(adoption,file.path(out_dir,"adoption.csv"),row.names=FALSE)
  print(summ,row.names=FALSE);print(adoption,row.names=FALSE)
  if(any(!is.na(res$error)))warning("Panel errors retained in results.csv; no incomplete group passes")
  invisible(summ)
}
if(sys.nframe()==0L) {
  a<-commandArgs(TRUE);run_k_selection_panel(if(length(a))as.integer(a[1])else 6L)
}
