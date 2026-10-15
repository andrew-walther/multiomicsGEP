# ============================================================
# Script: run_k_selection_panel.R
# Purpose: Apply the fixed protected-pruning adoption rule on calibrated beds
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, yaml, parallel; existing benchmark fitters
# Run: Rscript results/benchmark_sim/run_k_selection_panel.R [cores]
# ============================================================
# Fixed before execution: seeds 20261015 + 1:20; noise {0.25,1}; K {5,7,10,15}.
# Single: selected canonical shared point-Laplace model, protected pruning;
# reference alpha_F=0, signed point-Laplace L/F, without pruning.
# Multimodal: M4 selected model with and without protection (reference).
# Beds: n_train=400, n_validation=120, p=50+50, beta3=3.3.
# Metrics: loading recovery >=.7 with correct sign; false survival-active,
# validation C, retained K, convergence and final output-change criterion.
# Adoption at EVERY K/noise: recovery >=16/20, false <= reference+.5,
# all ELBO-converged, median K_final=3. Never tune after results.
# The multimodal fitter stops on outputs, not ELBO: report ELBO as NA;
# an unavailable mandatory convergence criterion cannot support adoption.
# Socket workers only; each worker is serial. Outputs include failed cells.

suppressPackageStartupMessages(library(parallel))
out_dir <- "results/benchmark_sim/outputs/k_selection"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE)
cores <- if (length(args)) as.integer(args[1]) else 6L

#' Load existing fitter and diagnostic definitions without executing their grids
#' @return Nothing; functions and their constants are loaded in the global environment.
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
  tryCatch({
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
        runaway_F = NA_real_, runaway_L = NA_real_, max_abs_beta = max(abs(fit$EBeta)))
    }
    cbind(id, error = NA_character_, runtime = proc.time()[["elapsed"]] - t0, metrics)
  }, error = function(e) cbind(id, error = conditionMessage(e),
                              runtime = proc.time()[["elapsed"]] - t0))
}

cfg <- yaml::read_yaml("config/globals.yml")
grid <- expand.grid(seed = 20261015L + 1:20, noise = c(.25, 1), K = c(5L,7L,10L,15L),
                    model = c("single", "multimodal"), arm = c("reference", "protected"),
                    stringsAsFactors = FALSE)
cl <- makeCluster(cores)
clusterExport(cl, c("load_panel_code", "panel_cell", "grid", "cfg"))
invisible(clusterEvalQ(cl, load_panel_code()))
# One batch per K, checkpointed before the next. All cells and arms remain fixed.
rows <- list()
for (k in unique(grid$K)) {
  idx <- which(grid$K == k)
  message("Starting K = ", k, " (", length(idx), " cells)")
  batch <- parLapplyLB(cl, idx, panel_cell)
  rows <- c(rows, batch)
  cols <- unique(unlist(lapply(rows, names)))
  res <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols,names(r))] <- NA; r[cols] }))
  write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)
}
stopCluster(cl)
summ <- do.call(rbind, lapply(split(res, list(res$model,res$arm,res$K,res$noise), drop=TRUE), function(d) {
  ok <- d[is.na(d$error), ]
  data.frame(model=d$model[1],arm=d$arm[1],K=d$K[1],noise=d$noise[1],n=nrow(d),
    errors=sum(!is.na(d$error)), recovered=sum(ok$recovered),false_active=mean(ok$false_active),
    median_K=median(ok$K_final),converged=sum(ok$converged),
    elbo_converged=if(all(is.na(ok$elbo_converged))) NA_integer_ else sum(ok$elbo_converged),
    outputs_converged=sum(ok$outputs_converged),C=mean(ok$c_validation))
}))
ref <- summ[summ$arm=="reference",c("model","K","noise","false_active")]
names(ref)[4] <- "reference_false_active"
summ <- merge(summ,ref,by=c("model","K","noise"))
summ$recovery_pass <- summ$recovered >= 16
summ$false_positive_pass <- summ$false_active <= summ$reference_false_active+.5
summ$convergence_pass <- !is.na(summ$elbo_converged) & summ$elbo_converged==20 & summ$errors==0
summ$rank_pass <- summ$median_K==3
write.csv(summ,file.path(out_dir,"summary.csv"),row.names=FALSE)
print(summ,row.names=FALSE)
if(any(!is.na(res$error))) warning("Panel contains errors; see results.csv; no adoption permitted.")
