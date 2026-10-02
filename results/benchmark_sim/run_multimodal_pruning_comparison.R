# ============================================================
# Script: run_multimodal_pruning_comparison.R
# Purpose: Compare the original multimodal YFB fit with the parsimony
#          framework (intercept + point-Laplace loadings + ebnm priors + ELBO
#          factor pruning) on simulations with known K (true K = 3).
#          Question: does the final factor count stop tracking K_init, and are
#          the variance-explaining and survival-explaining programs recovered?
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: survival, ebnm, yaml, parallel; code/fit_multimodal_yfb.R,
#               code/simulate_multimodal_yfb.R
# Run: Rscript results/benchmark_sim/run_multimodal_pruning_comparison.R [n_seeds] [cores]
# ============================================================
#
# Data: every variant sees the same data. Gaussian data are generated without
# truncation and shifted by a per-feature baseline of 4 noise SDs, so they
# are (almost surely) nonnegative and the original point-exponential model can
# be fit; a model with an intercept is unaffected by the shift.
# noise_scale = 1 is the September design (each true factor explains < 1% of
# variance); 0.25 is a moderate signal-to-noise ratio (~25%).
#
# Variants (control passed to fit_multimodal_yfb()):
#   original   EM prior updates, point-exponential F, no intercept, no pruning
#   noprune    intercept + point-Laplace F + ebnm priors, no pruning (ablation)
#   pruned     noprune + ELBO nullcheck, survival = partial log-lik at E[eta]
#   pruned_vc  pruned with the variance-corrected survival term
#
# Parallelization: one socket-cluster task per simulated dataset (parLapplyLB);
# each fit is single-threaded and small (n = 120, p = 50 + 50).

suppressPackageStartupMessages(library(survival))
for (f in c("code/multimodal_yfb_helpers.R", "code/preprocess_multimodal_yfb.R",
            "code/multimodal_yfb_updates.R", "code/fit_multimodal_yfb.R",
            "code/predict_multimodal_yfb.R", "code/simulate_multimodal_yfb.R")) source(f)

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 3L
cores <- if (length(args) >= 2) as.integer(args[2]) else 6L
out_dir <- "results/benchmark_sim/outputs/multimodal_yfb_pruning"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

laplace <- list(expression = "point_laplace", methylation = "point_laplace")
variants <- list(
  original  = list(),
  noprune   = list(intercept = TRUE, prior_update = "ebnm", prior_F = laplace),
  pruned    = list(intercept = TRUE, prior_update = "ebnm", prior_F = laplace,
                   prune = TRUE, survival_elbo = "plugin"),
  pruned_vc = list(intercept = TRUE, prior_update = "ebnm", prior_F = laplace,
                   prune = TRUE, survival_elbo = "corrected")
)
grid <- expand.grid(
  scenario = c("both_informative", "adverse_protective", "low_variance_prognostic_factor"),
  noise_scale = c(0.25, 1), K_init = c(5L, 8L, 12L, 16L),
  seed = 20261001L + seq_len(n_seeds), stringsAsFactors = FALSE)

#' Add a 4-noise-SD baseline to each feature (and truncate at 0, which
#' almost never binds)
with_baseline <- function(Y, truth, noise_scale) {
  Map(function(x, tau) pmax(sweep(x, 2, 4 * noise_scale / sqrt(tau), "+"), 0), Y, truth$Tau)
}

#' Factor-level metrics against the truth
#' Each true factor is matched to the estimated factor with the largest
#' |cor| of stacked loadings. A survival-active estimated factor is a false
#' positive if it is not the match of any true prognostic factor.
score_fit <- function(fit, data, Yva) {
  true_F <- do.call(rbind, data$truth$F)
  est_F <- do.call(rbind, fit$EF)
  r <- abs(suppressWarnings(stats::cor(true_F, est_F)))
  r[!is.finite(r)] <- 0
  best <- apply(r, 1, which.max)
  best_cor <- apply(r, 1, max)
  prognostic <- data$truth$beta != 0
  active <- fit$diagnostics$survival_active
  hit <- prognostic & best_cor >= 0.7 & active[best]
  direction <- hit & sign(fit$EBeta[best]) == sign(data$truth$beta)
  true_prognostic_matches <- unique(best[prognostic & best_cor >= 0.7])
  risk <- predict_multimodal_yfb(fit, Yva)$risk_scores
  data.frame(
    K_final = ncol(fit$EL),
    K_survival = sum(active),
    false_survival = sum(active[setdiff(seq_along(active), true_prognostic_matches)]),
    recovery_f1 = best_cor[1], recovery_f2 = best_cor[2], recovery_f3 = best_cor[3],
    factors_recovered = sum(best_cor >= 0.7),
    prognostic_recall = mean(hit[prognostic]),
    direction_recall = mean(direction[prognostic]),
    c_validation = multimodal_yfb_survival_metrics(risk, data$validation$time,
                                                   data$validation$event)$c_index,
    converged = fit$diagnostics$converged,
    n_pruned = NROW(fit$diagnostics$pruning))
}

run_one <- function(g) {
  data <- simulate_multimodal_yfb_data(grid$scenario[g], n_train = 120L, n_validation = 120L,
                                        p_expression = 50L, p_methylation = 50L,
                                        seed = grid$seed[g], truncate = FALSE,
                                        noise_scale = grid$noise_scale[g])
  Ytr <- with_baseline(data$training$Y, data$truth, grid$noise_scale[g])
  Yva <- with_baseline(data$validation$Y, data$truth, grid$noise_scale[g])
  do.call(rbind, lapply(names(variants), function(v) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(fit_multimodal_yfb(Ytr, data$training$time, data$training$event,
                                       grid$K_init[g], control = variants[[v]]),
                    error = function(e) e)
    runtime <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      return(data.frame(grid[g, ], variant = v, error = conditionMessage(fit),
                        runtime_seconds = runtime))
    }
    cbind(grid[g, ], variant = v, error = NA_character_, score_fit(fit, data, Yva),
          runtime_seconds = runtime)
  }))
}

# Socket cluster (fresh R processes): forked workers segfault on macOS when
# the BLAS library has been used in the parent process
cl <- parallel::makeCluster(cores)
parallel::clusterExport(cl, c("grid", "variants", "with_baseline", "score_fit", "run_one"))
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(survival))
  for (f in c("code/multimodal_yfb_helpers.R", "code/preprocess_multimodal_yfb.R",
              "code/multimodal_yfb_updates.R", "code/fit_multimodal_yfb.R",
              "code/predict_multimodal_yfb.R", "code/simulate_multimodal_yfb.R")) source(f)
}))
results <- parallel::parLapplyLB(cl, seq_len(nrow(grid)), function(g) {
  tryCatch(run_one(g), error = function(e) structure(conditionMessage(e), class = "try-error"))
})
parallel::stopCluster(cl)
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) warning(sum(failed), " datasets failed in a worker: ",
                         paste(vapply(results[failed], as.character, ""), collapse = " | "))
all_cols <- unique(unlist(lapply(results[!failed], names)))
results <- do.call(rbind, lapply(results[!failed], function(x) {
  x[setdiff(all_cols, names(x))] <- NA
  x[all_cols]
}))
write.csv(results, file.path(out_dir, "pruning_comparison.csv"), row.names = FALSE)

errors <- results[!is.na(results$error), ]
if (nrow(errors)) {
  cat("Fits that errored:\n"); print(errors[, c("scenario", "noise_scale", "K_init", "seed", "variant", "error")])
}
ok <- results[is.na(results$error), ]
summary_tab <- aggregate(
  cbind(K_final, K_survival, false_survival, factors_recovered, prognostic_recall,
        direction_recall, c_validation, converged, runtime_seconds) ~
    noise_scale + scenario + K_init + variant, data = ok, FUN = mean)
summary_tab <- summary_tab[order(summary_tab$noise_scale, summary_tab$scenario,
                                 summary_tab$variant, summary_tab$K_init), ]
print(summary_tab, digits = 3, row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "pruning_comparison_summary.csv"), row.names = FALSE)
