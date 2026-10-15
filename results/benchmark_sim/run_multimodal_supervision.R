# ============================================================
# Script: run_multimodal_supervision.R
# Purpose: Workstream C6 (M4): does the single-modality survival supervision
#          (supervised start + shared point-Laplace beta prior with warm-up)
#          carry into the multimodal model? Planted low-variance prognostic
#          program, both modalities, calibrated test bed.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, parallel; code/fit_multimodal_yfb.R and helpers,
#               code/simulate_multimodal_yfb.R
# Run: Rscript results/benchmark_sim/run_multimodal_supervision.R [n_seeds] [cores] [K_init]
# ============================================================
#
# Data: simulate_multimodal_yfb_data("low_variance_prognostic_factor",
#   n_train = 400, n_validation = 120, p = 50 + 50, beta = c(0, 0, 3.3),
#   truncate = FALSE, noise_scale 0.25 and 1) -- the calibrated bed
#   (DECISIONS.md 2026-10-15). Raw Gaussian data; both arms fit an intercept.
#
# Arms (fixed before running):
#   pruned      intercept + point-Laplace F + ebnm priors + ELBO pruning (10/1)
#   pruned_sup  pruned + init_supervised + shared point-Laplace beta prior,
#               10-iteration warm-up
#
# Metrics (as single modality): recovered = best |cor| between stacked
# planted loadings (expression + methylation) and an estimated factor >= 0.7
# with the correct coefficient sign; survival-active = fit's
# |E beta| >= 1.96 SD; false = survival-active factors other than the
# planted program's match; converged; K_final; validation C.
#
# Parallelization: socket cluster, one task per (seed, noise, arm).

suppressPackageStartupMessages({ library(survival); library(parallel) })
args    <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 20L
cores   <- if (length(args) >= 2) as.integer(args[2]) else 8L
K_INIT  <- if (length(args) >= 3) as.integer(args[3]) else 7L
out_dir <- sprintf("results/benchmark_sim/outputs/multimodal_supervision/K%d", K_INIT)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

laplace <- list(expression = "point_laplace", methylation = "point_laplace")
base_ctrl <- list(intercept = TRUE, prior_update = "ebnm", prior_F = laplace,
                  prune = TRUE, survival_elbo = "plugin")
arms <- list(
  pruned     = base_ctrl,
  pruned_sup = c(base_ctrl, list(init_supervised = TRUE,
                                 beta_prior = "shared_point_laplace", beta_warmup = 10L))
)

load_code <- function() {
  suppressPackageStartupMessages(library(survival))
  for (f in c("code/multimodal_yfb_helpers.R", "code/preprocess_multimodal_yfb.R",
              "code/multimodal_yfb_updates.R", "code/fit_multimodal_yfb.R",
              "code/predict_multimodal_yfb.R", "code/simulate_multimodal_yfb.R")) source(f)
}

run_cell <- function(seed, noise, arm) {
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor", n_train = 400L,
                                    n_validation = 120L, p_expression = 50L, p_methylation = 50L,
                                    seed = seed, truncate = FALSE, noise_scale = noise,
                                    beta = c(0, 0, 3.3))
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(fit_multimodal_yfb(d$training$Y, d$training$time, d$training$event,
                                     K_INIT, control = arms[[arm]]), error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error"))
    return(data.frame(seed = seed, noise = noise, arm = arm, error = conditionMessage(fit),
                      runtime = secs))
  true_F <- do.call(rbind, d$truth$F); est_F <- do.call(rbind, fit$EF)
  r <- suppressWarnings(stats::cor(true_F, est_F)); r[!is.finite(r)] <- 0
  j <- which(d$truth$beta != 0)
  best <- apply(abs(r), 1, which.max); best_cor <- apply(abs(r), 1, max)
  k <- best[j]
  active <- fit$diagnostics$survival_active
  matched <- if (best_cor[j] >= 0.7) k else integer(0)
  risk <- predict_multimodal_yfb(fit, d$validation$Y)$risk_scores
  pve_planted <- {
    S <- lapply(names(d$truth$F), function(m)
      scale(outer(d$training$L[, j], d$truth$F[[m]][, j]), scale = FALSE))
    sum(sapply(S, function(x) sum(x^2))) /
      sum(sapply(d$training$Y, function(x) sum(scale(x, scale = FALSE)^2)))
  }
  data.frame(seed = seed, noise = noise, arm = arm, error = NA_character_, runtime = secs,
             pve_planted = pve_planted, cor_f1 = best_cor[1], cor_f2 = best_cor[2],
             cor_planted = best_cor[j],
             recovered = best_cor[j] >= 0.7 && sign(fit$EBeta[k] * r[j, k]) == sign(d$truth$beta[j]),
             match_active = active[k], n_active = sum(active),
             false_active = sum(active[setdiff(seq_along(active), matched)]),
             K_final = ncol(fit$EL), K_reconstruction = fit$diagnostics$K_eff_reconstruction,
             K_survival = fit$diagnostics$K_eff_survival,
             converged = fit$diagnostics$converged, iterations = fit$diagnostics$iterations,
             max_abs_beta = max(abs(fit$EBeta)),
             c_validation = multimodal_yfb_survival_metrics(risk, d$validation$time,
                                                            d$validation$event)$c_index)
}

grid <- expand.grid(seed = 20261015L + seq_len(n_seeds), noise = c(0.25, 1),
                    arm = names(arms), stringsAsFactors = FALSE)
cl <- makeCluster(cores)
clusterExport(cl, c("arms", "K_INIT", "load_code", "run_cell", "grid"))
invisible(clusterEvalQ(cl, load_code()))
rows <- parLapplyLB(cl, seq_len(nrow(grid)),
                    function(i) run_cell(grid$seed[i], grid$noise[i], grid$arm[i]))
stopCluster(cl)
cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(res, list(res$arm, res$noise), drop = TRUE), function(d0) {
  d <- d0[is.na(d0$error), ]
  data.frame(arm = d0$arm[1], noise = d0$noise[1], n = nrow(d0), errors = sum(!is.na(d0$error)),
             recovered = sum(d$recovered), cor_planted_median = median(d$cor_planted),
             cor_f1_median = median(d$cor_f1), cor_f2_median = median(d$cor_f2),
             false_active_mean = mean(d$false_active), K_final_median = median(d$K_final),
             K_survival_median = median(d$K_survival), converged = sum(d$converged),
             max_abs_beta = max(d$max_abs_beta), c_validation_mean = mean(d$c_validation, na.rm = TRUE),
             pve_planted_max = max(d$pve_planted), runtime_median = median(d$runtime))
}))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
