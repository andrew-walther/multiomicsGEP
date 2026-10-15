# ============================================================
# Script: run_supervision_diagnostics.R
# Purpose: Test bed for survival supervision of the single-modality loadings
#          (Workstream C): does fit_cox_on_yf() recover a planted
#          low-variance prognostic program, and how much does survival enter
#          the F update?
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, parallel; code/fit_cox_on_yf.R,
#               code/predict_cox_on_yf.R, code/concordance_ci.R,
#               code/simulate_multimodal_yfb.R
# Run: Rscript results/benchmark_sim/run_supervision_diagnostics.R <arm_set> [n_seeds] [cores]
# ============================================================
#
# Data (fixed for every arm set): simulate_multimodal_yfb_data(
#   "low_variance_prognostic_factor", noise_scale = 0.25, truncate = FALSE),
# expression block only (n = 120 training, 120 validation, p = 50, 3 true
# programs). Program 3 is the planted prognostic program: score scale 0.15,
# true coefficient 1.10; programs 1-2 explain the variance and have
# coefficient 0. Each cohort is z-standardized by its own feature means and
# SDs, as the real-data pipeline standardizes each platform.
#
# Metrics (pre-specified, Workstream C success criteria):
#   recovered     best |cor(F_true[,3], EF_k)| >= 0.7 AND the matched
#                 program's coefficient has the true sign, i.e.
#                 sign(EBeta_k * cor) == sign(1.10)
#   active_z      program k is survival-active: |E[beta_k]| >= 1.96 sd(beta_k)
#   false_z       active_z programs other than the planted program's match
#                 (counted only when the match has |cor| >= 0.7)
#   false_thresh  same with classify_factors()'s |E[beta_k]| > 0.001 rule
#   c_validation  C on the validation cohort, orientation frozen at fit time
#   n_programs    programs with nonzero loadings (||EF_k||^2 > 1e-8)
#   converged     fit_cox_on_yf() convergence flag; runaway: max over
#                 iterations of ||F_k|| or ||L_k|| relative to iteration 1
#
# Parallelization: socket cluster, one task per (seed, arm). Fits are small
# (n = 120, p = 50, K = 7) and single-threaded.

suppressPackageStartupMessages({ library(survival); library(parallel) })

ARM_SETS <- list(
  # C1: the three update arms as they exist on main, plus the recommended
  # real-data prior (point-exponential L and F) as a reference
  c1 = list(
    alpha_F_0   = list(alpha_F = 0),
    alpha_F_0.5 = list(alpha_F = 0.5),
    equal_wt    = list(w_S = 1),
    ref_point_exponential = list(alpha_F = 0, prior_LF = "point_exponential")
  )
)

args    <- commandArgs(trailingOnly = TRUE)
arm_set <- if (length(args) >= 1) args[1] else "c1"
n_seeds <- if (length(args) >= 2) as.integer(args[2]) else 10L
cores   <- if (length(args) >= 3) as.integer(args[3]) else 6L
if (!arm_set %in% names(ARM_SETS)) stop("Unknown arm set: ", arm_set)
arms    <- ARM_SETS[[arm_set]]
K_FIT   <- 7L
out_dir <- file.path("results/benchmark_sim/outputs/supervision_diagnostics", arm_set)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# prior_LF = point_laplace: on z-standardized (signed, centered) data the
# recommended point-exponential L/F prior collapses to zero programs in most
# seeds of this test bed (see the ref_point_exponential arm; DECISIONS.md
# 2026-10-15), so the arms use the signed point-Laplace prior.
FIT_DEFAULTS <- list(K = K_FIT, max_iter = 300L, tol = 1e-5, prior_beta = "normal",
                     prior_LF = "point_laplace",
                     alpha = 0.5, init_method = "svd", verbose = FALSE,
                     trace_supervision = TRUE)

load_code <- function() {
  suppressPackageStartupMessages({ library(survival); library(ebnm) })
  suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
  source("code/predict_cox_on_yf.R"); source("code/concordance_ci.R")
  source("code/simulate_multimodal_yfb.R")
}

#' Build the expression-only test-bed dataset for one seed
#' @param seed Integer RNG seed.
#' @return list(Ytr, Yva, time/status for both cohorts, F_true, beta_true, pve3)
make_testbed <- function(seed) {
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
                                    n_train = 120L, n_validation = 120L,
                                    p_expression = 50L, p_methylation = 50L,
                                    seed = seed, truncate = FALSE, noise_scale = 0.25)
  Ytr_raw <- d$training$Y$expression
  F_true  <- d$truth$F$expression
  # Planted program's share of the centered training variance: rank-one signal
  # L_3 F_3' (centered) over the centered total sum of squares
  S3   <- scale(outer(d$training$L[, 3], F_true[, 3]), scale = FALSE)
  pve3 <- sum(S3^2) / sum(scale(Ytr_raw, scale = FALSE)^2)
  list(Ytr = scale(Ytr_raw), Yva = scale(d$validation$Y$expression),
       ttr = d$training$time, str = d$training$event,
       tva = d$validation$time, sva = d$validation$event,
       F_true = F_true, beta_true = d$truth$beta, pve3 = pve3,
       c_oracle = frozen_reverse_cindex(d$validation$eta, d$validation$time,
                                        d$validation$event))
}

#' Score one fit against the truth (metrics defined in the header)
score_fit <- function(fit, tb) {
  r <- suppressWarnings(stats::cor(tb$F_true, fit$EF))
  r[!is.finite(r)] <- 0
  best     <- apply(abs(r), 1, which.max)
  best_cor <- apply(abs(r), 1, max)
  j <- which(tb$beta_true != 0)                       # planted program (3)
  k <- best[j]
  beta_sd  <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  active_z <- abs(fit$EBeta) >= 1.96 * beta_sd & beta_sd > 0
  matched  <- if (best_cor[j] >= 0.7) k else integer(0)
  tr <- fit$history$trace
  rel <- function(m) max(apply(m, 2, function(v) max(v, na.rm = TRUE) / max(v[1], 1e-12)))
  risk <- predict_cox_on_yf(tb$Yva, fit$EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
  elbo <- fit$history$elbo_full
  data.frame(
    pve_planted   = tb$pve3,
    cor_f1 = best_cor[1], cor_f2 = best_cor[2], cor_planted = best_cor[j],
    match_k       = k,
    n_programs    = sum(colSums(fit$EF^2) > 1e-8),
    beta_match    = fit$EBeta[k],
    beta_z_match  = fit$EBeta[k] / max(beta_sd[k], 1e-12),
    recovered     = best_cor[j] >= 0.7 &&
                    sign(fit$EBeta[k] * r[j, k]) == sign(tb$beta_true[j]),
    match_active_z = active_z[k],
    n_active_z    = sum(active_z),
    false_z       = sum(active_z[setdiff(seq_along(active_z), matched)]),
    false_thresh  = sum((abs(fit$EBeta) > 0.001)[setdiff(seq_along(fit$EBeta), matched)]),
    c_validation  = frozen_reverse_cindex(risk, tb$tva, tb$sva),
    c_oracle      = tb$c_oracle,
    converged     = fit$history$converged,
    n_iter        = fit$history$n_iter,
    elbo_decreases = sum(diff(elbo) < -1e-8 * abs(head(elbo, -1))),
    runaway_F     = rel(tr$norm_F),
    runaway_L     = rel(tr$norm_L),
    shareA_final_match = tail(tr$shareA_median[, k], 1),
    shareA_unweighted_final_match = tail(tr$shareA_unweighted_median[, k], 1),
    row.names = NULL
  )
}

#' Fit one (seed, arm) cell, save its trace, return the score row
run_cell <- function(seed, arm) {
  tb  <- make_testbed(seed)
  t0  <- proc.time()[["elapsed"]]
  fit <- tryCatch(do.call(fit_cox_on_yf, c(list(tb$Ytr, tb$ttr, tb$str),
                                           modifyList(FIT_DEFAULTS, arms[[arm]]))),
                  error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error"))
    return(data.frame(seed = seed, arm = arm, error = conditionMessage(fit), runtime = secs))
  saveRDS(list(trace = fit$history$trace, elbo = fit$history$elbo_full,
               EBeta = fit$EBeta, EBeta2 = fit$EBeta2, EF = fit$EF,
               F_true = tb$F_true),
          file.path(out_dir, sprintf("trace_seed%d_%s.rds", seed, arm)))
  cbind(data.frame(seed = seed, arm = arm, error = NA_character_, runtime = secs),
        score_fit(fit, tb))
}

seeds <- 20261015L + seq_len(n_seeds)
grid  <- expand.grid(seed = seeds, arm = names(arms), stringsAsFactors = FALSE)

load_code()
cat(sprintf("Arm set %s: %d seeds x %d arms = %d fits on %d cores\n",
            arm_set, n_seeds, length(arms), nrow(grid), cores))
cl <- makeCluster(cores)
clusterExport(cl, c("arms", "FIT_DEFAULTS", "out_dir", "load_code", "make_testbed",
                    "score_fit", "run_cell", "grid"))
invisible(clusterEvalQ(cl, load_code()))
rows <- parLapplyLB(cl, seq_len(nrow(grid)), function(i) run_cell(grid$seed[i], grid$arm[i]))
stopCluster(cl)

all_cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(r) { r[setdiff(all_cols, names(r))] <- NA; r[all_cols] }))
n_err <- sum(!is.na(res$error))
if (n_err > 0) warning(n_err, " fits failed: ", paste(unique(res$error[!is.na(res$error)]), collapse = " | "))
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)

ok <- is.na(res$error)
summ <- do.call(rbind, lapply(split(res[ok, ], res$arm[ok]), function(d) data.frame(
  arm = d$arm[1], n = nrow(d), recovered = sum(d$recovered),
  cor_planted_median = median(d$cor_planted), cor_f1_median = median(d$cor_f1),
  cor_f2_median = median(d$cor_f2), n_programs_median = median(d$n_programs), match_active_z = sum(d$match_active_z),
  false_z_mean = mean(d$false_z), false_thresh_mean = mean(d$false_thresh),
  c_validation_mean = mean(d$c_validation), c_oracle_mean = mean(d$c_oracle),
  converged = sum(d$converged), elbo_decreases_total = sum(d$elbo_decreases),
  runaway_F_max = max(d$runaway_F), runaway_L_max = max(d$runaway_L),
  shareA_match_median = median(d$shareA_final_match),
  shareA_unweighted_median = median(d$shareA_unweighted_final_match),
  pve_planted_max = max(d$pve_planted), runtime_median = median(d$runtime))))
write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summ, digits = 3)
