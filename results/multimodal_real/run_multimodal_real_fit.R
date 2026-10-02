# ============================================================
# Script: run_multimodal_real_fit.R
# Purpose: First real-data application of the multimodal YFB model: train on
#          TCGA (expression + methylation), validate on ICGC, and compare with
#          expression-only YFB and two-step EBMF -> Cox baselines.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: survival, flashier, ebnm, yaml; code/fit_multimodal_yfb.R,
#               code/fit_cox_on_yf.R, code/run_multimodal_yfb_simulation.R,
#               code/load_multiomics_data.R
# Run: Rscript results/multimodal_real/run_multimodal_real_fit.R [K ...]
#      (default K = 5 7 10). Each K is cached in outputs/fits/.
# ============================================================
#
# Inputs are on the scales the current nonnegative point-exponential model
# accepts: log2 expression and methylation beta values (decision 2, 10/1).
# All arms see the same training subjects and screened features. The
# expression-only YFB arm uses the same raw log2 inputs, not the per-platform
# z-standardized inputs of the recommended single-modality configuration, so
# the comparison isolates the effect of adding methylation.
#
# C-index: Harrell's C with higher risk = shorter survival (reverse = TRUE),
# risk orientation frozen from the fit (no re-orientation on validation data).
#
# Memory/time: TCGA is 144 x (3,000 + 10,000). One joint fit holds a few
# n x p copies (~15 MB each), so memory is small; runtime is dominated by the
# per-feature tau and F updates over 13k features x up to 500 sweeps.

suppressPackageStartupMessages(library(survival))
source("code/multimodal_yfb_helpers.R")
source("code/preprocess_multimodal_yfb.R")
source("code/multimodal_yfb_updates.R")
source("code/fit_multimodal_yfb.R")
source("code/predict_multimodal_yfb.R")
source("code/simulate_multimodal_yfb.R")          # multimodal_yfb_survival_metrics()
source("code/update_beta.R"); source("code/update_L.R"); source("code/update_F.R")
source("code/update_tau.R"); source("code/compute_elbo.R")
# fit_cox_on_yf.R ends with a runner block that errors when real_Y is NULL;
# the function definitions are complete before it fires (same as run_tests.R)
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL)))
source("code/predict_cox_on_yf.R")
source("code/run_multimodal_yfb_simulation.R")    # baseline fitters
source("code/concordance_ci.R")
source("code/load_multiomics_data.R")

args <- commandArgs(trailingOnly = TRUE)
K_values <- if (length(args)) as.integer(args) else c(5L, 7L, 10L)
out_dir <- "results/multimodal_real/outputs"
dir.create(file.path(out_dir, "fits"), recursive = TRUE, showWarnings = FALSE)

cache <- "data/processed_multiomics_tcga_icgc.rds"
d <- if (file.exists(cache)) readRDS(cache) else build_multiomics_cohorts()
tr <- d$training
validation <- list(icgc_primary = d$validation_primary, icgc_all = d$validation_all)

# Evaluation ----

#' C-index (and bootstrap 95% CI) of a frozen risk score on each cohort
#' @param risk_fn Function mapping a cohort's `Y` to a risk score vector.
#' @return data.frame with one row per evaluation cohort.
score_cohorts <- function(risk_fn) {
  cohorts <- c(list(tcga_train = tr), validation)
  do.call(rbind, lapply(names(cohorts), function(nm) {
    x <- cohorts[[nm]]
    risk <- risk_fn(x$Y)
    ci <- bootstrap_concordance_ci(risk, x$time, x$event, B = 1000, seed = 1)
    data.frame(cohort = nm, n = length(risk), events = sum(x$event),
               c_index = multimodal_yfb_survival_metrics(risk, x$time, x$event)$c_index,
               ci_lower = ci$lower, ci_upper = ci$upper)
  }))
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  attr(value, "runtime_seconds") <- proc.time()[["elapsed"]] - t0
  value
}

# Fits ----

rows <- list()
for (K in K_values) {
  fit_path <- file.path(out_dir, "fits", sprintf("fits_K%02d.rds", K))
  if (file.exists(fit_path)) {
    fits <- readRDS(fit_path)
  } else {
    message("K = ", K, ": fitting joint multimodal YFB")
    joint <- timed(fit_multimodal_yfb(tr$Y, tr$time, tr$event, K))
    message("K = ", K, ": fitting expression-only YFB")
    expr_only <- timed(fit_multimodal_yfb_single_modality(
      tr$Y, tr$time, tr$event, "expression", K))
    message("K = ", K, ": fitting two-step EBMF -> Cox (stacked modalities)")
    two_step <- timed(fit_multimodal_yfb_ebmf_cox(tr$Y, tr$time, tr$event, K))
    fits <- list(joint = joint, expr_only = expr_only, two_step = two_step)
    saveRDS(fits, fit_path)
  }

  dj <- fits$joint$diagnostics
  arms <- list(
    joint_multimodal_yfb = list(
      risk = function(Y) predict_multimodal_yfb(fits$joint, Y)$risk_scores,
      info = data.frame(converged = dj$converged, iterations = dj$iterations,
                        K_eff_reconstruction = dj$K_eff_reconstruction,
                        K_eff_survival = dj$K_eff_survival)),
    expression_only_yfb = list(
      risk = function(Y) predict_multimodal_yfb_single_modality(fits$expr_only, Y)$risk_scores,
      info = data.frame(converged = fits$expr_only$history$converged,
                        iterations = fits$expr_only$history$n_iter,
                        K_eff_reconstruction = sum(fits$expr_only$history$factor_pve[
                          fits$expr_only$history$n_iter, ] >= 0.01),
                        K_eff_survival = NA_integer_)),
    two_step_ebmf_cox = list(
      risk = function(Y) predict_multimodal_yfb_ebmf_cox(fits$two_step, Y),
      info = data.frame(converged = NA, iterations = NA,
                        K_eff_reconstruction = fits$two_step$n_factors,
                        K_eff_survival = NA_integer_))
  )
  for (arm in names(arms)) {
    sc <- score_cohorts(arms[[arm]]$risk)
    rows[[length(rows) + 1L]] <- cbind(
      K_init = K, method = arm, sc, arms[[arm]]$info,
      runtime_seconds = attr(fits[[c(joint_multimodal_yfb = "joint",
        expression_only_yfb = "expr_only", two_step_ebmf_cox = "two_step")[[arm]]]],
        "runtime_seconds") %||% NA_real_)
  }
}
results <- do.call(rbind, rows)
rownames(results) <- NULL
print(results, digits = 3)
write.csv(results, file.path(out_dir, "real_fit_cindex.csv"), row.names = FALSE)
