# ============================================================
# Script:  run_init_constraint_sweep.R
# Purpose: Phase 2 synthetic sweep — evaluate DeSurv-style init constraints
#          (init_y_normalize × init_constraint) on synthetic benchmark data.
#          Tests whether any combination of Y column normalization and column
#          sum-to-1 constraints on EL or EF at initialization improves ELBO
#          or external C-index over the default (none/FALSE) baseline.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Usage:   Rscript results/sandbox/init_constraints/run_init_constraint_sweep.R
# Output:  results/sandbox/init_constraints/sandbox_init_summary.csv
# ============================================================

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R")) {
  if (file.exists("../../code/fit_modular.R")) setwd("../..")
}

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
  library(ebnm)
})

cfg <- yaml::read_yaml("config/globals.yml")

source("code/update_beta.R")
source("code/update_L.R")
source("code/update_F.R")
source("code/update_tau.R")
source("code/compute_elbo.R")
source("code/fit_modular.R")
source("code/predict.R")
source("code/train_test_split.R")
source("results/benchmark_sim/benchmark_helpers.R")

OUT_DIR  <- "results/sandbox/init_constraints"
OUT_FILE <- file.path(OUT_DIR, "sandbox_init_summary.csv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Generate synthetic benchmark data (same DGP as the main benchmark)
# --------------------------------------------------------------------------
cat("Generating synthetic benchmark data ...\n")
syn <- generate_synthetic_benchmark_data(
  n = cfg$synthetic$n, p = cfg$synthetic$p,
  K_true = cfg$synthetic$k_true, seed = cfg$synthetic$seed
)
cat(sprintf("  n=%d, p=%d, K_true=%d, censoring=%.1f%%\n",
            syn$n, syn$p, syn$K_true, 100 * syn$censoring_rate))

set.seed(cfg$evaluation$holdout_seed)
split_idx <- stratified_split(syn$status,
                              test_frac = cfg$evaluation$holdout_frac,
                              seed      = cfg$evaluation$holdout_seed)
Y_tr  <- syn$Y[split_idx$train, ]
Y_te  <- syn$Y[split_idx$test,  ]
t_tr  <- syn$time[split_idx$train];   s_tr <- syn$status[split_idx$train]
t_te  <- syn$time[split_idx$test];    s_te <- syn$status[split_idx$test]

# --------------------------------------------------------------------------
# Sweep design: 2 × 3 = 6 conditions
# --------------------------------------------------------------------------
init_y_norm_grid  <- c(FALSE, TRUE)
init_constraint_grid <- c("none", "col_l_sum1", "col_f_sum1")
K     <- cfg$benchmark$k_pdac_synthetic
PRIOR <- "point_normal"

cat(sprintf("\nSweeping %d conditions (init_y_normalize × init_constraint) ...\n",
            length(init_y_norm_grid) * length(init_constraint_grid)))
cat(sprintf("K=%d | prior=%s | max_iter=%d\n\n", K, PRIOR, cfg$cavi$max_iter))

results <- list()
cond_i  <- 0

for (yn in init_y_norm_grid) {
  for (ic in init_constraint_grid) {
    cond_i <- cond_i + 1
    label  <- sprintf("y_norm=%s | constraint=%s", yn, ic)
    cat(sprintf("[%d/6] %s ...\n", cond_i, label))

    fit <- fit_supervised_mf_modular(
      Y_tr, t_tr, s_tr,
      K               = K,
      max_iter        = cfg$cavi$max_iter,
      tol             = cfg$cavi$tol,
      prior_LF        = "point_exponential",
      prior_beta      = PRIOR,
      alpha           = cfg$benchmark$alpha,
      lambda          = cfg$benchmark$lambda,
      N_burnin        = cfg$benchmark$n_burnin,
      normalize_AB    = cfg$benchmark$normalize_ab,
      sign_correction = TRUE,
      init_y_normalize  = yn,
      init_constraint   = ic,
      verbose         = FALSE
    )

    pred  <- predict_supervised_mf(Y_te, fit$EF, fit$EBeta)
    c_idx <- as.numeric(
      concordance(Surv(t_te, s_te) ~ pred$risk_scores)$concordance
    )
    k_eff   <- sum(abs(fit$EBeta) > cfg$k_selection$beta_threshold)
    elbo    <- tail(fit$history$elbo_full, 1)

    cat(sprintf("  → iter=%d | ELBO=%.1f | K_eff=%d | C=%.4f\n",
                fit$history$n_iter, elbo, k_eff, c_idx))

    results[[cond_i]] <- data.frame(
      init_y_normalize = yn,
      init_constraint  = ic,
      final_elbo       = round(elbo, 2),
      k_eff            = k_eff,
      beta_max         = round(max(abs(fit$EBeta)), 4),
      c_index          = round(c_idx, 4),
      n_iters          = fit$history$n_iter,
      converged        = fit$history$converged,
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, results)
write.csv(df, OUT_FILE, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n\n", OUT_FILE))

# --------------------------------------------------------------------------
# Summary table
# --------------------------------------------------------------------------
baseline_elbo <- df$final_elbo[df$init_y_normalize == FALSE & df$init_constraint == "none"]
baseline_c    <- df$c_index[df$init_y_normalize == FALSE & df$init_constraint == "none"]

cat("=== Phase 2 Init Constraint Sweep — Synthetic Results ===\n")
cat(sprintf("  Baseline (none/FALSE): ELBO=%.1f | C=%.4f\n\n", baseline_elbo, baseline_c))
cat(sprintf("  %-10s  %-15s  %8s  %8s  %6s  %6s\n",
            "y_norm", "constraint", "ELBO", "dELBO", "K_eff", "C"))
cat(sprintf("  %s\n", paste(rep("-", 58), collapse="")))
for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  cat(sprintf("  %-10s  %-15s  %8.1f  %8.1f  %6d  %6.4f\n",
              r$init_y_normalize, r$init_constraint,
              r$final_elbo, r$final_elbo - baseline_elbo,
              r$k_eff, r$c_index))
}

best_idx  <- which.max(df$final_elbo)
best_c_idx <- which.max(df$c_index)
cat(sprintf("\n  Best ELBO: %s + %s (ELBO=%.1f, delta=%.1f, C=%.4f)\n",
            df$init_y_normalize[best_idx], df$init_constraint[best_idx],
            df$final_elbo[best_idx],
            df$final_elbo[best_idx] - baseline_elbo,
            df$c_index[best_idx]))
cat(sprintf("  Best C:    %s + %s (C=%.4f, ELBO=%.1f)\n",
            df$init_y_normalize[best_c_idx], df$init_constraint[best_c_idx],
            df$c_index[best_c_idx], df$final_elbo[best_c_idx]))

cat("\nSweep complete.\n")
