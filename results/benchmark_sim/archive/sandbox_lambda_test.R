# ============================================================
# Script:       sandbox_lambda_test.R
# Purpose:      Test three lambda conditions on the synthetic
#               benchmark to assess whether lambda = p/n survival
#               scaling improves beta recovery and hold-out C-index.
#               Sandbox only — no production changes.
# Author:       Claude Code (reviewed by Andrew Walther)
# Created:      2026-04-24
# Run from repo root:
#   Rscript results/benchmark_sim/sandbox_lambda_test.R
# ============================================================

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages({
  library(survival)
  library(yaml)
})

source("code/train_test_split.R")
source("code/predict.R")
suppressMessages(tryCatch(source("code/fit_modular.R"), error = function(e) invisible(NULL)))
source("code/select_alpha_cv.R")
source("results/benchmark_sim/run_ssbmf_benchmark.R")

cfg <- yaml::read_yaml("config/globals.yml")
set.seed(cfg$synthetic$seed)

# ---- Data generation ----------------------------------------
cat("Generating synthetic benchmark data ...\n")
data <- generate_synthetic_benchmark_data(
  n        = cfg$synthetic$n,
  p        = cfg$synthetic$p,
  K_true   = cfg$synthetic$k_true,
  seed     = cfg$synthetic$seed
)
beta_true <- data$beta_true
cat(sprintf("  n=%d  p=%d  K_true=%d\n", data$n, data$p, length(beta_true)))
cat(sprintf("  beta_true: [%s]\n",
            paste(sprintf("%.2f", beta_true), collapse=", ")))

# ---- Train/test split ---------------------------------------
splits <- stratified_split(data$status,
                           test_frac = 0.2, seed = cfg$synthetic$seed)
Y_train      <- data$Y[splits$train_idx, ]
time_train   <- data$time[splits$train_idx]
status_train <- data$status[splits$train_idx]
Y_test       <- data$Y[splits$test_idx, ]
time_test    <- data$time[splits$test_idx]
status_test  <- data$status[splits$test_idx]
n_train      <- nrow(Y_train)
p            <- ncol(Y_train)

# ---- Alpha CV (shared across all lambda conditions) ---------
cat("\nRunning alpha CV (shared, lambda=1 for selection) ...\n")
cv_res <- select_alpha_cv(
  Y_train, time_train, status_train,
  alpha_grid = cfg$cavi$alpha_grid,
  n_folds    = 5,
  K_max      = cfg$cavi$k_max,
  use_1se    = TRUE,
  seed       = cfg$synthetic$seed,
  max_iter   = cfg$cavi$max_iter,
  tol        = cfg$cavi$tol,
  prior_beta = "point_normal",
  verbose    = FALSE
)
alpha_opt <- cv_res$alpha_opt
cat(sprintf("  alpha_opt = %.2f\n", alpha_opt))

# ---- Lambda conditions to test ------------------------------
lambda_conditions <- c(
  "lambda=1 (current)"     = 1.0,
  "lambda=p/n (principled)" = p / n_train,
  "lambda=2*(p/n)"         = 2 * p / n_train
)
cat(sprintf("\nLambda values: 1.0 | %.2f | %.2f\n",
            p / n_train, 2 * p / n_train))

# ---- Run all three conditions --------------------------------
results <- lapply(names(lambda_conditions), function(label) {
  lam <- lambda_conditions[label]
  cat(sprintf("\n--- %s ---\n", label))

  fit <- fit_supervised_mf_modular(
    Y_train, time_train, status_train,
    K          = cfg$cavi$k_max,
    alpha      = alpha_opt,
    lambda     = lam,
    max_iter   = cfg$cavi$max_iter,
    tol        = cfg$cavi$tol,
    prior_beta = "point_normal",
    verbose    = FALSE
  )

  # K_eff
  pve_final  <- as.numeric(tail(fit$history$factor_pve, 1))
  active     <- (abs(fit$EBeta) > 0.05) | (pve_final > 0.01)
  K_eff      <- sum(active)

  # Beta recovery: RMSE vs truth (match estimated to true factors by loading correlation)
  est_beta <- fit$EBeta
  K_min  <- min(length(beta_true), length(est_beta))
  cormat <- cor(data$L[splits$train_idx, seq_len(K_min)],
                fit$EL[, seq_len(K_min)])
  cormat[is.na(cormat)] <- 0   # replace NA (zero-variance columns) with 0

  # Greedy matching: each true factor gets the best available estimated factor
  est_order <- integer(K_min)
  used      <- rep(FALSE, K_min)
  for (ki in seq_len(K_min)) {
    scores     <- abs(cormat[ki, ])
    scores[used] <- -Inf
    best_j       <- which.max(scores)
    est_order[ki] <- best_j
    used[best_j]  <- TRUE
  }
  sign_flip        <- sign(diag(cormat[seq_len(K_min), est_order]))
  sign_flip[sign_flip == 0] <- 1
  beta_est_matched <- sign_flip * est_beta[est_order]
  beta_rmse <- sqrt(mean((beta_est_matched[seq_len(K_min)] - beta_true[seq_len(K_min)])^2))

  # Hold-out C-index
  pred   <- predict_supervised_mf(Y_test, fit$EF, fit$EBeta)
  c_test <- as.numeric(concordance(
    Surv(time_test, status_test) ~ I(-pred$risk_scores)
  )$concordance)

  cat(sprintf("  K_eff=%d | beta RMSE=%.4f | hold-out C=%.4f | iters=%d\n",
              K_eff, beta_rmse, c_test, fit$history$n_iter))
  cat(sprintf("  EBeta (matched): [%s]\n",
              paste(sprintf("%+.3f", beta_est_matched[seq_len(K_min)]), collapse=", ")))
  cat(sprintf("  beta_true:       [%s]\n",
              paste(sprintf("%+.3f", beta_true[seq_len(K_min)]), collapse=", ")))

  data.frame(
    Label      = label,
    Lambda     = lam,
    K_eff      = K_eff,
    Beta_RMSE  = round(beta_rmse, 4),
    C_holdout  = round(c_test, 4),
    N_iter     = fit$history$n_iter,
    Converged  = fit$history$converged,
    stringsAsFactors = FALSE
  )
})

# ---- Summary table ------------------------------------------
summary_df <- do.call(rbind, results)
cat("\n\n=== LAMBDA SANDBOX RESULTS ===\n")
print(summary_df, row.names = FALSE)
cat(sprintf("\nbaseline (lambda=1) beta RMSE: %.4f | C: %.4f\n",
            summary_df$Beta_RMSE[1], summary_df$C_holdout[1]))
cat(sprintf("principled (p/n=%.1f) beta RMSE: %.4f | C: %.4f  [delta RMSE: %+.4f, delta C: %+.4f]\n",
            p / n_train,
            summary_df$Beta_RMSE[2], summary_df$C_holdout[2],
            summary_df$Beta_RMSE[2] - summary_df$Beta_RMSE[1],
            summary_df$C_holdout[2] - summary_df$C_holdout[1]))
cat(sprintf("over-scaled (2p/n=%.1f) beta RMSE: %.4f | C: %.4f  [delta RMSE: %+.4f, delta C: %+.4f]\n",
            2 * p / n_train,
            summary_df$Beta_RMSE[3], summary_df$C_holdout[3],
            summary_df$Beta_RMSE[3] - summary_df$Beta_RMSE[1],
            summary_df$C_holdout[3] - summary_df$C_holdout[1]))
