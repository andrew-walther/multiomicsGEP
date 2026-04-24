# ==============================================================================
# Script:       run_ssbmf_benchmark.R
# Purpose:      Phase 2 synthetic benchmark for the SBMF model.
#               Generates a corrected synthetic DGP, selects alpha via CV,
#               fits the final model, and writes benchmark tables/figures.
# Author:       Codex
# Created:      2026-04-23
# Dependencies: code/fit_modular.R, code/predict.R, code/train_test_split.R,
#               code/select_alpha_cv.R
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(survival)
})

if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (!file.exists("code/fit_modular.R") && file.exists("../code/fit_modular.R")) {
  setwd("..")
}

cfg <- yaml::read_yaml("config/globals.yml")

source("code/train_test_split.R")
source("code/predict.R")
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
source("code/select_alpha_cv.R")

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

save_plot_pair <- function(path_stub, width_pdf, height_pdf, width_png, height_png, draw_fun) {
  for (ext in c("pdf", "png")) {
    path <- paste0(path_stub, ".", ext)
    if (ext == "pdf") {
      pdf(path, width = width_pdf, height = height_pdf)
    } else {
      png(path, width = width_png, height = height_png, res = 140)
    }
    draw_fun()
    dev.off()
  }
}

calibrate_censor_scale <- function(event_times, base_censor, target = 0.30, n_iter = 40) {
  lo <- 1e-3
  hi <- max(event_times) * 10
  for (i in seq_len(n_iter)) {
    mid <- sqrt(lo * hi)
    censor_rate <- mean(base_censor * mid < event_times)
    if (censor_rate > target) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  hi
}

generate_synthetic_benchmark_data <- function(n = 300, p = 1000, K_true = 5,
                                              seed = 222, target_censoring = 0.30) {
  set.seed(seed)

  beta_true <- c(1.0, -0.8, 0, 0, 0)
  if (length(beta_true) != K_true) {
    stop("beta_true length must equal K_true.")
  }

  L_true <- matrix(rexp(n * K_true, rate = 1), n, K_true)
  signal_scale <- 0.25
  F_true <- matrix(0, p, K_true)
  for (k in seq_len(K_true)) {
    active <- sample.int(p, size = max(1L, round(0.05 * p)))
    F_true[active, k] <- if (k <= 2) signal_scale else 5 * signal_scale
  }
  tau_true <- rgamma(p, shape = 2, rate = 2)
  E <- sweep(matrix(rnorm(n * p), n, p), 2, sqrt(tau_true), "/") * signal_scale
  Y <- L_true %*% t(F_true) + E

  eta_true <- as.vector(L_true %*% beta_true)
  shape <- 1.5
  scale0 <- 0.01
  event_times <- (-log(runif(n)) / (scale0 * exp(eta_true)))^(1 / shape)
  censor_base  <- rexp(n, rate = 1)
  censor_scale <- calibrate_censor_scale(event_times, censor_base,
                                         target = target_censoring)
  censor_times <- censor_base * censor_scale
  time   <- pmin(event_times, censor_times)
  status <- as.integer(event_times <= censor_times)

  list(
    Y = Y,
    time = time,
    status = status,
    L_true = L_true,
    F_true = F_true,
    beta_true = beta_true,
    tau_true = tau_true,
    censoring_rate = mean(status == 0),
    shape = shape,
    scale0 = scale0,
    K_true = K_true
  )
}

greedy_rectangular_match <- function(cors) {
  score <- abs(cors)
  n_true <- nrow(cors)
  n_est  <- ncol(cors)
  n_match <- min(n_true, n_est)
  out <- vector("list", n_match)

  for (i in seq_len(n_match)) {
    idx <- which(score == max(score), arr.ind = TRUE)[1, ]
    r <- idx[1]
    c <- idx[2]
    out[[i]] <- data.frame(
      true_factor = r,
      est_factor = c,
      cor = cors[r, c],
      abs_cor = abs(cors[r, c]),
      sign = sign(cors[r, c])
    )
    score[r, ] <- -Inf
    score[, c] <- -Inf
  }

  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

compute_k_effective <- function(fit, beta_thresh = 0.05, pve_thresh = 0.01) {
  final_pve <- fit$history$factor_pve[nrow(fit$history$factor_pve), ]
  active <- (abs(fit$EBeta) > beta_thresh) | (final_pve > pve_thresh)
  list(
    K_effective = sum(active),
    final_pve = final_pve,
    active = active
  )
}

plot_factor_recovery <- function(cors, match_df, title_text) {
  cols <- colorRampPalette(c("#08306B", "#6BAED6", "#F7FBFF", "#FDBE85", "#A50F15"))(100)
  image(1:ncol(cors), 1:nrow(cors), t(abs(cors)),
        axes = FALSE, col = cols, zlim = c(0, 1),
        xlab = "Estimated factor", ylab = "True factor",
        main = title_text)
  axis(1, at = seq_len(ncol(cors)), labels = paste0("E", seq_len(ncol(cors))), las = 2, cex.axis = 0.6)
  axis(2, at = seq_len(nrow(cors)), labels = paste0("T", seq_len(nrow(cors))), las = 1, cex.axis = 0.8)
  box()
  for (i in seq_len(nrow(match_df))) {
    rect(match_df$est_factor[i] - 0.5, match_df$true_factor[i] - 0.5,
         match_df$est_factor[i] + 0.5, match_df$true_factor[i] + 0.5,
         border = ifelse(match_df$true_factor[i] <= 2, "#08519C", "#636363"),
         lwd = 2)
    text(match_df$est_factor[i], match_df$true_factor[i],
         sprintf("%.2f", match_df$cor[i]), cex = 0.7, col = "black")
  }
}

plot_ard_shrinkage <- function(factor_pve, match_df, active, title_text) {
  n_iter <- nrow(factor_pve)
  k_est  <- ncol(factor_pve)
  order_est <- c(match_df$est_factor[order(match_df$true_factor)],
                 setdiff(seq_len(k_est), match_df$est_factor))
  ordered <- factor_pve[, order_est, drop = FALSE]
  log_pve <- log10(ordered + 1e-8)
  cols <- rep("#BDBDBD", k_est)
  cols[seq_along(order_est) <= nrow(match_df)] <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD")[seq_len(nrow(match_df))]
  matplot(seq_len(n_iter), log_pve, type = "l", lty = 1, lwd = 2,
          col = cols, xlab = "Iteration", ylab = "log10(PVE + 1e-8)",
          main = title_text, bty = "n")
  abline(h = log10(0.01), col = "#444444", lty = 2)
  legend("topright",
         legend = c("Matched true factors", "Unmatched factors", "1% threshold"),
         col = c("#1F77B4", "#BDBDBD", "#444444"),
         lty = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n", cex = 0.8)
}

plot_beta_recovery <- function(beta_df, title_text) {
  x <- seq_len(nrow(beta_df))
  ylim <- range(c(beta_df$beta_true, beta_df$beta_est + 1.96 * beta_df$beta_sd,
                  beta_df$beta_est - 1.96 * beta_df$beta_sd)) * 1.15
  plot(x, beta_df$beta_est, pch = 16, col = "#1F77B4", ylim = ylim,
       xaxt = "n", xlab = "Aligned factor", ylab = expression(beta),
       main = title_text, bty = "n")
  axis(1, at = x, labels = beta_df$label, las = 2, cex.axis = 0.7)
  arrow_idx <- which(beta_df$ci_high > beta_df$ci_low)
  if (length(arrow_idx) > 0) {
    arrows(x[arrow_idx], beta_df$ci_low[arrow_idx],
           x[arrow_idx], beta_df$ci_high[arrow_idx],
           angle = 90, code = 3, length = 0.04, col = "#1F77B4", lwd = 1.5)
  }
  points(x, beta_df$beta_true, pch = 4, col = "#D62728", lwd = 2)
  abline(h = 0, col = "#666666", lty = 2)
  legend("topright", legend = c("Estimated (95% CI)", "True"),
         col = c("#1F77B4", "#D62728"), pch = c(16, 4), lwd = c(1, 2), bty = "n")
}

safe_cor_matrix <- function(true_mat, est_mat) {
  out <- matrix(0, nrow = ncol(true_mat), ncol = ncol(est_mat))
  for (i in seq_len(ncol(true_mat))) {
    x <- true_mat[, i]
    if (stats::sd(x) == 0) next
    for (j in seq_len(ncol(est_mat))) {
      y <- est_mat[, j]
      if (stats::sd(y) == 0) next
      out[i, j] <- suppressWarnings(stats::cor(x, y))
      if (!is.finite(out[i, j])) out[i, j] <- 0
    }
  }
  out
}

plot_elbo_trace <- function(history, title_text) {
  elbo_full <- if (!is.null(history$elbo_full)) history$elbo_full else history$ELBO_Full
  elbo_proxy <- if (!is.null(history$elbo_proxy)) history$elbo_proxy else history$ELBO_Proxy
  plot(elbo_full, type = "l", lwd = 2, col = "#1F77B4",
       xlab = "Iteration", ylab = "ELBO", main = title_text, bty = "n")
  lines(elbo_proxy, col = "#2CA02C", lwd = 2, lty = 2)
  legend("bottomright", legend = c("Full ELBO", "Genomics proxy"),
         col = c("#1F77B4", "#2CA02C"), lwd = 2, lty = c(1, 2), bty = "n")
}

plot_holdout_cindex <- function(cindex_df, title_text) {
  bar_cols <- c("#1F77B4", "#FF7F0E")
  barplot(cindex_df$C_Index, names.arg = cindex_df$Method, col = bar_cols,
          ylim = c(0, 1), main = title_text, ylab = "C-index", las = 2, bty = "n")
  abline(h = 0.5, col = "#666666", lty = 2)
}

run_ssbmf_benchmark <- function(output_root = "results/modular_sim_factor/ssbmf_benchmark",
                                synthetic_n = 300,
                                synthetic_p = 1000,
                                K_true = 5,
                                K_max = 15,
                                alpha_grid = c(0.1, 0.3, 0.5),
                                n_folds = 3,
                                max_iter = 80,
                                tol = 1e-4,
                                seed = cfg$synthetic$seed,
                                holdout_frac = 0.2,
                                quick = FALSE) {

  if (quick) {
    synthetic_n <- 120
    synthetic_p <- 300
    K_true <- 5
    K_max <- 8
    alpha_grid <- c(0.1, 0.3, 0.5)
    n_folds <- 3
    max_iter <- 35
    tol <- 1e-4
  }

  benchmark_root <- file.path(output_root, "synthetic")
  table_dir  <- file.path(benchmark_root, "tables")
  figure_dir  <- file.path(benchmark_root, "figures")
  ensure_dir(table_dir)
  ensure_dir(figure_dir)

  data <- generate_synthetic_benchmark_data(
    n = synthetic_n, p = synthetic_p, K_true = K_true, seed = seed
  )

  cat("=== SSBMF Synthetic Benchmark ===\n")
  cat(sprintf("  n=%d p=%d K_true=%d K_max=%d censoring=%.1f%%\n",
              synthetic_n, synthetic_p, K_true, K_max, 100 * data$censoring_rate))

  split <- stratified_split(data$status, test_frac = holdout_frac, seed = seed)
  train_idx <- split$train_idx
  test_idx  <- split$test_idx

  cv_res <- select_alpha_cv(
    data$Y[train_idx, , drop = FALSE],
    data$time[train_idx],
    data$status[train_idx],
    alpha_grid = alpha_grid,
    n_folds = n_folds,
    K_max = K_max,
    use_1se = TRUE,
    seed = seed,
    max_iter = max_iter,
    tol = tol
  )
  alpha_opt <- cv_res$alpha_opt

  holdout_fit <- fit_supervised_mf_modular(
    data$Y[train_idx, , drop = FALSE],
    data$time[train_idx],
    data$status[train_idx],
    K = K_max,
    alpha = alpha_opt,
    max_iter = max_iter,
    tol = tol,
    verbose = FALSE
  )
  holdout_pred <- predict_supervised_mf(data$Y[test_idx, , drop = FALSE],
                                        holdout_fit$EF, holdout_fit$EBeta)
  cindex_supervised <- as.numeric(concordance(
    Surv(data$time[test_idx], data$status[test_idx]) ~ I(-holdout_pred$risk_scores)
  )$concordance)

  pca_train <- prcomp(data$Y[train_idx, , drop = FALSE], rank. = min(5, K_max))
  pca_train_scores <- predict(pca_train, newdata = data$Y[train_idx, , drop = FALSE])
  pca_cox <- coxph(Surv(data$time[train_idx], data$status[train_idx]) ~ pca_train_scores)
  pca_test_scores <- predict(pca_train, newdata = data$Y[test_idx, , drop = FALSE])
  pca_lp_test <- as.vector(pca_test_scores %*% coef(pca_cox))
  cindex_pca <- as.numeric(concordance(
    Surv(data$time[test_idx], data$status[test_idx]) ~ I(-pca_lp_test)
  )$concordance)

  full_fit <- fit_supervised_mf_modular(
    data$Y,
    data$time,
    data$status,
    K = K_max,
    alpha = alpha_opt,
    max_iter = max_iter,
    tol = tol,
    verbose = FALSE
  )

  cors <- safe_cor_matrix(data$L_true, full_fit$EL)
  match_df <- greedy_rectangular_match(cors)
  est_order <- c(match_df$est_factor[order(match_df$true_factor)],
                 setdiff(seq_len(K_max), match_df$est_factor))
  aligned_est <- full_fit$EBeta[est_order]
  aligned_sd  <- sqrt(pmax(full_fit$EBeta2[est_order] - full_fit$EBeta[est_order]^2, 0))
  aligned_true <- c(data$beta_true[match_df$true_factor], rep(0, K_max - nrow(match_df)))
  beta_labels <- c(paste0("T", match_df$true_factor[order(match_df$true_factor)]),
                   paste0("N", seq_len(K_max - nrow(match_df))))
  beta_df <- data.frame(
    label = beta_labels,
    true_factor = c(match_df$true_factor[order(match_df$true_factor)], rep(NA_integer_, K_max - nrow(match_df))),
    est_factor  = est_order,
    beta_true   = aligned_true,
    beta_est    = aligned_est,
    beta_sd     = aligned_sd,
    ci_low      = aligned_est - 1.96 * aligned_sd,
    ci_high     = aligned_est + 1.96 * aligned_sd,
    ci_cover    = (aligned_true >= aligned_est - 1.96 * aligned_sd) &
                  (aligned_true <= aligned_est + 1.96 * aligned_sd)
  )

  final_k <- compute_k_effective(full_fit)
  final_pve_pct <- round(final_k$final_pve * 100, 2)
  ard_df <- data.frame(
    Factor = seq_len(K_max),
    Abs_Beta = round(abs(full_fit$EBeta), 4),
    PVE_Pct = round(final_pve_pct, 2),
    Active = final_k$active,
    Matched_True = {
      lookup <- rep(NA_integer_, K_max)
      lookup[match_df$est_factor] <- match_df$true_factor
      lookup
    }
  )

  holdout_df <- data.frame(
    Method = c("PCA", "Supervised SBMF"),
    C_Index = round(c(cindex_pca, cindex_supervised), 4)
  )

  elbo_df <- data.frame(
    Iteration = seq_along(full_fit$history$elbo_full),
    ELBO_Full = full_fit$history$elbo_full,
    ELBO_Proxy = full_fit$history$elbo_proxy,
    RMSE = full_fit$history$rmse
  )

  cors_df <- data.frame(
    true_factor = rep(seq_len(nrow(cors)), each = ncol(cors)),
    est_factor = rep(seq_len(ncol(cors)), times = nrow(cors)),
    cor = as.vector(cors),
    abs_cor = abs(as.vector(cors))
  )

  summary_df <- data.frame(
    Seed = seed,
    n = synthetic_n,
    p = synthetic_p,
    K_true = K_true,
    K_max = K_max,
    alpha_opt = alpha_opt,
    cv_rule = cv_res$selection_rule,
    cv_mean_cindex = cv_res$cv_table$mean_cindex[cv_res$cv_table$alpha == alpha_opt],
    holdout_supervised_cindex = cindex_supervised,
    holdout_pca_cindex = cindex_pca,
    K_effective = final_k$K_effective,
    censoring_rate = data$censoring_rate,
    converged = full_fit$history$converged,
    n_iter = full_fit$history$n_iter,
    final_elbo = tail(full_fit$history$elbo_full, 1),
    final_rmse = tail(full_fit$history$rmse, 1)
  )

  write.csv(cv_res$cv_table, file.path(table_dir, "alpha_cv_table.csv"), row.names = FALSE)
  write.csv(cv_res$fold_results, file.path(table_dir, "alpha_cv_fold_results.csv"), row.names = FALSE)
  write.csv(summary_df, file.path(table_dir, "benchmark_summary.csv"), row.names = FALSE)
  write.csv(cors_df, file.path(table_dir, "factor_recovery_correlations.csv"), row.names = FALSE)
  write.csv(match_df, file.path(table_dir, "factor_recovery_matches.csv"), row.names = FALSE)
  write.csv(ard_df, file.path(table_dir, "ard_summary.csv"), row.names = FALSE)
  write.csv(beta_df, file.path(table_dir, "beta_recovery_table.csv"), row.names = FALSE)
  write.csv(elbo_df, file.path(table_dir, "elbo_trace.csv"), row.names = FALSE)
  write.csv(holdout_df, file.path(table_dir, "holdout_cindex.csv"), row.names = FALSE)

  save_plot_pair(file.path(figure_dir, "factor_recovery_heatmap"), 8, 6, 900, 700, function() {
    par(mar = c(5, 5, 4, 2))
    plot_factor_recovery(cors, match_df,
                         sprintf("Factor Recovery Heatmap\nalpha=%.2f | K_eff=%d", alpha_opt, final_k$K_effective))
  })
  save_plot_pair(file.path(figure_dir, "ard_shrinkage"), 8, 5, 900, 600, function() {
    par(mar = c(5, 5, 4, 2))
    plot_ard_shrinkage(full_fit$history$factor_pve, match_df, final_k$active,
                       sprintf("ARD Shrinkage Trajectories\nalpha=%.2f | K_max=%d", alpha_opt, K_max))
  })
  save_plot_pair(file.path(figure_dir, "beta_recovery"), 8, 5.5, 900, 650, function() {
    par(mar = c(7, 5, 4, 2))
    plot_beta_recovery(beta_df,
                       sprintf("Beta Recovery with 95%% CIs\ncoverage=%.1f%%",
                               100 * mean(beta_df$ci_cover)))
  })
  save_plot_pair(file.path(figure_dir, "elbo_trace"), 8, 5, 900, 600, function() {
    par(mar = c(5, 5, 4, 2))
    plot_elbo_trace(elbo_df, "ELBO Convergence Trace")
  })
  save_plot_pair(file.path(figure_dir, "holdout_cindex"), 6.5, 4.5, 800, 500, function() {
    par(mar = c(6, 5, 4, 1))
    plot_holdout_cindex(holdout_df, "Held-Out C-index")
  })

  cat(sprintf("  Outputs written to %s\n", benchmark_root))

  invisible(list(
    summary = summary_df,
    cv = cv_res,
    full_fit = full_fit,
    holdout_fit = holdout_fit,
    factor_matches = match_df,
    beta = beta_df
  ))
}

if (sys.nframe() == 0) {
  run_ssbmf_benchmark()
}
