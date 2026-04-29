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
} else if (file.exists("code/fit_modular.R")) {
  # already at repo root (e.g. invoked via Rscript from project root)
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")   # running from results/benchmark_sim/
} else if (file.exists("../code/fit_modular.R")) {
  setwd("..")      # running from results/
}

cfg <- yaml::read_yaml("config/globals.yml")

source("code/train_test_split.R")
source("code/predict.R")
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
source("code/select_alpha_cv.R")
source("code/preprocess_desurv.R")

# Platform-specific log2+1 transform flag:
#   RNA-seq (raw counts)  → TRUE  (log2(x+1) is appropriate)
#   Proteomics / Microarray (already on log / normalized scale) → FALSE
PLATFORM_LOG_TRANSFORM <- c(
  TCGA_PAAD         = TRUE,
  CPTAC             = FALSE,
  Dijk              = TRUE,
  Moffitt_GEO_array = FALSE,
  PACA_AU_array     = FALSE,
  PACA_AU_seq       = TRUE,
  Puleo_array       = FALSE
)

# PDAC data root — same default as run_factor_modular_simulation.R
PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT", unset = path.expand(
  paste0("~/Library/CloudStorage/",
         "OneDrive-UniversityofNorthCarolinaatChapelHill/",
         "UNC Dissertation (Liu)/PDAC_data")
))

TRAINING_COHORTS  <- c("TCGA_PAAD", "CPTAC")
EXTERNAL_COHORTS  <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

PLATFORM_MAP <- c(
  TCGA_PAAD         = "RNA-seq",
  CPTAC             = "Proteomics",
  Dijk              = "RNA-seq",
  Moffitt_GEO_array = "Microarray",
  PACA_AU_array     = "Microarray",
  PACA_AU_seq       = "RNA-seq",
  Puleo_array       = "Microarray"
)

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

  beta_true <- as.numeric(cfg$synthetic$b_true)
  if (length(beta_true) != K_true) {
    stop(sprintf(
      "cfg$synthetic$b_true has length %d but K_true=%d. Update globals.yml.",
      length(beta_true), K_true
    ))
  }

  L_true <- matrix(rexp(n * K_true, rate = 1), n, K_true)
  signal_scale <- 0.25
  F_true <- matrix(0, p, K_true)
  for (k in seq_len(K_true)) {
    active <- sample.int(p, size = max(1L, round(0.05 * p)))
    F_true[active, k] <- signal_scale
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
  legend("bottomleft",
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

run_ssbmf_benchmark <- function(output_root = "results/benchmark_sim/outputs",
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
                                prior_beta = "point_normal",
                                lambda = 1.0,
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

  benchmark_root <- file.path(output_root, "synthetic", prior_beta)
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
    tol = tol,
    prior_beta = prior_beta
  )
  alpha_opt <- cv_res$alpha_opt

  holdout_fit <- fit_supervised_mf_modular(
    data$Y[train_idx, , drop = FALSE],
    data$time[train_idx],
    data$status[train_idx],
    K = K_max,
    alpha = alpha_opt,
    lambda = lambda,
    max_iter = max_iter,
    tol = tol,
    prior_beta = prior_beta,
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
    lambda = lambda,
    max_iter = max_iter,
    tol = tol,
    prior_beta = prior_beta,
    verbose = FALSE
  )

  cors <- safe_cor_matrix(data$L_true, full_fit$EL)
  match_df <- greedy_rectangular_match(cors)
  est_order <- c(match_df$est_factor[order(match_df$true_factor)],
                 setdiff(seq_len(K_max), match_df$est_factor))
  aligned_est <- full_fit$EBeta[est_order]
  aligned_sd  <- sqrt(pmax(full_fit$EBeta2[est_order] - full_fit$EBeta[est_order]^2, 0))
  aligned_true <- c(data$beta_true[match_df$true_factor[order(match_df$true_factor)]], rep(0, K_max - nrow(match_df)))
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
    prior_beta = prior_beta,
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
    par(mar = c(4, 5, 4, 1))
    bar_cols <- c("#1F77B4", "#FF7F0E")
    bp <- barplot(holdout_df$C_Index, names.arg = holdout_df$Method, col = bar_cols,
                  ylim = c(0, 1), main = "Held-Out C-index (80/20 split)",
                  ylab = "C-index", las = 1, cex.names = 0.9, bty = "n")
    abline(h = 0.5, col = "#666666", lty = 2)
    text(bp, holdout_df$C_Index + 0.03, sprintf("%.3f", holdout_df$C_Index),
         cex = 0.9, font = 2)
  })

  # --- RMSE trace ---
  save_plot_pair(file.path(figure_dir, "rmse_trace"), 8, 5, 900, 600, function() {
    par(mar = c(5, 5, 4, 2))
    plot(elbo_df$RMSE, type = "l", lwd = 2, col = "#1F77B4",
         xlab = "Iteration", ylab = "RMSE",
         main = "Genomics Reconstruction RMSE", bty = "n")
  })

  # --- PVE scree ---
  save_plot_pair(file.path(figure_dir, "pve_scree"), 7, 5, 800, 570, function() {
    par(mar = c(5, 5, 4, 2))
    pve_final <- as.numeric(full_fit$history$factor_pve[nrow(full_fit$history$factor_pve), ])
    ord <- order(pve_final, decreasing = TRUE)
    active_ord <- final_k$active[ord]
    factor_labels <- paste0("F", ord)
    bp <- barplot(pve_final[ord] * 100,
                  col = ifelse(active_ord, "#1F77B4", "#BDBDBD"),
                  names.arg = factor_labels,
                  xlab = "Factor (sorted by PVE)", ylab = "PVE (%)",
                  main = sprintf("PVE Scree  (K_eff = %d / %d)", final_k$K_effective, K_max),
                  las = 1, cex.names = 0.75, bty = "n")
    abline(h = 1, col = "#D62728", lty = 2)
    legend("topright", legend = c("Active (|b|>0.05 or PVE>1%)", "Pruned", "1% threshold"),
           fill = c("#1F77B4", "#BDBDBD", NA), border = NA,
           lty = c(NA, NA, 2), col = c(NA, NA, "#D62728"), bty = "n", cex = 0.75)
  })

  # --- Alpha CV curve ---
  cv_tbl_synth <- cv_res$cv_table
  save_plot_pair(file.path(figure_dir, "alpha_cv_curve"), 7, 5, 800, 550, function() {
    par(mar = c(5, 5, 4, 2))
    has_se <- !is.null(cv_tbl_synth$se_cindex) && any(!is.na(cv_tbl_synth$se_cindex))
    y_lim  <- if (has_se)
      range(c(cv_tbl_synth$mean_cindex - cv_tbl_synth$se_cindex,
              cv_tbl_synth$mean_cindex + cv_tbl_synth$se_cindex), na.rm = TRUE)
    else range(cv_tbl_synth$mean_cindex, na.rm = TRUE)
    plot(cv_tbl_synth$alpha, cv_tbl_synth$mean_cindex,
         type = "b", pch = 16, col = "#1F77B4", lwd = 2,
         xlab = expression(alpha), ylab = "Mean CV C-index",
         ylim = y_lim,
         main = sprintf("Alpha CV - Synthetic (n=%d, p=%d)\nalpha_opt=%.2f",
                        synthetic_n, synthetic_p, alpha_opt),
         bty = "n")
    if (has_se) {
      lo <- cv_tbl_synth$mean_cindex - cv_tbl_synth$se_cindex
      hi <- cv_tbl_synth$mean_cindex + cv_tbl_synth$se_cindex
      draw_idx <- which(is.finite(lo) & is.finite(hi) & (hi - lo) > 1e-10)
      if (length(draw_idx) > 0) {
        arrows(cv_tbl_synth$alpha[draw_idx], lo[draw_idx],
               cv_tbl_synth$alpha[draw_idx], hi[draw_idx],
               angle = 90, code = 3, length = 0.04, col = "#1F77B4", lwd = 1.2)
      }
    }
    abline(v = alpha_opt, col = "#D62728", lty = 2, lwd = 1.5)
    abline(h = 0.5, col = "#666666", lty = 3)
  })

  # --- GEP loading heatmap (top genes per matched factor) ---
  save_plot_pair(file.path(figure_dir, "gep_loading_heatmap"), 9, 6, 1000, 700, function() {
    layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
    par(mar = c(5, 6, 4, 1))
    n_match  <- nrow(match_df)
    ef_cols  <- est_order[seq_len(n_match)]
    ef_mat   <- full_fit$EF[, ef_cols, drop = FALSE]
    top_per  <- 20L
    top_idx  <- unique(unlist(lapply(seq_len(ncol(ef_mat)), function(k)
      order(abs(ef_mat[, k]), decreasing = TRUE)[seq_len(top_per)])))
    mat      <- t(ef_mat[top_idx, , drop = FALSE])
    pal      <- colorRampPalette(c("#08306B", "#F7FBFF", "#A50F15"))(100)
    image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat),
          col = pal, axes = FALSE,
          xlab = sprintf("Gene index (top %d per factor)", top_per),
          ylab = "",
          main = "GEP Loading Heatmap (matched factors)")
    axis(1, at = c(1, ncol(mat)), labels = c("1", as.character(ncol(mat))),
         cex.axis = 0.8)
    axis(2, at = seq_len(nrow(mat)),
         labels = paste0("GEP", match_df$true_factor[order(match_df$true_factor)]),
         las = 1, cex.axis = 0.85)
    box()
    # Color scale panel
    par(mar = c(5, 1, 4, 3))
    scale_vals <- seq(min(mat), max(mat), length.out = 100)
    image(1, scale_vals, matrix(scale_vals, nrow = 1),
          col = pal, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = pretty(scale_vals, n = 4),
         labels = formatC(pretty(scale_vals, n = 4), format = "f", digits = 2),
         las = 1, cex.axis = 0.8)
    mtext("Loading", side = 4, line = 2.5, cex = 0.8)
    box()
  })

  # --- KM 3-group on synthetic training set ---
  lp_train <- as.vector(full_fit$EL %*% full_fit$EBeta)
  save_plot_pair(file.path(figure_dir, "km_3group_training"), 7, 5, 800, 560, function() {
    par(mar = c(5, 5, 4, 1))
    plot_km_2group(lp_train, data$time, data$status,
                   title = sprintf("KM Stratification - Synthetic Training\n(n=%d, alpha=%.2f)",
                                   synthetic_n, alpha_opt))
  })

  # --- Tau precision distribution ---
  save_plot_pair(file.path(figure_dir, "tau_distribution"), 7, 5, 800, 550, function() {
    par(mar = c(5, 5, 4, 2))
    tau_est <- if (!is.null(full_fit$Tau)) full_fit$Tau else rep(NA_real_, synthetic_p)
    if (any(!is.na(tau_est))) {
      hist(log10(tau_est + 1e-8), breaks = 40,
           col = "#9ECAE1", border = "white",
           xlab = expression(log[10](hat(tau)[j])),
           main = "Gene Precision Distribution  (estimated vs. true)",
           bty = "n")
      abline(v = mean(log10(data$tau_true + 1e-8)),
             col = "#D62728", lwd = 2, lty = 2)
      legend("topright", legend = "True tau mean",
             col = "#D62728", lty = 2, lwd = 2, bty = "n")
    } else {
      plot.new()
      title("Tau not available in this fit")
    }
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

# ==============================================================================
# Real-data helpers
# ==============================================================================

#' Load a PDAC cohort without any preprocessing.
#'
#' Replicates the symlink trick from run_factor_modular_simulation.R::load_real_data()
#' but returns the raw expression matrix so that preprocess_desurv_cohort() can
#' be applied separately.
#'
#' Gene names are deduplicated (first occurrence kept) to avoid intersection
#' failures caused by duplicate SYMBOL entries in some platforms.
#'
#' @param dataset_name string: cohort name (must exist in PDAC data root)
#' @param pdac_root    string: path to the PDAC_data directory
#' @return list(Y, gene_names, time, status, n, p, dataset_name)
load_pdac_raw <- function(dataset_name, pdac_root) {
  if (!dir.exists(pdac_root))
    stop(sprintf("PDAC data root not found: %s\nSet PDAC_DATA_ROOT env var.", pdac_root))

  tmp_wd    <- tempfile("pdac_wd_")
  dir.create(tmp_wd, showWarnings = FALSE)
  data_link <- file.path(tmp_wd, "data")
  file.symlink(pdac_root, data_link)

  old_wd <- getwd()
  on.exit({
    setwd(old_wd)
    unlink(data_link)
  }, add = TRUE)

  setwd(tmp_wd)
  source(file.path(pdac_root, "load_data_internal.R"), local = TRUE)
  result <- load_data_internal(dataset_name)
  setwd(old_wd)

  keeps <- which(result$sampInfo$keep == 1)
  if (length(keeps) == 0)
    stop(sprintf("No valid samples for dataset '%s' after keep filter.", dataset_name))

  # Transpose: genes × samples → patients × genes (n × p)
  Y <- t(result$ex[, keeps])

  # Extract gene names
  fi <- result$featInfo
  if (is.data.frame(fi) && "SYMBOL" %in% names(fi)) {
    gene_names <- fi$SYMBOL
  } else if (is.character(fi)) {
    gene_names <- fi
  } else {
    gene_names <- rownames(result$ex)
  }
  if (length(gene_names) != ncol(Y))
    gene_names <- paste0("Gene", seq_len(ncol(Y)))

  # Deduplicate gene names: keep first occurrence to avoid intersection failures
  dup_mask  <- duplicated(gene_names)
  if (any(dup_mask)) {
    keep_cols  <- !dup_mask
    Y          <- Y[, keep_cols, drop = FALSE]
    gene_names <- gene_names[keep_cols]
  }

  colnames(Y) <- gene_names
  rownames(Y) <- NULL

  time   <- result$sampInfo$time[keeps]
  status <- as.integer(result$sampInfo$event[keeps])

  list(
    Y = Y, gene_names = gene_names,
    time = time, status = status,
    n = nrow(Y), p = ncol(Y),
    dataset_name = dataset_name
  )
}

#' KM survival curve for a 3-group stratification.
#'
#' Splits patients into tertiles (Low / Mid / High) by linear predictor lp.
#' Plots Kaplan–Meier curves and reports log-rank p-value.
#'
#' @param lp      numeric vector: linear predictor (L_test %*% beta)
#' @param time    numeric vector: follow-up time
#' @param status  integer vector: event indicator
#' @param title   plot title
plot_km_2group <- function(lp, time, status, title = "") {
  if (diff(range(lp)) < 1e-10) {
    plot.new()
    title(sprintf("%s\n(constant risk score — no stratification possible)", title))
    return(invisible(NULL))
  }
  med_cut <- quantile(lp, probs = 0.5)
  grp <- cut(lp, breaks = c(-Inf, med_cut, Inf),
             labels = c("Low risk", "High risk"), include.lowest = TRUE)

  km_fit  <- survfit(Surv(time, status) ~ grp)
  sd_test <- survdiff(Surv(time, status) ~ grp)
  p_val   <- 1 - pchisq(sd_test$chisq, df = 1)

  km_cols <- c("#2166AC", "#D6604D")
  plot(km_fit, col = km_cols, lwd = 2, lty = 1,
       xlab = "Time (days)", ylab = "Survival probability",
       main = sprintf("%s\nLog-rank p = %.4f", title, p_val),
       bty = "n", mark.time = TRUE)
  legend("topright",
         legend = c("Low risk", "High risk"),
         col = km_cols, lwd = 2, bty = "n", cex = 0.85)
  invisible(list(km_fit = km_fit, p_val = p_val, groups = grp))
}

#' Run the PDAC cross-cohort benchmark (Phase 3B).
#'
#' Training: TCGA-PAAD + CPTAC (DeSurv-preprocessed, intersected).
#' External validation: Dijk, Moffitt, PACA-AU array/seq, Puleo.
#'
#' @param output_root  directory for outputs (tables/ + figures/ subdirs)
#' @param pdac_root    path to PDAC_data directory
#' @param alpha_grid   candidate alpha values for CV
#' @param n_folds      number of CV folds
#' @param K_max        number of factors
#' @param max_iter     max CAVI iterations per fit
#' @param tol          ELBO convergence tolerance
#' @param top_n        number of top-variable genes per cohort (DeSurv spec: 2000)
#' @return invisibly: list with all results
run_real_data_benchmark <- function(
    training_mode          = "merged",
    prior_beta             = "point_normal",
    preprocessing_version  = "v1",
    output_root            = file.path("results/benchmark_sim/outputs/real_data",
                                       training_mode,
                                       if (preprocessing_version == "v1") prior_beta
                                       else paste0("v2_", prior_beta)),
    pdac_root              = PDAC_DATA_ROOT,
    alpha_grid             = c(0.1, 0.3, 0.5, 0.7, 0.9),
    n_folds                = 5,
    K_max                  = 10,
    max_iter               = 300,
    tol                    = 1e-5,
    lambda                 = 1.0,
    top_n                  = 2000,
    gene_intersection_only = FALSE) {

  preprocessing_version <- match.arg(preprocessing_version, c("v1", "v2"))

  if (!dir.exists(pdac_root)) {
    message(sprintf(
      "PDAC data not found at '%s' — skipping real-data benchmark.\n",
      pdac_root))
    return(invisible(NULL))
  }

  # Determine training cohorts for this mode
  active_train_cohorts <- switch(training_mode,
    merged    = c("TCGA_PAAD", "CPTAC"),
    tcga_only  = "TCGA_PAAD",
    cptac_only = "CPTAC",
    stop(sprintf("Unknown training_mode: '%s'. Use 'merged', 'tcga_only', or 'cptac_only'.",
                 training_mode))
  )

  table_dir  <- file.path(output_root, "tables")
  figure_dir <- file.path(output_root, "figures")
  ensure_dir(table_dir)
  ensure_dir(figure_dir)

  # ============================================================
  # 1. Load + preprocess training cohorts
  # ============================================================
  cat(sprintf("=== Phase 3B: PDAC Cross-Cohort Benchmark [mode: %s] ===\n", training_mode))
  cat("  Loading training cohorts:", paste(active_train_cohorts, collapse = " + "), "\n")

  train_raw <- lapply(setNames(active_train_cohorts, active_train_cohorts), function(ds) {
    cat(sprintf("    Loading %s ...\n", ds))
    load_pdac_raw(ds, pdac_root)
  })

  # ============================================================
  # v1 vs v2 preprocessing branch
  #   v1 (default): per-cohort log2 → top-N → rank, then intersect
  #   v2 (reordered): intersect raw → log2 → QN → merged variance →
  #                   top-N → rank  (preprocess_merged_cohorts)
  # Single-cohort modes always use the v1 per-cohort path regardless of
  # preprocessing_version (QN across a single cohort is a no-op and the
  # bug being fixed only affects the merged intersection step).
  # ============================================================

  if (preprocessing_version == "v2" && length(active_train_cohorts) > 1) {
    # v2 merged path ------------------------------------------------
    cat(sprintf("  Preprocessing version: v2 (intersect-first + QN)\n"))
    log_flags <- PLATFORM_LOG_TRANSFORM[active_train_cohorts]
    merged_v2 <- preprocess_merged_cohorts(
      cohort_raw_list    = train_raw,
      log_transform_flags = log_flags,
      top_n              = top_n
    )
    Y_train          <- merged_v2$Y
    training_gene_names <- merged_v2$gene_names
    n_train_genes    <- merged_v2$p
    time_train       <- unlist(lapply(active_train_cohorts, function(ds) train_raw[[ds]]$time))
    status_train     <- unlist(lapply(active_train_cohorts, function(ds) train_raw[[ds]]$status))
    cohort_labels    <- merged_v2$dataset_labels

    # expose intersection count for the checkpoint below
    train_preproc    <- NULL  # not used in v2 merged path
    n_raw_intersect  <- merged_v2$n_raw_intersect

    if (gene_intersection_only) {
      cat(sprintf(
        "\n=== INTERSECTION CHECKPOINT (v2) ===\n  Raw gene intersection: %d genes\n  After top-%d selection: %d genes\n=== Stopping here (gene_intersection_only=TRUE). ===\n",
        n_raw_intersect, top_n, n_train_genes
      ))
      return(invisible(list(n_raw_intersect = n_raw_intersect,
                            n_genes_selected = n_train_genes)))
    }

  } else {
    # v1 path (per-cohort preprocess then intersect) ----------------
    if (preprocessing_version == "v2")
      cat("  preprocessing_version=v2 ignored for single-cohort mode; using v1.\n")

    train_preproc <- lapply(active_train_cohorts, function(ds) {
      raw <- train_raw[[ds]]
      cat(sprintf("    Preprocessing %s (n=%d, p_raw=%d) ...\n", ds, raw$n, raw$p))
      preprocess_desurv_cohort(
        Y             = raw$Y,
        gene_names    = raw$gene_names,
        top_n         = top_n,
        log_transform = PLATFORM_LOG_TRANSFORM[[ds]],
        cohort_name   = ds
      )
    })
    names(train_preproc) <- active_train_cohorts

    if (length(active_train_cohorts) > 1) {
      train_intersected   <- intersect_preprocessed_cohorts(train_preproc, reference = 1)
      names(train_intersected) <- active_train_cohorts
      n_train_genes       <- train_intersected[[1]]$p
      training_gene_names <- train_intersected[[1]]$gene_names
      cat(sprintf("  Training gene intersection: %d genes\n", n_train_genes))
    } else {
      ds_single           <- active_train_cohorts[1]
      n_train_genes       <- train_preproc[[ds_single]]$p
      training_gene_names <- train_preproc[[ds_single]]$gene_names
      cat(sprintf("  Single-cohort training: %d genes\n", n_train_genes))
    }

    if (gene_intersection_only) {
      if (training_mode != "merged") {
        cat(sprintf("  gene_intersection_only only applies to 'merged' mode (current: '%s').\n",
                    training_mode))
        return(invisible(NULL))
      }
      cat(sprintf(
        "\n=== INTERSECTION CHECKPOINT ===\n  TCGA_PAAD: n=%d, p_raw=%d, p_after_top%d=%d\n  CPTAC:     n=%d, p_raw=%d, p_after_top%d=%d\n  Shared genes after intersection: %d\n  Event rates: TCGA=%.1f%%, CPTAC=%.1f%%\n=== Stopping here. Set gene_intersection_only=FALSE to proceed. ===\n",
        train_raw[["TCGA_PAAD"]]$n, train_raw[["TCGA_PAAD"]]$p, top_n, train_preproc[["TCGA_PAAD"]]$p,
        train_raw[["CPTAC"]]$n,     train_raw[["CPTAC"]]$p,     top_n, train_preproc[["CPTAC"]]$p,
        n_train_genes,
        100 * mean(train_raw[["TCGA_PAAD"]]$status),
        100 * mean(train_raw[["CPTAC"]]$status)
      ))
      return(invisible(list(
        n_intersect     = n_train_genes,
        n_TCGA          = train_raw[["TCGA_PAAD"]]$n,
        n_CPTAC         = train_raw[["CPTAC"]]$n,
        p_TCGA_preproc  = train_preproc[["TCGA_PAAD"]]$p,
        p_CPTAC_preproc = train_preproc[["CPTAC"]]$p
      )))
    }

    # Assemble final training matrices
    # cohort_labels persisted to final_model.rds for cohort-stratified heatmaps.
    if (length(active_train_cohorts) > 1) {
      merged_train  <- merge_preprocessed_cohorts(train_intersected,
                                                   dataset_labels = active_train_cohorts)
      Y_train       <- merged_train$Y
      time_train    <- unlist(lapply(active_train_cohorts, function(ds) train_raw[[ds]]$time))
      status_train  <- unlist(lapply(active_train_cohorts, function(ds) train_raw[[ds]]$status))
      cohort_labels <- merged_train$dataset_labels
    } else {
      ds_single     <- active_train_cohorts[1]
      Y_train       <- train_preproc[[ds_single]]$Y
      time_train    <- train_raw[[ds_single]]$time
      status_train  <- train_raw[[ds_single]]$status
      cohort_labels <- factor(rep(ds_single, nrow(Y_train)), levels = ds_single)
    }
  }
  n_train <- nrow(Y_train)
  cat(sprintf("  Training set [%s]: n=%d, p=%d, event_rate=%.1f%%\n",
              training_mode, n_train, ncol(Y_train), 100 * mean(status_train)))

  # Checkpoint 1: save training set summary
  # p_preproc: for v2, all cohorts share the merged p (no per-cohort preproc object).
  p_preproc_vals <- if (!is.null(train_preproc)) {
    sapply(active_train_cohorts, function(ds) train_preproc[[ds]]$p)
  } else {
    rep(n_train_genes, length(active_train_cohorts))
  }
  write.csv(
    data.frame(
      Cohort        = active_train_cohorts,
      n             = sapply(active_train_cohorts, function(ds) train_raw[[ds]]$n),
      p_raw         = sapply(active_train_cohorts, function(ds) train_raw[[ds]]$p),
      p_preproc     = p_preproc_vals,
      log_transform = PLATFORM_LOG_TRANSFORM[active_train_cohorts],
      platform      = PLATFORM_MAP[active_train_cohorts],
      event_rate    = sapply(active_train_cohorts, function(ds) mean(train_raw[[ds]]$status))
    ),
    file.path(table_dir, "training_cohort_summary.csv"),
    row.names = FALSE
  )
  cat("  CHECKPOINT 1: Training cohort summary saved.\n")

  # ============================================================
  # 2. Alpha cross-validation on merged training set
  # ============================================================
  cat(sprintf("  Running %d-fold CV over alpha grid: %s\n",
              n_folds, paste(alpha_grid, collapse = ", ")))
  cat("  (This may take several minutes — one fit per alpha × fold)\n")

  cv_res <- select_alpha_cv(
    Y_train, time_train, status_train,
    alpha_grid = alpha_grid,
    n_folds    = n_folds,
    K_max      = K_max,
    use_1se    = TRUE,
    seed       = 42,
    max_iter   = max_iter,
    tol        = tol,
    prior_beta = prior_beta,
    verbose    = FALSE
  )
  alpha_opt <- cv_res$alpha_opt
  cat(sprintf("  Alpha CV complete: alpha_opt = %.2f (rule: %s)\n",
              alpha_opt, cv_res$selection_rule))

  write.csv(cv_res$cv_table,
            file.path(table_dir, "alpha_cv_table.csv"), row.names = FALSE)
  write.csv(cv_res$fold_results,
            file.path(table_dir, "alpha_cv_fold_results.csv"), row.names = FALSE)

  # Alpha CV curve figure
  cv_tbl <- cv_res$cv_table
  save_plot_pair(file.path(figure_dir, "alpha_cv_curve"), 7, 5, 800, 550, function() {
    par(mar = c(5, 5, 4, 2))
    y_lim <- range(c(cv_tbl$mean_cindex - cv_tbl$se_cindex,
                     cv_tbl$mean_cindex + cv_tbl$se_cindex), na.rm = TRUE)
    plot(cv_tbl$alpha, cv_tbl$mean_cindex,
         type = "b", pch = 16, col = "#1F77B4", lwd = 2,
         xlab = expression(alpha), ylab = "Mean C-index (CV)",
         ylim = y_lim,
         main = sprintf("Alpha CV [%s, %d-fold]\nalpha_opt=%.2f",
                        training_mode, n_folds, alpha_opt),
         bty = "n")
    arrows(cv_tbl$alpha,
           cv_tbl$mean_cindex - cv_tbl$se_cindex,
           cv_tbl$alpha,
           cv_tbl$mean_cindex + cv_tbl$se_cindex,
           angle = 90, code = 3, length = 0.04, col = "#1F77B4", lwd = 1.2)
    abline(v = alpha_opt, col = "#D62728", lty = 2, lwd = 1.5)
    abline(h = 0.5, col = "#666666", lty = 3)
    legend("bottomright",
           legend = c(sprintf("alpha_opt = %.2f", alpha_opt), "C-index = 0.5"),
           col = c("#D62728", "#666666"), lty = c(2, 3), lwd = 1.5, bty = "n")
  })
  cat("  CHECKPOINT 2: Alpha CV complete + curve saved.\n")

  # ============================================================
  # 3. Final model fit on full training set
  # ============================================================
  cat(sprintf("  Fitting final model: K=%d, alpha=%.2f ...\n", K_max, alpha_opt))
  final_fit <- fit_supervised_mf_modular(
    Y_train, time_train, status_train,
    K          = K_max,
    alpha      = alpha_opt,
    lambda     = lambda,
    max_iter   = max_iter,
    tol        = tol,
    prior_beta = prior_beta,
    verbose    = FALSE
  )
  cat(sprintf("  Final fit: converged=%s, iter=%d, ELBO=%.4f\n",
              final_fit$history$converged, final_fit$history$n_iter,
              tail(final_fit$history$elbo_full, 1)))

  # Effective K (auto-prune thresholds)
  final_pve   <- as.numeric(tail(final_fit$history$factor_pve, 1))
  active_mask <- (abs(final_fit$EBeta) > 0.05) | (final_pve > 0.01)
  K_eff       <- sum(active_mask)

  # Training beta table
  beta_train_df <- data.frame(
    Factor    = seq_len(K_max),
    EBeta     = round(final_fit$EBeta, 4),
    EBeta_SD  = round(sqrt(pmax(final_fit$EBeta2 - final_fit$EBeta^2, 0)), 4),
    PVE_pct   = round(as.numeric(final_pve) * 100, 2),
    Active    = active_mask
  )
  write.csv(beta_train_df, file.path(table_dir, "training_beta_summary.csv"),
            row.names = FALSE)

  # ELBO trace for training fit
  elbo_train_df <- data.frame(
    Iteration  = seq_along(final_fit$history$elbo_full),
    ELBO_Full  = final_fit$history$elbo_full,
    ELBO_Proxy = final_fit$history$elbo_proxy
  )
  write.csv(elbo_train_df, file.path(table_dir, "training_elbo_trace.csv"),
            row.names = FALSE)

  # Save fitted factor matrices for downstream use (PH diagnostics, gene enrichment,
  # cohort-stratified loading heatmap). EL and cohort_labels added so Phase 1
  # diagnostic plots can be regenerated from the RDS without re-fitting.
  saveRDS(
    list(EF = final_fit$EF, EBeta = final_fit$EBeta,
         EL = final_fit$EL, cohort_labels = cohort_labels,
         time_train = time_train, status_train = status_train,
         alpha_opt = alpha_opt, training_gene_names = training_gene_names,
         training_mode = training_mode, prior_beta = prior_beta),
    file.path(table_dir, "final_model.rds")
  )
  cat(sprintf("  CHECKPOINT 3: Final model fit complete. K_eff=%d\n", K_eff))

  # ============================================================
  # 4. External cohort validation
  # ============================================================
  cat("  Running external cohort projections ...\n")

  external_results <- lapply(EXTERNAL_COHORTS, function(ds) {
    cat(sprintf("    Projecting %s ...\n", ds))
    tryCatch({
      raw <- load_pdac_raw(ds, pdac_root)
      preproc <- preprocess_desurv_cohort(
        Y           = raw$Y,
        gene_names  = raw$gene_names,
        top_n       = top_n,
        log_transform = PLATFORM_LOG_TRANSFORM[[ds]],
        cohort_name = ds
      )

      # Subset to training gene intersection (genes in training set)
      common_idx <- match(training_gene_names, preproc$gene_names)
      missing    <- sum(is.na(common_idx))
      if (missing > 0)
        warning(sprintf("%s: %d training genes not found in external cohort (will be zero-filled).",
                        ds, missing))

      # Align external Y to training gene order (fill missing with column mean = 0 after rank-transform)
      Y_ext <- matrix(0.0, nrow = raw$n, ncol = length(training_gene_names))
      present <- !is.na(common_idx)
      Y_ext[, present] <- preproc$Y[, common_idx[present], drop = FALSE]

      pred <- predict_supervised_mf(Y_ext, final_fit$EF, final_fit$EBeta)
      c_idx <- as.numeric(concordance(
        Surv(raw$time, raw$status) ~ I(-pred$risk_scores)
      )$concordance)

      # Proportional hazards diagnostic (Grambsch-Therneau test on risk score)
      ph_result <- tryCatch({
        lp_vec <- as.vector(pred$risk_scores)
        coxfit  <- coxph(Surv(raw$time, raw$status) ~ lp_vec, ties = "efron")
        zph     <- cox.zph(coxfit)
        tbl     <- zph$table
        list(
          ph_chisq = round(tbl["lp_vec", "chisq"], 3),
          ph_df    = as.integer(tbl["lp_vec", "df"]),
          ph_p     = round(tbl["lp_vec", "p"], 4)
        )
      }, error = function(e) list(ph_chisq = NA_real_, ph_df = NA_integer_, ph_p = NA_real_))

      km_res <- NULL
      if (sum(raw$status) >= 10) {
        km_stub <- file.path(figure_dir, sprintf("km_3group_%s", ds))
        save_plot_pair(km_stub, 7, 5, 800, 560, function() {
          par(mar = c(5, 5, 4, 1))
          plot_km_2group(
            lp     = pred$risk_scores,
            time   = raw$time,
            status = raw$status,
            title  = sprintf("%s (n=%d, %s)", ds, raw$n, PLATFORM_MAP[[ds]])
          )
        })
        km_res <- tryCatch({
          tertile_cut <- quantile(pred$risk_scores, probs = c(1/3, 2/3))
          grp <- cut(pred$risk_scores, breaks = c(-Inf, tertile_cut, Inf),
                     labels = c("Low", "Mid", "High"), include.lowest = TRUE)
          sd_test <- survdiff(Surv(raw$time, raw$status) ~ grp)
          1 - pchisq(sd_test$chisq, df = 2)
        }, error = function(e) NA_real_)
      }

      list(
        dataset      = ds,
        n            = raw$n,
        events       = sum(raw$status),
        p_intersect  = sum(!is.na(common_idx)),
        platform     = PLATFORM_MAP[[ds]],
        c_index      = round(c_idx, 4),
        km_logrank_p = km_res,
        event_rate   = mean(raw$status),
        ph_chisq     = ph_result$ph_chisq,
        ph_df        = ph_result$ph_df,
        ph_p         = ph_result$ph_p,
        status       = "ok"
      )
    }, error = function(e) {
      warning(sprintf("External cohort %s failed: %s", ds, conditionMessage(e)))
      list(dataset = ds, status = "error", error_msg = conditionMessage(e),
           c_index = NA_real_, km_logrank_p = NA_real_,
           ph_chisq = NA_real_, ph_df = NA_integer_, ph_p = NA_real_)
    })
  })
  names(external_results) <- EXTERNAL_COHORTS

  # External results table
  ext_df <- do.call(rbind, lapply(external_results, function(r) {
    data.frame(
      Cohort      = r$dataset,
      n           = if (is.null(r$n)) NA_integer_ else r$n,
      Platform    = if (is.null(r$platform)) NA_character_ else r$platform,
      p_intersect = if (is.null(r$p_intersect)) NA_integer_ else r$p_intersect,
      C_index     = r$c_index,
      KM_logrank_p = if (is.null(r$km_logrank_p)) NA_real_ else r$km_logrank_p,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(ext_df, file.path(table_dir, "external_cindex_table.csv"), row.names = FALSE)

  # PH diagnostics table
  ph_diag_df <- do.call(rbind, lapply(external_results, function(r) {
    ph_p_val <- if (is.null(r$ph_p)) NA_real_ else r$ph_p
    data.frame(
      Cohort   = r$dataset,
      n        = if (is.null(r$n))       NA_integer_ else r$n,
      Events   = if (is.null(r$events))  NA_integer_ else r$events,
      Platform = if (is.null(r$platform)) NA_character_ else r$platform,
      C_index  = if (is.null(r$c_index)) NA_real_    else r$c_index,
      PH_chisq = if (is.null(r$ph_chisq)) NA_real_   else r$ph_chisq,
      PH_df    = if (is.null(r$ph_df))    NA_integer_ else r$ph_df,
      PH_p     = ph_p_val,
      PH_flag  = if (is.na(ph_p_val)) "ERROR"
                 else if (ph_p_val < 0.05) "FLAG (p<0.05)" else "PASS",
      stringsAsFactors = FALSE
    )
  }))
  write.csv(ph_diag_df, file.path(table_dir, "ph_diagnostics_table.csv"), row.names = FALSE)

  # External C-index bar chart
  ext_ok <- ext_df[!is.na(ext_df$C_index), ]
  if (nrow(ext_ok) > 0) {
    save_plot_pair(file.path(figure_dir, "external_cindex_barchart"), 8, 5, 900, 570, function() {
      par(mar = c(9, 5, 4, 2))
      bp <- barplot(ext_ok$C_index,
                    names.arg = ext_ok$Cohort,
                    col = "#4292C6", ylim = c(0, 1),
                    main = "External cohort C-index (SSBMF)",
                    ylab = "C-index", las = 2, cex.names = 0.75, bty = "n")
      abline(h = 0.5, col = "#666666", lty = 2)
      # DeSurv reference band (~0.60–0.65 per paper Fig 2B)
      rect(bp[1] - 0.6, 0.60, bp[nrow(bp)] + 0.6, 0.65,
           col = adjustcolor("#D62728", 0.15), border = NA)
      legend("topright",
             legend = c("SSBMF", "DeSurv CV range (0.60–0.65)"),
             fill = c("#4292C6", adjustcolor("#D62728", 0.15)),
             border = NA, bty = "n", cex = 0.85)
    })
  }
  cat("  CHECKPOINT 4: External cohort projections complete.\n")

  # ============================================================
  # 5. DeSurv comparison table (reading-based)
  # ============================================================
  desurv_compare <- data.frame(
    Metric           = c("Optimal K", "Optimal alpha", "Training strategy",
                         "CV C-index", "Factor 1 direction", "Factor 2 direction",
                         "Factor 3 direction", "External KM log-rank"),
    DeSurv_published = c("3", "0.7",
                         "TCGA+CPTAC (n=273)",
                         "~0.60–0.65 (Fig 2B)",
                         "Protective (HR=0.37, immune/iCAF)",
                         "Non-prognostic (HR=0.99, exocrine)",
                         "Risky (HR=1.43, basal-like)",
                         "p < 0.0001"),
    SSBMF_ours       = c(
      as.character(K_eff),
      sprintf("%.2f (1SE rule)", alpha_opt),
      sprintf("%s (n=%d)", training_mode, n_train),
      sprintf("%.4f (alpha_opt fold mean)",
              cv_tbl$mean_cindex[cv_tbl$alpha == alpha_opt]),
      sprintf("Factor w/ most negative beta: %.3f",
              min(final_fit$EBeta)),
      sprintf("Factor closest to 0: %.4f",
              final_fit$EBeta[which.min(abs(final_fit$EBeta))]),
      sprintf("Factor w/ most positive beta: %.3f",
              max(final_fit$EBeta)),
      if (all(is.na(ext_df$KM_logrank_p))) "Not computed"
      else sprintf("%.4f (median across cohorts)",
                   median(ext_df$KM_logrank_p, na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(desurv_compare,
            file.path(table_dir, "desurv_comparison_table.csv"), row.names = FALSE)

  # Overall summary row
  summary_df <- data.frame(
    training_mode   = training_mode,
    prior_beta      = prior_beta,
    n_train         = n_train,
    p_genes         = ncol(Y_train),
    alpha_opt       = alpha_opt,
    cv_rule         = cv_res$selection_rule,
    K_max           = K_max,
    K_eff           = K_eff,
    converged       = final_fit$history$converged,
    n_iter          = final_fit$history$n_iter,
    final_elbo      = round(tail(final_fit$history$elbo_full, 1), 4),
    median_ext_cindex = round(median(ext_df$C_index, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
  write.csv(summary_df, file.path(table_dir, "realdata_benchmark_summary.csv"),
            row.names = FALSE)

  cat(sprintf("  Outputs written to %s\n", output_root))
  cat("=== Phase 3B complete ===\n")

  invisible(list(
    Y_train          = Y_train,
    cv               = cv_res,
    alpha_opt        = alpha_opt,
    final_fit        = final_fit,
    K_eff            = K_eff,
    external_results = external_results,
    desurv_compare   = desurv_compare
  ))
}

# ==============================================================================
# Entry point
# ==============================================================================

if (sys.nframe() == 0) {
  PRIOR_GRID <- c("point_normal", "point_laplace")

  # --- Synthetic benchmark: both priors ---
  for (pb in PRIOR_GRID) {
    cat(sprintf("\n\n===== Synthetic benchmark: prior_beta = %s =====\n\n", pb))
    run_ssbmf_benchmark(
      prior_beta = pb,
      max_iter   = cfg$cavi$max_iter,
      tol        = cfg$cavi$tol,
      alpha_grid = cfg$cavi$alpha_grid,
      K_max      = cfg$cavi$k_max,
      seed       = cfg$synthetic$seed,
      lambda     = cfg$cavi$lambda
    )
  }

  # --- Real-data benchmark: tcga_only + cptac_only × both priors ---
  # merged is kept for the multi-modal failure documentation (point_normal only)
  run_real_data_benchmark(
    training_mode = "merged",
    prior_beta    = "point_normal",
    max_iter      = cfg$cavi$max_iter,
    tol           = cfg$cavi$tol,
    alpha_grid    = cfg$cavi$alpha_grid,
    K_max         = cfg$cavi$k_max,
    lambda        = cfg$cavi$lambda
  )

  for (mode in c("tcga_only", "cptac_only")) {
    for (pb in PRIOR_GRID) {
      cat(sprintf("\n\n===== Real-data benchmark: mode=%s prior=%s =====\n\n", mode, pb))
      run_real_data_benchmark(
        training_mode = mode,
        prior_beta    = pb,
        max_iter      = cfg$cavi$max_iter,
        tol           = cfg$cavi$tol,
        alpha_grid    = cfg$cavi$alpha_grid,
        K_max         = cfg$cavi$k_max,
        lambda        = cfg$cavi$lambda
      )
    }
  }

  # --- Cross-mode × cross-prior comparison table ---
  all_combos <- rbind(
    data.frame(mode = "merged",    pb = "point_normal", stringsAsFactors = FALSE),
    expand.grid(mode = c("tcga_only", "cptac_only"),
                pb   = PRIOR_GRID,
                stringsAsFactors = FALSE)
  )
  cmp_rows <- lapply(seq_len(nrow(all_combos)), function(i) {
    f <- file.path("results/benchmark_sim/outputs/real_data",
                   all_combos$mode[i], all_combos$pb[i],
                   "tables", "realdata_benchmark_summary.csv")
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
  })
  cmp_rows <- do.call(rbind, Filter(Negate(is.null), cmp_rows))
  if (!is.null(cmp_rows) && nrow(cmp_rows) > 0) {
    write.csv(cmp_rows,
              "results/benchmark_sim/outputs/real_data/training_mode_prior_comparison.csv",
              row.names = FALSE)
    cat("\nCross-mode × prior comparison table written.\n")
  }
}
