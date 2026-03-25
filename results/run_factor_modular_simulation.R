# ==============================================================================
# Factor-Wise Modular Simulation — Export All Results
#
# Uses fit_supervised_mf_modular() from code/fit_modular.R, which implements
# V3 Algorithm 1: factor-wise Gauss-Seidel CAVI.
#
# Factor-wise updates (_k variants): for each k, L_k -> F_k -> beta_k.
# This is the canonical coordinate ascent order: F_k sees the updated L_k
# immediately, rather than waiting until the next outer iteration (block-wise).
#
# Compare with results/run_modular_simulation.R (DEPRECATED — block-wise):
#   Block-wise: all-L -> all-F -> all-beta -> tau (one sweep per iteration)
#   Factor-wise: for k in 1:K { L_k, F_k, beta_k } (Gauss-Seidel per V3)
#
# Outputs:
#   results/tables/factor_modular_sim/  -- 7 CSV tables
#   results/figures/factor_modular_sim/ -- 8 figures (PDF + PNG)
#
# Run from repo root:
#   Rscript results/run_factor_modular_simulation.R
# ==============================================================================

# Set working directory to repo root (portable: works locally and on Longleaf)
if (Sys.getenv("REPO_ROOT") != "") {
  setwd(Sys.getenv("REPO_ROOT"))
} else if (file.exists("code/update_L.R")) {
  # Already at repo root (e.g., local RStudio)
} else if (file.exists("../code/update_L.R")) {
  setwd("..")
} else {
  stop("Cannot find repo root. Run from project root or set REPO_ROOT env var.")
}

# ==============================================================================
# Source fit_modular.R via tryCatch
#
# fit_modular.R has a runner block at the bottom that stopifnot()s when
# DATA_MODE="real" and real_Y is NULL.  We wrap the source() in tryCatch so
# the stopifnot error is silently caught.  By the time the error fires,
# fit_supervised_mf_modular() and all its dependencies (update_L.R,
# update_F.R, update_beta.R, update_tau.R, library(survival), library(ebnm),
# calc_cox_taylor) are already defined in the global environment.
# ==============================================================================

suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
# Now fit_supervised_mf_modular() is available, along with all its dependencies
# (update_L.R, update_F.R, update_beta.R, update_tau.R, calc_cox_taylor, library(survival), library(ebnm))

# ==============================================================================
# Create Output Directories
# ==============================================================================

dir.create("results/tables/factor_modular_sim",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/factor_modular_sim", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Analytics Helpers
# (Functionally identical to helpers in results/run_modular_simulation.R)
# ==============================================================================

# C-index comparison: Supervised loadings vs. top-5 PCA components
get_cindex_comparison <- function(EL, data) {
  pca_y  <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  fit_l  <- coxph(Surv(data$time, data$status) ~ EL)
  list(
    c_original = round(summary(fit_pc)$concordance[1], 3),
    c_latent   = round(summary(fit_l)$concordance[1], 3)
  )
}

# Top n_top influential features per factor by |weight|
get_top_features <- function(EF, n_top = 10) {
  lapply(1:ncol(EF), function(k) {
    weights   <- EF[, k]
    order_idx <- order(abs(weights), decreasing = TRUE)
    data.frame(
      FeatureID = order_idx[1:n_top],
      Weight    = round(weights[order_idx[1:n_top]], 4)
    )
  })
}

# Per-factor summary: Beta, log-rank p, sparsity %, PVE %
get_factor_summary_table <- function(EL, EF, EBeta, data) {
  K <- ncol(EL)
  p_vals <- sapply(1:K, function(k) {
    grp     <- ifelse(EL[, k] > median(EL[, k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ grp)
    1 - pchisq(sd_test$chisq, 1)
  })
  total_var <- sum(data$Y^2)
  pve <- sapply(1:K, function(k) {
    sum((outer(EL[, k], EF[, k]))^2) / total_var * 100
  })
  nonzero_pct <- colMeans(EF != 0) * 100
  data.frame(
    Factor      = 1:K,
    Beta        = round(EBeta, 3),
    LogRank_P   = round(p_vals, 4),
    NonZero_Pct = round(nonzero_pct, 2),
    PVE_Pct     = round(pve, 2)
  )
}

# ==============================================================================
# Data Generation (n=250, p=1000, K=5)
# Same DGP as run_modular_simulation.R
# ==============================================================================

set.seed(42)

n <- 250; p <- 1000; K <- 5

L_true <- matrix(rnorm(n * K), n, K)
F_true <- matrix(0, p, K)
for (k in 1:K) {
  active <- sample(1:p, round(p * 0.05))        # 5% sparse factor structure
  F_true[active, k] <- rnorm(length(active), 0, 5)
}
Y <- L_true %*% t(F_true) + matrix(rnorm(n * p), n, p)

B_true   <- c(1.5, -1.2, 0.8, -0.5, 0.0)
eta_true <- as.vector(L_true %*% B_true)
raw_times  <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
cens_times <- rexp(n, rate = 1 / 50)
time   <- pmin(raw_times, cens_times)
status <- as.integer(raw_times <= cens_times)

data <- list(Y = Y, time = time, status = status,
             L_true = L_true, F_true = F_true, B_true = B_true)

# ==============================================================================
# Fit Model
# ==============================================================================

cat("=== Factor-Wise Modular Supervised MF — Simulation ===\n")
cat(sprintf("  n=%d  p=%d  K=%d  seed=42  (factor-wise CAVI, V3 Algorithm 1)\n", n, p, K))
cat(sprintf("  Censoring rate: %.1f%%\n\n", 100 * mean(status == 0)))

res <- fit_supervised_mf_modular(Y, time, status, K = K,
                                  max_iter = 100, tol = 1e-3, verbose = TRUE)

# Unpack results (E-prefix naming from fit_modular.R)
EL     <- res$EL
EL2    <- res$EL2
EF     <- res$EF
EF2    <- res$EF2
EBeta  <- res$EBeta
EBeta2 <- res$EBeta2
Tau    <- res$Tau
history <- res$history

beta_sd <- sqrt(pmax(EBeta2 - EBeta^2, 0))

# ==============================================================================
# Summary Printout
# ==============================================================================

cat("\n=== FACTOR SUMMARY TABLE ===\n")
summary_tab <- get_factor_summary_table(EL, EF, EBeta, data)
print(summary_tab)

cat("\n=== ESTIMATED vs TRUE BETA ===\n")
print(data.frame(
  Factor    = 1:K,
  Beta_true = B_true,
  Beta_est  = round(EBeta, 3),
  Abs_Error = round(abs(EBeta - B_true), 3),
  Sign_Match = sign(EBeta) == sign(B_true) | B_true == 0
))

cat("\n=== MODEL PERFORMANCE (C-INDEX) ===\n")
perf <- get_cindex_comparison(EL, data)
cat(sprintf("  Top-5 PCA  C-index: %.3f\n", perf$c_original))
cat(sprintf("  Supervised C-index: %.3f\n", perf$c_latent))

cat("\n=== PROPORTIONAL HAZARDS TEST ===\n")
print(cox.zph(coxph(Surv(time, status) ~ EL)))

# ==============================================================================
# Save CSV Tables  -->  results/tables/factor_modular_sim/
# ==============================================================================

# Factor summary
write.csv(summary_tab,
          "results/tables/factor_modular_sim/factor_summary_table.csv",
          row.names = FALSE)

# Beta comparison
beta_df <- data.frame(
  Factor       = 1:K,
  Beta_true    = B_true,
  Beta_est     = round(EBeta, 4),
  Beta2_est    = round(EBeta2, 6),
  Posterior_SD = round(beta_sd, 4),
  Abs_Error    = round(abs(EBeta - B_true), 4),
  Sign_Match   = sign(EBeta) == sign(B_true) | B_true == 0
)
write.csv(beta_df,
          "results/tables/factor_modular_sim/beta_comparison_table.csv",
          row.names = FALSE)

# C-index comparison
cindex_df <- data.frame(
  Method  = c("Top-5 PCA", "Supervised Latent L"),
  C_Index = c(perf$c_original, perf$c_latent)
)
write.csv(cindex_df,
          "results/tables/factor_modular_sim/cindex_comparison.csv",
          row.names = FALSE)

# Convergence history
history_df <- data.frame(
  Iteration  = seq_along(history$rmse),
  RMSE       = history$rmse,
  ELBO_Proxy = history$elbo_proxy
)
write.csv(history_df,
          "results/tables/factor_modular_sim/convergence_history.csv",
          row.names = FALSE)

# Top features per GEP
top_feats <- get_top_features(EF, 10)
for (k in seq_along(top_feats)) {
  write.csv(top_feats[[k]],
            sprintf("results/tables/factor_modular_sim/top_features_GEP%d.csv", k),
            row.names = FALSE)
}

# Loading correlations (true vs estimated)
cors <- cor(L_true, EL)
write.csv(round(cors, 4),
          "results/tables/factor_modular_sim/loading_correlation_matrix.csv")

# PH test
ph_test <- cox.zph(coxph(Surv(time, status) ~ EL))
ph_df <- data.frame(
  Factor  = rownames(ph_test$table),
  Chisq   = round(ph_test$table[, 1], 4),
  DF      = ph_test$table[, 2],
  P_Value = round(ph_test$table[, 3], 4)
)
write.csv(ph_df,
          "results/tables/factor_modular_sim/ph_test_results.csv",
          row.names = FALSE)

cat("\nCSV tables saved to results/tables/factor_modular_sim/\n")

# ==============================================================================
# Generate Figures  -->  results/figures/factor_modular_sim/
# ==============================================================================

# --- Figure 1: RMSE Convergence Trace ---
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig1_rmse_trace.pdf",
                        width = 8, height = 5)
  else              png("results/figures/factor_modular_sim/fig1_rmse_trace.png",
                        width = 800, height = 500, res = 120)
  par(mar = c(5, 5, 4, 2))
  plot(history$rmse, type = "l", lwd = 2.5, col = "#1f77b4",
       main = "Figure 1: Reconstruction RMSE Across CAVI Iterations\n(Factor-Wise Modular CAVI)",
       xlab = "Iteration", ylab = "RMSE", bty = "n",
       cex.lab = 1.2, cex.main = 1.3)
  abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
  legend("topright", legend = c("RMSE", "True Noise SD = 1.0"),
         col = c("#1f77b4", "#d62728"), lty = c(1, 2), lwd = c(2.5, 1.5), bty = "n")
  grid(col = "lightgray", lty = "dotted")
  dev.off()
}

# --- Figure 2: ELBO Proxy Trace ---
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig2_elbo_proxy.pdf",
                        width = 8, height = 5)
  else              png("results/figures/factor_modular_sim/fig2_elbo_proxy.png",
                        width = 800, height = 500, res = 120)
  par(mar = c(5, 5, 4, 2))
  plot(history$elbo_proxy, type = "l", lwd = 2.5, col = "#2ca02c",
       main = "Figure 2: Genomics ELBO Proxy Across Iterations\n(Factor-Wise Modular CAVI)",
       xlab = "Iteration",
       ylab = expression(E[q]*"[log P(Y | L, F, "*tau*")]"),
       bty = "n", cex.lab = 1.2, cex.main = 1.3)
  grid(col = "lightgray", lty = "dotted")
  dev.off()
}

# --- Figure 3: Beta Comparison (True vs Estimated) ---
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig3_beta_comparison.pdf",
                        width = 8, height = 5.5)
  else              png("results/figures/factor_modular_sim/fig3_beta_comparison.png",
                        width = 800, height = 550, res = 120)
  par(mar = c(5, 5, 4, 2))
  x_pos <- 1:K
  plot(x_pos, EBeta, pch = 16, cex = 1.8, col = "#1f77b4",
       ylim = range(c(B_true, EBeta + 1.96 * beta_sd, EBeta - 1.96 * beta_sd)) * 1.2,
       xlab = "Factor", ylab = expression(beta[k]),
       main = "Figure 3: Estimated vs. True Survival Coefficients\n(Factor-Wise Modular CAVI)",
       bty = "n", cex.lab = 1.2, cex.main = 1.3, xaxt = "n")
  axis(1, at = 1:K)
  arrows(x_pos, EBeta - 1.96 * beta_sd, x_pos, EBeta + 1.96 * beta_sd,
         angle = 90, code = 3, length = 0.08, col = "#1f77b4", lwd = 1.5)
  points(x_pos, B_true, pch = 4, cex = 2, col = "#d62728", lwd = 2.5)
  abline(h = 0, col = "gray50", lty = 3)
  legend("bottomleft", legend = c("Estimated (95% CI)", "True"),
         col = c("#1f77b4", "#d62728"), pch = c(16, 4),
         pt.cex = c(1.8, 2), pt.lwd = c(1, 2.5), bty = "n")
  grid(col = "lightgray", lty = "dotted")
  dev.off()
}

# --- Figure 4: GEP Heatmap ---
n_features    <- 50
top_var_genes <- order(rowSums(abs(EF)), decreasing = TRUE)[1:n_features]
F_sub <- EF[top_var_genes, ]
palette <- colorRampPalette(c("blue", "white", "red"))(100)
max_val <- max(abs(F_sub))

for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig4_gep_heatmap.pdf",
                        width = 10, height = 6)
  else              png("results/figures/factor_modular_sim/fig4_gep_heatmap.png",
                        width = 1000, height = 600, res = 120)
  layout(matrix(1:2, ncol = 2), widths = c(5, 1))
  par(mar = c(6, 4, 4, 1))
  image(1:nrow(F_sub), 1:ncol(F_sub), F_sub,
        main = "Figure 4: GEP Feature Weights (Top 50 Features)\n(Factor-Wise Modular CAVI)",
        xlab = "Feature ID", ylab = "Latent Factor",
        col = palette, axes = FALSE, zlim = c(-max_val, max_val))
  axis(1, at = 1:nrow(F_sub), labels = top_var_genes, las = 2, cex.axis = 0.55)
  axis(2, at = 1:ncol(F_sub), labels = paste0("F", 1:ncol(F_sub)), las = 1)
  box()
  par(mar = c(6, 1, 4, 3))
  legend_image <- as.matrix(seq(-max_val, max_val, length.out = 100))
  image(1, seq(-max_val, max_val, length.out = 100), t(legend_image),
        col = palette, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.8)
  mtext("Weight", side = 4, line = 2, cex = 0.8)
  layout(1)
  dev.off()
}

# --- Figure 5: Kaplan-Meier Survival Curves per Factor ---
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig5_kaplan_meier.pdf",
                        width = 12, height = 8)
  else              png("results/figures/factor_modular_sim/fig5_kaplan_meier.png",
                        width = 1200, height = 800, res = 120)
  par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
  for (k in 1:K) {
    groups  <- ifelse(EL[, k] > median(EL[, k]), "High Score", "Low Score")
    km      <- survfit(Surv(time, status) ~ groups)
    p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
    plot(km, col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n",
         main = paste0("Factor ", k, " (beta = ", round(EBeta[k], 2), ")\n", p_label),
         xlab = "Time", ylab = "Survival Probability")
    legend("bottomleft", legend = c("High", "Low"),
           col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n", cex = 0.8)
  }
  dev.off()
}

# --- Figure 6: Signal Recovery (Best-matched factor) ---
cors_mat    <- cor(L_true, EL)
best_match  <- apply(abs(cors_mat), 2, which.max)
target_est  <- which.max(abs(EBeta))
target_true <- best_match[target_est]
sign_corr   <- sign(cors_mat[target_true, target_est])
r_val       <- cors_mat[target_true, target_est]

for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig6_signal_recovery.pdf",
                        width = 7, height = 7)
  else              png("results/figures/factor_modular_sim/fig6_signal_recovery.png",
                        width = 700, height = 700, res = 120)
  par(mar = c(5, 5, 4, 2))
  plot(L_true[, target_true], EL[, target_est] * sign_corr,
       main = sprintf("Figure 6: Signal Recovery (Factor-Wise Modular CAVI)\n(Est Factor %d vs True Factor %d, r = %.3f)",
                      target_est, target_true, abs(r_val)),
       xlab = "Ground Truth Loading",
       ylab = "Estimated Loading (sign-corrected)",
       col = rgb(0, 0, 0, 0.35), pch = 16, cex = 1.2, bty = "n",
       cex.lab = 1.2, cex.main = 1.2)
  abline(0, 1, col = "#d62728", lwd = 2, lty = 2)
  abline(lm(I(EL[, target_est] * sign_corr) ~ L_true[, target_true]),
         col = "#1f77b4", lwd = 2)
  legend("topleft", legend = c("Identity line", "Best fit"),
         col = c("#d62728", "#1f77b4"), lty = c(2, 1), lwd = 2, bty = "n")
  grid(col = "lightgray", lty = "dotted")
  dev.off()
}

# --- Figure 7: Loading Correlation Heatmap (True vs Estimated) ---
cor_mat  <- cor(L_true, EL)
palette2 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig7_loading_correlations.pdf",
                        width = 7, height = 6)
  else              png("results/figures/factor_modular_sim/fig7_loading_correlations.png",
                        width = 700, height = 600, res = 120)
  par(mar = c(5, 5, 4, 5))
  image(1:K, 1:K, abs(cor_mat),
        col = palette2, zlim = c(0, 1),
        main = "Figure 7: |Correlation| Between True and Estimated Loadings\n(Factor-Wise Modular CAVI)",
        xlab = "True Factor", ylab = "Estimated Factor",
        axes = FALSE, cex.lab = 1.2, cex.main = 1.2)
  axis(1, at = 1:K); axis(2, at = 1:K)
  for (i in 1:K) for (j in 1:K) {
    text(i, j, sprintf("%.2f", cor_mat[i, j]), cex = 1.0,
         col = if (abs(cor_mat[i, j]) > 0.5) "white" else "black")
  }
  box()
  dev.off()
}

# --- Figure 8: Tau (Noise Precision) Distribution ---
for (ext in c("pdf", "png")) {
  if (ext == "pdf") pdf("results/figures/factor_modular_sim/fig8_tau_distribution.pdf",
                        width = 8, height = 5)
  else              png("results/figures/factor_modular_sim/fig8_tau_distribution.png",
                        width = 800, height = 500, res = 120)
  par(mar = c(5, 5, 4, 2))
  hist(Tau, breaks = 50, col = "#1f77b4AA", border = "white",
       main = "Figure 8: Estimated Feature-Specific Noise Precision\n(Factor-Wise Modular CAVI)",
       xlab = expression(hat(tau)[j]), ylab = "Count",
       cex.lab = 1.2, cex.main = 1.3)
  abline(v = 1.0, col = "#d62728", lwd = 2, lty = 2)
  abline(v = median(Tau), col = "#2ca02c", lwd = 2, lty = 3)
  legend("topright",
         legend = c("True (tau = 1.0)",
                    sprintf("Median est. (%.3f)", median(Tau))),
         col = c("#d62728", "#2ca02c"), lty = c(2, 3), lwd = 2, bty = "n")
  dev.off()
}

cat(sprintf("\nAll figures saved to results/figures/factor_modular_sim/ (PDF + PNG)\n"))
cat(sprintf("Total: %d CSVs + %d figure pairs\n",
            length(list.files("results/tables/factor_modular_sim", pattern = "\\.csv$")),
            length(list.files("results/figures/factor_modular_sim", pattern = "\\.pdf$"))))
