# ==============================================================================
# Run V2 Simulation and Export All Results
# Saves: text output, CSV tables, and PDF/PNG figures to results/
# ==============================================================================

# Set working directory to repo root
setwd("/Users/ajwalther/GithubProjects/multiomicsGEP")

# Source the V2 script in simulated mode — this runs the full pipeline
# We capture the console output
sink("results/simulation_console_output.txt")
cat("================================================================\n")
cat("Supervised Bayesian MF V2 — Full Simulation Run\n")
cat(sprintf("Timestamp: %s\n", Sys.time()))
cat("================================================================\n\n")

source("code/Supervised_Bayesian_MF_V2.R")

sink()  # stop capturing

# At this point, 'data' and 'res' are in the global environment from the source()

# ==============================================================================
# Save Tables as CSV
# ==============================================================================

# Factor summary table
factor_summary <- get_factor_summary_table(res, data)
write.csv(factor_summary, "results/factor_summary_table.csv", row.names = FALSE)

# Beta comparison table
beta_comparison <- data.frame(
  Factor    = 1:5,
  Beta_true = data$B_true,
  Beta_est  = round(res$Beta, 4),
  Beta2_est = round(res$Beta2, 6),
  Posterior_SD = round(sqrt(pmax(res$Beta2 - res$Beta^2, 0)), 4),
  Abs_Error = round(abs(res$Beta - data$B_true), 4),
  Sign_Match = sign(res$Beta) == sign(data$B_true) | data$B_true == 0
)
write.csv(beta_comparison, "results/beta_comparison_table.csv", row.names = FALSE)

# C-index comparison
perf <- get_cindex_comparison(res$L, data)
cindex_df <- data.frame(
  Method = c("Top-5 PCA", "Supervised Latent L"),
  C_Index = c(perf$c_original, perf$c_latent)
)
write.csv(cindex_df, "results/cindex_comparison.csv", row.names = FALSE)

# Convergence history
history_df <- data.frame(
  Iteration  = 1:length(res$history$rmse),
  RMSE       = res$history$rmse,
  ELBO_Proxy = res$history$elbo_proxy
)
write.csv(history_df, "results/convergence_history.csv", row.names = FALSE)

# Top features per GEP
top_feats <- get_top_features(res$F, 10)
for (k in 1:length(top_feats)) {
  write.csv(top_feats[[k]],
            sprintf("results/top_features_GEP%d.csv", k),
            row.names = FALSE)
}

# Correlation between true and estimated loadings
cors <- cor(data$L_true, res$L)
write.csv(round(cors, 4), "results/loading_correlation_matrix.csv")

# PH test
ph_test <- cox.zph(coxph(Surv(data$time, data$status) ~ res$L))
ph_df <- data.frame(
  Factor = rownames(ph_test$table),
  Chisq  = round(ph_test$table[, 1], 4),
  DF     = ph_test$table[, 2],
  P_Value = round(ph_test$table[, 3], 4)
)
write.csv(ph_df, "results/ph_test_results.csv", row.names = FALSE)

cat("CSV tables saved to results/\n")

# ==============================================================================
# Generate Figures
# ==============================================================================

K <- ncol(res$L)

# --- Figure 1: RMSE Convergence Trace ---
pdf("results/figures/fig1_rmse_trace.pdf", width = 8, height = 5)
par(mar = c(5, 5, 4, 2))
rmse_vals <- res$history$rmse
plot(rmse_vals, type = "l", lwd = 2.5, col = "#1f77b4",
     main = "Figure 1: Reconstruction RMSE Across CAVI Iterations",
     xlab = "Iteration", ylab = "RMSE", bty = "n", cex.lab = 1.2, cex.main = 1.3)
abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
legend("topright", legend = c("RMSE", "True Noise SD = 1.0"),
       col = c("#1f77b4", "#d62728"), lty = c(1, 2), lwd = c(2.5, 1.5), bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

png("results/figures/fig1_rmse_trace.png", width = 800, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
plot(rmse_vals, type = "l", lwd = 2.5, col = "#1f77b4",
     main = "Figure 1: Reconstruction RMSE Across CAVI Iterations",
     xlab = "Iteration", ylab = "RMSE", bty = "n", cex.lab = 1.2, cex.main = 1.3)
abline(h = 1.0, col = "#d62728", lty = 2, lwd = 1.5)
legend("topright", legend = c("RMSE", "True Noise SD = 1.0"),
       col = c("#1f77b4", "#d62728"), lty = c(1, 2), lwd = c(2.5, 1.5), bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

# --- Figure 2: ELBO Proxy Trace ---
pdf("results/figures/fig2_elbo_proxy.pdf", width = 8, height = 5)
par(mar = c(5, 5, 4, 2))
elbo_vals <- res$history$elbo_proxy
plot(elbo_vals, type = "l", lwd = 2.5, col = "#2ca02c",
     main = "Figure 2: Genomics ELBO Proxy Across Iterations",
     xlab = "Iteration", ylab = expression(E[q]*"[log P(Y | L, F, "*tau*")]"),
     bty = "n", cex.lab = 1.2, cex.main = 1.3)
grid(col = "lightgray", lty = "dotted")
dev.off()

png("results/figures/fig2_elbo_proxy.png", width = 800, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
plot(elbo_vals, type = "l", lwd = 2.5, col = "#2ca02c",
     main = "Figure 2: Genomics ELBO Proxy Across Iterations",
     xlab = "Iteration", ylab = expression(E[q]*"[log P(Y | L, F, "*tau*")]"),
     bty = "n", cex.lab = 1.2, cex.main = 1.3)
grid(col = "lightgray", lty = "dotted")
dev.off()

# --- Figure 3: Beta Comparison (True vs Estimated) ---
pdf("results/figures/fig3_beta_comparison.pdf", width = 8, height = 5.5)
par(mar = c(5, 5, 4, 2))
beta_true <- data$B_true
beta_est  <- res$Beta
beta_sd   <- sqrt(pmax(res$Beta2 - res$Beta^2, 0))
x_pos <- 1:K

plot(x_pos, beta_est, pch = 16, cex = 1.8, col = "#1f77b4",
     ylim = range(c(beta_true, beta_est + beta_sd, beta_est - beta_sd)) * 1.2,
     xlab = "Factor", ylab = expression(beta[k]),
     main = "Figure 3: Estimated vs. True Survival Coefficients",
     bty = "n", cex.lab = 1.2, cex.main = 1.3, xaxt = "n")
axis(1, at = 1:K)
# Error bars (posterior SD)
arrows(x_pos, beta_est - 1.96 * beta_sd, x_pos, beta_est + 1.96 * beta_sd,
       angle = 90, code = 3, length = 0.08, col = "#1f77b4", lwd = 1.5)
# True values
points(x_pos, beta_true, pch = 4, cex = 2, col = "#d62728", lwd = 2.5)
abline(h = 0, col = "gray50", lty = 3)
legend("bottomleft",
       legend = c("Estimated (95% CI)", "True"),
       col = c("#1f77b4", "#d62728"), pch = c(16, 4), pt.cex = c(1.8, 2),
       pt.lwd = c(1, 2.5), bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

png("results/figures/fig3_beta_comparison.png", width = 800, height = 550, res = 120)
par(mar = c(5, 5, 4, 2))
plot(x_pos, beta_est, pch = 16, cex = 1.8, col = "#1f77b4",
     ylim = range(c(beta_true, beta_est + beta_sd, beta_est - beta_sd)) * 1.2,
     xlab = "Factor", ylab = expression(beta[k]),
     main = "Figure 3: Estimated vs. True Survival Coefficients",
     bty = "n", cex.lab = 1.2, cex.main = 1.3, xaxt = "n")
axis(1, at = 1:K)
arrows(x_pos, beta_est - 1.96 * beta_sd, x_pos, beta_est + 1.96 * beta_sd,
       angle = 90, code = 3, length = 0.08, col = "#1f77b4", lwd = 1.5)
points(x_pos, beta_true, pch = 4, cex = 2, col = "#d62728", lwd = 2.5)
abline(h = 0, col = "gray50", lty = 3)
legend("bottomleft",
       legend = c("Estimated (95% CI)", "True"),
       col = c("#1f77b4", "#d62728"), pch = c(16, 4), pt.cex = c(1.8, 2),
       pt.lwd = c(1, 2.5), bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

# --- Figure 4: GEP Heatmap ---
pdf("results/figures/fig4_gep_heatmap.pdf", width = 10, height = 6)
n_features <- 50
top_var_genes <- order(rowSums(abs(res$F)), decreasing = TRUE)[1:n_features]
F_sub   <- res$F[top_var_genes, ]
palette <- colorRampPalette(c("blue", "white", "red"))(100)
max_val <- max(abs(F_sub))
layout(matrix(1:2, ncol = 2), widths = c(5, 1))
par(mar = c(6, 4, 4, 1))
image(1:nrow(F_sub), 1:ncol(F_sub), F_sub,
      main = "Figure 4: GEP Feature Weights (Top 50 Features)",
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

png("results/figures/fig4_gep_heatmap.png", width = 1000, height = 600, res = 120)
layout(matrix(1:2, ncol = 2), widths = c(5, 1))
par(mar = c(6, 4, 4, 1))
image(1:nrow(F_sub), 1:ncol(F_sub), F_sub,
      main = "Figure 4: GEP Feature Weights (Top 50 Features)",
      xlab = "Feature ID", ylab = "Latent Factor",
      col = palette, axes = FALSE, zlim = c(-max_val, max_val))
axis(1, at = 1:nrow(F_sub), labels = top_var_genes, las = 2, cex.axis = 0.55)
axis(2, at = 1:ncol(F_sub), labels = paste0("F", 1:ncol(F_sub)), las = 1)
box()
par(mar = c(6, 1, 4, 3))
image(1, seq(-max_val, max_val, length.out = 100), t(legend_image),
      col = palette, axes = FALSE, xlab = "", ylab = "")
axis(4, las = 1, cex.axis = 0.8)
mtext("Weight", side = 4, line = 2, cex = 0.8)
layout(1)
dev.off()

# --- Figure 5: Kaplan-Meier Survival Curves per Factor ---
pdf("results/figures/fig5_kaplan_meier.pdf", width = 12, height = 8)
summary_tab <- get_factor_summary_table(res, data)
par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
for (k in 1:K) {
  groups  <- ifelse(res$L[, k] > median(res$L[, k]), "High Score", "Low Score")
  km      <- survfit(Surv(data$time, data$status) ~ groups)
  p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
  plot(km, col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n",
       main = paste0("Factor ", k, " (beta = ", round(res$Beta[k], 2), ")\n", p_label),
       xlab = "Time", ylab = "Survival Probability")
  legend("bottomleft", legend = c("High", "Low"),
         col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()

png("results/figures/fig5_kaplan_meier.png", width = 1200, height = 800, res = 120)
par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
for (k in 1:K) {
  groups  <- ifelse(res$L[, k] > median(res$L[, k]), "High Score", "Low Score")
  km      <- survfit(Surv(data$time, data$status) ~ groups)
  p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
  plot(km, col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n",
       main = paste0("Factor ", k, " (beta = ", round(res$Beta[k], 2), ")\n", p_label),
       xlab = "Time", ylab = "Survival Probability")
  legend("bottomleft", legend = c("High", "Low"),
         col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()

# --- Figure 6: Signal Recovery (Best-matched factor) ---
pdf("results/figures/fig6_signal_recovery.pdf", width = 7, height = 7)
par(mar = c(5, 5, 4, 2))
cors_mat    <- cor(data$L_true, res$L)
best_match  <- apply(abs(cors_mat), 2, which.max)
target_est  <- which.max(abs(res$Beta))
target_true <- best_match[target_est]
sign_corr   <- sign(cors_mat[target_true, target_est])
r_val       <- round(cors_mat[target_true, target_est], 3)
plot(data$L_true[, target_true], res$L[, target_est] * sign_corr,
     main = sprintf("Figure 6: Signal Recovery\n(Est Factor %d vs True Factor %d, r = %.3f)",
                    target_est, target_true, abs(r_val)),
     xlab = "Ground Truth Loading", ylab = "Estimated Loading (sign-corrected)",
     col = rgb(0, 0, 0, 0.35), pch = 16, cex = 1.2, bty = "n",
     cex.lab = 1.2, cex.main = 1.2)
abline(0, 1, col = "#d62728", lwd = 2, lty = 2)
abline(lm(I(res$L[, target_est] * sign_corr) ~ data$L_true[, target_true]),
       col = "#1f77b4", lwd = 2)
legend("topleft", legend = c("Identity line", "Best fit"),
       col = c("#d62728", "#1f77b4"), lty = c(2, 1), lwd = 2, bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

png("results/figures/fig6_signal_recovery.png", width = 700, height = 700, res = 120)
par(mar = c(5, 5, 4, 2))
plot(data$L_true[, target_true], res$L[, target_est] * sign_corr,
     main = sprintf("Figure 6: Signal Recovery\n(Est Factor %d vs True Factor %d, r = %.3f)",
                    target_est, target_true, abs(r_val)),
     xlab = "Ground Truth Loading", ylab = "Estimated Loading (sign-corrected)",
     col = rgb(0, 0, 0, 0.35), pch = 16, cex = 1.2, bty = "n",
     cex.lab = 1.2, cex.main = 1.2)
abline(0, 1, col = "#d62728", lwd = 2, lty = 2)
abline(lm(I(res$L[, target_est] * sign_corr) ~ data$L_true[, target_true]),
       col = "#1f77b4", lwd = 2)
legend("topleft", legend = c("Identity line", "Best fit"),
       col = c("#d62728", "#1f77b4"), lty = c(2, 1), lwd = 2, bty = "n")
grid(col = "lightgray", lty = "dotted")
dev.off()

# --- Figure 7: Loading Correlation Heatmap (True vs Estimated) ---
pdf("results/figures/fig7_loading_correlations.pdf", width = 7, height = 6)
par(mar = c(5, 5, 4, 5))
cor_mat <- cor(data$L_true, res$L)
palette2 <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
image(1:K, 1:K, abs(cor_mat),
      col = palette2, zlim = c(0, 1),
      main = "Figure 7: |Correlation| Between True and Estimated Loadings",
      xlab = "True Factor", ylab = "Estimated Factor",
      axes = FALSE, cex.lab = 1.2, cex.main = 1.2)
axis(1, at = 1:K); axis(2, at = 1:K)
# Add text labels
for (i in 1:K) for (j in 1:K) {
  text(i, j, sprintf("%.2f", cor_mat[i, j]), cex = 1.0,
       col = if (abs(cor_mat[i, j]) > 0.5) "white" else "black")
}
box()
dev.off()

png("results/figures/fig7_loading_correlations.png", width = 700, height = 600, res = 120)
par(mar = c(5, 5, 4, 5))
image(1:K, 1:K, abs(cor_mat),
      col = palette2, zlim = c(0, 1),
      main = "Figure 7: |Correlation| Between True and Estimated Loadings",
      xlab = "True Factor", ylab = "Estimated Factor",
      axes = FALSE, cex.lab = 1.2, cex.main = 1.2)
axis(1, at = 1:K); axis(2, at = 1:K)
for (i in 1:K) for (j in 1:K) {
  text(i, j, sprintf("%.2f", cor_mat[i, j]), cex = 1.0,
       col = if (abs(cor_mat[i, j]) > 0.5) "white" else "black")
}
box()
dev.off()

# --- Figure 8: Noise Precision Recovery ---
pdf("results/figures/fig8_tau_distribution.pdf", width = 8, height = 5)
par(mar = c(5, 5, 4, 2))
hist(res$Tau, breaks = 50, col = "#1f77b4AA", border = "white",
     main = "Figure 8: Estimated Feature-Specific Noise Precision",
     xlab = expression(hat(tau)[j]), ylab = "Count",
     cex.lab = 1.2, cex.main = 1.3)
abline(v = 1.0, col = "#d62728", lwd = 2, lty = 2)
abline(v = median(res$Tau), col = "#2ca02c", lwd = 2, lty = 3)
legend("topright",
       legend = c(sprintf("True (tau = 1.0)"),
                  sprintf("Median est. (%.3f)", median(res$Tau))),
       col = c("#d62728", "#2ca02c"), lty = c(2, 3), lwd = 2, bty = "n")
dev.off()

png("results/figures/fig8_tau_distribution.png", width = 800, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
hist(res$Tau, breaks = 50, col = "#1f77b4AA", border = "white",
     main = "Figure 8: Estimated Feature-Specific Noise Precision",
     xlab = expression(hat(tau)[j]), ylab = "Count",
     cex.lab = 1.2, cex.main = 1.3)
abline(v = 1.0, col = "#d62728", lwd = 2, lty = 2)
abline(v = median(res$Tau), col = "#2ca02c", lwd = 2, lty = 3)
legend("topright",
       legend = c(sprintf("True (tau = 1.0)"),
                  sprintf("Median est. (%.3f)", median(res$Tau))),
       col = c("#d62728", "#2ca02c"), lty = c(2, 3), lwd = 2, bty = "n")
dev.off()

cat("\nAll figures saved to results/figures/ (PDF + PNG)\n")
cat(sprintf("Total files: %d CSVs + %d figure pairs\n",
            length(list.files("results", pattern = "\\.csv$")),
            length(list.files("results/figures", pattern = "\\.pdf$"))))
