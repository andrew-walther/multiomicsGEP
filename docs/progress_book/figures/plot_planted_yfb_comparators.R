# ============================================================
# Script: plot_planted_yfb_comparators.R
# Purpose: Plot planted loading matches and the matched simulation comparators.
# Created: 2026-11-03
# Dependencies: Base R
# ============================================================
# Run from the repository root. Inputs are saved synthetic replicate tables.
root <- "results/benchmark_sim/outputs"
figures <- "docs/progress_book/figures/2026-11-06"
runs <- read.csv(file.path(root, "supervision_diagnostics/c3c4_n400_b3.3/results.csv"))
arms <- c("alpha_F_0", "c4_point_laplace")
noise <- c(0.25, 1)
colors <- c("#704E9B", "#12766F")
set.seed(42)

png(file.path(figures, "planted-loading-recovery.png"), 1500, 770, res = 150)
par(mfrow = c(1, 2), mar = c(5, 4.8, 3, 1), oma = c(2, 0, 0, 0), las = 1)
for (level in noise) {
  plot(NA, xlim = c(0.5, 2.5), ylim = c(0, 1), xaxt = "n", xlab = "",
       ylab = "Best absolute loading correlation", main = paste("Noise", level))
  abline(h = 0.7, lty = 2, col = "gray50")
  labels <- character(2)
  for (i in seq_along(arms)) {
    d <- subset(runs, arm == arms[i] & noise == level)
    expected <- if (i == 1) 1 else if (level == 0.25) 17 else 19
    stopifnot(nrow(d) == 20L, !anyDuplicated(d$seed), all(is.finite(d$cor_planted)),
              sum(d$recovered) == expected, all(d$cor_planted[d$recovered] >= 0.7))
    points(jitter(rep(i, nrow(d)), amount = 0.16), d$cor_planted,
           pch = ifelse(d$recovered, 19, 1), col = colors[i], cex = 1.1)
    labels[i] <- sprintf("%s\n%d/%d recovered", c("Feedback off", "Joint YFB")[i],
                         sum(d$recovered), nrow(d))
  }
  axis(1, 1:2, labels, tick = FALSE)
}
mtext("Filled points meet the loading and survival-direction criteria. Dashed line = 0.70.",
      side = 1, outer = TRUE, line = 0.5, cex = 0.86)
dev.off()

path <- file.path(root, "supervision_comparators/round2_retry1")
panel <- read.csv(file.path(path, "results.csv"))
summary <- read.csv(file.path(path, "summary.csv"))
methods <- c("supervised", "ebmf_alone", "ebmf_cox", "ebmf_lasso_cox", "pca_cox")
labels <- c("Joint YFB", "EBMF alone", "EBMF + Cox", "EBMF + lasso Cox", "PCA + Cox")
stopifnot(nrow(panel) == 200L, !anyDuplicated(panel[c("seed", "noise", "method")]))

png(file.path(figures, "planted-matched-comparators.png"), 1800, 1260, res = 150)
par(mfrow = c(2, 2), mar = c(4.2, 10.6, 3, 2.5), oma = c(1, 0, 0, 0), las = 1)
for (level in noise) {
  stopifnot(sum(panel$noise == level) == 100L,
            all(table(panel$seed[panel$noise == level]) == length(methods)))
  for (metric in c("cor_planted", "C")) {
    correlation <- metric == "cor_planted"
    title <- if (correlation) "Loading match" else "Held-out concordance"
    xlab <- if (correlation) "Absolute loading correlation of selected program" else "C-index"
    plot(NA, xlim = c(0, 1.19), ylim = c(0.45, 5.55), xaxt = "n", yaxt = "n",
         xlab = xlab, ylab = "", main = paste0(title, " · noise ", level))
    axis(1, seq(0, 1, 0.2))
    axis(2, 5:1, labels, tick = FALSE, cex.axis = 0.85)
    abline(v = if (correlation) 0.7 else 0.5, lty = 2, col = "gray60")
    for (i in seq_along(methods)) {
      d <- subset(panel, method == methods[i] & noise == level)
      s <- subset(summary, method == methods[i] & noise == level)
      values <- d[[metric]][is.finite(d[[metric]])]
      expected <- if (correlation) s$n_cor_defined else s$n_C_defined
      stopifnot(nrow(d) == 20L, length(values) == expected, sum(d$recovered) == s$recovered)
      y <- 6 - i
      color <- if (i == 1) colors[2] else colors[1]
      if (length(values)) {
        points(values, jitter(rep(y, length(values)), amount = 0.16),
               pch = 19, cex = 0.65, col = adjustcolor(color, alpha.f = 0.55))
        segments(mean(values), y - 0.23, mean(values), y + 0.23, col = color, lwd = 2)
      }
      count <- if (!correlation && i == 2) "No risk" else sprintf("%d/20", length(values))
      text(1.09, y, count, cex = 0.72)
    }
  }
}
mtext("Vertical marks show means. Counts show defined values out of 20 datasets.",
      side = 1, outer = TRUE, line = 0, cex = 0.85)
dev.off()
