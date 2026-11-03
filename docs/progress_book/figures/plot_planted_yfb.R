# ============================================================
# Script: plot_planted_yfb.R
# Purpose: Plot recovery and validation C from saved simulation replicates.
# Created: 2026-11-03
# Dependencies: Base R
# ============================================================
# Run from the repository root. No model fitting or patient data is required.
input <- "results/benchmark_sim/outputs/supervision_diagnostics/c3c4_n400_b3.3/results.csv"
output <- "docs/progress_book/figures/2026-11-06/planted-yfb-performance.png"
results <- read.csv(input)
arms <- c("alpha_F_0", "c4_point_laplace")
noise <- c(0.25, 1)
cells <- expand.grid(arm = arms, noise = noise, stringsAsFactors = FALSE)
recovery <- n <- numeric(nrow(cells))
validation <- vector("list", nrow(cells))

for (i in seq_len(nrow(cells))) {
  rows <- results[results$arm == cells$arm[i] & results$noise == cells$noise[i], ]
  stopifnot(nrow(rows) == 20L, !anyDuplicated(rows$seed), !anyNA(rows$recovered))
  n[i] <- nrow(rows)
  recovery[i] <- sum(rows$recovered)
  validation[[i]] <- rows$c_validation[is.finite(rows$c_validation)]
}

colors <- c("#704E9B", "#12766F")
positions <- c(1, 2, 3.5, 4.5)
png(output, width = 1620, height = 790, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.4, 1), las = 1, oma = c(2, 0, 0, 0))
plot(NA, xlim = c(0.4, 5.1), ylim = c(0, 23), xaxt = "n", yaxt = "n",
     xlab = "Simulation noise", ylab = "Replicates recovering the planted program",
     main = "Recovery of the prognostic program")
axis(1, c(1.5, 4), noise)
axis(2, seq(0, 20, 5))
abline(h = seq(0, 20, 5), col = "gray90")
rect(positions - 0.28, 0, positions + 0.28, recovery, col = colors, border = NA)
text(positions, recovery + 1, paste0(recovery, "/", n), cex = 0.9)

plot(NA, xlim = c(0.4, 5.1), ylim = c(0.43, 0.96), xaxt = "n",
     xlab = "Simulation noise", ylab = "Held-out concordance index",
     main = "Validation across replicates")
axis(1, c(1.5, 4), noise)
abline(h = 0.5, col = "gray60", lty = 2)
for (i in seq_along(positions)) {
  values <- sort(validation[[i]])
  color <- colors[match(cells$arm[i], arms)]
  offsets <- seq(-0.16, 0.16, length.out = length(values))
  points(positions[i] + offsets, values, pch = 16, cex = 0.65,
         col = adjustcolor(color, alpha.f = 0.58))
  segments(positions[i] - 0.23, mean(values), positions[i] + 0.23,
           mean(values), col = color, lwd = 3)
  text(positions[i], 0.445, paste0("n = ", length(values)), cex = 0.7)
}
mtext("Purple  Feedback off       Teal  Joint YFB       Horizontal marks  Group means",
      side = 1, outer = TRUE, line = 0.5, cex = 0.9)
dev.off()
