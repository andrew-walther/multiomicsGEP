# ============================================================
# Script: plot_weight_curves_cv.R
# Purpose: Redraw the cross-validated-weight simulation curves for the progress book with plain labels.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: base R, grDevices
# Run from the repository root: Rscript docs/progress_book/figures/plot_weight_curves_cv.R
# ============================================================
# Reads the saved aggregate weight_curves.csv written by
# results/benchmark_sim/summarize_yfb_supervision_strength.R and draws the same
# panels with reader-facing legend labels. It refits and rewrites nothing else.

OUT <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
FIG <- "docs/progress_book/figures/2026-11-06"
curve <- read.csv(file.path(OUT, "weight_curves.csv"))
stopifnot(all(c("A1", "A2") %in% curve$arm))

png(file.path(FIG, "s13-weight-curves.png"), width = 1600, height = 1000, res = 120)
op <- par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
mets <- c(recovered = "Recovered (of 10)", mean_false_active = "Mean false-active",
          slope_median = "Calibration slope (median)", converged = "Converged (of 10)",
          c_validation = "Validation C", cv_c = "Training CV C",
          share_match = "Survival share, matched", mean_survival_only = "Survival-only directions")
# Simulation cells: planted variance share (0.25%, 1%, 5%) or the weaker-survival cell.
cols <- c(weak = "#d95f02", v0.0025 = "#1b9e77", v0.01 = "#7570b3", v0.05 = "#e7298a")
labs <- c(weak = "Weaker survival", v0.0025 = "0.25% share", v0.01 = "1% share", v0.05 = "5% share")
for (mt in names(mets)) {
  plot(NA, xlim = c(1, 50), ylim = range(curve[[mt]], na.rm = TRUE), log = "x", xlab = "w_S",
       ylab = mets[[mt]], main = mets[[mt]])
  if (mt == "slope_median") abline(h = c(.8, 1, 1.25), lty = c(3, 2, 3), col = "grey50")
  for (a in c("A1", "A2")) for (cl in names(cols)) {
    z <- curve[curve$arm == a & curve$cell == cl, ]
    lines(z$w_S, z[[mt]], col = cols[[cl]], lty = if (a == "A1") 1 else 2, type = "b", pch = 19, cex = .6)
  }
}
legend("topright", c(labs[names(cols)], "Block-weighted fit (solid)", "Joint-posterior fit (dashed)"),
       col = c(cols, "black", "black"), lty = c(1, 1, 1, 1, 1, 2), bty = "n", cex = .8)
par(op); dev.off()
