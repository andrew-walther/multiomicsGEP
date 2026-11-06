# ============================================================
# Script: plot_supervision_strength.R
# Purpose: November 6 chapter figures for S-13 (weight curves) and S-14 (supervised minus
#          two-step validation C near the detection threshold).
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: base R, grDevices
# Run (from the repository root): Rscript docs/progress_book/figures/plot_supervision_strength.R
# ============================================================
# Reads saved aggregate CSVs only (synthetic data); refits nothing.

S13 <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
S14 <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05"
FIG <- "docs/progress_book/figures/2026-11-06"
dir.create(FIG, showWarnings = FALSE)

# S-13: copy the summarizer's weight-curve panel ----
file.copy(file.path(S13, "weight_curves.png"), file.path(FIG, "s13-weight-curves.png"), overwrite = TRUE)

# S-14: per-dataset validation C, supervised (w_S = 1) against two-step ----
m <- read.csv(file.path(S14, "metrics.csv"))
w <- reshape(m[m$arm %in% c("twostep", "A1_w1"), c("strength", "share", "seed", "arm", "c_validation", "oracle_c")],
             idvar = c("strength", "share", "seed", "oracle_c"), timevar = "arm", direction = "wide")
png(file.path(FIG, "s14-paired-c.png"), width = 1500, height = 700, res = 140)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
for (s in c("strong", "realistic")) {
  x <- w[w$strength == s, ]
  lim <- range(c(x$c_validation.twostep, x$c_validation.A1_w1), na.rm = TRUE)
  cols <- c("0.05%" = "#d95f02", "0.1%" = "#7570b3", "0.2%" = "#1b9e77")
  lab <- paste0(x$share * 100, "%")
  plot(x$c_validation.twostep, x$c_validation.A1_w1, pch = 19, col = cols[lab], xlim = lim, ylim = lim,
       xlab = "Validation C, two-step fit", ylab = "Validation C, supervised fit (w_S = 1)",
       main = sprintf("%s survival (oracle C about %.2f)",
                      if (s == "strong") "Strong" else "Realistic", mean(x$oracle_c)))
  abline(0, 1, lty = 2, col = "grey50")
  legend("bottomright", names(cols), col = cols, pch = 19, title = "Planted share", bty = "n")
}
par(op); dev.off()
cat("Wrote", file.path(FIG, c("s13-weight-curves.png", "s14-paired-c.png")), sep = "\n")
