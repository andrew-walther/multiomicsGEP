# ============================================================
# Script:  presentation/walther_lab_meeting_06_18_2026/figs/make_external_cindex_fig.R
# Purpose: Build the 3-way external C-index comparison figure for the 6/18 lab
#          meeting deck: per-cohort grouped bars for the three model arms
#            - Loadings predictor      Lβ      (supervised, joint)
#            - Projection predictor    (YF)β   (supervised, recommended)
#            - Unsupervised two-step   EBMF → Cox
#          across the 5 held-out PDAC cohorts (+ a "Mean" group).
#
#          Legend labels models by CHARACTERISTICS (linear predictor), never by
#          internal code names. Numbers are read from the canonical CSVs.
#
#   Inputs:
#     results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv
#       (D3 = LB DeSurv-aligned, D4 = YFB DeSurv-aligned)
#     results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_results.csv
#   Output:
#     presentation/walther_lab_meeting_06_18_2026/assets/external_cindex.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-15
# Dependencies: base R only (grDevices, graphics)
# Usage:   Rscript presentation/walther_lab_meeting_06_18_2026/figs/make_external_cindex_fig.R
# ============================================================

# Resolve project root so the script runs from anywhere.
if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else {
  stop("Cannot locate project root (code/fit_modular.R).")
}

DESURV_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv")
EBMF_CSV   <- file.path(ROOT, "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_results.csv")
OUT_PNG    <- file.path(ROOT, "presentation/walther_lab_meeting_06_18_2026/assets/external_cindex.png")

if (!file.exists(DESURV_CSV)) stop("Missing desurv CSV: ", DESURV_CSV)
if (!file.exists(EBMF_CSV))   stop("Missing EBMF→Cox CSV: ", EBMF_CSV,
                                   "\nRun results/benchmark_sim/run_ebmf_cox_external.R first.")

desurv <- read.csv(DESURV_CSV, stringsAsFactors = FALSE)
ebmf   <- read.csv(EBMF_CSV,   stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# 1. Assemble a cohort × arm matrix of external C-index values.
# --------------------------------------------------------------------------

# Canonical cohort order + display labels.
cohort_levels <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array", "PACA_AU_seq", "Puleo_array")
cohort_pretty <- c("Dijk", "Moffitt", "PACA-AU\n(array)", "PACA-AU\n(seq)", "Puleo")

lb  <- desurv[desurv$model == "D3", ]  # Loadings predictor (Lβ), DeSurv-aligned
yfb <- desurv[desurv$model == "D4", ]  # Projection predictor ((YF)β), DeSurv-aligned

pick <- function(df, ch) {
  v <- df$c_index[df$cohort == ch]
  if (length(v) != 1) stop(sprintf("Expected exactly one C for cohort %s, got %d.", ch, length(v)))
  v
}

c_lb   <- vapply(cohort_levels, function(ch) pick(lb,   ch), numeric(1))
c_yfb  <- vapply(cohort_levels, function(ch) pick(yfb,  ch), numeric(1))
c_ebmf <- vapply(cohort_levels, function(ch) pick(ebmf, ch), numeric(1))

# Append the mean group.
mat <- rbind(
  Loadings   = c(c_lb,   mean(c_lb)),
  Projection = c(c_yfb,  mean(c_yfb)),
  EBMF       = c(c_ebmf, mean(c_ebmf))
)
colnames(mat) <- c(cohort_pretty, "Mean")

# --------------------------------------------------------------------------
# 2. Plot — grouped bar chart, UNC palette, legend by characteristics.
# --------------------------------------------------------------------------

# UNC brand colors: navy (loadings), Carolina blue (projection/recommended), gray (two-step).
col_lb   <- "#13294B"   # navy
col_yfb  <- "#4B9CD3"   # Carolina blue (recommended — emphasized)
col_ebmf <- "#9AA7B4"   # muted gray (unsupervised baseline)
bar_cols <- c(col_lb, col_yfb, col_ebmf)

leg_labels <- c(
  expression("Loadings " * eta == L * beta * " (supervised)"),
  expression("Projection " * eta == (YF) * beta * "  (recommended)"),
  "Unsupervised  EBMF → Cox"
)

png(OUT_PNG, width = 1450, height = 1020, res = 150)
op <- par(mar = c(5.4, 5.4, 3.4, 1.6), xpd = NA)

bp <- barplot(
  mat, beside = TRUE, col = bar_cols, border = NA,
  ylim = c(0.5, 0.72), xpd = FALSE,
  ylab = "External C-index", cex.lab = 1.35, cex.axis = 1.1,
  main = "External validation across 5 held-out PDAC cohorts",
  cex.main = 1.45, names.arg = colnames(mat), cex.names = 1.12,
  legend.text = FALSE
)

# Reference line at chance.
abline(h = 0.5, lty = 2, col = "gray55", lwd = 1.6)

# Value labels above each bar.
vals <- as.vector(mat)
text(as.vector(bp), vals + 0.006, labels = sprintf("%.2f", vals),
     cex = 0.82, col = "#13294B", srt = 0)

# Separator before the Mean group for visual emphasis.
n_groups <- ncol(mat)
sep_x <- (bp[3, n_groups - 1] + bp[1, n_groups]) / 2
abline(v = sep_x, lty = 3, col = "gray70", lwd = 1.2)

legend("topright", legend = leg_labels, fill = bar_cols, border = NA,
       bty = "n", cex = 1.05, inset = c(0.0, -0.02))

par(op)
dev.off()

cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
cat(sprintf("  Mean external C-index — Loadings: %.3f | Projection: %.3f | EBMF->Cox: %.3f\n",
            mean(c_lb), mean(c_yfb), mean(c_ebmf)))
