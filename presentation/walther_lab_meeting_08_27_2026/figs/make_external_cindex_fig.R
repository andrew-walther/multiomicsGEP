# ============================================================
# Script:  presentation/walther_lab_meeting_08_27_2026/figs/make_external_cindex_fig.R
# Purpose: Build the 2-way external C-index comparison figure for the 8/27 lab
#          meeting deck: per-cohort grouped bars for the two model arms
#            - Projection predictor    (YF)β   (supervised, recommended)
#            - Unsupervised two-step   EBMF -> Cox, K=40, LASSO stage 2
#                (the fairest independent baseline: large, YFB-uninformed K,
#                 with a regularized stage-2 Cox rather than a plain coxph()
#                 that overfits at K=40; DECISIONS.md 2026-08-20)
#          across the 5 held-out PDAC cohorts (+ a "Mean" group), plus the
#          pooled bootstrap difference CI.
#
#          Legend labels models by CHARACTERISTICS (linear predictor), never by
#          internal code names. Numbers are read from the canonical CSVs.
#
#   Inputs:
#     results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv (D4 = YFB)
#     results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_results.csv (k40 = LASSO stage 2)
#     results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds (D4 per-cohort risk)
#     results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_riskscores.rds (k40 per-cohort risk)
#     code/concordance_ci.R (bootstrap_concordance_diff_ci)
#   Output:
#     presentation/walther_lab_meeting_08_27_2026/assets/external_cindex.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Dependencies: base R (grDevices, graphics), survival (via concordance_ci.R)
# Usage:   Rscript presentation/walther_lab_meeting_08_27_2026/figs/make_external_cindex_fig.R
# ============================================================

# Resolve project root so the script runs from anywhere.
if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else {
  stop("Cannot locate project root (code/fit_modular.R).")
}

DESURV_CSV     <- file.path(ROOT, "results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv")
EBMF_REG_CSV   <- file.path(ROOT, "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_results.csv")
D4_RISK_RDS    <- file.path(ROOT, "results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds")
EBMF_RISK_RDS  <- file.path(ROOT, "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_riskscores.rds")
CI_HELPERS     <- file.path(ROOT, "code/concordance_ci.R")
OUT_PNG        <- file.path(ROOT, "presentation/walther_lab_meeting_08_27_2026/assets/external_cindex.png")

for (f in c(DESURV_CSV, EBMF_REG_CSV, D4_RISK_RDS, EBMF_RISK_RDS, CI_HELPERS)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}
source(CI_HELPERS)

desurv   <- read.csv(DESURV_CSV, stringsAsFactors = FALSE)
ebmf_reg <- read.csv(EBMF_REG_CSV, stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# 1. Assemble a cohort x arm matrix of external C-index values.
# --------------------------------------------------------------------------

cohort_levels <- c("Dijk", "Moffitt_GEO_array", "PACA_AU_array", "PACA_AU_seq", "Puleo_array")
cohort_pretty <- c("Dijk", "Moffitt", "PACA-AU\n(array)", "PACA-AU\n(seq)", "Puleo")

yfb  <- desurv[desurv$model == "D4", ]         # Projection predictor ((YF)beta), DeSurv-aligned
ebmf <- ebmf_reg[ebmf_reg$K == 40, ]           # EBMF -> Cox, K=40, LASSO stage 2 (fairest baseline)

pick <- function(df, ch) {
  v <- df$c_index[df$cohort == ch]
  if (length(v) != 1) stop(sprintf("Expected exactly one C for cohort %s, got %d.", ch, length(v)))
  v
}

c_yfb  <- vapply(cohort_levels, function(ch) pick(yfb,  ch), numeric(1))
c_ebmf <- vapply(cohort_levels, function(ch) pick(ebmf, ch), numeric(1))

mat <- rbind(
  Projection = c(c_yfb,  mean(c_yfb)),
  EBMF       = c(c_ebmf, mean(c_ebmf))
)
colnames(mat) <- c(cohort_pretty, "Mean")

# --------------------------------------------------------------------------
# 2. Pooled bootstrap difference CI (YFB - EBMF k40 LASSO), for the caption.
#    Reproduces DECISIONS.md 2026-08-20's "+0.026, 95% CI [0.0002, 0.0498]".
# --------------------------------------------------------------------------

d4_risk   <- readRDS(D4_RISK_RDS)$D4
ebmf_risk <- readRDS(EBMF_RISK_RDS)$k40

pooled_yfb <- pooled_ebmf <- pooled_time <- pooled_status <- c()
for (coh in cohort_levels) {
  if (!identical(d4_risk[[coh]]$time, ebmf_risk[[coh]]$time) ||
      !identical(d4_risk[[coh]]$status, ebmf_risk[[coh]]$status)) {
    stop("D4 and EBMF k40 risk-score lists are not patient-aligned for cohort: ", coh)
  }
  pooled_yfb    <- c(pooled_yfb, d4_risk[[coh]]$risk)
  pooled_ebmf   <- c(pooled_ebmf, ebmf_risk[[coh]]$risk)
  pooled_time   <- c(pooled_time, d4_risk[[coh]]$time)
  pooled_status <- c(pooled_status, d4_risk[[coh]]$status)
}

ci_pooled <- bootstrap_concordance_diff_ci(pooled_yfb, pooled_ebmf, pooled_time, pooled_status,
                                            B = 2000, seed = 1)
cat(sprintf("Pooled diff (YFB - EBMF k40 LASSO): %.4f, 95%% CI [%.4f, %.4f], n=%d, significant=%s\n",
            ci_pooled$estimate, ci_pooled$lower, ci_pooled$upper, length(pooled_time), ci_pooled$significant))

# Persist for the .qmd to read without recomputing the bootstrap at render time.
ci_out <- data.frame(
  diff_estimate = round(ci_pooled$estimate, 4),
  diff_lower    = round(ci_pooled$lower, 4),
  diff_upper    = round(ci_pooled$upper, 4),
  significant   = ci_pooled$significant,
  n             = length(pooled_time),
  n_events      = sum(pooled_status)
)
CI_OUT_CSV <- file.path(ROOT, "presentation/walther_lab_meeting_08_27_2026/assets/yfb_vs_ebmf_k40_pooled_ci.csv")
write.csv(ci_out, CI_OUT_CSV, row.names = FALSE)
cat(sprintf("Wrote %s\n", normalizePath(CI_OUT_CSV)))

# --------------------------------------------------------------------------
# 3. Plot - grouped bar chart, UNC palette, legend by characteristics.
# --------------------------------------------------------------------------

col_yfb  <- "#4B9CD3"   # Carolina blue (recommended, emphasized)
col_ebmf <- "#9AA7B4"   # muted gray (unsupervised baseline)
bar_cols <- c(col_yfb, col_ebmf)

leg_labels <- c(
  expression("Projection " * eta == (YF) * beta * "  (recommended)"),
  "Unsupervised EBMF -> Cox (K=40, LASSO stage 2)"
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

abline(h = 0.5, lty = 2, col = "gray55", lwd = 1.6)

vals <- as.vector(mat)
text(as.vector(bp), vals + 0.006, labels = sprintf("%.2f", vals),
     cex = 0.85, col = "#13294B", srt = 0)

n_groups <- ncol(mat)
sep_x <- (bp[2, n_groups - 1] + bp[1, n_groups]) / 2
abline(v = sep_x, lty = 3, col = "gray70", lwd = 1.2)

legend("topright", legend = leg_labels, fill = bar_cols, border = NA,
       bty = "n", cex = 1.05, inset = c(0.0, -0.02))

mtext(sprintf("Pooled bootstrap difference: +%.3f, 95%% CI [%.4f, %.4f]%s",
              ci_pooled$estimate, ci_pooled$lower, ci_pooled$upper,
              if (ci_pooled$significant) " (significant)" else ""),
      side = 1, line = 4, cex = 0.95, col = "#13294B")

par(op)
dev.off()

cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
cat(sprintf("  Mean external C-index — Projection: %.3f | EBMF->Cox (K=40, LASSO): %.3f\n",
            mean(c_yfb), mean(c_ebmf)))
