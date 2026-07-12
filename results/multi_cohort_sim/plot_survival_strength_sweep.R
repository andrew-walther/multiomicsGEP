# ============================================================
# Script:  results/multi_cohort_sim/plot_survival_strength_sweep.R
# Purpose: Generate the C-index-vs-strength comparison figure for
#          docs/reports/joint_vs_twostep_sweep_07_12_2026.qmd from the
#          results CSV produced by run_survival_strength_sweep.R.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-12
# Usage:   Rscript results/multi_cohort_sim/plot_survival_strength_sweep.R
# ============================================================

if (file.exists("code/fit_modular.R")) {
  # already at repo root
} else if (file.exists("../../code/fit_modular.R")) {
  setwd("../..")
}

suppressPackageStartupMessages(library(ggplot2))

IN_CSV  <- "results/multi_cohort_sim/outputs/survival_strength_sweep_results.csv"
OUT_PNG <- "results/multi_cohort_sim/outputs/survival_strength_sweep_plot.png"

if (!file.exists(IN_CSV))
  stop(sprintf("%s not found -- run run_survival_strength_sweep.R first.", IN_CSV))

r <- read.csv(IN_CSV)
agg <- aggregate(c_index ~ strength + arm, data = r,
                 FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
agg <- do.call(data.frame, agg)
names(agg)[3:4] <- c("mean_c", "se_c")

# YFB_alpha0 is C=0.5 by construction (beta is forced to the floor at
# alpha=0) -- informative for the internal-control check (see
# survival_strength_sweep_alpha_invariance.csv) but not for this plot, which
# compares predictive performance.
agg <- agg[agg$arm != "YFB_alpha0", ]
agg$arm <- factor(agg$arm,
                  levels = c("YFB_tuned", "PCA_Cox", "EBMF", "LB_tuned"),
                  labels = c("YFB (joint, tuned alpha)", "PCA + Cox (2-step)",
                             "EBMF + Cox (2-step)", "LB (joint, tuned alpha)"))

p <- ggplot(agg, aes(x = strength, y = mean_c, color = arm, group = arm)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_c - se_c, ymax = mean_c + se_c), width = 0.05, alpha = 0.6) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 0.25, 0.5, 1, 2, 4)) +
  labs(x = "Prognostic effect strength (beta_shared multiplier)",
       y = "Mean held-out C-index (oriented)", color = NULL,
       title = "Joint model vs. two-step baselines as survival signal strength grows",
       subtitle = "Error bars: +/- 1 SE across 5 seeds. Dashed line: chance (C=0.5).") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(OUT_PNG, p, width = 8, height = 5.5, dpi = 150)
cat(sprintf("Saved: %s\n", OUT_PNG))
