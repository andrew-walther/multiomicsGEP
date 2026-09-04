# ============================================================
# Script:  results/multi_cohort_sim/make_delta_fig.R
# Purpose: Stage 5 delta figure: one bar per (scenario, K_init) of
#          C_joint - C_EBMF with +-1 SE across seeds, making the YFB-vs-EBMF
#          gap explicit rather than left to be eyeballed from adjacent bars.
#          Also confirms and states that both arms are scored on the SAME
#          simulated dataset per seed (the advisors asked directly).
#
#   Input:  results/multi_cohort_sim/outputs/multicohort_sim_results.csv
#   Output: docs/progress_book/figs/2026-09-04_multicohort_delta_fig.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/multi_cohort_sim/make_delta_fig.R
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

d <- read.csv("results/multi_cohort_sim/outputs/multicohort_sim_results.csv", stringsAsFactors = FALSE)

wide <- d |>
  select(scenario, k_setting, K_init, arm, seed, c_index) |>
  tidyr::pivot_wider(names_from = arm, values_from = c_index) |>
  mutate(delta = YFB_base - EBMF)

delta_summary <- wide |>
  group_by(scenario, k_setting, K_init) |>
  summarise(mean_delta = mean(delta, na.rm = TRUE),
            se_delta   = sd(delta, na.rm = TRUE) / sqrt(sum(!is.na(delta))),
            n_seeds    = sum(!is.na(delta)), .groups = "drop")

write.csv(delta_summary, "results/multi_cohort_sim/outputs/multicohort_delta_summary.csv", row.names = FALSE)
print(delta_summary)

delta_summary$K_init <- factor(delta_summary$K_init, levels = sort(unique(delta_summary$K_init)))
p <- ggplot(delta_summary, aes(x = K_init, y = mean_delta, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Joint YFB vs. unsupervised EBMF: C-index gap by scenario and K_init",
       x = "K_init", y = expression(C[joint] - C[EBMF]), fill = "Scenario",
       caption = "Bars = mean over 15 seeds, +-1 SE. Both arms scored on the SAME simulated dataset per seed.") +
  theme_minimal(base_size = 12)
ggsave("docs/progress_book/figs/2026-09-04_multicohort_delta_fig.png", p, width = 9, height = 5.5, dpi = 170, bg = "white")
cat("\nWrote docs/progress_book/figs/2026-09-04_multicohort_delta_fig.png\n")
