# ============================================================
# Script: make_pruning_rank_figure.R
# Purpose: Redraw the multimodal simulation's final-rank figure for the
#          writeup from the saved sweep output (no fitting).
# Author: Andrew Walther
# Created: 2026-10-14
# Dependencies: ggplot2; results/benchmark_sim/outputs/multimodal_yfb_pruning/
#               pruning_comparison.csv
# Run (from the repository root):
#   Rscript paper/writeup/figures/make_pruning_rank_figure.R
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

# Inputs ----
sweep <- read.csv("results/benchmark_sim/outputs/multimodal_yfb_pruning/pruning_comparison.csv")
stopifnot(nrow(sweep) == 486, all(is.na(sweep$error)))

# The original model has no removal step, so its final rank equals the
# starting rank by construction. It is checked here and stated in the
# caption rather than drawn (a y = x line carries no information).
original <- sweep[sweep$variant == "original", ]
stopifnot(all(original$K_final == original$K_init))

# The variance-corrected variant is reported in the caption, not drawn.
pruned <- sweep[sweep$variant == "pruned", ]

## Summaries over the three seeds per setting ----
summ <- aggregate(K_final ~ scenario + noise_scale + K_init, data = pruned,
                  FUN = function(x) c(mean = mean(x), lo = min(x), hi = max(x)))
summ <- data.frame(summ[, 1:3], summ$K_final)

scenario_labels <- c(both_informative = "Two prognostic programs",
                     adverse_protective = "Adverse + protective",
                     low_variance_prognostic_factor = "Low-variance prognostic program")
summ$scenario <- factor(scenario_labels[summ$scenario], levels = scenario_labels)
noise_labels <- c("0.25" = "Moderate noise (noise SD x 0.25)",
                  "1" = "September design noise (noise SD x 1)")
summ$noise <- factor(noise_labels[as.character(summ$noise_scale)], levels = noise_labels)

k_grid <- sort(unique(sweep$K_init))

# Plot ----
# Categorical slots 1-3 of the reference palette (validated all-pairs);
# marker shape is a second encoding because slot 3 is below 3:1 contrast.
cols <- c("#2a78d6", "#eb6834", "#1baf7a")
ink <- "#52514e"

p <- ggplot(summ, aes(K_init, mean, colour = scenario, shape = scenario)) +
  geom_hline(yintercept = 3, linetype = 2, colour = ink, linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6),
                 linewidth = 0.5, alpha = 0.8) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.6), size = 2.3) +
  facet_wrap(~ noise, nrow = 1) +
  scale_colour_manual(values = cols, name = "Pruned framework, scenario:") +
  scale_shape_manual(values = c(16, 15, 17), name = "Pruned framework, scenario:") +
  scale_x_continuous(breaks = k_grid) +
  scale_y_continuous(breaks = 0:5, labels = c("0", "1", "2", "3 = true", "4", "5"),
                     limits = c(0, 5.2), expand = c(0, 0)) +
  labs(x = "Starting rank", y = "Final number of programs",
       title = "Pruned fits stay near the true rank at moderate noise; at September noise they under-count",
       subtitle = "Mean final rank over 3 seeds (bars: seed range). Dashed line: true rank 3.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top", legend.justification = "left",
        plot.title = element_text(size = 12, face = "bold", colour = "#0b0b0b"),
        plot.subtitle = element_text(size = 10, colour = ink),
        strip.text = element_text(size = 10, face = "bold", colour = "#0b0b0b"),
        axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
        panel.spacing = unit(1.2, "lines"),
        plot.background = element_rect(fill = "white", colour = NA))

out <- "paper/writeup/figures/multimodal_pruning_final_rank.png"
ggsave(out, p, width = 9, height = 4.4, dpi = 200)
cat("Wrote", out, "\n")
