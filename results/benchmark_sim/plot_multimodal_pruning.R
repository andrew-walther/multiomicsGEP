# ============================================================
# Script: plot_multimodal_pruning.R
# Purpose: Figures for the multimodal pruning comparison: final factor count
#          vs K_init, and factor / prognostic-program recovery by variant.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: ggplot2
# Run: Rscript results/benchmark_sim/plot_multimodal_pruning.R
# ============================================================

suppressPackageStartupMessages(library(ggplot2))
out_dir <- "results/benchmark_sim/outputs/multimodal_yfb_pruning"
res <- read.csv(file.path(out_dir, "pruning_comparison.csv"))
res <- res[is.na(res$error), ]

variant_labels <- c(original = "Original (no intercept, no pruning)",
                    pruned = "Pruned (partial log-lik)",
                    pruned_vc = "Pruned (variance-corrected)")
scenario_labels <- c(both_informative = "Two prognostic programs",
                     adverse_protective = "Adverse + protective",
                     low_variance_prognostic_factor = "Low-variance prognostic")
# Categorical slots 1-4 of the validated reference palette (fixed order);
# shapes give a second encoding because two hues are below 3:1 contrast
cols <- c(original = "#2a78d6", pruned = "#eb6834", pruned_vc = "#1baf7a")
shapes <- c(original = 16, pruned = 15, pruned_vc = 17)
res$variant <- factor(res$variant, levels = names(variant_labels))
res$scenario <- factor(res$scenario, levels = names(scenario_labels), labels = scenario_labels)
res$noise <- factor(ifelse(res$noise_scale == 1, "Original noise (factor PVE < 1%)",
                           "Moderate noise (factor PVE ~ 10%)"),
                    levels = c("Moderate noise (factor PVE ~ 10%)", "Original noise (factor PVE < 1%)"))

agg <- function(y) {
  a <- aggregate(res[[y]], res[c("noise", "scenario", "variant", "K_init")],
                 function(v) c(mean = mean(v), lo = min(v), hi = max(v)))
  cbind(a[1:4], as.data.frame(a$x))
}
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                  legend.direction = "vertical", strip.text = element_text(face = "bold")))
scales_v <- list(scale_colour_manual(values = cols, labels = variant_labels, name = NULL),
                 scale_shape_manual(values = shapes, labels = variant_labels, name = NULL))

k <- agg("K_final")
p1 <- ggplot(k, aes(K_init, mean, colour = variant, shape = variant)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  geom_hline(yintercept = 3, linetype = 2, colour = "grey40") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6,
                 position = position_dodge(width = 1.2)) +
  geom_line(linewidth = 0.7, position = position_dodge(width = 1.2)) +
  geom_point(size = 2.6, position = position_dodge(width = 1.2)) +
  facet_grid(noise ~ scenario) + scales_v +
  scale_x_continuous(breaks = c(3:10, 15)) +
  labs(x = expression(K[init]), y = "Final number of factors (mean, range over seeds)",
       title = "Final factor count vs starting K (true K = 3, dashed; y = x, dotted)")
ggsave(file.path(out_dir, "pruning_K_final.png"), p1, width = 9.5, height = 6.5, dpi = 150)

long <- rbind(
  transform(agg("factors_recovered"), metric = "True factors recovered (of 3; |cor| >= 0.7)"),
  transform(agg("prognostic_recall"), metric = "Prognostic programs found survival-active (share)"),
  transform(agg("false_survival"), metric = "False survival-active factors (count)"))
mod <- long[long$noise == levels(long$noise)[1], ]
p2 <- ggplot(mod, aes(K_init, mean, colour = variant, shape = variant)) +
  geom_line(linewidth = 0.7, position = position_dodge(width = 1.2)) +
  geom_point(size = 2.6, position = position_dodge(width = 1.2)) +
  facet_grid(metric ~ scenario, scales = "free_y", labeller = label_wrap_gen(28)) + scales_v +
  scale_x_continuous(breaks = c(3:10, 15)) +
  labs(x = expression(K[init]), y = NULL,
       title = "Recovery at moderate noise (mean over seeds)")
ggsave(file.path(out_dir, "pruning_recovery.png"), p2, width = 9.5, height = 7.5, dpi = 150)
cat("Figures written to", out_dir, "\n")
