# ============================================================
# Script: docs/progress_book/figures/plot_2026-11-20.R
# Purpose: Build the November 20 progress-book figure comparing validation C
#          of the supervised and two-step fits in the near-threshold simulation.
# Author: Andrew Walther
# Created: 2026-11-20
# Dependencies: ggplot2
# Run (from the repo root): Rscript docs/progress_book/figures/plot_2026-11-20.R
# ============================================================
# No model fitting and no patient-level data. The input is the committed
# per-dataset metrics table of the near-threshold simulation (2026-11-05).
# This reworks the November 6 scatter (figures/2026-11-06/s14-paired-c.png):
# both panels share one axis range so the two survival strengths are
# comparable, the region where supervision wins is shaded, and each panel
# states the mean C gain and planted-program recovery for each variance share.

suppressPackageStartupMessages(library(ggplot2))

# Setup ----

in_csv <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05/metrics.csv"
out_dir <- "docs/progress_book/figures/2026-11-20"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Categorical slots 1-3 of the validated default palette (light surface).
# Aqua is below 3:1 contrast, so share is also encoded by marker shape and
# named in the panel text.
share_cols <- c("0.05%" = "#2a78d6", "0.1%" = "#eb6834", "0.2%" = "#1baf7a")
share_shapes <- c("0.05%" = 16, "0.1%" = 17, "0.2%" = 15)
ink_muted <- "#52514e"

# Data: one row per dataset, supervised (w_S = 1) beside two-step ----

m <- read.csv(in_csv)
m <- m[m$arm %in% c("twostep", "A1_w1"), ]
w <- reshape(m[, c("strength", "share", "seed", "arm", "c_validation", "recovered", "oracle_c")],
             idvar = c("strength", "share", "seed", "oracle_c"), timevar = "arm", direction = "wide")
stopifnot(nrow(w) == 60, !anyNA(w$c_validation.A1_w1), !anyNA(w$c_validation.twostep))
w$share_lab <- factor(paste0(w$share * 100, "%"), levels = names(share_cols))

# Panel labels give the true-risk C: the concordance of the simulated true
# risk score, the best C any model could reach in that setting. The caption
# in the chapter defines the two survival settings.
oracle <- tapply(w$oracle_c, w$strength, mean)
panel_lab <- c(
  strong = sprintf("Very strong survival signal (true-risk C %.2f)", oracle[["strong"]]),
  realistic = sprintf("Moderate survival signal (true-risk C %.2f)", oracle[["realistic"]]))
w$panel <- factor(panel_lab[w$strength], levels = panel_lab)

# Per-cell summary: mean paired C gain and recovery counts out of 10 ----

cells <- split(w, list(w$strength, w$share_lab), drop = TRUE)
summ <- do.call(rbind, lapply(cells, function(x) data.frame(
  strength = x$strength[1], share_lab = x$share_lab[1],
  gain = mean(x$c_validation.A1_w1 - x$c_validation.twostep),
  rec_sup = sum(x$recovered.A1_w1), rec_two = sum(x$recovered.twostep), n = nrow(x))))
summ <- summ[order(summ$strength, summ$share_lab), ]
txt <- do.call(rbind, lapply(split(summ, summ$strength), function(s) data.frame(
  panel = factor(panel_lab[s$strength[1]], levels = panel_lab),
  label = paste0("Share   C gain   Recovered\n",
                 paste(sprintf("%-5s   %+.3f   %2d vs %-2d", s$share_lab, s$gain,
                               s$rec_sup, s$rec_two), collapse = "\n")))))
print(summ, row.names = FALSE, digits = 3)

# Figure ----

lim <- c(0.5, 0.97)
upper <- data.frame(x = c(lim[1], lim[1], lim[2]), y = c(lim[1], lim[2], lim[2]))

p <- ggplot(w, aes(c_validation.twostep, c_validation.A1_w1)) +
  geom_polygon(data = upper, aes(x, y), fill = "#f1f0ec", inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = ink_muted, linewidth = 0.4) +
  annotate("text", x = 0.515, y = 0.955, label = "Supervised fit higher", hjust = 0, vjust = 1,
           size = 3.3, colour = ink_muted, fontface = "italic") +
  annotate("text", x = 0.955, y = 0.515, label = "Two-step fit higher", hjust = 1, vjust = 0,
           size = 3.3, colour = ink_muted, fontface = "italic") +
  geom_point(aes(colour = share_lab, shape = share_lab), size = 2.4, alpha = 0.9) +
  geom_text(data = txt, aes(x = 0.955, y = 0.545, label = label), inherit.aes = FALSE,
            hjust = 1, vjust = 0, size = 3, colour = "#0b0b0b", lineheight = 1.05, family = "mono") +
  facet_wrap(~panel) +
  scale_colour_manual(values = share_cols, name = "Prognostic program's share\nof expression variance") +
  scale_shape_manual(values = share_shapes, name = "Prognostic program's share\nof expression variance") +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = "Validation C, two-step fit (unsupervised factorization, then Cox)",
       y = "Validation C, supervised SSBMF (w_S = 1)",
       title = "Supervision helps most when the prognostic program is weak and the survival signal is strong",
       subtitle = "One point per simulated dataset. Corner tables: mean gain in validation C and datasets (of 10) where each fit recovered the prognostic program.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e1"),
        strip.text = element_text(face = "bold", size = 10.5), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(colour = ink_muted, size = 10),
        plot.background = element_rect(fill = "white", colour = NA))

out <- file.path(out_dir, "near-threshold-paired-c.png")
ggsave(out, p, width = 11, height = 6.4, dpi = 150)
cat("Wrote", out, "\n")
