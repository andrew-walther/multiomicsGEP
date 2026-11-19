# ============================================================
# Script: docs/progress_book/figures/plot_2026-11-20.R
# Purpose: Build the November 20 progress-book figure comparing validation C
#          of the supervised and two-step fits in the near-threshold simulation
#          at K = 5 across four survival-signal settings.
# Author: Andrew Walther
# Created: 2026-11-20
# Dependencies: ggplot2
# Run (from the repo root): Rscript docs/progress_book/figures/plot_2026-11-20.R
# ============================================================
# No model fitting and no patient-level data. The inputs are the committed
# aggregate tables of the 2026-11-19 rerun
# (docs/plans/2026-11-19-near-threshold-k5-rerun-plan.md): per-dataset metrics
# and the loading similarity between the two fits of each dataset.
# Each panel is one survival setting. A dot is one simulated dataset: its
# two-step validation C (x) against its supervised validation C (y). The
# corner table gives, for each variance share, the mean paired C gain, the
# recovery counts and the median loading similarity between the two fits.

suppressPackageStartupMessages(library(ggplot2))

# Setup ----

run_dir <- "results/benchmark_sim/outputs/yfb_near_threshold_k5_2026-11-19"
out_dir <- "docs/progress_book/figures/2026-11-20"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Categorical slots 1-3 of the validated default palette (light surface).
# Aqua is below 3:1 contrast, so share is also encoded by marker shape and
# named in the panel tables.
share_cols <- c("0.05%" = "#2a78d6", "0.1%" = "#eb6834", "0.2%" = "#1baf7a")
share_shapes <- c("0.05%" = 16, "0.1%" = 17, "0.2%" = 15)
ink_muted <- "#52514e"

# Survival settings, weakest to strongest, with the true risk score SD that
# defines each (the very strong setting keeps the generator's coefficient).
settings <- c(real_like = "Real-data-like survival signal (SD 0.8",
              realistic = "Moderate survival signal (SD 1.1",
              middle = "Strong survival signal (SD 2.5",
              strong = "Very strong survival signal (SD 6 to 10")

# Data: one row per dataset, supervised (w_S = 1) beside two-step ----

m <- read.csv(file.path(run_dir, "metrics.csv"))
sim <- read.csv(file.path(run_dir, "loading_similarity.csv"))
m <- m[m$arm %in% c("twostep", "A1_w1"), ]
w <- reshape(m[, c("strength", "share", "seed", "arm", "c_validation", "recovered", "oracle_c")],
             idvar = c("strength", "share", "seed", "oracle_c"), timevar = "arm", direction = "wide")
w <- merge(w, sim, by = c("seed", "share", "strength"))
# Fail loud: every dataset must have both fits, a defined C and a similarity.
stopifnot(nrow(w) == 120, !anyNA(w$c_validation.A1_w1), !anyNA(w$c_validation.twostep),
          !anyNA(w$sim_mean))
w$share_lab <- factor(paste0(w$share * 100, "%"), levels = names(share_cols))

# Panel labels add the true-risk C: the concordance of the simulated true risk
# score on the validation patients, the best C any model could reach.
oracle <- tapply(w$oracle_c, w$strength, mean)
panel_lab <- setNames(sprintf("%s, true-risk C %.2f)", settings, oracle[names(settings)]), names(settings))
w$panel <- factor(panel_lab[w$strength], levels = panel_lab)

# Per-cell summary ----
# gain: mean paired C (supervised minus two-step); rec_*: datasets of 10 where
# each fit recovered the prognostic program; sim: median over datasets of the
# mean best loading correlation between supervised and two-step programs.

cells <- split(w, list(w$strength, w$share_lab), drop = TRUE)
summ <- do.call(rbind, lapply(cells, function(x) data.frame(
  strength = x$strength[1], share_lab = x$share_lab[1],
  gain = mean(x$c_validation.A1_w1 - x$c_validation.twostep),
  rec_sup = sum(x$recovered.A1_w1), rec_two = sum(x$recovered.twostep),
  sim = median(x$sim_mean), sim_lo = min(x$sim_min), n = nrow(x))))
summ <- summ[order(match(summ$strength, names(settings)), summ$share_lab), ]
write.csv(summ, file.path(run_dir, "cell_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)
txt <- do.call(rbind, lapply(split(summ, summ$strength), function(s) data.frame(
  panel = factor(panel_lab[s$strength[1]], levels = panel_lab),
  label = paste0("Share   C gain   Recovered   Loading r\n",
                 paste(sprintf("%-5s   %+.3f   %2d vs %-2d    %.3f", s$share_lab, s$gain,
                               s$rec_sup, s$rec_two, s$sim), collapse = "\n")))))

# Figure ----

lim <- c(0.45, 0.97)
upper <- data.frame(x = c(lim[1], lim[1], lim[2]), y = c(lim[1], lim[2], lim[2]))

p <- ggplot(w, aes(c_validation.twostep, c_validation.A1_w1)) +
  geom_polygon(data = upper, aes(x, y), fill = "#f1f0ec", inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = ink_muted, linewidth = 0.4) +
  annotate("text", x = lim[1] + 0.015, y = lim[2] - 0.015, label = "Supervised fit higher",
           hjust = 0, vjust = 1, size = 3.3, colour = ink_muted, fontface = "italic") +
  annotate("text", x = lim[2] - 0.015, y = lim[1] + 0.015, label = "Two-step fit higher",
           hjust = 1, vjust = 0, size = 3.3, colour = ink_muted, fontface = "italic") +
  geom_point(aes(colour = share_lab, shape = share_lab), size = 2.4, alpha = 0.9) +
  geom_text(data = txt, aes(x = lim[2] - 0.015, y = lim[1] + 0.055, label = label),
            inherit.aes = FALSE, hjust = 1, vjust = 0, size = 3, colour = "#0b0b0b",
            lineheight = 1.05, family = "mono") +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(values = share_cols, name = "Prognostic program's share\nof expression variance") +
  scale_shape_manual(values = share_shapes, name = "Prognostic program's share\nof expression variance") +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = "Validation C, two-step fit (unsupervised factorization, then Cox)",
       y = "Validation C, supervised SSBMF (w_S = 1)",
       title = "Supervision changes the fit only when the survival signal is strong",
       subtitle = paste("One point per simulated dataset, K = 5. Corner tables: mean gain in validation C, datasets (of 10)",
                        "where each fit recovered the prognostic program, and median loading correlation between the two fits.",
                        sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e1"),
        strip.text = element_text(face = "bold", size = 10.5), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = ink_muted, size = 10),
        plot.background = element_rect(fill = "white", colour = NA))

out <- file.path(out_dir, "near-threshold-paired-c.png")
ggsave(out, p, width = 11, height = 12, dpi = 150)
cat("Wrote", out, "\n")
