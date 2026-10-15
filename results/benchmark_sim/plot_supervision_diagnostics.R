# ============================================================
# Script: plot_supervision_diagnostics.R
# Purpose: Figure for the C1 diagnostic note: survival's share of the F-update
#          precision and E[beta] across iterations, by arm, for the estimated
#          program closest to the planted low-variance prognostic program.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: ggplot2; outputs of run_supervision_diagnostics.R
# Run: Rscript results/benchmark_sim/plot_supervision_diagnostics.R <arm_set> [arm,arm,...]
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

args    <- commandArgs(trailingOnly = TRUE)
arm_set <- if (length(args) >= 1) args[1] else "c1"
in_dir  <- file.path("results/benchmark_sim/outputs/supervision_diagnostics", arm_set)
res     <- read.csv(file.path(in_dir, "results.csv"))
arms    <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
             setdiff(unique(res$arm), c("ref_point_exponential", "canon_wS_0"))

# One row per (seed, arm, iteration) for the program matched to the planted one
long <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
  r <- res[i, ]
  if (!r$arm %in% arms || !is.na(r$error)) return(NULL)
  tr <- readRDS(file.path(in_dir, sprintf("trace_seed%d_noise%s_%s.rds", r$seed, r$noise, r$arm)))$trace
  k  <- r$match_k
  it <- seq_len(nrow(tr$beta_mean))
  rbind(
    data.frame(seed = r$seed, noise = r$noise, arm = r$arm, iter = it,
               metric = "Survival share of F precision (median over genes)",
               value = pmax(tr$shareA_median[, k], 1e-6)),
    data.frame(seed = r$seed, noise = r$noise, arm = r$arm, iter = it,
               metric = "|E[beta]| of the closest program",
               value = pmax(abs(tr$beta_mean[, k]), 1e-6)))
}))
long$arm <- factor(long$arm, levels = arms)

# Categorical slots 1-3 of the reference palette (validated: CVD dE 9.2)
cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")[seq_along(arms)], arms)

p <- ggplot(long, aes(iter, value, colour = arm, group = interaction(seed, arm))) +
  geom_line(linewidth = 0.5, alpha = 0.75) +
  scale_y_log10() +
  scale_colour_manual(values = cols, name = NULL) +
  facet_grid(metric ~ paste("noise", noise), scales = "free_y",
             labeller = label_wrap_gen(width = 28)) +
  labs(x = "CAVI iteration", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold"))

out <- file.path(in_dir, "supervision_trace.png")
ggsave(out, p, width = if (length(unique(long$noise)) > 1) 9 else 7, height = 5.5,
       dpi = 150, bg = "white")
cat("Wrote", out, "\n")
