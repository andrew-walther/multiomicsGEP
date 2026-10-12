# ============================================================
# Script: plot_real_cindex_by_kinit.R
# Purpose: External concordance (ICGC primary PDAC) of the pruned multimodal
#          fit across starting ranks, against the established single-modality
#          model on the same donors.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: ggplot2; outputs/real_sweep_summary.csv
# Run: Rscript results/multimodal_real/plot_real_cindex_by_kinit.R
# ============================================================
suppressPackageStartupMessages(library(ggplot2))
x <- read.csv("results/multimodal_real/outputs/real_sweep_summary.csv")
j <- x[x$method == "joint_multimodal_yfb" & x$config == "pruned" & x$cohort == "icgc_primary", ]
single_modality_c <- 0.685  # D4 on the same 50 donors (DECISIONS.md 2026-10-02)
p <- ggplot(j, aes(K_init, c_index)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey55") +
  geom_hline(yintercept = single_modality_c, linetype = 2, colour = "#eb6834") +
  annotate("text", x = 15.4, y = 0.715, hjust = 1, size = 3, colour = "#a8441c",
           label = "established single-modality model,\nsame 50 donors (0.685)") +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = "#2a78d6", linewidth = 0.6) +
  geom_point(colour = "#2a78d6", size = 2.6) +
  geom_text(aes(label = paste0(n_survival_active, " SA")), nudge_x = 0.4, size = 2.8, colour = "grey30") +
  scale_x_continuous(breaks = sort(unique(j$K_init))) + coord_cartesian(ylim = c(0.35, 0.82)) +
  labs(x = expression(K[init] ~ "(= final K; nothing was pruned)"),
       y = "External C on ICGC primary PDAC (95% CI)",
       title = "Multimodal model on real data: external concordance by starting rank",
       subtitle = "SA = number of survival-active programs") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("results/multimodal_real/outputs/real_cindex_by_kinit.png", p, width = 8, height = 4.5, dpi = 150)
