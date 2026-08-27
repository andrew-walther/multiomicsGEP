# ============================================================
# Script:  presentation/walther_lab_meeting_06_18_2026/figs/make_synthetic_figs.R
# Purpose: Build the multi-cohort SIMULATION figures for the 6/18 lab-meeting
#          deck directly from the canonical simulation outputs, so every panel
#          is reproducible and the unsupervised EBMF -> Cox arm is shown
#          alongside the supervised arms wherever a comparative claim is made.
#
#          Figures written into assets/:
#            syn_recovery.png         - factor-recovery correlation by arm/scenario
#            syn_specacc.png          - shared-vs-specific classification accuracy
#            syn_cindex.png           - held-out C-index by arm/scenario (incl. EBMF)
#            syn_recovery_heatmap.png - true vs estimated patient-loading matrix L
#
#          Source figure logic mirrors docs/reports/multicohort_sim_proposal_06_14_26.qmd
#          (the true-vs-estimated L heatmap and the C-index panel are lifted from it),
#          recast as standalone, audience-facing PNGs with no internal config labels.
#
#   Inputs:
#     results/multi_cohort_sim/outputs/multicohort_sim_results.csv
#     results/multi_cohort_sim/outputs/multicohort_sim_example.rds
#   Author:  Claude Code (reviewed by Andrew Walther)
#   Created: 2026-06-15
#   Dependencies: ggplot2, dplyr, tidyr, pheatmap, gridExtra, grid
#   Usage:   Rscript presentation/walther_lab_meeting_06_18_2026/figs/make_synthetic_figs.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

# Resolve project root so the script runs from anywhere.
if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

ASSETS  <- file.path(ROOT, "presentation/walther_lab_meeting_06_18_2026/assets")
OUT_MC  <- file.path(ROOT, "results/multi_cohort_sim/outputs")
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

res <- read.csv(file.path(OUT_MC, "multicohort_sim_results.csv"), stringsAsFactors = FALSE)
ex  <- readRDS(file.path(OUT_MC, "multicohort_sim_example.rds"))

# Audience-facing arm labels (characteristics, no internal codes) and a fixed
# colour for each arm shared across all synthetic figures.
ARM_DISP <- c(YFB_base = "Projection (YF)β", LB_base = "Loadings Lβ",
              YFB_cohort = "Projection +cohort", LB_cohort = "Loadings +cohort",
              EBMF = "EBMF → Cox (unsupervised)")
ARM_LEVELS <- unname(ARM_DISP)
ARM_COLS <- setNames(
  c("#4B9CD3", "#13294B", "#9ecae1", "#6b7d8f", "#d62728"), ARM_LEVELS)
SCEN_LEVELS <- c("all_shared", "hybrid", "nothing_shared")
SCEN_DISP   <- c(all_shared = "All shared",
                 hybrid = "Shared & cohort-specific",
                 nothing_shared = "Nothing shared (control)")

res <- res %>%
  filter(arm %in% names(ARM_DISP)) %>%
  mutate(arm_lab  = factor(ARM_DISP[arm], levels = ARM_LEVELS),
         scen_lab = factor(SCEN_DISP[scenario], levels = unname(SCEN_DISP[SCEN_LEVELS])))

# --------------------------------------------------------------------------
# 1. Factor-recovery correlation: shared and study-specific factors by arm.
#    rec_shared / rec_specific = mean best |cor| between estimated and true
#    gene programs over the shared / study-specific true factors (5 seeds).
# --------------------------------------------------------------------------
rec_long <- res %>%
  select(scen_lab, arm_lab, rec_shared, rec_specific) %>%
  pivot_longer(c(rec_shared, rec_specific),
               names_to = "ftype", values_to = "rec") %>%
  filter(!is.na(rec)) %>%
  mutate(ftype = recode(ftype, rec_shared = "Shared factors",
                        rec_specific = "Study-specific factors")) %>%
  group_by(scen_lab, arm_lab, ftype) %>%
  summarise(mean = mean(rec), se = sd(rec) / sqrt(n()), .groups = "drop")

p_rec <- ggplot(rec_long, aes(arm_lab, mean, fill = arm_lab)) +
  geom_col(alpha = 0.9, width = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25) +
  facet_grid(ftype ~ scen_lab) +
  scale_fill_manual(values = ARM_COLS, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Factor recovery  (mean best |correlation|)") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        strip.text = element_text(face = "bold"))
ggsave(file.path(ASSETS, "syn_recovery.png"), p_rec, width = 10, height = 6.2, dpi = 130)
cat("Wrote syn_recovery.png\n")

# --------------------------------------------------------------------------
# 2. Specificity classification accuracy: fraction of factors correctly typed
#    as shared vs. study-specific (purely from per-cohort loading energy).
# --------------------------------------------------------------------------
spec_df <- res %>%
  filter(!is.na(spec_acc)) %>%
  group_by(scen_lab, arm_lab) %>%
  summarise(mean = mean(spec_acc), se = sd(spec_acc) / sqrt(n()), .groups = "drop")

p_spec <- ggplot(spec_df, aes(arm_lab, mean, fill = arm_lab)) +
  geom_col(alpha = 0.9, width = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = pmin(1, mean + se)), width = 0.25) +
  facet_wrap(~ scen_lab) +
  scale_fill_manual(values = ARM_COLS, guide = "none") +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(x = NULL, y = "Specificity-classification accuracy") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        strip.text = element_text(face = "bold"))
ggsave(file.path(ASSETS, "syn_specacc.png"), p_spec, width = 10, height = 4.6, dpi = 130)
cat("Wrote syn_specacc.png\n")

# --------------------------------------------------------------------------
# 3. Held-out C-index by arm and scenario, INCLUDING the EBMF -> Cox baseline
#    (mean +/- SE over 5 seeds; dashed line = chance 0.5).
# --------------------------------------------------------------------------
cidx <- res %>%
  group_by(scen_lab, arm_lab) %>%
  summarise(C = mean(c_index, na.rm = TRUE),
            se = sd(c_index, na.rm = TRUE) / sqrt(sum(!is.na(c_index))), .groups = "drop")

p_cidx <- ggplot(cidx, aes(arm_lab, C, fill = arm_lab)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_col(alpha = 0.9, width = 0.8) +
  geom_errorbar(aes(ymin = C - se, ymax = C + se), width = 0.25) +
  geom_text(aes(label = sprintf("%.2f", C)), vjust = -0.6, size = 3.4) +
  facet_wrap(~ scen_lab) +
  scale_fill_manual(values = ARM_COLS, guide = "none") +
  coord_cartesian(ylim = c(0.4, 0.95)) +
  labs(x = NULL, y = "Held-out C-index") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        strip.text = element_text(face = "bold"))
ggsave(file.path(ASSETS, "syn_cindex.png"), p_cidx, width = 11, height = 5.0, dpi = 130)
cat("Wrote syn_cindex.png\n")

# --------------------------------------------------------------------------
# 4. True vs estimated patient-loading matrix L (hybrid scenario), showing the
#    model recovers the block-sparse shared/specific structure with no indicator.
#    (Lifted from the proposal report's top panel; saved standalone.)
# --------------------------------------------------------------------------
suppressPackageStartupMessages({ library(pheatmap); library(gridExtra); library(grid) })

d   <- ex$hybrid$data
spl <- ex$hybrid$split
cid <- droplevels(d$cohort_id[spl$train_idx])
fit <- ex$hybrid$YFB_base
col_idx <- fit$match$match           # estimated column for each true factor
labs <- d$factor_labels; flab <- paste0("F", seq_along(labs))
ord  <- order(cid); cid_ord <- cid[ord]
gap_row <- which(diff(as.integer(cid_ord)) != 0)
Ltrue <- d$L_true[spl$train_idx, , drop = FALSE][ord, , drop = FALSE]
Lhat  <- fit$EL[ord, col_idx, drop = FALSE]
rownames(Ltrue) <- rownames(Lhat) <- paste0("P", seq_len(nrow(Ltrue)))
colnames(Ltrue) <- colnames(Lhat) <- flab
row_anno <- data.frame(Cohort = cid_ord, row.names = rownames(Ltrue))
cohort_pal <- setNames(c("#2166AC", "#D6604D"), levels(cid_ord))
col_anno <- data.frame(
  Type = factor(ifelse(labs == "shared", "Shared (prognostic)", "Study-specific"),
                levels = c("Shared (prognostic)", "Study-specific")),
  `Abs.beta` = abs(fit$EBeta[col_idx]), row.names = flab, check.names = FALSE)
anno_colors <- list(Cohort = cohort_pal,
                    Type = c("Shared (prognostic)" = "#D62728", "Study-specific" = "#BDBDBD"),
                    `Abs.beta` = c("#FFFFFF", "#08306B"))
heat_pal <- colorRampPalette(c("#FFFFFF", "#FCBBA1", "#FB6A4A", "#A50F15"))(100)
ph_args <- list(cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE,
                show_colnames = TRUE, color = heat_pal, annotation_row = row_anno,
                annotation_col = col_anno, annotation_colors = anno_colors,
                gaps_row = gap_row, border_color = NA, fontsize_col = 11, fontsize = 11,
                legend = TRUE, annotation_legend = TRUE, silent = TRUE)
p_true <- do.call(pheatmap, c(list(Ltrue, main = "True patient loadings L"), ph_args))
p_hat  <- do.call(pheatmap, c(list(Lhat,  main = "Estimated L  (projection model, no cohort indicator)"), ph_args))

png(file.path(ASSETS, "syn_recovery_heatmap.png"), width = 1700, height = 950, res = 135)
grid.arrange(p_true$gtable, p_hat$gtable, ncol = 2)
dev.off()
cat("Wrote syn_recovery_heatmap.png\n")

cat("\nDone — synthetic figures regenerated (EBMF arm included; no internal labels).\n")
