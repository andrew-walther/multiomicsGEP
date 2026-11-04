# ============================================================
# Script: make_figures.R
# Purpose: Collect and draw every PNG used by the Part III prelim deck
#          (presentation/prelim-slides/figures/), from saved aggregate results only.
# Author: Andrew Walther
# Created: 2026-11-03
# Dependencies: ggplot2
# ============================================================
#
# Run from the multiomicsGEP root:
#   Rscript presentation/prelim-slides/tools/make_figures.R
#
# Two kinds of figure:
#   1. Copies of figures that already have their own scripts elsewhere in the
#      repository (writeup, prelim chapter, progress book). The source path is
#      listed next to each copy so the provenance stays visible; rerun the
#      source script first if its inputs changed.
#   2. Figures drawn here for the deck (the external concordance forest plot).
# No patient-level data is read.

suppressPackageStartupMessages(library(ggplot2))

# Paths ----
out_dir <- "presentation/prelim-slides/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Copied figures ----
# Names on the left are the deck's file names; paths on the right are the
# scripted sources (each made by the script named in the comment).
copies <- c(
  # paper/writeup/figures/make_writeup_figures.R
  "single_sim_delta.png"        = "paper/writeup/figures/single_sim_delta.png",
  "km_risk_external.png"        = "paper/writeup/figures/km_risk_external.png",
  "km_programs_external.png"    = "paper/writeup/figures/km_programs_external.png",
  # paper/writeup/figures/make_real_data_workflows.R
  "real_data_workflow_single.png" = "paper/writeup/figures/real_data_workflow_single.png",
  # paper/prelim/ figure scripts (Chapter 4)
  "programs_heatmap.png"        = "paper/prelim/figures/single_programs_heatmap.png",
  "percohort_beta.png"          = "paper/prelim/figures/single_percohort_beta.png",
  # docs/progress_book/figures/plot_planted_yfb.R, plot_planted_yfb_comparators.R,
  # verify_october_yfb.R (backup slides on supervised fits, 2026-11-03)
  "planted_matched_comparators.png" = "docs/progress_book/figures/2026-11-06/planted-matched-comparators.png",
  "planted_yfb_performance.png" = "docs/progress_book/figures/2026-11-06/planted-yfb-performance.png",
  "joint_external_concordance.png" = "docs/progress_book/figures/2026-11-06/external-concordance.png",
  "joint_program_correspondence.png" = "docs/progress_book/figures/2026-11-06/program-correspondence.png"
)
missing <- copies[!file.exists(copies)]
if (length(missing)) stop("Missing source figures: ", paste(missing, collapse = ", "))
ok <- file.copy(copies, file.path(out_dir, names(copies)), overwrite = TRUE)
if (!all(ok)) stop("Copy failed for: ", paste(names(copies)[!ok], collapse = ", "))

# External concordance forest plot ----
## Inputs ----
# Model: frozen single-modality fit (benchmark snapshot, 2026-10-14), with
# independent-record bootstrap intervals. Comparators: two-step fits
# (factorize, then Cox on projection scores) at rank 7 (unpenalized) and
# 40 factors (LASSO Cox). Comparator intervals were not saved, so they are
# drawn as points only.
cohort_labels <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt",
                   PACA_AU_array = "PACA-AU array", PACA_AU_seq = "PACA-AU RNA-seq",
                   Puleo_array = "Puleo")

model <- read.csv("results/benchmark_snapshot_2026-10-14/single_modality/external_cindex_ci.csv")
model <- subset(model, model == "YFB (recommended)")
k7 <- read.csv("results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_results_k7.csv")
k40 <- subset(read.csv("results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_regularized_results.csv"),
              k_tag == "k40")

pts <- rbind(
  data.frame(cohort = model$cohort, c = model$c_index, lo = model$ci_lo, hi = model$ci_hi,
             arm = "Joint model (survival does not shape loadings)"),
  data.frame(cohort = k7$cohort, c = k7$c_index, lo = NA, hi = NA,
             arm = "Two-step, 7 programs"),
  data.frame(cohort = k40$cohort, c = k40$c_index, lo = NA, hi = NA,
             arm = "Two-step, 40 programs, LASSO Cox")
)
# Fail loud if any cohort is missing from an arm
stopifnot(all(table(pts$arm) == 5), all(pts$cohort %in% names(cohort_labels)))

# Mean C across the five cohorts, per arm (matches writeup Table: 0.627, 0.623, 0.601)
means <- aggregate(c ~ arm, pts, mean)
pts <- rbind(pts, data.frame(cohort = "Mean", c = means$c, lo = NA, hi = NA, arm = means$arm))
pts$label <- factor(ifelse(pts$cohort == "Mean", "Mean of five", cohort_labels[pts$cohort]),
                    levels = rev(c(unname(cohort_labels), "Mean of five")))
pts$arm <- factor(pts$arm, levels = c("Joint model (survival does not shape loadings)",
                                      "Two-step, 7 programs", "Two-step, 40 programs, LASSO Cox"))

## Plot ----
dodge <- position_dodge(width = 0.6)
p <- ggplot(pts, aes(x = c, y = label, colour = arm, shape = arm)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.9,
                 position = dodge, na.rm = TRUE) +
  geom_point(size = 3.6, position = dodge) +
  scale_colour_manual(values = c("#13294B", "#4B9CD3", "#C9A227"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "External concordance (C-index)", y = NULL) +
  coord_cartesian(xlim = c(0.45, 0.8)) +
  theme_minimal(base_size = 26) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
ggsave(file.path(out_dir, "external_c_forest.png"), p, width = 10, height = 8, dpi = 150, bg = "white")

print(round(setNames(means$c, means$arm), 3))
cat("Wrote", length(copies) + 1, "figures to", out_dir, "\n")
