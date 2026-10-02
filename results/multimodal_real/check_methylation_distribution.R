# ============================================================
# Script: check_methylation_distribution.R
# Purpose: Compare the distribution of screened TCGA methylation values on
#          the beta scale and after asin(2*beta - 1), to inform which scale
#          the Gaussian residual model should use.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: ggplot2; code/load_multiomics_data.R
# Run: Rscript results/multimodal_real/check_methylation_distribution.R
# ============================================================

suppressPackageStartupMessages(library(ggplot2))
source("code/load_multiomics_data.R")

out_dir <- "results/multimodal_real/outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache <- "data/processed_multiomics_tcga_icgc.rds"
d <- if (file.exists(cache)) readRDS(cache) else build_multiomics_cohorts()

beta <- d$training$Y$methylation              # n x p, values in [0, 1]
asin_t <- asin(2 * beta - 1)                  # range [-pi/2, pi/2]

# Per-CpG moments: the Gaussian model is fit per feature (feature-wise tau),
# so normality matters within each CpG across subjects.
#   skewness = E[(x - mu)^3] / sd^3; excess kurtosis = E[(x - mu)^4] / sd^4 - 3
moments <- function(Y) {
  z <- scale(Y)
  data.frame(skew = colMeans(z^3), kurt = colMeans(z^4) - 3)
}
m_beta <- moments(beta); m_asin <- moments(asin_t)
summary_tab <- data.frame(
  scale = c("beta", "asin(2*beta - 1)"),
  median_abs_skew = c(median(abs(m_beta$skew)), median(abs(m_asin$skew))),
  frac_abs_skew_gt_1 = c(mean(abs(m_beta$skew) > 1), mean(abs(m_asin$skew) > 1)),
  median_excess_kurtosis = c(median(m_beta$kurt), median(m_asin$kurt)),
  # Shapiro-Wilk on each CpG (n = 144); fraction rejected at 0.05
  frac_shapiro_reject = c(
    mean(apply(beta, 2, function(x) stats::shapiro.test(x)$p.value) < 0.05),
    mean(apply(asin_t, 2, function(x) stats::shapiro.test(x)$p.value) < 0.05))
)
print(summary_tab, digits = 3)
write.csv(summary_tab, file.path(out_dir, "methylation_distribution_summary.csv"),
          row.names = FALSE)

# Figure: pooled values (all screened CpGs) and per-CpG skewness on each scale
pooled <- rbind(
  data.frame(scale = "beta", value = as.vector(beta)),
  data.frame(scale = "asin(2*beta - 1)", value = as.vector(asin_t)))
skews <- rbind(data.frame(scale = "beta", skew = m_beta$skew),
               data.frame(scale = "asin(2*beta - 1)", skew = m_asin$skew))
p1 <- ggplot(pooled, aes(value)) + geom_histogram(bins = 80) +
  facet_wrap(~scale, scales = "free_x") +
  labs(title = "Screened TCGA CpGs (10k, n = 144): pooled values", x = NULL, y = "count")
p2 <- ggplot(skews, aes(skew)) + geom_histogram(bins = 80) +
  facet_wrap(~scale) + geom_vline(xintercept = c(-1, 1), linetype = 2) +
  labs(title = "Per-CpG skewness across subjects", x = "skewness", y = "CpGs")
ggsave(file.path(out_dir, "methylation_distribution_pooled.png"), p1, width = 8, height = 3.5, dpi = 150)
ggsave(file.path(out_dir, "methylation_distribution_skew.png"), p2, width = 8, height = 3.5, dpi = 150)
