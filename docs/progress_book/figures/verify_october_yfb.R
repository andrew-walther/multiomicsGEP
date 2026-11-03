# ============================================================
# Script: verify_october_yfb.R
# Purpose: Compare saved October YFB and rank-7 EBMF fits on matched PDAC inputs.
# Created: 2026-11-03
# Dependencies: survival, yaml, local PDAC training inputs, saved fits
# ============================================================
# Run from the repository root. Only aggregate output is written.
# Historical log flags reproduce the benchmark gene set for this comparison.
source("code/concordance_ci.R")
source("code/preprocess_desurv.R")
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
PLATFORM_LOG_TRANSFORM[c("TCGA_PAAD", "Dijk", "PACA_AU_seq")] <- TRUE
train <- cfg$pdac$training_cohorts
raw <- lapply(setNames(train, train), load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
pp <- preprocess_merged_cohorts(
  cohort_raw_list = raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[train],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank")
d4_genes <- readRDS("results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
stopifnot(identical(pp$gene_names, d4_genes))
Y_train <- pp$Y
figures <- "docs/progress_book/figures/2026-11-06"
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
library(survival)
out <- "results/benchmark_sim/outputs"
joint <- readRDS(file.path(out, "supervised_real_check/fits.rds"))$supervised
ebmf <- readRDS(file.path(out, "ebmf_cox_external/ebmf_cox_external_fit_k7.rds"))
joint_risk <- readRDS(file.path(out, "supervised_real_check/riskscores.rds"))
ebmf_risk <- readRDS(file.path(out, "ebmf_cox_external/ebmf_cox_external_riskscores_k7.rds"))
audit <- read.csv(file.path(out, "ebmf_cox_external/2026-10-30-orientation-rescore.csv"))
cohorts <- names(joint_risk)
stopifnot(identical(d4_genes, ebmf$train_genes), setequal(cohorts, names(ebmf_risk)))

# Compare paired predictions with the archived orientations.
differences <- matrix(NA_real_, 2000, length(cohorts))
performance <- matrix(NA_real_, length(cohorts), 2, dimnames = list(cohorts, c("Joint YFB", "EBMF → Cox")))
set.seed(42)
for (j in seq_along(cohorts)) {
  co <- cohorts[j]
  a <- joint_risk[[co]]
  b <- ebmf_risk[[co]]
  stopifnot(identical(a$time, b$time), identical(a$status, b$status))
  performance[j, ] <- c(frozen_reverse_cindex(a$supervised, a$time, a$status),
                       frozen_reverse_cindex(b$risk, b$time, b$status))
  stopifnot(abs(performance[j, 2] - audit$c_frozen_training[audit$arm == "ebmf_cox_k7" & audit$cohort == co]) < 1e-10)
  differences[, j] <- replicate(nrow(differences), {
    i <- sample.int(length(a$time), replace = TRUE)
    frozen_reverse_cindex(a$supervised[i], a$time[i], a$status[i]) -
      frozen_reverse_cindex(b$risk[i], a$time[i], a$status[i])
  })
}
print(performance)
print(colMeans(performance))
print(quantile(rowMeans(differences), c(.025, .975)))
print(apply(differences, 2, quantile, c(.025, .975)))

# Compare centered projected scores on identical genes and training patients.
joint_scores <- Y_train %*% joint$EF
ebmf_scores <- Y_train %*% ebmf$F_ebmf
correspondence <- abs(cor(joint_scores, ebmf_scores))
basis <- qr(scale(ebmf_scores, scale = FALSE))
centered <- scale(joint_scores, scale = FALSE)
span_r2 <- 1 - colSums(qr.resid(basis, centered)^2) / colSums(centered^2)
programs <- subset(read.csv(file.path(out, "supervised_real_check/programs.csv")), arm == "supervised")
active <- programs$program[programs$survival_active]
closest <- max.col(correspondence, ties.method = "first")
loading_match <- cor(joint$EF, ebmf$F_ebmf)[cbind(seq_along(closest), closest)]
print(data.frame(programs, closest_ebmf = closest,
                 best_score_match = apply(correspondence, 1, max), loading_match, span_r2))

# Exploratory overrepresentation in the risk-oriented upper tail.
sets <- readRDS(file.path(out, "pathway_enrichment/pdac_genesets.rds"))
enrichment <- do.call(rbind, lapply(active, function(k) {
  weights <- joint$EF[, k] * sign(joint$EBeta[k])
  top <- d4_genes[order(weights, decreasing = TRUE)[1:100]]
  cat("Program", k, "risk-associated markers:", paste(head(top, 12), collapse = ", "), "\n")
  do.call(rbind, lapply(names(sets), function(set) {
    members <- intersect(sets[[set]], d4_genes)
    overlap <- length(intersect(top, members))
    data.frame(program = k, set, overlap, available = length(members),
               p = phyper(overlap - 1, length(members), length(d4_genes) - length(members),
                          length(top), lower.tail = FALSE))
  }))
}))
enrichment$q <- p.adjust(enrichment$p, "BH")
print(enrichment[order(enrichment$p), ], row.names = FALSE)

# Figures contain aggregate values only.
colors <- c("#12766F", "#704E9B")
png(file.path(figures, "external-concordance.png"), 1450, 850, res = 145)
par(mar = c(5, 10, 4, 2), las = 1)
plot(NA, xlim = c(.50, .72), ylim = c(.4, 5.7), yaxt = "n", xlab = "External concordance index",
     ylab = "", main = "External concordance by cohort")
axis(2, 5:1, c("Dijk", "Moffitt", "PACA-AU array", "PACA-AU RNA-seq", "Puleo"), tick = FALSE)
abline(v = seq(.50, .70, .05), col = "gray90")
segments(performance[, 1], 5:1, performance[, 2], 5:1, col = "gray65", lwd = 2)
points(performance[, 1], 5:1, col = colors[1], pch = 19, cex = 1.4)
points(performance[, 2], 5:1, col = colors[2], pch = 17, cex = 1.4)
legend("bottomright", colnames(performance), pch = c(19, 17), col = colors, bty = "n")
dev.off()

png(file.path(figures, "program-correspondence.png"), 1500, 1050, res = 150)
par(mar = c(6.5, 6, 4, 9), las = 1)
k <- nrow(correspondence)
image(seq_len(k), seq_len(k), t(correspondence[k:1, ]), axes = FALSE, xlab = "Unsupervised EBMF program",
      ylab = "", zlim = c(0, 1), col = hcl.colors(25, "Blues 3", rev = TRUE),
      main = "How closely EBMF matches the joint program scores")
axis(1, seq_len(k), tick = FALSE)
axis(2, seq_len(k), paste0("YFB ", k:1, ifelse(k:1 %in% active, " *", "")), tick = FALSE)
axis(4, seq_len(k), sprintf("%.4f%%", 100 * span_r2[k:1]), tick = FALSE)
text(k + 1.06, k + 0.67, "EBMF span R²", xpd = NA, cex = .88)
for (i in seq_len(k)) for (j in seq_len(k)) {
  value <- correspondence[i, j]
  text(j, k + 1 - i, sprintf("%.2f", value), col = if (value > .6) "white" else "#172C44")
}
mtext("Cells show absolute score correlations. * Survival-active in the joint fit.", side = 1, line = 4.6, cex = .78)
dev.off()
