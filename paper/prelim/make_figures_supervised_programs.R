# ============================================================
# Script: paper/prelim/make_figures_supervised_programs.R
# Purpose: Program figure, program table, gene-set enrichment, subtype/stroma
#          concordance and correspondence results for the survival-supervised
#          YFB fit on the training cohorts and its fixed external projections.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: ggplot2, patchwork, fgsea, msigdbr, BiocParallel, yaml
# Run (from the repo root): Rscript paper/prelim/make_figures_supervised_programs.R
# Inputs : saved all-patient fits (w_S = 1 and w_S = 0), the training matrix in
#          the git-ignored data/ folder, the saved PDAC gene sets, and the
#          historical alpha_F = 0 benchmark fit (read only). The subtype analysis
#          uses the local PDAC_DATA_ROOT and saved per-series subtype references.
# Outputs: aggregate CSVs in OUT_DIR, figures in paper/prelim/figures,
#          tables in paper/prelim/tables. No patient-level values are written.
# ============================================================

.libPaths(c("data/round2_R_library", .libPaths()))  # fgsea / msigdbr fall back here
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(fgsea)
  library(msigdbr)
})
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")
source("code/predict_ssbmf_yfb.R")
set.seed(20261106)

## cairo_pdf embeds fonts but needs the cairo shared library (it fails to load without XQuartz's X11 libs).
## capabilities("cairo") can report TRUE when the library cannot load, so test by opening a device.
CAIRO_OK <- tryCatch({ f <- tempfile(fileext = ".pdf"); suppressWarnings(cairo_pdf(f)); dev.off(); file.exists(f) },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (!CAIRO_OK) warning("cairo_pdf unavailable: figures are written with pdf() and their fonts are NOT embedded",
                       call. = FALSE)
PDF_DEVICE <- if (CAIRO_OK) cairo_pdf else function(filename, ...) grDevices::pdf(file = filename, ..., useDingbats = FALSE)

# Settings ----
FIT_SUP   <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits/size_D4_N2218_w1_f0.rds"
FIT_UNSUP <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits/size_D4_N2218_w0_f0.rds"
TRAIN_RDS <- "data/yfb_supervised_real_2026-11-04.rds"          # not in git
GENESETS  <- "results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds"
BENCHMARK <- "results/benchmark_snapshot_2026-10-14/single_modality"
CMB_RDATA <- file.path(PDAC_DATA_ROOT, "original", "cmbSubtypes.RData")
OUT_DIR   <- "results/benchmark_sim/outputs/prelim_ch4_2026-11-06"
FIG_DIR   <- "paper/prelim/figures"
TAB_DIR   <- "paper/prelim/tables"
N_POLE    <- 4L      # genes per pole in the table
N_HEAT    <- 3L      # genes per pole in the print-size heatmap
Z_ACTIVE  <- 1.96    # survival-active: |E[beta]| >= 1.96 posterior SD
for (d in c(OUT_DIR, FIG_DIR, TAB_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Program labels ----
# Training gene-set and signature results motivate these descriptive labels.
LABELS <- c("P1", "Classical / basal-like", "Proliferation / classical",
            "P4", "Mitotic / stromal signals", "P6", "P7")

# Helpers ----

#' Escape LaTeX special characters in a string.
#' @param x Character vector.
#' @return Escaped character vector.
#' @examples
#' tex_esc("10%")
tex_esc <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#$])", "\\\\\\1", x)
  gsub("−", "$-$", x)
}

#' Format a signed number for LaTeX with a true minus sign.
#' @param x Numeric values.
#' @param d Digits after the decimal point.
#' @return Formatted character vector.
#' @examples
#' fmt_signed(-0.25)
fmt_signed <- function(x, d = 2) {
  s <- formatC(x, format = "f", digits = d)
  sub("^-", "$-$", s)
}

#' Top genes at each pole of a loading vector.
#' @param w Named numeric loading vector.
#' @param n Genes per pole.
#' @return List with positive and negative gene names.
#' @examples
#' poles(c(A = 1, B = -1), 1)
poles <- function(w, n) {
  o <- order(w, decreasing = TRUE)
  list(pos = names(w)[o[seq_len(n)]], neg = names(w)[rev(o)[seq_len(n)]])
}

# Load fits and training data ----
fit <- readRDS(FIT_SUP)
fit0 <- readRDS(FIT_UNSUP)
dat <- readRDS(TRAIN_RDS)
genes <- fit$genes
stopifnot(identical(genes, fit0$genes), all(genes %in% dat$genes), length(fit$te) == 0)
EF <- fit$state$EF
EF0 <- fit0$state$EF
rownames(EF) <- rownames(EF0) <- genes
K <- ncol(EF)
stopifnot(K == length(LABELS), !anyDuplicated(genes))
Yfull <- dat$Y
colnames(Yfull) <- dat$genes
Y <- Yfull[, genes, drop = FALSE]
Z <- predict_ssbmf_yfb(Y, EF, fit$state$EBeta, fit$state$EF_norms)$ZF_test
stopifnot(all(abs(apply(Z, 2, sd) - 1) < 1e-8))
cat(sprintf("Fit: %d genes, %d programs, n = %d patients\n", length(genes), K, nrow(Y)))

# Program table (T4 inputs) ----
# The compact fit stores both moments of beta, but does not store VBeta.
beta <- fit$state$EBeta
bsd <- sqrt(pmax(fit$state$EBeta2 - beta^2, 0))
stopifnot(all(bsd > 0))
prog <- data.frame(program = seq_len(K), beta = beta, beta_sd = bsd, z = beta / bsd,
                   short = LABELS)
prog$active <- abs(prog$z) >= Z_ACTIVE
ref <- read.csv("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/programs_D4.csv")
ref <- ref[ref$arm == "w1", ]
stopifnot(identical(ref$program, prog$program), all(abs(ref$beta - prog$beta) < 1e-8),
          all(abs(ref$beta_sd - prog$beta_sd) < 1e-8), identical(ref$survival_active, prog$active))
pl <- lapply(seq_len(K), function(k) poles(EF[, k], N_POLE))
prog$top_pos <- vapply(pl, function(p) paste(p$pos, collapse = ", "), "")
prog$top_neg <- vapply(pl, function(p) paste(p$neg, collapse = ", "), "")
write.csv(prog, file.path(OUT_DIR, "programs_table.csv"), row.names = FALSE)
print(prog[, c("program", "beta", "beta_sd", "z", "active")])

# Enrichment (B2, T5) ----
msig <- rbind(msigdbr(species = "Homo sapiens", collection = "H"),
              msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
              msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
              msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"))
msig_collection <- c(H = "Hallmark", "CP:REACTOME" = "Reactome",
                     "CP:KEGG_MEDICUS" = "KEGG Medicus", "CP:KEGG_LEGACY" = "KEGG Legacy")
msig$collection <- unname(msig_collection[ifelse(msig$gs_collection == "H", "H", msig$gs_subcollection)])
pdac <- readRDS(GENESETS)
sets <- lapply(c(split(msig$gene_symbol, msig$gs_name), pdac), unique)
msig_map <- unique(msig[c("gs_name", "collection")])
stopifnot(!anyDuplicated(msig_map$gs_name), length(unique(msig$db_version)) == 1L)
coll <- c(setNames(msig_map$collection, msig_map$gs_name), setNames(rep("PDAC", length(pdac)), names(pdac)))
stopifnot(!anyDuplicated(names(sets)), !anyNA(coll))
enr <- do.call(rbind, lapply(which(prog$active), function(k) {
  r <- fgseaMultilevel(sets, EF[, k], minSize = 10, maxSize = 500, scoreType = "std",
                       nPermSimple = 10000, BPPARAM = BiocParallel::SerialParam())
  r$padj <- p.adjust(r$pval, "BH", n = nrow(r))
  data.frame(program = k, collection = coll[r$pathway], set = r$pathway, size = r$size,
             NES = r$NES, pval = r$pval, padj = r$padj, log2err = r$log2err,
             msigdb_version = ifelse(coll[r$pathway] == "PDAC", NA_character_, unique(msig$db_version)),
             leading_edge = vapply(r$leadingEdge, function(x) paste(head(x, 8), collapse = ";"), ""))
}))
enr <- enr[order(enr$program, enr$padj), ]
write.csv(enr, file.path(OUT_DIR, "programs_enrichment.csv"), row.names = FALSE)
print(with(enr, table(program, collection, estimable = is.finite(padj))))

# Correspondence (B4, F8) ----
# (a) w_S = 1 versus w_S = 0 loadings on the same genes
Ra <- cor(EF, EF0)
# (b) w_S = 1 loadings versus the alpha_F = 0 benchmark programs on shared genes
old <- readRDS(file.path(BENCHMARK, "desurv_comparison_fits.rds"))$D4$EF
rownames(old) <- readRDS(file.path(BENCHMARK, "d4_gene_names.rds"))
shared <- intersect(genes, rownames(old))
# The pruned benchmark has zero columns. Their correlations are undefined.
live_old <- apply(old[shared, , drop = FALSE], 2, sd) > 0
Rb <- matrix(NA_real_, nrow = K, ncol = ncol(old))
Rb[, live_old] <- cor(EF[shared, , drop = FALSE], old[shared, live_old, drop = FALSE])
cat(sprintf("Shared genes with the alpha_F = 0 benchmark: %d of %d new, %d old\n",
            length(shared), length(genes), nrow(old)))
corr <- rbind(
  data.frame(comparison = "w_S=1 vs w_S=0", new_program = rep(seq_len(K), K),
             other_program = rep(seq_len(K), each = K),
             r = as.vector(Ra), abs_r = abs(as.vector(Ra)), n_genes = length(genes)),
  data.frame(comparison = "w_S=1 vs alpha_F=0 benchmark", new_program = rep(seq_len(K), ncol(old)),
             other_program = rep(seq_len(ncol(old)), each = K),
             r = as.vector(Rb), abs_r = abs(as.vector(Rb)), n_genes = length(shared)))
write.csv(corr, file.path(OUT_DIR, "programs_correspondence.csv"), row.names = FALSE)
# Overlap of the top-25 genes at the matching pole (sign of r) for each best old match
topov <- do.call(rbind, lapply(1:K, function(k) {
  j <- which.max(abs(Rb[k, ])); s <- sign(Rb[k, j])
  wn <- EF[shared, k] * s; wo <- old[shared, j]
  tn <- names(sort(wn, decreasing = TRUE))[1:25]; to <- names(sort(wo, decreasing = TRUE))[1:25]
  data.frame(new_program = k, best_old_program = j, r = Rb[k, j], overlap_top25 = length(intersect(tn, to)),
             shared_top_genes = paste(intersect(tn, to), collapse = ", "))
}))
write.csv(topov, file.path(OUT_DIR, "programs_correspondence_topgenes.csv"), row.names = FALSE)
cat("Best old-program match per new program:\n"); print(topov)

# Figure F7: loadings heatmap ----
act_order <- c(which(prog$active), which(!prog$active))
col_lab <- sprintf("P%d%s", act_order, ifelse(prog$active[act_order], "*", ""))
rows <- do.call(rbind, lapply(act_order, function(k) {
  g <- c(head(pl[[k]]$pos, N_HEAT), head(pl[[k]]$neg, N_HEAT))
  data.frame(block = k, key = paste0(k, ":", g), gene = g)
}))
rows$block <- factor(rows$block, levels = act_order, labels = paste0("P", act_order))
long <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  data.frame(key = rows$key[i], block = rows$block[i], gene = rows$gene[i],
             program = factor(act_order, levels = act_order),
             w = vapply(act_order, function(k) EF[rows$gene[i], k] / max(abs(EF[, k])), 0))
}))
long$key <- factor(long$key, levels = rev(rows$key))
lab_fun <- setNames(rows$gene, rows$key)
p1 <- ggplot(long, aes(program, key, fill = w)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top", labels = setNames(col_lab, act_order)) +
  scale_y_discrete(labels = lab_fun) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-1, 1), name = "Loading / largest |loading|\nof the program") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8.5, colour = "black"),
        axis.text.x.top = element_text(size = 9, colour = "black"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 9, face = "bold"),
        panel.spacing.y = unit(4, "pt"), legend.position = "bottom",
        legend.title = element_text(size = 8.5), legend.text = element_text(size = 8.5),
        legend.key.width = unit(0.35, "in"), legend.key.height = unit(0.12, "in"))
ggsave(file.path(FIG_DIR, "sup_programs_heatmap.pdf"), p1, width = 6.5, height = 7.7,
       device = PDF_DEVICE)

# Figure F8: correspondence ----
#' Plot a matrix of loading correlations against another fit.
#' @param R Numeric correlation matrix, new programs in rows and comparison programs in columns.
#' @param xl,yl Axis titles.
#' @param lab_sign Whether printed labels retain the sign of the correlation.
#' @return A ggplot with absolute correlations mapped to tile fill.
#' @examples
#' tile(diag(2), "Comparison", "New", TRUE)
tile <- function(R, xl, yl, lab_sign) {
  d <- expand.grid(new = seq_len(nrow(R)), other = seq_len(ncol(R)))
  d$r <- as.vector(R)
  d$label <- ifelse(is.na(d$r), "NA", sprintf("%.2f", if (lab_sign) d$r else abs(d$r)))
  ggplot(d, aes(factor(other), factor(new, levels = nrow(R):1), fill = abs(r))) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label, colour = !is.na(r) & abs(r) > 0.6), size = 2.9) +
    scale_colour_manual(values = c("black", "white"), guide = "none") +
    scale_fill_gradient(low = "white", high = "#104281", limits = c(0, 1),
                        na.value = "#e9ecef", name = "|r|") +
    scale_x_discrete(position = "top") +
    labs(x = xl, y = yl) + theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 8.5, colour = "black"),
          axis.title = element_text(size = 9), legend.position = "none")
}
pa <- tile(Ra, "Two-step program", "Supervised program", FALSE)
pb <- tile(Rb, "Historical benchmark program", NULL, TRUE)
ggsave(file.path(FIG_DIR, "sup_programs_correspondence.pdf"),
       (pa | pb) + plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = "A"),
       width = 6.5, height = 3.2, device = PDF_DEVICE)

# Table T4 ----
act_txt <- ifelse(prog$active, "yes", "no")
rowsT <- vapply(act_order, function(k) {
  p <- prog[prog$program == k, ]
  sprintf("%d & %s & %s (%s) & %s & %s & %s & %s \\\\",
          k, tex_esc(p$short),
          fmt_signed(p$beta), formatC(p$beta_sd, format = "f", digits = 2), fmt_signed(p$z, 1),
          act_txt[k],
          tex_esc(paste(pl[[k]]$pos[1:min(N_POLE, 4)], collapse = ", ")),
          tex_esc(paste(pl[[k]]$neg[1:min(N_POLE, 4)], collapse = ", ")))
}, "")
tex <- c("\\begin{table}[tbp]", "\\centering", "\\small", "\\setlength{\\tabcolsep}{2.5pt}",
  "\\caption{Programs of the survival-supervised fit on the training cohorts ($n=273$). $E[\\beta]$ is the posterior mean coefficient per SD of the program score $Z_k = Y\\,E[F_k]$ (posterior SD in parentheses), $z = E[\\beta]/\\mathrm{SD}$, and a program is survival-active when $|z|\\ge 1.96$. A positive $\\beta$ indicates an association between a higher score and a higher hazard. Top genes are those with the largest positive and most negative loadings. The positive pole is the direction in which the score $Z_k$ increases.}",
  "\\label{tab:programs}",
  "\\begin{tabular}{@{}r p{1.12in} r r c p{1.3in} p{1.3in}@{}}", "\\toprule",
  "\\# & Label & $E[\\beta]$ (SD) & $z$ & Active & Positive pole genes & Negative pole genes \\\\",
  "\\midrule", rowsT, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(TAB_DIR, "tab_programs.tex"))

# Table T5 ----
topE <- do.call(rbind, lapply(which(prog$active), function(k) {
  x <- enr[enr$program == k & is.finite(enr$padj), ]
  rbind(transform(head(x[x$NES > 0, ], 2), pole = "+"),
        transform(head(x[x$NES < 0, ], 2), pole = "-"))
}))
set_name <- gsub("_", " ", sub("^(HALLMARK|REACTOME|KEGG_MEDICUS_REFERENCE|KEGG)_", "", topE$set))
adj_p <- ifelse(topE$padj < 0.001, "$<$0.001", formatC(topE$padj, format = "f", digits = 3))
rowsE <- sprintf("%d & %s & %s & %s & %s & %s \\\\", topE$program, ifelse(topE$pole == "+", "$+$", "$-$"),
                 tex_esc(set_name), topE$collection, fmt_signed(topE$NES), adj_p)
texE <- c("\\begin{table}[tbp]", "\\centering", "\\small",
  "\\caption{Two gene sets with the smallest adjusted $p$-values at each pole of the survival-active programs. Enrichment was computed on signed loadings over the 2,218 modeled genes. The test included MSigDB Hallmark, Reactome, KEGG and published PDAC sets. Benjamini--Hochberg adjustment was applied within program. Positive normalized enrichment scores (NES) indicate the positive pole.}",
  "\\label{tab:enrichment}", "\\begin{tabular}{@{}r c p{2.6in} l r r@{}}", "\\toprule",
  "Program & Pole & Gene set & Collection & NES & adj.\\ $p$ \\\\", "\\midrule", rowsE,
  "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(texE, file.path(TAB_DIR, "tab_enrichment.tex"))
# Subtype and stroma concordance (B3) ----
e <- new.env()
load(CMB_RDATA, envir = e)
signature_groups <- setNames(e$subtypeGeneList, e$schemaList)
MT <- signature_groups$MT
MS <- signature_groups$MS
CO <- signature_groups$Collisson
BA <- signature_groups$Bailey
sigs <- list(
  "Moffitt basal-like"    = MT$geneSymbol[MT[["Basal-like"]]],
  "Moffitt classical"     = MT$geneSymbol[MT$Classical],
  "Moffitt activated stroma" = MS$geneSymbol[MS$Activated],
  "Moffitt normal stroma" = MS$geneSymbol[MS$Normal],
  "Collisson classical"   = CO$geneSymbol[CO$Classical],
  "Collisson exocrine-like" = CO$geneSymbol[CO$Exocrine],
  "Collisson quasi-mesenchymal" = CO$geneSymbol[CO$QM],
  "Bailey squamous"       = BA$geneSymbol[BA$Squamous & !BA$NotUnique],
  "Bailey immunogenic"    = BA$geneSymbol[BA$Immunogenic & !BA$NotUnique],
  "Bailey pancreatic progenitor" = BA$geneSymbol[BA$PancreaticProgenitor & !BA$NotUnique],
  "Bailey ADEX"           = BA$geneSymbol[BA$ADEX & !BA$NotUnique],
  "DeSurv classical tumor" = pdac$DeSurv_D1_ClassicalTumor,
  "DeSurv stromal-immune"  = pdac$DeSurv_D2_StromalImmune,
  "DeSurv basal-like tumor" = pdac$DeSurv_D3_BasalLikeTumor)
#' Correlate each fixed program score with an outcome-free subtype score.
#' @param scores Patient by program projection matrix.
#' @param signature Numeric vector with one score per row of scores.
#' @param cohort Series name.
#' @param name Signature name.
#' @param n_sig,n_present Signature gene count and count available for scoring.
#' @return Aggregate Spearman correlations and sample counts by program.
#' @examples
#' correlate_signature(matrix(1:8, ncol = 2), 4:1, "Example", "Signature", 2, 2)
correlate_signature <- function(scores, signature, cohort, name, n_sig, n_present) {
  use <- is.finite(signature)
  if (sum(use) >= 3 && sd(signature[use]) > 0) {
    rho <- apply(scores[use, , drop = FALSE], 2, function(z) {
      if (sd(z) == 0) NA_real_ else cor(z, signature[use], method = "spearman")
    })
  } else {
    rho <- rep(NA_real_, ncol(scores))
  }
  data.frame(cohort = cohort, signature = name, n_sig = n_sig, n_present = n_present,
             n = sum(use), program = seq_len(ncol(scores)), spearman_r = rho)
}

#' Compare projections with mean gene z-scores over the measured signature genes.
#' @param expression Patient by gene matrix, already standardized within series.
#' @param scores Patient by program projection matrix.
#' @param cohort Series name.
#' @return Aggregate correlations and coverage for all signatures.
#' @examples
#' # score_signatures(Y, Z, "TCGA_PAAD")
score_signatures <- function(expression, scores, cohort) {
  stopifnot(nrow(expression) == nrow(scores))
  do.call(rbind, lapply(names(sigs), function(name) {
    genes_present <- intersect(unique(sigs[[name]]), colnames(expression))
    signature <- if (length(genes_present) >= 2) {
      rowMeans(expression[, genes_present, drop = FALSE])
    } else {
      rep(NA_real_, nrow(expression))
    }
    correlate_signature(scores, signature, cohort, name,
                        length(unique(sigs[[name]])), length(genes_present))
  }))
}

#' Compare projections with the published PurIST basal probability by sample ID.
#' @param scores Patient by program projection matrix.
#' @param sample_id Patient sample identifiers in score matrix row order.
#' @param cohort Series name, used to locate the saved PurIST probabilities.
#' @return Aggregate Spearman correlations and matched sample counts.
#' @examples
#' # score_purist(Z, tcga_ids, "TCGA_PAAD")
score_purist <- function(scores, sample_id, cohort) {
  path <- file.path(PDAC_DATA_ROOT, "original", paste0(cohort, ".caf_subtype.rds"))
  subtype <- readRDS(path)$Subtype
  stopifnot(length(sample_id) == nrow(scores), !anyNA(sample_id), !anyDuplicated(sample_id),
            !anyDuplicated(subtype$sampID))
  # The expression loader adds _seq to distinguish array and RNA-seq IDs.
  if (cohort == "PACA_AU_seq") {
    stopifnot(all(endsWith(sample_id, "_seq")))
    sample_id <- sub("_seq$", "", sample_id)
  }
  probability <- subtype$PurIST.prob[match(sample_id, subtype$sampID)]
  stopifnot(all(is.finite(probability)))
  correlate_signature(scores, probability, cohort, "PurIST basal probability", NA_integer_, NA_integer_)
}

training <- lapply(cfg$pdac$training_cohorts, function(cohort) {
  raw <- load_pdac_raw(cohort, PDAC_DATA_ROOT)
  idx <- which(dat$cohort == cohort)
  x <- raw$Y[, match(dat$genes, raw$gene_names), drop = FALSE]
  prep <- preprocess_desurv_cohort(x, dat$genes, top_n = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[cohort]], cohort_name = cohort,
    rank_transform = FALSE, per_platform_standardize = TRUE)
  # Check the prepared expression and follow-up rows before matching sample IDs.
  stopifnot(identical(dat$time[idx], raw$time), identical(dat$status[idx], raw$status),
            isTRUE(all.equal(unname(prep$Y[, match(dat$genes, prep$gene_names)]),
                             unname(Yfull[idx, ]), tolerance = 1e-12)))
  rbind(score_signatures(Yfull[idx, ], Z[idx, ], cohort),
        score_purist(Z[idx, ], raw$sampID, cohort))
})
external <- lapply(cfg$pdac$external_cohorts, function(cohort) {
  cat("Subtype concordance:", cohort, "\n")
  raw <- load_pdac_raw(cohort, PDAC_DATA_ROOT)
  prep <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[cohort]], cohort_name = cohort,
    rank_transform = FALSE, per_platform_standardize = TRUE)
  colnames(prep$Y) <- prep$gene_names
  index <- match(genes, prep$gene_names)
  hit <- !is.na(index)
  aligned <- matrix(0, raw$n, length(genes))
  aligned[, hit] <- prep$Y[, index[hit], drop = FALSE]
  # Project with the frozen loadings and training normalization. Only expression
  # and published subtype probabilities enter the external correlations.
  projection <- sweep(aligned %*% EF, 2, fit$state$EF_norms, "/")
  rbind(score_signatures(prep$Y, projection, cohort),
        score_purist(projection, raw$sampID, cohort))
})
conc <- do.call(rbind, c(list(score_signatures(Yfull, Z, "Pooled training")), training, external))
write.csv(conc, file.path(OUT_DIR, "programs_subtype_concordance.csv"), row.names = FALSE)

cat("Done.\n")
