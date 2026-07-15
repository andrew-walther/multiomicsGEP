# ============================================================
# Script:  results/benchmark_sim/run_pathway_enrichment.R
# Purpose: Pathway/gene-set enrichment on the D4 (recommended) model's 7 gene
#          expression programs -- fgsea ranked-by-weight (primary) + ORA on
#          top-N genes (confirmatory) across MSigDB collections (Hallmark,
#          Reactome, KEGG, GO:BP) and a custom PDAC gene-set collection
#          (Moffitt basal/classical, Bailey 2016, DeSurv D1-D3).
#
#          Enriches all 7 programs (C1: active-vs-inactive sanity check) but
#          focuses reported output on Programs 3 (Protective) and 7 (Adverse).
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     T1_enrichment_active.csv   -- fgsea results, Programs 3 & 7, all collections
#     T2_top_genes.csv           -- top-100 weighted genes, Programs 3 & 7
#     C1_all_programs_summary.csv -- coherent-hit count per program (all 7)
#     fgsea_results_all.rds     -- full fgsea result data.frame (all 7 programs)
#     ora_results_active.rds   -- full ORA result data.frame (Programs 3 & 7, N in {50,100,150})
#     pdac_genesets.rds / genesets_manifest.txt -- from build_pdac_genesets()
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-14
# Usage:   Rscript results/benchmark_sim/run_pathway_enrichment.R
#          Rscript results/benchmark_sim/run_pathway_enrichment.R --quick
#          (--quick restricts to the Hallmark collection only, for a fast dry run)
# ============================================================

args       <- commandArgs(trailingOnly = TRUE)
QUICK_MODE <- "--quick" %in% args

# Navigate to project root if invoked from a subdirectory
if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Section: Load weights & gene-set collections ----

d4 <- load_d4_weights()
ALL_PROGRAMS <- 1:7
ACTIVE_PROGRAMS <- c(3, 7)

msigdb_collections <- if (QUICK_MODE) {
  get_msigdb_collections(collections = list(Hallmark = list(collection = "H")))
} else {
  get_msigdb_collections()
}

pdac_root <- Sys.getenv("PDAC_DATA_ROOT")
desurv_si_pdf <- if (nchar(pdac_root) > 0) {
  file.path(dirname(pdac_root), "papers", "DeSurv", "si_appendix.pdf")
} else {
  ""
}
pdac_custom_available <- nchar(pdac_root) > 0 && dir.exists(pdac_root) &&
  file.exists(desurv_si_pdf) && nchar(Sys.which("pdftotext")) > 0

if (pdac_custom_available) {
  pdac_genesets <- build_pdac_genesets(pdac_root, desurv_si_pdf, OUT_DIR)
  collections <- c(msigdb_collections, list(PDAC_custom = pdac_genesets))
} else {
  message("run_pathway_enrichment: PDAC_DATA_ROOT and/or DeSurv SI appendix PDF ",
          "unavailable -- skipping the custom PDAC gene-set collection (Moffitt/Bailey/DeSurv). ",
          "MSigDB collections only.")
  collections <- msigdb_collections
}

cat(sprintf("Collections: %s\n", paste(names(collections), collapse = ", ")))
cat(sprintf("Programs enriched: %s\n", paste(ALL_PROGRAMS, collapse = ", ")))

# Section: fgsea -- all 7 programs, all collections ----

fgsea_results <- do.call(rbind, lapply(ALL_PROGRAMS, function(k) {
  weights_k <- d4$EF[, k]
  per_collection <- lapply(names(collections), function(coll_name) {
    tryCatch({
      res <- run_fgsea_program(weights_k, collections[[coll_name]], seed = 1,
                                collection = coll_name)
      res$program <- k
      res$label <- d4$program_labels[[as.character(k)]]
      res
    }, error = function(e) {
      message(sprintf("run_pathway_enrichment: fgsea skipped for program %d, collection %s (%s)",
                       k, coll_name, conditionMessage(e)))
      NULL
    })
  })
  do.call(rbind, per_collection)
}))

saveRDS(fgsea_results, file.path(OUT_DIR, "fgsea_results_all.rds"))

t1 <- fgsea_results[fgsea_results$program %in% ACTIVE_PROGRAMS & fgsea_results$padj < 0.10, ]
t1 <- t1[order(t1$program, t1$padj),
         c("program", "label", "collection", "set", "size", "NES", "pval", "padj", "leading_edge")]
write.csv(t1, file.path(OUT_DIR, "T1_enrichment_active.csv"), row.names = FALSE)

n_hits_active <- sapply(ACTIVE_PROGRAMS, function(k) sum(t1$program == k))
if (any(n_hits_active == 0)) {
  warning(sprintf("run_pathway_enrichment: program(s) %s returned ZERO sets at padj<0.10",
                   paste(ACTIVE_PROGRAMS[n_hits_active == 0], collapse = ", ")))
}
cat(sprintf("T1: %d enriched sets at padj<0.10 across Programs %s\n",
            nrow(t1), paste(ACTIVE_PROGRAMS, collapse = ",")))

# Section: C1 -- active vs inactive program summary ----

c1 <- do.call(rbind, lapply(ALL_PROGRAMS, function(k) {
  sub <- fgsea_results[fgsea_results$program == k & fgsea_results$padj < 0.10, ]
  data.frame(program = k, label = d4$program_labels[[as.character(k)]],
             n_coherent_hits = nrow(sub))
}))
write.csv(c1, file.path(OUT_DIR, "C1_all_programs_summary.csv"), row.names = FALSE)
cat("C1 (coherent-hit count per program):\n")
print(c1)

# Section: ORA -- Programs 3 & 7, top-N in {50, 100, 150} ----

background <- d4$gene_names
ora_results <- do.call(rbind, lapply(ACTIVE_PROGRAMS, function(k) {
  do.call(rbind, lapply(c(50, 100, 150), function(n) {
    top_genes <- top_n_genes_table(d4$EF, k, d4$program_labels, n = n)$gene
    do.call(rbind, lapply(names(collections), function(coll_name) {
      tryCatch({
        res <- run_ora_program(top_genes, background, collections[[coll_name]],
                                collection = coll_name)
        res$program <- k
        res$label <- d4$program_labels[[as.character(k)]]
        res$top_n <- n
        res
      }, error = function(e) {
        message(sprintf("run_pathway_enrichment: ORA skipped for program %d, N=%d, collection %s (%s)",
                         k, n, coll_name, conditionMessage(e)))
        NULL
      })
    }))
  }))
}))
saveRDS(ora_results, file.path(OUT_DIR, "ora_results_active.rds"))
cat(sprintf("ORA: %d total enriched-set rows across Programs %s x N in {50,100,150}\n",
            nrow(ora_results), paste(ACTIVE_PROGRAMS, collapse = ",")))

# Section: T2 -- top-100 weighted genes, Programs 3 & 7 ----

t2 <- top_n_genes_table(d4$EF, ACTIVE_PROGRAMS, d4$program_labels, n = 100)
write.csv(t2, file.path(OUT_DIR, "T2_top_genes.csv"), row.names = FALSE)

# Section: F1-F3 -- figures ----

suppressPackageStartupMessages(library(ggplot2))

f1_data <- prepare_dotplot_data(fgsea_results, programs = ACTIVE_PROGRAMS, top_n = 10)
if (nrow(f1_data) > 0) {
  f1 <- plot_enrichment_dotplot(f1_data)
  ggsave(file.path(OUT_DIR, "F1_enrichment_dotplot.png"), f1, width = 16, height = 8, dpi = 150)
  ggsave(file.path(OUT_DIR, "F1_enrichment_dotplot.pdf"), f1, width = 16, height = 8)
  cat("F1 written.\n")
} else {
  message("run_pathway_enrichment: no sets available for F1 (t1 was empty) -- skipping F1")
}

# F2: running-ES plot for the single top (lowest padj) set per active program.
for (k in ACTIVE_PROGRAMS) {
  top_set_row <- t1[t1$program == k, ][1, ]
  if (nrow(top_set_row) == 0 || is.na(top_set_row$set)) {
    message(sprintf("run_pathway_enrichment: no headline set for program %d -- skipping F2", k))
    next
  }
  geneset <- collections[[top_set_row$collection]][[top_set_row$set]]
  f2 <- plot_running_es(d4$EF[, k], geneset,
                         title = sprintf("%s (Program %d, %s)", top_set_row$set, k, top_set_row$label))
  fname <- gsub("[^A-Za-z0-9_-]", "_", sprintf("F2_running_es_program%d_%s", k, top_set_row$set))
  ggsave(file.path(OUT_DIR, paste0(fname, ".png")), f2, width = 7, height = 4.5, dpi = 120)
  cat(sprintf("F2 written for program %d (%s).\n", k, top_set_row$set))
}

# F3: gene-weight heatmap over the union of headline leading-edge genes (Programs 3 & 7).
f3_genes <- unique(unlist(lapply(strsplit(t1$leading_edge[t1$program %in% ACTIVE_PROGRAMS], ";"),
                                  function(g) g[seq_len(min(20, length(g)))])))
if (length(f3_genes) > 0) {
  plot_geneweight_heatmap(d4$EF, f3_genes, d4$program_labels,
                           filename = file.path(OUT_DIR, "F3_geneweight_heatmap.png"))
  cat(sprintf("F3 written (%d genes).\n", length(f3_genes)))
} else {
  message("run_pathway_enrichment: no leading-edge genes available for F3 -- skipping F3")
}

cat("\nrun_pathway_enrichment.R complete. Outputs written to", OUT_DIR, "\n")
