# ============================================================
# Script:  results/benchmark_sim/run_ebmf_pathway_concordance.R
# Purpose: ROADMAP.md "A/B comparison: SSBMF vs unsupervised EBMF" (line 359),
#          Part 3 -- does the unsupervised EBMF factor that best matches each
#          survival-active program (run_ebmf_factor_correspondence.R: EBMF_F1
#          <-> Program 3, EBMF_F2 <-> Program 7) show the SAME pathway/gene-set
#          enrichment as that program's own fgsea results, or does supervision
#          change not just which genes are selected but what biology they
#          represent?
#
#          EBMF's factors are signed (flashier's default point_normal prior on
#          F), unlike SBMF's point-exponential F (non-negative by construction,
#          hence the one-sided scoreType="pos" used for Programs 3/7). Each
#          EBMF factor is sign-aligned to its matched program (multiplied by
#          the sign of their correlation, from run_ebmf_factor_correspondence.R)
#          before a two-sided fgsea (scoreType="std") -- an unsupervised factor
#          has no inherent "positive" direction, so alignment must be done
#          explicitly rather than assumed.
#
#          Enrichment is run against the SAME collections already used for
#          Programs 3/7 (MSigDB Hallmark/Reactome/KEGG/GO:BP + the existing
#          custom PDAC collection), then concordance is read off directly by
#          checking whether each of Program 3/7's own significant sets
#          (padj<0.10, from the existing fgsea_results_all.rds) also comes up
#          significant for the matched, sign-aligned EBMF factor.
#
#   Output: results/benchmark_sim/outputs/pathway_enrichment/
#     T6_ebmf_pathway_concordance.csv   (every Program-3/7 significant set,
#                                         with the matched EBMF factor's NES/padj)
#     ebmf_fgsea_results.rds            (full fgsea output for EBMF_F1/F2, all collections)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-03
# Usage:   Rscript results/benchmark_sim/run_ebmf_pathway_concordance.R
#          (requires prior runs of run_ebmf_cox_external.R, run_pathway_enrichment.R,
#           and run_ebmf_factor_correspondence.R)
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

source("code/pathway_enrichment.R")

OUT_DIR <- "results/benchmark_sim/outputs/pathway_enrichment"
ACTIVE_PROGRAMS <- c(3, 7)

ebmf_path <- "results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit.rds"
correspondence_path <- file.path(OUT_DIR, "T5_ebmf_factor_correspondence.csv")
pdac_genesets_path <- file.path(OUT_DIR, "pdac_genesets.rds")
fgsea_all_path <- file.path(OUT_DIR, "fgsea_results_all.rds")

for (p in c(ebmf_path, correspondence_path, pdac_genesets_path, fgsea_all_path)) {
  if (!file.exists(p)) {
    stop("run_ebmf_pathway_concordance: ", p, " not found -- run the prerequisite ",
         "scripts listed in this script's header first")
  }
}

ebmf <- readRDS(ebmf_path)
d4 <- load_d4_weights()
correspondence <- read.csv(correspondence_path, stringsAsFactors = FALSE)
pdac_genesets <- readRDS(pdac_genesets_path)
fgsea_all <- readRDS(fgsea_all_path)

if (!identical(ebmf$train_genes, d4$gene_names)) {
  stop("run_ebmf_pathway_concordance: EBMF's train_genes and D4's gene_names ",
       "differ -- re-derive the gene alignment before proceeding (same check as ",
       "run_ebmf_factor_correspondence.R).")
}

# Section: sign-align each matched EBMF factor to its program ----

aligned_factors <- lapply(seq_len(nrow(correspondence)), function(i) {
  row <- correspondence[i, ]
  j <- row$best_ebmf_factor
  w <- ebmf$F_ebmf[, j] * sign(row$r)
  names(w) <- d4$gene_names
  w
})
names(aligned_factors) <- sprintf("EBMF_F%d_aligned_to_P%d", correspondence$best_ebmf_factor, correspondence$program)

cat("Sign-aligned EBMF factors (multiplied by sign of correlation with matched program):\n")
print(data.frame(factor = names(aligned_factors), r = correspondence$r,
                  sign_applied = sign(correspondence$r)))

# Section: gene-set collections -- same ones used for Programs 3/7 ----

msigdb_collections <- get_msigdb_collections()
collections <- c(msigdb_collections, list(PDAC_custom = pdac_genesets))
cat(sprintf("Collections: %s\n", paste(names(collections), collapse = ", ")))

# Section: fgsea on the two sign-aligned EBMF factors, two-sided ----

ebmf_fgsea <- do.call(rbind, lapply(seq_along(aligned_factors), function(i) {
  fname <- names(aligned_factors)[i]
  program_k <- correspondence$program[i]
  per_collection <- lapply(names(collections), function(coll_name) {
    tryCatch({
      res <- run_fgsea_program(aligned_factors[[i]], collections[[coll_name]], seed = 1,
                                collection = coll_name, scoreType = "std")
      res$ebmf_factor <- fname
      res$matched_program <- program_k
      res
    }, error = function(e) {
      message(sprintf("run_ebmf_pathway_concordance: fgsea skipped for %s, collection %s (%s)",
                       fname, coll_name, conditionMessage(e)))
      NULL
    })
  })
  do.call(rbind, per_collection)
}))
saveRDS(ebmf_fgsea, file.path(OUT_DIR, "ebmf_fgsea_results.rds"))

# Section: T6 -- concordance table, Program 3/7's own significant sets vs. the matched EBMF factor ----

t6 <- do.call(rbind, lapply(ACTIVE_PROGRAMS, function(k) {
  prog_sig <- fgsea_all[fgsea_all$program == k & fgsea_all$padj < 0.10, ]
  ebmf_fname <- names(aligned_factors)[correspondence$program == k]
  ebmf_res <- ebmf_fgsea[ebmf_fgsea$ebmf_factor == ebmf_fname, ]

  do.call(rbind, lapply(seq_len(nrow(prog_sig)), function(i) {
    set_name <- prog_sig$set[i]
    coll_name <- prog_sig$collection[i]
    match_row <- ebmf_res[ebmf_res$set == set_name & ebmf_res$collection == coll_name, ]
    data.frame(
      program = k, label = d4$program_labels[[as.character(k)]],
      collection = coll_name, set = set_name,
      program_NES = prog_sig$NES[i], program_padj = prog_sig$padj[i],
      ebmf_factor = ebmf_fname,
      ebmf_NES = if (nrow(match_row) == 1) match_row$NES else NA,
      ebmf_padj = if (nrow(match_row) == 1) match_row$padj else NA,
      concordant = if (nrow(match_row) == 1) !is.na(match_row$padj) && match_row$padj < 0.10 &&
        sign(match_row$NES) == sign(prog_sig$NES[i]) else FALSE,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(t6, file.path(OUT_DIR, "T6_ebmf_pathway_concordance.csv"), row.names = FALSE)

cat(sprintf("\nT6: %d/%d of Programs 3+7's significant sets are also significant (padj<0.10, same direction) ",
            sum(t6$concordant), nrow(t6)))
cat("for their matched, sign-aligned unsupervised EBMF factor:\n")
print(t6[, c("program", "label", "collection", "set", "program_padj", "ebmf_padj", "concordant")])

cat("\nrun_ebmf_pathway_concordance.R complete. Outputs written to", OUT_DIR, "\n")
