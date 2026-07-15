# ============================================================
# tests/test_pathway_enrichment.R
# Tests for code/pathway_enrichment.R
#
# Usage: sourced by tests/run_tests.R
# ============================================================

source("code/pathway_enrichment.R")

cat("\n========================================\n")
cat("  Tests: pathway_enrichment.R\n")
cat("========================================\n\n")

cat("=== T1: load_d4_weights() ===\n")

run_test("T1.1: EF has 2064 rows (genes) and 7 columns (programs)", {
  d4 <- load_d4_weights()
  assert_true(nrow(d4$EF) == 2064, msg = "EF should have 2064 rows")
  assert_true(ncol(d4$EF) == 7, msg = "EF should have 7 columns")
})

run_test("T1.2: gene_names length matches EF rows, no duplicates", {
  d4 <- load_d4_weights()
  assert_length(d4$gene_names, 2064)
  assert_true(sum(duplicated(d4$gene_names)) == 0, msg = "gene names must be unique (or de-dup logged)")
})

run_test("T1.3: EF rows are named with gene symbols", {
  d4 <- load_d4_weights()
  assert_equal(rownames(d4$EF), d4$gene_names)
})

run_test("T1.4: EBeta has length 7", {
  d4 <- load_d4_weights()
  assert_length(d4$EBeta, 7)
})

run_test("T1.5: program_labels guards against the stale (pre-06-16) labeling", {
  d4 <- load_d4_weights()
  assert_equal(d4$program_labels[["7"]], "Adverse")
  assert_equal(d4$program_labels[["3"]], "Protective")
})

run_test("T1.6: inactive programs are labeled Inactive", {
  d4 <- load_d4_weights()
  for (k in c("1", "2", "4", "5", "6")) {
    assert_equal(d4$program_labels[[k]], "Inactive", msg = paste("program", k, "should be Inactive"))
  }
})

cat("=== T2: run_fgsea_program() ===\n")

# Synthetic ranking: 100 genes, weights strictly decreasing gene1 (highest) -> gene100 (lowest).
.fgsea_test_weights <- local({
  w <- setNames(as.numeric(100:1), paste0("gene", 1:100))
  w
})
.fgsea_matched_set   <- paste0("gene", 1:10)     # concentrated at the top of the ranking
.fgsea_scrambled_set <- paste0("gene", 91:100)   # concentrated at the bottom of the ranking

run_test("T2.1: matched (top-ranked) gene set enriches with NES > 0", {
  res <- run_fgsea_program(
    .fgsea_test_weights,
    list(matched = .fgsea_matched_set, scrambled = .fgsea_scrambled_set),
    seed = 1, minSize = 5, maxSize = 50
  )
  nes_matched <- res$NES[res$set == "matched"]
  assert_true(length(nes_matched) == 1 && nes_matched > 0, msg = "matched set should have NES > 0")
})

run_test("T2.2: matched and scrambled sets give different NES (not a constant-output bug)", {
  res <- run_fgsea_program(
    .fgsea_test_weights,
    list(matched = .fgsea_matched_set, scrambled = .fgsea_scrambled_set),
    seed = 1, minSize = 5, maxSize = 50
  )
  nes_matched   <- res$NES[res$set == "matched"]
  nes_scrambled <- res$NES[res$set == "scrambled"]
  assert_true(abs(nes_matched - nes_scrambled) > 0.5,
              msg = "matched vs scrambled NES should differ substantially")
})

run_test("T2.3: matched set has smaller p-value than scrambled set", {
  res <- run_fgsea_program(
    .fgsea_test_weights,
    list(matched = .fgsea_matched_set, scrambled = .fgsea_scrambled_set),
    seed = 1, minSize = 5, maxSize = 50
  )
  p_matched   <- res$pval[res$set == "matched"]
  p_scrambled <- res$pval[res$set == "scrambled"]
  assert_true(p_matched < p_scrambled, msg = "matched set should be more significant than scrambled")
})

run_test("T2.4: output columns include set, size, NES, pval, padj, leading_edge", {
  res <- run_fgsea_program(
    .fgsea_test_weights,
    list(matched = .fgsea_matched_set, scrambled = .fgsea_scrambled_set),
    seed = 1, minSize = 5, maxSize = 50
  )
  expected_cols <- c("set", "size", "NES", "pval", "padj", "leading_edge")
  assert_true(all(expected_cols %in% names(res)), msg = "missing expected fgsea output columns")
})

run_test("T2.5: zero-overlap gene sets fail loud rather than silently returning nothing", {
  err <- tryCatch({
    run_fgsea_program(
      .fgsea_test_weights,
      list(no_overlap = paste0("notagene", 1:10)),
      seed = 1, minSize = 5, maxSize = 50
    )
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", msg = "run_fgsea_program should fail loud on zero-overlap gene sets")
})

cat("=== T3: run_ora_program() ===\n")

.ora_background     <- paste0("gene", 1:100)
.ora_top_genes       <- paste0("gene", 1:10)      # top-10 weighted genes
.ora_matched_set    <- paste0("gene", 1:15)       # 10/15 overlap with top_genes
.ora_scrambled_set  <- paste0("gene", 50:65)      # 0 overlap with top_genes

run_test("T3.1: matched set (high overlap) is significant; scrambled (zero overlap) is absent", {
  res <- run_ora_program(
    .ora_top_genes, .ora_background,
    list(matched = .ora_matched_set, scrambled = .ora_scrambled_set)
  )
  p_matched <- res$pval[res$set == "matched"]
  assert_true(length(p_matched) == 1 && p_matched < 0.05,
              msg = "matched set should appear and be significant")
  # clusterProfiler::enricher() correctly omits zero-overlap sets from its result
  # entirely (no evidence of enrichment, not a "high p-value" row) -- this is
  # expected ORA behavior, not a bug.
  assert_true(!("scrambled" %in% res$set), msg = "zero-overlap set should not appear in ORA output")
})

run_test("T3.2: output columns include set, size, pval, padj, leading_edge", {
  res <- run_ora_program(
    .ora_top_genes, .ora_background,
    list(matched = .ora_matched_set, scrambled = .ora_scrambled_set)
  )
  expected_cols <- c("set", "size", "pval", "padj", "leading_edge")
  assert_true(all(expected_cols %in% names(res)), msg = "missing expected ORA output columns")
})

run_test("T3.3: zero-overlap-with-everything fails loud", {
  err <- tryCatch({
    run_ora_program(
      paste0("notagene", 1:10), .ora_background,
      list(matched = .ora_matched_set)
    )
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", msg = "run_ora_program should fail loud when top_genes overlap nothing")
})

cat("=== T4: parse_desurv_si_text() ===\n")

# A miniature synthetic mimic of the SI appendix's pdftotext -layout output:
# two factors, 3-gene rows x 2 columns = 6 genes per factor (instead of the
# real 45 rows x 6 columns = 270), same structural quirks as the real PDF:
# a form-feed before "Table", an en-dash column-range header row, a lone
# page-number footer line, and a mixed-case gene symbol (like C19orf33).
.mini_si_text <- c(
  "\fTable S7: Top 6 genes for DeSurv factor D1 (Factor D1 (Classical tumor)), ranked by x",
  "",
  "                1–63         4–6",
  "                GENE1        c19orf99",
  "                GENE2        GENE5",
  "                GENE3        GENE6",
  "",
  "                                   12",
  "\fTable S8: Top 6 genes for DeSurv factor D2 (Factor D2 (stromal/immune)), ranked by x",
  "",
  "                1–63         4–6",
  "                GENEA        GENED",
  "                GENEB        GENEE",
  "                GENEC        GENEF",
  "",
  "                                   13",
  "Bailey, Peter, David K Chang, et al. 2016. citation text follows here."
)

run_test("T4.1: parse_desurv_si_text extracts exactly 6 genes for D1 and D2", {
  res <- parse_desurv_si_text(.mini_si_text,
                               factor_patterns = c(D1 = "Table S7:", D2 = "Table S8:"),
                               end_pattern = "^Bailey, Peter")
  assert_length(res$D1, 6)
  assert_length(res$D2, 6)
})

run_test("T4.2: parse_desurv_si_text drops the column-range header and page-number footer", {
  res <- parse_desurv_si_text(.mini_si_text,
                               factor_patterns = c(D1 = "Table S7:", D2 = "Table S8:"),
                               end_pattern = "^Bailey, Peter")
  assert_false("12" %in% res$D1, msg = "page-number footer leaked into D1 gene list")
  assert_false(any(grepl("–", res$D1)), msg = "column-range header leaked into D1 gene list")
})

run_test("T4.3: parse_desurv_si_text keeps mixed-case gene symbols (e.g. C19orf-style)", {
  res <- parse_desurv_si_text(.mini_si_text,
                               factor_patterns = c(D1 = "Table S7:", D2 = "Table S8:"),
                               end_pattern = "^Bailey, Peter")
  assert_true("c19orf99" %in% res$D1, msg = "mixed-case gene symbol was incorrectly dropped")
})

run_test("T4.4: parse_desurv_si_text fails loud if a factor doesn't hit the expected gene count", {
  err <- tryCatch({
    parse_desurv_si_text(.mini_si_text,
                          factor_patterns = c(D1 = "Table S7:", D2 = "Table S8:"),
                          end_pattern = "^Bailey, Peter",
                          expected_n = 270)
    "no error"
  }, error = function(e) "error")
  assert_equal(err, "error", msg = "should fail loud when gene count != expected_n")
})

cat("=== T5: get_msigdb_collections() ===\n")

run_test("T5.1: Hallmark collection returns ~50 non-empty gene sets", {
  res <- get_msigdb_collections(collections = list(Hallmark = list(collection = "H")))
  assert_true("Hallmark" %in% names(res))
  assert_true(length(res$Hallmark) >= 45 && length(res$Hallmark) <= 55,
              msg = sprintf("expected ~50 Hallmark sets, got %d", length(res$Hallmark)))
  assert_true(all(vapply(res$Hallmark, length, integer(1)) > 0), msg = "found an empty Hallmark set")
})

run_test("T5.2: gene sets are character vectors of gene symbols", {
  res <- get_msigdb_collections(collections = list(Hallmark = list(collection = "H")))
  one_set <- res$Hallmark[[1]]
  assert_true(is.character(one_set), msg = "gene set should be a character vector")
  assert_true(all(grepl("^[A-Za-z0-9.-]+$", one_set)), msg = "gene set contains non-symbol tokens")
})

cat("=== T6: top_n_genes_table() ===\n")

run_test("T6.1: returns top-N rows per program, ordered by descending weight", {
  d4 <- load_d4_weights()
  tbl <- top_n_genes_table(d4$EF, programs = c(3, 7), program_labels = d4$program_labels, n = 10)
  assert_true(nrow(tbl) == 20, msg = "expected 10 rows x 2 programs = 20")
  t3 <- tbl[tbl$program == 3, ]
  assert_true(all(diff(t3$weight) <= 0), msg = "weights should be sorted descending within a program")
})

run_test("T6.2: program 7's top gene matches the known top-weighted gene (ITGA3)", {
  d4 <- load_d4_weights()
  tbl <- top_n_genes_table(d4$EF, programs = 7, program_labels = d4$program_labels, n = 5)
  assert_equal(tbl$gene[1], "ITGA3")
  assert_equal(tbl$label[1], "Adverse")
})

report_results("pathway_enrichment.R")
