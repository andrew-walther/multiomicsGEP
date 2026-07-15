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

report_results("pathway_enrichment.R")
