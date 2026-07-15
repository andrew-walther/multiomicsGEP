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

report_results("pathway_enrichment.R")
