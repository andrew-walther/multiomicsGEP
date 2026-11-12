# ============================================================
# Script: test_methylation_factor_characterization_k5.R
# Purpose: Check K = 5 characterization classification, annotation and permutation.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R
# ============================================================
source("results/multimodal_real/characterize_methylation_factors_k5.R")
set.seed(17)
shared <- rnorm(144)
independent <- qr.resid(qr(cbind(1, shared)), rnorm(144))
methylation <- cbind(shared, independent)
expression <- cbind(shared + rnorm(144, sd = 0.05), rnorm(144))
r <- abs(safe_cor(methylation, expression))
best <- apply(r, 1, max)
stopifnot(identical(unname(classify_specific(best, c(4L, 3L), c(0.9, 0.8))), c(FALSE, TRUE)),
          identical(classify_specific(0.1, 2L, 0.9), FALSE),
          identical(classify_specific(0.1, 3L, 0.69), FALSE),
          identical(classify_specific(0.5, 3L, 0.7), FALSE),
          identical(classify_specific(0.49, 3L, 0.7), TRUE),
          is.na(classify_specific(NA_real_, 3L, 0.8)))
a <- data.frame(UCSC_RefGene_Name = c("A;B", "C"), UCSC_RefGene_Group = c("TSS200;Body", "5'UTR"),
                Relation_to_Island = c("N_Shore", ""), row.names = c("cg1", "cg2"))
b <- expand_annotation(a)
stopifnot(nrow(b) == 3L, b$gene[b$cpg == "cg1"] == c("A", "B"),
          b$region[b$cpg == "cg1"] == c("TSS200", "Body"),
          identical(island_class(c("Island", "N_Shore", "S_Shelf", "")), c("island", "shore", "shelf", "open sea")))
p1 <- overlap_permutation(letters, letters[1:10], letters[1:5], 42L, 100L)
p2 <- overlap_permutation(letters, letters[1:10], letters[1:5], 42L, 100L)
p3 <- overlap_permutation(letters, letters[1:10], letters[20:24], 42L, 100L)
stopifnot(identical(p1, p2), p1["observed"] == 5, p3["observed"] == 0,
          p1["p_upper"] < p3["p_upper"], p1["null_sd"] > 0)
cat("K = 5 methylation characterization: 15 assertions passed.\n")
