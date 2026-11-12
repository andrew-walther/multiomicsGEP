# ============================================================
# Script: k5_program_correspondence.R
# Purpose: Match the programs of the K = 5 working fit to those of the reported K = 7 fit
#          by training score correlation and top-gene overlap.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: base R; code/predict_ssbmf_yfb.R
# Run (repo root): Rscript results/benchmark_sim/k5_program_correspondence.R
# ============================================================
# Both fits are the S-11 all-patient w_S = 1 fits on the same 2,218 genes and 273 training
# patients (per-cohort z-scores), differing only in rank. The program score is the
# training-normalized projection z_k = Y E[f_k] / ||Y E[f_k]||-scale (unit SD), the quantity
# multiplied by beta_k in the YFB linear predictor eta = sum_k z_k beta_k.
# Outputs are aggregate: all 5 x 7 score correlations, and for each K = 5 program its best
# K = 7 match with the overlap of the 25 largest loadings at the matched pole. The training
# matrix is read from the git-ignored data/ folder; no patient-level value is written.

suppressPackageStartupMessages(source("code/predict_ssbmf_yfb.R"))
K7_FIT <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits/size_D4_N2218_w1_f0.rds"
K5_FIT <- "results/benchmark_sim/outputs/k5_working_2026-11-12/s11/fits/size_D4_N2218_w1_f0.rds"
OUT    <- "results/benchmark_sim/outputs/k5_working_2026-11-12/program_correspondence"
N_TOP  <- 25L     # genes per pole, as in the Chapter 4 historical correspondence table

#' Posterior coefficient summary and survival activity of one fit
#' @param s Checkpoint state (EBeta, EBeta2, EF). @return Data frame: program, beta, sd, z, active.
#' @examples
#' # beta_summary(fit$state)
beta_summary <- function(s) {
  sdb <- sqrt(pmax(s$EBeta2 - s$EBeta^2, 0))
  alive <- colSums(s$EF^2) > 1e-8
  data.frame(program = seq_along(s$EBeta), beta = s$EBeta, beta_sd = sdb,
             z = ifelse(sdb > 0, s$EBeta / sdb, NA_real_),
             active = alive & sdb > 0 & abs(s$EBeta) >= 1.96 * sdb)
}

#' Genes with the N largest loadings at one pole
#' @param w Named loading vector. @param sgn +1 for the positive pole, -1 for the negative.
#' @param n Number of genes. @return Character vector.
#' @examples
#' # top_pole(c(a = 1, b = -2), -1, 1)  # "b"
top_pole <- function(w, sgn, n) names(sort(sgn * w, decreasing = TRUE))[seq_len(n)]

if (sys.nframe() == 0L) {
  z  <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  f7 <- readRDS(K7_FIT); f5 <- readRDS(K5_FIT)
  stopifnot(identical(f7$genes, f5$genes), length(f5$genes) == 2218L,
            ncol(f7$state$EF) == 7L, ncol(f5$state$EF) == 5L)
  Y  <- z$Y[, match(f5$genes, z$genes), drop = FALSE]
  Z7 <- predict_ssbmf_yfb(Y, f7$state$EF, f7$state$EBeta, f7$state$EF_norms)$ZF_test
  Z5 <- predict_ssbmf_yfb(Y, f5$state$EF, f5$state$EBeta, f5$state$EF_norms)$ZF_test
  R  <- cor(Z5, Z7)                          # 5 x 7 training score correlations
  b5 <- beta_summary(f5$state); b7 <- beta_summary(f7$state)
  EF5 <- f5$state$EF; EF7 <- f7$state$EF; rownames(EF5) <- rownames(EF7) <- f5$genes

  all_pairs <- data.frame(k5_program = rep(1:5, 7), k7_program = rep(1:7, each = 5),
                          score_r = as.vector(R))
  best <- do.call(rbind, lapply(1:5, function(k) {
    j <- which.max(abs(R[k, ])); s <- sign(R[k, j])
    pos <- length(intersect(top_pole(EF5[, k], 1, N_TOP), top_pole(EF7[, j], s, N_TOP)))
    neg <- length(intersect(top_pole(EF5[, k], -1, N_TOP), top_pole(EF7[, j], -s, N_TOP)))
    second <- sort(abs(R[k, -j]), decreasing = TRUE)[1]
    data.frame(k5_program = k, k7_match = j, score_r = R[k, j], next_best_abs_r = second,
               loading_r = cor(EF5[, k], EF7[, j]),
               top25_overlap_pos = pos, top25_overlap_neg = neg,
               k5_beta = b5$beta[k], k5_z = b5$z[k], k5_active = b5$active[k],
               k7_beta = b7$beta[j], k7_z = b7$z[j], k7_active = b7$active[j])
  }))
  unmatched <- setdiff(1:7, best$k7_match)
  stopifnot(!anyDuplicated(best$k7_match))
  k7_unmatched <- data.frame(k7_program = unmatched,
    best_abs_r_with_k5 = apply(abs(R[, unmatched, drop = FALSE]), 2, max),
    k7_beta = b7$beta[unmatched], k7_z = b7$z[unmatched], k7_active = b7$active[unmatched])
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_pairs, file.path(OUT, "score_correlations.csv"), row.names = FALSE)
  write.csv(best, file.path(OUT, "matches.csv"), row.names = FALSE)
  write.csv(k7_unmatched, file.path(OUT, "k7_unmatched.csv"), row.names = FALSE)
  print(round(R, 3)); print(best, digits = 3); print(k7_unmatched, digits = 3)
}
