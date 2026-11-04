# ============================================================
# Script: summarize_yfb_supervised_real.R
# Purpose: Attempt S-11 tables: S3 programs and CV C, feature-set-size curve with the
#          registered preferred-size rule and elbow plot, survival-screening comparison.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: survival, yaml; code/feature_sensitivity.R, code/concordance_ci.R,
#               code/predict_ssbmf_yfb.R
# Run: Rscript results/benchmark_sim/summarize_yfb_supervised_real.R
# ============================================================
# Reads checkpoints from run_yfb_supervised_real.R; refits nothing. Needs the git-ignored
# prepared matrix (data/) only for program agreement and pooled out-of-fold C. Writes
# aggregate tables and gene-level program summaries; no patient-level values.

suppressPackageStartupMessages(library(survival))
source("code/concordance_ci.R"); source("code/predict_ssbmf_yfb.R"); source("code/feature_sensitivity.R")
OUT   <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04"
LOCAL <- "data/yfb_supervised_real_2026-11-04.rds"
z <- readRDS(LOCAL)
ck <- lapply(list.files(file.path(OUT, "fits"), full.names = TRUE), readRDS)
m <- do.call(rbind, lapply(ck, function(x) { r <- x$metrics
  for (nm in c("converged", "n_iter", "n_alive", "n_active", "max_abs_beta", "share_active_median",
               "share_alive_median", "recon_train", "fit_sec", "C", "cvpl", "risk_sd"))
    if (is.null(r[[nm]])) r[[nm]] <- NA
  r }))
key <- function(part, set, N, arm, fold) sprintf("%s_%s_N%d_%s_f%d", part, set, N, arm, fold)
names(ck) <- key(m$part, m$set, m$N, m$arm, m$fold)
write.csv(m, file.path(OUT, "fits.csv"), row.names = FALSE)
cat("Fits read:", nrow(m), " errors:", sum(!is.na(m$error)), "\n")

# Per (part, set, N, arm) summary ----

#' Cross-validated summary of one group of fits
#' @param d Rows of one (part, set, N, arm) group, folds 0..5.
#' @return One-row data frame.
#' @examples
#' # cv_row(d)
cv_row <- function(d) {
  cvd <- d[d$fold > 0, ]; full <- d[d$fold == 0, ]
  oof <- rep(NA_real_, length(z$time))
  for (i in seq_len(nrow(cvd))) {
    x <- ck[[key(cvd$part[i], cvd$set[i], cvd$N[i], cvd$arm[i], cvd$fold[i])]]
    if (!is.null(x$risk)) oof[x$te] <- x$risk
  }
  data.frame(part = d$part[1], set = d$set[1], N = d$N[1], arm = d$arm[1],
             n_genes = median(d$n_genes), folds = nrow(cvd), errors = sum(!is.na(d$error)),
             mean_C = mean(cvd$C, na.rm = TRUE), se_C = sd(cvd$C, na.rm = TRUE) / sqrt(sum(!is.na(cvd$C))),
             n_C_na = sum(is.na(cvd$C)),
             pooled_oof_C = if (anyNA(oof)) NA_real_ else frozen_reverse_cindex(oof, z$time, z$status),
             cv_partial_loglik = sum(cvd$cvpl), folds_converged = sum(cvd$converged, na.rm = TRUE),
             mean_active_cv = mean(cvd$n_active), full_converged = if (nrow(full)) full$converged else NA,
             full_active = if (nrow(full)) full$n_active else NA,
             full_alive = if (nrow(full)) full$n_alive else NA,
             full_max_beta = if (nrow(full)) full$max_abs_beta else NA,
             share_alive_median = median(d$share_alive_median, na.rm = TRUE),
             fit_sec_median = median(d$fit_sec, na.rm = TRUE))
}
grp <- split(m, list(m$part, m$set, m$N, m$arm), drop = TRUE)
summ <- do.call(rbind, lapply(grp, cv_row))
summ <- summ[order(summ$part, summ$arm, summ$n_genes, summ$set), ]
write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)

# Feature-set size ----
sz <- summ[summ$part == "size", ]
if (nrow(sz)) {
  ## Program agreement with the all-genes fit (same arm), patient-level projections ----
  big <- names(z$sets)[which.max(lengths(z$sets))]
  agree <- do.call(rbind, lapply(unique(sz$arm), function(a) {
    ref <- ck[[key("size", big, length(z$sets[[big]]), a, 0L)]]
    if (is.null(ref$state)) return(NULL)
    refY <- z$Y[, match(ref$genes, z$genes)]
    rs <- ref$state; sdb <- sqrt(pmax(rs$EBeta2 - rs$EBeta^2, 0))
    ref_active <- which(abs(rs$EBeta) >= 1.96 * sdb & sdb > 0 & colSums(rs$EF^2) > 1e-8)
    do.call(rbind, lapply(names(z$sets), function(s) {
      x <- ck[[key("size", s, length(z$sets[[s]]), a, 0L)]]
      if (is.null(x$state)) return(NULL)
      newY <- z$Y[, match(x$genes, z$genes)]
      ag <- feature_program_agreement(rs, ref$genes, x$state, x$genes, refY, newY)
      pc <- abs(ag$projection_correlation)
      best <- if (length(ref_active)) apply(pc[ref_active, , drop = FALSE], 1, max, na.rm = TRUE) else numeric()
      rr <- predict_ssbmf_yfb(refY, rs$EF, rs$EBeta, EF_norms = rs$EF_norms)$risk_scores
      rn <- predict_ssbmf_yfb(newY, x$state$EF, x$state$EBeta, EF_norms = x$state$EF_norms)$risk_scores
      data.frame(arm = a, set = s, n_genes = length(x$genes), ref_active = length(ref_active),
                 min_match_active = if (length(best)) min(best) else NA_real_,
                 matches = paste(sprintf("%.2f", best), collapse = " "),
                 risk_cor = if (sd(rn) > 0 && sd(rr) > 0) cor(rr, rn) else NA_real_)
    }))
  }))
  write.csv(agree, file.path(OUT, "size_agreement.csv"), row.names = FALSE)

  ## Registered preferred-size rule ----
  pref <- do.call(rbind, lapply(unique(sz$arm), function(a) {
    s <- sz[sz$arm == a, ]; s <- s[order(s$n_genes), ]
    ag <- agree[agree$arm == a, ]
    ok_prog <- function(set) { r <- ag[ag$set == set, ]
      nrow(r) == 1 && (r$ref_active == 0 || isTRUE(r$min_match_active >= 0.7)) }
    lg <- s[nrow(s), ]; bs <- s[which.max(s$mean_C), ]
    rule <- function(target) { c1 <- s$mean_C >= target$mean_C - target$se_C
      hit <- which(c1 & vapply(s$set, ok_prog, TRUE)); if (length(hit)) s$set[min(hit)] else NA }
    data.frame(arm = a, all_genes_C = lg$mean_C, all_genes_se = lg$se_C,
               preferred_vs_all_genes = rule(lg), best_set = bs$set, best_C = bs$mean_C,
               preferred_vs_best = rule(bs))
  }))
  write.csv(pref, file.path(OUT, "size_preferred.csv"), row.names = FALSE)

  ## Elbow plot ----
  png(file.path(OUT, "size_elbow.png"), width = 1100, height = 700, res = 130)
  cols <- c(w1 = "#d95f02", w0 = "#7570b3")
  se0 <- ifelse(is.na(sz$se_C), 0, sz$se_C)
  plot(NA, xlim = range(sz$n_genes), ylim = range(c(sz$mean_C - se0, sz$mean_C + se0), na.rm = TRUE),
       log = "x", xlab = "Genes in the input (log scale)", ylab = "Mean held-out C over 5 folds (± SE)",
       main = "Training CV C by feature-set size")
  for (a in names(cols)) { s <- sz[sz$arm == a, ]; s <- s[order(s$n_genes), ]
    off <- if (a == "w1") 1.02 else 0.98
    arrows(s$n_genes * off, s$mean_C - s$se_C, s$n_genes * off, s$mean_C + s$se_C,
           angle = 90, code = 3, length = .03, col = cols[a])
    lines(s$n_genes * off, s$mean_C, type = "b", pch = 19, col = cols[a]) }
  legend("bottomright", c("w_S = 1 (survival in the loadings)", "w_S = 0 (unsupervised loadings)"),
         col = cols, lty = 1, pch = 19, bty = "n")
  dev.off()
}

# S3: programs of the D4 all-patient fits ----
progs <- do.call(rbind, lapply(c("w1", "w0"), function(a) {
  x <- ck[[key("size", "D4", length(z$sets$D4), a, 0L)]]
  if (is.null(x$state)) return(NULL)
  s <- x$state; sdb <- sqrt(pmax(s$EBeta2 - s$EBeta^2, 0)); alive <- colSums(s$EF^2) > 1e-8
  do.call(rbind, lapply(seq_len(ncol(s$EF)), function(k) {
    f <- s$EF[, k]; o <- order(f)
    data.frame(arm = a, program = k, alive = alive[k], beta = s$EBeta[k], beta_sd = sdb[k],
               z = if (sdb[k] > 0) s$EBeta[k] / sdb[k] else NA_real_,
               survival_active = alive[k] && sdb[k] > 0 && abs(s$EBeta[k]) >= 1.96 * sdb[k],
               top_positive = paste(rev(tail(x$genes[o], 10)), collapse = " "),
               top_negative = paste(head(x$genes[o], 10), collapse = " "))
  }))
}))
if (!is.null(progs)) write.csv(progs, file.path(OUT, "programs_D4.csv"), row.names = FALSE)

# Survival screening ----
sc <- summ[summ$part == "screen", ]
if (nrow(sc)) {
  paired <- do.call(rbind, lapply(split(sc, list(sc$N, sc$arm), drop = TRUE), function(g) {
    N <- g$N[1]; a <- g$arm[1]
    f <- function(set) m[m$part == "screen" & m$set == set & m$N == N & m$arm == a & m$fold > 0, c("fold", "C")]
    dd <- merge(f("S"), f("U"), by = "fold", suffixes = c("_S", "_U")); dif <- dd$C_S - dd$C_U
    ov <- vapply(1:5, function(fo) { s <- ck[[key("screen", "S", N, a, fo)]]$genes
      u <- ck[[key("screen", "U", N, a, fo)]]$genes; if (is.null(s) || is.null(u)) NA_real_ else length(intersect(s, u)) / length(u) }, 1)
    data.frame(N = N, arm = a, n_genes = g$n_genes[1],
               C_U = g$mean_C[g$set == "U"], se_U = g$se_C[g$set == "U"],
               C_S = g$mean_C[g$set == "S"], se_S = g$se_C[g$set == "S"],
               diff_S_minus_U = mean(dif, na.rm = TRUE),
               diff_se = sd(dif, na.rm = TRUE) / sqrt(sum(!is.na(dif))), paired_folds = sum(!is.na(dif)),
               mean_overlap_S_in_U = mean(ov, na.rm = TRUE),
               active_U_full = g$full_active[g$set == "U"], active_S_full = g$full_active[g$set == "S"])
  }))
  write.csv(paired[order(paired$arm, paired$N), ], file.path(OUT, "screen_paired.csv"), row.names = FALSE)
  print(paired, digits = 3, row.names = FALSE)
}
print(summ[, c("part", "set", "arm", "n_genes", "mean_C", "se_C", "pooled_oof_C", "folds_converged",
               "mean_active_cv", "full_active", "share_alive_median", "fit_sec_median")], digits = 3, row.names = FALSE)
if (exists("pref")) print(pref, row.names = FALSE)
if (exists("agree")) print(agree, digits = 3, row.names = FALSE)
