# ============================================================
# Script: summarize_balance.R
# Purpose: Compare the balance options of rebuild Part 4 with w = 1 and the
#          two-step comparator on the low-variance prognostic cells (O4),
#          and check that O1-O3 and O7 still hold under each option.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: base R
# Usage: Rscript results/ssbmf/summarize_balance.R   (from the repository root)
# ============================================================

base <- "results/ssbmf/outputs"
read_jobs <- function(dir) {
  f <- list.files(file.path(base, dir, "jobs"), full.names = TRUE)
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, function(x) { r <- readRDS(x); r[, intersect(names(r), cols)] }))
}
w1 <- read.csv(file.path(base, "sim_panel", "sim_panel_rows.csv"))
w1 <- w1[w1$cell %in% c("C", "C2", "C3") & w1$arm %in% c("ssbmf_point_laplace", "two_step"), ]
w1$w_surv <- ifelse(w1$arm == "two_step", 0, 1); w1$block_weight <- 1
cols <- names(w1)
rows <- rbind(w1, read_jobs("balance_eff"), read_jobs("balance_cvw"))
# compare like with like: cell C2 has only K_init = 6 for the CV arm
rows$fit_id <- paste(rows$cell, rows$seed, rows$K_init, rows$arm)

summ <- do.call(rbind, lapply(split(rows, list(rows$cell, rows$arm), drop = TRUE), function(x) {
  f <- x[!duplicated(x$fit_id), ]
  lv <- x[x$beta_true != 0, ]
  hv <- x[x$pve_true >= 0.02, ]
  cov_rows <- x[(x$pve_true >= 0.02 & x$best_r >= 0.9) | (x$pve_true < 0.02 & x$best_r >= 0.7), ]
  data.frame(cell = x$cell[1], arm = x$arm[1], fits = nrow(f),
             K_init = paste(sort(unique(f$K_init)), collapse = ","),
             O4_recovered = sprintf("%d/%d", sum(lv$best_r >= 0.7 & lv$active & lv$sign_ok), nrow(lv)),
             lowvar_matched = mean(lv$best_r >= 0.7), lowvar_best_r = mean(lv$best_r),
             lowvar_active = mean(lv$active),
             O1_no_decrease = mean(f$step_decreases == 0), O2_converged = mean(f$converged),
             O3_highvar = mean(hv$best_r >= 0.9),
             O7_coverage90 = if (nrow(cov_rows)) mean(cov_rows$covered) else NA,
             mean_val_C = mean(f$val_C),
             w_surv_chosen = paste(names(table(f$w_surv)), table(f$w_surv), sep = ":", collapse = " "),
             median_block_weight = median(f$block_weight),
             rejections_median = median(f$rejections),
             runtime_median_s = median(f$runtime_s))
}))
summ <- summ[order(summ$cell, summ$arm), ]
write.csv(summ, file.path(base, "balance_summary.csv"), row.names = FALSE)
options(width = 220)
print(summ, digits = 3, row.names = FALSE)

# C2 restricted to K_init = 6 for a like-for-like comparison with the CV arm
c2 <- rows[rows$cell == "C2" & rows$K_init == 6 & rows$beta_true != 0, ]
if (nrow(c2)) {
  cat("\nCell C2, K_init = 6 only:\n")
  print(aggregate(cbind(recovered = best_r >= 0.7 & active & sign_ok, best_r = best_r) ~ arm, data = c2, FUN = mean))
}
