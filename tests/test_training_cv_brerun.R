# ============================================================
# Script: test_training_cv_brerun.R
# Purpose: Verify beta-prior dispatch and paired K=5 rerun summaries without patient data.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml; run_training_cv_redo.R helpers
# ============================================================
source('results/benchmark_sim/run_training_cv_redo.R')

# Prior dispatch ----
# Capture the fitter call to verify defaults and the opt-in shared prior.
#' Capture fitter arguments without running a model
#' @param ... Arguments passed to the fitter.
#' @return Named argument list.
#' @examples
#' fit_ssbmf_yfb(K = 5)
fit_ssbmf_yfb <- function(...) list(...)
cfg <- yaml::read_yaml('config/globals.yml')
x <- matrix(1:12, 4)
a <- f_fit(x, 1:4, c(1, 0, 1, 0), 5, FALSE, cfg)
stopifnot(identical(a$prior_beta, 'normal'), identical(a$beta_prior_mode, 'per_program'))
b <- f_fit(x, 1:4, c(1, 0, 1, 0), 5, FALSE, cfg,
           prior_beta = 'point_normal', beta_prior_mode = 'shared')
stopifnot(identical(b$prior_beta, 'point_normal'), identical(b$beta_prior_mode, 'shared'),
          b$K == 5, b$prior_LF == 'point_laplace', b$canonicalize, b$w_S == 1,
          b$max_iter == cfg$cavi$max_iter, b$tol == cfg$cavi$tol)

# Paired summaries ----
F_PRIV <- tempfile('brerun-fits-'); F_OUT <- tempfile('brerun-output-')
dir.create(F_PRIV); dir.create(F_OUT)
parts <- c('rank', 'b_rerun_PLAP', 'b_rerun_PN', 'b_rerun_NSH')
delta <- c(-0.03, -0.01, 0, 0.01, 0.02)
for (i in seq_along(parts)) for (fold in 0:5) {
  offset <- if (fold == 0 || i == 1) 0 else delta[fold] - (if (i == 2) 0.05 else 0)
  row <- data.frame(K = 5L, fold = fold, error = NA_character_, warnings = '',
    C_within = 0.6 + fold / 100 + offset, C_pooled = 0.65 + fold / 100 + 2 * offset,
    converged = fold != 4, n_active = fold %% 3)
  saveRDS(list(metrics = row, state = list(EBeta = if (fold == 2) rep(0, 5) else 1:5),
    pip = if (i == 3) seq(0.1, 0.9, length.out = 5) else NULL),
    file.path(F_PRIV, sprintf('%s_D4_Z_K05_f%d.rds', parts[i], fold)))
}
out <- f_b_summary(5)
stopifnot(nrow(out) == 4, all(out$n_folds == 5), all(out$folds_all_beta_zero == 1),
          all(out$folds_converged == 4), !out$criterion_c_pass[2], out$criterion_c_pass[3],
          abs(out$diff_within[2] + 0.052) < 1e-12,
          abs(out$se_within[3] - sd(delta) / sqrt(5)) < 1e-12,
          abs(out$diff_pooled[3] - 2 * mean(delta)) < 1e-12,
          out$full_pip_1[3] == 0.1, out$full_pip_5[3] == 0.9,
          is.na(out$full_pip_1[1]))
# A missing fold must fail instead of silently reducing the paired denominator.
unlink(file.path(F_PRIV, 'b_rerun_PN_D4_Z_K05_f3.rds'))
stopifnot(inherits(suppressWarnings(try(f_b_summary(5), silent = TRUE)), 'try-error'))
unlink(c(F_PRIV, F_OUT), recursive = TRUE)
cat('K=5 rerun regression checks passed (dispatch, paired estimates, rule, PIPs, missing fold).\n')
