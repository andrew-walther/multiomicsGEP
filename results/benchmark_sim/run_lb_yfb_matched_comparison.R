# ============================================================
# Script: run_lb_yfb_matched_comparison.R
# Purpose: Compare LB and YFB on the same planted-program simulations at
#          matched survival-to-genomics precision weights, with each model's
#          two-step (unsupervised fit, then Cox) reference.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, ebnm, yaml; code/fit_ssbmf_lb.R,
#               code/predict_ssbmf_lb.R, code/fit_ssbmf_yfb.R,
#               code/predict_ssbmf_yfb.R
# Usage:  Rscript results/benchmark_sim/run_lb_yfb_matched_comparison.R
#         (about 10 minutes on a laptop; single process, no HPC needed)
# ============================================================
#
# Matching. The LB L update has precision (1 - alpha) A_gen + alpha A_surv =
# (1 - alpha) [A_gen + w A_surv] with w = alpha / (1 - alpha). YFB with
# survival_block_weight = TRUE multiplies the Cox term by w_S. So LB at
# alpha = w / (1 + w) and YFB at w_S = w weight survival against genomics
# equally (LB also tempers everything by 1 - alpha). w in {1, 10}.
#
# Data. Planted programs: Y = L F' + noise, eta = L beta with
# beta = (1, -0.8, 0.5, 0), n = 200 train / 100 test, p = 400, K = 4, two
# signal amplitudes, six seeds each (the LB review design). Each model gets
# its own preprocessing: LB the raw nonnegative Y with point-exponential
# priors (its design), YFB training-standardized Y with point-Laplace priors
# and score standardization (the S-11/S-13 configuration).
#
# Scores. LB test scores are fold-in loadings (predict_ssbmf_lb_foldin());
# YFB test scores are the program scores Y* F (exact, so score_converged is
# TRUE). Recovery is the mean |r| between true test loadings of the three
# survival programs and one-to-one matched estimated scores. LB's training C
# uses its in-sample, outcome-informed loadings and its test C uses fold-in
# loadings, so train-test gaps are not comparable across models. This is an exploratory comparison, not a
# registered study; it does not change either model's configuration.

suppressPackageStartupMessages({ library(survival); library(ebnm); library(yaml) })
suppressMessages(tryCatch(source("code/fit_ssbmf_lb.R"), error = function(e) NULL))
suppressMessages(tryCatch(source("code/fit_ssbmf_yfb.R"), error = function(e) NULL))
source("code/predict_ssbmf_lb.R"); source("code/predict_ssbmf_yfb.R")
cfg <- yaml::read_yaml("config/globals.yml")
OUT <- "results/benchmark_sim/outputs/lb_yfb_matched_2026-11-06"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Helpers ----

#' Simulate planted-program LB data (review design)
#' @param seed RNG seed. @param amp Planted loading amplitude in F.
#' @return list(Y, time, status, L, beta, tr, te).
#' @examples
#' # d <- sim_planted(1, 0.5)
sim_planted <- function(seed, amp, ntr = 200, nte = 100, p = 400, K = 4,
                        beta = c(1.0, -0.8, 0.5, 0)) {
  set.seed(seed); n <- ntr + nte
  L <- matrix(rexp(n * K), n, K); Fm <- matrix(0, p, K)
  for (k in 1:K) Fm[sample.int(p, round(0.1 * p)), k] <- amp
  tau <- rgamma(p, 2, 2)
  Y <- L %*% t(Fm) + sweep(matrix(rnorm(n * p), n, p), 2, sqrt(tau), "/") * 0.5
  eta <- as.vector(L %*% beta)
  ev <- (-log(runif(n)) / (0.01 * exp(eta)))^(1 / 1.5)
  cen <- rexp(n, 1) * median(ev) * 1.5
  list(Y = Y, time = pmin(ev, cen), status = as.integer(ev <= cen),
       L = L, beta = beta, tr = 1:ntr, te = (ntr + 1):n)
}

#' Cox-direction concordance (larger score = higher hazard)
#' @param eta Risk score. @param time,status Outcomes. @return C-index or NA.
#' @examples
#' # cidx(rnorm(10), rexp(10), rep(1, 10))
cidx <- function(eta, time, status) {
  if (stats::sd(eta) == 0) return(0.5)
  as.numeric(concordance(Surv(time, status) ~ eta, reverse = TRUE)$concordance)
}

#' Mean |r| between true survival-program loadings and one-to-one matched scores
#'
#' Greedy one-to-one matching on |r| over all K true programs and K estimated
#' columns, so one estimated column cannot stand in for several programs; a
#' program matched to a dead column scores 0. Averaged over `progs`.
#' @param Ltrue n x K true loadings. @param S n x K estimated scores.
#' @param progs Columns of Ltrue to score (nonzero beta). @return Mean |r|.
#' @examples
#' # recovery(matrix(rnorm(20), 10), matrix(rnorm(20), 10), 1:2)
recovery <- function(Ltrue, S, progs) {
  R <- suppressWarnings(abs(cor(Ltrue, S)))
  R[!is.finite(R)] <- 0
  best <- numeric(ncol(Ltrue)); M <- R
  for (i in seq_len(min(dim(M)))) {
    ij <- which(M == max(M), arr.ind = TRUE)[1, ]
    best[ij[1]] <- R[ij[1], ij[2]]
    M[ij[1], ] <- -1; M[, ij[2]] <- -1
  }
  mean(best[progs])
}

#' Evaluate an expression, counting (not hiding) its warnings
#' @param expr Expression. @return list(value, n_warn).
#' @examples
#' # with_warn_count(log(-1))
with_warn_count <- function(expr) {
  n <- 0L
  value <- withCallingHandlers(expr, warning = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") })
  list(value = value, n_warn = n)
}

#' Fit a Cox model on training scores and return test risk (two-step)
#' @param Str,Ste Training/test scores. @param time,status Training outcomes.
#' @return Test risk scores.
#' @examples
#' # twostep_risk(Str, Ste, time, status)
twostep_risk <- function(Str, Ste, time, status) {
  live <- which(apply(Str, 2, stats::sd) > 1e-8)
  if (!length(live)) return(rep(0, nrow(Ste)))
  cf <- coef(coxph(Surv(time, status) ~ Str[, live, drop = FALSE],
                   control = coxph.control(timefix = FALSE)))
  cf[!is.finite(cf)] <- 0
  as.vector(Ste[, live, drop = FALSE] %*% cf)
}

# Runs ----
W <- c(1, 10); K <- 4; progs <- 1:3
rows <- list()
for (amp in c(0.5, 1.0)) for (seed in 1:6) {
  d <- sim_planted(seed, amp); tr <- d$tr; te <- d$te
  ttr <- d$time[tr]; str_ <- d$status[tr]
  # YFB preprocessing: standardize with training mean/sd
  mu <- colMeans(d$Y[tr, ]); sdv <- apply(d$Y[tr, ], 2, sd)
  Zs <- sweep(sweep(d$Y, 2, mu, "-"), 2, sdv, "/")

  record <- function(model, w, fit, Ste, Str, risk_te, risk_tr, secs, n_warn, score_conv) {
    h <- fit$history
    data.frame(amp, seed, model, w, converged = isTRUE(h$converged), n_iter = h$n_iter,
               n_warn = n_warn, score_converged = score_conv,
               secs = round(secs, 1), live = sum(colSums(fit$EF^2) > 1e-8),
               rec_test = recovery(d$L[te, ], Ste, progs),
               C_train = cidx(risk_tr, ttr, str_),
               C_test = cidx(risk_te, d$time[te], d$status[te]))
  }

  ## LB: two-step (alpha = 0) and supervised at alpha = w / (1 + w) ----
  for (w in c(0, W)) {
    t0 <- proc.time()[3]
    fw <- with_warn_count(fit_ssbmf_lb(d$Y[tr, ], ttr, str_, K = K, alpha = w / (1 + w),
                                       max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                                       verbose = FALSE))
    fit <- fw$value
    secs <- proc.time()[3] - t0
    pw <- with_warn_count(predict_ssbmf_lb_foldin(d$Y[te, ], fit)); fi <- pw$value
    nw <- fw$n_warn + pw$n_warn
    if (w == 0) {
      risk_te <- twostep_risk(fit$EL, fi$L_test, ttr, str_)
      risk_tr <- twostep_risk(fit$EL, fit$EL, ttr, str_)
      rows[[length(rows) + 1]] <- record("LB", "two-step", fit, fi$L_test, fit$EL, risk_te, risk_tr,
                                         secs, nw, fi$converged)
    } else {
      rows[[length(rows) + 1]] <- record("LB", w, fit, fi$L_test, fit$EL, fi$risk_scores,
                                         as.vector(fit$EL %*% fit$EBeta), secs, nw, fi$converged)
    }
  }

  ## YFB: two-step (w_S = 0) and supervised at w_S = w ----
  for (w in c(0, W)) {
    t0 <- proc.time()[3]
    fw <- with_warn_count(fit_ssbmf_yfb(Zs[tr, ], ttr, str_, K = K,
                                          max_iter = cfg$cavi$max_iter, tol = cfg$cavi$tol,
                                          prior_LF = "point_laplace", prior_beta = "normal",
                                          canonicalize = TRUE, w_S = w,
                                          survival_block_weight = w > 0,
                                          alpha = cfg$benchmark$alpha, verbose = FALSE))
    fit <- fw$value
    secs <- proc.time()[3] - t0
    Ste <- Zs[te, ] %*% fit$EF; Str <- Zs[tr, ] %*% fit$EF
    if (w == 0) {
      rows[[length(rows) + 1]] <- record("YFB", "two-step", fit, Ste, Str,
                                         twostep_risk(Str, Ste, ttr, str_),
                                         twostep_risk(Str, Str, ttr, str_), secs, fw$n_warn, TRUE)
    } else {
      p_te <- predict_ssbmf_yfb(Zs[te, ], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
      p_tr <- predict_ssbmf_yfb(Zs[tr, ], fit$EF, fit$EBeta, fit$EF_norms)$risk_scores
      rows[[length(rows) + 1]] <- record("YFB", w, fit, Ste, Str, p_te, p_tr, secs, fw$n_warn, TRUE)
    }
  }
  cat(sprintf("amp %.1f seed %d done\n", amp, seed))
}

res <- do.call(rbind, rows)
res$w <- as.character(res$w)
write.csv(res, file.path(OUT, "matched_comparison_fits.csv"), row.names = FALSE)

# Summary by signal level, model and weight ----
sm <- aggregate(cbind(converged, score_converged, n_warn, n_iter, live, rec_test, C_train, C_test, secs) ~ amp + model + w,
                data = res, FUN = mean)
sm$gap_train_test <- sm$C_train - sm$C_test
sm <- sm[order(sm$amp, sm$model, sm$w), ]
write.csv(sm, file.path(OUT, "matched_comparison_summary.csv"), row.names = FALSE)
print(format(sm, digits = 3), row.names = FALSE)
cat("True-score test C (oracle):\n")
for (amp in c(0.5, 1.0)) {
  cs <- sapply(1:6, function(s) { d <- sim_planted(s, amp); cidx(as.vector(d$L[d$te, ] %*% d$beta), d$time[d$te], d$status[d$te]) })
  cat(sprintf("  amp %.1f: %.3f\n", amp, mean(cs)))
}
