# ============================================================
# Script: run_yfb_supervised_real_refit.R
# Purpose: Attempt S-15. Real-data YFB fit with the survival-block weight (A1), w_S chosen by
#          five-fold CV partial log-likelihood, with a nested-CV estimate of the procedure.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, ebnm, yaml, parallel; results/benchmark_sim/run_yfb_supervised_real.R
#               (s11_load, prepared input), code/fit_ssbmf_yfb.R, code/predict_ssbmf_yfb.R
# Run: Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R fit [cores]
#      Rscript results/benchmark_sim/run_yfb_supervised_real_refit.R summary
# Registration: docs/plans/2026-11-06-supervised-real-refit-plan.md
# ============================================================
# Training patients only (TCGA_PAAD, CPTAC); no external cohort is loaded. Patient-level inputs
# stay in data/ (git-ignored); checkpoints are .rds (git-ignored); only aggregates are written.
#
# Indexing. Outer level o = 0 uses all patients with the S-11 folds as its inner folds; o = 1..5
# uses the patients outside S-11 fold o, with inner folds stratified by cohort x event
# (set.seed(2026 + o)). Inner level i = 0 fits the whole outer set; i = 1..5 leaves inner fold i
# out. A fit (o, i, w) is scored on: inner fold i (i > 0), outer fold o (i = 0, o > 0), nothing
# (o = i = 0, the final all-patient fit).
#
# Parallelization: PSOCK, one BLAS thread per worker, one checkpoint per fit; a rerun skips
# finished fits. Memory: 273 x 2,218 doubles (~5 MB) per worker.

Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source("results/benchmark_sim/run_yfb_supervised_real.R")   # s11_load, LOCAL, NFOLD (main guarded)
# EXHIBIT_K (exhibit_rank.R, sourced above) reruns at another rank into a separate folder.
S15_OUT <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_refit_2026-11-06", "s15")
S15_W   <- c(1, 2, 5, 10, 20, 50)
S15_S11 <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04", "s11")

#' Load code as s11_load(), then re-source the fitter and predictor from the committed snapshot
#'
#' Another session may edit code/fit_ssbmf_yfb.R and code/predict_ssbmf_yfb.R in this folder
#' while the grid runs; the snapshot (written by the driver from `git show HEAD:`) keeps every
#' worker on the commit recorded in provenance.txt.
#' @param snap Directory holding the two snapshot files. @return NULL, invisibly.
#' @examples
#' # s15_load(tempdir())
s15_load <- function(snap) {
  s11_load()
  suppressMessages(tryCatch(source(file.path(snap, "fit_ssbmf_yfb.R")), error = function(e) {
    if (conditionMessage(e) != "!is.null(real_Y) is not TRUE") stop(e)
  }))
  source(file.path(snap, "predict_ssbmf_yfb.R"))
  invisible(NULL)
}

#' Prepared D4 data and the nested fold structure
#' @return List: Y (n x 2,218), time, status, cohort, outer folds, inner fold list (index o + 1).
#' @examples
#' # z <- s15_data()
s15_data <- function() {
  z <- readRDS(LOCAL)
  Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
  inner <- vector("list", NFOLD + 1L)
  inner[[1]] <- z$fold
  for (o in seq_len(NFOLD)) {
    f <- integer(length(z$time)); U <- which(z$fold != o)
    set.seed(2026L + o)
    for (co in unique(z$cohort)) for (st in unique(z$status)) {
      idx <- U[z$cohort[U] == co & z$status[U] == st]
      if (length(idx)) f[idx] <- sample(rep(seq_len(NFOLD), length.out = length(idx)))
    }
    inner[[o + 1L]] <- f            # 0 for patients outside the outer training set
  }
  list(Y = Y, time = z$time, status = z$status, cohort = z$cohort, outer = z$fold, inner = inner)
}

#' Fit and score one (outer, inner, weight) task, saving a checkpoint
#'
#' CV partial log-likelihood (Verweij-van Houwelingen) for i > 0:
#' l_U(eta) - l_T(eta), with eta = Y_U E[F] E[beta] from the fit on T = U minus inner fold i.
#' @param task Row with o, i, w. @param d s15_data() output. @param cfg Config.
#' @return Status string.
#' @examples
#' # s15_task(data.frame(o = 0, i = 1, w = 1), d, cfg)
s15_task <- function(task, d, cfg) {
  key <- sprintf("o%d_i%d_w%g", task$o, task$i, task$w)
  path <- file.path(S15_OUT, "fits", paste0(key, ".rds"))
  if (file.exists(path)) return(paste(key, "done earlier"))
  U <- if (task$o == 0L) seq_along(d$time) else which(d$outer != task$o)
  fi <- d$inner[[task$o + 1L]]
  Tr <- if (task$i == 0L) U else U[fi[U] != task$i]
  E <- if (task$i > 0L) U[fi[U] == task$i] else if (task$o > 0L) which(d$outer == task$o) else integer(0)
  set.seed(42L); t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(fit_ssbmf_yfb(d$Y[Tr, , drop = FALSE], d$time[Tr], d$status[Tr],
                    K = exhibit_fit_k(cfg), max_iter = cfg$cavi$max_iter,
                    tol = cfg$cavi$tol, prior_LF = "point_laplace", prior_beta = "normal",
                    canonicalize = TRUE, survival_block_weight = TRUE, w_S = task$w,
                    alpha = cfg$benchmark$alpha, verbose = FALSE, trace_supervision = TRUE),
                  error = function(e) e)
  row <- data.frame(o = task$o, i = task$i, w = task$w, n_train = length(Tr), n_test = length(E),
                    error = NA_character_)
  if (inherits(fit, "error")) {
    row$error <- conditionMessage(fit); saveRDS(list(metrics = row), path)
    return(paste(key, "ERROR", row$error))
  }
  Ytr <- d$Y[Tr, , drop = FALSE]
  sdb <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0)); alive <- colSums(fit$EF^2) > 1e-8
  active <- abs(fit$EBeta) >= 1.96 * sdb & sdb > 0 & alive
  shareA <- tail(fit$history$trace$shareA_median, 1)[1, ]
  el <- fit$history$elbo_full[seq_len(fit$history$n_iter)]
  dec <- if (length(el) > 1) max(c(0, -diff(el) / pmax(abs(head(el, -1)), 1))) else 0
  row <- cbind(row, data.frame(converged = fit$history$converged, n_iter = fit$history$n_iter,
    max_rel_elbo_decrease = dec, n_alive = sum(alive), n_active = sum(active),
    share_active_median = if (any(active)) median(shareA[active]) else NA_real_,
    share_alive_median = median(shareA[alive]), max_abs_beta = max(abs(fit$EBeta)),
    fit_sec = proc.time()[["elapsed"]] - t0, C = NA_real_, cvpl = NA_real_))
  risk <- NULL
  if (length(E)) {
    risk <- as.numeric(predict_ssbmf_yfb(d$Y[E, , drop = FALSE], fit$EF, fit$EBeta,
                                         EF_norms = fit$EF_norms)$risk_scores)
    row$C <- if (sd(risk) > 1e-8) frozen_reverse_cindex(risk, d$time[E], d$status[E]) else NA_real_
  }
  if (task$i > 0L) {
    eta <- as.numeric(predict_ssbmf_yfb(d$Y[U, , drop = FALSE], fit$EF, fit$EBeta,
                                        EF_norms = fit$EF_norms)$risk_scores)
    inT <- U %in% Tr
    row$cvpl <- calc_cox_taylor_yf(eta, d$time[U], d$status[U])$logPL -
                calc_cox_taylor_yf(eta[inT], d$time[Tr], d$status[Tr])$logPL
  }
  saveRDS(list(metrics = row, E = E, risk = risk,
               state = if (task$i == 0L) list(EF = fit$EF, EL = fit$EL, EBeta = fit$EBeta,
                 EBeta2 = fit$EBeta2, EF_norms = fit$EF_norms, shareA = shareA, Tr = Tr) else NULL), path)
  sprintf("%s ok %.0fs", key, row$fit_sec)
}

#' Chosen weight for one outer level: largest summed CV partial log-likelihood among weights
#' whose five inner fits ran without error; ties to the smaller weight
#' @param m Metrics data frame. @param o Outer level. @return The chosen weight (NA if none).
#' @examples
#' # s15_choose(m, 0)
s15_choose <- function(m, o) {
  x <- m[m$o == o & m$i > 0, ]
  ok <- tapply(is.na(x$error) & is.finite(x$cvpl), x$w, sum)
  pl <- tapply(x$cvpl, x$w, sum)
  cand <- as.numeric(names(ok)[ok == NFOLD])
  if (!length(cand)) return(NA_real_)
  pl <- pl[as.character(cand)]
  min(cand[pl == max(pl)])
}

# Summary ----

#' Write the S-15 aggregate tables
#' @return NULL; writes CSVs to S15_OUT and prints the verdict.
#' @examples
#' # s15_summary()
s15_summary <- function() {
  m <- do.call(rbind, lapply(list.files(file.path(S15_OUT, "fits"), full.names = TRUE), function(f) {
    r <- readRDS(f)$metrics
    for (nm in c("converged", "n_iter", "max_rel_elbo_decrease", "n_alive", "n_active",
                 "share_active_median", "share_alive_median", "max_abs_beta", "fit_sec", "C", "cvpl"))
      if (is.null(r[[nm]])) r[[nm]] <- NA
    r }))
  cat(sprintf("%d of %d fits; %d errors\n", nrow(m), 6 * 6 * length(S15_W), sum(!is.na(m$error))))
  write.csv(m, file.path(S15_OUT, "metrics.csv"), row.names = FALSE)
  # Weight curve on the S-11 folds (o = 0)
  x <- m[m$o == 0 & m$i > 0, ]
  curve <- do.call(rbind, lapply(split(x, x$w), function(z) data.frame(w = z$w[1],
    folds_ok = sum(is.na(z$error)), folds_converged = sum(z$converged, na.rm = TRUE),
    cvpl = sum(z$cvpl), mean_C = mean(z$C, na.rm = TRUE), se_C = sd(z$C, na.rm = TRUE) / sqrt(NFOLD),
    share_active_median = median(z$share_active_median, na.rm = TRUE),
    share_alive_median = median(z$share_alive_median, na.rm = TRUE),
    mean_active = mean(z$n_active, na.rm = TRUE), median_iter = median(z$n_iter, na.rm = TRUE))))
  write.csv(curve, file.path(S15_OUT, "weight_curve.csv"), row.names = FALSE)
  print(curve, digits = 3, row.names = FALSE)
  # Reproducibility check against S-11 w_S = 1
  s11 <- read.csv(file.path(S15_S11, "fits.csv"))
  s11 <- s11[s11$set == "D4" & s11$arm == "w1" & s11$fold > 0, c("fold", "C")]
  rep <- merge(s11, x[x$w == 1, c("i", "C")], by.x = "fold", by.y = "i", suffixes = c("_s11", "_s15"))
  write.csv(rep, file.path(S15_OUT, "reproduce_s11.csv"), row.names = FALSE)
  cat(sprintf("Max |C difference| vs S-11 at w = 1: %.2e\n", max(abs(rep$C_s11 - rep$C_s15))))
  # Nested CV
  nest <- do.call(rbind, lapply(seq_len(NFOLD), function(o) {
    w <- s15_choose(m, o)
    r <- if (is.na(w)) NULL else m[m$o == o & m$i == 0 & m$w == w, ]
    data.frame(outer = o, chosen_w = w, C = if (length(r) && nrow(r)) r$C else NA_real_)
  }))
  write.csv(nest, file.path(S15_OUT, "nested_cv.csv"), row.names = FALSE)
  print(nest)
  wstar <- s15_choose(m, 0)
  full <- m[m$o == 0 & m$i == 0 & m$w == wstar, ]
  nested_C <- mean(nest$C, na.rm = TRUE); nested_se <- sd(nest$C, na.rm = TRUE) / sqrt(NFOLD)
  # Registered S-11 references (K = 7); at another rank they are read from that rank's S-11
  # summary (mean fold C and its SE of the D4 w_S = 1 and two-step arms).
  ref <- c(w1 = 0.663, twostep = 0.652, se = 0.024)
  if (!is.na(exhibit_k())) {
    sm <- read.csv(file.path(S15_S11, "summary.csv"))
    sm <- sm[sm$part == "size" & sm$set == "D4", ]
    ref <- c(w1 = sm$mean_C[sm$arm == "w1"], twostep = sm$mean_C[sm$arm == "w0cox"],
             se = sm$se_C[sm$arm == "w1"])
    stopifnot(length(ref) == 3L, all(is.finite(ref)))
  }
  verdict <- data.frame(w_star = wstar, nested_mean_C = nested_C, nested_se_C = nested_se,
    s11_w1_C = ref[["w1"]], s11_twostep_C = ref[["twostep"]], threshold = ref[["w1"]] - ref[["se"]],
    full_fit_ok = nrow(full) == 1 && is.na(full$error) &&
      (isTRUE(full$converged) || full$max_rel_elbo_decrease <= 1e-6),
    nested_ok = nested_C >= ref[["w1"]] - ref[["se"]])
  verdict$adopt <- verdict$w_star > 1 && verdict$full_fit_ok && verdict$nested_ok
  write.csv(verdict, file.path(S15_OUT, "verdict.csv"), row.names = FALSE)
  print(verdict)
  # Programs of the full fit at w* against the S-11 unsupervised loadings
  ck <- readRDS(file.path(S15_OUT, "fits", sprintf("o0_i0_w%g.rds", wstar)))$state
  d <- s15_data(); Ytr <- d$Y[ck$Tr, , drop = FALSE]
  w0 <- readRDS(file.path(S15_S11, "fits", "size_D4_N2218_w0_f0.rds"))$state$EF
  Z <- Ytr %*% ck$EF; eta <- drop(Z %*% ck$EBeta); sdb <- sqrt(pmax(ck$EBeta2 - ck$EBeta^2, 0))
  pve <- colSums(ck$EL^2) * colSums(ck$EF^2) / sum(Ytr^2)
  act <- abs(ck$EBeta) >= 1.96 * sdb & sdb > 0 & colSums(ck$EF^2) > 1e-8
  rr <- suppressWarnings(abs(cor(ck$EF, w0))); rr[!is.finite(rr)] <- 0
  prog <- data.frame(k = seq_along(pve), pve = pve, beta = ck$EBeta, beta_sd = sdb, active = act,
    role = ifelse(act, ifelse(pve >= .01, "prognostic_program", "survival_only_direction"),
                  ifelse(pve >= .01, "genomics_only_program", "inactive")),
    std_beta = abs(ck$EBeta) * apply(Z, 2, sd),
    risk_var_share = apply(Z, 2, sd)^2 * ck$EBeta^2 / var(eta),
    survival_precision_share = ck$shareA, best_r_unsupervised = apply(rr, 1, max))
  write.csv(prog, file.path(S15_OUT, "programs_full_fit.csv"), row.names = FALSE)
  print(prog, digits = 3)
  invisible(NULL)
}

# Main ----
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE); mode <- args[1]
  cores <- if (length(args) > 1) as.integer(args[2]) else 9L
  dir.create(file.path(S15_OUT, "fits"), recursive = TRUE, showWarnings = FALSE)
  if (mode == "fit") {
    snap <- file.path(tempdir(), "s15_code"); dir.create(snap, showWarnings = FALSE)
    for (f in c("fit_ssbmf_yfb.R", "predict_ssbmf_yfb.R"))
      system2("git", c("show", paste0("HEAD:code/", f)), stdout = file.path(snap, f))
    writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()),
                 "fit_ssbmf_yfb.R and predict_ssbmf_yfb.R sourced from this commit"),
               file.path(S15_OUT, "provenance.txt"))
    design <- expand.grid(i = 0:NFOLD, o = 0:NFOLD, w = S15_W)
    design <- design[order(-design$w, design$o, design$i), ]   # slow large weights first
    cl <- parallel::makePSOCKcluster(cores, outfile = file.path(S15_OUT, "workers.log"))
    parallel::clusterExport(cl, c("S15_OUT", "LOCAL", "NFOLD", "s11_load", "s15_load", "s15_data",
                                  "s15_task", "snap", "exhibit_k", "exhibit_fit_k"), envir = environment())
    parallel::clusterEvalQ(cl, { Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1); s15_load(snap); d <- s15_data()
      cfg <- yaml::read_yaml("config/globals.yml"); NULL })
    st <- parallel::parLapplyLB(cl, split(design, seq_len(nrow(design))), function(t)
      tryCatch(s15_task(t, d, cfg), error = function(e) paste("TASK ERROR", conditionMessage(e))),
      chunk.size = 1L)
    parallel::stopCluster(cl)
    writeLines(unlist(st), file.path(S15_OUT, "task_status.txt"))
  }
  if (mode == "summary") { s11_load(); s15_summary() }
}
