# ============================================================
# Script: run_methylation_single_modality.R
# Purpose: Compare fixed methylation-only YFB fits with matched expression.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: survival, ebnm, yaml, Rcpp, parallel, impute; existing YFB/data helpers
# ============================================================
# Run from the repository root: Rscript <this script> prepare|smoke|run|summary
# Working-rank expression first: Rscript <this script> rerun --K=5
# Then all CpG counts: Rscript <this script> sweep --K=5
# Registration: docs/prompts/2026-11-12-codex-c-d-e-prompts.md, Prompt 1.
# Patient matrices, checkpoints and scores remain under ignored data/.
# At most three concurrent fit processes, each with one BLAS/OpenMP thread.
# 100k/all-CpG fits run serially. For Longleaf, request one CPU per fit and
# allow several GB per process for input copies, posterior moments and residuals.
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
C_LOCAL <- "data/methylation_single_modality_2026-11-12"
C_OUT <- "results/multimodal_real/outputs/methylation_single_modality_2026-11-12"
C_SCRIPT <- "results/multimodal_real/run_methylation_single_modality.R"
C_BASE_LOCAL <- C_LOCAL
C_BASE_OUT <- C_OUT
C_K <- NULL

#' Select an isolated rank without changing the historical config key
#' @param K Positive integer rank; NULL keeps the historical directory/default.
#' @return NULL invisibly; sets rank-specific private and aggregate paths.
#' @examples
#' # c_set_rank(5)
c_set_rank <- function(K = NULL) {
  historical <- yaml::read_yaml("config/globals.yml")$benchmark$k_merged_yfb_desurv
  if (is.null(K)) K <- historical
  stopifnot(length(K) == 1L, is.finite(K), K > 0L, K == as.integer(K))
  C_K <<- as.integer(K)
  C_LOCAL <<- if (K == historical) C_BASE_LOCAL else file.path(C_BASE_LOCAL, sprintf("K%02d", K))
  C_OUT <<- if (K == historical) C_BASE_OUT else file.path(C_BASE_OUT, sprintf("K%02d", K))
  dir.create(C_LOCAL, recursive = TRUE, showWarnings = FALSE)
  dir.create(C_OUT, recursive = TRUE, showWarnings = FALSE)
  if (K != historical) {
    design <- file.path(C_LOCAL, "design.rds")
    original <- normalizePath(file.path(C_BASE_LOCAL, "design.rds"), mustWork = TRUE)
    if (!file.exists(design) && !file.symlink(original, design)) stop("Could not link frozen design.")
    if (!identical(unname(tools::md5sum(design)), unname(tools::md5sum(original))))
      stop("Rank-specific design differs from the registered design.")
  }
  invisible(NULL)
}

# Pure analysis helpers ----

#' Estimate feature scaling using only fitting patients
#' PL uses (Y_ij - mean_tr,j) / SD_tr,j. PE uses Y_ij / SD_tr,j.
#' @param Y Numeric patients by features matrix.
#' @param arm PL for z-scores or PE for division by sample SD.
#' @return List with transformed Y, center and sd vectors.
#' @examples
#' c_scale(matrix(1:12, 4), "PL")
c_scale <- function(Y, arm) {
  stopifnot(arm %in% c("PL", "PE"), all(is.finite(Y)))
  s <- apply(Y, 2, stats::sd)
  if (any(!is.finite(s) | s <= 0)) stop("Constant/nonfinite feature in fitting patients.")
  center <- if (arm == "PL") colMeans(Y) else rep(0, ncol(Y))
  list(Y = sweep(sweep(Y, 2, center), 2, s, "/"), center = center, sd = s)
}

#' Apply frozen fitting-patient scaling to new patients
#' @param Y Numeric patients by features matrix in training feature order.
#' @param scaling Result from c_scale.
#' @return Numeric transformed matrix.
#' @examples
#' c_transform(matrix(1:12, 4), c_scale(matrix(1:12, 4), "PL"))
c_transform <- function(Y, scaling) {
  sweep(sweep(Y, 2, scaling$center), 2, scaling$sd, "/")
}

#' Allocate one reproducible event-balanced five-fold split
#' @param event Binary event vector.
#' @param seed Integer random seed.
#' @return Integer fold vector in patient order.
#' @examples
#' c_folds(rep(0:1, each = 10))
c_folds <- function(event, seed = 2026L) {
  set.seed(seed)
  fold <- integer(length(event))
  for (e in unique(event)) {
    i <- which(event == e)
    fold[i] <- sample(rep(1:5, length.out = length(i)))
  }
  fold
}

#' Classify exactly nonzero loadings and credible survival coefficients
#' Activity requires |E[β_k]| ≥ 1.96 sqrt(E[β_k²] - E[β_k]²), positive
#' posterior variance and nonzero feature loadings.
#' @param fit YFB fit containing EF, EBeta and EBeta2.
#' @return List of logical alive and active vectors.
#' @examples
#' # c_activity(fit)
c_activity <- function(fit) {
  alive <- colSums(fit$EF^2) > 0
  s <- sqrt(pmax(fit$EBeta2 - fit$EBeta^2, 0))
  list(alive = alive, active = alive & s > 0 & abs(fit$EBeta) >= 1.96 * s)
}

#' Match programs by maximum total absolute score correlation
#' Exhaustive one-to-one assignment is inexpensive at the registered K = 7.
#' Constant programs contribute zero. Correlations use common fitting patients.
#' @param full_scores,fold_scores Patient by K score matrices.
#' @return Integer fold column for each full-fit column.
#' @examples
#' c_match(cbind(1:5, (1:5)^2), cbind((1:5)^2, -(1:5)))
c_match <- function(full_scores, fold_scores) {
  r <- abs(suppressWarnings(stats::cor(full_scores, fold_scores)))
  r[!is.finite(r)] <- 0
  k <- ncol(r)
  permutations <- matrix(integer(0), 1, 0)
  for (j in seq_len(k)) {
    permutations <- do.call(rbind, lapply(seq_len(nrow(permutations)), function(i) {
      available <- setdiff(seq_len(k), permutations[i, ])
      cbind(matrix(rep(permutations[i, ], each = length(available)),
                   nrow = length(available)), available)
    }))
  }
  score <- apply(permutations, 1, function(p) sum(r[cbind(seq_len(k), p)]))
  as.integer(permutations[which.max(score), ])
}

#' Summarize five paired fold differences without dropping failed folds
#' For d_f = C_methylation,f - C_expression,f, SE = sd(d_f) / sqrt(5).
#' @param meth,expression Numeric concordance vectors ordered by folds 1 to 5.
#' @return Named mean difference and SE, both NA unless all five are defined.
#' @examples
#' c_paired(rep(.6, 5), c(.5, .6, .5, .6, .5))
c_paired <- function(meth, expression) {
  stopifnot(length(meth) == 5, length(expression) == 5)
  d <- meth - expression
  if (any(!is.finite(d))) return(c(mean_difference = NA_real_, se = NA_real_))
  c(mean_difference = mean(d), se = stats::sd(d) / sqrt(5))
}

# Preparation ----

#' Freeze outcome-free CpG screening, D4 genes and shared folds
#' @return NULL invisibly, with private design and aggregate screening/fold counts.
#' @examples
#' # c_prepare()
c_prepare <- function() {
  path <- file.path(C_LOCAL, "design.rds")
  if (file.exists(path)) stop("Registered design already exists.")
  u <- readRDS("data/round2_tcga_feature_universe.rds")
  source("code/cpg_feature_sensitivity.R")
  validate_cpg_reconstruction(u$methylation, u$baseline$Y$methylation)
  obs <- u$observed_summary
  stopifnot(identical(obs$probe, colnames(u$methylation)))
  # The saved mask uses observed β means before imputation and reproduces
  # the registered 182,934 count. Rank eligible imputed β by sample variance,
  # following cpg_count_manifest's decreasing variance/sort convention.
  eligible <- obs$probe[obs$keep_mean]
  stopifnot(length(eligible) == 182934L, nrow(u$methylation) == 144L)
  v <- apply(u$methylation[, eligible, drop = FALSE], 2, stats::var)
  stopifnot(all(is.finite(v) & v > 0))
  ranked <- names(sort(v, decreasing = TRUE))
  d4 <- readRDS("data/yfb_supervised_real_2026-11-04.rds")$d4
  stopifnot(all(d4 %in% colnames(u$expression_tcga_raw)),
            identical(rownames(u$methylation), rownames(u$expression_tcga_raw)),
            identical(rownames(u$methylation), names(u$time)))
  counts <- c(10000L, 25000L, 50000L, 100000L, length(ranked))
  z <- list(methylation = u$methylation, expression = u$expression_tcga_raw[, d4, drop = FALSE],
            ranked = ranked, counts = counts, time = u$time, event = u$event,
            fold = c_folds(u$event), seed = 2026L)
  saveRDS(z, path)
  write.csv(data.frame(universe = ncol(u$methylation), eligible_observed_mean = length(eligible),
    eligible_imputed_mean = sum(colMeans(u$methylation) >= .1 & colMeans(u$methylation) <= .9),
    n = length(z$time), events = sum(z$event), D4_genes = length(d4)),
    file.path(C_OUT, "screen.csv"), row.names = FALSE)
  write.csv(as.data.frame(table(event = z$event, fold = z$fold)),
            file.path(C_OUT, "fold_counts.csv"), row.names = FALSE)
  write.csv(data.frame(path = c("data/round2_tcga_feature_universe.rds",
    "data/yfb_supervised_real_2026-11-04.rds", "config/globals.yml"),
    md5 = unname(tools::md5sum(c("data/round2_tcga_feature_universe.rds",
    "data/yfb_supervised_real_2026-11-04.rds", "config/globals.yml")))),
    file.path(C_OUT, "input_hashes.csv"), row.names = FALSE)
  invisible(NULL)
}

#' Freeze committed dependencies and the reviewed runner before fitting
#' The private snapshot excludes concurrent uncommitted changes in other tasks.
#' @return NULL invisibly, with immutable source/config files and source hashes.
#' @examples
#' # c_freeze()
c_freeze <- function() {
  snapshot <- file.path(C_LOCAL, "source")
  if (dir.exists(snapshot)) {
    manifest <- file.path(C_OUT, "source_hashes.csv")
    if (!file.exists(manifest)) stop("Incomplete private source snapshot: missing manifest.")
    hashes <- read.csv(manifest, stringsAsFactors = FALSE)
    required <- c("code/fit_ssbmf_yfb.R", "code/multimodal_yfb_sweep.cpp",
                  "config/globals.yml", "runner.R")
    if (!all(required %in% hashes$path)) stop("Incomplete private source manifest.")
    if (!identical(unname(tools::md5sum(file.path(snapshot, hashes$path))), hashes$md5))
      stop("Private source snapshot is incomplete or has changed.")
    return(invisible(NULL))
  }
  dir.create(snapshot, recursive = TRUE)
  if (!is.null(C_K) && C_LOCAL != C_BASE_LOCAL) {
    baseline <- file.path(C_BASE_LOCAL, "source")
    hashes <- read.csv(file.path(C_BASE_OUT, "source_hashes.csv"), stringsAsFactors = FALSE)
    if (!identical(unname(tools::md5sum(file.path(baseline, hashes$path))), hashes$md5))
      stop("Historical source snapshot differs from its manifest.")
    if (!all(file.copy(file.path(baseline, c("code", "config")), snapshot, recursive = TRUE)))
      stop("Could not copy historical model dependencies.")
  } else {
  archive <- file.path(snapshot, "dependencies.tar")
  status <- system2("git", c("archive", "--format=tar",
    paste0("--output=", shQuote(archive)), "HEAD", "code", "config"))
  if (status != 0L) stop("Could not freeze committed model dependencies.")
  status <- system2("tar", c("-xf", shQuote(archive), "-C", shQuote(snapshot)))
  if (status != 0L) stop("Could not extract frozen model dependencies.")
  unlink(archive)
  }
  if (!file.copy(C_SCRIPT, file.path(snapshot, "runner.R"))) stop("Could not freeze runner.")
  revision <- if (!is.null(C_K) && C_LOCAL != C_BASE_LOCAL)
    read.csv(file.path(C_BASE_OUT, "source_revision.csv"))$revision else
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  write.csv(data.frame(revision = revision), file.path(C_OUT, "source_revision.csv"), row.names = FALSE)
  packages <- c("survival", "ebnm", "yaml", "Rcpp", "impute")
  versions <- sapply(packages, utils::packageVersion, simplify = FALSE)
  write.csv(data.frame(package = packages, version = unlist(lapply(versions, as.character))),
            file.path(C_OUT, "package_versions.csv"), row.names = FALSE)
  files <- list.files(snapshot, recursive = TRUE, full.names = TRUE)
  write.csv(data.frame(path = substring(files, nchar(snapshot) + 2),
    md5 = unname(tools::md5sum(files))), file.path(C_OUT, "source_hashes.csv"), row.names = FALSE)
  invisible(NULL)
}

# Individual isolated fits ----

#' Load fitting and scoring helpers in a child process
#' @return NULL invisibly.
#' @examples
#' # c_load()
c_load <- function() {
  snapshot <- normalizePath(file.path(C_LOCAL, "source"), mustWork = TRUE)
  previous <- getwd()
  setwd(snapshot)
  on.exit(setwd(previous))
  source("code/fit_ssbmf_yfb.R")
  source("code/predict_ssbmf_yfb.R")
  source("code/concordance_ci.R")
  if (!multimodal_yfb_use_cpp()) stop("Compiled sweep required for the registered large study.")
  invisible(NULL)
}

#' Fit and checkpoint one registered matrix, setup and fold
#' @param modality methylation or expression.
#' @param count Number of CpGs, ignored for expression.
#' @param arm PL or PE.
#' @param fold Zero for full fit, otherwise held-out fold 1 to 5.
#' @return NULL invisibly. Private checkpoint includes failures and warnings.
#' @examples
#' # c_one("methylation", 10000, "PL", 0)
c_one <- function(modality, count, arm, fold) {
  key <- sprintf("%s_%d_%s_f%d", modality, count, arm, fold)
  path <- file.path(C_LOCAL, paste0(key, ".rds"))
  if (file.exists(path)) return(invisible(NULL))
  c_load()
  z <- readRDS(file.path(C_LOCAL, "design.rds"))
  cfg <- yaml::read_yaml(file.path(C_LOCAL, "source/config/globals.yml"))
  raw <- if (modality == "expression") z$expression else
    z$methylation[, z$ranked[seq_len(count)], drop = FALSE]
  tr <- if (fold == 0L) seq_along(z$time) else which(z$fold != fold)
  te <- which(z$fold == fold)
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  answer <- tryCatch(withCallingHandlers({
    scaling <- c_scale(raw[tr, , drop = FALSE], arm)
    set.seed(42L)
    fit <- fit_ssbmf_yfb(scaling$Y, z$time[tr], z$event[tr],
      K = if (is.null(C_K)) cfg$benchmark$k_merged_yfb_desurv else C_K, max_iter = cfg$cavi$max_iter,
      tol = cfg$cavi$tol, prior_LF = if (arm == "PL") "point_laplace" else "point_exponential",
      prior_beta = "normal", canonicalize = TRUE, w_S = 1,
      alpha = cfg$benchmark$alpha, center_projection = arm == "PE", prune = FALSE,
      verbose = FALSE)
    pred <- predict_ssbmf_yfb(c_transform(raw, scaling), fit$EF, fit$EBeta,
                             fit$EF_norms, projection_center = fit$projection_center)
    act <- c_activity(fit)
    row <- data.frame(modality, count, arm, fold, n_features = ncol(raw),
      status = if (fit$history$converged) "converged" else "capped",
      converged = fit$history$converged, iterations = fit$history$n_iter,
      seconds = proc.time()[["elapsed"]] - start,
      nonzero = sum(act$alive), survival_active = sum(act$active),
      C = if (fold > 0 && stats::sd(pred$risk_scores[te]) > 1e-8)
        frozen_reverse_cindex(pred$risk_scores[te], z$time[te], z$event[te]) else NA_real_,
      error = NA_character_)
    # Keep full posterior and full-row scores privately for reproducibility and
    # permutation-aware fold matching. External outcomes never enter this stage.
    list(row = row, fit = fit, scaling = scaling[c("center", "sd")],
         scores = pred$ZF_test, risk = pred$risk_scores, tr = tr, te = te,
         features = colnames(raw), activity = act)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    message("FIT WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(row = data.frame(modality, count, arm, fold,
    n_features = ncol(raw), status = "failed", converged = FALSE, iterations = NA_integer_,
    seconds = proc.time()[["elapsed"]] - start, nonzero = NA_integer_,
    survival_active = NA_integer_, C = NA_real_, error = conditionMessage(e))))
  answer$warnings <- warnings
  saveRDS(answer, path)
  print(answer$row)
  invisible(NULL)
}

#' Launch one child with operating-system peak memory measurement
#' @param modality,count,arm,fold Registered job identifiers.
#' @return One checkpoint's metrics, including peak resident memory in MiB.
#' @examples
#' # c_launch("methylation", 10000, "PL", 0)
c_launch <- function(modality, count, arm, fold) {
  key <- sprintf("%s_%d_%s_f%d", modality, count, arm, fold)
  path <- file.path(C_LOCAL, paste0(key, ".rds"))
  log <- file.path(C_LOCAL, paste0(key, ".log"))
  if (!file.exists(path)) {
    # macOS time -l reports maximum resident set size in bytes.
    launch_start <- proc.time()[["elapsed"]]
    status <- system2("/usr/bin/time", c("-l", shQuote(file.path(R.home("bin"), "Rscript")),
      shQuote(file.path(C_LOCAL, "source/runner.R")), "one", modality, count, arm, fold,
      if (!is.null(C_K)) paste0("--K=", C_K)), stdout = log, stderr = log)
    if (status != 0L || !file.exists(path)) {
      error <- paste("Child exit", status, "for", key, "see private log", log)
      message(error)
      # Preserve an existing fit if only resource measurement failed. Otherwise
      # record an explicit failed row so subsequent registered counts still run.
      if (!file.exists(path)) saveRDS(list(row = data.frame(modality, count, arm, fold,
        n_features = count, status = "failed", converged = FALSE, iterations = NA_integer_,
        seconds = proc.time()[["elapsed"]] - launch_start, nonzero = NA_integer_,
        survival_active = NA_integer_, C = NA_real_, error = error), warnings = character()), path)
      else stop(error)
    }
  }
  a <- readRDS(path)
  memory <- if (file.exists(log)) grep("maximum resident set size", readLines(log), value = TRUE) else character()
  a$row$peak_rss_MiB <- if (length(memory)) as.numeric(trimws(sub("maximum.*", "", memory[1]))) / 1024^2 else NA_real_
  lines <- if (file.exists(log)) readLines(log) else character()
  wall <- grep("^[[:space:]]*[0-9.]+[[:space:]]+real", lines, value = TRUE)
  a$row$process_seconds <- if (length(wall)) as.numeric(strsplit(trimws(wall[1]), "[[:space:]]+")[[1]][1]) else NA_real_
  a$row$memory_measured <- is.finite(a$row$peak_rss_MiB)
  a$row$warning_count <- length(a$warnings)
  config_path <- file.path(C_LOCAL, "source/config/globals.yml")
  a$row$K <- if (!is.null(a$fit)) ncol(a$fit$EF) else if (!is.null(C_K)) C_K else if (file.exists(config_path))
    yaml::read_yaml(config_path)$benchmark$k_merged_yfb_desurv else NA_integer_
  a$row
}

#' Run full fits in increasing size, then only eligible five-fold comparisons
#' @return NULL invisibly, with aggregate metrics and CV eligibility table.
#' @examples
#' # c_run()
c_run <- function() {
  z <- readRDS(file.path(C_LOCAL, "design.rds"))
  counts <- z$counts
  rm(z); gc()
  full <- list(); eligibility <- list()
  for (n in counts) {
    rows <- lapply(c("PL", "PE"), function(a) c_launch("methylation", n, a, 0L))
    full <- c(full, rows)
    projected <- 5 * sum(vapply(rows, function(r) if (is.finite(r$process_seconds)) r$process_seconds else r$seconds, 0.0))
    eligibility[[length(eligibility) + 1L]] <- data.frame(count = n,
      projected_both_setups_CV_seconds = projected,
      PL_converged = rows[[1]]$converged, PE_converged = rows[[2]]$converged,
      timing_eligible = projected <= 3 * 3600)
    write.csv(do.call(rbind, full), file.path(C_OUT, "full_fits.csv"), row.names = FALSE)
    write.csv(do.call(rbind, eligibility), file.path(C_OUT, "cv_eligibility.csv"), row.names = FALSE)
    cat("COUNT", n, "projected five-fold both setups:", projected / 3600, "hours\n")
  }
  for (a in c("PL", "PE")) full[[length(full) + 1L]] <- c_launch("expression", 2218L, a, 0L)
  write.csv(do.call(rbind, full), file.path(C_OUT, "full_fits.csv"), row.names = FALSE)
  tasks <- expand.grid(modality = "expression", count = 2218L, arm = c("PL", "PE"),
                       fold = 1:5, stringsAsFactors = FALSE)
  for (e in eligibility) if (e$timing_eligible) for (a in c("PL", "PE")) {
    if (e[[paste0(a, "_converged")]]) tasks <- rbind(tasks,
      data.frame(modality = "methylation", count = e$count, arm = a, fold = 1:5))
  }
  # Forking only orchestrates fresh Rscript children. No fitter/BLAS is loaded
  # in this parent, avoiding the historical macOS fork-after-Rcpp failure.
  small <- split(tasks[tasks$count < 100000, ], seq_len(sum(tasks$count < 100000)))
  run_task <- function(t) c_launch(t$modality, t$count, t$arm, t$fold)
  cv <- parallel::mclapply(small, run_task, mc.cores = 3L, mc.preschedule = FALSE)
  if (any(!vapply(cv, function(x) is.data.frame(x) && nrow(x) == 1L, logical(1))))
    stop("Missing/failed CV worker result; inspect private logs.")
  large <- tasks[tasks$count >= 100000, ]
  for (i in seq_len(nrow(large))) cv[[length(cv) + 1L]] <- run_task(large[i, ])
  cv <- do.call(rbind, cv)
  planned <- with(tasks, paste(modality, count, arm, fold))
  observed <- with(cv, paste(modality, count, arm, fold))
  if (anyDuplicated(observed) || !setequal(planned, observed) || nrow(cv) != nrow(tasks))
    stop("CV result coverage differs from the registered task list.")
  write.csv(cv, file.path(C_OUT, "cv_fits.csv"), row.names = FALSE)
  invisible(NULL)
}

# Descriptive summaries and fixed external scoring ----

#' Read one private fit checkpoint
#' @param modality,count,arm,fold Registered job identifiers.
#' @return Saved fit result or NULL when CV was deferred.
#' @examples
#' # c_read("methylation", 10000, "PL", 0)
c_read <- function(modality, count, arm, fold) {
  path <- file.path(C_LOCAL, sprintf("%s_%d_%s_f%d.rds", modality, count, arm, fold))
  if (file.exists(path)) readRDS(path) else NULL
}

#' Describe full programs, fold stability and paired training CV
#' No tuning or selection uses these summaries or external outcomes.
#' @return NULL invisibly, writing aggregate CV and program tables.
#' @examples
#' # c_summary()
c_summary <- function() {
  z <- readRDS(file.path(C_LOCAL, "design.rds"))
  cv <- read.csv(file.path(C_OUT, "cv_fits.csv"))
  eligibility <- read.csv(file.path(C_OUT, "cv_eligibility.csv"))
  K <- if (is.null(C_K)) yaml::read_yaml(file.path(C_LOCAL, "source/config/globals.yml"))$benchmark$k_merged_yfb_desurv else C_K
  for (name in c("full_fits.csv", "cv_fits.csv", "smoke.csv")) {
    path <- file.path(C_OUT, name)
    if (file.exists(path)) {
      metrics <- read.csv(path)
      if ("K" %in% names(metrics) && any(metrics$K != K, na.rm = TRUE))
        stop("Saved metric rank differs from the frozen model configuration.")
      metrics$K <- K
      write.csv(metrics, path, row.names = FALSE)
    }
  }
  programs <- matches <- summaries <- list()
  for (a in c("PL", "PE")) {
    ex <- c_read("expression", 2218L, a, 0L)
    ex_cv <- cv[cv$modality == "expression" & cv$arm == a, ]
    for (n in z$counts) {
      m <- c_read("methylation", n, a, 0L)
      rows <- cv[cv$modality == "methylation" & cv$count == n & cv$arm == a, ]
      c_values <- rows$C[match(1:5, rows$fold)]
      paired <- c_paired(c_values, ex_cv$C[match(1:5, ex_cv$fold)])
      converged_folds <- sum(rows$converged)
      timing_ok <- eligibility$timing_eligible[match(n, eligibility$count)]
      cv_status <- if (nrow(rows) == 5) {
        if (converged_folds == 5) "complete" else "complete_with_nonconverged_folds"
      } else if (!m$row$converged) "not_run_full_fit_nonconverged"
      else if (!timing_ok) "deferred_over_3_hours" else "incomplete"
      summaries[[length(summaries) + 1L]] <- data.frame(K = K, count = n, arm = a,
        cv_status = cv_status, completed_folds = nrow(rows), converged_folds = converged_folds,
        defined_C_folds = sum(is.finite(c_values)),
        mean_C = if (all(is.finite(c_values))) mean(c_values) else NA_real_,
        expression_mean_C = if (nrow(ex_cv) == 5 && all(is.finite(ex_cv$C))) mean(ex_cv$C) else NA_real_,
        mean_difference = paired[1], paired_SE = paired[2],
        meets_convergence_criterion = m$row$converged && converged_folds >= 4)
      if (is.null(m$fit) || !m$row$converged) next
      k <- ncol(m$scores)
      stable <- integer(k)
      for (f in 1:5) {
        mf <- c_read("methylation", n, a, f)
        if (is.null(mf) || is.null(mf$fit)) next
        mapping <- c_match(m$scores[mf$tr, , drop = FALSE], mf$scores[mf$tr, , drop = FALSE])
        r <- suppressWarnings(stats::cor(m$scores[mf$tr, , drop = FALSE],
                                        mf$scores[mf$tr, , drop = FALSE]))
        matched_r <- abs(r[cbind(seq_len(k), mapping)])
        # A constant fold score cannot substantiate a corresponding program.
        active <- mf$activity$active[mapping] & is.finite(matched_r) & matched_r > 0 & mf$row$converged
        stable <- stable + active
        matches[[length(matches) + 1L]] <- data.frame(K = K, count = n, arm = a, fold = f,
          full_program = seq_len(k), fold_program = mapping, abs_r = matched_r,
          fold_converged = mf$row$converged, matched_active = active)
      }
      r <- if (!is.null(ex$fit) && ex$row$converged)
        abs(suppressWarnings(stats::cor(m$scores, ex$scores))) else matrix(NA_real_, k, k)
      max_r <- apply(r, 1, function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_)
      assessable <- nrow(rows) == 5 && converged_folds >= 4 && ex$row$converged
      programs[[length(programs) + 1L]] <- data.frame(K = K, count = n, arm = a,
        program = seq_len(k), nonzero_loadings = m$activity$alive,
        beta_mean = m$fit$EBeta,
        beta_sd = sqrt(pmax(m$fit$EBeta2 - m$fit$EBeta^2, 0)),
        full_survival_active = m$activity$active, active_converged_folds = stable,
        max_abs_expression_r = max_r, methylation_only_assessable = assessable,
        methylation_only = if (assessable) m$activity$active & stable >= 3 & max_r < .5 else NA)
    }
  }
  write.csv(do.call(rbind, summaries), file.path(C_OUT, "cv_summary.csv"), row.names = FALSE)
  if (length(programs)) write.csv(do.call(rbind, programs), file.path(C_OUT, "programs.csv"), row.names = FALSE)
  if (length(matches)) write.csv(do.call(rbind, matches), file.path(C_OUT, "fold_program_matching.csv"), row.names = FALSE)
  invisible(NULL)
}

#' Score saved fixed full fits on matched ICGC donors for description only
#' Uses the existing donor-level clinical reconciliation and whole-universe KNN
#' helper. Features, scaling, fit and Cox score direction are frozen on TCGA.
#' @return NULL invisibly. Private external matrices/scores and aggregate C table.
#' @examples
#' # c_external()
c_external <- function() {
  c_load()
  previous <- getwd()
  setwd(normalizePath(file.path(C_LOCAL, "source")))
  tryCatch(source("code/load_multiomics_data.R"), finally = setwd(previous))
  cfg <- yaml::read_yaml(file.path(C_LOCAL, "source/config/globals.yml"))
  cache <- file.path(C_LOCAL, "external.rds")
  if (!file.exists(cache)) {
    u <- readRDS(file.path(C_LOCAL, "design.rds"))
    # Read every prepared CpG before subsetting to preserve the loader's KNN
    # imputation neighborhood. Clinical rows are collapsed to one per donor.
    ext <- read_icgc_multiomics(cfg$multiomics_data$root, colnames(u$methylation),
      cfg$multiomics_data$icgc_primary_histology, cfg$multiomics_data$icgc_primary_sample_type)
    saveRDS(ext, cache)
    rm(u); gc()
  } else ext <- readRDS(cache)
  full <- read.csv(file.path(C_OUT, "full_fits.csv"))
  rows <- list()
  for (i in seq_len(nrow(full))) {
    task <- full[i, ]
    m <- c_read(task$modality, task$count, task$arm, 0L)
    raw <- if (task$modality == "expression") t(ext$expression) else t(ext$methylation)
    ids <- match(m$features, colnames(raw))
    if (anyNA(ids)) stop("External cohort is missing frozen features.")
    risk <- NULL
    score_error <- if (is.null(m$fit)) "No valid saved full fit" else NA_character_
    if (!is.null(m$fit)) {
      pred <- predict_ssbmf_yfb(c_transform(raw[, ids, drop = FALSE], m$scaling),
        m$fit$EF, m$fit$EBeta, m$fit$EF_norms, projection_center = m$fit$projection_center)
      if (any(!is.finite(pred$risk_scores))) stop("Nonfinite fixed external scores.")
      risk <- pred$risk_scores
      saveRDS(list(risk = risk, ids = rownames(raw)),
        file.path(C_LOCAL, sprintf("external_scores_%s_%d_%s.rds", task$modality, task$count, task$arm)))
    }
    for (population in c("primary_PDAC", "all_matched_sensitivity")) {
      keep <- if (population == "primary_PDAC") which(ext$clinical$primary_pdac) else seq_along(ext$time)
      rows[[length(rows) + 1L]] <- data.frame(K = if (is.null(C_K)) cfg$benchmark$k_merged_yfb_desurv else C_K, modality = task$modality,
        count = task$count, arm = task$arm, population = population,
        n = length(keep), events = sum(ext$event[keep]), full_status = task$status,
        full_converged = task$converged, full_iterations = task$iterations,
        score_error = score_error,
        C_external_description_only = if (!is.null(risk))
          frozen_reverse_cindex(risk[keep], ext$time[keep], ext$event[keep]) else NA_real_)
    }
  }
  write.csv(do.call(rbind, rows), file.path(C_OUT, "external_description.csv"), row.names = FALSE)
  invisible(NULL)
}

#' Complete the matched expression reference before the methylation sweep
#' @return NULL invisibly; full/CV fit rows and aggregate expression summaries.
#' @examples
#' # c_set_rank(5); c_rerun()
c_rerun <- function() {
  stopifnot(!is.null(C_K), C_LOCAL != C_BASE_LOCAL)
  c_freeze()
  full <- lapply(c("PL", "PE"), function(a) c_launch("expression", 2218L, a, 0L))
  tasks <- expand.grid(modality = "expression", count = 2218L, arm = c("PL", "PE"),
                       fold = 1:5, stringsAsFactors = FALSE)
  write.csv(do.call(rbind, full), file.path(C_OUT, "full_fits.csv"), row.names = FALSE)
  cv <- parallel::mclapply(split(tasks, seq_len(nrow(tasks))), function(t)
    c_launch(t$modality, t$count, t$arm, t$fold), mc.cores = 3L, mc.preschedule = FALSE)
  if (any(!vapply(cv, function(x) is.data.frame(x) && nrow(x) == 1L, logical(1))))
    stop("Missing/failed CV worker result; inspect private logs.")
  cv <- do.call(rbind, cv)
  keys <- with(cv, paste(modality, count, arm, fold))
  if (anyDuplicated(keys) || !setequal(keys, with(tasks, paste(modality, count, arm, fold))))
    stop("CV coverage differs from the rerun task list.")
  write.csv(cv, file.path(C_OUT, "cv_fits.csv"), row.names = FALSE)
  summary <- do.call(rbind, lapply(c("PL", "PE"), function(a) {
    rows <- cv[cv$arm == a, ]; f <- full[[match(a, c("PL", "PE"))]]
    data.frame(K = C_K, arm = a, mean_C = if (all(is.finite(rows$C))) mean(rows$C) else NA_real_,
      defined_folds = sum(is.finite(rows$C)), converged_folds = sum(rows$converged),
      mean_nonzero = mean(rows$nonzero), mean_survival_active = mean(rows$survival_active),
      full_status = f$status, full_nonzero = f$nonzero, full_survival_active = f$survival_active)
  }))
  write.csv(summary, file.path(C_OUT, "expression_summary.csv"), row.names = FALSE)
  if (any(c(vapply(full, function(r) r$status == "failed", logical(1)), cv$status == "failed")))
    stop("Rerun contains failed fits; aggregate diagnostics were saved.")
  invisible(NULL)
}

# Main ----
if (sys.nframe() == 0L) {
  dir.create(C_LOCAL, recursive = TRUE, showWarnings = FALSE)
  dir.create(C_OUT, recursive = TRUE, showWarnings = FALSE)
  args <- commandArgs(TRUE)
  rank_arg <- grep("^--K=", args, value = TRUE)
  stopifnot(length(rank_arg) <= 1L)
  if (length(rank_arg)) c_set_rank(as.numeric(sub("^--K=", "", rank_arg)))
  args <- args[!grepl("^--K=", args)]
  mode <- args[1]
  if (mode == "prepare") c_prepare()
  else if (mode == "one") c_one(args[2], as.integer(args[3]), args[4], as.integer(args[5]))
  else if (mode == "smoke") {
    c_freeze()
    rows <- lapply(c("PL", "PE"), function(a) c_launch("methylation", 2000L, a, 1L))
    write.csv(do.call(rbind, rows), file.path(C_OUT, "smoke.csv"), row.names = FALSE)
    if (any(vapply(rows, function(r) r$status == "failed", logical(1)))) stop("Smoke fit failed.")
  } else if (mode == "run") { c_freeze(); c_run(); c_summary(); c_external() }
  else if (mode == "rerun") c_rerun()
  else if (mode == "sweep") {
    # Explicit rank required. The historical stopped study cannot resume here.
    stopifnot(!is.null(C_K), C_LOCAL != C_BASE_LOCAL)
    c_freeze(); c_run(); c_summary()
    rows <- rbind(read.csv(file.path(C_OUT, "full_fits.csv")), read.csv(file.path(C_OUT, "cv_fits.csv")))
    if (any(rows$status == "failed")) stop("Sweep contains failed fits; inspect aggregate diagnostics.")
  }
  else if (mode == "summary") { c_summary(); c_external() }
  else stop("Use prepare, smoke, run, summary, rerun --K=5, or sweep --K=5.")
}
