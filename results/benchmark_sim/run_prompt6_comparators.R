# ============================================================
# Script: run_prompt6_comparators.R
# Purpose: Compare frozen S-11 risks and training-only PDAC comparators.
# Created: 2026-11-06
# Dependencies: survival, yaml, glmnet, DeSurv
# ============================================================

suppressPackageStartupMessages(library(survival))
source("results/benchmark_sim/run_prompt6_desurv.R")
source("code/preprocess_desurv.R")
source("code/predict_ssbmf_yfb.R")
source("code/concordance_ci.R")

#' Load the fixed S-11 models and their training design.
#' @param fold 0 for the full fit, or the S-11 validation fold number.
#' @return List with Y, outcomes, fold labels, and the two fixed model states.
#' @examples
#' # prompt6_reference()
prompt6_reference <- function(fold = 0L) {
  d <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
  fits <- lapply(1:0, function(w) readRDS(file.path(
    "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/fits",
    sprintf("size_D4_N2218_w%d_f%d.rds", w, fold))))
  stopifnot(identical(fits[[1]]$genes, fits[[2]]$genes))
  genes <- fits[[1]]$genes
  stopifnot(all(genes %in% d$genes))
  d$Y <- d$Y[, match(genes, d$genes), drop = FALSE]
  d$genes <- genes
  d$w1 <- fits[[1]]$state
  d$w0 <- fits[[2]]$state
  live <- which(colSums(d$w0$EF^2) > 1e-8)
  train <- if (fold == 0L) rep(TRUE, nrow(d$Y)) else d$fold != fold
  score <- d$Y[train, , drop = FALSE] %*% d$w0$EF[, live, drop = FALSE]
  d$cox0 <- coxph(Surv(d$time[train], d$status[train]) ~ score)
  stopifnot(all(is.finite(coef(d$cox0))))
  d$b0 <- numeric(ncol(d$w0$EF))
  d$b0[live] <- coef(d$cox0)
  d
}

#' Score the frozen S-12 reference arms.
#' @param Y A cohort matrix aligned to the D4 genes.
#' @param d Output of prompt6_reference().
#' @return List of Cox-oriented risk vectors.
#' @examples
#' # prompt6_risk(d$Y, d)
prompt6_risk <- function(Y, d) {
  list(w1 = as.numeric(predict_ssbmf_yfb(Y, d$w1$EF, d$w1$EBeta,
    EF_norms = d$w1$EF_norms)$risk_scores),
    twostep = drop(Y %*% d$w0$EF %*% d$b0))
}

#' Link eligible sample IDs to subtype calls and donor IDs.
#' @param co Cohort name. @param raw Output of load_pdac_raw().
#' @return Data frame aligned to raw rows, with missing calls retained.
#' @examples
#' # prompt6_metadata("CPTAC", raw)
prompt6_metadata <- function(co, raw) {
  root <- file.path(PDAC_DATA_ROOT, "original")
  ids <- as.character(raw$sampID)
  stopifnot(length(ids) == raw$n, !anyNA(ids), all(nzchar(ids)), !anyDuplicated(ids))
  source_ids <- if (co == "PACA_AU_seq") sub("_seq$", "", ids) else ids
  subtype <- readRDS(file.path(root, paste0(co, ".caf_subtype.rds")))$Subtype
  calls <- read.csv(file.path(root, paste0(co, "_subtype.csv")))
  stopifnot(!anyDuplicated(subtype$sampID), !anyDuplicated(calls$ID))
  donor <- paste0(co, ":", ids)
  if (co %in% c("PACA_AU_array", "PACA_AU_seq")) {
    meta <- readRDS(file.path(root, paste0(co, ".rds")))$sampInfo
    stopifnot(!anyDuplicated(meta$icgc_sample_id))
    donor_ids <- as.character(meta$icgc_donor_id[match(source_ids, meta$icgc_sample_id)])
    stopifnot(!anyNA(donor_ids), all(nzchar(donor_ids)))
    donor <- paste0("PACA_AU:", donor_ids)
  }
  j <- match(source_ids, subtype$sampID)
  data.frame(prob = subtype$PurIST.prob[j], PurIST = subtype$PurIST[j],
             DeCAF = calls$DeCAF[match(source_ids, calls$ID)], donor = donor)
}

#' Apply S-12 preprocessing and zero-fill missing frozen training genes.
#' @param co External cohort name. @param genes Frozen D4 genes in model order.
#' @return List with prepared matrix, outcome vectors, IDs and subtype calls.
#' @examples
#' # prompt6_external("Dijk", d$genes)
prompt6_external <- function(co, genes) {
  raw <- load_pdac_raw(co, PDAC_DATA_ROOT)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
    log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
    rank_transform = FALSE, per_platform_standardize = TRUE)
  Y <- matrix(0, raw$n, length(genes))
  j <- match(genes, pre$gene_names)
  Y[, !is.na(j)] <- pre$Y[, j[!is.na(j)], drop = FALSE]
  list(Y = Y, time = raw$time, status = raw$status, ids = raw$sampID,
       meta = prompt6_metadata(co, raw))
}

#' Compute external C, calibration, and paired and pooled comparisons.
#' @param external Named cohort list with time, status, meta and risks.
#' @return Named aggregate data frames for model performance and coverage.
#' @examples
#' # prompt6_evaluate(external)
prompt6_evaluate <- function(external) {
  models <- names(external[[1]]$risks)
  by_cohort <- list(); pairs <- list(); hr <- list()
  for (model in models) {
    paired <- list(); pooled <- list()
    for (co in names(external)) {
      x <- external[[co]]
      all <- data.frame(time = x$time, event = x$status, risk = x$risks[[model]],
                        reference = x$risks$w1, donor = x$meta$donor, cohort = co)
      ok <- is.finite(all$risk)
      a <- all[ok, ]
      row <- data.frame(cohort = co, model = model, eligible = nrow(all), n = nrow(a),
                        dropped = sum(!ok), events = sum(a$event), C = NA_real_,
                        lo = NA_real_, hi = NA_real_, cal_slope = NA_real_, cal_se = NA_real_)
      if (nrow(a) >= 10L && sum(a$event) >= 2L && sd(a$risk) > 0) {
        row$C <- frozen_reverse_cindex(a$risk, a$time, a$event)
        set.seed(1)
        boot <- replicate(2000, {
          j <- sample.int(nrow(a), replace = TRUE)
          frozen_reverse_cindex(a$risk[j], a$time[j], a$event[j])
        })
        stopifnot(all(is.finite(boot)))
        row[1, c("lo", "hi")] <- quantile(boot, c(.025, .975))
        cal <- coxph(Surv(time, event) ~ risk, data = a)
        row$cal_slope <- unname(coef(cal))
        row$cal_se <- sqrt(vcov(cal)[1, 1])
        pooled[[co]] <- a
        if (model != "w1" && all(is.finite(a$reference))) paired[[co]] <- a
      }
      by_cohort[[length(by_cohort) + 1L]] <- row
    }
    if (length(paired)) {
      col <- function(nm) lapply(paired, `[[`, nm)
      ci <- bootstrap_concordance_diff_ci_stratified(col("reference"), col("risk"),
        col("time"), col("event"), donor_list = col("donor"), frozen = TRUE)
      pairs[[model]] <- data.frame(model = model, strata = paste(names(paired), collapse = "|"),
        n = sum(lengths(col("risk"))), donors = length(unique(unlist(col("donor")))),
        diff = ci$estimate, lo = ci$lower, hi = ci$upper, B = ci$B)
    }
    if (length(pooled)) {
      a <- do.call(rbind, pooled)
      a$risk_z <- a$risk / sd(a$risk)
      m <- summary(coxph(Surv(time, event) ~ risk_z + strata(cohort), data = a))
      # Review addition: robust (sandwich) SE clustered on donor, because the
      # 47 PACA-AU donors measured on both assays enter the pooled fit twice.
      rc <- coxph(Surv(time, event) ~ risk_z + strata(cohort) + cluster(donor), data = a)
      rse <- sqrt(vcov(rc)[1, 1])
      hr[[model]] <- data.frame(model = model, n = nrow(a), events = sum(a$event),
        hr = m$conf.int[1, 1], lo = m$conf.int[1, 3], hi = m$conf.int[1, 4],
        p = m$coefficients[1, 5], donors = length(unique(a$donor)),
        donor_lo = exp(coef(rc) - qnorm(.975) * rse), donor_hi = exp(coef(rc) + qnorm(.975) * rse))
    }
  }
  list(external = do.call(rbind, by_cohort), paired = do.call(rbind, pairs),
       pooled_hr = do.call(rbind, hr))
}

if (sys.nframe() == 0L) {
  options(warn = 1)
  d <- prompt6_reference()
  train <- lapply(cfg$pdac$training_cohorts, load_pdac_raw, pdac_root = PDAC_DATA_ROOT)
  stopifnot(identical(as.numeric(unlist(lapply(train, `[[`, "time"))), as.numeric(d$time)),
            identical(as.numeric(unlist(lapply(train, `[[`, "status"))), as.numeric(d$status)))
  prob <- unlist(Map(function(co, raw) prompt6_metadata(co, raw)$prob,
                     cfg$pdac$training_cohorts, train))
  purist <- coxph(Surv(d$time, d$status) ~ prob)
  stopifnot(is.finite(coef(purist)))
  alpha <- seq(0, 1, by = .25)
  cv_warnings <- rep("", length(alpha))
  cv <- lapply(seq_along(alpha), function(i) withCallingHandlers(
    glmnet::cv.glmnet(d$Y, Surv(d$time, d$status), alpha = alpha[i], family = "cox",
      foldid = d$fold, type.measure = "deviance", cox.ties = "breslow"),
    warning = function(w) {
      cv_warnings[i] <<- paste(cv_warnings[i], conditionMessage(w), sep = " | ")
    }))
  tune <- data.frame(alpha = alpha, lambda = vapply(cv, `[[`, numeric(1), "lambda.min"),
    cv_deviance = vapply(cv, function(f) min(f$cvm, na.rm = TRUE), numeric(1)),
    path_warning = trimws(cv_warnings))
  best <- order(tune$cv_deviance, -alpha)[1]
  tune$selected <- seq_along(alpha) == best
  work <- "data/prompt6_desurv"
  final_path <- file.path(work, "final.rds")
  partial <- "--without-desurv" %in% commandArgs(TRUE)
  if (!partial) desurv <- readRDS(final_path)
  external <- setNames(lapply(cfg$pdac$external_cohorts, function(co) {
    x <- prompt6_external(co, d$genes)
    x$risks <- c(prompt6_risk(x$Y, d), list(
      PurIST = as.numeric(coef(purist)) * x$meta$prob,
      elastic_net = as.numeric(predict(cv[[best]], x$Y, s = "lambda.min", type = "link")),
      clinical = rep(NA_real_, length(x$time))))
    if (!partial) {
      a <- desurv_author_prepare(setNames(list(prompt6_author_data(co)), co),
                                 genes = rownames(desurv$W))
      stopifnot(identical(as.character(colnames(a$ex)), as.character(x$ids)),
                identical(as.numeric(a$sampInfo$time), as.numeric(x$time)),
                identical(as.numeric(a$sampInfo$event), as.numeric(x$status)))
      x$risks$DeSurv <- desurv_author_risk(desurv, a$ex, 270L)
    }
    x
  }), cfg$pdac$external_cohorts)
  stopifnot(length(intersect(external$PACA_AU_array$meta$donor,
                             external$PACA_AU_seq$meta$donor)) == 47L)
  result <- prompt6_evaluate(external)
  prior <- read.csv(paste0("results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/",
                           "external_cindex.csv"))
  actual <- result$external[result$external$model %in% c("w1", "twostep"), ]
  check <- merge(actual, prior, by = c("cohort", "model"))
  stopifnot(nrow(check) == 10L, all(check$n.x == check$n.y),
            max(abs(check$C.x - check$C.y)) < 1e-10)
  clinical <- data.frame(cohort = c(cfg$pdac$training_cohorts, cfg$pdac$external_cohorts),
    stage = c("AJCC.pathologic.tumor.stage", "tumor_stage_pathological", NA,
              "T; N; M", "AJCC.Pathology.Stage", "AJCC.Pathology.Stage", NA),
    grade = c("Grade", NA, "differentiation grade", "Diff", "Tumour.Grade", "Tumour.Grade",
              "Characteristics.TNM.tumour.grading."),
    age = c("Age.at.initial.pathologic.diagnosis", "age", "age (year)", NA,
            "Age.at.Diagnosis.in.Years", "Age.at.Diagnosis.in.Years", NA),
    margin = c("Residual.tumor", "residual_tumor", NA, "Margin", NA, NA,
               "Characteristics.resection.margin."))
  for (i in seq_len(nrow(clinical))) {
    fields <- unlist(clinical[i, -1], use.names = FALSE)
    fields <- unlist(strsplit(fields[!is.na(fields)], "; ", fixed = TRUE))
    meta <- readRDS(file.path(PDAC_DATA_ROOT, "original",
                             paste0(clinical$cohort[i], ".rds")))$sampInfo
    stopifnot(all(fields %in% names(meta)))
  }
  stopifnot(all(colSums(is.na(clinical[, -1])) > 0))
  result$clinical_fields <- clinical
  result$tuning <- tune
  result$training <- data.frame(model = c("w1", "twostep", "PurIST", "elastic_net", "clinical"),
    eligible = nrow(d$Y), n = c(nrow(d$Y), nrow(d$Y), purist$n, nrow(d$Y), 0L),
    n_cv = c(0L, 0L, 0L, length(alpha) * length(unique(d$fold)), 0L),
    final_new = c(0L, 1L, 1L, length(alpha), 0L))
  if (!partial) {
    starts <- do.call(rbind, lapply(Sys.glob(file.path(work, "seed_*_status.csv")), read.csv))
    result$desurv_starts <- starts
    result$desurv_final <- data.frame(n = nrow(d$Y), converged = desurv$convergence,
      iterations = length(desurv$lossit), nan_flag = desurv$nan_flag)
  }
  out <- if (partial) "data/prompt6_partial" else "results/benchmark_sim/outputs/prompt6_comparators"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(result))
    write.csv(result[[nm]], file.path(out, paste0(nm, ".csv")), row.names = FALSE)
}
