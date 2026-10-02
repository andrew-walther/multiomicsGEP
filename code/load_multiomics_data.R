# ============================================================
# Script: load_multiomics_data.R
# Purpose: Build matched expression + methylation + survival cohorts (TCGA
#          training, ICGC validation) in the multimodal YFB data contract.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: impute (Bioconductor), yaml; code/preprocess_desurv.R
# ============================================================
#
# Output contract (matches code/preprocess_multimodal_yfb.R):
#   Y     = list(expression = n x p_expr, methylation = n x p_meth), rows are
#           subjects, columns are features, finite and nonnegative
#   time  = named positive follow-up times (months), names = subject IDs
#   event = named 0/1 event indicators (1 = death)
#
# Source files (see docs/plans/Working_Plan_10_1_26.md, Workstream 1):
#   TCGA expression  tcga/expr_log_normalized.rds        log2(TMM-CPM + 1), genes x 150
#   TCGA methylation tcga/tcga_meth_filtered_sex.rds$meth beta values, CpGs x 150
#   TCGA survival    PDAC_data/original/TCGA_PAAD.survival_data.rds
#   ICGC expression  icgc/icgc_matched_data.rds$expr      unlogged, genes x 79 donors
#   ICGC methylation icgc/icgc_matched_data.rds$meth      beta values, CpGs x 79 donors
#   ICGC survival    icgc/icgc_matched_data.rds$info.expr one row per RNA sample
#   ICGC CpG filter  icgc/icgc_meth_filtered_sex.rds      rownames only
#
# Memory: the TCGA methylation object is ~0.8 GB in RAM and the ICGC matched
# object ~0.6 GB. Building both cohorts peaks at a few GB; fine on a 16 GB
# laptop. The full 305k-CpG fit (no screening) is a Longleaf job.

if (!exists("select_top_variable_genes")) source("code/preprocess_desurv.R")

# Survival and clinical parsing ----

#' Keep subjects with usable survival and report the ones dropped
#'
#' Fails loudly: every dropped subject is listed by reason in the returned
#' `dropped` table and summarized in a message.
#'
#' @param ids Character subject IDs.
#' @param time Numeric follow-up times (months), aligned with `ids`.
#' @param event Numeric 0/1 event indicators, aligned with `ids`.
#' @param label Cohort label used in the message.
#' @return List with `keep` (logical), `time`, `event` (named, kept subjects
#'   only), and `dropped` (data.frame of id and reason).
#' @examples
#' filter_usable_survival(c("a", "b"), c(5, NA), c(1, 0), "toy")
#' @family multiomics_data
filter_usable_survival <- function(ids, time, event, label) {
  reason <- rep(NA_character_, length(ids))
  reason[is.na(time)] <- "missing time"
  reason[is.na(reason) & time <= 0] <- "time <= 0"
  reason[is.na(reason) & !(event %in% c(0, 1))] <- "missing or invalid event"
  keep <- is.na(reason)
  if (any(!keep)) {
    message(label, ": dropped ", sum(!keep), " of ", length(ids),
            " subjects without usable survival (",
            paste(names(table(reason)), table(reason), sep = " = ", collapse = "; "),
            ").")
  }
  list(
    keep = keep,
    time = stats::setNames(as.numeric(time[keep]), ids[keep]),
    event = stats::setNames(as.numeric(event[keep]), ids[keep]),
    dropped = data.frame(id = ids[!keep], reason = reason[!keep],
                         stringsAsFactors = FALSE)
  )
}

#' Collapse ICGC per-RNA-sample clinical rows to one row per donor
#'
#' `info.expr` has one row per RNA sample; one donor (DO33168) has two. The
#' survival fields must agree across a donor's rows, otherwise this stops.
#' A donor counts as primary PDAC only if every one of its samples is a
#' primary tumour with PDAC histology.
#'
#' @param info `info.expr` data.frame from `icgc_matched_data.rds`.
#' @param histology HistoSubtype value that defines PDAC.
#' @param sample_type Sample.type value that defines a primary tumour.
#' @return data.frame with one row per donor: icgc_donor_id,
#'   submitted_donor_id, time, event, primary_pdac, HistoSubtype, Sample.type.
#' @examples
#' info <- data.frame(icgc_donor_id = c("D1", "D1"), submitted_donor_id = "S1",
#'   survival_months = 10, censored = "death",
#'   HistoSubtype = "Pancreatic Ductal Adenocarcinoma",
#'   Sample.type = "Primary tumour")
#' collapse_icgc_clinical(info, "Pancreatic Ductal Adenocarcinoma", "Primary tumour")
#' @family multiomics_data
collapse_icgc_clinical <- function(info, histology, sample_type) {
  donors <- unique(as.character(info$icgc_donor_id))
  rows <- lapply(donors, function(d) {
    x <- info[as.character(info$icgc_donor_id) == d, , drop = FALSE]
    time <- unique(x$survival_months)
    status <- unique(as.character(x$censored))
    if (length(time) != 1L || length(status) != 1L) {
      stop("ICGC donor ", d, " has conflicting survival across its RNA samples.")
    }
    data.frame(
      icgc_donor_id = d,
      submitted_donor_id = as.character(x$submitted_donor_id[1]),
      time = as.numeric(time),
      # censored is "death" (event), "censor" (censored) or NA
      event = ifelse(is.na(status), NA_real_, as.numeric(status == "death")),
      primary_pdac = all(!is.na(x$HistoSubtype) & x$HistoSubtype == histology &
                           !is.na(x$Sample.type) & x$Sample.type == sample_type),
      HistoSubtype = paste(unique(as.character(x$HistoSubtype)), collapse = "; "),
      Sample.type = paste(unique(as.character(x$Sample.type)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Methylation imputation ----

#' KNN-impute missing methylation beta values (CpGs x subjects)
#'
#' Follows the code Yusha provided (9/18): only CpGs with at least one missing
#' value are passed to `impute::impute.knn()`, which fills each missing entry
#' from the k = 10 nearest CpGs (Euclidean distance over observed subjects).
#' Imputation is done separately within each cohort, so no information moves
#' between training and validation. `impute.knn` uses a fixed internal seed,
#' so results are reproducible.
#'
#' @param meth Numeric matrix, CpGs in rows and subjects in columns.
#' @param label Cohort label used in the message.
#' @return The matrix with missing values filled.
#' @examples
#' m <- matrix(runif(200), 20, 10); m[3, 2] <- NA
#' impute_methylation_knn(m, "toy")
#' @family multiomics_data
impute_methylation_knn <- function(meth, label) {
  idx <- which(rowSums(is.na(meth)) > 0)
  if (length(idx) == 0L) return(meth)
  message(label, ": imputing ", sum(is.na(meth)), " missing beta values (",
          signif(100 * mean(is.na(meth)), 2), "%) in ", length(idx), " CpGs by KNN.")
  # impute.knn prints its cluster-splitting progress with cat(); silence it
  utils::capture.output(
    imputed <- suppressWarnings(impute::impute.knn(meth[idx, , drop = FALSE]))$data
  )
  meth[idx, ] <- imputed
  if (anyNA(meth)) stop(label, ": missing methylation values remain after imputation.")
  meth
}

# Cohort readers ----

#' Read the TCGA-PAAD matched expression, methylation and survival data
#'
#' @param root Directory holding `tcga/` (globals: multiomics_data$root).
#' @param pdac_root Directory holding `original/TCGA_PAAD.survival_data.rds`.
#' @param cpgs CpG IDs to keep (rows of the methylation matrix).
#' @return List with `expression` (genes x subjects), `methylation`
#'   (CpGs x subjects, imputed), `time`, `event` (named), and `dropped`.
#' @family multiomics_data
read_tcga_multiomics <- function(root, pdac_root, cpgs) {
  expr <- readRDS(file.path(root, "tcga", "expr_log_normalized.rds"))
  meth <- readRDS(file.path(root, "tcga", "tcga_meth_filtered_sex.rds"))$meth
  surv <- readRDS(file.path(pdac_root, "original", "TCGA_PAAD.survival_data.rds"))

  ids <- colnames(expr)
  if (!identical(sort(ids), sort(colnames(meth)))) {
    stop("TCGA expression and methylation sample barcodes differ.")
  }
  missing_surv <- setdiff(ids, surv$sampID)
  if (length(missing_surv) > 0L) {
    stop("TCGA samples without a survival record: ", paste(missing_surv, collapse = ", "))
  }
  s <- surv[match(ids, surv$sampID), ]
  usable <- filter_usable_survival(ids, s$time, s$event, "TCGA")
  keep_ids <- ids[usable$keep]

  meth <- meth[cpgs, keep_ids, drop = FALSE]
  list(
    expression = expr[, keep_ids, drop = FALSE],
    methylation = impute_methylation_knn(meth, "TCGA"),
    time = usable$time, event = usable$event, dropped = usable$dropped,
    n_whitelisted = sum(s$whitelist[usable$keep], na.rm = TRUE)
  )
}

#' Read the ICGC (PACA-AU) matched expression, methylation and survival data
#'
#' Expression is unlogged in the source file, so log2(x + 1) is applied here.
#' Expression and methylation columns are both ICGC donor IDs (DO...).
#'
#' @param root Directory holding `icgc/` (globals: multiomics_data$root).
#' @param cpgs CpG IDs to keep.
#' @param histology,sample_type Definition of a primary PDAC donor.
#' @return List with `expression`, `methylation` (features x subjects),
#'   `time`, `event` (named), `clinical` (one row per kept donor, including
#'   `primary_pdac`), and `dropped`.
#' @family multiomics_data
read_icgc_multiomics <- function(root, cpgs, histology, sample_type) {
  x <- readRDS(file.path(root, "icgc", "icgc_matched_data.rds"))
  ids <- colnames(x$expr)
  if (!identical(ids, colnames(x$meth))) {
    stop("ICGC expression and methylation donor IDs are not identically ordered.")
  }
  clinical <- collapse_icgc_clinical(x$info.expr, histology, sample_type)
  if (!setequal(clinical$icgc_donor_id, ids)) {
    stop("ICGC clinical donors do not match the expression/methylation donors.")
  }
  clinical <- clinical[match(ids, clinical$icgc_donor_id), ]
  usable <- filter_usable_survival(ids, clinical$time, clinical$event, "ICGC")
  keep_ids <- ids[usable$keep]

  meth <- x$meth[cpgs, keep_ids, drop = FALSE]
  list(
    expression = log2_plus1_transform(x$expr[, keep_ids, drop = FALSE]),
    methylation = impute_methylation_knn(meth, "ICGC"),
    time = usable$time, event = usable$event,
    clinical = clinical[usable$keep, , drop = FALSE],
    dropped = usable$dropped
  )
}

# Screening and assembly ----

#' Screen features on the training cohort
#'
#' Genes: DeSurv combined mean + variance rank (`select_top_variable_genes(...,
#' method = "combined_rank")`). CpGs: largest variance of the beta values.
#' Zero-variance features are removed first. Both rules see only the training
#' cohort, so the validation cohort plays no part in feature selection.
#'
#' @param train List with `expression` and `methylation`, features x subjects.
#' @param n_genes,n_cpgs Number of genes and CpGs to keep.
#' @return List with character vectors `genes` and `cpgs`.
#' @family multiomics_data
screen_multiomics_features <- function(train, n_genes, n_cpgs) {
  expr <- t(train$expression)
  expr <- expr[, apply(expr, 2, stats::var) > 0, drop = FALSE]
  genes <- select_top_variable_genes(expr, colnames(expr), top_n = n_genes,
                                     method = "combined_rank")$gene_names
  cpg_var <- apply(train$methylation, 1, stats::var)
  cpgs <- names(sort(cpg_var[cpg_var > 0], decreasing = TRUE))[seq_len(min(n_cpgs, sum(cpg_var > 0)))]
  list(genes = genes, cpgs = cpgs)
}

#' Convert a cohort to the multimodal YFB contract on fixed features
#'
#' @param cohort Output of `read_tcga_multiomics()` or `read_icgc_multiomics()`.
#' @param features Output of `screen_multiomics_features()`.
#' @param subjects Optional subset of subject IDs (default: all).
#' @return List with `Y`, `time`, `event` in the contract described at the top.
#' @family multiomics_data
as_multimodal_yfb_cohort <- function(cohort, features, subjects = names(cohort$time)) {
  Y <- list(
    expression = t(cohort$expression[features$genes, subjects, drop = FALSE]),
    methylation = t(cohort$methylation[features$cpgs, subjects, drop = FALSE])
  )
  list(Y = Y, time = cohort$time[subjects], event = cohort$event[subjects])
}

#' Build the TCGA training and ICGC validation cohorts
#'
#' Steps: (1) CpGs = intersection of the TCGA and ICGC sex-chromosome-filtered
#' probe sets; genes = symbols present in both cohorts' expression. (2) Read
#' each cohort, drop subjects without usable survival, KNN-impute methylation
#' within cohort. (3) Screen features on TCGA only. (4) Return TCGA training,
#' ICGC primary-PDAC validation (main) and all-donor ICGC validation
#' (sensitivity), all on the same features.
#'
#' @param root,pdac_root Data directories (default from config/globals.yml).
#' @param n_genes,n_cpgs Screening sizes (default from config/globals.yml).
#' @return List with `training`, `validation_primary`, `validation_all`
#'   (each `Y`, `time`, `event`), `features`, `icgc_clinical`, and `summary`.
#' @examples
#' \dontrun{
#' data <- build_multiomics_cohorts()
#' fit <- fit_multimodal_yfb(data$training$Y, data$training$time,
#'                           data$training$event, K = 7)
#' }
#' @family multiomics_data
build_multiomics_cohorts <- function(root = NULL, pdac_root = NULL,
                                     n_genes = NULL, n_cpgs = NULL) {
  g <- yaml::read_yaml("config/globals.yml")
  cfg <- g$multiomics_data
  root <- root %||% cfg$root
  pdac_root <- pdac_root %||% cfg$pdac_root
  n_genes <- n_genes %||% g$preprocessing$top_n_genes_desurv
  n_cpgs <- n_cpgs %||% cfg$n_cpgs

  # (1) Shared feature universe -- probe names, never positions
  tcga_cpgs <- rownames(readRDS(file.path(root, "tcga", "tcga_meth_filtered_sex.rds"))$meth)
  icgc_cpgs <- rownames(readRDS(file.path(root, "icgc", "icgc_meth_filtered_sex.rds"))$meth)
  cpgs <- intersect(tcga_cpgs, icgc_cpgs)
  message("CpGs in both filtered probe sets: ", length(cpgs))

  # (2) Read cohorts
  tcga <- read_tcga_multiomics(root, pdac_root, cpgs)
  icgc <- read_icgc_multiomics(root, cpgs, cfg$icgc_primary_histology,
                               cfg$icgc_primary_sample_type)
  shared_genes <- intersect(rownames(tcga$expression), rownames(icgc$expression))
  tcga$expression <- tcga$expression[shared_genes, , drop = FALSE]
  message("Genes in both cohorts: ", length(shared_genes))

  # (3) Screen on TCGA only
  features <- screen_multiomics_features(tcga, n_genes, n_cpgs)

  # (4) Assemble
  primary_ids <- icgc$clinical$icgc_donor_id[icgc$clinical$primary_pdac]
  out <- list(
    training = as_multimodal_yfb_cohort(tcga, features),
    validation_primary = as_multimodal_yfb_cohort(icgc, features, primary_ids),
    validation_all = as_multimodal_yfb_cohort(icgc, features),
    features = features,
    icgc_clinical = icgc$clinical,
    dropped = list(tcga = tcga$dropped, icgc = icgc$dropped)
  )
  out$summary <- data.frame(
    cohort = c("TCGA (training)", "ICGC primary PDAC (validation)", "ICGC all (sensitivity)"),
    n = c(length(out$training$time), length(out$validation_primary$time),
          length(out$validation_all$time)),
    events = c(sum(out$training$event), sum(out$validation_primary$event),
               sum(out$validation_all$event)),
    genes = length(features$genes), cpgs = length(features$cpgs)
  )
  out
}

