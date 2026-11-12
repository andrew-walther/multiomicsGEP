# ============================================================
# Script: characterize_methylation_factors_k5.R
# Purpose: Describe saved K = 5 methylation and expression factors.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: IlluminaHumanMethylation450kanno.ilmn12.hg19, fgsea, msigdbr, yaml
# ============================================================

# Fixed descriptive thresholds ----
# Specificity requires expression |r| < 0.5, at least three converged
# matched folds with |r| >= 0.5, and recurrence at another count at |r| >= 0.7.

#' Expand aligned 450K gene and region annotations
#' @param annotation Annotation data frame with row names and UCSC fields.
#' @return One row per CpG and gene-region entry, retaining unannotated CpGs.
#' @examples
#' # expand_annotation(annotation)
expand_annotation <- function(annotation) {
  out <- lapply(seq_len(nrow(annotation)), function(i) {
    #' Split one semicolon-delimited annotation field
    #' @param x A single annotation string.
    #' @return Character vector retaining an empty field.
    #' @examples
    #' # split_field("A;B")
    split_field <- function(x) if (is.na(x) || !nzchar(x)) "" else strsplit(x, ";", fixed = TRUE)[[1]]
    genes <- split_field(annotation$UCSC_RefGene_Name[i])
    regions <- split_field(annotation$UCSC_RefGene_Group[i])
    if (length(genes) != length(regions)) stop("Gene and region annotation lengths differ for ", rownames(annotation)[i])
    data.frame(cpg = rownames(annotation)[i], gene = genes, region = regions,
               island_relation = annotation$Relation_to_Island[i], stringsAsFactors = FALSE)
  })
  unique(do.call(rbind, out))
}

#' Evaluate prespecified factor-specificity conditions
#' @param max_expression Largest absolute expression-factor correlation, or NA.
#' @param stable_folds Number of converged folds with matched |r| at least 0.5.
#' @param recurrence Largest absolute correlation with another CpG count, or NA.
#' @return Logical vector, with NA when the comparison is unavailable.
#' @examples
#' classify_specific(c(0.9, 0.2), c(4, 3), c(0.8, 0.8))
classify_specific <- function(max_expression, stable_folds, recurrence) {
  max_expression < 0.5 & stable_folds >= 3 & recurrence >= 0.7
}

#' Calculate a seeded CpG-overlap permutation distribution
#' @param eligible CpG IDs eligible for sampling.
#' @param hits CpG IDs annotated to expression-factor top genes.
#' @param observed_top CpG IDs in the methylation factor top set.
#' @param seed Random seed.
#' @param nperm Number of random samples.
#' @return Observed count, null mean and upper-tail Monte Carlo p-value.
#' @examples
#' overlap_permutation(letters, letters[1:3], letters[1:5], seed = 42, nperm = 10)
overlap_permutation <- function(eligible, hits, observed_top, seed = 2026L, nperm = 1000L) {
  stopifnot(!anyDuplicated(eligible), all(observed_top %in% eligible), length(observed_top) <= length(eligible))
  hit <- eligible %in% hits
  observed <- sum(observed_top %in% hits)
  set.seed(seed)
  null <- replicate(nperm, sum(sample(hit, length(observed_top), replace = FALSE)))
  c(observed = observed, null_mean = mean(null), null_sd = stats::sd(null),
    p_upper = (1 + sum(null >= observed)) / (nperm + 1))
}

#' Describe CpG location in relation to islands
#' @param x Illumina island relationship labels.
#' @return Island, shore, shelf or open sea.
#' @examples
#' island_class(c("Island", "N_Shore", ""))
island_class <- function(x) {
  x[is.na(x)] <- ""
  ifelse(x == "Island", "island", ifelse(grepl("Shore", x), "shore",
    ifelse(grepl("Shelf", x), "shelf", "open sea")))
}

#' Safely correlate columns, allowing constant score vectors
#' @param x,y Numeric matrices sharing rows.
#' @return Matrix of Pearson correlations; constant columns yield NA.
#' @examples
#' safe_cor(matrix(1:6, 3), matrix(6:1, 3))
safe_cor <- function(x, y) suppressWarnings(stats::cor(x, y))

#' Run the descriptive analysis using saved checkpoints
#' @return Invisibly NULL; writes aggregate K = 5 CSVs only.
#' @examples
#' # characterize_k5()
characterize_k5 <- function() {
  local <- "data/methylation_single_modality_2026-11-12/K05"
  base <- "results/multimodal_real/outputs/methylation_single_modality_2026-11-12/K05"
  out <- file.path(base, "factor_characterization")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(normalizePath(file.path(local, "R-library"), mustWork = TRUE), .libPaths()))
  source("code/pathway_enrichment.R")
  z <- readRDS(file.path(local, "design.rds"))
  stopifnot(nrow(z$expression) == 144L, ncol(z$expression) == 2218L,
    identical(rownames(z$expression), rownames(z$methylation)), length(z$ranked) == 182934L,
    !anyDuplicated(colnames(z$expression)), !anyDuplicated(z$ranked))
  #' Read and validate one converged full-fit checkpoint
  #' @param modality,count,arm Checkpoint identifiers.
  #' @param extension Whether to use the 1,000-sweep directory.
  #' @return Checkpoint list, or NULL if unavailable or unconverged.
  #' @examples
  #' # read_fit("expression", 2218, "PL")
  read_fit <- function(modality, count, arm, extension = FALSE) {
    path <- file.path(local, if (extension) "iteration1000" else "",
                      sprintf("%s_%d_%s_f0.rds", modality, count, arm))
    if (!file.exists(path)) return(NULL)
    x <- readRDS(path)
    if (is.null(x$fit) || !isTRUE(x$row$converged)) return(NULL)
    stopifnot(nrow(x$scores) == 144L, ncol(x$scores) == 5L,
              ncol(x$fit$EF) == 5L, length(x$features) == nrow(x$fit$EF))
    expected_features <- if (modality == "expression") colnames(z$expression) else z$ranked[seq_len(count)]
    stopifnot(identical(x$features, expected_features))
    if (!is.null(rownames(x$scores))) stopifnot(identical(rownames(x$scores), rownames(z$expression)))
    if (!is.null(rownames(x$fit$EF))) stopifnot(identical(rownames(x$fit$EF), x$features))
    x
  }
  specs <- data.frame(arm = c(rep("PL", 5), rep("PE", 2)),
                      count = c(z$counts, tail(z$counts, 2)))
  fits <- lapply(seq_len(nrow(specs)), function(i) {
    f <- read_fit("methylation", specs$count[i], specs$arm[i])
    if (is.null(f)) stop("Missing converged methylation checkpoint: ", specs$arm[i], " ", specs$count[i])
    f
  })
  expr <- list(PL = read_fit("expression", 2218L, "PL"),
               PE = read_fit("expression", 2218L, "PE", TRUE))
  if (is.null(expr$PL)) stop("Missing converged PL expression checkpoint")
  if (is.null(expr$PE)) message("K = 5 PE expression full fit unavailable or not converged; comparisons omitted.")
  matching <- read.csv(file.path(base, "fold_program_matching.csv"))
  stopifnot(all(matching$K == 5L))
  score_rows <- cca_rows <- stability_rows <- list()
  for (i in seq_along(fits)) {
    m <- fits[[i]]; arm <- specs$arm[i]; count <- specs$count[i]; ex <- expr[[arm]]
    r <- if (is.null(ex)) matrix(NA_real_, 5, 5) else abs(safe_cor(m$scores, ex$scores))
    ij <- expand.grid(methylation_factor = 1:5, expression_factor = 1:5)
    score_rows[[i]] <- data.frame(K = 5L, arm, count, ij, abs_r = as.vector(r))
    cc <- rep(NA_real_, 5)
    if (!is.null(ex) && qr(scale(m$scores, scale = FALSE))$rank == 5L &&
        qr(scale(ex$scores, scale = FALSE))$rank == 5L)
      cc <- stats::cancor(m$scores, ex$scores)$cor
    cca_rows[[i]] <- data.frame(K = 5L, arm, count, dimension = 1:5, canonical_correlation = cc)
    recurrence <- vapply(1:5, function(k) {
      others <- which(specs$arm == arm & specs$count != count)
      v <- unlist(lapply(others, function(j) abs(safe_cor(m$scores[, k, drop = FALSE], fits[[j]]$scores))))
      if (any(is.finite(v))) max(v, na.rm = TRUE) else NA_real_
    }, numeric(1))
    best <- apply(r, 1, function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_)
    best_factor <- apply(r, 1, function(x) if (any(is.finite(x))) which.max(replace(x, !is.finite(x), -Inf)) else NA_integer_)
    nfold <- vapply(1:5, function(k) {
      a <- matching[matching$arm == arm & matching$count == count & matching$full_program == k, ]
      sum(a$fold_converged & is.finite(a$abs_r) & a$abs_r >= 0.5)
    }, integer(1))
    stability_rows[[i]] <- data.frame(K = 5L, arm, count, factor = 1:5,
      closest_expression_factor = best_factor, max_expression_abs_r = best,
      stable_converged_folds = nfold, max_other_count_abs_r = recurrence,
      methylation_specific = classify_specific(best, nfold, recurrence))
  }
  write.csv(do.call(rbind, score_rows), file.path(out, "score_correlations.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cca_rows), file.path(out, "canonical_correlations.csv"), row.names = FALSE)
  stability <- do.call(rbind, stability_rows)
  write.csv(stability, file.path(out, "factor_stability.csv"), row.names = FALSE)

  # The 450K gene and group fields are parallel semicolon-separated lists.
  annotation_data <- new.env(parent = emptyenv())
  annotation_package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  utils::data("Other", package = annotation_package, envir = annotation_data)
  utils::data("Islands.UCSC", package = annotation_package, envir = annotation_data)
  annotation <- as.data.frame(annotation_data$Other)
  islands <- as.data.frame(annotation_data$Islands.UCSC)
  stopifnot(all(rownames(annotation) %in% rownames(islands)))
  annotation$Relation_to_Island <- islands[rownames(annotation), "Relation_to_Island"]
  if (!all(z$ranked %in% rownames(annotation))) stop("Eligible CpGs missing from 450K annotation")
  annotation <- annotation[z$ranked, , drop = FALSE]
  pairs <- expand_annotation(annotation)
  pairs$island_class <- island_class(pairs$island_relation)
  promoter <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
  cpg_island <- pairs[!duplicated(pairs$cpg), c("cpg", "island_class")]
  promoter_cpg <- unique(pairs$cpg[pairs$region %in% promoter])
  body_cpg <- unique(pairs$cpg[pairs$region == "Body"])
  top_rows <- enrich_rows <- list()
  top_ids <- vector("list", length(fits))
  for (i in seq_along(fits)) {
    m <- fits[[i]]; arm <- specs$arm[i]; count <- specs$count[i]
    top_ids[[i]] <- vector("list", 5)
    for (k in 1:5) {
      ids <- m$features[order(abs(m$fit$EF[, k]), decreasing = TRUE)[1:100]]
      top_ids[[i]][[k]] <- ids
      a <- pairs[pairs$cpg %in% ids, , drop = FALSE]
      a$loading <- m$fit$EF[match(a$cpg, m$features), k]
      a$loading_sign <- sign(a$loading)
      top_rows[[length(top_rows) + 1L]] <- data.frame(K = 5L, arm, count, factor = k, a)
      for (category in c("island", "shore", "shelf", "open sea", "promoter", "gene body")) {
        hits <- switch(category, promoter = promoter_cpg, `gene body` = body_cpg,
                       cpg_island$cpg[cpg_island$island_class == category])
        enrich_rows[[length(enrich_rows) + 1L]] <- data.frame(K = 5L, arm, count, factor = k,
          category, top_fraction = mean(ids %in% hits), eligible_fraction = mean(z$ranked %in% hits))
      }
    }
  }
  write.csv(do.call(rbind, top_rows), file.path(out, "top_cpgs_annotated.csv"), row.names = FALSE)
  write.csv(do.call(rbind, enrich_rows), file.path(out, "annotation_enrichment.csv"), row.names = FALSE)
  p <- unique(pairs[pairs$region %in% promoter & pairs$gene %in% colnames(z$expression), c("cpg", "gene")])
  # Rank once per feature, then calculate each requested Spearman correlation.
  mr <- apply(z$methylation[, unique(p$cpg), drop = FALSE], 2, rank)
  er <- apply(z$expression, 2, rank)
  p$spearman_r <- vapply(seq_len(nrow(p)), function(i) stats::cor(mr[, p$cpg[i]], er[, p$gene[i]]), numeric(1))
  if (any(!is.finite(p$spearman_r))) stop("Undefined promoter CpG-gene correlations")
  p$K <- 5L
  write.csv(p[, c("K", "cpg", "gene", "spearman_r")], file.path(out, "promoter_cpg_gene_correlations.csv"), row.names = FALSE)
  overlap <- list()
  for (i in seq_along(fits)) for (k in 1:5) {
    ex <- expr[[specs$arm[i]]]; best <- stability$closest_expression_factor[(i - 1L) * 5L + k]
    v <- c(observed = NA_real_, null_mean = NA_real_, null_sd = NA_real_, p_upper = NA_real_)
    if (!is.null(ex) && !is.na(best)) {
      genes <- ex$features[order(abs(ex$fit$EF[, best]), decreasing = TRUE)[1:100]]
      hits <- unique(pairs$cpg[pairs$gene %in% genes])
      v <- overlap_permutation(z$ranked, hits, top_ids[[i]][[k]])
    }
    subset_r <- p$spearman_r[p$cpg %in% top_ids[[i]][[k]]]
    overlap[[length(overlap) + 1L]] <- data.frame(K = 5L, arm = specs$arm[i], count = specs$count[i],
      methylation_factor = k, expression_factor = best, t(v), promoter_pairs = length(subset_r),
      promoter_fraction_negative = mean(subset_r < 0), promoter_below_minus_0.3 = sum(subset_r < -0.3),
      promoter_below_minus_0.5 = sum(subset_r < -0.5), all_promoter_pairs = nrow(p),
      all_fraction_negative = mean(p$spearman_r < 0), all_below_minus_0.3 = sum(p$spearman_r < -0.3),
      all_below_minus_0.5 = sum(p$spearman_r < -0.5))
  }
  write.csv(do.call(rbind, overlap), file.path(out, "cross_modality_overlap.csv"), row.names = FALSE)
  sets <- list(Hallmark = get_msigdb_collections(collections = list(Hallmark = list(collection = "H")))$Hallmark,
               PDAC = readRDS("results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds"))
  fg <- list()
  warning_rows <- list()
  factor_sets <- c(lapply(seq_along(fits), function(i) list(modality = "methylation", arm = specs$arm[i], count = specs$count[i], fit = fits[[i]])),
                   lapply(names(expr)[!vapply(expr, is.null, logical(1))], function(a) list(modality = "expression", arm = a, count = 2218L, fit = expr[[a]])))
  for (f in factor_sets) for (k in 1:5) {
    weights <- as.numeric(safe_cor(z$expression, f$fit$scores[, k, drop = FALSE]))
    names(weights) <- colnames(z$expression)
    if (any(!is.finite(weights))) stop("Undefined gene correlations in ", f$modality, " ", f$arm, " factor ", k)
    for (collection in names(sets)) {
      result <- withCallingHandlers(
        run_fgsea_program(weights, sets[[collection]], seed = 2026L, collection = collection, scoreType = "std"),
        warning = function(w) {
          warning_rows[[length(warning_rows) + 1L]] <<- data.frame(
            K = 5L, modality = f$modality, arm = f$arm, count = f$count,
            factor = k, collection = collection, message = conditionMessage(w))
          message("fgsea warning [", f$modality, " ", f$arm, " count=", f$count,
                  " factor=", k, " ", collection, "]: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      fg[[length(fg) + 1L]] <- data.frame(K = 5L, modality = f$modality, arm = f$arm,
                                          count = f$count, factor = k, gene_universe = 2218L, result)
    }
  }
  write.csv(do.call(rbind, fg), file.path(out, "fgsea_by_factor.csv"), row.names = FALSE)
  warning_table <- if (length(warning_rows)) do.call(rbind, warning_rows) else
    data.frame(K = integer(), modality = character(), arm = character(), count = integer(),
               factor = integer(), collection = character(), message = character())
  write.csv(warning_table, file.path(out, "analysis_warnings.csv"), row.names = FALSE)
  # Align ABSOLUTE purity by the exact saved TCGA sample identifier. Only
  # factor-level correlations and complete-case counts leave the private data.
  root <- Sys.getenv("PDAC_DATA_ROOT")
  if (!nzchar(root)) root <- yaml::read_yaml("config/globals.yml")$multiomics_data$pdac_root
  purity_path <- file.path(root, "original/TCGA_PAAD.rds")
  if (file.exists(purity_path)) {
    sample_info <- readRDS(purity_path)$sampInfo
    stopifnot(all(c("Tumor.Sample.ID", "ABSOLUTE.Purity") %in% names(sample_info)),
              !anyDuplicated(sample_info$Tumor.Sample.ID),
              all(rownames(z$expression) %in% sample_info$Tumor.Sample.ID),
              is.numeric(sample_info$ABSOLUTE.Purity))
    purity <- sample_info$ABSOLUTE.Purity[match(rownames(z$expression), sample_info$Tumor.Sample.ID)]
    purity_rows <- list()
    for (f in factor_sets) for (k in 1:5) {
      good <- is.finite(purity) & is.finite(f$fit$scores[, k])
      r <- if (sum(good) >= 3L && stats::sd(purity[good]) > 0 &&
               stats::sd(f$fit$scores[good, k]) > 0)
        stats::cor(purity[good], f$fit$scores[good, k]) else NA_real_
      purity_rows[[length(purity_rows) + 1L]] <- data.frame(K = 5L, modality = f$modality,
        arm = f$arm, count = f$count, factor = k, n_complete = sum(good),
        correlation = r, abs_correlation = abs(r))
    }
    write.csv(do.call(rbind, purity_rows), file.path(out, "purity_correlations.csv"), row.names = FALSE)
  } else if (dir.exists(root)) {
    candidates <- list.files(root, pattern = "purity|absolute", recursive = TRUE,
                             full.names = TRUE, ignore.case = TRUE)
    text_files <- list.files(root, pattern = "\\.(csv|tsv|txt)$", recursive = TRUE,
                             full.names = TRUE, ignore.case = TRUE)
    header_candidates <- text_files[vapply(text_files, function(path) {
      header <- tryCatch(readLines(path, n = 2L, warn = FALSE), error = function(e) character())
      any(grepl("purity|absolute", header, ignore.case = TRUE))
    }, logical(1))]
    candidates <- unique(c(candidates, header_candidates))
    if (length(candidates)) stop("Potential local purity files require inspection: ", paste(candidates, collapse = ", "))
    message("No local purity estimates found under PDAC_DATA_ROOT. No purity correlations calculated.")
  } else {
    message("Purity data root is missing. No purity correlations calculated.")
  }
  # Hash private input files without writing patient data into the aggregate output.
  input_paths <- c(design = file.path(local, "design.rds"), matching = file.path(base, "fold_program_matching.csv"),
    pdac_sets = "results/benchmark_sim/outputs/pathway_enrichment/pdac_genesets.rds",
    vapply(seq_len(nrow(specs)), function(i) file.path(local, sprintf("methylation_%d_%s_f0.rds", specs$count[i], specs$arm[i])), character(1)),
    expression_PL = file.path(local, "expression_2218_PL_f0.rds"))
  if (!is.null(expr$PE)) input_paths <- c(input_paths, expression_PE = file.path(local, "iteration1000/expression_2218_PE_f0.rds"))
  if (file.exists(purity_path)) input_paths <- c(input_paths, purity = purity_path)
  annotation_path <- system.file("data", package = annotation_package)
  annotation_files <- list.files(annotation_path, full.names = TRUE)
  input_paths <- c(input_paths, annotation_files)
  hallmark_file <- tempfile(fileext = ".rds")
  saveRDS(sets$Hallmark, hallmark_file, version = 3)
  provenance <- data.frame(K = 5L, input = c(basename(input_paths), "Hallmark gene sets (serialized)"),
    md5 = unname(tools::md5sum(c(input_paths, hallmark_file))))
  unlink(hallmark_file)
  write.csv(provenance, file.path(out, "input_hashes.csv"), row.names = FALSE)
  versions <- data.frame(K = 5L, package = c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "fgsea", "msigdbr", "yaml"),
    version = vapply(c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "fgsea", "msigdbr", "yaml"), function(x) as.character(utils::packageVersion(x)), character(1)))
  write.csv(versions, file.path(out, "package_versions.csv"), row.names = FALSE)
  invisible(NULL)
}

if (sys.nframe() == 0L) characterize_k5()
