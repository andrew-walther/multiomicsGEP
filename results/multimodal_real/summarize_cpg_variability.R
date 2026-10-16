# ============================================================
# Script: summarize_cpg_variability.R
# Purpose: Describe full-TCGA probe variability without selecting a model or filter.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: base R, ggplot2, yaml (existing project dependencies)
# ============================================================
# Prospective Task5-QC1: observed raw beta values on the unchanged TCGA training
# patients; all available TCGA probes, with the frozen screened10k reference.
# No imputation, model fit, survival-outcome use, SD threshold or new default.
# Mean rule: retain .1<=mean(beta)<=.9. Proportion comparison: remove if >=.95
# of observed beta values are <.1, or >=.95 are >.9 (not combined extremes).
# Display both beta and asin(2*beta-1) sample SD; missing/constant probes stay
# in the inventory and summary denominators. Histograms use80 bins, as the
# earlier distribution checker did; the bins are presentation, not selection.
# Serial local run, typically minutes. Reading the raw four-matrix RDS may
# peak at several GB; block summaries avoid an extra full transformed matrix.
# No nested parallelism or Longleaf submission. For HPC, request8GB/one CPU.
# Run ONLY after the design/code commit: Rscript .../summarize_cpg_variability.R run

# Training-patient alignment ----

#' Restrict raw CpGs to the established training patients, by identifier
#' @param raw_meth Raw beta matrix, probes x patients.
#' @param screened_meth Cached training beta matrix, patients x screened probes.
#' @return Raw matrix on every original probe and the cached patient order.
#' @examples
#' # align_cpg_training(raw_meth, training$Y$methylation)
#' @family cpg_variability
#' @seealso cpg_variability_table
align_cpg_training <- function(raw_meth, screened_meth) {
  ids <- rownames(screened_meth)
  if (!is.matrix(raw_meth) || !is.matrix(screened_meth) ||
      is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      is.null(colnames(raw_meth)) || anyDuplicated(colnames(raw_meth)) ||
      !all(ids %in% colnames(raw_meth)) ||
      !all(colnames(screened_meth) %in% rownames(raw_meth)))
    stop("Raw and cached probe/patient identities must be unique and compatible.")
  raw_meth[, ids, drop=FALSE]
}

# Observed moments ----

#' Calculate observed-value CpG moments in bounded blocks
#' Sample variance is sum_i (beta_i - mean(beta))^2 / (n_observed - 1).
#' Undefined SD at n<2 stays NA. All-missing probes remain in the output.
#' @param beta Numeric raw-beta matrix, probes x subjects, values in [0,1] or NA.
#' @param block_size Number of probes summarized at once; affects memory only.
#' @return Data frame of probe identity, counts, mean, variance/SD, transformed
#'   SD, extreme proportions and descriptive filter membership.
#' @examples
#' cpg_variability_table(rbind(cg1=c(.1,.3,.5),cg2=c(.2,NA,.8)))
#' @family cpg_variability
#' @seealso summarize_cpg_filters
cpg_variability_table <- function(beta, block_size=10000L) {
  if (!is.matrix(beta) || !is.numeric(beta) || !nrow(beta) || !ncol(beta) ||
      is.null(rownames(beta)) || anyNA(rownames(beta)) ||
      any(!nzchar(rownames(beta))) || anyDuplicated(rownames(beta)))
    stop("A nonempty numeric beta matrix with unique probe names is required.")
  if (length(block_size)!=1L || !is.finite(block_size) || block_size<1 ||
      block_size!=as.integer(block_size)) stop("Positive integer block size required.")
  parts <- vector("list",ceiling(nrow(beta)/block_size))
  for (b in seq_along(parts)) {
    index <- seq.int((b-1L)*block_size+1L,min(b*block_size,nrow(beta)))
    x <- beta[index,,drop=FALSE]
    if (any(!is.na(x) & (!is.finite(x) | x<0 | x>1)))
      stop("Observed raw beta values must be finite and within [0,1].")
    n <- as.integer(rowSums(!is.na(x)))
    mu <- rowSums(x,na.rm=TRUE)/n;mu[n==0L] <- NA_real_
    centered <- sweep(x,1,mu,"-")
    variance <- rowSums(centered^2,na.rm=TRUE)/(n-1L)
    variance[n<2L] <- NA_real_
    # Transform each observed value, never the mean or the already-computed SD.
    a <- asin(2*x-1);ma <- rowSums(a,na.rm=TRUE)/n
    va <- rowSums(sweep(a,1,ma,"-")^2,na.rm=TRUE)/(n-1L)
    va[n<2L] <- NA_real_
    low <- rowSums(x<.1,na.rm=TRUE)/n
    high <- rowSums(x>.9,na.rm=TRUE)/n
    low[n==0L] <- high[n==0L] <- NA_real_
    parts[[b]] <- data.frame(probe=rownames(x),n_observed=n,
      n_missing=ncol(x)-n,mean_beta=mu,var_beta=variance,sd_beta=sqrt(variance),
      sd_asin=sqrt(va),prop_below_01=low,prop_above_09=high,
      prop_outside=low+high,keep_mean=mu>=.1 & mu<=.9,
      keep_proportion=low<.95 & high<.95,constant=variance==0,
      stringsAsFactors=FALSE,row.names=NULL)
  }
  ans <- do.call(rbind,parts);rownames(ans) <- NULL;ans
}

# Descriptive subsets ----

#' Summarize filter counts and finite-variability denominators
#' @param tab Output of cpg_variability_table.
#' @param screened_ids Identifiers of the current screened reference probes.
#' @return Data frame including all, each filter, intersections and screened set.
#' @examples
#' # summarize_cpg_filters(tab, c("cg1","cg2"))
#' @family cpg_variability
#' @seealso cpg_variability_plots
summarize_cpg_filters <- function(tab, screened_ids) {
  if (anyNA(screened_ids) || anyDuplicated(screened_ids) ||
      !all(screened_ids %in% tab$probe))stop("Screened probes must be unique and present.")
  mean_ok <- tab$keep_mean %in% TRUE;prop_ok <- tab$keep_proportion %in% TRUE
  screened <- tab$probe %in% screened_ids
  groups <- list("All available"=rep(TRUE,nrow(tab)),"Mean filter"=mean_ok,
    "Proportion filter"=prop_ok,"Both filters"=mean_ok & prop_ok,
    "Mean only"=mean_ok & !prop_ok,"Proportion only"=prop_ok & !mean_ok,
    "Current screened set"=screened)
  rows <- vector("list",length(groups))
  for (i in seq_along(groups)) {
    g <- names(groups)[i];take <- groups[[g]]
    rows[[i]] <- data.frame(subset=g,n_probes=sum(take),
      n_sd_beta=sum(take & is.finite(tab$sd_beta)),
      n_sd_asin=sum(take & is.finite(tab$sd_asin)),
      n_sd_missing=sum(take & !is.finite(tab$sd_beta)),
      n_constant=sum(take & tab$constant %in% TRUE),n_screened=sum(take & screened),
      n_mean_unknown=sum(take & is.na(tab$keep_mean)),
      n_proportion_unknown=sum(take & is.na(tab$keep_proportion)))
  }
  do.call(rbind,rows)
}

# Plots ----

#' Plot observed per-CpG SD distributions and mean-versus-SD density
#' @param tab Output of cpg_variability_table.
#' @param screened_ids Current screened reference probe identifiers.
#' @param n_patients Number of TCGA training patients represented.
#' @return Named list of ggplot objects sd_histogram and mean_sd; no files written.
#' @examples
#' # plots <- cpg_variability_plots(tab, screened_ids,144)
#' @family cpg_variability
#' @seealso summarize_cpg_filters
cpg_variability_plots <- function(tab, screened_ids, n_patients) {
  groups <- list("All available"=rep(TRUE,nrow(tab)),
    "Mean filter"=tab$keep_mean %in% TRUE,
    "Proportion filter"=tab$keep_proportion %in% TRUE,
    "Current screened set"=tab$probe %in% screened_ids)
  long <- list();i <- 0L
  for (g in names(groups))for (scale in c("Beta SD","Arcsine SD")) {
    v <- if(scale=="Beta SD")tab$sd_beta else tab$sd_asin
    take <- groups[[g]] & is.finite(v);i <- i+1L
    long[[i]] <- data.frame(subset=g,scale=scale,SD=v[take])
  }
  d <- do.call(rbind,long);d$subset <- factor(d$subset,levels=names(groups))
  d$scale <- factor(d$scale,levels=c("Beta SD","Arcsine SD"))
  h <- ggplot2::ggplot(d,ggplot2::aes(x=SD))+
    ggplot2::geom_histogram(bins=80,fill="#315a79",colour="white",linewidth=.1)+
    ggplot2::facet_wrap(~scale+subset,scales="free",ncol=4,drop=FALSE)+
    ggplot2::labs(title=paste("Per-CpG variability across",n_patients,"TCGA training patients"),
      subtitle="Observed values before imputation; no SD cutoff applied",
      x="Sample SD across observed patients",y="Number of CpGs")+
    ggplot2::theme_bw(base_size=11)
  take <- is.finite(tab$mean_beta) & is.finite(tab$sd_beta)
  m <- ggplot2::ggplot(tab[take,],ggplot2::aes(x=mean_beta,y=sd_beta))+
    ggplot2::geom_bin_2d(bins=80)+ggplot2::scale_fill_continuous(trans="log10")+
    ggplot2::labs(title="TCGA probe mean and variability (observed beta scale)",
      x="Mean beta across observed patients",y="Sample SD",fill="CpGs")+
    ggplot2::theme_bw(base_size=11)
  list(sd_histogram=h,mean_sd=m)
}

# Descriptive driver ----

#' Run the registered TCGA-only descriptive analysis into a fresh directory
#' @param raw_path RDS with raw TCGA probe-by-patient beta matrix in $meth.
#' @param training_path Existing TCGA-only cache with $Y$methylation.
#' @param benchmark_path Frozen RDS with $joint$training_spec feature/patient IDs.
#' @param out_dir Fresh output directory; never overwritten.
#' @return Aggregate summaries invisibly; local probe inventory RDS, metadata,
#'   warning/error records and two aggregate PNG figures are saved.
#' @examples
#' # run_cpg_variability(raw_path,training_path,benchmark_path,out_dir)
#' @family cpg_variability
#' @seealso cpg_variability_table
run_cpg_variability <- function(raw_path,training_path,benchmark_path,out_dir) {
  if (dir.exists(out_dir))stop("Fresh CpG descriptive output directory required.")
  dir.create(out_dir,recursive=TRUE)
  warnings <- character()
  #' Retain each warning without suppressing its display
  #' @param w Warning condition.
  #' @return Updated process-local warning vector.
  #' @examples
  #' # record_warning(simpleWarning("fixture"))
  record_warning <- function(w) warnings <<- c(warnings,conditionMessage(w))
  #' Preserve the error before stopping this descriptive run
  #' @param e Error condition.
  #' @return Does not return; propagates the original error.
  #' @examples
  #' # record_error(simpleError("fixture"))
  record_error <- function(e) {
    writeLines(conditionMessage(e),file.path(out_dir,"error.txt"));stop(e)
  }
  tryCatch(withCallingHandlers({
    before <- tools::md5sum(c(raw_path,training_path,benchmark_path))
    raw <- readRDS(raw_path)$meth
    raw_n <- ncol(raw)
    # Extract only the molecular matrix and its identities; do not use outcomes.
    screened <- readRDS(training_path)$Y$methylation
    spec <- readRDS(benchmark_path)$joint$training_spec
    frozen <- spec$feature_names$methylation
    if (is.null(frozen) || !setequal(frozen,colnames(screened)) ||
        !identical(spec$subject_ids,rownames(screened)))
      stop("Current probe/patient identities do not match the frozen training manifest.")
    raw <- align_cpg_training(raw,screened);gc()
    tab <- cpg_variability_table(raw)
    cached_tab <- cpg_variability_table(t(screened))
    summary <- summarize_cpg_filters(tab,frozen)
    saveRDS(tab,file.path(out_dir,"probe_inventory.rds"))
    write.csv(summary,file.path(out_dir,"filter_summary.csv"),row.names=FALSE)
    # Reference rows distinguish observed raw values from the unchanged cached
    # imputation. No imputation is rerun, and no missing value is silently filled.
    raw_screened <- tab[match(frozen,tab$probe),]
    reference <- list(observed_all=tab,observed_frozen_screened=raw_screened,
      cached_imputed_frozen_screened=cached_tab)
    ref_rows <- vector("list",length(reference))
    for (i in seq_along(reference)) {
      z <- reference[[i]]
      ref_rows[[i]] <- data.frame(source=names(reference)[i],n_probes=nrow(z),
        n_missing=sum(z$n_missing),n_undefined_sd=sum(!is.finite(z$sd_beta)),
        median_mean_beta=median(z$mean_beta,na.rm=TRUE),
        median_sd_beta=median(z$sd_beta,na.rm=TRUE),
        median_sd_asin=median(z$sd_asin,na.rm=TRUE))
    }
    write.csv(do.call(rbind,ref_rows),file.path(out_dir,"input_reference_summary.csv"),row.names=FALSE)
    metadata <- data.frame(raw_samples=raw_n,training_patients=nrow(screened),
      raw_tcga_probes=nrow(raw),frozen_screened_probes=length(frozen),
      raw_missing_values=sum(tab$n_missing),cached_missing_values=sum(cached_tab$n_missing),
      imputation_rerun=FALSE,mean_low=.1,mean_high=.9,extreme_fraction=.95,
      histogram_bins=80,screened_matches_frozen=TRUE,
      cross_reactive_snp_filter="unverified upstream; no new QC applied")
    write.csv(metadata,file.path(out_dir,"metadata.csv"),row.names=FALSE)
    # Hash patient identities locally; never write patient names to committed CSV.
    ids <- tempfile("cpg_patient_ids_");writeLines(rownames(screened),ids)
    id_hash <- unname(tools::md5sum(ids));unlink(ids)
    provenance <- data.frame(input=c("raw_tcga","tcga_training_cache","frozen_fit"),
      path=c(raw_path,training_path,benchmark_path),md5=unname(before))
    write.csv(provenance,file.path(out_dir,"source_manifest.csv"),row.names=FALSE)
    writeLines(c(paste("training_patient_ids_md5",id_hash),
      "All available TCGA probes are described; this is not a new shared-feature eligibility rule.",
      "Frozen screened IDs are verified against the benchmark fit; its scores are not recomputed.",
      "Raw observed summaries precede imputation; the existing screened cache is reported separately.",
      "No outcomes, survival scores, new filter/default, or SD cutoff are used.",
      paste("R",getRversion(),"ggplot2",utils::packageVersion("ggplot2"))),
      file.path(out_dir,"provenance.txt"))
    plots <- cpg_variability_plots(tab,frozen,nrow(screened))
    ggplot2::ggsave(file.path(out_dir,"cpg_sd_histogram.png"),plots$sd_histogram,
      width=14,height=7,dpi=150)
    ggplot2::ggsave(file.path(out_dir,"cpg_mean_sd.png"),plots$mean_sd,
      width=8,height=5,dpi=150)
    if (!identical(before,tools::md5sum(c(raw_path,training_path,benchmark_path))))
      stop("An input hash changed during the descriptive analysis.")
    print(metadata);print(summary)
    invisible(summary)
  },warning=record_warning),error=record_error,finally={
    write.csv(data.frame(message=warnings),file.path(out_dir,"warnings.csv"),row.names=FALSE)
  })
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (!identical(args,"run"))stop("Use the explicitly preregistered 'run' entry point.")
  cfg <- yaml::read_yaml("config/globals.yml")
  run_cpg_variability(
    file.path(cfg$multiomics_data$root,"tcga","tcga_meth_filtered_sex.rds"),
    "data/round2_tcga_training.rds",
    "results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds",
    "results/multimodal_real/outputs/feature_sensitivity/cpg_qc")
}
