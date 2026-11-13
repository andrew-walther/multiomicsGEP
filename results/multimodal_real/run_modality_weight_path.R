# ============================================================
# Script: run_modality_weight_path.R
# Purpose: Describe four fixed TCGA reconstruction-weight fits at 10000 CpGs.
# Author: Andrew Walther
# Created: 2026-11-12
# Dependencies: yaml, ebnm, survival, IlluminaHumanMethylation450kanno.ilmn12.hg19
# ============================================================
# One sequential worker and one numerical thread. No CV or external scoring.
# Patient-level fits and scores remain under ignored data/.
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
source("results/multimodal_real/run_multimodal_tcga_comparison.R")
source("results/benchmark_sim/run_modality_weight_simulation.R")

#' Classify live programs and match their patient scores to fixed single-modality fits
#' @param fit Fitted posterior moments, precision and optional modality powers.
#' @param groups Feature-order modality labels.
#' @param scores Patient-by-program joint projection scores.
#' @param expression,methylation Patient-by-program reference projection scores.
#' @return One aggregate row per program, retaining undefined score correlations.
#' @examples
#' # mwpath_programs(fit, groups, Y %*% fit$EF, expression_scores, methylation_scores)
mwpath_programs <- function(fit, groups, scores, expression, methylation) {
  p <- mwsim_precision(fit,groups)
  e <- p[p$modality=="expression",]; m <- p[p$modality=="methylation",]
  live <- e$precision+m$precision>0
  classification <- rep("dead",nrow(e))
  classification[live & (m$empty | (!is.na(m$share) & m$share<.05))] <- "expression-exclusive"
  classification[live & (e$empty | (!is.na(e$share) & e$share<.05))] <- "methylation-exclusive"
  classification[live & !is.na(e$share) & !is.na(m$share) & e$share>=.05 & m$share>=.05] <- "joint"
  #' Find the largest defined absolute score correlation for each joint program
  #' @param reference Patient-by-program single-modality score matrix.
  #' @return Numeric vector with NA if all reference comparisons are undefined.
  #' @examples
  #' # best(expression)
  best <- function(reference) apply(suppressWarnings(abs(cor(scores,reference))),1,function(r)
    if (any(is.finite(r))) max(r[is.finite(r)]) else NA_real_)
  er <- best(expression); mr <- best(methylation)
  data.frame(program=e$program,live=live,expression_share=e$share,methylation_share=m$share,
    expression_empty=e$empty,methylation_empty=m$empty,class=classification,
    expression_abs_r=er,methylation_abs_r=mr,
    new_program=live & is.finite(er) & is.finite(mr) & er<.5 & mr<.5)
}

#' Extract top loadings for every exclusive or new program with 450K annotations
#' @param fit Fitter result.
#' @param programs Program classification table.
#' @param features Feature identifiers in fitted order.
#' @param groups Modality labels in fitted order.
#' @param annotation CpG annotation data frame indexed by CpG identifier.
#' @return Aggregate top-20 loading table; expression genes have no CpG location.
#' @examples
#' # mwpath_top(fit, programs, features, groups, annotation)
mwpath_top <- function(fit,programs,features,groups,annotation) {
  result <- list()
  for (k in programs$program[programs$live & (programs$class!="joint" | programs$new_program)]) {
    cls <- programs$class[programs$program==k]
    mods <- if (cls=="expression-exclusive") "expression" else if(cls=="methylation-exclusive") "methylation" else unique(groups)
    for (m in mods) {
      ix <- which(groups==m)
      ix <- head(ix[order(-abs(fit$EF[ix,k]),features[ix])],20)
      a <- annotation[match(features[ix],rownames(annotation)),,drop=FALSE]
      result[[length(result)+1L]] <- data.frame(program=k,modality=m,rank=seq_along(ix),
        feature=features[ix],loading=fit$EF[ix,k],
        annotated_gene=if(m=="expression") features[ix] else a$UCSC_RefGene_Name,
        gene_region=if(m=="expression") NA_character_ else a$UCSC_RefGene_Group,
        island_relation=if(m=="expression") NA_character_ else a$Relation_to_Island,
        promoter=if(m=="expression") NA else grepl("(^|;)(TSS1500|TSS200|5'UTR|1stExon)(;|$)",a$UCSC_RefGene_Group))
    }
  }
  if (!length(result)) return(data.frame(program=integer(),modality=character(),rank=integer(),
    feature=character(),loading=numeric(),annotated_gene=character(),gene_region=character(),
    island_relation=character(),promoter=logical()))
  do.call(rbind,result)
}

#' Run the four prospectively fixed descriptive TCGA fits and summarize the reused baseline
#' @return Invisible list of fit and program summaries; writes aggregate CSVs only.
#' @examples
#' # mwpath_real_run()
mwpath_real_run <- function() {
  out <- "results/multimodal_real/outputs/modality_weight_path_2026-11-12"
  private <- "data/modality_weight_path_2026-11-12"
  if(dir.exists(out)||dir.exists(private)) stop("Refusing to overwrite a weight-path run.")
  input <- c(design="data/methylation_single_modality_2026-11-12/K05/design.rds",
    expression="data/methylation_single_modality_2026-11-12/K05/expression_2218_PL_f0.rds",
    methylation="data/methylation_single_modality_2026-11-12/K05/methylation_10000_PL_f0.rds",
    unweighted="data/multimodal_tcga_comparison_2026-11-12/joint_grouped_fold0.rds")
  hashes <- tools::md5sum(input)
  cfg <- yaml::read_yaml("config/globals.yml")
  design <- readRDS(input[["design"]]); expression <- readRDS(input[["expression"]]); methylation <- readRDS(input[["methylation"]])
  baseline <- readRDS(input[["unweighted"]])
  stopifnot(length(design$time)==144L,ncol(design$expression)==2218L,cfg$benchmark$k_working_yfb==5L,
    isTRUE(expression$row$converged),isTRUE(methylation$row$converged),isTRUE(baseline$row$converged))
  x <- mtcga_prepare_fold(design,0L)
  Y <- do.call(cbind,x$training); groups <- rep(names(x$training),vapply(x$training,ncol,integer(1)))
  stopifnot(identical(colnames(Y),baseline$features),identical(groups,baseline$groups))
  # Read the annotation installed for task C without modifying its library or data.
  .libPaths(c(normalizePath("data/methylation_single_modality_2026-11-12/K05/R-library",mustWork=TRUE),.libPaths()))
  ann <- new.env(parent=emptyenv())
  utils::data("Other",package="IlluminaHumanMethylation450kanno.ilmn12.hg19",envir=ann)
  utils::data("Islands.UCSC",package="IlluminaHumanMethylation450kanno.ilmn12.hg19",envir=ann)
  annotation <- as.data.frame(ann$Other)
  islands <- as.data.frame(ann$Islands.UCSC)
  annotation$Relation_to_Island <- islands[rownames(annotation),"Relation_to_Island"]
  stopifnot(all(colnames(x$training$methylation)%in%rownames(annotation)))
  dir.create(out,recursive=TRUE);dir.create(private,recursive=TRUE)
  selection_path <- "results/benchmark_sim/outputs/feature_group_simulation_2026-11-12_weight_path/selection.csv"
  if(!file.exists(selection_path)) stop("Part A selection must be recorded before descriptive TCGA fitting.")
  selection <- read.csv(selection_path)
  weights <- c(unweighted=1,w070=.7,w_sqrt5=1/sqrt(5),w030=.3,w_count=2218/10000)
  writeLines(c("registration=docs/prompts/2026-11-12-codex-d-weight-path-prompt.md",
    "scope=four new full-data fits; reused unweighted baseline; no CV; no external scoring",
    "count_weight_exact=2218/10000=0.2218; prompt rounded display 0.22",
    "data=144 matched patients;2218 D4 genes;10000 ranked CpGs;all-patient feature z-scores",
    "fit=grouped point-Laplace L/F;Normal beta;K5;canonicalize TRUE;wS1;common SVD starts",
    "classification=exclusive if other block empty or precision<5%;joint if both>=5%",
    "new_program=best absolute patient-score correlation<0.5 with both fixed single-modality fits",
    "promoter=any 450K UCSC region TSS1500,TSS200,5'UTR,1stExon;island=450K Relation_to_Island",
    "workers=1;numerical_threads=1"),file.path(out,"settings.txt"))
  sources <- c("results/multimodal_real/run_modality_weight_path.R","results/multimodal_real/run_multimodal_tcga_comparison.R",
    "results/benchmark_sim/run_modality_weight_simulation.R","code/fit_ssbmf_yfb.R","code/modality_weights_yfb.R","config/globals.yml")
  write.csv(data.frame(path=sources,md5=unname(tools::md5sum(sources))),file.path(out,"source_hashes.csv"),row.names=FALSE)
  summaries <- programs <- tops <- warning_rows <- list()
  init <- mtcga_init(Y,5L)
  for (arm in names(weights)) {
    w <- weights[[arm]];warn <- character();started <- proc.time()[["elapsed"]]
    if(arm=="unweighted") { fit <- baseline$fit; seconds <- baseline$row$seconds } else {
      fit <- withCallingHandlers(fit_ssbmf_yfb(Y,design$time,design$event,K=5L,
        max_iter=cfg$cavi$max_iter,tol=cfg$cavi$tol,alpha=cfg$benchmark$alpha,
        prior_LF="point_laplace",prior_beta="normal",canonicalize=TRUE,w_S=1,
        feature_group=groups,modality_weight=c(expression=1,methylation=w),
        EL_init=init$EL_init,EF_init=init$EF_init,verbose=FALSE),
        warning=function(e){warn<<-c(warn,conditionMessage(e));invokeRestart("muffleWarning")})
      seconds <- proc.time()[["elapsed"]]-started
      saveRDS(list(fit=fit,features=colnames(Y),groups=groups,scaling=list(center=x$center,sd=x$sd),
        scores=Y%*%fit$EF,weight=w,warnings=warn),file.path(private,paste0(arm,".rds")))
    }
    scores <- Y%*%fit$EF
    p <- mwpath_programs(fit,groups,scores,expression$scores,methylation$scores)
    p$arm<-arm;p$weight<-w
    # A weight is validated only if Part A explicitly selected that same value.
    selected <- if("selected_weight"%in%names(selection)) selection$selected_weight[[1L]] else NA_real_
    validated <- is.finite(selected) && abs(w-selected)<1e-12
    p$simulation_status <- if(validated) "usable in registered simulation" else "not validated by simulation"
    programs[[arm]]<-p
    t <- mwpath_top(fit,p,colnames(Y),groups,annotation); t$arm<-rep(arm,nrow(t));t$weight<-rep(w,nrow(t));tops[[arm]]<-t
    summaries[[arm]]<-data.frame(arm,weight=w,reused=arm=="unweighted",converged=fit$history$converged,
      iterations=fit$history$n_iter,seconds,live=sum(p$live),empty_expression=sum(p$expression_empty),
      empty_methylation=sum(p$methylation_empty),expression_exclusive=sum(p$class=="expression-exclusive"),
      methylation_exclusive=sum(p$class=="methylation-exclusive"),joint=sum(p$class=="joint"),
      new_programs=sum(p$new_program),warnings=length(warn),simulation_status=p$simulation_status[1L])
    if(length(warn))warning_rows[[arm]]<-data.frame(arm,warning=warn)
    for(nm in c("summaries","programs","tops"))write.csv(do.call(rbind,get(nm)),file.path(out,paste0(nm,".csv")),row.names=FALSE,na="NA")
    if(length(warning_rows))write.csv(do.call(rbind,warning_rows),file.path(out,"warnings.csv"),row.names=FALSE)
    message("Completed ",arm," w=",w," iterations=",fit$history$n_iter," seconds=",round(seconds,1))
  }
  if(!identical(hashes,tools::md5sum(input)))stop("Read-only input hash changed.")
  write.csv(data.frame(path=unname(input),md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
  invisible(list(summaries=do.call(rbind,summaries),programs=do.call(rbind,programs)))
}
if(sys.nframe()==0L)mwpath_real_run()
