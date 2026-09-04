suppressPackageStartupMessages({library(yaml); library(survival)})
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R"); source("code/update_tau.R")
source("code/compute_elbo.R"); source("code/update_F_cohort.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error=function(e) invisible(NULL)))
source("code/preprocess_desurv.R")

fits <- readRDS("results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_fits.rds")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
pp <- preprocess_merged_cohorts(train_raw, PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
        top_n=cfg$preprocessing$top_n_genes_desurv, rank_transform=FALSE,
        per_platform_standardize=TRUE, normalize_method="none",
        selection_per_cohort=TRUE, selection_method="combined_rank")
train_genes <- pp$gene_names
EXTERNAL_COHORTS <- cfg$pdac$external_cohorts

oriented_cindex <- function(risk, time, status) {
  if (sd(risk)==0) return(NA_real_)
  c <- as.numeric(concordance(Surv(time,status)~risk)$concordance); max(c,1-c)
}

ext_data <- list()
for (ec in EXTERNAL_COHORTS) {
  raw_ext <- load_pdac_raw(ec, PDAC_DATA_ROOT)
  pre_ext <- preprocess_desurv_cohort(raw_ext$Y, raw_ext$gene_names, top_n=NULL,
               log_transform=PLATFORM_LOG_TRANSFORM[[ec]], cohort_name=ec,
               rank_transform=FALSE, per_platform_standardize=TRUE)
  common <- intersect(train_genes, pre_ext$gene_names)
  ext_data[[ec]] <- list(Y_ext=pre_ext$Y[,match(common,pre_ext$gene_names),drop=FALSE],
                          train_idx=match(common,train_genes), time=raw_ext$time, status=raw_ext$status)
}

arms <- c("joint_yfb","joint_yfb_cohort_L","joint_yfb_beta_c","joint_yfb_all_c")
per_cohort_c <- list()
for (arm in arms) {
  fit <- fits[[arm]]
  use_pooled <- !is.null(fit$beta_cohort_id)
  cs <- sapply(EXTERNAL_COHORTS, function(ec) {
    d <- ext_data[[ec]]; EF_sub <- fit$EF[d$train_idx,,drop=FALSE]
    beta <- if (use_pooled) fit$EBeta_pooled else fit$EBeta
    pred <- predict_cox_on_yf(d$Y_ext, EF_sub, beta, EF_norms=fit$EF_norms)
    oriented_cindex(pred$risk_scores, d$time, d$status)
  })
  per_cohort_c[[arm]] <- cs
}
pc <- as.data.frame(per_cohort_c)
cat("=== Per-cohort C by arm ===\n"); print(round(pc,4))

cat("\n=== Leave-one-study-out: mean external C by arm, excluding each cohort in turn ===\n")
loo_rows <- list()
for (excl in EXTERNAL_COHORTS) {
  means <- sapply(arms, function(a) mean(pc[[a]][rownames(pc)!=excl]))
  ranking <- names(sort(means, decreasing=TRUE))
  loo_rows[[excl]] <- data.frame(excluded=excl, t(round(means,4)), ranking_1st=ranking[1], stringsAsFactors=FALSE)
}
loo <- do.call(rbind, loo_rows)
print(loo)
write.csv(loo, "results/benchmark_sim/outputs/cohort_beta_comparison/leave_one_study_out.csv", row.names=FALSE)
write.csv(cbind(cohort=EXTERNAL_COHORTS, pc), "results/benchmark_sim/outputs/cohort_beta_comparison/per_cohort_c_by_arm.csv", row.names=FALSE)
cat("\nFull-sample ranking:", paste(names(sort(colMeans(pc), decreasing=TRUE)), collapse=" > "), "\n")
