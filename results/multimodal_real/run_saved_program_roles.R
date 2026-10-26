# ============================================================
# Script: run_saved_program_roles.R
# Purpose: Characterize saved joint programs against the frozen benchmark without refitting.
# Created: 2026-10-15
# Dependencies: yaml and tested multimodal program-agreement helpers
# Run only after Task5-Roles1 commit: Rscript results/multimodal_real/run_saved_program_roles.R
# ============================================================
# Descriptive audit of two existing full-TCGA fits, both unconverged. No model,
# new survival score, external outcome, biological label transfer or selection.
# Norm shares are within-modality rank-one norm fractions, not additive PVE
# or survival-explained fractions. Coefficients retain fitted orientation.
source("code/multimodal_program_agreement.R")
cfg <- yaml::read_yaml("config/globals.yml")
out <- "results/multimodal_real/outputs/feature_sensitivity/saved_program_roles"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths <- c(benchmark="results/benchmark_snapshot_2026-10-14/multimodal/fits_K07_point_laplace-point_laplace_pruned.rds",
  profiled="results/benchmark_sim/outputs/round2/profiled/cells/cell_102.rds",
  training="data/round2_tcga_training.rds")
before <- tools::md5sum(paths)
old <- readRDS(paths[["benchmark"]])$joint
new <- readRDS(paths[["profiled"]])
Y <- readRDS(paths[["training"]])$Y
stopifnot(identical(old$training_spec,new$training_spec),
  !old$diagnostics$converged,!new$diagnostics$converged,
  old$diagnostics$iterations==500L,new$diagnostics$iterations==500L)
roles_old <- multimodal_feature_roles(old,Y,cfg$multimodal_yfb$prognostic_z_threshold)
roles_new <- multimodal_feature_roles(new,Y,cfg$multimodal_yfb$prognostic_z_threshold)
roles <- rbind(cbind(source="benchmark_snapshot",roles_old),
  cbind(source="profiled_unadopted",roles_new))
roles$converged <- FALSE;roles$iterations <- 500L
write.csv(roles,file.path(out,"program_roles.csv"),row.names=FALSE)
a <- multimodal_feature_agreement(old,new,Y,Y)
p <- a$pairs
p$old_source <- "benchmark_snapshot";p$new_source <- "profiled_unadopted"
p$old_norm_share <- roles_old$norm_share[match(paste(p$old_factor,p$modality),
  paste(roles_old$factor,roles_old$modality))]
p$new_norm_share <- roles_new$norm_share[match(paste(p$new_factor,p$modality),
  paste(roles_new$factor,roles_new$modality))]
write.csv(p,file.path(out,"program_matches.csv"),row.names=FALSE)
write.csv(a$activity_correlation,file.path(out,"activity_correlation.csv"))
write.csv(a$projection_correlation,file.path(out,"projection_correlation.csv"))
for(m in names(a$loading_correlation))write.csv(a$loading_correlation[[m]],
  file.path(out,paste0(m,"_loading_correlation.csv")))
write.csv(data.frame(source=c("benchmark_snapshot","profiled_unadopted"),
  risk_correlation_to_benchmark=c(1,a$risk_correlation),
  credible_programs=c(sum(roles_old$credible[roles_old$modality=="expression"]),
    sum(roles_new$credible[roles_new$modality=="expression"])),
  converged=FALSE,iterations=500L),file.path(out,"comparison_summary.csv"),row.names=FALSE)
stopifnot(identical(before,tools::md5sum(paths)))
write.csv(data.frame(path=unname(paths),md5=unname(before)),file.path(out,"input_hashes.csv"),row.names=FALSE)
writeLines(c("Both saved fits hit500 iterations; roles are descriptive and not a model adoption.",
  "No new survival scoring or biological enrichment/label transfer was performed.",
  "Shares normalize rank-one molecular norms separately in each modality; not PVE or survival contributions.",
  "Greedy activity matching preserves undefined and ambiguous cases; coefficients retain fitted orientation.",
  capture.output(sessionInfo())),file.path(out,"interpretation.txt"))
print(read.csv(file.path(out,"comparison_summary.csv")))
