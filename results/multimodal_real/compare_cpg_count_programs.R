# ============================================================
# Script: compare_cpg_count_programs.R
# Purpose: Report provisional loading agreement between saved 10k and 25k methylation fits.
# Created: 2026-10-15
# Dependencies: existing tested program-agreement helpers; base R
# Run after Task5-A1 prospective commit and terminal M3 fit:
# Rscript results/multimodal_real/compare_cpg_count_programs.R
# ============================================================
# Descriptive comparison only. Same statistical formulation and original
# controls; 25k uses log-space arithmetic after the original calculation failed.
# Numerical implementation differs, and the 10k fit capped: this is not a
# clean converged count-sensitivity curve or a biological-label transfer.
source("code/multimodal_program_agreement.R")
out<-"results/multimodal_real/outputs/feature_sensitivity/cpg_count_agreement"
if(file.exists(file.path(out,"program_pairs.csv")))stop("Preserve the completed comparison.")
old_path<-"results/multimodal_real/outputs/feature_sensitivity/methylation_count/cell_01.rds"
new_path<-"results/multimodal_real/outputs/feature_sensitivity/logspace_case/fit.rds"
status_path<-"results/multimodal_real/outputs/feature_sensitivity/logspace_case/case.csv"
if(!file.exists(status_path))stop("Wait for the existing M3 case to terminate; never start another fit here.")
status<-feature_read_terminal(status_path)
if(nzchar(status$error))stop("M3 terminal diagnostics report an error; no comparison is implied.")
paths<-c(old_path,new_path,status_path,"data/round2_tcga_feature_universe.rds",
 "code/multimodal_program_agreement.R","code/feature_sensitivity.R")
hashes<-tools::md5sum(paths)
old<-readRDS(old_path);new<-readRDS(new_path)
z<-readRDS("data/round2_tcga_feature_universe.rds")
# Adapter validates exact patient order and selects each fit's frozen features.
blocks<-list(methylation=z$methylation)
a<-multimodal_feature_adapter(old,blocks);b<-multimodal_feature_adapter(new,blocks)
stopifnot(length(a$genes)==10000L,length(b$genes)==25000L,
 identical(old$training_spec$subject_ids,new$training_spec$subject_ids))
comparison<-multimodal_feature_agreement(old,new,blocks,blocks)
pairs<-comparison$pairs
pairs$old_features<-10000L;pairs$new_features<-25000L
pairs$old_arithmetic<-"original";pairs$new_arithmetic<-"logspace"
pairs$old_converged<-old$diagnostics$converged;pairs$new_converged<-new$diagnostics$converged
pairs$old_iterations<-old$diagnostics$iterations;pairs$new_iterations<-new$diagnostics$iterations
# Keep unmatched or undefined factors; never replace NA by zero correlation.
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(pairs,file.path(out,"program_pairs.csv"),row.names=FALSE)
write.csv(comparison$activity_correlation,file.path(out,"activity_correlation.csv"))
write.csv(comparison$projection_correlation,file.path(out,"projection_correlation.csv"))
write.csv(comparison$loading_correlation$methylation,file.path(out,"common_cpg_loading_correlation.csv"))
write.csv(data.frame(old_features=10000L,new_features=25000L,
 risk_correlation=comparison$risk_correlation,old_converged=old$diagnostics$converged,
 new_converged=new$diagnostics$converged),file.path(out,"risk_agreement.csv"),row.names=FALSE)
stopifnot(file.copy("results/multimodal_real/outputs/feature_sensitivity/logspace_case/results_with_benchmark.csv",
 file.path(out,"results_with_benchmark.csv"),overwrite=FALSE))
stopifnot(identical(hashes,tools::md5sum(paths)))
write.csv(data.frame(path=paths,md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
writeLines(c("Provisional saved-fit agreement, not a complete converged CpG count curve.",
 "The 25k statistical model is unchanged but its numerical Cox implementation is logspace.",
 "Match by absolute activity correlation without reuse; report signed loading r on common CpGs.",
 "Keep undefined/unmatched factors and matching margins; no threshold or transferred biological label.",
 "Companion results_with_benchmark.csv preserves frozen contextual scores without rescoring.",
 capture.output(sessionInfo())),file.path(out,"interpretation.txt"))
print(pairs[,c("old_factor","new_factor","activity_correlation","loading_correlation")])
