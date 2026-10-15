# ============================================================
# Script: run_round2_profiled_cox.R
# Purpose: Evaluate the approved opt-in profiled Cox curvature without external selection.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: existing fixed-design runner and opt-in weighted-centering helper
# Run: Rscript results/benchmark_sim/run_round2_profiled_cox.R [cores|prepare|finalize]
# ============================================================
# Task2-4, committed BEFORE execution, after exact Task2-3 saved-fold identity.
# Andrew explicitly approved this additional causal remedy. One change from
# each model's historical-framework control list: profile_cox_intercept=TRUE.
# Retain per-program/per-factor Normal beta, NO supervised initialization,
# point-Laplace molecular priors, w_S1 for Single, K7, Single no pruning and
# multimodal ordinary pruning. Do not combine a shared prior with this test.
# Profile only the arbitrary constant in the frozen Cox working quadratic:
# Hc=diag(w)-ww'/sum(w). Weighted-center the survival design and old predictor;
# keep raw genomic residuals, raw risk prediction and the original reported ELBOs,
# pruning scores, convergence monitors/rules, caps and tolerances unchanged.
# Multimodal fixed-working objective reports its actual profiled quadratic.
# Single retains its original ELBO monitor; monotonicity is not guaranteed
# by refreshing either working approximation. No criterion is redefined.
# Default/off fits must remain identical to feature/supervised-F.
#
# Simulation:80 NEW cells, unchanged n400/validation120,p50 per modality,
# beta3=3.3,20 seeds20261016:20261035 at noises.25/1 in both models. Same
# original80 K7 simulation references (Single alpha_F0/multimodal10/1).
# Require recovery>=16/20 each noise, null-active<=reference+.5 and every
# fit converged under its original model rule. No K-rank gate until adopted.
#
# TCGA: same approved144 patients/3000genes/10000 imputed CpGs and original
# event-stratified5folds seed42; Single scaling estimated within training fold.
# Screening/imputation cache remains fixed: conditional, not fully nested CV.
# Fit10 new candidate CV cells, reuse exactly10 unchanged R1 framework
# references only after its100-cell aggregate, matching cache/fold hash and
# cell-key/error/provenance checks. No refit, deletion or score substitution.
# Five finite matched C pairs required; mean delta>=-sd(delta)/sqrt5.
# The saved multimodal reference fold1 has C=NA: this comparison therefore
# cannot pass unchanged, even if the new candidate avoids collapse.
# Report the candidate's own mean only if all its5 folds are finite, separately
# from paired delta/SE; never average only the reference's four finite folds.
#
# Also fit both candidates on full TCGA:2 NEW fits. No-collapse requires>=1
# credible program (SD>0, |mean beta|>=1.96SD) AND nonconstant risk. Report
# all caps against the convergence goal; no adoption by ignoring failed gates.
# Total102 records:92 new fits plus10 cached framework CV references.
# Every real-data result table includes labelled frozen benchmark rows.
# No external data/outcomes/predictions for diagnosis or selection. Only a
# final configuration selected without external values can get the specified
# once-only external report. No automatic benchmark replacement or goal claim.
#
# First validate the exact profiled scalar updates/uncertainty using dense
# quadratic calculations, genomic residual preservation and default identity.
# Mathematical derivation: docs/plans/Round2_profiled_curvature_proposal.md.
# Longleaf: round2_profiled_cox.slurm,92 new serial cells,8GB/cell,12h;
# prepare freezes reference metrics then finalize aggregates all102 records.
# Local six sockets, no nested fitter parallelism; runtime/BLAS recorded.
# Dense TCGA fits may exceed1hour. No jobs submitted automatically.
source("results/benchmark_sim/run_round2_fixed_prior.R")
if(sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  out <- "results/benchmark_sim/outputs/round2/profiled"
  if(length(a) && a[1]=="cell") {
    load_panel_code();source("code/concordance_ci.R")
    fixed_prior_cell(as.integer(a[2]),out,candidate="profiled")
  } else if(length(a) && a[1] %in% c("prepare","finalize")) {
    run_fixed_prior_panel(out_dir=out,mode=a[1],candidate="profiled")
  } else run_fixed_prior_panel(if(length(a))as.integer(a[1]) else 6L,
    out_dir=out,candidate="profiled")
}
