# ============================================================
# Script: run_real_multimodal.R
# Purpose: Expression + methylation real-data analysis with ssbmf (rebuild
#          Task 17): on TCGA (n = 144) fit the joint two-block model and the
#          expression-only model through the same ssbmf_fit(), K_init by
#          training CV for each, freeze, score ICGC primary PDAC (n = 50) once,
#          and compare the two with a paired donor-level bootstrap (O6).
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: ssbmf (installed in SSBMF_LIB for parallel CV), survival
# Usage: SSBMF_LIB=<lib> Rscript results/ssbmf/run_real_multimodal.R [cores]   (repository root)
# Data: the matched cache data/processed_multiomics_tcga_icgc.rds in the main
#       checkout (git-ignored), used as built: expression log2 scale (top
#       3,000 genes by combined rank on TCGA), methylation beta-values (top
#       10,000 CpGs by variance on TCGA), KNN-imputed, no scaling. The model's
#       feature intercepts absorb the location of each feature.
# Rules: ICGC is scored once per frozen model and chooses nothing.
# ============================================================

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
source("results/ssbmf/real_helpers.R")
load_ssbmf()
library(survival)
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[1]) else 8L
out_dir <- "results/ssbmf/outputs/real_multimodal"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache <- "/Users/ajwalther/GithubProjects/multiomicsGEP/data/processed_multiomics_tcga_icgc.rds"
d <- readRDS(cache)
tr <- d$training; va <- d$validation_primary
stopifnot(nrow(tr$Y$expression) == 144, nrow(va$Y$expression) == 50)

# Drop features that are constant in training (ssbmf_data refuses them) ----
drop_const <- function(Y) lapply(Y, function(A) A[, apply(A, 2, sd) > 0, drop = FALSE])
Ytr <- drop_const(tr$Y)
dropped <- vapply(names(Ytr), function(m) ncol(tr$Y[[m]]) - ncol(Ytr[[m]]), 0)
message("Constant training features dropped: ", paste(names(dropped), dropped, collapse = ", "))

models <- list(joint = Ytr, expression_only = Ytr["expression"])
K_grid <- c(3L, 5L, 7L, 10L)
ctl <- ssbmf_control()
results <- list(); scores <- list(); diag <- list()
for (mod in names(models)) {
  dd <- ssbmf_data(models[[mod]], tr$time, tr$event, orientation = "samples_by_features")
  t0 <- Sys.time()
  cv <- ssbmf_cv(dd, K_grid = K_grid, nfolds = 5, control = ctl, seed = 1, cores = cores)
  fit <- cv$fit
  message(sprintf("%s: CV and refit %.1f min; K_init = %d", mod, as.numeric(difftime(Sys.time(), t0, units = "mins")), cv$best))
  saveRDS(list(cv = cv, fit = fit), file.path(out_dir, sprintf("frozen_fit_%s.rds", mod)))   # git-ignored
  write.csv(cbind(model = mod, cv$cv), file.path(out_dir, sprintf("cv_by_K_%s.csv", mod)), row.names = FALSE)
  roles <- cbind(model = mod, fit$roles, coef(fit, level = 0.95)[, c("lower", "upper")],
                 surv_share = apply(fit$surv_share, 1, max))
  write.csv(roles, file.path(out_dir, sprintf("frozen_fit_roles_%s.csv", mod)), row.names = FALSE)
  diag[[mod]] <- data.frame(model = mod, K_init = cv$best, converged = fit$converged, sweeps = fit$sweeps,
                            step_decreases = fit$step_decreases, rejections = sum(fit$rejections),
                            n_survival_active = sum(fit$roles$role == "survival_active"),
                            n_reconstruction_only = sum(fit$roles$role == "reconstruction_only"),
                            training_C = cindex(fit$eta, tr$time, tr$event),
                            cv_mean_C_selected = cv$cv$mean_C[cv$cv$K == cv$best],
                            cv_fits_converged = sum(cv$fold_C$converged), cv_fits = nrow(cv$fold_C),
                            cv_fits_with_decrease = sum(cv$fold_C$step_decreases > 0))
  # Freeze, then score ICGC primary once
  Yva <- lapply(names(models[[mod]]), function(m) va$Y[[m]][, colnames(models[[mod]][[m]]), drop = FALSE])
  names(Yva) <- names(models[[mod]])
  scores[[mod]] <- predict(fit, newdata = Yva, orientation = "samples_by_features")
}
diag <- do.call(rbind, diag)
write.csv(diag, file.path(out_dir, "frozen_fit_diagnostics.csv"), row.names = FALSE)

ext <- data.frame(model = names(scores), n = 50, events = sum(va$event),
                  C_icgc = vapply(scores, cindex, 0, time = va$time, status = va$event))
write.csv(ext, file.path(out_dir, "icgc_C.csv"), row.names = FALSE)
bt <- paired_donor_bootstrap(list(ICGC_primary = list(risk_new = scores$joint, risk_ref = scores$expression_only,
                                                      time = va$time, status = va$event, donor = names(va$time))),
                             B = 2000, seed = 1)
write.csv(cbind(comparison = "joint - expression only", bt$per_cohort),
          file.path(out_dir, "icgc_paired_bootstrap.csv"), row.names = FALSE)
print(diag); print(ext); print(bt$per_cohort)
