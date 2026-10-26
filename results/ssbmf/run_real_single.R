# ============================================================
# Script: run_real_single.R
# Purpose: Expression-only real-data analysis with ssbmf (rebuild Task 16):
#          choose K_init by training CV on TCGA + CPTAC, freeze the fit, score
#          the five external cohorts once, and compare with the frozen D4
#          benchmark (O5) and a two-step comparator.
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: ssbmf (installed in SSBMF_LIB for parallel CV), survival
# Usage: SSBMF_LIB=<lib> Rscript results/ssbmf/run_real_single.R [cores]   (repository root)
#        Run prepare_real_single.R first.
# Rules: external cohorts are scored once, after the configuration is frozen
#        on training data; they choose nothing. Risk = eta (log hazard), no
#        direction chosen on evaluation data.
# ============================================================

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
source("results/ssbmf/real_helpers.R")
load_ssbmf()
library(survival)
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[1]) else 8L
out_dir <- "results/ssbmf/outputs/real_single"
dat <- readRDS(file.path(out_dir, "data_expression.rds"))
tr <- dat$train

# 1. Training only: K_init by 5-fold CV with Cox strata by cohort ----
dd <- ssbmf_data(list(expression = tr$Y), tr$time, tr$status, strata = tr$cohort,
                 orientation = "samples_by_features")
ctl <- ssbmf_control()                      # provisional defaults, w_S = 1
K_grid <- c(5L, 7L, 10L, 15L)
t0 <- Sys.time()
cv <- ssbmf_cv(dd, K_grid = K_grid, nfolds = 5, control = ctl, seed = 1, cores = cores)
fit <- cv$fit
message(sprintf("CV and refit: %.1f min; selected K_init = %d", as.numeric(difftime(Sys.time(), t0, units = "mins")), cv$best))
saveRDS(list(cv = cv, fit = fit, control = ctl), file.path(out_dir, "frozen_fit_expression.rds"))  # git-ignored
write.csv(cv$cv, file.path(out_dir, "cv_by_K.csv"), row.names = FALSE)
write.csv(cv$fold_C, file.path(out_dir, "cv_folds.csv"), row.names = FALSE)

# Diagnostics of the frozen fit (training only) ----
roles <- cbind(fit$roles, coef(fit, level = 0.95)[, c("lower", "upper")],
               surv_share = fit$surv_share[, 1])
write.csv(roles, file.path(out_dir, "frozen_fit_roles.csv"), row.names = FALSE)
diag <- data.frame(K_init = cv$best, converged = fit$converged, sweeps = fit$sweeps,
                   step_decreases = fit$step_decreases, rejections = sum(fit$rejections),
                   objective = fit$objective$final[["total"]],
                   n_survival_active = sum(fit$roles$role == "survival_active"),
                   n_reconstruction_only = sum(fit$roles$role == "reconstruction_only"),
                   training_C = cindex(fit$eta, tr$time, tr$status),
                   cv_fits_converged = sum(cv$fold_C$converged), cv_fits = nrow(cv$fold_C),
                   cv_fits_with_decrease = sum(cv$fold_C$step_decreases > 0))
write.csv(diag, file.path(out_dir, "frozen_fit_diagnostics.csv"), row.names = FALSE)

# Two-step comparator at the same K_init (EBMF with w_S = 0, then stratified Cox) ----
fit2 <- ssbmf_fit(dd, cv$best, ssbmf_control(w_surv = 0))
keep <- fit2$roles$role != "inactive"
Z2 <- predict(fit2, type = "scores")[, keep, drop = FALSE]
cox2 <- coxph(Surv(tr$time, tr$status) ~ Z2 + strata(tr$cohort))
b2 <- numeric(ncol(fit2$scores)); b2[keep] <- coef(cox2)

# 2. Freeze, then score the five external cohorts once ----
ext <- dat$external
res <- do.call(rbind, lapply(names(ext), function(co) {
  e <- ext[[co]]
  r_new <- score_external(fit, e$Y, "expression")
  f2 <- fit2; f2$beta$mean[] <- b2
  r_two <- score_external(f2, e$Y, "expression")
  ext[[co]]$risk_new <<- r_new; ext[[co]]$risk_two <<- r_two
  data.frame(cohort = co, n = nrow(e$Y), events = sum(e$status), genes_present = e$genes_present,
             C_ssbmf = cindex(r_new, e$time, e$status), C_D4 = cindex(e$frozen_d4_risk, e$time, e$status),
             C_two_step = cindex(r_two, e$time, e$status))
}))
res <- rbind(res, data.frame(cohort = "mean", n = sum(res$n), events = sum(res$events), genes_present = NA,
                             C_ssbmf = mean(res$C_ssbmf), C_D4 = mean(res$C_D4), C_two_step = mean(res$C_two_step)))
write.csv(res, file.path(out_dir, "external_C.csv"), row.names = FALSE)
print(res, digits = 3)

# Donor-level paired bootstrap: ssbmf vs D4, ssbmf vs two-step ----
mk <- function(ref) lapply(ext, function(e) list(risk_new = e$risk_new, risk_ref = e[[ref]], time = e$time,
                                                 status = e$status, donor = e$donor))
bt_d4 <- paired_donor_bootstrap(mk("frozen_d4_risk"), B = 2000, seed = 1)
bt_two <- paired_donor_bootstrap(mk("risk_two"), B = 2000, seed = 1)
boot <- rbind(cbind(comparison = "ssbmf - D4", rbind(bt_d4$per_cohort, bt_d4$mean)),
              cbind(comparison = "ssbmf - two-step", rbind(bt_two$per_cohort, bt_two$mean)))
write.csv(boot, file.path(out_dir, "external_paired_bootstrap.csv"), row.names = FALSE)
print(boot, digits = 3)
o5 <- res$C_ssbmf[res$cohort == "mean"] >= 0.627 - 0.01
cat(sprintf("\nO5 (mean external C >= 0.617): %s (%.4f)\n", if (o5) "met" else "not met", res$C_ssbmf[res$cohort == "mean"]))
