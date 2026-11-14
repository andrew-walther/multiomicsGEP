# ============================================================
# Script: test_share_targeted_structure.R
# Purpose: Verify structure calculations and signature loading.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: survival, code/simulate_multimodal_yfb.R
# ============================================================
local({
source("code/simulate_multimodal_yfb.R", local = TRUE)
source("results/benchmark_sim/share_targeted_structure.R", local = TRUE)
set.seed(31)
Y <- scale(matrix(rnorm(360), 60, 6), scale = FALSE)
F <- cbind(c(1, 0, 0, 0, 0, 0), c(0, 1, 0, 0, 0, 0))
L <- Y %*% F
fit <- list(EF = F, EL = L, EBeta = c(1, .5), EBeta2 = c(1.01, .26), EF_norms = c(1, 1))
ref <- fit; ref$EF <- cbind(-F[, 2], F[, 1])
full <- list(center = rep(0, 6), arms = list(w1 = list(fit = fit), twostep = list(fit = ref)))
fold <- list(arms = list(w1 = list(fit = ref)))
d <- list(Y = Y, time = rexp(60), status = rep(c(1, 1, 0), 20))
a <- s16_structure_arm(full, rep(list(fold), 5), d, "w1")
recon <- L %*% t(F)
expected <- (sum((Y - recon + tcrossprod(L[, 1], F[, 1]))^2) - sum((Y - recon)^2)) / sum(Y^2)
stopifnot(isTRUE(all.equal(a$program$semipartial_R2[1], expected)),
          all(a$program$active_folds == 5), all(a$program$main_program),
          all(a$matches$estimated_factor[order(a$matches$true_factor)] == c(2, 1)))
z <- Y %*% F
mfull <- survival::coxph(survival::Surv(d$time, d$status) ~ z, control = survival::coxph.control(timefix = FALSE))
zdrop <- z[, 2, drop = FALSE]
mdrop <- survival::coxph(survival::Surv(d$time, d$status) ~ zdrop, control = survival::coxph.control(timefix = FALSE))
stopifnot(isTRUE(all.equal(a$program$delta_logPL[1], as.numeric(logLik(mfull) - logLik(mdrop)))))
missing <- rep(list(fold), 5); missing[[2]] <- list(arms = list())
b <- s16_structure_arm(full, missing, d, "w1")
stopifnot(all(b$program$main_program), all(b$program$active_folds == 4),
          all(b$program$unknown_folds == 1))
inactive <- fold
inactive$arms$w1$fit$EBeta <- c(0, 0)
inactive$arms$w1$fit$EBeta2 <- c(1, 1)
missing[[1]] <- inactive
c <- s16_structure_arm(full, missing, d, "w1")
stopifnot(all(is.na(c$program$main_program)), all(c$program$active_folds == 3))
missing[[3]] <- inactive
q <- s16_structure_arm(full, missing, d, "w1")
stopifnot(!any(q$program$main_program), all(q$program$active_folds == 2))
cache <- tempfile(fileext = ".rds"); subtype <- tempfile(fileext = ".RData")
on.exit(unlink(c(cache, subtype)), add = TRUE)
saveRDS(list(Moffitt_Classical = c("C1", "C2"), Moffitt_BasalLike = "B1"), cache)
schemaList <- c("other", "Elyada_CAF")
subtypeGeneList <- list(data.frame(), data.frame(geneSymbol = c("A", "B", "C"),
  iCAF = c(TRUE, FALSE, TRUE), myCAF = c(FALSE, TRUE, TRUE)))
save(schemaList, subtypeGeneList, file = subtype)
sets <- s16_structure_genesets(cache, subtype)
stopifnot(identical(sets$Moffitt_Classical, c("C1", "C2")),
          identical(sets$Moffitt_BasalLike, "B1"),
          identical(sets$Elyada_iCAF, c("A", "C")),
          identical(sets$Elyada_myCAF, c("B", "C")))
cat("share-targeted structure checks passed\n")
})
