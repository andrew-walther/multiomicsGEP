# ============================================================
# Script: test_feature_selection_helpers.R
# Purpose: Tests for the survival-screen experiment helpers (review focus items 1, 2, 5).
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: survival
# Run from repository root: Rscript results/feature_selection/test_feature_selection_helpers.R
# ============================================================
library(survival)
source("results/feature_selection/feature_selection_helpers.R")
n_pass <- 0L; n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("[PASS]", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("[FAIL]", label, "\n") }
}

set.seed(1)
n <- 300
time <- rexp(n); status <- rbinom(n, 1, 0.7)
risk <- -log(time) + rnorm(n, sd = 0.5)

# harrell_c ----
check("harrell_c matches survival::concordance(reverse = TRUE)",
  isTRUE(all.equal(harrell_c(risk, time, status),
    as.numeric(concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance))))
check("harrell_c is above 0.5 for a true risk score", harrell_c(risk, time, status) > 0.6)
check("harrell_c does not flip an inverted score", harrell_c(-risk, time, status) < 0.4)
check("harrell_c is NA for constant risk", is.na(harrell_c(rep(1, n), time, status)))

# training_direction ----
check("training_direction is 1 for a risk score", training_direction(risk, time, status) == 1L)
check("training_direction is -1 for a good-prognosis score",
  training_direction(-risk, time, status) == -1L)
tr <- 1:200; te <- 201:300
d1 <- training_direction(risk[tr], time[tr], status[tr])
perm <- sample(te)
time2 <- time; status2 <- status; time2[te] <- time[perm]; status2[te] <- status[perm]
check("permuting held-out outcomes leaves the direction unchanged",
  d1 == training_direction(risk[tr], time2[tr], status2[tr]))

# survival_screen_genes ----
p <- 40
Y <- matrix(rnorm(n * p), n, p)
lp <- rowSums(Y[, 1:5]) * 0.6
t_ev <- rexp(n, exp(lp)); cens <- rexp(n, 0.3)
tt <- pmin(t_ev, cens); ss <- as.integer(t_ev <= cens)
cohort <- rep(c("A", "B"), each = n / 2)
Y[, 10] <- 0                                       # constant column
# Column 11 tracks cohort only; cohort B has much higher baseline hazard.
tt[cohort == "B"] <- tt[cohort == "B"] / 5
Y[, 11] <- ifelse(cohort == "B", 3, 0) + rnorm(n, sd = 0.01)
sel <- survival_screen_genes(Y, tt, ss, cohort, 5)
check("screen recovers the five survival genes", setequal(sel, 1:5))
check("screen returns exactly n_keep genes", length(sel) == 5L)
check("constant column gets p = 1 and is counted", attr(sel, "p_values")[10] == 1 &&
  attr(sel, "n_failed") >= 1L)
check("cohort-only column is not in the top 5 when stratified", !(11 %in% sel))
check("n_keep above the candidate count stops",
  inherits(tryCatch(survival_screen_genes(Y, tt, ss, cohort, p + 1), error = identity), "error"))

# unsupervised_gene_set ----
entry <- setNames(sample.int(100), paste0("g", 1:100))
small <- unsupervised_gene_set(entry, 20); big <- unsupervised_gene_set(entry, 50)
check("unsupervised sets are nested", all(small$genes %in% big$genes))
check("unsupervised set reaches the requested size", small$n_shared_actual == 20L)

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0L) quit(status = 1)
