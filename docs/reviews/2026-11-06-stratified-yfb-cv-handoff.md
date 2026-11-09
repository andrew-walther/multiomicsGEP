# Stratified YFB cross-validation review handoff

## Scope

This follows commit `a463eea`, which added within-stratum projection centering
to the canonical YFB fitter. The CV changes use its returned training means.
Stratified concordance uses Harrell's within-stratum comparable pairs, pooled
by pair counts. Held-out partial likelihood continues to use each held-out
stratum's risk sets. The rank and likelihood helpers still assign folds by
event status. Centered CV requires each held-out stratum to be represented in
its training fold. A rare stratum missing from training raises the predictor's
existing unknown-label error. This avoids silently choosing a different
centering policy in CV.

A new external stratum can be scored on its own with one fixed center vector,
such as `rep(0, ncol(Y_ext))`. All patients in that stratum receive the same
additive score shift for any change in that center. The ranking and
within-stratum concordance are invariant to the shift. Comparisons between
different strata are outside this interpretation.

## Implementation and verification

| File | Change or check |
|---|---|
| [code/select_K.R](../../code/select_K.R) | Validate and subset YFB stratum labels, project test rows with the fold's fitted center, and score within-stratum pairs. |
| [code/compute_cv_loglik.R](../../code/compute_cv_loglik.R) | Use the predictor for both training diagnostics and test likelihood. Both receive the training fold's center, and the training concordance uses within-stratum pairs. |
| [code/predict_ssbmf_yfb.R](../../code/predict_ssbmf_yfb.R) | Add an example for a previously unseen stratum using a fixed zero center. |
| [tests/test_select_K_cv.R](../../tests/test_select_K_cv.R) | KCV-T20 reconstructs test risk and checks concordance against the stratified formula interface in every fold. Before the fix, the fitter rejected the full-length stratum labels. |
| [tests/test_compute_cv_loglik.R](../../tests/test_compute_cv_loglik.R) | CVLL-T20 checks training concordance, held-out log partial likelihood and the actual held-out risk vector. Before the fix, training C differed from the within-stratum calculation. |
| [tests/test_projection_centering.R](../../tests/test_projection_centering.R) | Compare the zero and within-cohort mean center for an external stratum. The scores differ by a constant and their order agrees. |
| [README.md](../../README.md), [AGENTS.md](../../AGENTS.md), [PROJECT_STATUS.qmd](../../PROJECT_STATUS.qmd), [DECISIONS.md](../../DECISIONS.md), [review index](README.md) | Record usage, verification scope and the outstanding numerical decision. |

`Rscript tests/run_tests.R` exited successfully with 783 passing tests
and zero failures in the shared checkout, including the concurrent flat-β
tests. The suite emitted the bounded-β profile-likelihood fallback warning.
Both unstratified CV helpers returned identical results to their versions at
`f5b91e0` on the synthetic rank-CV fixture. Both centered helpers also
raised the unknown-stratum error in a singleton-stratum check. The tests
compare predictions and concordance values because partial likelihood alone
cannot detect omission of a within-stratum shift.

The CV changes are in the shared working tree. The flat-β code and tests
and the progress-book changes belong to separate work. In the shared
README and AGENTS files, the sections titled “Stratified YFB
cross-validation” and “Stratified YFB scoring” are this handoff's changes.
Legacy and study-specific runners were outside this follow-up.

## Scaling diagnostic

This is a deterministic calculation on a supplied factor state, rather than
a model fit or an adoption study. There are two strata with four rows each.
Gene 1 is the stratum mean (1 or 5) plus ε times the repeated vector
(−1.5, −0.5, 0.5, 1.5). Gene 2 alternates 2 and 3. Set the one loading
column to (1, 0), patient activities to gene 1, and both β moments to zero.

| ε | Within-stratum score SD | ‖F‖ after within-stratum canonicalization | A_gen at τ = 1 | Classification |
|---:|---:|---:|---:|---|
| 0 | 0 | 1 | 104 | genomics only |
| 5 × 10⁻⁷ | 5.976 × 10⁻⁷ | 1 | 104 | genomics only |
| 2 × 10⁻⁶ | 2.390 × 10⁻⁶ | 4.183 × 10⁵ | 5.943 × 10⁻¹⁰ | genomics only |
| 10⁻⁴ | 1.195 × 10⁻⁴ | 8.367 × 10³ | 1.486 × 10⁻⁶ | genomics only |

The corresponding global-centered score SD was about 2.138 throughout.
Canonicalization preserves the outer product and the reported PVE, about
0.667 here. The numerical floors are separate. Below the canonicalization
floor of 10⁻⁶ the rescale is skipped. Just above it the loading can grow
sharply, while activity second moments and therefore the genomics component
of loading precision shrink. At ε = 2 × 10⁻⁶, the constructed state with
τ = 1 and E[β²] = 0 would satisfy the fitter's separate dead-update
condition, max(A_gen + w_S A_surv) < 10⁻⁸. This is a counterexample to
treating that condition as scale-invariant. The run does not show that
actual fits take that trajectory. It also does not show that a factor
explaining genuine expression variation should be removed simply because
it lacks within-stratum survival information.

Reproduce the table from the repository root:

```r
tryCatch(source("code/fit_ssbmf_yfb.R"), error = function(e) {
  if (!exists("canonicalize_yf_factors")) stop(e)
})
source("code/select_K.R")

g <- rep(c("A", "B"), each = 4)
within <- rep(c(-1.5, -0.5, 0.5, 1.5), 2)
for (epsilon in c(0, 5e-7, 2e-6, 1e-4)) {
  Y <- cbind(rep(c(1, 5), each = 4) + epsilon * within, rep(c(2, 3), 4))
  Yc <- Y - (rowsum(Y, g) / as.vector(table(g)))[g, ]
  L <- Y[, 1, drop = FALSE]
  F <- matrix(c(1, 0), ncol = 1)
  for (X in list(scale(Y, center = TRUE, scale = FALSE), Yc)) {
    z <- canonicalize_yf_factors(X, L, L^2, F, F^2, 0, 0, 1)
    stopifnot(isTRUE(all.equal(z$EL %*% t(z$EF), L %*% t(F))))
    cl <- classify_factors(z, Y)
    print(data.frame(epsilon, sd = sd(drop(X %*% F)), divisor = z$scale,
                     F_norm = sqrt(sum(z$EF^2)), A_gen_at_tau1 = sum(z$EL2),
                     pve = cl$PVE, category = cl$category))
  }
}
```

## Open decision for review

Keep the fitter and pruning rule as written for this follow-up. Review the CV
integration and the deliberate missing-stratum error first. For the numerical
issue, the next scientific decision is how to separate an estimable
genomics-only program from a collapsed program in a way that is robust to
factor scaling. An eventual change needs a registered evaluation of planted
recovery, false-active programs and convergence, and should preserve genuine
genomics-only structure. This diagnostic supplies a reproducible edge case
for that evaluation.
