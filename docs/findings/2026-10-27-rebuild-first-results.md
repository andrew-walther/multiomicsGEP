# Rebuilt survival-supervised factorization (`ssbmf`): first results

Written 2026-10-26/27 for readers outside the project. Branch
`rebuild/unified-ssbmf` (not merged, not pushed). Every choice described
here as provisional was made without review and is listed in Section 5 for
Andrew to accept or change.

## 1. Summary

The method was rebuilt from a written derivation as an R package, `ssbmf`,
with one fitting routine for one or more data blocks. The rebuilt fit does
what the current code did not: it maximizes one stated objective, never
lowered that objective in any simulation or real-data fit, and converged in
every simulation fit, including 2,000 expression plus 10,000 correlated
methylation features. It recovers high-variance programs and, on the
expression-only PDAC data, predicts the five external cohorts as well as the
frozen benchmark (mean C 0.628 against 0.627).

It does not yet do the thing a supervised factorization is for. At the
unweighted setting ($w_S = 1$) the Cox term supplies about $10^{-5}$ to
$4\times10^{-4}$ of the loading-update precision for a prognostic program (about
$5\times10^{-7}$ on the real data), so survival hardly moves the loadings: a
program explaining 1% of the variance is found equally well without
supervision, and one explaining 0.25% is found less often than by the
two-step comparator. A survival weight chosen by training cross-validation
(the last-resort option) makes that program survival-active and raises
validation C, but the recovered loadings still match the planted program
only partly (|r| 0.40 against 0.31 at $w_S = 1$; the threshold is 0.7), and
its coefficient intervals then cover the truth only 27% of the time. On
TCGA expression plus methylation, the joint fit converges but its shared β
prior shrinks to zero, so it predicts ICGC worse than the expression-only
fit (C 0.554 vs 0.707).

## 2. Criteria for a working model (O1–O8)

| # | Criterion | Threshold | Result | Status |
|---|---|---|---|---|
| O1 | Monotone objective | no step lowers it by > 1e-8 × \|obj\| | 180/180 simulation fits and every real-data fit; worst relative step −9.9e-13; 0 rejected steps at $w_S = 1$ | met |
| O2 | Convergence | rel. change < 1e-6 within 1,000 sweeps in ≥ 95% of simulation fits and every real fit at the selected K | simulation 180/180 (median 17–97 sweeps, max 279); expression real data: selected fit converged (297 sweeps), 19/20 CV fold fits; multimodal: both selected fits converged (147 and 343 sweeps), 17/20 and 19/20 CV fold fits; one CV fold fit at a selected K_init (expression-only, K_init = 7) stopped at 1,000 sweeps without converging | met for simulation and every frozen fit; one CV fold fit at a selected K did not converge |
| O3 | High-variance program recovery | best-match \|r\| ≥ 0.9 in ≥ 90% | 630/630 (one and two blocks, two noise levels, K_init 4, 6, 8) | met |
| O4 | Low-variance prognostic recovery | matched (\|r\| ≥ 0.7), survival-active, correct sign in ≥ 80%, and more often than two-step | 1% of variance: 100% for ssbmf and for two-step (no advantage); 0.25%: ssbmf 13–20% vs two-step 27%; with CV-chosen $w_S$: 0.25% cell 0/15 at \|r\| ≥ 0.7 (mean \|r\| 0.40), survival-active 12/15 | not met |
| O5 | Expression-only real data | mean external C ≥ 0.617 | 0.628 (D4 0.627; paired donor-level difference +0.002, 95% interval −0.021 to 0.024) | met |
| O6 | Multimodal real data | both fits converge; ICGC C with paired donor-level interval | both converged; ICGC primary C: joint 0.554, expression-only 0.707; paired donor-level difference −0.153 (95% interval −0.266 to −0.051) | met as defined (the joint fit is worse) |
| O7 | Coefficient uncertainty | 90% intervals cover truth in 80–95% | 82% overall (577/701); one block 85–90%, two blocks 70–75% | met overall; not met in two-block cells |
| O8 | Package | `R CMD check` with no errors or warnings; quickstart and vignettes run | check Status OK (0 errors, 0 warnings, 0 notes); quickstart runs (about 9 s); vignettes not written | partly met |

## 3. What the rebuilt model does differently from the current code

| Topic | Current code (`code/fit_cox_on_yf.R`, `code/fit_multimodal_yfb.R`) | `ssbmf` |
|---|---|---|
| Implementations | separate single-modality and multimodal fitters with many option flags | one `ssbmf_fit()`; one block is the M = 1 case |
| Objective | working quadratic expansion of the Cox term, refreshed during fitting (not a bound; the refreshed objective fell in 18 of 50 refreshes in an earlier audit) | Gaussian-moment (Jensen) bound $\sum_i\delta_i[\bar\eta_i-\log\sum_{j\in R_i}e^{\bar\eta_j+v_j/2}]$, a deterministic function of the variational moments; checked below Monte Carlo values of $E_q\log PL$ in 24/24 random states (median gap 0.49 nats) |
| Monotonicity | not guaranteed; not checked per step | every step is checked against the exact objective; steps that would lower it are shortened or rejected and counted |
| Survival in the loadings | single-modality benchmark used $\alpha_F = 0$ (survival never enters F); multimodal update dropped the within-factor coupling term | survival enters every loading update through its exact gradient in the mean parameters (conjugate-computation messages), including the coupling term |
| KL bookkeeping | KL evaluated against stale priors in the multimodal code; β prior variance not carried through rescaling | every prior is refitted inside the same `ebnm` call that returns the posterior, and the KL is computed against that prior |
| Scale | per-iteration canonicalization and sign correction inside the loop | no rescaling inside the loop; a canonical scale (unit SD of each program score) only for reported output |
| β prior | per-program or shared priors that could absorb at β = 0 | shared Normal prior with variance fitted by marginal likelihood of the score-based pseudo-observations; zero is a fixed point only when no program carries survival signal (tested: β leaves zero within 5 sweeps); on the joint real-data fit the variance did fall to about 2e-7, because no program carried signal at the converged loadings |
| Intercept, cohorts | z-scoring per cohort; optional cohort offsets in L | feature intercepts in every block; cohorts as Cox strata only |
| Program roles | delete-and-refit pruning (multimodal) or `classify_factors()` | `ssbmf_classify()` on one fit: survival-active (lfsr ≤ 0.05), reconstruction-only (≥ 0.1% of a block's variance), inactive |
| Uncertainty in β | point estimates in reports | posterior mean, SD, intervals and lfsr per program from `coef()` |
| Numerics | — | log-space risk sets (finite at \|η\| = 1e4); η scale stable from 500 to 5,000 correlated features; `ebnm` solved on rescaled data, impossible `ebnm` results (from programs switching off) replaced by the no-information solution |

The derivation is `derivations/unified_ssbmf/ssbmf_derivation.{qmd,pdf}`;
`verify_derivation.R` reproduces the bound check and a finite-difference
stationarity check (all gradients after an update ≤ 7.7e-7 per posterior SD).

## 4. Results

### 4.1 Simulation at $w_S = 1$

Registered design: `results/ssbmf/sim_design.md` (written before any run).
n = 200 training and 200 validation samples, 4 true programs on disjoint
sets of 10% of the features, features in groups of 20 sharing a nuisance
factor, Weibull survival with 30% censoring, K_init ∈ {4, 6, 8}, 5 datasets
per cell. Arms: ssbmf with point-Laplace loadings (provisional default),
ssbmf with point-exponential loadings, and the two-step comparator (the same
model at $w_S = 0$, then Cox on the program scores).

| Cell | Design | Converged | Max sweeps | O3 | Low-variance program recovered | O7 (90% coverage) | Validation C |
|---|---|---|---|---|---|---|---|
| A | 2,000 features, noise SD 1 | 45/45 | 143 | 100% | — | 85–92% | 0.668 (all arms) |
| A2 | 2,000 features, noise SD 2 | 45/45 | 50 | 100% | — | 85–90% | 0.668–0.669 |
| B | 2,000 + 10,000 features | 45/45 | 279 | 100% | — | 70–80% | 0.700 |
| C | prognostic program at 1% of variance | 45/45 | 63 | 100% | 15/15 every arm | 87–95% | 0.699–0.701 |
| C2 | as C, two blocks | 45/45 | 101 | 100% | 15/15 every arm | 75–85% | 0.712–0.713 |
| C3 | prognostic program at 0.25% of variance | 45/45 | 75 | 100% | ssbmf 2/15 and 3/15; two-step 4/15 | 86–94% | ssbmf 0.543–0.557; two-step 0.584 |

Ranges are over arms (point-Laplace, point-exponential, two-step). Cell C3
was added after one pilot fit showed that every arm recovers the 1% program;
no method setting was changed. Point-exponential loadings converged in fewer
sweeps (median 17 vs 92 in the two-block cell). The survival share of the
loading-update precision for the prognostic program was a median 4e-4 (C),
1e-5 (C2) and 1e-13 (C3, where it was usually not found, so its β stayed
near 0 and the program received no survival information).

### 4.2 Balancing survival and reconstruction (rebuild Part 4)

Because O4 failed at $w_S = 1$, the options were tried in the plan's order on
the low-variance cells, with the same datasets:

1. **Effective-number-of-features block weight** (no tuning parameter;
   $w_m = p_{\mathrm{eff},m}/p_m$ from the residual covariance). It set
   $w_m \approx 0.23$ (one block) and 0.06 (two blocks), weakening the
   reconstruction information the low-variance program needs: recovery fell to
   3/15 (C), 0/15 (C2) and 0/15 (C3). **Dropped.**
2. **Non-absorbing β prior.** Already the default.
3. **Survival weight chosen by training cross-validation** (last resort),
   grid $w_S \in \{1, 10, 100\}$, one-standard-error rule preferring the
   smallest weight. Fits at $w_S \ge 1{,}000$ were dropped from the grid after
   one pilot dataset (Section 5, item 9).

| Cell | Arm | Low-variance program recovered | Mean best \|r\| | Survival-active | $w_S$ chosen | O7 coverage | Validation C |
|---|---|---|---|---|---|---|---|
| C (1%) | $w_S = 1$ | 15/15 | 0.965 | 15/15 | — | 87% | 0.701 |
| C | effective-feature weight | 3/15 | 0.53 | 4/15 | — | 89% | 0.548 |
| C | CV-chosen $w_S$ | 15/15 | 0.965 | 15/15 | 1 (12), 10 (3) | 67% | 0.693 |
| C | two-step | 15/15 | 0.965 | 15/15 | — | 95% | 0.701 |
| C2 (1%, two blocks; K_init = 6 for the CV arm) | $w_S = 1$ | 15/15 | 0.964 | 15/15 | — | 75% | 0.713 |
| C2 | effective-feature weight | 0/15 | 0.12 | 4/15 | — | 93% | 0.550 |
| C2 | CV-chosen $w_S$ | 5/5 | 0.964 | 5/5 | 1 (5) | 63% | 0.723 |
| C2 | two-step | 15/15 | 0.964 | 15/15 | — | 85% | 0.713 |
| C3 (0.25%) | $w_S = 1$ | 2/15 | 0.31 | 2/15 | — | 88% | 0.543 |
| C3 | effective-feature weight | 0/15 | 0.05 | 0/15 | — | 100% | 0.500 |
| C3 | CV-chosen $w_S$ | 0/15 | 0.40 | 12/15 | 1 (2), 10 (7), 100 (6) | 27% | 0.603 |
| C3 | two-step | 4/15 | 0.31 | 6/15 | — | 94% | 0.584 |

With the survival weight chosen by training CV, survival does move the
loadings toward the prognostic direction in the hardest cell: the program
becomes survival-active in most fits and validation C rises above both the
unweighted model and the two-step comparator. The loadings found are a
mixture, not the planted program (|r| 0.40), and the power likelihood makes
the β intervals far too narrow (27% coverage of nominal 90%). When the
program is easy to find (1% cells), CV keeps $w_S = 1$ or 10 and nothing
changes. Large weights also slow the fit: the column-wise message step
overshoots when many correlated features move together, so many steps are
damped or rejected (median 22 rejected steps per fit in C3).

### 4.3 Real data, expression only (TCGA + CPTAC → five external cohorts)

Benchmark preprocessing unchanged (per-cohort z-scoring, DeSurv gene
selection, 273 × 2,064, identical gene list to the frozen D4 fit; the
double log₂ on TCGA, Dijk and PACA-AU RNA-seq is kept so the comparison is on
identical inputs). Cox strata by cohort, $w_S = 1$, point-Laplace loadings.

- Training 5-fold CV over K_init ∈ {5, 7, 10, 15}: mean C 0.640, 0.652, 0.649,
  0.656 (SE about 0.04). The one-standard-error rule selected K_init = 5.
- Frozen fit: converged in 297 sweeps, no objective decreases. Two
  survival-active programs (β per SD of score 0.47, 95% interval 0.32 to
  0.62; −0.19, −0.35 to −0.04) and three reconstruction-only programs.
  Survival share of the loading precision about 5e-7.
- External cohorts, scored once (risk = η, direction fixed by the model):

| Cohort | n | ssbmf | D4 (frozen) | Two-step, same K |
|---|---|---|---|---|
| Dijk | 90 | 0.616 | 0.634 | 0.598 |
| Moffitt | 123 | 0.533 | 0.549 | 0.531 |
| PACA-AU array | 63 | 0.655 | 0.648 | 0.667 |
| PACA-AU RNA-seq | 52 | 0.692 | 0.657 | 0.701 |
| Puleo | 288 | 0.646 | 0.645 | 0.649 |
| Mean | 616 | **0.628** | 0.627 | 0.629 |

Donor-level paired bootstrap (2,000 replicates; PACA-AU donors with records
in both series resampled together): ssbmf − D4 mean difference +0.002 (95%
interval −0.021 to 0.024); ssbmf − two-step −0.001 (−0.010 to 0.009). No
cohort's interval for ssbmf − D4 excludes zero.

### 4.4 Real data, expression + methylation (TCGA → ICGC primary PDAC)

Matched TCGA training data (n = 144, 75 events; 3,000 genes on the log
scale, 10,000 CpG β-values; the cache as built, no scaling), $w_S = 1$,
point-Laplace loadings, K_init by 5-fold training CV over {3, 5, 7, 10};
ICGC primary PDAC (n = 50, 30 events) scored once per frozen fit.

| Fit | CV mean C by K_init (3, 5, 7, 10) | Selected K_init | Converged (sweeps) | CV folds converged | Survival-active programs | ICGC C |
|---|---|---|---|---|---|---|
| Joint (expression + methylation) | 0.566, 0.561, 0.548, 0.586 | 3 | yes (147) | 17/20 | 0 | 0.554 |
| Expression only | 0.557, 0.574, 0.616, 0.597 | 7 | yes (343) | 19/20 | 2 | 0.707 |

Paired donor-level bootstrap, joint − expression only: −0.153 (95% interval
−0.266 to −0.051). In the joint fit the shared Normal prior on β was fitted
to a variance of about 2e-7: at the converged loadings no program's score
carried survival signal above its noise level, so all three programs are
reconstruction-only and the joint risk score is essentially uninformative.
The first joint program explains 35% of the methylation variance and 8% of
the expression variance. The expression-only fit has two survival-active
programs (β per SD −0.36, 95% interval −0.59 to −0.13; +0.26, 0.06 to 0.45).

This repeats, with a correct objective, the earlier observation that
methylation dominates the joint factorization: with 10,000 CpGs against
3,000 genes and one Gaussian observation per entry, the programs follow the
methylation variance. The comparison is on 50 donors and the K_init choice
by the one-standard-error rule fell on the smallest value of the grid.

## 5. Provisional choices to review

1. **Cox treatment:** the Gaussian-moment Jensen bound (option (c) of the
   derivation) instead of the fixed-curvature bound (exact but inflates the
   survival precision by one to two orders of magnitude) or the working
   quadratic (not a bound).
2. **Non-conjugate updates:** gradient messages in the mean parameters with
   damping and an exact-objective check; up to 8 halvings before rejection;
   up to 25 inner iterations for β.
3. **Loading prior:** point-Laplace (signed) by default; point-exponential
   gave the same recovery and converged faster in simulation.
4. **β prior:** shared Normal, variance by marginal likelihood (no
   hyperprior, no lower bound on the inclusion probability).
5. **Program roles:** survival-active if lfsr ≤ 0.05; reconstruction-only if
   ≥ 0.1% of a block's centered variance.
6. **Numerics:** precision cap $\tau_j \le 10^6/\widehat{\mathrm{Var}}(y_j)$;
   tolerance 1e-6 relative over a sweep, 1,000 sweeps; `ebnm` results that
   violate their bounds treated as no information.
7. **Initialization:** SVD of the blocks scaled by their root mean square
   and concatenated; survival-score start available but off.
8. **Reporting scale:** β per SD of the training program score.
9. **Simulation design:** the cells, thresholds and the definition of
   "survival-active" for the two-step arm (Wald p ≤ 0.10, matching a
   one-sided 0.05 lfsr) in `results/ssbmf/sim_design.md`; cell C3 added
   after a pilot; the $w_S$ grid for the CV arm reduced to {1, 10, 100}
   after a pilot showed fits at $w_S \ge 1{,}000$ were slow and no better.
10. **Real data:** current benchmark preprocessing (including the double log);
    external genes missing from a cohort contribute nothing (filled with the
    training mean); multimodal inputs used as in the cache (log expression,
    methylation β-values, no scaling) with K_init ∈ {3, 5, 7, 10}.
11. **Balance:** the effective-feature rule is implemented but off and
    dropped on this evidence; the CV-chosen weight is available in
    `ssbmf_cv(w_grid = ...)` but not the default.

## 6. Tests and checks

- Package tests: 216 passing (`devtools::test("ssbmf")`), covering the Cox
  primitives against `survival::coxph` (ties and strata, to 1e-8), data
  validation, closed-form KL checks, monotonicity of each update on 20 random
  states, stationarity, agreement with `flashier` at $w_S = 0$ (fitted means
  within 0.34%), scale stability in p, prediction alignment and direction,
  roles, canonical scaling, the simulator, CV, and the balance rules.
- `R CMD check --no-manual`: Status OK.
- `verify_derivation.R`: bound below Monte Carlo in 24/24 states; stationarity
  gradients ≤ 7.7e-7.
- The old repository test suite (`tests/run_tests.R`) was not run: no file in
  `code/` or `tests/` was changed.

## 7. What this means and what comes next

The rebuild settles the algorithmic questions: with a bound in place of the
working quadratic, exact KL bookkeeping and per-step checks, the fit is
monotone and converges at realistic sizes, and real-data prediction is
unchanged. The supervision question is a modelling one. Counting every
matrix entry as an independent Gaussian observation gives reconstruction
about $10^3$ to $10^6$ times the precision the Cox term contributes to a
prognostic program's loadings, so at $w_S = 1$ the model is effectively
unsupervised in the loadings. The same counting lets 10,000 CpGs dominate
3,000 genes in the joint fit. The effective-feature
correction moves weight in the wrong direction for this purpose. A survival
weight chosen by cross-validation recovers prognostic signal but not a clean
program. Open questions for review: whether a CV-chosen $w_S$ is acceptable given
that it breaks the calibration of the β intervals; whether supervision
should act on a smaller, survival-screened feature set; whether recovery
should be judged on loadings or on the prognostic score; how the two blocks
should be balanced in the joint model (per-block weights by training CV, or
a methylation set reduced to the scale of the expression set); and whether
the column-wise message step should be replaced by feature-wise
(Gauss–Seidel) updates, which would remove most damped and rejected steps at
large weights.

Files: package `ssbmf/`; derivation `derivations/unified_ssbmf/`; runners
and committed summaries `results/ssbmf/` (fits are git-ignored under
`results/ssbmf/outputs/`).
