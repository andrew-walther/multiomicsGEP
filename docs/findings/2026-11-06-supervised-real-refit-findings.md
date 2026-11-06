# Attempt S-15: real-data refit with the CV-chosen survival weight — findings

Registration: [`docs/plans/2026-11-06-supervised-real-refit-plan.md`](../plans/2026-11-06-supervised-real-refit-plan.md).
Runner: `results/benchmark_sim/run_yfb_supervised_real_refit.R`, with the fitter and predictor
sourced from commit `5ef54b3`. Output: `results/benchmark_sim/outputs/yfb_supervised_refit_2026-11-06/`.
There were 216 of 216 fits and no errors. Training patients only (TCGA_PAAD and CPTAC, n = 273,
2,218 D4 genes).

## Verdict

**w* = 1, so nothing is adopted.** The S-11 fit (w_S = 1) and its S-12 external scoring stand,
and no new external scoring is run. The reproducibility check passes: the w = 1 fold fits match
S-11 to within 3 × 10⁻¹⁶ in C.

## Weight curve on the S-11 folds

| w_S | CV partial log-likelihood | Mean CV C (SE) | Folds converged | Survival share of loading precision, active programs (median) |
|---|---|---|---|---|
| 1 | −775.168 | 0.663 (0.024) | 3 of 5 | 1.2 × 10⁻⁶ |
| 2 | −775.179 | 0.663 (0.025) | 3 of 5 | 2.3 × 10⁻⁶ |
| 5 | −775.213 | 0.661 (0.025) | 4 of 5 | 5.9 × 10⁻⁶ |
| 10 | −775.244 | 0.661 (0.025) | 4 of 5 | 1.2 × 10⁻⁵ |
| 20 | −775.380 | 0.659 (0.026) | 4 of 5 | 2.5 × 10⁻⁵ |
| 50 | −775.429 | 0.653 (0.028) | 4 of 5 | 7.3 × 10⁻⁵ |

Folds that did not converge reached the 300-sweep cap. Every fit ran without error.

**Nested CV.** The weights chosen inside the five outer training sets were 5, 2, 50, 1 and 10.
The nested mean C was 0.656 (SE 0.027), against 0.663 for w_S = 1 and 0.652 for the two-step
fit. It is above the registered threshold of 0.639.

## What this shows

1. **On the real data the weight does not matter in this range.** The CV partial
   log-likelihood changes by 0.26 across a fifty-fold range of weights, out of about 775. That
   is why the weights chosen in the nested folds jump from 1 to 50.
2. **Survival still cannot move the real-data loadings.** Even at w_S = 50, survival supplies
   about 7 × 10⁻⁵ of the loading precision of the active programs. Every program's loadings
   correlate 1.000 with the unsupervised (w_S = 0) loadings. The share grows linearly with
   w_S, so reaching a share near 0.1 would take a weight in the tens of thousands. That is the
   same order as DeSurv's normalized weighting (about 4,000 on our scale). In the S-9
   simulation, weights of 40 or more stopped converging, so a larger weight would need the
   stability work first.
3. **The supervised fit's real-data programs are a small, clear set.** At w_S = 1 three programs
   are active, with standardized coefficients of 0.41, 0.22 and 0.23 (log hazard ratio per SD of
   the program score). Their ratios to the largest are 1, 0.54 and 0.57. Each carries 21% to 71%
   of the risk-score variance. All three pass the proposed ratio rule (≥ 0.25 of the largest
   and ≥ 5% of the risk variance). The fold-stability part of that rule has not yet been
   computed. One inactive program has a standardized coefficient of 0.15 and 10% of the risk
   variance. It is the borderline case for that rule.
4. **The CV choice behaves as it did in simulation.** It picks the smallest weight when the data
   do not support more. On these data it gives no reason to move away from w_S = 1.

The full-fit program table is `programs_full_fit.csv`, and the weight curve is
`weight_curve.csv`.

## Limits

- Five folds of 273 patients, so CV C differences below about 0.03 are within noise.
- The weight grid stops at 50. The precision-share result says a much larger weight would be
  needed to change the loadings, but this grid cannot test whether that would help.
