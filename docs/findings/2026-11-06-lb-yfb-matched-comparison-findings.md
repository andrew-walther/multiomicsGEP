# LB and YFB at matched survival weights (2026-11-06)

Exploratory comparison, not a registered study. It changes neither model's
configuration.

## Question

How do the two parameterizations compare on the same data when survival gets
the same weight against genomics, and how does each compare with its own
two-step version (unsupervised fit, then Cox)?

## Design

- **Data.** Planted programs: Y = LFᵀ + noise and η = Lβ, with β = (1, −0.8, 0.5, 0).
  - n = 200 training and 100 test, p = 400, K = 4.
  - Each program loads on 10% of genes, at amplitude 0.5 or 1.0.
  - Six seeds per amplitude, so 12 datasets in all.
  - This is the design used in the 2026-11-06 LB review.
- **Matched weights.** In the LB L update the survival-to-genomics precision ratio
  is w = α/(1 − α). YFB with the survival-block weight multiplies the Cox term
  by w_S. So LB at α = w/(1 + w) and YFB at w_S = w give survival the same weight
  against genomics in the objective. The comparison uses w = 1 (LB α = 0.5) and
  w = 10 (LB α = 0.91). Equal weight is not equal influence:
  - LB also tempers its L and F updates by 1 − α. At w = 10 their noise
    variances are 11 times larger.
  - LB fits β with Cox weight α (0.5 or 0.91). YFB fits β with weight 1.
  - The survival precision is per patient in LB and per gene in YFB, on
    different scales.
- **Configurations.** Each model ran as it was designed to be used.
  - LB: raw nonnegative Y with point-exponential priors.
  - YFB: training-standardized Y, point-Laplace priors, score standardization
    and the survival-block weight. This is the S-11/S-13 configuration.
  - Two-step references: LB at α = 0 and YFB at w_S = 0, each followed by a Cox
    fit on its training scores.
- **Test scores.** LB uses fold-in loadings (`predict_ssbmf_lb_foldin()`). YFB
  uses its program scores Y\*F.
- **Recovery.** The mean |r| between the true test loadings of the three
  survival programs and their estimated scores, matched one to one. A program
  matched to a dead column scores 0.
- **Oracle.** The test C of the true risk score is 0.790 at both amplitudes.

Runner: `results/benchmark_sim/run_lb_yfb_matched_comparison.R`. Outputs:
`results/benchmark_sim/outputs/lb_yfb_matched_2026-11-06/`.

## Results

The table gives means over six seeds. All 72 fits converged, and so did every
LB fold-in. Nine LB fits, from three datasets at every weight, produced one
warning each. It came from the warm-start Cox fit on the initial SVD loadings
("coefficient may be infinite"). No YFB fit warned.

| Amplitude | Model | Weight | Live programs | Test recovery | Train C | Test C (range) |
|---|---|---|---|---|---|---|
| 0.5 | LB | two-step | 3.0 | 0.664 | 0.712 | 0.710 (0.54–0.79) |
| 0.5 | LB | w = 1 | 2.8 | 0.828 | 0.719 | 0.724 (0.57–0.80) |
| 0.5 | LB | w = 10 | 1.8 | 0.366 | 0.720 | 0.598 (0.50–0.70) |
| 0.5 | YFB | two-step | 4.0 | 0.843 | 0.762 | 0.782 (0.76–0.80) |
| 0.5 | YFB | w_S = 1 | 4.0 | 0.843 | 0.763 | 0.783 (0.76–0.80) |
| 0.5 | YFB | w_S = 10 | 4.0 | 0.848 | 0.771 | 0.780 (0.75–0.80) |
| 1.0 | LB | two-step | 3.0 | 0.772 | 0.735 | 0.739 (0.63–0.80) |
| 1.0 | LB | w = 1 | 3.0 | 0.769 | 0.742 | 0.757 (0.72–0.80) |
| 1.0 | LB | w = 10 | 3.0 | 0.705 | 0.757 | 0.727 (0.65–0.79) |
| 1.0 | YFB | two-step | 4.0 | 0.793 | 0.766 | 0.787 (0.76–0.81) |
| 1.0 | YFB | w_S = 1 | 4.0 | 0.809 | 0.766 | 0.787 (0.76–0.80) |
| 1.0 | YFB | w_S = 10 | 4.0 | 0.809 | 0.769 | 0.788 (0.76–0.81) |

"Live programs" counts columns of F with any nonzero gene weight.

## What this shows

1. **YFB predicts better and more consistently.**
   - YFB's test C is 0.78–0.79, close to the oracle's 0.790, and varies by
     about 0.05 across seeds.
   - LB's test C at w = 1 is 0.72–0.76 and ranges from 0.57 to 0.80 across
     seeds.
2. **Supervision helps LB at w = 1, mostly at the weaker signal.**
   - At the weaker signal, LB's test C rises over its two-step version by 0.015
     on average (five of six seeds higher). Recovery rises from 0.664 to 0.828.
   - At the stronger signal, test C rises by 0.018 on average. One seed accounts
     for nearly all of it (+0.121); the other five change by less than 0.012.
     Recovery is unchanged.
3. **LB at w = 10 fails at the weaker signal.**
   - Only 1.8 of 4 programs remain live on average.
   - Two of six fits shrank every β to zero, giving test C 0.500.
   - One fit separated the training data (train C 1.00, test C 0.51).
   - At the stronger signal, w = 10 keeps three programs, with test C 0.727
     against 0.739 for the two-step fit.

   This repeats the LB collapse at strong supervision recorded in DECISIONS.md
   (2026-07-12). Here the cause cannot be separated from the 11-fold tempering of
   the L and F noise that comes with α = 0.91.
4. **YFB's predictions are stable across weights, but its recovered programs
   are not stable for every seed.**
   - Mean test C and mean recovery at w_S = 1 and 10 match the two-step fit.
   - The one mean recovery gain (0.016 at the stronger signal) comes from a
     single seed, which gained 0.095.
   - At the weaker signal and w_S = 10, the per-seed recovery changes from the
     two-step fit range from −0.137 to +0.169, which cancel in the mean.
   - Survival carries little weight in the F update here because these planted
     programs explain a large share of the variance, as the precision-share
     analysis predicts. S-14 found its gains near the detection threshold, not in
     data like these.
5. **LB loses a planted program in every fit.** Every LB fit, two-step or
   supervised, ended with at least one of the four planted programs dead. This
   matches the dead programs from the pmax(SVD) start seen in the LB review. YFB
   kept all four in every fit.

## Limits

- **Small synthetic sample.** These are 12 small synthetic datasets with LB-style
  planted structure. YFB scores estimate L(FᵀF), so the design does not favor
  YFB.
- **Methods as used, not parameterizations alone.** The models differ in
  preprocessing, priors and the tempering described above. The comparison is
  between the two methods as used, not between the parameterizations with
  everything else held fixed.
- **Recovery measures different quantities.** Recovery compares fold-in loadings
  for LB with program scores for YFB.
- **Train–test gaps are not comparable across models.** LB's training C uses its
  in-sample, outcome-informed loadings, while its test C uses fold-in.

## Decision

None. The results support the current plan:

- YFB stays the primary model.
- LB stays an alternative with fold-in and is not run at strong supervision.
