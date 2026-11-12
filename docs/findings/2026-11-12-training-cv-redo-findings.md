# Task F findings: gene screen, input scaling and rank chosen by training CV alone

Date: 2026-11-12. Author: Andrew Walther.

**Result.** Using only the 273 TCGA and CPTAC training patients, the D4 gene screen and
per-cohort z-score input are kept, and the rank changes from K = 7 to **K = 5**. No
alternative screen or scaling beats the current one, and both lose by 0.012 to 0.020 in
within-cohort CV concordance. For the rank, CV concordance is flat from K = 5 to 12
(0.649 to 0.666) and clearly lower at K = 2 to 4 (0.611 to 0.620). The registered one-SE
rule therefore picks K = 5. K = 5 and K = 7 differ by −0.003 (paired SE 0.007), and K = 5
converges in all five folds where K = 7 converges in three. The external results that once
supported keeping K = 7 are no longer used for selection. Task A's prior choice holds at
K = 5: point-exponential on X/sd reaches 0.519 against 0.655 for point-Laplace.

Registration: the 2026-11-12 "F" entry in the
[attempts log](../plans/2026-10-16-round2-attempts-log.md), task F of the
[week plan](../plans/2026-11-09-week-plan.md). Runner:
`results/benchmark_sim/run_training_cv_redo.R`. Aggregate outputs:
`results/benchmark_sim/outputs/training_cv_redo_2026-11-12/`.

## Design

The settled model from tasks A and B: point-Laplace L and F, per-program Normal β,
`fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)`, no pruning, pooled training without strata.
Every arm uses the S-11 five folds (balanced by cohort and event), so each is paired with
task A's point-Laplace fits, which serve as the D4, z-score, K = 7 reference. The primary
metric is within-cohort concordance (pairs only within TCGA or within CPTAC). External
cohorts were not used.

- **Gene screen.** The 2026-05-27 comparison D1 to D5 crossed two screens with two models
  and a cohort indicator. Under YFB without a cohort indicator, two screens remain. D2: per-cohort
  z-score, quantile normalization of the merged matrix, top 2,000 genes by variance after
  normalization. D4 (current): per-cohort combined mean and variance rank, top 3,000 per
  cohort, intersection (2,218 genes), per-cohort z-score. The two sets share 258 genes.
  D1 and D3 are the LB model (task E), and D5 adds a cohort indicator, which the cohort
  decision replaced. Both screens are outcome-free and computed once on all 273 patients.
- **Input scaling.** On D4, under point-Laplace: per-cohort z-score (current), X/sd within
  cohort with a centered projection, and within-patient ranks followed by per-cohort
  z-score. The rank-transform point-exponential arm is reported in task A.
- **Rank.** K = 2 to 12 on D4 with the z-score input. `select_K_cv()` was not used. It
  scores pooled concordance, would need `strata_id` in the fit to score within-cohort
  pairs, and returns only the chosen K.

**Rules, fixed before fitting.** For the screen and the scaling, the current choice is kept
unless an alternative is better by more than one SE of the paired fold difference. For K,
the one-SE rule: the smallest K whose mean within-cohort C is at least the best mean minus
the SE of the fold concordances at the best K.

All 78 new fits returned without errors. With the six point-exponential fits of the
task A rerun, that makes 84.

## Gene screen and input scaling (K = 7)

| Input | Genes | Mean within-cohort C (SE) | Difference from current (paired SE) | Folds converged | Survival-active per fold (mean) |
|---|---|---|---|---|---|
| D4, per-cohort z-score (current) | 2,218 | 0.658 (0.026) | | 3/5 | 2.6 |
| D2 screen | 2,000 | 0.646 (0.030) | −0.012 (0.014) | 5/5 | 1.4 |
| D4, X/sd, centered projection | 2,218 | 0.638 (0.023) | −0.020 (0.007) | 5/5 | 3.0 |
| D4, rank then z-score | 2,218 | 0.641 (0.022) | −0.017 (0.009) | 5/5 | 2.8 |

None of the alternatives is better than the current input, so D4 and the per-cohort
z-score stay. The point-Laplace prior handles the X/sd input reasonably (0.638), unlike
point-exponential in task A (0.573), because its signed programs need not encode the
cohort mean. The alternatives converge in more folds than the current input.

## Rank

| K | Mean within-cohort C (SE) | Difference from K = 7 (paired SE) | Folds converged | Survival-active per fold (mean) | Survival-active, full fit |
|---|---|---|---|---|---|
| 2 | 0.611 (0.022) | −0.047 (0.009) | 5/5 | 1.0 | 1 |
| 3 | 0.617 (0.019) | −0.040 (0.010) | 5/5 | 1.0 | 1 |
| 4 | 0.620 (0.024) | −0.038 (0.009) | 5/5 | 1.6 | 2 |
| **5** | **0.655 (0.029)** | −0.003 (0.007) | 5/5 | 2.6 | 3 |
| 6 | 0.659 (0.022) | +0.001 (0.005) | 4/5 | 2.4 | 3 |
| 7 | 0.658 (0.026) | | 3/5 | 2.6 | 3 |
| 8 | 0.656 (0.027) | −0.001 (0.003) | 3/5 | 3.0 | 3 |
| 9 | 0.649 (0.022) | −0.009 (0.005) | 2/5 | 2.8 | 4 |
| 10 | 0.666 (0.027) | +0.008 (0.004) | 2/5 | 3.2 | 4 |
| 11 | 0.661 (0.021) | +0.003 (0.007) | 1/5 | 3.4 | 4 |
| 12 | 0.654 (0.020) | −0.004 (0.007) | 0/5 | 3.4 | 4 |

The best mean is at K = 10 (0.666, SE 0.027), so the threshold is 0.639 and the smallest K
above it is 5. Under point-Laplace every column keeps nonzero loadings at every K, so the
number of programs with nonzero loadings equals K. The empirical-Bayes priors do not shrink
whole columns to zero here. The number of survival-active programs rises slowly with K, from
3 at K = 5 to 7 to 4 at K = 9 to 12. Convergence within 300 sweeps worsens as K grows. The
held-out partial log-likelihood follows the same shape (−155.9 at K = 2 and 3, −154.7 at
K = 5, −155.0 at K = 7, −153.9 at K = 11). No program at any K correlates above 0.09 with
the cohort label.

**How decisive is K = 5?** The one-SE rule uses the SE of the fold concordances (0.027),
which includes fold-to-fold variation shared by every K. Paired differences are much more
precise (SE about 0.005), and on that scale K = 5 to 8 are indistinguishable while K = 2
to 4 are clearly worse. K = 10 is 0.008 above K = 7 (paired SE 0.004). DeSurv used the
same unpaired one-SE rule. The rule picks the smallest rank on the plateau. It does not
show that K = 5 predicts better than K = 7.

## Task A rerun at K = 5

| Arm | Mean within-cohort C | Difference (SE) | Folds converged | Survival-active per fold | Max \|cor\| with cohort |
|---|---|---|---|---|---|
| point-Laplace, z-score | 0.655 | | 5/5 | 2.6 | 0.05 |
| point-exponential, X/sd | 0.519 | −0.135 (0.021) | 0/5 | 0 | 0.998 |

Task A's criterion 4 fails again, so point-Laplace stays at the new rank.

## Consequences

- K = 5 is adopted for new work under the registered rule. The external C screens of
  July and August that supported K = 7 are no longer used for selection.
- Results already reported at K = 7 (S-11 to S-15, the S-12 external scores, and the
  task A, A2, B, C and E work this week) stay on record and are labeled with K = 7.
- `config/globals.yml` (`benchmark$k_merged_yfb_desurv: 7`) is unchanged here. Many
  historical runners read that key, so changing it would alter their reruns. Whether to
  add a new key or change this one is a decision for Andrew.
- Prelim Chapter 4 reports the fixed K = 7 fit, and its numbers stay. K = 5 is described as
  the training-CV choice for the next fit, as the week plan specifies.

## Limits

- One fold split (seed 2026). The plan's repetition with four more seeds when comparisons
  are close was not run. The K = 5 to 12 plateau is exactly that case.
- K = 7 to 12 do not converge in every fold within 300 sweeps. Their concordance is taken
  from the state at the cap.
- The screen comparison covers the two screens left after removing model and cohort-indicator
  variants. The gene-set sizes and survival screens compared in S-11 were not redone.

## Text for the shared status files

Recorded by this session in DECISIONS.md, PROJECT_STATUS.qmd, the attempts log and the
week-plan checklist.
