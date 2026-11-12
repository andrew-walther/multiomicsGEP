# Codex prompt: task D next steps (matched-TCGA comparison and a recoverable scenario)

Date: 2026-11-12. Author: Andrew Walther.

This continues [the D rerun prompt](2026-11-12-codex-d-rerun-prompt.md). The rerun
(`041a3e5d` to `b56cbeac`, [findings](../findings/2026-11-12-feature-group-fitter-findings.md))
fixed the grouped crash: all 120 fits converged. No arm passed recovery in any cell,
including expression alone, while the true-risk C was about 0.78. The
`low_variance_prognostic_factor` scenario at 120 training patients is too weak for the gate
to separate the arms. The "Rules for both parts" section at the end of
[the D/E follow-up prompt](2026-11-12-codex-d-e-followup-prompt.md) applies unchanged.

Andrew decided on 2026-11-12 to do two things at once (recorded in DECISIONS.md). They are
independent. Use one worker for each, at most 2 in total. Do Part A first if you must
choose, because the prelim deadline depends on it.

## Part A: matched-TCGA comparison, descriptive, without a passing simulation gate

Andrew overrides the registered stop. The rule assumed the simulation could separate
the arms, and the rerun showed it cannot, since the single-modality model also fails.
This comparison is descriptive. Label it in every table and in the findings as **not
validated by simulation**. Nothing in it is adopted.

The design is step 4 of [the finishing prompt](2026-11-12-codex-d-finish-prompt.md),
unchanged:

- the 144 matched TCGA patients and the five folds from task C;
- task C's 10,000-CpG point-Laplace setup (CpG z-scores) and D4 genes for expression,
  each z-scored on training folds only;
- K = 5 and the working configuration (point-Laplace L and F, per-program Normal β,
  w_S = 1, `canonicalize = TRUE`, no strata, no pruning, limits from `config/globals.yml`);
- four fits per fold and on the full data: (1) expression alone, (2) methylation alone,
  (3) both with `feature_group` by modality (the multimodal model; assert that
  `fit$prior_F` holds K expression and K methylation priors), (4) both with a pooled prior.

Reuse task C's prepared inputs and fold labels. Read them from task C's saved preparation,
and do not edit `results/multimodal_real/run_methylation_single_modality.R` or its outputs.
Write a new runner, for example `results/multimodal_real/run_multimodal_tcga_comparison.R`.

Report:

- per fold and on average, the TCGA C for each fit, and the paired differences of fits 3
  and 4 against fit 1 with SE = SD of the five differences / √5;
- convergence per fit, with any iteration-cap fits named. Task C's expression fits hit
  the cap in most folds. Report caps, do not rerun to avoid them;
- each modality's L-precision share for every live program in fits 3 and 4;
- the empty-block counts in fit 3, by modality and program;
- in the full-data fit 3, each program's fitted point-mass weight and Laplace scale for
  each modality, next to the pooled prior from fit 4;
- programs in the full-data fits 3 and 4 that match no program in fits 1 and 2
  (|loading correlation| < 0.5 on the shared features), with their top-20 genes and CpGs.

Patient-level inputs, fits and scores stay under ignored `data/`. Commit only aggregates.

## Part B: simulation rerun on a scenario the single-modality model can recover

Same runner, arms, cells, seeds and configuration as the rerun, with one change: the
scenario is `both_informative` instead of `low_variance_prognostic_factor`. In
`both_informative`, programs 1 and 2 are prognostic (β = 0.70 and 0.35) with the usual
activity scale, and program 3 is not prognostic. Keep 120 training patients, close to the
144 matched TCGA patients.

- **Survival regeneration.** As before, but with two prognostic directions. For k = 1, 2,
  g_k is the standardized loading D⁻¹F[, k] stacked across modalities in the fitted column
  order and scaled to unit norm. η = a · X_c (g_1 β_1 + g_2 β_2), with a set per seed so
  that SD(η) on training rows equals the SD of the original generator's η. Validation
  uses the same g, a and training centers.
- **Recovery (gate).** Scored on program 1, the stronger one, with the existing rule
  (|cor(g_1, F_k)| ≥ 0.7, survival-active, correct sign).
- **Program 2.** Report its recovery the same way, but descriptively.
- **False activity.** Survival-active programs other than the fitted programs matched to
  true programs 1 and 2.
- **Expression-only arm.** Scored against the expression parts of g_1 and g_2, each
  renormalized to unit norm.
- **Other columns.** True-risk C, empty-block counts and precision shares, as in the rerun.

4 cells × 3 arms × 10 seeds = 120 fits, in new output and private directories (for
example `feature_group_simulation_2026-11-12_informative`). Do not touch earlier runs.

**Registered rules, set before any fit.**

- Gate unchanged: recovery ≥ 8/10, mean false activity ≤ 1, convergence ≥ 9/10.
- The scenario counts as able to separate the arms only in cells where expression-only
  passes the gate. Compare grouped and pooled only in those cells. If grouped passes
  where pooled fails, or pooled fails false activity where grouped passes, report that
  separate priors helped. If pooled passes where grouped fails, report that they hurt.
- If expression-only fails in every cell, report that this sample size cannot test the
  priors. Do not change the scenario or the sample size in this task.
- Part B's outcome does not change Part A, which stays descriptive.

## Order and commits

1. Tests: the two-direction survival regeneration (unit-norm g_k, SD(η) matches,
   validation uses training statistics) and the two-program scoring. Part A's runner
   checks: fold labels match task C's, scaling uses training folds only, the `prior_F`
   assertion. Commit.
2. Runner changes and a smoke run (one fold for A, one seed per cell for B). Full suite.
   Commit.
3. Part A full run. Commit aggregates. Part B full run, concurrently if possible. Commit
   aggregates.
4. Findings. Part A goes in a new note,
   `docs/findings/2026-11-12-multimodal-tcga-comparison-findings.md`. Part B goes in a new
   dated section of the feature-group findings note. Both get a "Text for the shared
   status files" section and rows in `docs/findings/README.md`. Independent review of the
   code and findings, then the final commit and push.

**Done when.** Part A is reported with the not-validated label. Part B is reported against
the gate and the separation rule. The full suite passes, and the review is recorded.
