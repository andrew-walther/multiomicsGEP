# Codex prompt: finish task D (feature-group simulation and matched-data comparison)

Date: 2026-11-12. Author: Andrew Walther.

This continues Part 2 of the
[D/E follow-up prompt](2026-11-12-codex-d-e-followup-prompt.md). The fitter
(`feature_group` in `fit_ssbmf_yfb()`), its tests and the paused runner are committed in
`8a3b5517`, `6f31a0d7` and `ae5da859`
([findings](../findings/2026-11-12-feature-group-fitter-findings.md)). Andrew answered
the open recovery-target question on 2026-11-12. This prompt fixes the rest of the design.
The "Rules for both parts" section at the end of the follow-up prompt applies unchanged.

## What the study has to show

The reason for the multimodal model is a separate F prior for each data type. Expression
and methylation programs can differ in sparsity and scale, and one pooled prior has to
compromise between them. The step-3 simulation as prepared fits only the grouped model,
on data where both modalities draw loadings with the same sparsity (15% nonzero). It can
show that the grouped fitter recovers the prognostic program, but not that separate priors
help. Two additions do that: a pooled-prior arm on the same data, and cells where the two
modalities differ in sparsity.

## Decision 1: one recovery target, by regenerating survival (Andrew, 2026-11-12)

The generator writes Y = LF′ + noise on the raw scale and draws survival from η = YFβ.
After training-row z-scoring with SDs D, the reconstruction loading is D⁻¹F and the
prognostic direction is DF. These differ because the generator's feature precisions vary
(τ ~ Gamma(2, 2)). YFB has one F for both roles.

Remove the difference instead of choosing between them. This is the `truth_mode = "shared"`
correction that S-9 adopted (`yfb_decision_data()` in
`results/benchmark_sim/run_yfb_supervision_decision.R`) and that task E's S-9 and S-14
cells used, so D's gate means what A's and E's gates mean.

For each seed:

1. Generate the data with the generator unchanged except for the sparsity argument below.
2. z-score each modality on training rows, and scale validation with the training means
   and SDs (`fgsim_scale()`).
3. Truth direction: g = the prognostic program's standardized loading D⁻¹F[, 3], stacked
   across modalities in the fitted column order and scaled to unit norm. Keep the
   per-modality pieces for block-level reporting.
4. η = X_c g · a, where X_c is the centered standardized training matrix and a is a single
   constant per seed. Set a so that SD(η) on training rows equals SD of the original
   generator's η on the same training rows. This keeps the survival signal strength of
   the registered scenario.
5. Redraw training and validation event and censoring times from η with the generator's
   own Weibull settings and censoring calibration (`multimodal_yfb_calibrate_censoring()`),
   using a fixed seed offset. Validation uses the same g, a and training centers.

Recovery is scored against g only. Do not report raw-F or DF correlations.

## Decision 2: the arms and cells

All fits use the working configuration: K = 5, point-Laplace L and F, per-program
Normal β, w_S = 1, `canonicalize = TRUE`, pooled, no strata, no pruning, tolerance and
iteration limit from `config/globals.yml`. Seeds 1–10 in every cell.

| Factor | Levels |
|---|---|
| Feature ratio, expression : methylation | 1:1 (50 : 50) and 1:5 (50 : 250) |
| Methylation loading sparsity | same as expression (15% nonzero), or different (3% nonzero) |
| Prior arm | grouped (`feature_group` by modality) and pooled (`feature_group = NULL`) |

That is 4 cells × 2 arms × 10 seeds = 80 fits. Both arms in a seed get the identical
scaled matrix, survival draw and initialization. Assert this as task E did.

To set methylation sparsity, add an argument `methylation_active_fraction = NULL` to
`simulate_multimodal_yfb_data()` and pass it to `simulate_multimodal_yfb_loadings()` for the
methylation block only. NULL keeps the current behavior. Add a test that NULL reproduces
the current output exactly for a fixed seed and that 0.03 gives the requested number of
nonzero methylation loadings per program.

## Decision 3: recovery, false activity and the adoption rule

Bring `fgsim_recovery()` in line with S-9's scoring (`yfb_decision_score()` in the S-9
runner):

- matched program = the fitted F column with largest |cor(g, F_k)| over all features;
- recovered = |cor| ≥ 0.7, survival-active, and E[β_k]·cor(g, F_k) > 0 (correct sign);
- survival-active = |E[β_k]| ≥ 1.96 SD(β_k), with SD(β_k) from `fit$VBeta` when present
  (AGENTS.md: use the variance carriers, not EBeta2 − EBeta²);
- false active = survival-active programs other than a recovered match.

Per cell and arm, report recovery successes out of 10, mean false-active programs,
converged out of 10 and mean validation C (with the defined count). The gate is the task D
criterion: recovery ≥ 8/10, mean false activity ≤ 1, convergence ≥ 9/10. Failed seeds stay
in the denominator.

Adoption rule for the real-data fit, registered here before any fit: the grouped prior is
used in step 4 unless pooled passes the gate in some cell where grouped fails. If grouped
fails the gate in every cell, stop and report before step 4.

Descriptive reporting, not part of the rule:

- the per-modality block correlations cor(g_m, F_m,k) for the matched program;
- each group's fitted F prior (point-mass weight and Laplace scale) in the grouped arm,
  against the true nonzero fraction of that modality, and the pooled prior;
- each modality's share of the L precision, Σ_j τ_mj E[F²_mjk] over the modality's
  features, for every live program. Methylation dominance means its share exceeds its
  feature share of the total in the matched program. Report it. Do not add modality
  weights; they need their own derivation first.

Save each warning message individually, as the prepared runner already does. Patient-level
or per-seed fits stay under ignored `data/`; commit only aggregate CSVs and `settings.txt`.

## Step 4: matched TCGA patients, descriptive

Run only if the step-3 gate passes for the arm the rule selects. Use the 144 matched TCGA
patients and the five folds from task C, and task C's 10,000-CpG point-Laplace setup
(CpG z-scores) with D4 genes for expression, both z-scored on training folds only. K = 5,
working configuration. Four fits per fold:

1. expression alone,
2. methylation alone,
3. both, grouped prior,
4. both, pooled prior.

Report per-fold and mean within-cohort C (one cohort here, so this is the TCGA C), the
paired differences of fits 3 and 4 against fit 1 with SE = SD of the five differences / √5,
each modality's L-precision share in fits 3 and 4, and the programs in the full-data
joint fit whose top-20 CpGs or genes do not match any program in the single-modality
full-data fits (|loading correlation| < 0.5 on the shared features). Nothing in step 4 is
adopted. It answers the task D question of whether the joint fit adds prognostic or
interpretive information.

## Order and commits

1. Tests: the sparsity argument, the survival regeneration (g is unit norm, SD(η) matches,
   validation uses training statistics), the identical-input assertion and the S-9 sign
   check in `fgsim_recovery()`. Commit.
2. Runner changes, smoke run on one seed per cell, full suite. Commit.
3. Full simulation (80 fits, at most 2 workers). Commit aggregates.
4. Step 4 if the rule allows. Commit aggregates.
5. Update `docs/findings/2026-11-12-feature-group-fitter-findings.md` with sections for the
   simulation and step 4 and a revised "Text for the shared status files". Independent
   review of the code and findings, then the final commit and push.

**Done when.** Step 3 is reported against the gate and the adoption rule for every cell
and arm; step 4 is reported or marked not run with the reason; the full suite passes; the
independent review is recorded in the findings note.
