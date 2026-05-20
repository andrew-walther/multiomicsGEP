# YFB Fit Improvements Plan — 2026-05-08

## Summary

This plan defines a later implementation pass to make YFB fitting, evaluation,
and K cross-validation internally consistent and interpretable. The work should
be staged incrementally so each change can be validated in isolation before it
is folded into the benchmark and default workflow.

The target outcomes are:

- single-cohort and merged YFB fits use aligned preprocessing logic and produce
  interpretable prognostic signals
- risk scores have one project-wide orientation convention
- YFB K-CV evaluates the same model definition that will be used in the final
  full-data fit
- K selection for YFB becomes interpretable enough to compare against LB and
  against shrinkage-based `K_eff` diagnostics

## Key Changes

### Phase 1 — Define and enforce one risk-score convention

- Make `risk_score` mean one thing everywhere: larger score = larger hazard =
  worse prognosis.
- Remove mixed evaluation conventions across benchmark runners and CV helpers;
  stop relying on some paths using `risk_score` and others using
  `I(-risk_score)`.
- Treat orientation as part of the fitted model, not an evaluation afterthought.
- Replace the current implicit beta-negation behavior with an explicit stored
  orientation field on fitted objects, such as `risk_sign` or `orientation`,
  learned at fit time and applied in prediction.
- Define the orientation rule once and reuse it for LB and YFB. Preferred rule:
  choose the orientation that maximizes the training Cox objective or yields
  concordance >= 0.5 on the training set, then store that orientation in the
  fit object.
- Update all prediction and evaluation helpers so full-data runs, external
  validation, synthetic holdout, and CV folds all consume the same oriented
  predictor.

### Phase 2 — Align YFB preprocessing across single-cohort and merged modes

- Refactor YFB preprocessing so single-cohort and merged training both flow
  through one explicit preprocessing policy layer rather than two unrelated code
  paths.
- Separate preprocessing options into named switches that can be logged and
  compared:
  - log transform
  - top-gene selection
  - per-platform standardization
  - quantile normalization
  - rank transform on or off
  - centering and scaling outputs used at fit and predict time
- Add a YFB-specific single-cohort preprocessing mode that can disable rank
  transform, matching the main lesson from the successful merged YFB
  configuration.
- Preserve the successful merged YFB path as a pinned reference configuration
  while testing single-cohort variants.
- Ensure external-cohort preprocessing mirrors the fitted training contract for
  each YFB mode, so prediction is not mixing incompatible transforms.

### Phase 3 — Make YFB CV evaluate the intended production model

- Update `select_K_cv()` so the YFB branch uses the same preprocessing,
  projection scaling, prior choice, and stored orientation logic as the
  benchmark fit.
- Fold orientation into each training fold as part of the learned model, then
  apply that fold-specific orientation to its held-out predictions.
- Add an explicit YFB CV configuration object or argument set so CV cannot
  silently drift from the benchmark runner.
- Revisit the YFB K grid after preprocessing is stabilized; include a smaller-K
  regime that reflects the current merged YFB result (`K=3`) and plausible
  single-cohort regimes.
- Preserve `K_eff` and shrinkage summaries as diagnostics, but do not use them
  as the sole YFB selection rule once CV is repaired.

### Phase 4 — Reframe K selection as “CV primary, shrinkage secondary”

- For YFB, adopt held-out prognostic performance as the primary K-selection
  criterion once CV is repaired.
- Keep ARD and shrinkage-based `K_eff` as secondary diagnostics used to
  interpret whether the chosen K is over-complete or under-activated.
- Define the later comparison workflow explicitly:
  - run larger `K_max` fits and inspect `K_eff`, beta magnitudes, and factor
    stability
  - run K-CV on a narrowed candidate grid
  - select the smallest K with stable held-out performance
- For documentation continuity, record that current YFB `K_eff` is useful for
  diagnosis but not yet sufficient as a final rank-selection rule on real PDAC
  data.

### Phase 5 — Add instrumentation before final benchmark reruns

- Add structured fit metadata to benchmark outputs for YFB:
  - preprocessing variant name
  - rank-transform flag
  - per-platform-standardization flag
  - orientation chosen
  - prior on beta
  - K, `K_eff`, `beta_max`
- Add intermediate diagnostics that explain why a run succeeds or fails:
  - distribution and scale of `ZF`
  - fraction of near-zero beta coefficients
  - training concordance before and after orientation assignment
  - fold-level CV concordance and orientation decisions
- Use these diagnostics to compare:
  - current single-cohort baseline
  - single-cohort no-rank variants
  - merged reference configuration
  - normal versus point-normal priors

## Test Plan

- Unit tests for fit and predict orientation:
  - fitted YFB object stores an explicit orientation field
  - prediction applies the stored orientation consistently
  - CV folds and full-data fits report concordance under the same sign
    convention
- Unit tests for preprocessing contracts:
  - single-cohort and merged YFB preprocessing return consistently structured
    outputs
  - disabling rank transform works in single-cohort mode
  - external preprocessing preserves training-gene ordering and compatible
    scaling
- Unit tests for CV consistency:
  - `select_K_cv(model="YFB")` uses the same preprocessing and orientation path
    as benchmark fitting
  - YFB CV output reports fold-level metadata needed to diagnose sign and
    collapse behavior
- Comparative validation runs:
  - single-cohort YFB with current preprocessing versus no-rank variant
  - merged YFB reference configuration remains reproducible
  - YFB normal-prior K-CV no longer produces systematically anti-concordant
    curves purely from sign convention mismatch
  - point-normal collapse, if still present, is measurable as a real shrinkage
    phenomenon rather than a preprocessing or orientation artifact
- Acceptance criteria:
  - one risk-score convention across code paths
  - YFB single-cohort runs have an interpretable preprocessing story and no
    hidden mismatch with merged mode
  - YFB K-CV curve can be interpreted directly without mentally converting `C`
    to `1-C`
  - final YFB K recommendation can be defended as a held-out result, with
    `K_eff` used as supporting evidence

## Assumptions

- This plan is implementation-focused and intentionally excludes advisor-facing
  report revisions, except where logs and metadata must be added to support
  later interpretation.
- Work should be executed incrementally, with each phase benchmarked before the
  next phase changes defaults.
- The successful merged YFB configuration remains the temporary reference
  baseline: per-platform standardization, no rank transform, and small-K
  fitting.
- The intended saved draft path for later execution is
  `docs/plans/YFB-fit-improvements-05_08_26.md`.
