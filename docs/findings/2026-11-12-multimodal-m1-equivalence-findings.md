# M = 1 multimodal equivalence and KL checks

Author: Andrew Walther  
Date: 2026-11-12

## Result and scope

The registered M = 1 acceptance criterion is **not met**. The pooled and stratified tests both stop because `fit_multimodal_yfb()` has no `strata_id` argument. The multimodal fitter also does not implement `center_projection` or the single-modality `w_S` control. Passing those names in `control` does not establish their use. The test must check that each requested option is implemented when the interface is extended. No two-modality simulation or real-data comparison was run.

The test was written and run before changes to the update module. Its baseline was 0/4 checks: the two equivalence fixtures stopped at the unsupported argument, the scalar-variance KL check returned a missing value, and the flat shared-β likelihood returned a second moment of zero instead of one. After the two narrow KL fixes, the focused standalone test has two expected equivalence failures and two passing KL checks. It is deliberately excluded from the general test runner until the acceptance comparison can run and pass. Both KL checks also run in the existing multimodal update test file.

## Closest expressible pooled comparison

For a diagnostic, both fitters received the same 32 by 12 nonnegative simulated matrix, outcomes, two factors, point-exponential L/F priors, Normal β, five sweeps and tolerance 10⁻⁵. The multimodal fit used one inner sweep, `prior_update = "ebnm"`, feature-specific noise and the corrected survival ELBO. Both used the uncentered projection because the multimodal fitter cannot center it. The single fit used `canonicalize = TRUE` and `w_S = 1`. Each fitter used its own initialization. This run was made after the KL fixes. It is not the registered centered comparison, and the ELBO difference includes differences in the objective and KL accounting.

| Quantity | Maximum absolute difference | Difference divided by max(1, magnitudes) |
|---|---:|---:|
| E[L] | 0.368667 | 0.0353423 |
| E[F] | 0.0760769 | 0.0760769 |
| E[β] | 0.582739 | 0.582739 |
| Final ELBO | 356.705 | 0.583668 |

The seed was 20261112. Reproduction uses the fixture in `tests/test_multimodal_m1_equivalence.R`, with `center_projection = FALSE` for the single fit and the unsupported multimodal argument and controls omitted. Different initializations and update sequences prevent attributing these gaps to one cause.

## Differences requiring a design decision

- The canonical single fitter updates β, L, then F for each factor. The multimodal loop updates L, F, then β and canonicalizes each factor immediately. The single fitter canonicalizes after the sweep. Their initialization and β warm starts also differ.
- The single fitter centers the survival projection globally or within strata when requested and uses strata in Cox risk sets. The multimodal fitter currently uses neither operation. Its optional weighted Cox-intercept profiling is a different operation.
- The single fitter applies `w_S` to survival updates and to its monitored objective. The multimodal loop has no corresponding control.
- The single fitter stops on relative full ELBO change by default and applies its dead-program rule. The multimodal fitter stops on relative reconstruction and predictor change and uses its own handling of degenerate updates. Both record an objective each outer iteration, but they do not monitor the same objective under all options.
- Prior fitting and KL sign conventions differ. The single fitter records negative KL contributions. The multimodal objective subtracts positive KL contributions. These signs must remain consistent if code is shared.

Matching the registered contract requires coordinated changes to the sweep, initialization, centering and strata, objective, and stopping behavior. This exceeds a small alignment, so no fitter rewrite was attempted. Two paths are available for review: extend the multimodal loop to mirror the canonical single sweep and prove the reduction, or build a shared M-modality sweep around the canonical single-modality algorithm. Either needs a separate design and regression plan before implementation. The latter may reduce long-term drift but still requires explicit validation of the M > 1 formulas.

## KL edge cases

`multimodal_yfb_coordinate_kl()` now repeats a scalar pseudo-observation variance to match the observation count before indexing. This includes one Gaussian normalizer per observation, matching the November 6 single-modality correction. `multimodal_yfb_update_beta_shared_k()` now uses the existing shared β update for an exactly flat likelihood, preserving the fixed prior's first two moments and zero KL. Nonfinite precision and an inconsistent zero-precision signal fail explicitly. A review also found that the shared-prior refit treated a flat likelihood as finite pseudo-noise. The refit now excludes flat coordinates and retains the current prior if all coordinates are flat. Tests cover scalar variance, Normal and point-Laplace shared priors, and the informative-only and all-flat prior refits. These fixes do not establish equivalence of the complete fitters or validate pruning generally.

Focused verification with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript tests/test_multimodal_yfb_updates.R` passed 13/13. The focused trace checks passed 7/7, including the prior-refit regression. The registered equivalence test remains failing by design pending implementation. The full regression suite passed 804 tests with zero failures after the final patch. Independent review passed 20 focused checks and approved the bounded changes. The suite emitted its existing bounded β optimizer profile-likelihood fallback warning.

## Text for the shared status files

The M = 1 multimodal equivalence gate remains open. The current multimodal fitter lacks projection centering and stratified risk sets and differs from the canonical single-modality sweep and stopping rule. Two KL edge cases were corrected and tested. A design decision is needed before aligning the fitters. These results do not authorize adoption of the multimodal fitter or a new model.
