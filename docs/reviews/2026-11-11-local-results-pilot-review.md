# Independent review: local results pilot

Date: 2026-11-11
Author: Andrew Walther

## Verdict

Approved as a bounded local software pilot. The implementation uses synthetic inputs and the canonical YFB fitter and predictor. It does not establish scientific performance, authorize model adoption, or migrate the production results pipeline.

## Scope and checks

The review compared the implementation with `docs/plans/2026-11-11-local-results-pilot-plan.md`. It inspected `results/pilot/functions.R`, `_targets.R`, `run.R`, `install.R`, `registry.yml`, `package-versions.csv`, and `verify.R`, as well as `tests/test_results_pilot.R` and the `results_pilot` settings in `config/globals.yml`.

The synthetic fixture uses training-set centers and scales for held-out rows. The second family summarizes the same inputs without fitting. The artifact workflow records source, configuration and package fingerprints, checks artifact checksums before presenting the catalog, and preserves the previous completed pointer when a build fails. The pilot is serial and does not require private data or cached fits.

`Rscript results/pilot/verify.R` passed all 14 acceptance checks in an isolated temporary source copy. These cover a complete initial build, an unchanged rebuild, reproducibility across independent stores, figure-only changes, model and fixture invalidation, a documentation-only change, declaration of the second family, missing-artifact detection and repair, pointer preservation and catalog refusal after a failed build, and the catalog command in a fresh R process without refitting. The full repository test suite is reported separately.

## Limits

The installer reuses an exact-version package already visible in an R library when available. It therefore does not necessarily copy every locked package into the project-local library. The pilot-specific settings inherit K, maximum iterations and tolerance from the existing `cavi` configuration. The fixture and its figures are software checks, not evidence of recovery or predictive performance.

During review, a graph variable collision caused the catalog command to report status instead. That collision was corrected. The acceptance harness also exposed differences between graph evaluation environments and fitter side effects. The fitter now runs serially in an isolated local R process, and the CLI and acceptance calls use the same graph environment. The independent 14-check rerun passed after these corrections.
