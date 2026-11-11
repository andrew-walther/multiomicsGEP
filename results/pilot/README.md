# Local results-rebuild pilot

This is a small implemented software pilot of the
[G design](../../docs/plans/2026-11-11-results-rebuild-design.md).
It uses synthetic data only. Its tables and figures demonstrate reproducibility
and dependency tracking, not evidence about SSBMF's scientific performance.
The existing result folders and production reports are unchanged.

## Run from the repository root

Tested with R 4.5.2. Install the exact locked package versions first:

```sh
Rscript results/pilot/install.R
Rscript results/pilot/run.R status
Rscript results/pilot/run.R build
Rscript results/pilot/run.R catalog
Rscript results/pilot/verify.R
```

Setup uses `package-versions.csv`, including `targets` 1.12.0 and its hard
dependencies. It reuses exact versions visible in existing R libraries and
installs missing or mismatched versions in ignored `results/pilot/library/`.
It changes neither the global library nor package versions there. Source-package
installation needs an R-compatible compiler/toolchain. Missing archived versions
or installation failures stop setup. No fallback to newer versions is allowed.

`status`, the default command, checks freshness without fitting. `build` runs
stale targets explicitly, serially with a local R controller and one isolated
fit subprocess at a time. No cluster account,
background worker, publishing step or patient data is used. `catalog` refuses
stale, incomplete or checksum-mismatched output. It does not silently show the
previous run as current after a failed build.

Generated runs, cached fit objects and the local dependency library are ignored.
Each run lives in `results/pilot/local/runs/<run-id>/`. It contains each family's
CSV table, PDF figure and adjacent provenance YAML, plus `manifest.yml` and
`catalog.csv`. The local `latest.yml` points to a completed **software pilot**.
It never changes `results/current.yml` or declares a scientific release.

## What the two families do

`prediction` generates the README's synthetic example from configured settings,
standardizes inputs using training rows only, calls the canonical
`fit_ssbmf_yfb()` and `predict_ssbmf_yfb()`, and reports held-out concordance,
training reconstruction RMSE and iteration count. Its figure displays the
synthetic held-out concordance. Scores retain the fitted risk orientation.
A capped fit, nonfinite state or nonfinite score fails the build.

`input_audit` reports the actual sample/feature counts and training event fraction
from the same fixture. Its figure displays the event fraction. It depends on
the synthetic inputs, not on fitting. This is the second independently registered
family used to demonstrate that the collection can grow.

The `results_pilot` section of `config/globals.yml` specifies the fixture and
explicit model options. Rank, iteration cap and tolerance inherit the existing
`cavi` settings. The pilot does not change scientific benchmark settings.

To add another small family, define a function returning a numeric `metric` and
`value` table in `functions.R`, then add its ID, title, single R `command` and
`figure_metric` to `registry.yml`. The command is trusted repository code.
The graph discovers its referenced targets and adds its table/figure/provenance
targets without changing either existing family's implementation. The current
figure template supports fraction/concordance metrics on a zero-to-one axis.
Other scientific plots and real-data families require a later implementation.

## Freshness, versions and failure handling

The graph tracks the fixture settings, model settings, verified package/runtime
state and R/C++ model files. The pilot conservatively tracks every R/C++ file in
`code/`, so an unrelated model-module edit can also invalidate the fit. The
future migration can narrow that boundary after auditing sourced dependencies.
Function bodies are tracked independently. A figure-only change leaves the fit
cached. An unrelated prose edit does not invalidate fitting.

The fit subprocess prevents the canonical fitter's compiled-helper and RNG
state from changing the controller's dependency graph. Prediction helpers load
in a private R environment. Both paths call the existing canonical code.

The manifest records the source commit and dirty flag, source/config content
hashes, explicit model/fixture settings, R/platform and exact package versions.
Content-derived run IDs preserve previous runs. Existing different artifacts
cannot be overwritten. Each CSV/PDF pair has a sidecar resolving its run manifest
and checksums. The catalog includes the family, title, run ID and artifact paths.
The PDF device suppresses timestamps to make rebuilding deleted figures stable.

Dirty or exported-source builds are labeled exploratory. Every run is explicitly
marked `scientific_release: false`, including clean-commit runs. This pilot does
not implement scientific review/promotion, patient data handling, donor-aware
external scoring, checkpoint archival or migration of historical runners.
Those remain requirements of the larger G design.

## Verification and function guide

`verify.R` copies source files into an owned temporary directory and uses
independent output stores. It tests unchanged-build caching, deterministic
aggregate reproduction, figure/model/config/prose invalidation, adding the
second family, missing-artifact repair and failure preservation. It deletes only
its temporary copies. No project source is modified by these checks.

Fast unit checks are included in `Rscript tests/run_tests.R` and can also be run
with `Rscript tests/test_results_pilot.R`. They require the existing core
dependencies, not `targets` or a model fit. The separate workflow acceptance
script requires setup because it runs the actual target graph.

| Function or script | Purpose, inputs and outputs |
|---|---|
| `install.R` | Read the lock, process hard dependencies in order and install exact missing versions locally. |
| `pilot_settings()` | Read YAML and return separate fixture/model lists, including inherited CAVI settings. |
| `pilot_runtime()` | Verify visible raw package versions and return package/R/platform provenance. Raw version strings preserve CRAN's hyphenated versions. |
| `pilot_source_state()` / `pilot_context()` | Read Git state and declared file hashes, then derive the run identity and manifest fields. No fit is run. |
| `pilot_simulate()` | Turn fixture settings into standardized train/test matrices and survival vectors. Only training means/SDs standardize test rows. |
| `pilot_fit()` | Pass explicit settings and training rows to the canonical fitter, rejecting capped/nonfinite fits. Returns the existing fit structure. |
| `pilot_prediction_summary()` / `pilot_input_summary()` | Convert the fit/input objects into aggregate metric/value data frames for the two families. |
| `pilot_finalize()` / `pilot_artifacts()` | Finalize immutable files, make a base-R PDF and write adjacent provenance. Fail on conflicting content. |
| `pilot_release()` / `pilot_catalog()` | Finalize a complete pilot manifest/catalog and local pointer, or refuse to read stale/missing/changed artifacts. |
| `_targets.R` / `run.R` | Declare dependencies and provide explicit local status/build/catalog commands. Registry entries generate each family's targets. |
| `verify.R` | Exercise rebuild/failure behavior in isolated temporary source copies. Returns a passing-check count or exits nonzero. |

The implementation uses base R for simulation, files and PDF plotting, `yaml`
for readable configuration/manifests, `digest` for a run identity, `targets`
for dependency tracking and its existing `callr` dependency to isolate the
canonical fit process. The canonical fitter supplies its existing statistical
packages. No new fitting algorithm or scientific comparison is introduced.
