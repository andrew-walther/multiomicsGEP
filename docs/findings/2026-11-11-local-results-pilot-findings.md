# Local results-rebuild pilot verification

Date: 2026-11-11. Author: Andrew Walther.

The authorized [small pilot](../plans/2026-11-11-local-results-pilot-plan.md)
is implemented in `cc11b34`, following the registration in `9613ab6`.
It is a software demonstration, not a new scientific result or an implementation
of the full production migration in G. All original six requested items were
already complete before Andrew authorized this extension.

## Delivered workflow

The [pilot guide](../../results/pilot/README.md) documents setup and the
`status`, `build`, `catalog` and verification commands. `targets` executes jobs
locally and serially. The canonical fitter runs in a fresh local R subprocess
to isolate its compiled helpers and RNG side effects from dependency tracking.
The canonical prediction and concordance helpers load in a private environment.
No fitting equations, priors or scientific defaults were changed.

The `prediction` family uses the existing README synthetic generator with
training-only standardization. It calls the canonical YFB fitter/predictor and
writes an aggregate table and PDF figure. A capped/nonfinite fit fails the build.
The independently registered `input_audit` family summarizes the same synthetic
inputs without fitting. Registry entries create family targets without changing
the first family's implementation.

Each family has a table, figure and adjacent provenance. An immutable run
manifest records source/configuration fingerprints, commit and dirty state,
runtime/package versions and artifact checksums. A catalog lists the two
families and their artifacts. A local pilot-only pointer advances after the
required outputs complete. Missing, stale or changed artifacts cannot be
presented as a fresh complete pilot. Previous completed runs are preserved.
Every output is explicitly marked as software-pilot evidence rather than a
scientific release. The nine pre-existing results folders remain in place.

## Verification evidence

| Check | Result |
|---|---|
| Fast unit checks, `Rscript tests/test_results_pilot.R` | 7 passed, 0 failed. Covers deterministic inputs, actual dimensions, raw version strings, immutable writes, manifest/catalog linkage, stale/tampered artifacts and incomplete-release preservation. |
| Workflow acceptance, `Rscript results/pilot/verify.R` | 14 passed. Independently rerun with the same result. About 21 seconds for the isolated acceptance run on this machine. |
| Final full suite, `Rscript tests/run_tests.R` | 801 passed, 0 failed, exit 0. The bounded-β profile-likelihood fallback warning remains visible. |
| Clean Git checkout at `cc11b34` | Documented `install.R`, `run.R build` and `run.R catalog` all exited 0. No ignored pilot library, fit cache, generated output or private input was copied into it. |
| Reproduction across the working and clean checkout | Both families' aggregate CSVs were byte-identical under the same tested R/package environment. |
| Unchanged clean-checkout rebuild | All 14 targets skipped, no dispatched fit. |
| Independent review | [Approved bounded software pilot](../reviews/2026-11-11-local-results-pilot-review.md), including the corrected fresh-process CLI check. |

The 14 acceptance checks cover initial completion, no-fit cached rebuilds,
independent output stores, figure-only changes without refitting, model and
fixture invalidation, prose edits without refitting, adding the second family,
missing-artifact refusal/repair, failed-family preservation and the fresh-process
catalog CLI. Test mutations occur only in owned temporary source copies.

A command-variable collision and differences in graph evaluation environments
were found during verification and corrected. Isolating the fitter removed its
compiled-state side effects from the controller. The PNG/Cairo device required
unavailable X11 libraries, so the implemented figures use base R PDFs with
timestamps disabled. No system graphics library was installed.

## Reproducibility scope and remaining work

The tested environment is R 4.5.2 with 57 locked hard-dependency packages,
including `targets` 1.12.0. The installer reuses exact versions already visible
and installs missing versions only in the local ignored library. The clean
checkout downloaded its missing pilot packages using that installer. Expected
warnings identifying initially missing packages preceded successful installation.
This verifies a clean checkout on this machine, not a clean operating-system
installation or cross-platform numerical identity. Source installation requires
the usual compiler/toolchain.

The pilot conservatively tracks all R/C++ files in `code/`. A model-module change
can therefore trigger more refitting than a future audited production graph
would need. Full-commit provenance can refresh cheap artifacts after a
documentation commit, while the model fit remains cached. Its fixed figure
template is for fraction/concordance metrics. Existing CAVI rank, cap and
tolerance are inherited rather than independently chosen for a research study.

The next owner can use this tested starting point to plan relevant real-data
families, richer figures, scientific review/promotion, resumable large sweeps,
historical-runner migration and checkpoint retention. Each needs its own agreed
scope. No task A to F, production publication, patient-data analysis, Longleaf job,
historical cleanup or external checkpoint archive was started here.


A fresh closeout review checked these counts and source commits against the
captured verification logs. Its stale-status findings were corrected in the
week checklist, plans index and full-workflow design, so the completed pilot is
now clearly distinguished from the future scientific implementation/migration.
