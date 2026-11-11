# Small local results-rebuild pilot

Date: 2026-11-11. Author: Andrew Walther.

Andrew authorized this implementation after completing items 00, 0a, 0b, 0c,
G and 0. It is a bounded software pilot of the
[results-rebuild design](2026-11-11-results-rebuild-design.md), not a new study
or a migration of the existing results tree.

## Deliverables and acceptance

1. Build a deterministic synthetic fixture, call the canonical YFB fitter and
   predictor, and produce an aggregate prediction/reconstruction table and figure.
   Use signed point-Laplace molecular priors, Normal coefficients and the existing
   canonicalization option. Fixed smoke settings are separate from scientific
   benchmark settings and live under `results_pilot` in `config/globals.yml`.
   A capped or nonfinite fit must fail the build, not produce a successful release.
2. Register a second family summarizing the same synthetic inputs without fitting.
   Its registry entry declares the producer, figure metric and dependencies.
   Add this family without changing the first family's implementation.
3. Use local serial `targets`, an ignored project-local dependency library and
   ignored local caches/runs. Every artifact resolves to a manifest with source,
   input/specification and package fingerprints. Produce a navigable artifact
   catalog. Label all outputs software-pilot evidence, not current scientific
   results. Do not alter the production current-results pointer.
4. Verify an unchanged build reruns no fits, a figure-only change does not rerun
   the fit, a model dependency/config change invalidates the fit, missing output
   is detected, and failed/stale builds cannot be presented through the catalog
   as a complete fresh pilot. Preserve previous completed runs.
5. Reproduce deterministic aggregate output from a clean checkout without private
   data, cached fits or personal paths. Record exact tested package versions and
   installation instructions. Run the full repository test suite and a bounded
   pilot acceptance script. A fresh independent review precedes completion.

## Limits

No external data, new model comparison, scientific acceptance decision, Longleaf
job, production publication, old-output deletion or folder relocation is allowed.
The pilot runs locally with one job at a time. The canonical fit uses a fresh
local R subprocess to isolate its compiled-helper and RNG side effects from the
dependency graph. This execution adjustment was made during verification,
without changing the fitting algorithm or settings. Larger releases, checkpoint archival and
migration remain later tasks. Its demonstration figures must not enter the
prelim or paper as method-performance evidence.

## Progress

- [ ] Implementation, acceptance checks, independent review and handoff update.
  - Implementation and independent 14-check acceptance run complete. Clean-checkout verification and final full-suite rerun are in progress.
