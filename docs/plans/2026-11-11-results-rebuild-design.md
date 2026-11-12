# Results that remain aligned with the model

Date: 2026-11-11. Author: Andrew Walther.

Status: full-workflow design. At G delivery no pipeline was implemented and no
result or runner was moved. Andrew subsequently authorized the
[small local software pilot](../../results/pilot/README.md), now implemented and
[verified](../findings/2026-11-11-local-results-pilot-findings.md). Full scientific
workflow implementation and broader migration remain separate tasks.
The later approved [267-file early-artifact migration](../findings/2026-11-11-early-benchmark-migration-findings.md)
is complete, with current November results, snapshot and checkpoints retained.
Completing G meant delivering the design, not executing the cleanup.
Revised November 11 after Andrew's feedback: local execution is the default,
the collection must accept new studies, and external reproduction is an explicit
acceptance criterion. These design clarifications did not themselves install or run a pipeline.

## Problem and recommendation

The [0b inventory](../findings/2026-11-11-results-inventory-findings.md) identifies
204 scripts, 251 result/exhibit locations and 82 report sources. Historical
diagnostics, reported evidence and current model code share directories. A figure
can therefore remain plausible after its generating model changes. Several
writers also have hardcoded destinations that differ from their current location.

Use one explicit registry of result families and a dependency-tracked rebuild,
with a small `targets` pilot before adopting it for long fits. Rebuilt analyses
must call the same fit and prediction functions, use explicit model settings and
write immutable run folders. A current-release index should advance only when
all required outputs pass verification and review. Frozen reports retain their
original release and never silently acquire new results.

This is a proposed architecture. Andrew's November 11 preference for a smaller
active tree and archival of past work informs it. His approved branch cleanup
does not authorize a results-folder migration.

The intended next deliverable is an implemented workflow, not an indefinite
advice document. That is a separate task with a bounded pilot, then a staged
migration. G supplies its specification. Andrew's review of this design does
not by itself move historical files or start the week's analysis experiments.

The implementation should meet four practical goals:

- Refresh affected results after a model change and identify stale results before reporting.
- Provide a browsable collection of tables and figures for papers, the prelim and presentations.
- Keep every result tied to its model, data and analysis version, including comparisons with several model variants.
- Let another researcher reproduce synthetic results from a clean checkout and reproduce real-data results after obtaining and preparing the required data.

## Result register and proposed scope

The inventory remains a dated audit. Its [script table](../findings/2026-11-11-results-scripts-inventory.csv),
[folder table](../findings/2026-11-11-results-folders-inventory.csv) and
[report-source table](../findings/2026-11-11-results-report-consumers.csv) supply
paths, source expressions, candidate writers and consumers. The companion
[script dispositions](2026-11-11-results-rebuild-script-dispositions.csv) propose
a role for every inventoried script. They are a migration discussion, not proof
that a script reproduces a result or is safe to remove. The proposal has 23
current candidates, 166 historical scripts and 15 retired references. A historical
simulation can supply a design to port without making its old runner a default
target. A current helper is a dependency, not an independently run analysis.

At implementation, keep the maintained result register in
`results/registry.yml`. Each family gets an identifier, scope, producer function,
input targets, outputs, consumer documents, model applicability, required checks
and resource class. Do not repeat scientific verdicts from PROJECT_STATUS.qmd.
The status document links to the release index and findings. Dated findings and
attempt logs remain the evidence behind scientific decisions.

| Family | Existing source candidates, relative to `results/benchmark_sim/` unless stated | Inputs → outputs → consumers | Rebuild policy |
|---|---|---|---|
| Simulation recovery and calibration | `run_yfb_survival_weight_p2000.R`, `run_yfb_supervision_strength.R`, `run_yfb_near_threshold.R` and summarizers | Registered DGP, seeds, fitter/options → recovery, false-active, convergence and prediction aggregates → simulation findings, progress book, writeup and prelim figure builders | Core for relevant fit/update/prediction changes. Port registered designs, preserving historical results separately. Failed fits stay in denominators. |
| Training CV and final expression fit | `run_yfb_supervised_real.R`, `summarize_yfb_supervised_real.R`, `run_yfb_supervised_real_refit.R` | Approved preprocessing, training inputs, folds, model specification → CV tables, final fit and program counts → real-data findings and current report exhibits | Core for model, input, fold, screening or prediction changes. The historical refit's weight search does not authorize new tuning. Current fixed settings come from the approved model specification. |
| External validation | `score_yfb_supervised_external.R`, `save_yfb_external_patient_scores.R` | Fixed selected model, frozen preprocessing/prediction, donor mapping → aggregate within-series and pooled-stratified scores, intervals and fixed-cutpoint curves → external findings and KM builders | Core downstream scoring after the model is frozen. External outcomes never feed selection. Keep donor assays together for intervals. Patient scores remain ignored. |
| Fixed comparators and program rules | `run_prompt6_*`, `run_prompt7_*` | Registered training protocol, method implementations and comparator versions → aggregate comparisons, fixed risk groups and rule audits → November comparator findings and exhibits | Core comparison tables refresh when the SSBMF arm changes. Cache unchanged comparator fits only when their own inputs/settings/code/version hashes agree. No external retuning. |
| Program biology and sensitivity | `results/benchmark_sim/` pathway/cohort runners, `results/multimodal_real/` interpretation helpers | Current program loadings, feature universe and versioned annotation sets → enrichment/concordance/sensitivity aggregates → biology findings and current reports | Conditional core for applicable model/data changes. Register the intended analysis first, with teach-back where the ROADMAP requires it. Historical program numbers cannot stand in for new programs. |
| Multimodal and LB comparison | Historical multimodal and matched LB/YFB runners listed in 0b | An approved future model/registration → parity/recovery/CV/scoring aggregates → future findings and reports | Disabled until the respective task is authorized and reviewed. Do not include these in the current expression release by filename alone. |
| Exhibits and reports | `paper/prelim/make_*`, `docs/progress_book/figures/`, dated meeting builders | One release's aggregate tables and annotations → figures/tables → report source and rendered document | Core exhibits rebuild when their inputs or builder change. Render into a new report version. Historical circulated versions stay frozen. |
| Method-development diagnostics | October mechanism screens, pruning reruns, old rank/preprocessing studies and legacy runners | Original registration, code revision and saved results → dated diagnostic findings | Historical, excluded from default rebuild. Replaying requires an explicit historical environment or a new registration. No pruning adoption. |
| Publication and sync | Book deploy and prelim sync/hook scripts | Approved rendered release → published site or dissertation mirror | Separate explicit action, never a dependency of a fit or default rebuild. |

“Core” means required to claim a complete current release. It does not mean all
families rerun after every edit. A plotting-only edit rebuilds its exhibits.
A changed fitter invalidates all applicable fits and descendants. An annotation
update invalidates biology and its exhibits. A changed external donor map
invalidates the affected external scoring only. These dependencies must be
explicit and verified, rather than inferred from dates or script names.

## Adding and selecting results

The analysis family is the unit of registration. One R script can contain several
experiments, provide a helper, summarize old outputs or publish a report. The nine
top-level folders are not nine current analyses, and the 204-script mapping is not
a list to execute wholesale. The proposed classifications need a scientific
relevance check before selecting the first active collection. A historical
negative result can remain valuable evidence without being rerun by default.

For each selected family, record its question, registered design, model variants,
modalities, input cohort set, settings, producer, expected output schema, exhibit
builders and intended claims. Add tags such as simulation recovery, prediction,
program biology or numerical diagnostic. An aggregate catalog should expose
the model/run identifier, status, table and figure links, supporting findings,
and whether an artifact is current, stale, historical or awaiting review.
This is the entry point for choosing presentation material. A current result
must still be scientifically appropriate for the claim being made.

Adding a family should require a registry entry, a function using the common
fit/predict interface, output/manifest definitions, validation checks and a small
synthetic smoke case. Existing families should not require editing for that
addition. If the family becomes required for a release, the release specification
must explicitly include it and remain incomplete until its outputs are reviewed.
Otherwise it can be an optional extension. A changed family invalidates its own
outputs and consumers without forcing unrelated fits to repeat.

Task E's YFB/LB comparison, task C's methylation-only fit and task D's multimodal
studies can therefore join the collection when their own work is authorized.
Keep model variant and modality in the registry and each arm's manifest. Shared
training folds, input fingerprints and paired evaluation definitions should be
recorded for comparisons. A comparison can legitimately include several model
specifications under one code release. Its catalog must identify each arm rather
than pretending every table was generated by one identical fitted model.

## Reproduction from another checkout

The first implemented release must supply documented setup commands, pinned
dependencies, a synthetic fixture or deterministic generator, fixed seeds and
the build instructions. Test that path in a clean checkout, with no untracked
local cache or personal absolute path available. Compare deterministic aggregate
checksums or predeclared numerical tolerances. Keep a machine-readable list of
the generated artifacts and failed or unavailable families.

Real-data reproduction additionally requires data access. A clone alone cannot
provide the ignored patient data. Document acquisition sources and access needs,
the expected input schema, cohort/donor matching, preprocessing, feature-selection
rules and non-identifying data-version fingerprints. Verify the build from those
inputs rather than an undocumented saved fit. A user without the data should be
able to run the synthetic workflow and render shipped aggregate exhibits, while
real-data targets clearly report that their inputs are unavailable. They must
not be silently skipped in a build claiming a complete real-data release.

Use the same contract here and later in the package/paper repositories. External
reproduction is an implementation acceptance test, not a claim already established
by this design or by the current scripts.

## Model state, provenance and release completeness

Every run has a machine-readable manifest. Record the repository commit, dirty
state and code-content hashes, explicit fitter/predictor and all options, priors,
K, survival weight, centering, strata and cohort definitions, preprocessing and
screen settings, input fingerprints, fold assignment fingerprint, random seeds,
package/compiler versions and execution environment. Use stable non-identifying
input identifiers in tracked manifests. Patient identifiers, per-patient folds
and private paths belong only in ignored manifests beside private checkpoints.

A package version alone is insufficient if a local build changes without a
version bump. Record the installed package source commit/build fingerprint.
An unrelated documentation commit should not force model refitting: store the
full commit for provenance, but use relevant dependency content hashes to decide
staleness. A dirty development run is exploratory and cannot advance the released
current index without a reproducible committed source state.

Use `YYYY-MM-DD-<source-commit>-<spec-hash>-<attempt-id>` as a run identifier.
The attempt suffix prevents collisions when repeating an identical specification.
Write into a temporary directory and finalize only on successful validation.
The manifest records required versus optional families, statuses, warnings,
convergence failures, expected and observed replicate counts, output checksums
and source-run identifiers for reused artifacts. Every aggregate table must
include a run identifier or have an adjacent provenance sidecar named for that
file, resolving its manifest and output checksum. Figures carry the same linkage
in a sidecar, and captions state the model state. Consumers verify that linkage
before reading an artifact. A reused comparator artifact
keeps its own provenance. It is not relabeled as newly fitted.

`results/current.yml` points to one reviewed release manifest, not independently
chosen “latest” folders for each table. Promotion must fail if a required family
is missing, stale, failed, unreviewed or incompatible with the model specification.
Technical validation and independent scientific review are separate requirements.
For results entering documents, record the independent reviewer, review date,
findings/disposition link and outstanding issues. Passing hash or schema checks
does not establish scientific acceptance. Prior current remains available,
explicitly labeled with its older model state.
An in-progress run can be inspected but cannot be presented as the updated model.
Current reports check the manifest before rendering and stop on stale inputs.
Historical reports instead resolve their pinned release and display its date and
model state. There is no silent fallback from a missing new table to an old one.

## Proposed layout and preservation of current production

```text
results/
  README.md                 # navigation, current release and historical index
  registry.yml              # maintained family/dependency register
  current.yml               # reviewed release pointer
  pipeline/                 # thin functions calling canonical code/package
  runs/<run-id>/            # aggregate outputs and public provenance only
  history/README.md         # index of preserved historical locations and states
  benchmark_snapshot_2026-10-14/  # immutable, unchanged
data/results/<run-id>/      # ignored fits, scores, folds and private manifests
_targets/                  # ignored local workflow store, may contain patient data
```

Start with `history/README.md` indexing the existing historical paths in place.
This reduces active navigation without immediately moving files. Later physical
relocation is a separate approved migration, with a source-to-destination map,
consumer checks and checksums. The repository's `docs/` rule still applies:
findings, plans and reviews remain in their type folders, with superseded status
in their READMEs. Do not create an archive folder for those documents.

Current production keeps its existing paths until parity is demonstrated for a
small new release, all direct and dynamic readers are checked, and reports can
select a manifest. Temporary compatibility exports may be necessary, but must
record their source release and have one writer. Never replace old results under
the same filename and assume downstream figures refreshed.

| Existing area | Proposed treatment before any physical move |
|---|---|
| `benchmark_sim/` | Separate core candidates from historical diagnostics using the companion mapping. Wrap or port only selected functions. Preserve dated outputs with their original registrations. |
| `multimodal_real/`, `cohort_lmm_sim/`, `multi_cohort_sim/` | Index historical studies and reusable helpers. Future multimodal core targets require their own accepted model and registration. |
| `feature_selection/` | Preserve old experiment. A new training-only screening target follows task F's approved design, not an automatic replay. |
| `tables/`, `figures/`, `legacy/` | Preserve as historical evidence. Audit legacy writers that still target top-level paths before consolidation. Do not run them in the default build. |
| Frozen benchmark | Keep manifest, outputs and tag immutable. It is a reference, never a rerun destination. |
| Report and presentation assets outside `results/` | Leave circulated assets in place. New builders consume a pinned release and write to a new report version. Fix the August meeting builder's June output destination before ever enabling it. |

## Checkpoint retention and archive verification

Keep the current and last reported model's required checkpoints locally, plus
any checkpoint required by an unresolved review or active diagnostic. Keep small
aggregate tables, manifests and dated findings in Git. Intermediate checkpoints
may be proposed for archival only after their consumers and recovery needs are
known. Do not remove a failed fit that supplies a failure denominator or explains
an unresolved numerical problem.

Older bulky checkpoints may move to durable institution-approved storage, with
the same access restrictions as the source data. Record checksums, file counts,
sizes, model/input/environment identifiers and a non-sensitive archive locator.
Verify the copy, restore a representative complete analysis dependency chain in
an isolated location, and regenerate its aggregate summary before considering
local removal. If restoration fails, retain the source and report the failure.
The storage destination, access and funding/retention period still need Andrew's
decision. An ignored folder or Git tag alone does not back up ignored checkpoints.

## Local execution and workflow choice

Recommend `targets` for the pilot because it tracks declared R dependencies and
file targets. Its [walkthrough](https://books.ropensci.org/targets/walkthrough.html)
and [target reference](https://docs.ropensci.org/targets/reference/tar_target.html)
describe updating downstream steps and tracking returned files. This only works
for declared inputs. A `system()` call running a sourced historical script does
not automatically expose every file, environment variable or compiled module it
uses. Wrap analysis functions with explicit arguments and declare code, config,
input and compiled-source fingerprints. Returned file targets must enumerate
the actual files. Missing or changed files invalidate the corresponding family.

The [package documentation](https://books.ropensci.org/targets/packages.html)
also distinguishes environment management from dependency tracking. Pin the
environment and explicitly depend on the fitter package fingerprint so a changed
installed binary invalidates results. The [local-data documentation](https://books.ropensci.org/targets/data.html)
explains that the workflow store serializes objects. Treat the entire `_targets/`
store as private and ignored, with no patient-bearing store metadata in Git.

| Choice | Fit for this project | Proposed role |
|---|---|---|
| `targets` | R objects, folds, fits, tables and reports can form an explicit graph. Long jobs and hidden dependencies require careful wrappers. | Preferred pilot, then scale only after invalidation tests pass. |
| Makefile | File dependencies can be explicit, but R function/config/package changes require additional dependency declarations. | Useful outer command wrapper, not a second freshness authority. |
| Numbered runners plus manifest checks | Low initial disruption, but dependencies must be maintained manually and can drift. | Fallback if the pilot fails. Require content-hash checks, not timestamp-only freshness. |
| Editing hook | Can remind a developer to rebuild. | Optional reminder only. It cannot establish freshness or scientific acceptance. |

Run the pilot and normal workflow locally. Longleaf access, SLURM knowledge and
a cluster backend are not prerequisites. Begin serially with a small synthetic
fixture and no patient data. Measure representative wall time and peak resident
memory before choosing local concurrency. More workers make extra copies of
matrices and fit objects, so start with one worker and avoid nested parallelism.
Give each task one output writer and deterministic seeds. Cap BLAS threads to
the intended local CPU budget if enabling several worker processes.

The historical [attempt log](2026-10-16-round2-attempts-log.md) provides context,
not forecasts for the changed model. Its Task5-F1 entry reports 24 expression
fits across four gene counts taking 553.487 seconds in total, with a maximum R
allocation estimate of 819.3 MB. That memory measure is not peak resident memory.
Its Task5-M1 entry records a historical-framework 10,000-feature methylation fit
taking 963.827 seconds
and reaching the 500-iteration cap without convergence. These examples support
trying local execution first. They do not establish the cost or convergence of
future supervised grids or multimodal fits.

Classify targets as quick exhibits, ordinary fits or long sweeps. Default status
and report commands should inspect existing artifacts without launching fits.
A requested refresh runs only stale prerequisites of the chosen families. Long
local sweeps should have an explicit estimated fit count, time/memory budget,
resumable per-seed or per-fold outputs and progress reporting. Preserve completed
valid work across interruptions, but never count an incomplete sweep as a
complete release. Do not launch a large fitting job automatically on every save
of model code. Detect staleness first, then use an explicit rebuild action.

If measured workloads exceed Andrew's available local time or memory, discuss
reducing an authorized workload or moving it before executing elsewhere.
Longleaf is an optional later backend requiring Andrew's approval and practical
setup/run guidance. The [distributed-computing documentation](https://books.ropensci.org/targets/crew.html)
provides a possible future route, not a dependency of the local pilot. The
analysis functions and manifest contract should remain the same across backends.
Finalized outputs are recorded only after success. Retry exhaustion and
convergence failures remain visible and prevent promotion where required.

## Package and paper pathway

The [November 2 pathway](2026-11-02-package-and-paper-pathway-plan.md) assigns
development evidence here, adopted methods to `SSBMF`, and reproducible paper
analyses to `SSBMF-paper`. Its instruction to leave this repository's historical
work in place differs from the later request to reduce the active results tree.
The staged proposal reconciles those aims: index history in place first, approve
any physical migration separately, and preserve all scientific evidence.

Until the package meets the pathway's acceptance/parity gates, thin targets here
call the canonical development fit and prediction functions. Once a package
version is accepted, port the registered analysis designs to `SSBMF-paper`, pin
the package fingerprint and verify parity. Do not maintain two independent
scientific implementations or adopt the standalone rebuild branch by default.
Historical diagnostics and findings remain here. The November 9 training-only
selection rule governs new work, superseding any older suggestion to use external
validation for model choice. Nothing in this design starts tasks A to F.

## Implementation acceptance and decisions for Andrew

The bounded pilot was authorized and completed in a separate task. Its
[verification findings](../findings/2026-11-11-local-results-pilot-findings.md)
record the implemented subset of the acceptance contract below. Real-data
families, scientific promotion and physical migration remain later work:

1. Register one synthetic simulation, a fit, an aggregate table and an exhibit.
   Run locally and reproduce from a clean checkout with the documented setup.
   Verify that an unchanged rebuild performs no fits and reproduces checksums
   for deterministic aggregate artifacts. For serialized or rendered artifacts
   with variable metadata, define numerical/content tolerances in advance and
   verify those rather than demanding byte identity.
2. Change a relevant fitter dependency and verify that affected fits and every
   downstream exhibit become stale. Change a plotting function and verify that
   no fit reruns. Change an input/config/package fingerprint and verify the
   expected invalidation. An unrelated prose edit must not invalidate a fit.
3. Delete an output and inject a failed branch. Verify that promotion and current
   report rendering fail, with the old release still available and labeled.
   Confirm the complete registered seed denominator is retained.
4. Add a small real-data family only after private-store handling, paired folds
   and donor-safe scoring are reviewed. Compare aggregates with the fixed
   historical reference using declared numerical tolerances, then audit readers.
5. Approve the manifest/index contract and any compatibility exports. Only then
   plan migration and a verified checkpoint archive. Deployment remains explicit.
6. Add a second small result family without changing the first. Verify its catalog
   entries, model-variant labels and freshness checks. This tests extensibility
   before adding this week's larger comparison or multimodal studies.

Andrew approved the bounded `targets` pilot, which is now complete. Remaining
decisions are the scope and manifest/layout contract for the first scientific
release, which conditional families it must include, and where durable private
archives will live. The default
recommendation is the expression model's simulation, training CV, final fit,
external scores, fixed comparators and associated exhibits. Biology is required
when making current biology claims. Multimodal/LB expansion and new sensitivity
studies remain gated by their own task approvals. Full scientific implementation
and migration are not started by completion of the software pilot.


## Design verification

Independent review checked task G's seven requirements, the package transition
and the 204-path mapping against the inventory. All paths appear exactly once.
The review's requests for per-output provenance links, explicit independent
scientific review and appropriate reproducibility tolerances were incorporated.
Local document links and mapped source paths were checked. No fit, pipeline,
archive migration or model change was executed for this design.


A second independent review checked the local-first and extensibility revision
against Andrew's feedback and the historical runtime entries. Its clarification
that methylation timing came from the historical framework was incorporated.
The full design is not implemented by the bounded pilot. Local execution remains the default, with
any later Longleaf use requiring a separate decision and practical guidance.
