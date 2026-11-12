# Results navigation and the first migration proposal

Date: 2026-11-11. Author: Andrew Walther.

Andrew requested an updated continuation prompt and action on results
organization after the original housekeeping, combined CV review and local
pilot were complete. Tasks A to F remain for his discussion with the next owner.

## Completed navigation stage

The [results index](../../results/README.md) now separates the tested local
pilot, current reported November evidence, the immutable benchmark and
historical work. Seven reported result families link to their saved aggregates
and findings, so readers need not infer their relevance from a runner name.
The [history index](../../results/history/README.md) organizes historical paths
in place and records important dependencies. It does not move or delete outputs.

The [continuation prompt](../prompts/2026-11-11-week-plan-continuation-prompt.md)
identifies all completed work and current verification evidence. It explicitly
requires discussion of the next scoped A-to-F or workflow item before executing
it. The existing handoff remains the detailed record of commits and limits.

This stage changes navigation, not analysis code, stored result contents or
scientific status. The dated 0b inventory retains its original snapshot rather
than being rewritten to imply a new analysis has run.

## First physical batch, subsequently approved

| Source | Proposed destination | Tracked files at this check |
|---|---|---:|
| `results/figures/` | `results/history/early-benchmarks/figures/` | 110 |
| `results/tables/` | `results/history/early-benchmarks/tables/` | 157 |

This would consolidate 267 early benchmark artifacts beneath indexed history.
It would not move `benchmark_sim/outputs/`, the current reported November
evidence, the frozen snapshot, the pilot or any ignored checkpoint.
Andrew subsequently approved this batch. The [migration record](2026-11-11-early-benchmark-migration-findings.md)
and [per-file manifest](2026-11-11-early-benchmark-migration-manifest.csv) record
its implementation and preservation checks. The table above retains the original
source/destination proposal.

At `94c7a74`, before these navigation documents were added, a literal source
scan for `results/figures/` and `results/tables/` found thirteen referencing source
documents/scripts:

- `PROJECT_STATUS.qmd`
- `longleaf_setup/README.md`, `longleaf_setup/run_modular_simulation.sl` and `longleaf_setup/run_v2_simulation.sl`
- The V1, V2 and V3 April 9 presentation sources under `presentation/walther_lab_meeting_04_09_2026/`
- `docs/findings/2026-11-11-results-inventory-findings.md`
- `docs/plans/2026-03-25-fit-modular-plan.md`
- `results/legacy/modular_sim_factor/run_factor_modular_simulation.R`
- `results/legacy/modular_sim_factor/run_prior_k_comparison.R`
- `results/legacy/modular_sim_block/run_modular_simulation.R`
- `results/legacy/full_sim/run_simulation.R`

These are literal references, not a complete dependency graph. Relative paths,
constructed paths and embedded report assets still need checking. The dated
inventory should retain its original paths and point to the migration record.
Historical circulated presentations should not be regenerated merely to change
a storage path. Decide each reference's treatment rather than replacing strings
across historical records indiscriminately.

If approved, the bounded implementation plan is:

1. Record source/destination and checksums for every moved file, including a
   check for ignored files or other local-only artifacts in the two sources.
2. Audit readers and writers, including dynamic and relative paths. Resolve
   legacy writer destinations without adopting or rerunning old analyses.
3. Move the approved batch, update current navigation and necessary executable
   paths, and preserve the historical study/code-revision interpretation.
4. Verify identical file contents, resolved current links and affected tests.
   Commit the move, reference changes and migration record together, then push.

The main test runner also sources
`results/multi_cohort_sim/fit_pca_cox.R`. This is a concrete reason not to move
every older-looking result folder in one batch. Larger consolidation of
`benchmark_sim/` needs the family-by-family dependency work in G. Checkpoint
archival additionally needs a storage destination and restoration verification.

## Verification

Tracked source counts were checked against Git. New local links and navigation
targets were checked against the filesystem. No model, runner, test, parameter,
patient input or saved result was changed, and no fit or publication was run.
The latest implementation suite remains 801 passing tests. No R suite rerun
was needed for this documentation-only stage.


Independent review checked navigation, handoff scope, indices and the proposed
batch counts. Its correction to include the root status document and two `.sl`
launchers in the pinned literal-reference scan was incorporated. The result is
13 referencing source files, with dynamic/relative dependencies still requiring
an audit before any move.
