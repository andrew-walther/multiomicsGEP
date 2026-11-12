# Historical results index

This index organizes historical work at its existing paths. It does not move,
delete or authorize rerunning a result. For reported November evidence and the
local software pilot, start at the [results index](../README.md).

Historical means preserve the study with its original model, registration and
code revision. It does not mean the result is scientifically worthless. A failed
method can be important evidence without belonging in the default rebuild.

| Existing location | Historical role and evidence |
|---|---|
| [October 14 benchmark snapshot](../benchmark_snapshot_2026-10-14/README.md) | Immutable benchmark reference with its own manifest and tag. Never a new-run destination. |
| [Early aggregate tables](../tables/) and [figures](../figures/) | Early benchmark and prior/rank outputs. Some retired runners still write to these paths. These are the proposed first physical migration batch, pending approval and reader/writer checks. |
| [Legacy runners and reports](../legacy/) | Original full, modular-block and factor-modular studies. Some scripts retain pre-archive destinations. Preserve provenance before any path change. |
| [Cohort-indicator simulation](../cohort_lmm_sim/) | Earlier synthetic cohort-indicator runner. Its source and output registration are listed in the inventory. |
| [Shared/cohort-specific simulations](../multi_cohort_sim/) | Earlier simulation studies and reusable helpers. The main test suite sources `fit_pca_cox.R` here, so the whole directory cannot be moved without updating that dependency. |
| [Feature-selection experiment](../feature_selection/) | Earlier training-only screening experiment. A future task F study is a new registered result, not an automatic overwrite of this one. |
| [Multimodal studies](../multimodal_real/) | Historical matched expression/methylation fits, feature-count sensitivity and collapse diagnostics, plus reusable tools. These are not the future task D fitter by default. |
| [Benchmark studies and diagnostics](../benchmark_sim/) | Mixed model states. Current reported families are linked separately in the results index. October mechanism screens, old rank/preprocessing studies and pruning reruns remain historical. |

The [0b findings and detailed tables](../../docs/findings/2026-11-11-results-inventory-findings.md)
record producers, saved locations, report consumers and unresolved mappings at
their pinned source snapshot. The [proposed script dispositions](../../docs/plans/2026-11-11-results-rebuild-script-dispositions.csv)
separate current candidates, historical studies and retired references. Neither
is a blanket deletion list or an executable dependency graph.

Before a physical move, approve a source/destination batch, check direct and
dynamic readers and writers, preserve file checksums and scientific records,
and verify the affected build/test commands. Private checkpoints require their
own retention and verified archive decision. Moving directories does not remove
large blobs from existing Git history.
