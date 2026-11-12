# Early benchmark artifacts

This approved archive contains 110 figure files and 157 table files moved from
the former top-level `results/figures/` and `results/tables/` directories.
The [migration manifest](../../../docs/findings/2026-11-11-early-benchmark-migration-manifest.csv)
records each original path, destination, byte count and SHA-256. Contents were
preserved, not regenerated. See the
[verification record](../../../docs/findings/2026-11-11-early-benchmark-migration-findings.md).

These are historical cohort benchmarks and prior/rank comparisons, not a new
current scientific release. The [reported November evidence](../../README.md#november-reported-result-families)
and frozen benchmark remain at their existing locations. Do not overwrite the
archive by rerunning an old script against changed model code.

| Location | Contents |
|---|---|
| [figures/](figures/) | Historical per-cohort figures and later rank-selection exhibits. |
| [tables/](tables/) | Historical per-cohort and cross-dataset aggregate outputs. |

The older baseline synthetic/full/blockwise outputs were already archived under
`results/legacy/figures/` and `results/legacy/tables/`; they were not part of this
move. In particular, baseline synthetic figures remain separate from the later
rank-selection artifacts in this archive's `figures/synthetic/` directory.
Their source report and circulated renderings remain under `results/legacy/`.

The archive is not a default build target. For the local versioned software
workflow use the [pilot guide](../../pilot/README.md). A scientific rerun needs
its own registered design, compatible model version and new output location.
