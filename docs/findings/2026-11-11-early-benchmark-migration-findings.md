# Early benchmark artifact migration

Date: 2026-11-11. Author: Andrew Walther.

Status: approved implementation in progress. Andrew approved the 267-file batch
in the [navigation proposal](2026-11-11-results-navigation-findings.md), provided
current results are preserved and old results are archived rather than lost.
The source checkout before this batch is `2a696b7`.

## Scope and preservation baseline

The [migration manifest](2026-11-11-early-benchmark-migration-manifest.csv)
records every tracked source/destination, byte count and SHA-256 before moving:
110 figure files and 157 table files. The sources contain no ignored checkpoint
or other local data object. Each contains one ignored Finder `.DS_Store`, which
will be moved unchanged alongside its directory without committing its contents.

Destinations are `results/history/early-benchmarks/figures/` and
`results/history/early-benchmarks/tables/`. No current November output,
frozen-benchmark file, local pilot, patient input or ignored checkpoint is in
this batch. The pre-move check recorded hashes for 5,009 tracked protected files
under benchmark outputs, the frozen snapshot and the pilot, plus size, mtime
and inode metadata for 3,180 ignored RDS files. Private RDS contents were not read.

## Reference changes and historical meaning

Only output path prefixes change in the four affected legacy R runners and two
historical SLURM launchers. The setup README's path examples follow those roots.
No model, fitting setting or scientific computation changes, and no old runner
or cluster job will be executed for this migration.

Three legacy report sources and the April 9 presentation sources contain
relative/constructed readers. Some were already broken by earlier directory
changes. This batch repairs the source artifact paths to the new archive roots,
without altering scientific text, rendering presentations or regenerating
circulated documents. The dated inventory keeps its original snapshot paths and
links to this migration record. No generated scientific output is rewritten.

## Verification and closeout

Pending completion: exact post-move checksums, protected-file checks,
reference/path validation, affected R/shell parsing, full repository suite,
independent review and updated continuation prompt/checklist.
