# Early benchmark artifact migration

Date: 2026-11-11. Author: Andrew Walther.

Status: complete and independently reviewed. Andrew approved the 267-file batch
in the [navigation proposal](2026-11-11-results-navigation-findings.md), provided
current results are preserved and old results are archived rather than lost.
The source checkout before this batch is `2a696b7`.
The approved migration and reference changes are committed in `872b92b0`;
`e32a8a1` records the scope and pre-move preservation manifest.

## Scope and preservation baseline

The [migration manifest](2026-11-11-early-benchmark-migration-manifest.csv)
records every tracked source/destination, byte count and SHA-256 before moving:
110 figure files and 157 table files. The sources contain no ignored checkpoint
or other local data object. Each contains one ignored Finder `.DS_Store`, which
were moved unchanged alongside their directories without committing their contents.

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

Two PDAC report sources and the April 9 presentation sources contain affected
relative/constructed readers. Their path roots were updated, including repair
of pre-existing wrong relative depth in the presentation sources. The separate
baseline synthetic report already correctly uses `results/legacy/figures/` and
`results/legacy/tables/`; it was left unchanged after validation distinguished
those artifacts from later rank-selection outputs in the moved directories.
The presentations' baseline synthetic image paths now resolve to that existing
legacy collection. Scientific text and rendered/circulated documents were not
changed or regenerated. The dated inventory keeps its original snapshot paths
and links to this record. No generated scientific output was rewritten.

## Verification and closeout

- All 267 destination files match their pre-move byte counts and SHA-256.
  Their former paths are absent. Both ignored Finder metadata files also match.
- All 5,009 tracked protected files match their baseline hashes. All 3,180
  ignored RDS files retain size, modification time and inode metadata.
- All 91 explicit archived Markdown image references in the six checked report/
  presentation sources resolve. Their prose is unchanged after path normalization.
- Four changed legacy R scripts parse successfully, and both modified shell
  launchers pass `bash -n`. No legacy analysis or SLURM job was executed.
- The tracked circulated HTML files have no external old-root artifact URLs.
  Those HTML files and all circulated PDFs remain unchanged.
- `Rscript tests/run_tests.R` passed 801 tests with zero failures and exit status 0.
  The bounded-β profile-likelihood fallback warning remains visible. The suite
  ran outside the restricted tool sandbox because of its known OpenMP limitation.

The continuation prompt, handoff, current navigation and status pointers now
identify the completed batch. Larger migration and tasks A to F remain separate
work for Andrew's discussion with the next owner. A fresh independent review
found no actionable issues. It verified all 267 files against both the manifest
and the original Git blobs, the protected tracked hashes and RDS metadata,
the runner/source path changes, preservation of circulated documents and the
801-test log. It reviewed references statically and did not rerun a study.
