# Round2 detailed traces

The 15 detailed iteration CSVs from `diagnosis1` and `framework_fold_trace`
are preserved at Git revision
`3e9fecc99c72065da91994986eb1eb8c4d214dd6`. The summary and comparison CSVs
remain at their existing paths. See `TRACE_SHA256SUMS` for the original
paths and content hashes.

From the repository root, restore the archive into a separate directory.

```bash
archive_dir=$(mktemp -d)
cut -c 67- results/multimodal_real/outputs/round2/TRACE_SHA256SUMS |
while IFS= read -r trace_file; do
  mkdir -p "$archive_dir/$(dirname "$trace_file")"
  git show "3e9fecc99c72065da91994986eb1eb8c4d214dd6:$trace_file" > "$archive_dir/$trace_file" || exit 1
done
manifest="$PWD/results/multimodal_real/outputs/round2/TRACE_SHA256SUMS"
(cd "$archive_dir" && shasum -a 256 -c "$manifest")
```

The diagnostic runners construct these trace paths when writing new runs.
The framework replay reads `diagnosis1/comparison_all.csv`, its saved
framework fit, the fixed-prior cell and the restricted local training cache.
The fixed-prior runner also reads `comparison_all.csv`. The October 15
notebook refers to the detailed trace names as provenance. It does not read
them when rendering. Full replay performs a new fit and is separate from
restoring historical trace bytes.
