# Round 2 attempt history

Authority: `docs/plans/Codex_Handoff_Round2_10_16_26.md` and Andrew's request.
Log each attempt here and commit its design before running it. Append outcomes;
never remove or rewrite a failed attempt. Defaults, thresholds and test beds stay
fixed. Selection uses simulation and training data only. New real-data tables
must include the immutable benchmark snapshot rows.

## Benchmark and integration record

Before Round 2, all ten SHA256 entries passed in
`results/benchmark_snapshot_2026-10-14/`. The annotated benchmark tag resolves to
`a352f17ca9cff523caa80756aa2af085e62551e1`; neither tag nor snapshot was changed.
Created `round2/supervised` from clean `feature/supervised-F` at `2bebbd4`.
Task 0 integration completed in the required order with `--no-ff` merges,
preserving all entries and tests. Successive full-suite counts were 553, 554,
556, 559 and 559, all with zero failures. No scientific attempt has run yet.
