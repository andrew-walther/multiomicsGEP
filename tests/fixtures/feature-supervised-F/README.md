# Baseline sources from `feature/supervised-F`

Copies of `code/fit_cox_on_yf.R` and `code/fit_multimodal_yfb.R` as they stood on the local
branch `feature/supervised-F` at commit 2bebbd4 (copied 2026-10-30). Four tests compare the
current fitters' default behaviour with these baselines bit for bit. They used to read the files with
`git show feature/supervised-F:...`, which fails for anyone without that local branch (the pushed
`origin/feature/supervised-F` predates `canonicalize`). Read them with
`read_baseline_source()` in `tests/test_helpers.R`. Do not edit or source these files directly.
