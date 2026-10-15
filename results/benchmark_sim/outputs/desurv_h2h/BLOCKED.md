# DeSurv head-to-head: stopped before fitting (2026-10-15)

No usable DeSurv implementation is installed or vendored in the searched locations.
`requireNamespace("DeSurv", quietly=TRUE)` and the lowercase spelling both return FALSE.
The only configured R library is `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library`.
The checkout (including ignored `data/`), its local `data/multiomicsGEP_code` scripts, and
PDAC_data at the configured OneDrive root contain no DeSurv package/fitting implementation.
`code/preprocess_desurv.R` supplies preprocessing, not the DeSurv method.

The authors' software page links to [rashidlab/DeSurv](https://github.com/rashidlab/DeSurv)
and [DeSurv-paper](https://github.com/rashidlab/DeSurv-paper). The latter documents the
companion package, installed using `remotes::install_github("rashidlab/DeSurv")`, and
Bayesian optimization for model tuning. The package repository could not be fetched through
the browser. No package was installed and no method was reimplemented, as required by Prompt A.

Needed to resume: a working, version-pinned DeSurv installation or local source, its documented
rank/penalty tuning and prediction interfaces, and confirmation that its nonnegative-input
requirements permit the requested D4 per-platform z-standardized matrices. Do not silently
shift or truncate those centered matrices. The citation in `paper/ssbmf-refs.bib` remains
`@unpublished{Young2026}`; the replication README identifies PNAS materials but this audit
does not independently establish publication metadata.

No DeSurv fit, supervised refit, external scoring, bootstrap, C table, or head-to-head runner
was created. The prompt stopped at its explicitly stated implementation prerequisite.

Configuration carried from Prompt E (`m5/k-selection`, ac572f6): protected pruning not adopted;
starting K remains 7 and the handoff single-modality configuration is retained without pruning.
