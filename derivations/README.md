# derivations/

Variational inference derivations for the supervised Bayesian matrix
factorization (SSBMF) models. Each model folder holds one full derivation with
its update algorithm. Supporting pieces sit in a subfolder beneath it.

## Models

| Folder | Model | Contents |
|--------|-------|----------|
| [`SSBMF-LB-model-derivations/`](SSBMF-LB-model-derivations/) | LB, η = Lβ | [`MF_Derivations_UpdateAlgo_REVISED`](SSBMF-LB-model-derivations/MF_Derivations_UpdateAlgo_REVISED.pdf) (.tex/.pdf): full CAVI derivation with the algorithm summary. `modular/` holds the same updates as longer self-contained documents (`qB/`, `qL/`, `qF/`, `qTau/`). `archive/` holds the dated January–March 2026 drafts and the companion for the earlier single-file implementation |
| [`SSBMF-YFB-model-derivations/`](SSBMF-YFB-model-derivations/) | YFB, η = (YF)β | `SSBMF-YFB-model-derivation-full` (.tex/.pdf): full derivation with the algorithm summary, and `verify_derivation.R`. `components/` holds the April 2026 per-update derivations (`qL_unsupervised/`, `qF_supervised/`, `qBeta/`, `ELBO/`), the October 15 canonical-scale F update (`canonical_ZF/`) and the November 5 survival-block weight note (`survival_block_weight/`). The full derivation supersedes them |
| [`SSBMF-YFB-multimodal-model-derivations/`](SSBMF-YFB-multimodal-model-derivations/) | Joint expression + methylation YFB | [derivation](SSBMF-YFB-multimodal-model-derivations/multimodal_YFB_derivation.pdf), implementation specification, extended technical notes, and `verify_derivation.R`, which checks the algebra numerically |
| [`SSBMF-LB-vs-YFB-comparison/`](SSBMF-LB-vs-YFB-comparison/) | Both | Side-by-side comparison of the two parameterizations: where survival enters each update, balance of the two likelihoods, and prediction for new patients |

## Early work

[`early-work/`](early-work/) holds background notes that preceded the
implemented models and were not used for the methodology code.

| Folder | Contents |
|--------|----------|
| `early-work/EBMF/` | Empirical Bayes matrix factorization background |
| `early-work/initial-SSBMF/` | First survival-supervised likelihood derivations, which preceded the CAVI derivations |

## Folder names before 2026-11-06

The circulated 2026-10-16 writeup and older rendered pages may use the old names.

| Old | New |
|-----|-----|
| `MF_UpdateDerivations/` | `SSBMF-LB-model-derivations/` |
| `Modular_UpdateDerivations/` | `SSBMF-LB-model-derivations/modular/` |
| `cox_on_YF/<part>/` | `SSBMF-YFB-model-derivations/components/<part>/` |
| `yfb_survival_block_weight/` | `SSBMF-YFB-model-derivations/components/survival_block_weight/` |
| `multimodal_YFB/` | `SSBMF-YFB-multimodal-model-derivations/` |
| `EBMF/` | `early-work/EBMF/` |
| `SurvivalMF/` | `early-work/initial-SSBMF/` |
