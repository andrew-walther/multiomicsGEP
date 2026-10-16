# derivations/

Variational inference derivations for the model. The complete current
derivations, with algorithms, are also in the appendices of the
[project write-up](../paper/writeup/writeup.pdf).

## Current

| Folder | Contents |
|--------|----------|
| [`MF_UpdateDerivations/`](MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf) | Corrected CAVI derivations for the LB model (η = Lβ): `MF_Derivations_UpdateAlgo_REVISED` (.tex/.pdf) |
| `Modular_UpdateDerivations/` | One derivation per LB update: `qB/`, `qL/`, `qF/`, `qTau/` |
| `cox_on_YF/` | YFB model (η = (YF)β): q(L), q(F) with survival supervision, q(β), the ELBO, and the canonical-scale F update |
| `multimodal_YFB/` | Joint expression + methylation YFB model: [derivation](multimodal_YFB/multimodal_YFB_derivation.pdf), implementation specification, extended technical notes, and `verify_derivation.R`, which checks the algebra numerically |

## Background and early work

| Folder | Contents |
|--------|----------|
| `EBMF/` | Empirical Bayes matrix factorization background |
| `SurvivalMF/` | Early survival-supervised likelihood derivations that preceded the CAVI derivations above |
| `MF_UpdateDerivations/archive/` | Dated drafts of the CAVI derivations (January–March 2026) and the companion that maps them to the earlier single-file implementation; superseded by `_REVISED` |
