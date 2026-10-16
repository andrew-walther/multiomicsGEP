# paper/writeup/

Project write-up: *Survival-Supervised Bayesian Matrix Factorization for
Prognostic Molecular Programs in Pancreatic Cancer*. It is organized as a
manuscript (Methods, Evaluation, Results, Biological findings, Discussion), with
complete derivations in the appendices.

**Read it:** [`writeup.pdf`](writeup.pdf) or [`writeup.html`](writeup.html).

These are the canonical circulation copies. The [accuracy-review
dispositions](_process/Accuracy_Review_Dispositions_10_16_26.md) distinguish
applied document corrections from donor-aware inference, preprocessing and
provenance work that remains scientifically unresolved.

## Research provenance

The reviewed manuscript and its saved figures are published on `main`.
Its experimental results and accuracy checks refer to research commit
[`d890f7a`](https://github.com/andrew-walther/multiomicsGEP/tree/d890f7abdfab37c388de8edb11e4b7dc312c9649)
on `round2/supervised`. The source inventory points to evidence at that
commit; some diagnostic code, aggregate outputs and plans are present there
but have not been integrated into `main`. The recorded 710-test result also
applies to that research revision. Transferring the manuscript does not
adopt those experimental model changes.

Rendering the manuscript uses the supplied sources and figures and does not
require those experiments. To reproduce historical figure calculations or
diagnostics, consult the pinned research revision, its registered designs
and its data requirements; these are distinct from rendering the writeup.

## Contents

| Path | Contents |
|------|----------|
| `writeup.qmd` | Source of the main text |
| `appendices/single_modality.qmd` | Appendix A: single-modality model derivations and algorithms |
| `appendices/multimodal.qmd` | Appendix B: multimodal (expression + methylation) derivations and algorithms |
| `appendices/diagnostics_and_provenance.qmd` | Appendix C: retained numerical diagnostics, participant flow, public study identifiers and software-provenance limits |
| `references.bib` | Manuscript bibliography; cited entries render with both formats |
| `figures/` | Figures as vector PDF and 300-dpi PNG, with scripts that draw them |
| `_process/` | Planning page, source inventory, review and verification records from drafting; kept as the record behind each revision, not needed to read or build the report |

## Rendering

Run from the repository root. Rendering reads only saved aggregates and figures;
it fits no models.

```bash
quarto render paper/writeup/writeup.qmd --to pdf
```

```bash
quarto render paper/writeup/writeup.qmd --to html
```

If the standalone Quarto reports `Bad CPU type in executable`, use the copy
bundled with RStudio:
`/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`.

## Rebuilding figures

The real-data training/validation flowcharts use audited aggregate cohort
counts and the documented historical pipeline. They read no patient data and
do not fit or score a model. Shape types distinguish cohort inputs, processing,
fitted models and evaluation; dashed callouts explain rank selection and
within-patient matching. Rebuild both diagrams with base R:

```bash
Rscript paper/writeup/figures/make_real_data_workflows.R
```

Figures are drawn from saved fits, scores and aggregates. The survival-curve
figures also need the held-out expression data.

```bash
Rscript paper/writeup/figures/make_pruning_rank_figure.R
```

```bash
PDAC_DATA_ROOT=data/PDAC_data Rscript paper/writeup/figures/make_writeup_figures.R
```
