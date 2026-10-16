# paper/writeup/

Project write-up: *Survival-Supervised Bayesian Matrix Factorization for
Prognostic Molecular Programs in Pancreatic Cancer*. It is organized as a
manuscript (Methods, Evaluation, Results, Biological findings, Discussion), with
complete derivations in the appendices.

**Read it:** [`writeup.pdf`](writeup.pdf) or [`writeup.html`](writeup.html).

## Contents

| Path | Contents |
|------|----------|
| `writeup.qmd` | Source of the main text |
| `appendices/single_modality.qmd` | Appendix A: single-modality model derivations and algorithms |
| `appendices/multimodal.qmd` | Appendix B: multimodal (expression + methylation) derivations and algorithms |
| `figures/` | Figures as vector PDF and 300-dpi PNG, with the two scripts that draw them |
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

Figures are drawn from saved fits, scores and aggregates. The survival-curve
figures also need the held-out expression data.

```bash
Rscript paper/writeup/figures/make_pruning_rank_figure.R
```

```bash
PDAC_DATA_ROOT=data/PDAC_data Rscript paper/writeup/figures/make_writeup_figures.R
```
