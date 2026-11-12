# paper/writeup/2026-11-13-writeup/

Version of the project writeup dated November 13, 2026. **Draft. Results of the November 9–13
tasks pending.** The previous version is [`../2026-11-06-writeup/`](../2026-11-06-writeup/README.md),
which is not edited.

This version presents the current model, the supervised YFB fit ($w_S$ = 1, 2,218 genes,
rank 7), and compares it only with external benchmarks: the rank-matched two-step fit (EBMF at
$w_S$ = 0, then Cox), EBMF at ranks 15 and 40 followed by Cox, DeSurv, elastic-net Cox and
PurIST. It states the rule that the external cohorts only score a model fixed on training data,
names every model by its configuration, and adds a section with marked placeholders for the
analyses of November 9–13. Comparisons with earlier versions of our own model (the October 15
supervision check, the $\alpha_F$ = 0 benchmark and its single-log refit, the 40-factor LASSO
comparator, historical loading matches, the table of past external-cohort uses and the earlier
multimodal appendix) are removed. No retained number from the November 6 version has changed.
The [initial plan](_process/plan.md) describes the November 6 version and is kept for reference.

## Reproducibility

The [source](writeup.qmd) evaluates the saved aggregate CSVs and stops if a
required input is missing. The aggregate PACA-AU donor audit
transcribes the October 16 [accuracy dispositions](../2026-10-16-writeup/_process/Accuracy_Review_Dispositions_10_16_26.md).
The [number and figure ledger](NUMBER-SOURCES.md) identifies the input files
for every table, figure and result in the text.

From the repository root, with Quarto, R, knitr and a working LuaLaTeX
installation, run:

```bash
quarto render paper/writeup/2026-11-13-writeup/writeup.qmd --to html
quarto render paper/writeup/2026-11-13-writeup/writeup.qmd --to pdf
```

The PDF reuses the Chapter 4 vector figures directly. The HTML embeds SVG
conversions of the same files. To regenerate the SVGs after a Chapter 4
figure changes, run the command below from the repository root.

```bash
for svg in paper/writeup/2026-11-13-writeup/figures/sup_*.svg; do
  stem=$(basename "$svg" .svg)
  pdftocairo -svg "paper/prelim/figures/$stem.pdf" "$svg"
done
```

The writeup needs the committed aggregate evidence and figures. Rendering
does not need patient records, fitted objects or a new analysis. The
[HTML](writeup.html) and [PDF](writeup.pdf) of this version are drafts. The
November 6 version remains the circulated record.

## Contents

| Part | Content |
|---|---|
| Summary | Aim, model, training and external comparison with the two-step fit, other comparators, simulation results |
| 1 Methods | Notation, model, fitting and prediction, the two-step comparator, choice of the survival weight, the model behind the external results |
| 2 Evaluation | Data, uses of the external cohorts, simulation designs, metrics |
| 3 Work of November 9–13 | Registered analyses of the week with result placeholders A1 to F4 and EXT |
| 4 Results | Supervised fit against the rank-7 two-step fit, real-data loadings, the three simulation studies, the weight refit, pruning, comparison with other prognostic models, risk groups, program stability, mixed-signal and cross-design simulations |
| 5 Programs | Programs of the supervised fit, enrichment and published signatures, individual programs, external associations and Kaplan–Meier curves |
| 6 Discussion | Interpretation, limitations, continuing work |
| Appendix A | Summary of the model derivation |
| Appendix B | Changes since the November 6 version |
| Appendix C | Comparator calibration by series and numbers at risk |

## Manuscript-planning extension

The expanded source reports the fixed prognostic comparisons, individual
program evidence, and the training-cutpoint risk groups. Appendix C gives
the per-series calibration and risk counts. The two new plots are generated
by `ggplot2` from saved aggregate CSVs at render time. The rebuild commands
above also apply to the expanded source.
