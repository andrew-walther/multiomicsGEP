# November 13 writeup

This draft presents the rank-5 supervised YFB working model and its training and external results, including the rank-matched two-step comparison, risk groups and program exhibits. The prior and input-scaling comparisons use a rank-7 reference fit of the same configuration, and the simulation studies fit rank 7. The writeup labels each rank-7 result. The draft also reports methylation alone and a descriptive matched-TCGA feature-group comparison, which is **not validated by simulation**. The previous [November 6 version](../2026-11-06-writeup/README.md) remains the circulated record.

The [source](writeup.qmd) reads saved aggregate CSVs. The [number and figure ledger](NUMBER-SOURCES.md) lists the evidence. Rendering does not require patient records, fitted objects or a new analysis. From the repository root, with Quarto, R, knitr and LuaLaTeX installed, run:

```bash
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/2026-11-13-writeup/writeup.qmd --to html
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render paper/writeup/2026-11-13-writeup/writeup.qmd --to pdf
```

The PDF uses Chapter 4 vector figures. HTML uses SVG conversions in `figures/`. To regenerate the SVGs after a Chapter 4 figure changes:

```bash
for svg in paper/writeup/2026-11-13-writeup/figures/sup_*.svg; do
  stem=$(basename "$svg" .svg)
  pdftocairo -svg "paper/prelim/figures/$stem.pdf" "$svg"
done
```

The source contains methods, evaluation, current analyses, results, program descriptions and discussion. Appendices document the derivation, changes from the prior version, comparator details and prior sensitivity. The HTML and PDF are draft outputs.

The PDF contains 47 pages. The source evaluates its aggregate checks
during both renders. Ten Chapter 4 figures are reused, alongside two plots
generated from aggregate results.
