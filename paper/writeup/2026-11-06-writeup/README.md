# paper/writeup/2026-11-06-writeup/

Version of the project writeup dated November 6, 2026. **In progress (Phase 1 done, paused).**
The previous version is [`../2026-10-16-writeup/`](../2026-10-16-writeup/README.md), which is not
edited.

This version starts from a new outline, not from the October 16 source. Its centerpiece is the
survival-weighted YFB fit (`fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1)`, 2,218 D4 genes, single
log). It leads with the comparison against the same factorization at rank 7 followed by Cox,
reports the pre-registered October 15 null result, and lists every use of the external cohorts
(GitHub issue #6 items 2, 4, 5 and 6). The section plan, with the evidence file for every number,
is [`_process/plan.md`](_process/plan.md).

## Reproducibility

`writeup.qmd` runs R. Every number is computed from saved CSVs under
`results/benchmark_sim/outputs/` (issue #6 item 8); a missing file stops the render. The
registered October 15 values are read from Git history (commit `65cf90a`). Exceptions, cited to
their source files in the text: the PACA-AU donor counts (47 shared, 568 people) and the
November 3 joint-versus-EBMF difference, which has no saved table.

```bash
quarto render paper/writeup/2026-11-06-writeup/writeup.qmd --to html
```

## Where it stopped (2026-11-06, about 04:40)

Done: section plan; `writeup.qmd` with CSV-driven chunks; drafts of Methods, Evaluation
(including the external-cohort uses table), Results 3.1 to 3.5 and 3.7 (tables only), Programs
table, Discussion draft, Appendices A and B (short). HTML renders cleanly. PDF not yet rendered.

Next (Phase 2, after the prelim Chapter 4 revamp lands):
1. Add Chapter 4 figures (simulation weight curves, S-14 paired C, external C by series, program
   heatmap and correspondence, KM if produced) and the program biology.
2. Add §3.6 from the S-15 findings and CSVs.
3. Finish the Summary and Limitations; one-page algorithm in Appendix A.
4. Render HTML and PDF; human-writer review by independent reviewer subagents; number-trace check.
5. Update `../README.md` status and PROJECT_STATUS.qmd.
