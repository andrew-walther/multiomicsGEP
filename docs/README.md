# docs/

| Path | Contents |
|------|----------|
| [`progress_book/`](https://andrew-walther.github.io/multiomicsGEP/) | Progress notebook, one chapter per advisor meeting (rendered to GitHub Pages) |
| [`plans/`](plans/README.md) | Current plans and status; superseded plans in `plans/archive/` |
| `reports/` | Dated analysis reports and advisor summaries (`.qmd` with rendered `.pdf`/`.html`) |
| `lb_walkthroughs/` | Walkthroughs of each LB update (`update_beta`, `update_L`, `update_F`, `update_tau`) and the LB fitting loop (`fit_modular`), as `.qmd`/`.pdf`/`.html`, plus [`update_L_fix.md`](lb_walkthroughs/update_L_fix.md), the debugging guide for the L update |
| `PDAC_data_audit` (`.qmd`/`.pdf`/`.html`) | Audit of the seven PDAC expression cohorts |
| `notes/` | YFB derivation notes and consistency audits (May 2026); external reviews in `notes/gemini/` |
| `progress_report/` | Earlier status updates (April–August 2026) |
| `reviews/` | Reviews of the progress notebook and prelim chapter drafts |

`make -C docs` re-renders the walkthroughs and the data audit.
