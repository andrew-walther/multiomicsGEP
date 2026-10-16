# Writeup publication on main — 2026-10-16

Andrew requested that the most current reviewed writeup live on `main`.
The complete manuscript bundle is transferred from research commit
`d890f7abdfab37c388de8edb11e4b7dc312c9649`: source, appendices, bibliography,
saved figures, figure scripts, rendered PDF/HTML and drafting/review records.
The root README, roadmap, decisions and project-status source receive only
the related reporting corrections and future risk-classification question.
Main's existing restart plans and model implementation are preserved.

Research code, diagnostic aggregate outputs, patient data and the separately
modified advisor-meeting PDF are outside this transfer. Source-inventory
paths and historical accuracy checks refer to the pinned research commit;
the writeup README documents how to locate them. No experiment is adopted
or executed by publishing this manuscript.

## Verification

- HTML and PDF render successfully from main's manuscript bundle.
- The main text and all three appendix sources match the pinned revision.
- The rebuilt 51-page PDF has exactly the same extracted text as that
  revision's approved PDF.
- All HTML fragment links resolve and all 24 bibliography entries appear.
- Both figure scripts parse successfully. They are copied unchanged from
  the tested research revision; no model or projection is rerun.
- The historical 710-test result belongs to the research revision, not a
  newly executed main-branch model test suite.
- Whitespace checks pass with the source inventory's existing CRLF line
  endings recognized.
