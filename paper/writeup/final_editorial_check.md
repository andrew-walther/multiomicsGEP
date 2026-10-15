# Final editorial check

**Checked:** 2026-10-14 on `writeup/draft` at `d80b724`; clean tree, 2 commits ahead of origin. Sources read: `writeup.qmd` (Summary through Continuing work, in full), both appendices where they changed, `verification_round2.md`, and the "Second verification" section of `review_dispositions.md`. Rendered PDF: 39 pages; main text ends on p.24, Appendix A starts on p.25. Nothing was edited, rendered, run, committed or pushed for this check.

## 1. N1 and N2

| Item | Verdict | Evidence |
|---|---|---|
| N1 Ridge in the first sweep-normalizer formulas | **Correctly applied** | `writeup.qmd` §1.2 and `appendices/single_modality.qmd` A.1 both give $d_k^{(r)} = (\lVert\bar F^{(r-1)}_{\cdot k}\rVert_2^2 + 10^{-10})^{1/2}$. This matches `code/fit_cox_on_yf.R:643–646` (`sqrt(colSums(EF_global^2) + 1e-10)`) and A.10. The random target's zero-direction convention is kept separate. Rendered: the formula appears in the PDF on p.3 and p.25, and three times in the HTML. |
| N2 Disposition-file pointer | **Correctly applied** | `PROJECT_STATUS.qmd:141` now points to `paper/writeup/review_dispositions.md`. |

## 2–3. Coherence, contradictions, redundancy, undefined terms, references

The report reads coherently in manuscript order (Methods → Evaluation → Results → Biology → Discussion). No numerical contradiction was found between text, tables, captions and Summary. All `@sec/@tbl/@fig/@eq` references and appendix links resolve; a scripted check found no dangling labels.

### Required corrections (small; no analysis needed)

| # | Location | Problem | Fix |
|---|---|---|---|
| E1 | Table 7 caption | The row is labelled "Single-block model", but the caption says "the single-modality model is a single-block fit…". The label and caption disagree. | Say "the single-block model is a fit of the model of @sec-single-model to unstandardized log2 expression with at most 100 sweeps". |
| E2 | §2.3 Comparators, "Multimodal comparators" | It lists the expression-only fit and "three further arms", but Table 7 also reports the base multimodal model (0.453), which §2.3 never introduces. It is defined only in §2.2's simulation design. | Add the base model: "Four further arms are reported: the base multimodal model without intercepts or factor removal, a single-block fit …, and two-step EBMF …". |
| E3 | Continuing work, item 1 | "Survival has little effect even when it enters the loading update" is stated as a finding, but the report contains no evaluated fit with survival in the single-modality loading update. The closest evidence is the multimodal framework's failure to recover the low-variance program. | Base it on reported evidence: "In the multimodal framework, where survival does enter the loading update, the planted low-variance prognostic program was not recovered. Plausible reasons are that the survival contribution is multiplied by the current coefficient…; that its precision is small relative to reconstruction; and that the single-modality loading scale is not fixed during the update." |
| E4 | Uncited exhibits | Table 1 is never cited in the text. Tables 10 (A.11) and 11 (B.13) are not cited in their appendices. A manuscript cites every exhibit. | Notation: "@tbl-models summarizes how the two models differ." A.11: "@tbl-math-vs-impl compares…". B.13: "@tbl-mm-reconcile lists…". |
| E5 | First prose use of abbreviations | ICGC (first in §2.1, "TCGA and ICGC patients"), CpG (§2.1) and PACA-AU (first in §3.2 text) are defined only inside table captions. Captions are self-contained, but the prose should also define a term at its first use. | §2.1: "International Cancer Genome Consortium (ICGC) patients…" and "CpG (cytosine–guanine dinucleotide) sites". §3.2 or §2.1: "the Australian ICGC cohort (PACA-AU)". |

### Optional editorial suggestions

- **O1 Characterization methods next to the findings.** §2.5 sits about ten pages before §4. Moving it to open §4 (as "Characterization approach", with §4.1 and §4.2 following) would put each method next to the result it supports. This conflicts with author item B11, which asked for it in Evaluation, so it is Andrew's call. It saves no space. A lighter alternative is a single line at the top of §4 pointing to §2.5.
- **O2 Cohort-specific coefficient simulation.** It appears only in §3.1 and is not described in §2.2. Either add one clause to §2.2 ("a separate ten-dataset variant gave one program a survival effect in one cohort only") or drop the sentence; it is peripheral to the evaluated configuration.
- **O3 Duplicate wording.**
  - The within-cohort standardization appears in both §1.4 and §2.1. Keep §1.4 and refer to it.
  - The Table 8 caption gives the validation sizes twice ("($n = 50$)…($n = 67$)" and "Validation sets: 50…(30 deaths)…"). Keep one.
  - §1.6 Approximations and the first Limitations paragraph overlap; see M4.
- **O4 Minor wording.**
  - §1.2 "onto its loading columns": "its" can be misread as the patient's. Use "onto the loading columns".
  - §3.1 "(@tbl-single-sim; rank 12 gave 0.98 and 1.00)": say which scenario each number belongs to (all shared; hybrid).

## 4. Methods cuts toward 15–20 pages

The main text is 24 pages, and Methods occupies pp.2–7. Part of that is white space before each unbreakable algorithm: p.3 and p.6 are about two-thirds full. The estimates below are approximate, because floats and table reserves will move; only a re-render can confirm them. Every cut preserves the scientific qualifications. The comparator rows of Table 7 are retained as instructed.

| Rank | Cut | Moves to | What stays in Methods | Approx. saving |
|---|---|---|---|---|
| M1 | Move Algorithms 1 and 2 out of the main text | Start of Appendix A and Appendix B, as "A.0/B.0 Algorithm", keeping their Revisit markers | A 3–4 sentence update order per model (initialize; per sweep: coefficients, activities, loadings, precisions; stop rule; risk direction). The evaluated settings sentence. | **1.5–2 pages**, including the two white-space gaps the minipages force |
| M2 | Condense the single-modality "Survival predictor" and "Cox working approximation" paragraphs | Sweep-fixed $D^{(r)}$, ridge, zero-column convention and $D_*$ → A.1/A.10. @eq-cox-working and the Poisson-style explanation → A.2 | @eq-zf. "Prediction uses the loading norms stored with the fit." "Coefficients are updated with a diagonal working approximation to the Cox partial likelihood, not its second-order Taylor expansion (Appendix A.2)." | 0.5 page |
| M3 | Condense the multimodal "Updates" and "Factor removal" paragraphs | The @eq-mm-elbo display and its timing detail → B.11; retarget the @eq-mm-elbo reference in §3.5 | In words: removal uses a score with a plug-in survival term and KL penalties cached against earlier priors; deletion only when the score does not decrease; runs on unconverged fits; $K > 1$; no guarantee for the exact evidence | 0.3–0.4 page |
| M4 | Fold §1.6 Approximations into Limitations and the appendices | The Jensen/plug-in sentence → B.11 (already there); "not exact coordinate ascent on one ELBO" → Limitations | Nothing separate; Limitations paragraph 1 already lists the same approximations | 0.3 page |
| M5 | Move the cross-validation numbers for the starting rank out of §1.3 | §3.3 (Results), which already discusses the plateau | Definitions of starting rank and classes; "K_init = 7 was chosen by cross-validation and external concordance (@sec-results-rank)" | 0.15 page |
| — | Keep Table 1 | — | Table 1 becomes the compact overview that lets M2 and M3 shorten the prose; cite it (E4) | — |

**Expected result:** about 2.7–3.3 pages saved, taking the main text to roughly 21 pages. Methods cuts alone are unlikely to reach 20. The next-lowest-loss options outside Methods are:
- condense §2.2 Simulation designs to one paragraph per model and move the remaining design detail into the appendices (about 0.5 page);
- move the starting-rank sensitivity sentence and the cohort-specific simulation (O2) into notes (about 0.2 page).

After any cut, re-measure the table reserves; Tables 7, 8 and 9 and Figures 3–5 are pinned or reserved, and their positions will change.

## Not checked

- Appendix algebra beyond the N1 change; it was verified in `verification_round2.md`.
- HTML at phone width.
- Numerical values against saved outputs; no numbers changed after `fbf8317`.
