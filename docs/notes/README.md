# Notes

Derivations and technical notes: mathematical background written to explain or
check the method, not plans or results. Code walkthroughs of the LB updates are in
[`../lb_walkthroughs/`](../lb_walkthroughs/); the current multimodal derivation is in the
repository's `derivations/` folder.

| Date | Document | Purpose | Status |
|---|---|---|---|
| 2026-05-08 | [yfb-derivation-notes](2026-05-08-yfb-derivation-notes.qmd) | Derives YFB (η = YFβ̃) from LB and how each variational update changes | complete (describes the May code, whose F update is genomics-only) |
| 2026-05-08 | [yfb-mathematical-consistency-audit-notes](2026-05-08-yfb-mathematical-consistency-audit-notes.qmd) | Checks nine theoretical YFB claims against the implemented algorithm | complete (May snapshot) |
| 2026-05-08 | [yfb-supplemental-proofs-notes](2026-05-08-yfb-supplemental-proofs-notes.qmd) | Proofs for YFB: train/test consistency, decoupling of L, the dual-source F update, sign invariance | complete |
| 2026-05-08 | [yfb-update-derivations-detailed-notes](2026-05-08-yfb-update-derivations-detailed-notes.qmd) | Algebra-level derivation of the YFB q(β̃), q(L) and q(F) updates, mapped to the code | complete |
| 2026-05-08 | [sign-orientation-issue-notes](2026-05-08-sign-orientation-issue-notes.qmd) | Separates factor-sign indeterminacy, the risk-score orientation convention and the code's sign flip | complete (orientation revisited in DECISIONS 2026-09-04 and 2026-10-02) |
| 2026-05-30 | [pathway-enrichment-overview-notes](2026-05-30-pathway-enrichment-overview-notes.qmd) | What pathway enrichment is and what to expect, illustrated with D4 weights, before any analysis ran | superseded by [the 2026-07-15 enrichment report](../findings/2026-07-15-pathway-enrichment-report.qmd) |
