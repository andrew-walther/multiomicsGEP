# Lab Meeting Notes — SSBMF Project (August 27, 2026)

> **Note for Claude Code:** This document mirrors `Meeting_Notes_Followup_8_21_26.md`'s structure —
> a planning-facing transcription of what was asked, not a narrative report (that's
> `docs/progress_book/chapters/2026-09-04.qmd` and `DECISIONS.md`). **Provenance note, unlike the
> 8/21 version:** this was reconstructed on 2026-09-04 from the asks as recorded in the
> post-meeting implementation plan (`we-recently-worked-on-splendid-puzzle.md`), not from a live
> transcript taken during the meeting itself. The direct quotes below are preserved verbatim from
> that plan; anything else is Andrew's own paraphrase of what was asked. If you have the actual
> meeting notes, replace or supplement this document with those — this is a best-effort
> reconstruction, not a primary source.

The meeting went well overall ("results & method are a great basis for the paper" — Naim). Two
items raised were substantive model/method gaps, not just requests for better write-ups, and both
were addressed as real implementation (Stage 1 and Stage 2 of the 9/4 plan) rather than
documentation-only fixes.

---

## ✅ Clarifications (resolved)

1. **"Is the log-likelihood/BIC used to select K_init cross-validated?"** — asked directly. It was
   not; Andrew had told the advisor it was ("You told your advisor [the LL/BIC] were
   cross-validated"), which was incorrect. **Resolved:** a genuine held-out survival log-likelihood
   and a bi-cross-validated genomics log-likelihood now exist alongside the in-sample criteria,
   reported separately (`code/compute_cv_loglik.R`; chapter §1; DECISIONS.md 2026-09-04).
2. **"Most parsimonious within 1 SE of the best"** — the advisors invoked this rule as the form
   K_init selection should take. **Resolved:** this is already implemented, not a new ask —
   `select_K_cv()`'s `use_1se` option, and the joint (K, α) Bayesian-optimization search run in July
   (DECISIONS.md 2026-07-13, Step 4). Andrew's own note at the time, "we already tried this out I
   think?", is correct.
3. **"Are multiple true factors merging into one estimated factor?"** (under-specified K_init) —
   the advisors' question. **Resolved:** yes, directly confirmed with the specific diagnostic
   requested (each estimated factor's top-2 correlations with true factors, not just the
   best-match-only summary) — 7 of 27 estimated factors show a merge signature when K_init is
   under-specified (2, 3, 4); 0 of 114 do at K_init=6, 12, or 20
   (`results/multi_cohort_sim/run_top2_match_diagnostic.R`; chapter §4).
4. **"Are blue/gray applied to the same sim dataset?"** — asked directly, re: the joint-vs-EBMF
   simulation comparison figure. **Resolved:** yes, confirmed directly from the runner (one
   simulated dataset generated per seed, both arms scored against it) and stated in the delta
   figure's caption (chapter §4).
5. **"What is different about K_init=13/14/15?"** (the external-C dip near K_init=11-13) — asked
   directly. **Resolved, but not as expected:** a falsifiable prediction was stated before testing
   it — if the dip is a single-init CAVI local-optimum artifact (the documented failure mode
   elsewhere in this project), a 15-restart best-ELBO multistart search should remove it. **It
   didn't** — the original single-init fit was already the best of all 15 restarts at both K_init=11
   and K_init=13. The dip is real and reproducible; what specifically differs about these two K
   values is still open (`ROADMAP.md`; chapter §1; DECISIONS.md 2026-09-04).
6. **Threshold sensitivity, raised 2026-08-27 (Andrew)**: "the K_init=11/13 dips in the K-sweep...
   are consistent with real boundary sensitivity, not just noise" — i.e., could the dip be an
   artifact of where `beta_thresh`/`pve_thresh` happen to sit, rather than a real property of the
   fit? **Resolved:** no — the retained-factor count stays flat at 2 around both dips under every
   threshold in the stable range (`results/benchmark_sim/run_threshold_sensitivity.R`; chapter §2).
7. **Beta comparability, raised at the meeting**: "does a fixed `beta_threshold` even compare
   factors fairly?" — since it tests raw `|β̂_k|` against one cutoff, assuming every factor's
   projection score is on a comparable scale. **Resolved:** checked directly — `sd(ZF_k)` does vary
   across factors, but a raw vs. variance-standardized `|β_k|·sd(ZF_k)` ranking gives the identical
   ordering for the two survival-active factors in the current fit. The fixed cutoff is defensible
   in practice here (chapter §2; would need rechecking if a future fit shows a much larger spread).
8. **"No elastic-net/L1 penalty anywhere in the model"** — a resolved clarification carried over
   from the plan, not a new implementation item: shrinkage on `beta` is purely the EBNM prior,
   whose variance `ebnm()` estimates from a single `(x_k, s_k)` pair per factor. No separate L1/
   elastic-net regularization exists anywhere in the CAVI loop.

---

## Cohort membership in the survival term (the second substantive gap)

**The ask**: cohort membership currently reaches only the genomics side of the model (`cohort_id`,
fixed indicator columns in L/F) and cannot express differential survival signal at all — can it?

**Resolved**: yes — new `beta_cohort_id` argument (`code/update_beta_cohort.R`), cohort-specific
survival coefficients, distinct from both `cohort_id` and `strata_id`. Full results in chapter §3
and DECISIONS.md 2026-09-04:
- Performance: statistically indistinguishable from the plain model (95% CI on the difference
  includes 0), and beats both the existing `cohort_id` extension and a two-step alternative.
- Interpretability: reveals that one of the two survival-active programs' effect is almost entirely
  specific to one training cohort — invisible to the shared-coefficient model or to external C-index
  alone. Confirmed a second way with simulated ground truth (a factor known by construction to be
  cohort-specific is correctly attributed to the right cohort in a clear majority of cases).
- The manuscript-relevant "using more samples improves factor estimation" claim (asserted in the
  8/21 notes, never tested) was tested directly this session too: pooling clearly helps (mean
  external C 0.627 pooled vs. 0.514/0.578 training on either cohort alone).

---

## Still open (carried to `ROADMAP.md`, not resolved this session)

- **Why K_init=11 and 13 specifically land on a worse solution** — confirmed real and not a
  multistart artifact, but the mechanism is unexplained. Deflation-init or warm-starting from the
  converged K_init=7 solution are the natural next things to try.
- **Literature grounding for the retention thresholds** — the sensitivity range is now resolved
  (item 6 above), but a principled (non-reverse-engineered) threshold value from the sparse-factor
  literature, for the manuscript's methods section, is still open.
- **Whether `alpha_F > 0` (survival feedback into the factorization itself) could be stabilized** —
  raised this session while explaining why this model's kept programs resemble an unsupervised
  method's own highest-variance factors rather than uncovering a lower-variance one: by design, the
  factorization step receives no survival signal at the current default (`alpha_F=0`). Discovering a
  genuinely low-variance, otherwise-overlooked prognostic program would require enabling
  `alpha_F > 0`, disabled by default for a documented CAVI instability. Worth revisiting with this
  session's own stabilization techniques (deflation-init, warm-start), not tested here.
- **Multi-modality integration** — still out of scope: no such data exists under `PDAC_DATA_ROOT`
  (all 7 datasets are gene-level expression/protein), and the current preprocessing row-binds
  patients into one matrix, structurally unable to represent per-modality feature blocks.

---

## Summary for Planning

1. Both substantive gaps raised at the meeting (in-sample vs. cross-validated selection criteria;
   cohort membership limited to the genomics side) were addressed as real model/analysis
   implementation, not documentation — see `docs/progress_book/chapters/2026-09-04.qmd` for the
   narrative and `DECISIONS.md` (six `2026-09-04` entries) for full technical detail.
2. Every direct question recorded from the meeting (items 1-7 above) has a stated resolution.
3. Three items remain genuinely open and are tracked in `ROADMAP.md`: the K_init=11/13 mechanism,
   threshold literature grounding, and whether `alpha_F > 0` can be stabilized.
