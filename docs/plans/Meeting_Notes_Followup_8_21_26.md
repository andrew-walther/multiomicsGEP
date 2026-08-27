# Lab Meeting Notes — SSBMF Project (August 2026)

> **Note for Claude Code:** This document is a reference for planning future improvements to the SSBMF / multiomicsGEP project. The **"Preparation for 8/27 Rashid Lab Meeting"** section below is **TOP PRIORITY** — it directly drives the slide deck due for the 8/27/2026 lab meeting and should be addressed first, ahead of the other backlog items in this document.

---

## 🔴 TOP PRIORITY: Preparation for 8/27 Rashid Lab Meeting

### Setup
- Use the `presentation/` directory. Make a copy of the `walther_lab_meeting_06_18_2026` presentation folder.
  - Review this presentation and identify relevant updates made since then.
  - Create a copy called `walther_lab_meeting_08_27_2026` to produce the new slide deck.
- Structure:
  - Brief introduction of the problem, method, and model.
  - Review the Quarto book created previously and pull figures from the 8/3 meeting update, specifically:
    1. Figure showing the genes included in each selected factor.
    2. Figure showing curves as survival signal increases from 0.

### New Content to Present

**1. K selection (ARD shrinkage prior)**
- Show the "curve" over a range of K, across multiple selection metrics.
- Make **4 plots**, one for each of: ELBO, BIC, Log-Likelihood, and external C-index. **Range: start at K=2 (fall back to K=4 if K=2/3 are unstable), through K=20** — see [Clarifications](#-clarifications-resolved).
- Frame these 4 plots as a **consensus recommendation** for K_init — ELBO/BIC/Log-Likelihood are model-fit metrics (largely interchangeable), C-index is prediction-focused (external cohort performance); the deck should highlight where they broadly agree, not treat each as an independent selector. See [K Selection & Active Factor Framework](#k-selection--active-factor-framework).

**2. Joint vs. 2-step method comparison**
- Show a figure illustrating the progression of joint vs. 2-step methods as survival strength increases.

**3. Selected factor gene composition**
- Show a figure illustrating key genes included in each of the selected factors (genomics- and survival-active factors).

### Next Steps Section (for the deck)
1. **Value-add focus** — importance of prioritizing either:
   - Integrating multiple cohorts (cohort indicator included in the L matrix), and/or
   - Integrating multiple modalities (modality indicator included in the L matrix).
2. **Biological characterization** — illustrate the procedure/methods for characterizing the biological makeup and relevance of the selected factors (the ~5 candidate methods are already documented in project context — pull from there).

### Key Takeaways / Framing for This Deck
- Minor revisions from the 6/18 deck: specify **ARD-based K selection** instead of K selection via cross-validation on the external C-index.
- Emphasize new results:
  - Gene lists for selected factors (adverse vs. protective).
  - Joint method and strict EBMF show **equivalent performance when there is no survival signal**.
  - Joint method **outperforms** the 2-step method (EBMF + Cox).
- Close with forward-looking direction:
  - Understanding biology of selected factors.
  - Integrating multiple cohorts and/or modalities.

---

## ✅ Clarifications (resolved)

1. **K sweep range**: Start at **K=2**, unless K=2 and/or K=3 show instability — in that case, start at **K=4**. Sweep up to **K=20** either way. (Instability at low K was observed previously, so K=4 is the likely practical starting point, but K=2 should be attempted first for completeness.)
2. **"YFB"** refers to a specific model parameterization — defined in the project's existing model/code documentation, not an undefined term. Claude Code should pull the exact definition from the project context (e.g., `CLAUDE.md` / model spec) rather than treat this as ambiguous.
3. **"5 methods" for biological characterization**: already documented in the project context (`CLAUDE.md` or prior notes) — Claude Code should reference that rather than re-derive or guess.
4. **Items marked "?"** (uniform prior on K, extending to multiple modalities, modality-specific priors, etc.) are **open methodology considerations** — options to evaluate empirically and adopt only if they demonstrably improve results, not committed decisions.
5. **Key Results Needed for Manuscript**: these are the results intended to demonstrate the "value" / contribution of the SSBMF method in the final manuscript — i.e., the core evidence for why the joint method matters, not incidental analyses.

**Still open:** *(resolved below — see "K Selection & Active Factor Framework")*

### K Selection & Active Factor Framework
*(Added based on 8/26 clarification — this explains how the metrics relate and fits together as a two-stage process.)*

**Stage 1 — Selecting K_init (consensus across metrics):**
- **BIC, ELBO, and Log-Likelihood** are **model-focused** metrics — they assess model fit/parsimony and are broadly interchangeable with one another (any one, or a combination, can serve this role).
- **External C-index** is **prediction-focused** — it measures the model's predictive performance on external validation cohorts, not model fit per se.
- These four metrics together are not meant to each independently "pick" K — rather, they form a **consensus recommendation**: look across all four curves and identify where they broadly agree on a plateau/optimal K_init, rather than relying on a single metric in isolation.

**Stage 2 — Selecting the active factor subset (after K_init is chosen):**
- Once K_init is fixed, the model still needs to determine which of the K factors are "active" (i.e., retained as meaningful).
- A factor is considered active if it satisfies **either** of:
  1. It explains sufficient variance in the genomics/expression matrix (unsupervised relevance), **or**
  2. It explains the survival signal (supervised relevance via the Cox/survival component).
- This means a factor doesn't need to contribute to survival prediction to be kept — explaining expression variance alone is sufficient, and vice versa.

This two-stage framing (consensus K_init selection → active factor subset selection) should carry through to both the 8/27 deck and any implementation/documentation of the K-selection procedure.

---

## Full Meeting Notes (Backlog / Future Improvements)

### K Selection
- Method update: cross-validation has been superseded by **ARD shrinkage** as the primary selection method (both methods can still be compared).
  - **CV**: C-index (external cohorts).
  - **ARD**: uses ELBO (also consider BIC).
- *(Open consideration)* Uniform prior on K — evaluate against literature precedent; adopt if it improves K selection.
  - Discrete — check literature for precedent.
- Choosing a starting value of K is a hard question.
  - Idea: run EBMF / consensus NMF first to get a reasonable starting K.
- Current setup: pivot from CV to ARD shrinkage prior to select "K effective" from a larger starting K value.
- **Follow-up task**: Starting K sweep should begin at **K=2** (fall back to **K=4** if K=2/3 are unstable), through **K=20** — see [Clarifications](#-clarifications-resolved). Compute an "objective" metric (ELBO / BIC / Log-Likelihood / C-index) for each starting K. Build a figure showing the "curve" to identify a performance "plateau" — select the smallest starting K value beyond which additional info leads to overfitting or no further improvement.

### Performance Metrics
Used together as a **consensus** to select K (see [K Selection & Active Factor Framework](#k-selection--active-factor-framework) for how these combine — BIC/ELBO/Log-Likelihood are model-fit metrics and largely interchangeable; C-index is prediction-focused):
- BIC
- Log-likelihood
- C-index

### Multiple Random Initializations
- Sweep K from **2 to 20** (true factors), falling back to starting at K=4 if K=2/3 prove unstable — see [Clarifications](#-clarifications-resolved). See where the model converges.
- For each K, evaluate BIC/CV performance of the survival model — as one input into the consensus K_init recommendation (see [K Selection & Active Factor Framework](#k-selection--active-factor-framework)).
- Note: increasing K plateaued after 3 for Amber's data.

### Factor Normalization
- Normalize factors before comparing Beta coefficients — need to be on the same scale to compare across factors. Refers to the **YFB model parameterization** (defined in project context/model docs).

### Output Interpretation
*(Related to Stage 2 of the [K Selection & Active Factor Framework](#k-selection--active-factor-framework) — a factor can be retained for explaining expression variance alone, survival signal alone, or both.)*
- If a non-survival factor still has a small survival regression coefficient, that's acceptable — its impact will be minimal if any.
- Trust factors with larger coefficients.
- Factors with small coefficients have small/no impact — not a concern.

### Next Focus Areas
- Streamline implementation of K selection.
- *(Open consideration)* Extend to multiple modalities — evaluate and adopt if it improves results.
- *(Open consideration)* Cohort membership dummy variable:
  - Indicator competes to explain variation in Y — needs a penalty.
  - One cohort: variability due to biology alone.
  - Two cohorts: variability due to biology, batch effect, and randomness.
- Multiple modalities *(open consideration, contingent on the above)*:
  - Stack methylation and RNA expression in the feature matrix.
  - *(Open consideration)* Different priors for methylation vs. gene expression components — evaluate empirically.
  - Add additional columns per modality.
- For integrating multiple cohorts — expected results to demonstrate:
  - Consensus learned factors are more reproducible in terms of association with survival.
  - Leave-one-study-out (study-aware) validation.
  - Using more samples improves factor estimation.
  - Train on cohort A or B → validate.
  - Train on cohort A & B → validate (expect better performance vs. single-cohort training).

### Key Results Needed for Manuscript
*(These are the core results intended to demonstrate the "value"/contribution of the SSBMF method in the final manuscript.)*
1. Comparison of SSBMF vs. EBMF for recovering gene programs.
2. Comparison of SSBMF vs. 2-step (EBMF + Cox) for external validation C-index.

---

## Summary for Planning
When building the improvement plan for this project:
1. **Immediate priority**: execute the 8/27 lab meeting presentation tasks above (deck copy, 4 K-selection plots, joint-vs-2-step figure, gene composition figure, next-steps slide).
2. **Secondary priority**: address the broader backlog — K selection methodology refinement, multi-cohort/multi-modality extensions, manuscript result generation (SSBMF vs. EBMF, SSBMF vs. 2-step C-index comparisons).
