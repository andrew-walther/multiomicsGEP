# Rashid Lab Meeting Notes — 6/18/2026
**Topic:** SSBMF Feedback & Next Steps

---

## Model Specification Questions

### Baseline Hazard — h₀(t)
- How is h₀(t) specified? It should be.
  - Question raised: isn't this canceled out in the Cox partial likelihood?

### Mixing Parameter (Genomics ↔ Survival)
- Mixing between genomics fit and survival fit should always default to 1/2 (α).

### Normalization w.r.t. n & p
- **Normalize with respect to n and p:**
  - Divide expression term by *np*
  - Divide survival term by *n*
- **Purpose:** addresses the expression component dominating the objective function relative to the survival component.

### Lambda (λ) Scaling
- Why isn't λ scaled at all?
  - Current scale for λ is written as (0, ∞) — should be constrained to (0, 1)
  - Add a (1 − λ) term in front of L_gen

### CAVI
- CAVI = block coordinate descent (definitional note, no action item)

---

## Preprocessing

- Z-standardize (Z-std) input data
- **Platform correction:** must correct for platform *after* Z-std
  - L is sample-level
  - Consider handling in the Cox model instead — different studies have different baseline survival rates
    - Add a strata variable to the linear predictor to fit a different baseline hazard per dataset
    - In R: `+ strata(study)` (study-specific variable)

---

## Factor Number (k)

- Current model needs 7 factors — why so many?
  - DeSurv was able to reduce down to 3 factors
- Question: if survival component is turned off, does k shrink?

---

## Synthetic Data / Simulation Diagnostics

- Differences in synthetic data results could stem from how the data is simulated (tunable).
- **Key diagnostic logic:**
  - If survival has no impact, there should be no difference between our joint model and a 2-step model
  - We should see the same C-index between joint and 2-step models if survival has no impact
- **Proposed checks:**
  - Constrain λ to zero
  - Pick 1–2 factors related to survival (range 0 → larger value) and check for differences in fit
  - Review how Amber simulated her data

---

## Factor Interpretation

- **Adverse factor:** basal-like
- **Protective factor:** see Amber's work (DeSurv paper)

---

## External Dataset Normalization

- Are external datasets normalized?
  - Dijk / Moffitt / PACA-AU (array/seq) / PULEO
- Consider rank transform?

---

## Amber's Work (DeSurv) — Reference Points

- Gene ID correlation with known signals
- Gene ontology used for factor characterization

---

## Repository / Manuscript Organization

- Organize scripts into the manuscript repo
- Running Quarto (.qmd) draft
- **Task:** set up a dedicated paper repo
  - Use SSBMF repo: https://github.com/andrew-walther/SSBMF-paper

**Action item:** Need to set up a follow-up meeting.

---

## Workflow / Tooling Notes

### Git Commit Review Automation
- On each commit, trigger a Codex review
  - Pass a fixed prompt for between-commit review
  - **Claude Code task:** create a global post-commit hook that:
    1. Triggers an external Codex review of the diff between the current and prior commit
    2. Launches the review automatically after each commit

### Literature Search Efficiency
- Current approach (web search each time) has accuracy/maintenance concerns
- **More efficient alternative:** RAG (Retrieval-Augmented Generation) database
  - Download PDFs as you go, or have Claude automate the downloading
  - Build a RAG database for later literature searching

---

## Next Up
- [ ] Gene ontology for genes in features
- [ ] Paper repo — ongoing updates
- [ ] Figure out HPC cluster workflow
