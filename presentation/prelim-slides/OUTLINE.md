# Talk outline: SSBMF (prelim Part III) — placeholder

Not started. This is a starting point for the Part III deck. See `README.md` for the
format and the prelim adaptation rules.

## Starting slot list

The slot list below is copied from the prelim deck outline
(`bios-dissertation/slides/prelim/OUTLINE.md`, "Part III · SSBMF"). Idea codes such as
[P3-M1] point to that file's ideas bank (Part B). These are slot ideas only.

Before using any slot, check its claim against the current project state: the method is
being rebuilt as the `ssbmf` package, so results and recommended settings may change.
Sources: `PROJECT_STATUS.qmd`, `DECISIONS.md` and `docs/plans/`.

| Section | Slot | Claim (provisional) | Ideas |
|---|---|---|---|
| Background | 1 | Pancreatic cancer is lethal, and its tumors carry prognostic expression programs | P3-B1 |
| Background | 2 | The usual two-step pipeline finds programs first and tests survival afterward | P3-B2 |
| Background | 3 | No method yet combines survival supervision, uncertainty and data-driven rank | P3-B3, P3-B4 |
| Methods | 4 | SSBMF fits a matrix factorization and a Cox model together (notation slide) | P3-M1 |
| Methods | 5 | A new patient's risk score needs only their molecular data | P3-M2 |
| Methods | 6 | Start with too many programs, then shrink | P3-M3 |
| Simulation Results | 7 | Under-specifying the starting rank merges true programs | P3-R1 |
| Application | 8 | Two survival-associated programs match the basal-like and classical subtypes | P3-R3 |
| Application | 9 | External concordance slightly exceeds a tuned two-step pipeline | P3-R4, P3-R5 |
| Application | 10 | Both programs keep their direction in held-out cohorts | P3-R6 |
| Methods | 11 | The multimodal extension adds methylation with its own priors | P3-M5, P3-M6 |
| Simulation Results | 12 | The multimodal model selects the correct number of programs in simulation | P3-R8 |
| Application | 13 | On real data, methylation adds little so far | P3-R9 |
| Discussion | 14 | Take-home | P3-T1 |
| Discussion | 15 | Limitations and remaining work | P3-T2, P3-T3 |
| Backup | B1–B14 | Backup ideas | P3-K1 to P3-K14 |

The prelim outline counts this as slightly over budget, so expect merges (5 into 4, 11
into 12). The sections follow Parts I and II. The multimodal slots (11–13) may fit
better as their own short block before the Discussion.

## Open questions for Andrew

- How much of the multimodal (methylation) work belongs in the main talk, given where
  the rebuild stands in December?
- Which results will be final enough to present by 2026-12-07?
- Target length: about 20 minutes, or shorter while the project is still in progress?
- Scope wording for the take-home, so that no claim implies survival-guided discovery or
  a methylation benefit (prelim outline open question 15).
