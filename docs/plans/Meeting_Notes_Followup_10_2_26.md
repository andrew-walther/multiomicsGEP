# Meeting Notes and Follow-up Plan — SSBMF Project (October 2, 2026)

Sources: Andrew's written notes and the Zoom summary from the 10/2 meeting with Yusha and Naim, plus
two post-meeting clarification sessions. The results discussed are in
`docs/progress_book/chapters/2026-10-02.qmd`. This plan supersedes the open-work sections of
`Working_Plan_10_1_26.md`; that file is still the record of what was done on 10/1–10/2.

**Timeline:** prelim chapter due 11/16, prelim exam 12/7, target final defense April 2027,
graduation May 2027. Naim asked for the paper writing to speed up: draft as complete a paper as
possible before the prelim.

---

## Deliverables

| # | Deliverable | Due | Workstream |
|---|---|---|---|
| D1 | Multimodal derivation PDF sent to Yusha and Naim | week of 10/5 | A0 |
| D2 | Unified write-up (methods, derivations, algorithms, simulation, real data, biology) sent to Yusha and Naim | first full draft by 11/2 | A |
| D3 | **Survival supervision of F**: a method in which survival shapes the program loadings, shown to recover a planted low-variance prognostic program in simulation | method fixed by 11/9 | C |
| D4 | Feature-count sensitivity: C and agreement between factors as the number of genes and CpGs grows | first results by 10/25 | B |
| D5 | Prelim chapter updated from D2–D4 | 11/16 | A |
| D6 | Journal choice and move to that journal's template | after D2 | A4 |

D3 is required to finalize the method. A method in which survival never informs F (α_F = 0) is
not the final method; see `DECISIONS.md` 2026-10-02.

---

## What the meeting concluded

1. **Write it up first.** A unified document covering methods (derivations and algorithms for the
   single- and multiple-modality models), results (simulation and real data) and biology, sent to
   Yusha and Naim. The target journal is chosen from that document, not before.
2. **Filter features less.** The Bayesian model has separate sparse priors per modality and a
   per-feature noise variance. It should be able to separate signal from noise without the
   stringent screening inherited from the NMF/DeSurv pipeline (3,000 genes; top 10,000 CpGs by
   variance, about 3% of the 305,352 available).
3. **Do a sensitivity analysis on the number of features.** Fix K, vary the number of features one
   modality at a time, and track the C-index and agreement between factors.
4. **Survival has to inform the factorization.** The advisors' view: the plain Bayesian formulation
   should let the priors balance the matrix-factorization and survival terms. If needed, add an
   estimated mixing parameter.
5. **Two acceptable framings for the multimodal model.** Either the joint model improves C over the
   single-modality model, or it finds biological factors that single-modality analysis does not,
   with no gain in C. First confirm the implementation works as intended, then decide how to frame
   it.

## Corrections and clarifications from the post-meeting sessions

- **L is still nonnegative** (point-exponential prior). The *loadings* F are now signed
  (point-Laplace or Normal prior, fit by `ebnm`). The intercept made signed F necessary: loadings
  measure deviation from each feature's baseline level. The meeting notes say "L is now not
  non-negative", which is wrong.
  - `config/globals.yml` still sets `prior_F` to `point_exponential` for both modalities. The
    real-data runner sets `point_laplace` when the pruned framework is used
    (`results/multimodal_real/run_multimodal_real_fit.R:65`). Make the defaults consistent before
    the sensitivity runs.
- **Why nothing is pruned on real data.** The pruning code behaves as designed. Every program's
  molecular-fit gain (−39,000 to −302,000 in ELBO) is far larger than its complexity cost (+9,000 to
  +29,000), while survival contributes at most −2.9 (chapter §D). On real tumors the ELBO is
  effectively asking whether a program explains molecular variance, and many sources of variation
  (purity, stroma, immune content) do. Most real-data fits also stopped at the 500-iteration cap,
  so the removal comparisons use unconverged fits. That is unlikely to change the conclusion given
  the margins, but it has not been checked.
- **The two methylation programs are different kinds.** The methylation-only program carries 59% of
  the methylation signal. The immune program is a shared program defined by immune genes that also
  carries 21% of the methylation signal. Neither is prognostic.
- **Top-variance CpG selection probably shapes what is found.** The highest-variance CpGs mostly
  track cell composition, so the current filter overrepresents the immune and purity signal.
  Composition is still the largest source of variation in methylation data and will produce
  programs under any filter.
- **The sparse priors on F select features within a program. They do not balance modalities or
  survival in the total objective.** Every feature contributes a reconstruction term to the ELBO,
  whether or not its loading is shrunk to zero. This matters for ELBO-based pruning. For the F
  update itself, the balance works per gene (see the mechanism section below).
- **The α_F history.** The meeting notes say "α_F > 0 led to instability". The 2026-08-20 entries in
  `DECISIONS.md` record that the instability does not reproduce on real data under the current
  preprocessing. A 16-configuration grid (α_F ∈ {0, 0.1, 0.3, 0.5} × freeze length) gave no gain in
  external C over α_F = 0. The grid never tested the case α_F exists for: recovering a low-variance
  prognostic program, which can only be measured in simulation against a known truth.

---

## How survival enters F, and why it currently has almost no effect

This section explains the mechanism so that the fixes in Workstream C can be judged.

### The update

The single-modality model has η = (YF)β̃, with program scores Z = YF. For gene j and program k, the
coordinate-ascent update for q(f_jk) combines two sources (`code/update_F_surv_YFB.R:27-31,136`):

```
A_F[j] = (1 − α_F) · τ_j Σ_i E[l_ik²]          +  α_F · E[β̃_k²] · Σ_i w_i y_ij²
B_F[j] = (1 − α_F) · τ_j (R_kᵀ E[L_k])[j]      +  α_F · E[β̃_k]  · Σ_i w_i z₋ₖ,i y_ij
          └──────── genomics ────────┘            └──────────── survival ────────────┘
posterior for f_jk: ebnm(x = B_F / A_F, s = 1 / √A_F)
```

Here w_i and z₋ₖ are the Cox model's quadratic-approximation weights and working response with
program k's contribution removed.

- **α_F = 0 drops survival from the F update entirely.** The plain Bayesian posterior adds both
  sources with weight 1. That gives the same pseudo-observation x = B_F / A_F as α_F = 0.5, with
  twice its precision. The advisors' "simple Bayesian formulation" is therefore the equal-weight
  update, not α_F = 0.
- The multimodal model already uses the equal-weight update with no α
  (`code/multimodal_yfb_updates.R:471-472`). It is stable there, and it still misses a planted
  low-variance prognostic program in simulation (chapter §B). This is direct evidence that equal
  weights alone are not enough.

### Three separate problems

1. **Gating by β.** The survival terms are multiplied by E[β̃_k] and E[β̃_k²]. A program whose
   coefficient is near 0 receives no survival information in its loadings. That is correct for a
   program that really is not prognostic. It also means a low-variance prognostic direction that no
   current program points at can never be found: no program carries a nonzero β toward it, so
   survival never pulls any program's loadings in that direction. Empirical-Bayes shrinkage of β
   makes this stronger.
2. **Scale imbalance.** At initialization the genomics precision was about 10⁴ times the survival
   precision (`DECISIONS.md` 2026-04-30; measured under the April preprocessing, so the ratio
   needs re-measuring). Under the model as written this is the correct Bayesian weighting: each
   gene's loading is informed by n expression values with residual precision τ_j, but by survival
   only through a 75-event Cox model on a sum over all genes. Equal weights are therefore not
   "balanced" in the sense the advisors intend. This is the argument for a weight or a structural
   change, and the reason it is justified rather than ad hoc: the low-rank Gaussian model for Y is
   misspecified, and weighting likelihood terms is the standard generalized-Bayes response to
   misspecification.
3. **Scale feedback (the April instability).** LFᵀ is unchanged if F_k is multiplied by c and L_k
   divided by c. The genomics term does not fix the scale of F, but the survival term depends on it
   through Z = YF. When survival pushed F_k up, L_k shrank in the next L update, which lowered the
   genomics precision for F_k, which raised survival's relative weight, and the loop ran away
   (`DECISIONS.md` 2026-04-30, steps 1–5). The multimodal model removes this freedom by rescaling
   each program every iteration so that sd(YF_k) = 1, adjusting L, β and the priors to match
   (`code/fit_multimodal_yfb.R:162-213`). The single-modality model has no such step.

So the update scheme does need revising. Problem 3 needs a fix to the scale (canonicalization).
Problem 1 needs a way to start a program in a survival-relevant direction. Problem 2 needs a weight
or a different prior on β. Each can be tested separately.

---

## Workstream A — Unified write-up (start of the manuscript)  `[Priority 1, runs throughout]`

Goal: one document that a reader can use to understand, reproduce and judge the method, and that
becomes the manuscript once a journal is chosen. Write the parts that will not change first.
Sections that depend on Workstream C get placeholders, not stale numbers.

**A0 (week of 10/5).** Send Yusha and Naim the multimodal derivation,
`derivations/multimodal_YFB/multimodal_YFB_derivation.pdf` (nine pages), with
`multimodal_YFB_technical_notes.pdf` as the long version. Check first that it reflects the 10/1
changes (intercept, signed F, `ebnm`, pruning); re-render if not.

**A1. Outline and location.** `paper/multiomicsGEP_manuscript.qmd` is a 26-line Quarto stub from
February. Replace it with the outline below, or start a new file and delete the stub (decision for
Andrew). Use the existing `paper/ssbmf-refs.bib`.

1. **Introduction:** PDAC subtypes, supervised factorization, DeSurv and how this model differs.
2. **Methods**
   - Single-modality model: η = (YF)β, priors, CAVI updates, ELBO, ARD-based K selection.
     Sources: `docs/notes/YFB_derivation_05_08_26.qmd`,
     `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.pdf`,
     `docs/update_*.qmd`, `code/SupervisedMF_Context.md`.
   - **Survival supervision of F:** placeholder until Workstream C is done. It will contain the
     update, the scale canonicalization and whatever C4 decides about weights or the β prior.
   - Multiple cohorts: per-platform z-standardization; the optional `cohort_id`, `strata_id` and
     `beta_cohort_id` terms (see A3).
   - Multimodal model: shared L, modality-specific F_m, intercept, `ebnm` priors, pruning.
     Source: `derivations/multimodal_YFB/`.
   - One algorithm box per model, with pseudocode.
3. **Simulation:** single-modality benchmark, multi-cohort simulation, multimodal pruning sweep,
   recovery of a low-variance prognostic program (from C).
4. **Real data:** single-modality external validation (5 cohorts, KM curves, two-step baselines);
   multimodal TCGA → ICGC; feature-count sensitivity (from B).
5. **Biology:** Program 7 (basal-like, adverse) and Program 3 (classical, protective); pathway
   enrichment; the multimodal programs (chapter §E). Rerun the characterization if C changes the
   final programs.
6. **Discussion:** limitations and open problems.

**A2. Fill each section from existing material** (progress-book chapters, `PROJECT_STATUS.qmd`,
the 7/15 executive summary, the prelim chapter). Mark every number with its source file.

**A3. Multiple cohorts: find out what is still open.** The notes say "need to address the multiple
cohorts method (standardized → looked more at pooled)" and "C-index collapse?". The record shows:
per-platform z-standardization is required for mixed RNA-seq and proteomics training; 10 of 12
other preprocessing choices collapse to β = 0; the pooled TCGA + CPTAC fit beats single-cohort
fits; the cohort-specific options are performance-neutral (`ROADMAP.md`, `DECISIONS.md` 2026-09-04
addendum). Ask Yusha and Naim whether the open issue is presentation or a method change before
starting new analysis.

**A4. Journal choice** (after A2): statistics journal versus subject-matter journal, depending on
whether the biology turns out to be the main finding.

Verify: D1 and D2 are sent; every number traces to a file in the repo.

## Workstream C — Survival supervision of F  `[Priority 1 for method work]`

Goal (D3): survival shapes the program loadings strongly enough to recover a low-variance prognostic
program, without instability and without losing external C.

**Test bed for every step:** the planted low-variance prognostic scenario. Use the multimodal
simulation's "low-variance prognostic" scenario run with one modality, and the single-modality
benchmark simulator with an added prognostic program that explains under 1% of variance. Use at
least 10 seeds per setting.

**Pre-specified success criteria** (thresholds as in the prelim chapter, `DECISIONS.md`
2026-10-02):
- the planted program is recovered, with absolute loading correlation ≥ 0.7 and the correct sign
  of β, in a large majority of seeds (Andrew to set the rate before C4 runs);
- false prognostic programs do not increase relative to α_F = 0;
- every fit converges, with no runaway in ‖F_k‖ or ‖L_k‖;
- on real data, the lower 95% limit of ΔC against the α_F = 0 model on the 5 external cohorts stays
  above −0.02.

**Selection metric for any weight or tuning choice:** cross-validated C on the training data, as
agreed. Report the cross-validated survival partial log-likelihood (`code/compute_cv_loglik.R`)
alongside it: with 75 events, C is coarse and may not separate nearby settings, and the
log-likelihood is a smoother check.

**C0. Branch. Done 2026-10-12.** `codex/multimodal-yfb` was merged to `main` (`2f5dd64`, 531/531
tests). Workstream C runs on `feature/supervised-F`, created from `main`. The write-up session shares
the same checkout, so check `git branch --show-current` before every commit.

**Findings from scoping C1** (from reading the code on 2026-10-12; C1 should confirm them):
1. **The monitored objective and the F update use different weightings.** The ELBO that
   `fit_cox_on_yf()` reports, and checks for convergence, is
   (1−α)·genomics/divisor + α·survival/divisor, with α = 0.5 and a normalizer on each term
   (`code/fit_cox_on_yf.R:780`). The F update instead mixes the raw, un-normalized genomics and
   survival precisions with α_F (`code/fit_cox_on_yf.R:727`, `code/update_F_surv_YFB.R:136`). The
   updates therefore do not climb the objective being monitored. Whatever weighting C4 chooses must
   enter both.
2. **F is the only route for survival into the factorization.** Under η = (YF)β the L update is
   genomics-only (`code/fit_cox_on_yf.R:701`). With α_F = 0 survival is absent from the
   factorization entirely.

**C1. Diagnostics, no method change (2–3 days; one step at a time, each taught back before the
next).**
1. **Trace.** Add an opt-in `trace_supervision = FALSE` argument to `fit_cox_on_yf()`. When on,
   `history` records per program k and iteration:
   - E[β̃_k] and E[β̃_k²];
   - survival's share of the loading precision, α_F·A_surv / A_F (median and 90th percentile over
     genes), and the same share for B_F;
   - sd(YF_k), ‖F_k‖ and ‖L_k‖.

   Verify:
   - with the trace off, output is bit-for-bit identical (new test);
   - with it on, the logged share matches a hand computation on a toy fit;
   - the full suite still passes (531).
2. **Test bed.**
   - Data: `simulate_multimodal_yfb_data(scenario = "low_variance_prognostic_factor")`, keeping
     the expression block only, moderate noise, 10 seeds.
   - Fit: `fit_cox_on_yf()` at K = 7 with the recommended settings, apart from the arms in step 3.
   - Verify: the planted program explains under 1% of variance.
3. **Arms.**
   - α_F = 0;
   - α_F = 0.5;
   - the exact equal-weight Bayesian update: α_F = 0.5 with both precision terms doubled, using a
     new flag that is off by default.

   For each arm record:
   - recovery of the planted program (absolute loading correlation ≥ 0.7 and the correct sign of β);
   - false prognostic programs;
   - held-out C;
   - convergence.

   The α_F = 0 arm should reproduce non-recovery. The runner is
   `results/benchmark_sim/run_supervision_diagnostics.R`. Time one fit before running the rest.
4. **Diagnostic note** with one figure: survival share and E[β̃] across iterations, by arm, for the
   estimated program closest to the planted one. It answers four questions:
   - Is the planted program's direction pulled at all?
   - What is the survival share now (the April ratio of about 10⁻⁴ is stale)?
   - Does the April scale runaway reproduce in simulation?
   - How far apart are the monitored objective and the F update (finding 1)?

   The note goes into the progress book, and any decision into `DECISIONS.md`. It also feeds the
   write-up's methods discussion.

**C2. Fix the scale (problem 3).**
- Port the multimodal canonicalization (rescale each program to sd(YF_k) = 1 with matching changes
  to L, β and the priors) into `fit_cox_on_yf()`.
- Tests:
  - LFᵀ, η and the ELBO are unchanged by one canonicalization step;
  - α_F = 0 fits are unchanged up to scale;
  - the existing suite still passes.
- Rerun C1. Expected: no runaway at any α_F, which makes C3–C4 interpretable.
- Then drop α_F as a separate knob: the update becomes equal-weight Bayes, plus at most one survival
  weight w_S from C4.

**C3. Start a survival-relevant program (problem 1).** Two options; implement the first, the
second only if needed.
- **Supervised initialization:** initialize one program's loadings from the Cox score direction
  Yᵀ(martingale residuals) of a null or α_F = 0 fit. Its β then starts away from 0.
- **Greedy supervised program:** after a fit converges, add a program initialized from the Cox
  score direction of the current residual risk, backfit, and keep it only if the selection metric
  improves. This follows the greedy-then-backfit pattern of `flashier`.

**C4. Weight survival or loosen the β prior (problem 2).** Compare on the test bed:
- equal weights (the plain Bayesian model) as the baseline;
- **prior on β:** empirical Bayes (current) versus a fixed wider variance versus a hierarchical or
  heavier-tailed prior. This is the "let the priors balance it" option the advisors proposed:
  E[β̃_k²] sets how much survival can move F_k;
- **one survival weight w_S** on the Cox term (a tempered likelihood), chosen by cross-validated C
  from a small grid. Note: w_S cannot be estimated by maximizing the ELBO, because the weighted
  ELBO changes with w_S and has no fixed target. If a data-driven estimate is wanted later,
  consider generalized-Bayes learning-rate methods (SafeBayes, Grünwald; Holmes and Walker 2017).
  These have not been checked in detail.

Pick the simplest option that meets the success criteria. A prior-based solution is preferred over
a weight if both meet them: it needs no extra tuning and fits the Bayesian framing.

**C5. Revise the update scheme, only if C2–C4 fall short.**
- Joint update of (F_k, β_k) in one block instead of separate coordinate steps. The mean-field
  split q(F)q(β) ignores their posterior correlation, which is part of why gating is so strong.
- Damped F updates.
- Several β updates per F update, or another update order.

**C6. Carry the result into the multimodal model.**
- Use the same β prior or weight, plus modality weights w_m ∝ 1/p_m if Workstream B shows that
  methylation dominates at larger p.
- Pruning: either the weighted ELBO, or prune and select K by cross-validated C (`code/select_K.R`).
  Compare the two on the multimodal simulation and on TCGA → ICGC.

**C7. Real-data confirmation.**
- Refit the recommended single-modality configuration with the chosen method (K = 7, 5 external
  cohorts).
- Report ΔC against the α_F = 0 model and the two-step baselines.
- Compare the programs against Programs 3 and 7. Rerun pathway enrichment if the programs change.

Record C2, C4 and C6 outcomes in `DECISIONS.md`. Add C1 and C4 results to the progress book.

## Workstream B — Feature-count sensitivity  `[Priority 2; runs in the background during A and C]`

Goal (D4): how C, the number of retained programs and the factors themselves change as more
features are included.

Workstream B starts with the current method. The B1 and B4 runs that matter for the paper are
repeated once C settles the method. The feature-screening pipeline built here is reused unchanged.

**B0. Preparation.**
- Make the `prior_F` defaults consistent (see Corrections).
- Add a feature-count argument for each screening rule to `screen_multiomics_features()`
  (`code/load_multiomics_data.R:229`). Screening stays TCGA-only.
- Decide the iteration cap and convergence check for large p.

**B1. Expression only** (the multimodal fitter with one modality, supported since `387062e`).
- Fixed K = 7, pruning off.
- Genes: 2,000 / 3,000 (current) / 5,000 / 10,000 / 15,000 / all 17,428 shared. Rank by the DeSurv
  rule so that each set contains the smaller ones. Use a lenient mean-expression floor for the
  larger sets.
- Outputs for each fit:
  - ICGC primary C (n = 50) and all-donor C (n = 67);
  - TCGA cross-validated C;
  - ELBO, iterations, converged or not, runtime;
  - which programs are prognostic;
  - agreement with the 3,000-gene fit: match factors by correlation of L, then correlate F over the
    shared genes.
- Then repeat at 2–3 feature counts with K_init over a range and pruning on.

**B2. Methylation only.**
- For all 305,352 CpGs on TCGA, compute mean beta, SD, and the proportion of samples outside
  [0.1, 0.9].
- Compare two filters by size and overlap: the mean rule from the meeting (drop mean > 0.9 or
  < 0.1), and a proportion rule (drop a CpG if ≥ 95% of samples are below 0.1 or above 0.9). The
  mean rule would also drop tumor-specific hypermethylated CpGs that are near 0 in most tumors.
- Plot a histogram of per-CpG SD for the retained CpGs (Yusha's request), on the beta and arcsine
  scales. Choose the SD cutoffs from it.
- Check whether the upstream probe files already removed cross-reactive and SNP-overlapping probes.
  If not, remove them.
- Fit methylation only at K = 7 for about 10,000 / 25,000 / 50,000 / 100,000 / all retained CpGs.
  Same outputs as B1.

**B3. Joint fits** at matched settings from B1 and B2, compared against expression-only at the same
gene count.

**B4. Single-modality SSBMF (`fit_cox_on_yf.R`, 7-cohort pipeline).**
- Refit with 2,000 / 3,000 / 5,000 / 10,000 genes.
- Compare against Amber's DeSurv gene set: overlap, correlation of F over the shared genes, and
  external C.

**B5. Simulation at large p.**
- Rerun the pruning simulation with many added noise features (e.g. p = 15,000).
- If pruning stops working there too, the real-data result is mainly caused by p.

**Compute.**
- At 305,352 CpGs × 144 patients, one matrix is about 350 MB, and the F update is about 30× slower
  per iteration than at 10,000 CpGs.
- Run the large B2/B3 fits on Longleaf: one fit per SLURM array task, about 4–8 GB of memory each.
- For local runs use socket clusters; forked workers segfault on macOS.

Verify: one table and one figure per modality (C and agreement versus number of features); the B5
result stated before interpreting B1–B3.

---

## Order of work

*Superseded on 2026-10-13 by "Sequencing and handoffs" below; kept as the original schedule.*

| When | Write-up (A) | Method (C) | Sensitivity (B) |
|---|---|---|---|
| 10/5–10/11 | A0 send derivation; A1 outline | C0 merge and branch (done 10/12); C1 diagnostics | B0 preparation; start B1, B5 |
| 10/12–10/18 | A2: single-modality methods, simulation, real data | C2 canonicalization and tests; C3 supervised initialization | B2 CpG filters and histogram |
| 10/19–10/25 | A2: multimodal methods, biology | C4 β prior and w_S comparison | B1 and B2 results; B4 |
| 10/26–11/1 | A2: Results placeholders filled | C4 decision; C6 multimodal; C5 only if needed | B3 joint fits |
| 11/2–11/8 | **D2 sent** | C7 real-data confirmation | B1/B4 rerun with the new method |
| 11/9–11/16 | **D5 prelim chapter** from D2 | write-up of C into Methods | — |
| 11/17–12/7 | Prelim preparation; A4 journal choice | — | — |

If C runs late, D2 goes out on 11/2 with the supervision section as a described plan and the C1–C4
results so far, rather than waiting.

## Decisions needed (Andrew and advisors)

1. **C4:** prior-based balance (β prior) or an explicit survival weight w_S; rate threshold for
   "recovered in a large majority of seeds".
2. ~~**C0:** merge `codex/multimodal-yfb` to `main` now.~~ Done 2026-10-12 (`2f5dd64`).
3. **A1:** replace the February manuscript stub, or start a new file.
4. **A3:** what the open multiple-cohorts issue is.
5. **B2:** mean-only CpG filter or proportion rule.
6. **A4:** statistics journal or subject-matter journal (after the write-up).

---

## Sequencing and handoffs (decided 2026-10-13)

This section replaces the earlier Workstream C kickoff prompt and the parallel schedule in "Order of
work". The write-up comes first because Andrew shares it on Friday 10/16 (draft targeted for
Wednesday 10/14). Method work follows, and its results then revise the write-up, which is the
foundation for the manuscript.

**One session at a time, one checkout.** Claude Code and Codex take turns in
`~/GithubProjects/multiomicsGEP`; no worktrees and no second clone. Each task has its own branch off
`main`. Run `git branch --show-current` before every commit, stage only the files you changed, and
don't push without asking.

| Step | Who | Branch | Prompt |
|---|---|---|---|
| 1. Review the planning page; Codex revises it and writes the drafting prompt | Andrew + Codex | `writeup/draft` (new, off `main`) | Prompt 1 |
| 2. Draft the write-up | Claude Code | `writeup/draft` | Prompt 2 (written by Codex in step 1) |
| 3. Review the draft | Codex | `writeup/draft` | Prompt 3 |
| 4. Fix the review findings; Andrew approves; merge to `main` | Claude Code | `writeup/draft` | Prompt 4 |
| 5. Method work (M0–M8 below) | Claude Code | `feature/supervised-F`, after merging `main` | Prompt 5 |
| 6. Review the method work, once at the end | Codex | `feature/supervised-F` | Prompt 6 |
| 7. Revise the write-up from its "Revisit after feature/supervised-F" list | Claude Code | `writeup/revision` (off `main` after step 6 merges) | Prompt 7 |

### Method work scope (step 5)

Goal: a strong single-modality YFB model (expression) and multimodal YFB model (expression +
methylation), with results in simulation and on real data that show value over EBMF/NMF alone, the
two-step EBMF → Cox method and DeSurv; a K-selection procedure whose K_eff accounts for both
survival and reconstruction; and an explanation of the biology of the final programs.

| # | Work | Covers |
|---|---|---|
| M0 | Merge `main` into `feature/supervised-F`; run the test suite and report the count | Setup |
| M1 | C1 diagnostics: `trace_supervision`, planted low-variance test bed, three α_F arms, diagnostic note | Why supervision currently has little effect |
| M2 | C2 scale canonicalization in `fit_cox_on_yf()`; drop α_F as a separate knob; then define ZF (normalized projection) and settle the subscript notation in `derivations/cox_on_YF/` (10/1 plan WS4.1) | α_F; ZF clarification |
| M3 | C3 supervised initialization; C4 β prior versus survival weight, including the β prior family comparison (Normal, point-Laplace, point-normal; 10/1 plan WS4.3); C5 only if needed | α_F; prior variation |
| M4 | C6: carry the chosen method into the multimodal model | Multimodal method |
| M5 | K selection for both models: compare the weighted ELBO with pruning by cross-validated C plus reconstruction; replicated K panels; 1-SE rule (10/1 plan WS4.7) | Parsimonious K_eff |
| M6 | Robustness: jittered SVD starts × 20–100 (10/1 plan WS4.4); F prior families; Workstream B feature counts | Multiple seeds; prior variation; feature count |
| M7 | Comparator matrix for both models, simulation and real data: EBMF/NMF alone, two-step EBMF → Cox, DeSurv. First confirm whether a head-to-head DeSurv run exists; the current DeSurv scripts compare preprocessing and gene overlap | Value over existing methods |
| M8 | C7 real-data confirmation; rerun the biology for the final programs (risk scores, KM curves, pathway enrichment, CpG interpretation for the multimodal programs) | Biology |

**Working mode for step 5.** Work through M0–M8 without stopping for step-by-step teach-back, and
give Andrew a summary at the end. Stop only for:
1. the recovery rate that counts as "a large majority of seeds", set before the C4 runs;
2. the C4 choice between a prior on β and a survival weight;
3. the K-selection criterion in M5;
4. any result that fails a pre-specified success criterion (report it; don't work around it).

Each stop comes with a short recommendation.

### Prompt 1 — Codex: planning-page review and drafting prompt

```text
Read paper/writeup/planning.qmd (rendered as planning.html),
paper/writeup/implementation_plan.md, paper/writeup/source_inventory.csv, and the
"Sequencing and handoffs" section of docs/plans/Meeting_Notes_Followup_10_2_26.md.

First create the branch writeup/draft from main (git switch -c writeup/draft main).
All write-up work happens there.

I will give you feedback on the planning page. Then:
1. Revise planning.qmd and implementation_plan.md to reflect my feedback, and
   re-render planning.html.
2. Add a section "Revisit after feature/supervised-F" to implementation_plan.md.
   List every passage of the planned write-up that depends on the method work
   (M0-M8 in the sequencing section): the single-modality F update and alpha_F,
   the definition of ZF, K selection, comparator results, and the final programs
   and their biology.
3. Add a task that produces an evidence-gap table in the write-up: each claim the
   manuscript needs to make, the existing evidence and its source file, and what
   is missing. Claims to cover: value over EBMF/NMF alone, over two-step
   EBMF -> Cox, and over DeSurv; K selection that accounts for survival and
   reconstruction; the multimodal model's added value; the biology.
4. Write the drafting prompt for Claude Code to paper/writeup/drafting_prompt.md.
   It must:
   - work on writeup/draft and execute Tasks 2, 3 and 5 of implementation_plan.md;
   - target a complete draft by Wednesday 10/14 for sharing on Friday 10/16, and
     state what is cut to make that feasible;
   - leave Task 4's independent review to Codex (step 3 of the sequence);
   - change no model code, defaults or results;
   - mark every passage on the "Revisit" list so it can be updated later.
Run git branch --show-current before committing, stage only your files, and
ask before committing or pushing.
```

### Prompt 2 — Claude Code: draft the write-up

```text
On branch writeup/draft, read paper/writeup/drafting_prompt.md and follow it.
```

### Prompt 3 — Codex: review the draft

*Superseded 2026-10-14 by `paper/writeup/review_prompt.md`, which adds Andrew's feedback on the initial draft. Use: "On branch writeup/draft, read paper/writeup/review_prompt.md and follow it."*

```text
On branch writeup/draft, review the write-up draft (paper/writeup/writeup.qmd and
the appendices) against the "Review focus" and "Global constraints" in
paper/writeup/implementation_plan.md and against source_inventory.csv.
Check:
- every quantitative claim against its recorded path and selector;
- equations and algorithm boxes against the code and the evaluated configuration;
- comparator identity, cohort sizes, orientation conventions and denominators;
- overstated claims, especially about survival shaping F under alpha_F = 0 and
  about the biology;
- that every passage on the "Revisit after feature/supervised-F" list is marked.
Write findings to paper/writeup/review.md, ranked by severity, each with
location, evidence and a suggested fix. Do not edit the draft. Ask before
committing review.md.
```

### Prompt 4 — Claude Code: fix review findings

*Amended 2026-10-14: also apply every item in "Author feedback on the initial draft" in `paper/writeup/implementation_plan.md`.*

```text
On branch writeup/draft, read paper/writeup/review.md. Verify each finding
against the source before acting on it; fix the confirmed ones, and record the
disposition of every finding (fixed, not confirmed with reason, or deferred to
the method work) in review.md. Re-render HTML and PDF, then summarize the changes
and ask me before merging writeup/draft to main.
```

### Prompt 5 — Claude Code: method work

```text
Resuming multiomicsGEP method work (Workstream C and the M0-M8 scope). The
write-up has been drafted, reviewed and merged to main.

Read first:
- the project CLAUDE.md;
- docs/plans/Meeting_Notes_Followup_10_2_26.md: "How survival enters F",
  Workstream C including "Findings from scoping C1", and "Sequencing and
  handoffs" (scope M0-M8 and the working mode);
- DECISIONS.md, entry "2026-10-02 — alpha_F = 0 cannot be the final method";
- paper/writeup/implementation_plan.md, section "Revisit after
  feature/supervised-F", and the write-up's evidence-gap table.

Then:
1. git fetch; confirm main is clean and report its head. Switch to
   feature/supervised-F and merge main into it (no rebase, no force-push).
   Run Rscript tests/run_tests.R and report the pass count.
2. Work through M0-M8 in order without stopping for teach-back. Stop only at the
   four decision points listed under "Working mode for step 5", each with a
   recommendation.
3. Keep the success criteria in Workstream C fixed; don't relax them after
   seeing results. Record decisions in DECISIONS.md and results in the progress
   book.
4. At the end, give me a summary: what changed, results against each success
   criterion, and what the write-up's "Revisit" list now needs.
Run git branch --show-current before every commit, stage only your files, and
don't push without asking.
```

### Prompt 6 — Codex: review the method work

```text
Review the method work on feature/supervised-F (git diff main...feature/supervised-F).
Read the "Sequencing and handoffs" section and Workstream C of
docs/plans/Meeting_Notes_Followup_10_2_26.md first.
Check:
- that the code matches the derivations and DECISIONS.md entries it cites;
- that the success criteria were applied as pre-specified, not adjusted after
  the results;
- that tests assert correct behavior, not just output shape; run
  Rscript tests/run_tests.R and report the count;
- simulation and comparator design (seeds, planted program, baselines,
  cross-validation without leakage);
- that reported numbers trace to saved outputs.
Write findings to docs/reviews/<date>_supervised_F_review.md, ranked by
severity, each with location, evidence and a suggested fix. Do not edit the
code. Ask before committing.
```

### Prompt 7 — Claude Code: revise the write-up

```text
feature/supervised-F has been reviewed and merged to main. Create
writeup/revision from main. Update the write-up from the "Revisit after
feature/supervised-F" list in paper/writeup/implementation_plan.md and from the
method-work summary, refresh the evidence-gap table, re-render HTML and PDF,
and ask me before committing.
```
