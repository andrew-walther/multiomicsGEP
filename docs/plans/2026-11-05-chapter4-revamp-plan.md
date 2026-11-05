# Chapter 4 revamp: supervised YFB model as the centerpiece

Written 2026-11-05 from a bios-dissertation planning session (Round 3 prose review of
Chapters 1–3). Status: **plan, not started**. Follow-up prompt:
[`2026-11-05-chapter4-revamp-prompt.md`](2026-11-05-chapter4-revamp-prompt.md).

Chapter source: `paper/prelim/project3-ssbmf.qmd`. The post-commit hook syncs it into
`bios-dissertation/prelim/project-proposals/project3-ssbmf/draft/` on any branch. Edits go on
`fix/yfb-supervision`; mirroring to `main` is optional, since that branch will be merged later
(Andrew, 2026-11-05).

**Deadline.** Revamp done by about 2026-11-07, so Andrew has about a week to review and edit before
the written prelim is due on 2026-11-16. Any checklist item that can't be produced by then is
written as planned evaluation, not left as a gap. The oral exam slides (2026-12-07) may show newer
results than the chapter.

## 1. Why the chapter changes

The current chapter (synced 2026-10-30, about 4,860 words) presents the frozen 2026-10-14 fit with
α_F = 0, in which survival never enters the loading update. DECISIONS 2026-10-02 says α_F = 0
can't be the final method. As of DECISIONS 2026-11-04, the preferred single-modality model is the
survival-weighted YFB fit:

- **Fitter:** `fit_ssbmf_yfb(canonicalize = TRUE, w_S = 1, prior_LF = "point_laplace",
  prior_beta = "normal")` (`code/fit_ssbmf_yfb.R`; runner
  `results/benchmark_sim/run_yfb_supervised_real.R`, lines 59–61). The code default for
  `prior_LF` is point-exponential, so the chapter must name **point-Laplace** for L and F and a
  **Normal** prior for β. The current chapter and the slide outline say otherwise.
- **Predictor:** η = (Y E[F])β. Under canonicalization each program is rescaled so that
  sd(Y F_k) = 1.
- **Objective:** E_q log p(Y | L, F, τ) + w_S × the Cox partial log-likelihood re-expanded at
  E[η], with Var(η) from q(β) and q(F), − KL. w_S multiplies the Cox term in the F update, the β
  update and the ELBO. Under canonicalization w_S = 0 also gives β = 0, so the unsupervised
  comparator is the w_S = 0 loadings followed by a separate Cox fit ("two-step", `w0cox`).
- **Rank and data:** K = 7 fixed, no pruning. D4 gene set (2,218 genes, single log, per-cohort
  z-scores), no survival screen. Training is TCGA-PAAD + CPTAC (273 patients, 139 events). w_S is
  chosen by 5-fold training CV C.

**Evidence available now** (every number from the findings files named):

| Result | Value | Source |
|---|---|---|
| S-9 simulation, 2,000 genes | planted program recovered in 10/10 datasets in every variance cell at w_S = 1; CV picks w_S = 1 in 55/60; every w_S ≥ 40 fails to converge | `docs/findings/2026-11-04-survival-weight-p2000-findings.md` |
| S-9 calibration (corrected 11-05) | slope 1.45–3.23 at 2,000 genes (attenuated, not overconfident); w_S = 10 closer to calibrated at larger variances; w_S ≥ 100 overfits | same, and DECISIONS 2026-11-04 (correction) |
| S-11 training CV, D4 | C 0.663 (SE 0.024) vs two-step 0.652 (SE 0.026); paired +0.011 (SE 0.005); 3 survival-active programs; 2 of 5 folds stop at the 300-sweep cap | `docs/findings/2026-11-04-supervised-real-data-findings.md` |
| S-11 caveat | survival supplies 10⁻⁸–10⁻⁶ of the loading precision; w_S = 1 and w_S = 0 loadings agree at \|r\| = 1.000 | same |
| S-11 size and screen | D4 preferred; survival screening adds no clear gain; null check passes (false-active 0.1–0.5 per cell) | same |
| S-12 external, 5 cohorts | C 0.584 / 0.537 / 0.656 / 0.681 / 0.644 (mean 0.620); paired gain over two-step +0.007 (95% CI −0.002 to 0.017); calibration slopes 0.36–1.57 (Puleo 0.99), closer to 1 than two-step in every cohort | `docs/findings/2026-11-05-supervised-external-scoring-findings.md` |

**What the chapter must not claim.**
- That survival shapes the real-data programs. It doesn't at w_S = 1. The gain comes from the joint
  β update with its shrinkage prior.
- That the external gain is established. The pooled interval includes zero.
- Any comparison with the 2026-10-14 benchmark (mean 0.627) as paired. That benchmark used an older
  pipeline.

Making survival move the loadings (a data-driven supervision strength) is the main open method
question. It belongs in the chapter as proposed methodology.

**Coordination.** The S-11 session (configuration review, external scoring) is still running.
Before producing anything, the follow-up reads the latest `DECISIONS.md`, `docs/findings/` and
`git log`. It doesn't rerun, retune or rescore what that session owns, and it doesn't touch its
output folders (`yfb_supervised_real_2026-11-04/`, `yfb_supervised_external_2026-11-05/`). New
outputs go to new dated folders.

## 2. Target outline

The outline matches Chapters 2 and 3: Introduction, Methods, Results (simulation, then the PDAC
application), Discussion. Proposed work sits inside the section it belongs to, not only in Future
Directions. Headings are descriptive, as in Chapter 3; the current Results headings are claims.

There are no bold run-in paragraph openers. The current chapter has 8, at lines 42, 74, 88, 129,
131, 141, 143 and 145. Use a `###` subheading where the content needs one, otherwise a plain topic
sentence. The bold "Summary." paragraph at line 42 goes: Chapters 2 and 3 open directly with the
Introduction, and its content moves into the Introduction's last paragraph.

Folder rule (`paper/prelim/AGENTS.md`): continuing work is described as additional methodology
and expected evidence, with no fixed-date timelines. Proposed methods and expected results are never
presented as findings.

1. **Introduction**
   - The problem: prognostic expression programs in PDAC, and the two-step pipeline (unsupervised
     factorization, then Cox on plug-in scores).
   - The gap: survival doesn't inform discovery, and score uncertainty is ignored.
   - What SSBMF does.
   - Aims and outline, stated as completed aims and proposed aims. *Proposed in this section:* the
     full project's aims (data-driven supervision strength, M-modality model, rank selection).
2. **Methods**
   1. Data and notation; the YFB model (Y = LFᵀ + E; η = (YF)β; Cox partial likelihood).
   2. Priors and empirical Bayes: point-Laplace on L and F, Normal on β, ebnm fits.
   3. Survival weight w_S and canonicalization: definition, and why canonicalization makes w_S
      comparable across programs. Also that w_S = 0 gives β = 0 and how the two-step comparator
      is built.
   4. Variational inference and the ELBO: sweeps, Cox re-expansion at E[η] with Var(η), stopping
      rule (tolerance 10⁻⁵, at most 300 sweeps).
   5. Choosing w_S by training CV; prediction in new cohorts (per-cohort z-scores, frozen
      orientation).
   6. Simulation design (S-9): planted low-variance prognostic program, 50 and 2,000 genes,
      variance cells, w_S grid, recovery, false activity, calibration, CV choice. Table from
      `yfb_survival_weight_2026-11-04/grid/design.csv`.
   7. PDAC data and preprocessing: training and external cohorts, single log, D4 gene set,
      evaluation metrics (C, paired bootstrap, calibration slope).
   8. *Proposed:* **Data-driven supervision strength.** Survival barely moves the real loadings at
      w_S = 1. Candidates come from DECISIONS 2026-11-04: a calibrated likelihood power
      (SafeBayes or Holmes–Walker tempering), a downweighted genomics term (S-10), or normalized
      loss terms as in DeSurv. This subsection also states the acceptance criteria.
   9. *Proposed:* **M-modality model with modality-specific priors on F.** Y_m = L F_mᵀ + E_m for
      m = 1..M, with a shared L, a program- and modality-specific point-Laplace prior g_{F_mk}
      fitted by empirical Bayes, and η from the stacked projections. M = 1 is the
      single-modality model. Expression and methylation are the planned case. Source:
      `derivations/multimodal_YFB/multimodal_YFB_derivation.qmd` (+ technical notes). It's
      written as derived but not yet implemented in this form.
   10. *Proposed:* **Rank selection.** K is fixed at 7 now. Evidence-based pruning (flashier-style)
       under the supervised objective, and the replicated K-selection criterion.
   11. *Proposed:* **Comparators still to run:** DeSurv head-to-head (`results/benchmark_sim/outputs/desurv_h2h/BLOCKED.md`;
       deferred 10-15), and the LB (η = Lβ) vs YFB comparison planned in DECISIONS 2026-11-04.
3. **Results**
   1. Simulation: recovery, false activity, calibration and convergence across w_S, and the CV
      choice of w_S (S-9). *Planned evaluation:* the same grid under the data-driven supervision
      variants.
   2. Training CV and feature-set size (S-11): CV C against the two-step fit, the size curve,
      screening.
   3. External validation (S-12): five-cohort C, paired difference, calibration slopes; KM by
      risk tertile.
   4. Programs and biology: the 7 programs, the 3 survival-active ones, top genes, enrichment,
      subtype and stroma concordance, and correspondence with the α_F = 0 programs.
      *Planned evaluation:* whether stronger supervision recovers programs the w_S = 1 fit misses.
   5. *Planned evaluation:* M-modality and rank selection. Expected evidence and decision rules
      only, with no October multimodal numbers.
4. **Discussion**
   - What the supervised fit does: better-calibrated β, a small gain in discrimination. What it
     doesn't do: move the loadings.
   - Limitations: K fixed; sweep cap; external gain interval includes zero; tertile risk groups
     are within-cohort.
   - Continuing work, tied back to Methods 2.8–2.11, plus transferable risk classifications
     (ROADMAP 2026-10-16) and the R package.

**Material leaving the chapter** goes to `paper/prelim/notes/trimmed-material-2026-11.md`, dated,
with the reason:
- the α_F = 0 framing and its results (tab:external 0.627 vs 0.601, fig:km-risk, fig:km-programs,
  fig:programs, fig:delta-c, fig:percohort), unless a supervised counterpart replaces them;
- the October multimodal parsimony results (fig:mm-pruning, tab:mm-external, the "Multimodal
  model" methods text). Multimodal is proposed work (Andrew, 2026-11-05).

The 2026-10-14 benchmark may still be cited once as the historical reference. It is labelled as
such, per DECISIONS 2026-11-04 "Presentation".

## 3. Results checklist (cross-reference with Andrew's in-progress work)

Status: **exists**, an output on disk is ready to use; **partial**, data exists but the
figure/table doesn't; **missing**. "Andrew started?" is for Andrew to fill in. Paths are relative
to the multiomicsGEP root and were checked on disk 2026-11-05.

### Tables

| # | Item | Replaces | Status | Data source | Script | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|---|---|
| T1 | S-9 simulation design | tab:sim-design | partial | `results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04/grid/design.csv` | typed from design.csv | write LaTeX table | | [ ] |
| T2 | Training CV: C, SE, paired difference, folds converged, survival-active count | (new) | partial | `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/summary.csv`, `fits.csv` | `results/benchmark_sim/summarize_yfb_supervised_real.R` | LaTeX table | | [ ] |
| T3 | External C in 5 cohorts, two-step C, paired difference, calibration slope | tab:external | partial | `results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05/external_cindex.csv`, `external_paired_diff.csv` | `results/benchmark_sim/score_yfb_supervised_external.R` (S-11 session; don't rerun) | LaTeX table; historical column labelled not paired | | [ ] |
| T4 | Programs: β, SD, z, survival-active, top genes at each pole, short label | (in-text list) | partial | `results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04/programs_D4.csv` | — | table + biological labels (needs B1–B3) | | [ ] |
| T5 | Enrichment summary for the survival-active programs | (new) | missing | all-patient fit `fits/size_D4_N2218_w1_f0.rds` (`state$EF`, `genes`) | adapt `code/pathway_enrichment.R` (written for the old D4 fit) | run, then table | | [ ] |

### Figures

| # | Item | Replaces | Status | Data source | Script | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|---|---|
| F1 | Simulation: recovery, false activity, survival share, calibration, convergence and C against w_S | fig:delta-c | exists (progress book) | `yfb_survival_weight_2026-11-04/grid/metrics.csv`, `cal_slope_exact.csv` | `docs/progress_book/figures/2026-11-06/survival-weight-curves.png` (written 11-05 04:00, after the slope correction; its plotting script wasn't found by name, so locate it) | confirm it uses `cal_slope_exact.csv`; restyle for the chapter | | [ ] |
| F2 | CV choice of w_S across datasets | (new) | partial | `grid/cv_chosen.csv`, `cv_chosen_cells.csv` | — | small figure or fold into T1/text | | [ ] |
| F3 | Real-data CV C against feature-set size | (new) | exists | `yfb_supervised_real_2026-11-04/size_elbow.png` (= progress book `real-data-size-curve.png`) | `summarize_yfb_supervised_real.R` | restyle for the chapter | | [ ] |
| F4 | External C by cohort with intervals, w_S = 1 vs two-step | tab:external visual | partial | `external_cindex.csv` | progress book `external-concordance.png` is the **October** fit | new plot from S-12 CSVs | | [ ] |
| F5 | KM by risk tertile in external cohorts | fig:km-risk | missing | needs per-patient external risk scores; S-12 saved only aggregates | extend `score_yfb_supervised_external.R` to save risks (same fixed model, no refit) | **ask Andrew/S-11 session first**: is saving per-patient scores within the "scored once" rule? | | [ ] |
| F6 | KM by program-score tertile for the survival-active programs | fig:km-programs | missing | same as F5 | same | same question | | [ ] |
| F7 | Program heatmap (loadings of top genes, 7 programs) | fig:programs | partial | `fits/size_D4_N2218_w1_f0.rds` | adapt `paper/prelim/make_figures.R` (uses `load_d4_weights()` for the old fit) | new figure | | [ ] |
| F8 | Correspondence of w_S = 1 programs with the α_F = 0 benchmark programs (and with w_S = 0) | (new) | partial | f0 fits for w1 and w0; `results/benchmark_snapshot_2026-10-14/single_modality/` | progress book `program-correspondence.png` is the October fit | new figure or table; links the new programs to the old Programs 3/7 | | [ ] |
| F9 | Calibration: external slopes, w_S = 1 vs two-step | (new) | partial | `external_cindex.csv` (`cal_slope`, `cal_se`) | — | small dot plot or a T3 column | | [ ] |
| F10 | Per-cohort β (TCGA vs CPTAC) for the survival-active programs | fig:percohort | missing | needs per-training-cohort fits | — | optional; drop if not produced by 11-07 | | [ ] |
| F11 | Model diagram (Y, L, F, β, w_S) | (new) | exists | `docs/progress_book/figures/two-step-vs-joint.svg`, `supervision-states.svg` | `build_diagrams.py` | optional; convert to PDF/PNG if used | | [ ] |

### Biology

| # | Item | Status | Source / method | Remaining | Andrew started? | Done |
|---|---|---|---|---|---|---|
| B1 | Top genes and pole labels for the 3 survival-active programs (hypoxia/glycolysis vs normal pancreas; ribosomal/proliferation vs classical epithelial; third adverse program) | partial | `programs_D4.csv`; the findings file's gene lists | confirm labels against B2–B3 | | [ ] |
| B2 | Pathway enrichment (MSigDB Hallmark, Reactome/KEGG via fgsea or clusterProfiler) | missing | adapt `code/pathway_enrichment.R` to `state$EF` of the f0 fit | run and summarize | | [ ] |
| B3 | Subtype and stroma concordance (Moffitt basal-like/classical, Moffitt normal/activated stroma, Collisson, Bailey, PurIST) | missing | correlate program scores with published signature scores in training and external cohorts | find or write the signature-scoring code; no repo code found by grep for "PurIST" | | [ ] |
| B4 | Correspondence with the α_F = 0 Programs 3 (protective) and 7 (adverse) | partial | F8 | report \|r\| and shared top genes | | [ ] |
| B5 | Direction of the survival-active programs out of sample | missing | per-cohort Cox on program scores (from F5/F6 data) | needs per-patient external scores | | [ ] |

## 4. Other documents the revamp affects

- **Chapter 1 (bios-dissertation `prelim/literature-review/draft/literature-review-draft.qmd`).**
  These passages describe Chapter 4 in α_F = 0 terms and must be rewritten with the revamp. Round
  3 (2026-11-05) fixed only their prose; exact lines are in
  `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/ch1/` ("Passages for the
  Chapter 4 revamp"):
  - preamble paragraph 3, "The third project, presented in Chapter 4 ...";
  - §1.3.4 "Joint Survival-Supervised Factorization and Its Gaps", closing paragraph;
  - §1.3.6 rank-selection paragraph citing Young2026;
  - §1.3 Summary, closing paragraph ("Chapter 4 presents this framework: a single-modality model
    with external validation ... and a multimodal extension ... with first results on real data").
    The multimodal part is now proposed work.
  - The file's header comment ("Section 1.3 does not describe SSBMF's current implementation
    state") still applies. The literature-review framing about recovering low-variance programs
    must match what the supervised fit shows: survival doesn't move the real loadings yet.
- **Prelim abstract** (`bios-dissertation/prelim/abstract/abstract-text.qmd`; claims sourced in
  `notes/abstract-fact-sheet.md`). The standalone abstract was uploaded 2026-11-02 and fills pages
  iii–iv exactly. Its Chapter 4 sentences (mean 0.627, "slightly higher" than the two-step,
  multimodal first results) no longer describe the chapter. **Andrew decides** whether the
  prelim's copy changes. Any change needs an equal cut and a fact-sheet update.
- **Slides** (`presentation/prelim-slides/OUTLINE.md`, Draft 1, 2026-11-03):
  - the backup slides B1–B9 on the supervised fit move into the main deck;
  - the β prior is Normal, not point-normal;
  - the slides may go beyond the chapter.
- **bios-dissertation full prelim:** `prelim/full-draft/assembly_manifest.json` lists the Chapter 4
  figures. Update it for the new figures, then rebuild (`build_prelim.py`,
  `test_prelim_assembly.py`).
- **ROADMAP / DECISIONS** in both repos: record the revamp, the trimmed material, and the
  supervised model as the chapter's centerpiece.
- **Round 3 human-writer review of Chapter 4** runs last, on the revamped chapter. Use the round 3
  brief in `bios-dissertation/prelim/full-draft/notes/human-writer-review/round3/WRITER-BRIEF.md`.

## 5. Decisions

**Settled by Andrew (2026-11-05):**
- **S-11 session:** it is still running. The revamp works around it (section 1).
- **Branch:** `fix/yfb-supervision`; mirroring to `main` is optional.
- **Multimodal:** proposed work, described as planned methodology with modality-specific priors on
  F. The October multimodal results leave the chapter.
- **Scope of the chapter:** it shows the state at writing time; the slides may go further.

**For Andrew during the follow-up:**
1. F5/F6/B5 need per-patient external risk scores from the fixed model. Is saving them (no refit,
   no retuning) consistent with the "score once" rule, and should the S-11 session produce them?
2. Should the prelim abstract copy change (section 4)?
3. Should F10 (per-cohort β) stay in scope, or move to planned evaluation?

**Checked 2026-11-05:**
- **Prior wording:** point-Laplace on L and F, Normal on β (S-11 runner lines 59–61).
- **Score uncertainty:** the old chapter's aim "score uncertainty in the single-modality
  coefficients" may be partly met already, since the YFB objective uses Var(η) from q(F). The
  follow-up confirms this in the β update before writing it either way.

## 6. Order of work (target: done by about 2026-11-07)

1. Read the latest S-11 session state (section 1). Settle decisions 1–3 with Andrew.
2. Produce T1–T5, F1–F9 and B1–B5 into `paper/prelim/figures/` and a new
   `paper/prelim/make_figures_supervised.R`. Use new dated output folders for anything computed.
   Tick the checklist as each item lands.
3. Move outgoing text to `paper/prelim/notes/trimmed-material-2026-11.md`.
4. Rewrite the chapter to section 2, using the human-writer skill while drafting. Every number
   traces to a file in section 1 or the checklist.
5. Commit the chapter on `fix/yfb-supervision`. Let the hook sync, then rebuild the full prelim
   in bios-dissertation and compare page counts.
6. Rewrite the Chapter 1 passages (section 4) in bios-dissertation.
7. Run the Round 3 review on Chapter 4 (writer, independent reviewer, apply, report).
8. Update ROADMAP/DECISIONS in both repos.
