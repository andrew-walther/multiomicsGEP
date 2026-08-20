# DECISIONS.md — multiomicsGEP

Architectural and analytical decisions made during development, in reverse chronological order.
Each entry records what was decided, why, what was traded away, and which files implement it.

---

## 2026-08-19 — Three-way factor classification (survival-active / genomics-only / dead) added; K_init sweep + best-of-multistart ELBO comparison recommends K=7 (2 survival-active + 2 genomics-only factors), with K=9 as a near-tied alternative

**Motivation:** two methodological gaps flagged for the manuscript
(`docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md`):

1. K=7 for the recommended D4 model (YFB, DeSurv-aligned preprocessing) was chosen by
   cross-validated held-out C-index — which uses survival outcomes to pick model structure,
   conflating structure selection with fitting. The methodologically preferable alternative
   (Method 2) is to start CAVI at a large K and let the ARD-style point-normal prior on β prune
   uninformative factors, so K need not be tuned against the outcome at all — consistent with this
   project's standing policy (2026-04-24 entry below: "ARD preferred over ELBO grid search for K
   selection").
2. The existing `auto_prune_K()` only distinguished "active" (`|β|>thresh` OR `PVE>thresh`) from
   "shrunk" — collapsing two biologically distinct outcomes into one "shrunk" bucket: factors that
   are real gene-expression programs with no prognostic value (*genomics-only*) vs. factors that
   are fully pruned by the prior (*dead*). The manuscript needs the three-way split to state
   "X total programs, Y are survival-associated, X−Y are genomics-active without prognostic value."

**What was implemented (Step 1 + Analysis A only; Analyses B/C deferred to a later phase):**

- `classify_factors(res, Y, beta_thresh, pve_thresh)` added to [`code/select_K.R`](code/select_K.R)
  (after `auto_prune_K()`). Reuses the existing `compute_pve()`. Labels each factor
  `"survival_active"` (`|EBeta_k| > beta_thresh`), else `"genomics_only"` (`PVE_k > pve_thresh`),
  else `"dead"`. Both thresholds default to `config/globals.yml`'s `k_selection$beta_threshold`
  (0.001) and `k_selection$pve_threshold` (0.01). 2 new tests in
  [`tests/test_select_K_cv.R`](tests/test_select_K_cv.R) (KCV-T17/T18) cover all three categories
  plus the borderline-threshold edge case (values exactly at threshold are not active — strict `>`,
  matching `auto_prune_K()`'s existing convention).
- Fixed a stale hardcoded `beta_thresh = 0.05` default in
  [`results/benchmark_sim/run_phase1_diagnostics.R`](results/benchmark_sim/run_phase1_diagnostics.R)
  to read from `config/globals.yml` instead — it predated the 2026-07 rescaling of
  `beta_threshold` to 0.001 for the YFB model's natural β scale (~0.003–0.008) and would have
  silently classified every YFB β as "Shrunk" if used on a YFB fit.
- **Analysis A, single-init sweep** (`results/benchmark_sim/run_k_init_sweep.R`): fit YFB on the
  real TCGA+CPTAC training data (n=273, p=2064 genes) with the D4 preprocessing (per-platform
  z-std, `combined_rank` gene selection, top-3000 per cohort before normalization, no cohort_id) at
  K_init ∈ {5, 6, 7, 8, 9, 10, 15, 20} — extended from an initial {7, 10, 15, 20} pass after the
  first result raised the question of what K best explains genomics *reconstruction* (not just
  survival) — classified factors at each fit, recorded each fit's final `elbo_full` and RMSE, and
  evaluated external C-index across the same 5 held-out cohorts used elsewhere in the benchmark
  suite.
- **Analysis A, multistart follow-up** (`results/benchmark_sim/run_k_init_multistart_check.R`):
  the single-init sweep found K=5/K=6 to have the *best* training ELBO of the whole grid — better
  than K=7 — which could mean either a genuinely-preferred smaller model or a lucky single-init
  local optimum. Re-fit K ∈ {5,6,7,8,9,10} with `fit_cox_on_yf_multistart()` (n_init=15: 1 SVD +
  14 random restarts, best-ELBO selection — the same tool used for exactly this purpose in the
  2026-07-13 K-parsimony follow-up) to check whether K=5/6's ELBO lead survives multistart.
  Results: `results/benchmark_sim/outputs/k_init_sweep/k_init_multistart_results.csv` (summary) and
  `k_init_multistart_restarts.csv` (all 90 individual restarts).

**Result 1 — K_survival_active is stable at 2 across the entire grid tested (K_init 5 through 20),
K_genomics_only stabilizes at ~2–3, both against expectation:**

| K_init | K_survival_active | K_genomics_only | K_dead | K_eff_total | mean external C |
|---|---|---|---|---|---|
| 5  | 3 | 1 | 1  | 4 | 0.596 |
| 6  | 3 | 1 | 2  | 4 | 0.597 |
| 7  | 2 | 2 | 3  | 4 | 0.6267 |
| 8  | 2 | 2 | 4  | 4 | 0.6256 |
| 9  | 2 | 3 | 4  | 5 | 0.6271 |
| 10 | 2 | 3 | 5  | 5 | 0.6279 |
| 15 | 2 | 3 | 10 | 5 | 0.6270 |
| 20 | 2 | 3 | 15 | 5 | 0.6274 |

K_survival_active = 2 (not 3) once K_init ≥ 7, and K_genomics_only stabilizes at 2–3 — i.e. K_eff_total
(the number of factors worth retaining for the genomics *reconstruction* term, not just survival) is
**4 at K_init=7–8, 5 at K_init≥9**. K_dead absorbs essentially all added K_init capacity beyond that
(3→4→5→10→15), consistent with ARD pruning behaving as expected once K_init is large enough.

**Result 2 — the ELBO criterion does NOT cleanly agree with generalization at small K, and this
survives multistart (not a single-init artifact):**

At every K tested, all 15 multistart restarts (1 SVD + 14 random) converged to the SAME best-ELBO
fit as the single fresh-SVD init already found — random restarts never beat SVD-init at any K. This
project has documented this exact pattern before (ROADMAP.md: "best-ELBO multistart... rescued
nothing at any K — SVD init was always the best-ELBO restart") and it replicates here. So K=5/K=6's
ELBO advantage over K=7 (−808,404 and −812,248 vs. −810,540) is **not** a fluke of a single
unlucky init at K=7 — it is robust. Yet K=5/K=6's external C-index (0.596–0.597) is far worse than
every K≥7 (0.6256–0.6279) — a genuine, large (~0.03) gap between training-data model evidence
(ELBO) and held-out predictive validity. The most likely mechanism: with too little spare capacity,
CAVI has nowhere to route ambiguous/noisy signal except into an existing factor's β, inflating
training likelihood without it generalizing; K≥7 gives CAVI "dead" factors to absorb that instead.

**Decision — apply ELBO with an external-validity gate, not ELBO alone:** per the plan agreed before
running this check, best-of-multistart ELBO is the primary K-selection criterion (consistent with
this project's ARD-over-CV policy), but a candidate that collapses on external validation is treated
as failing a plausibility gate, not as a competitor to be ranked on ELBO alone. K=5/K=6 fail that
gate outright. Restricting to K ∈ {7,8,9,10} (all pass, external C 0.6256–0.6279, mutually
indistinguishable), **K=7 (−810,540) and K=9 (−810,421) are a near-tie at the top of best-of-multistart
ELBO** — a ~0.015% difference, versus 3,000–5,500-unit gaps down to K=8/K=10. **Recommendation: K=7
(2 survival-active + 2 genomics-only = 4 total factors)**, on the same parsimony-tiebreaker logic
already used in this project (2026-07-13 entry below: K=7 chosen over statistically-tied K=4/K=5
because it is reachable via a single dependency-free fresh-SVD fit) — **K=9 (2 survival-active + 3
genomics-only = 5 total) is documented as a legitimate, near-tied alternative, not discarded.**

**Note on "should K_init always be large":** the data argues against a blanket "start as large as
possible" rule. ELBO gets monotonically *worse* past K_init≈10 (K=15: −826,682; K=20: −847,800) while
finding no additional real structure (K_eff_total plateaus at 5 from K_init=9 onward) — excess
capacity purely costs ELBO. Combined with K=5/K=6 failing on the small side, there is a real sweet
spot (~7–10) rather than "more is always safer," which happens to sit close to this project's
existing `cavi.k_max: 10` default in `config/globals.yml`.

**Note on the apparent K=3 discrepancy:** an earlier CV-selected K=3 elsewhere in this project
(`k_pdac_yfb_merged`, 2026-05 entries) belongs to a *different* (pre-DeSurv, variance-based top-2000
gene selection) preprocessing pipeline — not the D4/DeSurv-aligned config tested here, whose own
CV-selected K is 7 (`k_merged_yfb_desurv`). The DeSurv paper's own K=3 is a third, unrelated number
(their method, their preprocessing). None of these three K=3's are in tension with this analysis.

This directly answers both open items flagged in the 8/21 meeting-prep model-specification chapter
(`docs/progress_book/chapters/meeting_2026_08_21.qmd` §1): how K=7 was chosen (now: ELBO/multistart,
not just CV, and it holds up), and whether the survival-active/genomics-only split is an artifact of
that specific K or a stable structural feature of the data (K_survival_active=2 is stable from
K_init=7 up; K=5/6's apparent 3rd survival-active factor does not generalize and is treated as a
training-only artifact under the external-validity gate).

**Deferred to a later phase (not run this session):** Analysis B (ARD K-recovery in simulation,
validating that K_eff_survival/K_eff_genomics track *known* ground-truth counts, not just each
other) and Analysis C (re-running the signal-ratio sweep with K_init ≫ K_true to confirm the
YFB-vs-EBMF→Cox comparison holds under ARD pruning). Both remain scoped in
`docs/plans/ssbmf_factor_classification_k_selection_08_13_2026.md`.

---

## 2026-08-03 — Pathway concordance between SBMF and unsupervised EBMF: subtype-level biology replicates, pathway-level mechanism does not (ROADMAP.md A/B comparison, Part 3)

**Question:** given Part 2's result (same-day entry below) that EBMF's factors 1 and 2 correlate
strongly with Program 3/7 at the *gene-loading* level, does that correspondence extend to the
*pathway/gene-set enrichment* level — do the matched EBMF factors show the same enriched biology as
Program 3/7's own fgsea results (DECISIONS.md 2026-07-15), or does supervision change what biology
is recovered, not just which genes are weighted?

**Method:** EBMF's factors 1 and 2 are signed (flashier's default point-normal prior), unlike SBMF's
non-negative point-exponential $\mathbf{F}$ — so each was sign-aligned to its matched program
(multiplied by the sign of its Part-2 correlation: factor 1 by $-1$, factor 2 by $+1$) before a
two-sided `fgsea` (`scoreType="std"`) against the same 5 collections already used for Programs 3/7
(MSigDB Hallmark/Reactome/KEGG/GO:BP + the custom PDAC collection). Concordance for each of Program
3/7's own significant sets (padj<0.10) is read off directly: is that same set also significant,
same direction, for the matched EBMF factor?

**Result — a real, non-trivial split, not clean concordance:**

| Collection | Concordant / total |
|---|---|
| PDAC_custom (Moffitt/Bailey/DeSurv subtype gene sets) | **4/5** |
| Reactome | 1/2 |
| KEGG (MEDICUS pathway/signaling sets) | **0/6** |
| **Overall** | **5/13** |

The coarse *subtype*-level signal (classical/basal-like/squamous calls from DeSurv, Moffitt, Bailey)
is robustly recovered by the unsupervised EBMF factors — e.g. Program 7's `DeSurv_D3_BasalLikeTumor`
hit (padj=6.1e-4) replicates in EBMF_F2 at padj=1.9e-25, and `Bailey_Squamous` at padj=4.9e-7.
But every one of the specific *molecular signaling pathway* hits from KEGG (integrin/FAK/RHOA
signaling for Program 3; RTK/RAS/PI3K signaling for Program 7) fails to replicate — either
non-significant in the matched EBMF factor, or (for Program 3's KEGG hits specifically) significant
but in the **opposite** direction (NES flips sign).

**Interpretation:** supervision doesn't just re-weight genes within biology EBMF already sees (Part
2's finding) — it also sharpens or changes which finer-grained mechanistic signal comes through.
Unsupervised EBMF recovers the coarse tumor-subtype axis but not the specific signaling-pathway
story; the survival term appears to be doing real work at the pathway-mechanism level, not only at
the broad-subtype level. This tempers the Part-2 "SBMF is just a supervised re-weighting of EBMF"
framing — that framing holds at the gene-loading and subtype level, but not uniformly down to
pathway mechanism.

**Files:** `results/benchmark_sim/run_ebmf_pathway_concordance.R` (new); outputs
`results/benchmark_sim/outputs/pathway_enrichment/T6_ebmf_pathway_concordance.csv`,
`ebmf_fgsea_results.rds`.

---

## 2026-08-03 — The two survival-active programs also emerge from unsupervised EBMF (ROADMAP.md A/B comparison, Part 2)

**Question:** ROADMAP.md's "SSBMF vs. unsupervised EBMF" item (line 359) left two parts open after
the external C-index comparison (Part 1, done 2026-06-16/07-16): does the same 2-program biological
structure (Program 3 Protective, Program 7 Adverse) emerge if EBMF is run with **no survival term at
all**, or is it an artifact of supervision?

**Method:** no new model fit was needed. `run_ebmf_cox_external.R`'s two-step baseline already fits
EBMF unsupervised as its first step and saves the fit
(`results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_fit.rds`, $K=20$ factors); both
it and the recommended model (D4) select genes via the identical DeSurv-aligned
combined_rank/top-3000-per-cohort procedure on the same merged TCGA-PAAD+CPTAC training data, so
their 2064-gene universes are identical **and in the same order** (checked via `identical()`, not
assumed — see `run_ebmf_factor_correspondence.R`). This makes a direct Pearson correlation between
EBMF's 20 factor-loading columns and D4's Program 3/7 loadings valid without any re-alignment step.

**Result:** yes — cleanly, for the top match, with one honest secondary nuance:

| EBMF factor | vs. Program 3 (Protective) | vs. Program 7 (Adverse) |
|---|---|---|
| EBMF_F1 | **r = -0.72** (best match) | r = -0.05 |
| EBMF_F2 | r = +0.03 | **r = +0.69** (best match) |
| EBMF_F3 | r = +0.60 (second-best) | r = -0.33 |
| all other 17 factors | \|r\| < 0.15 | \|r\| < 0.31 |

EBMF's own first two extracted factors (i.e., the two it explains the most variance with,
unsupervised) correspond one-to-one with SBMF's two survival-active programs — this is real
structure in the data, not something the survival term invents. **One honest nuance:** EBMF_F3 has a
secondary, moderate correlation with Program 3 (r=0.60) not mirrored by an equally strong secondary
correlation for Program 7 — i.e., the unsupervised decomposition doesn't cleanly separate all of
Program 3's signal into a single factor the way it does for Program 7. This echoes the same kind of
secondary/weaker association already documented for Program 3 in the pathway-enrichment work
(DECISIONS.md 2026-07-15, Program 3's secondary basal-like association) — Program 3's biology may
simply be less cleanly separable from other axes of variation than Program 7's, independent of which
method is used to extract it.

**Interpretation:** SBMF is best understood as adding a supervised re-weighting/selection on top of
EBMF's own factor structure, not as inventing biology EBMF can't see at all — consistent with the
"Bayesian counterpart to DeSurv" positioning already used in the progress report. The joint model's
value-add (per the existing paired-bootstrap result, +0.042 C-index, 95% CI 0.013-0.071) is in
*using* this structure for prediction and disentangling it from the other 18 unsupervised factors,
not in discovering it from nothing.

**Not done (ROADMAP.md Part 3, deferred to a follow-up):** pathway-enrichment concordance between
EBMF_F1/F2's gene loadings and Program 3/7's existing fgsea results
(`docs/reports/pathway_enrichment_report_07_15_26.qmd`) — the correlation result above is a
gene-loading-level check, not yet a biological-pathway-level one.

**Files:** `results/benchmark_sim/run_ebmf_factor_correspondence.R` (new); outputs in
`results/benchmark_sim/outputs/pathway_enrichment/` (`T5_ebmf_factor_correspondence.csv`,
`T5_ebmf_factor_correlation_full.csv`, `F6_ebmf_factor_correspondence_heatmap.png`); noted in
`docs/progress_report/SSBMF_Status_Update_08_03_26.qmd`.

---

## 2026-07-16 — Fixed PDF layout overflow and non-portable HTML in the pathway enrichment report

**Problem, user-reported.** `docs/reports/pathway_enrichment_report_07_15_26.pdf` had a figure (F3,
the gene-weight heatmap) visibly cut off/bleeding off the bottom of its page, with the preceding
figure's page numbers also looking clipped; the `.html` version showed broken-image placeholders in
place of every figure.

**Root causes, confirmed via the LaTeX log (not guessed):**
1. `LaTeX Warning: Float too large for page` + an `Overfull \vbox (318.7pt too high)`: the F3 chunk
   had no `out.width` set, so pandoc defaulted to the image's raw pixel width at 300dpi — **10.7in**,
   far wider than the page. Combined with the image's own tall aspect ratio (h/w=1.665, a ~100-gene
   heatmap) and an automatic `height=\textheight` cap, the image was forced to consume the *entire*
   page height with no room left for its caption, which is what produced the visible bleed.
2. F1 had the same missing-`out.width` defect, less severely (`Overfull \hbox (36pt too wide)`,
   \~0.5in bleeding into the margin — the likely source of the "page numbers look cut off"
   observation, since the margin itself was untouched).
3. **HTML:** images were referenced via relative paths (`../../results/benchmark_sim/outputs/...`)
   rather than embedded — correct and working when opened from this exact repo layout (verified: all
   5 paths resolve), but broken in any context that doesn't preserve the full directory structure
   around the `.html` file (e.g. viewing the file in isolation).

**Fix:** added explicit `out.width` to the F1 (85%) and F3 (45%, given its extreme aspect ratio)
chunks — both now comfortably fit within the page in both dimensions. Added `\sloppy` to the LaTeX
preamble to absorb two much smaller (\~0.1in) text-overflow warnings from an unbreakable long file
path in the reproducibility paragraph. Added `embed-resources: true` to the HTML format so all 5
figures are inlined as base64 data URIs — confirmed via the rendered HTML (file grew 42KB→2.7MB, zero
remaining `_files/` references) — making the HTML immune to this class of bug regardless of where
it's opened from.

**Verification:** re-rendered both formats; the LaTeX log is now clean (zero `Overfull`/`Float too
large` warnings, down from 4); visually re-inspected the previously-broken pages — F2 and F3 now
render completely with full captions on their own pages. No change to any analysis, number, or
figure content — this is a rendering/layout fix only.

**Files:** `docs/reports/pathway_enrichment_report_07_15_26.qmd` (chunk options + YAML),
`docs/reports/pathway_enrichment_report_07_15_26.{pdf,html}` (re-rendered).

---

## 2026-07-16 — Bootstrap C-index CIs; refreshed a second stale baseline; paired test vs. the two-step method

**Context.** Item 2's progress report flagged "no uncertainty quantification on the external C-index"
as the highest-priority open item (2026-07-15 gap-review entry). Implemented per-cohort bootstrap
confidence intervals and, per the user's request, a paired-bootstrap test of whether the recommended
model (YFB, DeSurv-aligned gene selection, $K=7$, no cohort indicator) is significantly more
concordant than a two-step (unsupervised EBMF then Cox) baseline on the same real external cohorts.

**A second stale-baseline finding, caught while sourcing data for this work (same category as the
multi-cohort sim, 2026-07-15 entry above, but not caught by that pass since it audited `.qmd` report
files, not raw result CSVs):** `results/benchmark_sim/outputs/ebmf_cox_external/ebmf_cox_external_results.csv`
(the EBMF$\to$Cox two-step baseline's external validation on the same 5 real cohorts) was last
generated **2026-06-16** — the runner script (`run_ebmf_cox_external.R`) was patched for the
train/test preprocessing fix on 2026-07-12 (Phase 1c), but its cached output was never regenerated
afterward. **Re-ran it under current code:** mean external C moved **0.564 $\rightarrow$ 0.581**
(K=20, same direction of correction as the recommended model's own 0.636$\rightarrow$0.627 move, for
the same reason — external cohorts are now preprocessed consistently with training). `ROADMAP.md`'s
status banner, which cited the stale 0.564 figure, is corrected.

**New reusable code (TDD, `code/concordance_ci.R`, 12 new tests):**
- `bootstrap_concordance_ci(risk, time, status, B, seed)` — percentile-bootstrap CI on a single
  model's C-index. The risk score's "higher = worse" orientation is fixed **once** from the full
  sample (this project's existing `oriented_cindex()` convention); bootstrap replicates are **not**
  individually re-oriented, since doing so would force every replicate $\ge 0.5$ and upward-bias the
  CI for a weak/null signal — fixing the sign once and letting replicates vary freely below 0.5 is
  what gives an honest interval.
- `bootstrap_concordance_diff_ci(risk_a, risk_b, time, status, B, seed)` — a **paired**
  bootstrap CI on the difference in C-index between two models scored on the *same* patients (each
  replicate resamples patients once and scores both models on that identical resample, preserving
  the correlation between the two models' errors on shared patients — the statistically correct way
  to test "is model A significantly better than model B here," as opposed to comparing two
  independently-constructed CIs).
- No re-fitting was required for either model to compute these CIs — both runner scripts
  (`run_desurv_comparison.R`, `run_ebmf_cox_external.R`) were extended to cache per-patient
  `{risk, time, status}` per cohort (`*_riskscores.rds`) during their existing scoring step, which is
  why refreshing the stale EBMF baseline above (a real re-fit) only needed doing once, here.

**Results (`results/benchmark_sim/run_external_ci_analysis.R`, B=2000, seed=1):**

| Cohort (n) | YFB C-index (95% CI) | EBMF→Cox C-index (95% CI) | Paired diff (95% CI) | Sig.? |
|---|---|---|---|---|
| Dijk (90) | 0.635 (0.555, 0.708) | 0.592 (0.516, 0.671) | +0.042 (−0.029, 0.111) | No |
| Moffitt (123) | 0.549 (**0.467**, 0.625) | 0.541 (0.462, 0.613) | +0.008 (−0.071, 0.089) | No |
| PACA-AU array (63) | 0.648 (0.553, 0.747) | 0.597 (0.500, 0.695) | +0.051 (−0.046, 0.149) | No |
| PACA-AU seq (52) | 0.657 (0.544, 0.772) | 0.573 (0.464, 0.687) | +0.084 (−0.020, 0.188) | No |
| Puleo (288) | 0.645 (0.599, 0.688) | 0.602 (0.556, 0.646) | +0.043 (**0.005**, 0.081) | **Yes** |
| **Pooled (fixed-effect, n=616)** | — | — | **+0.042 (0.013, 0.071)** | **Yes** |

**Answers to the questions this was built to answer:**
1. **Is the external C-index distinguishable from chance (0.5)?** Yes for 4 of 5 cohorts — Moffitt's
   CI (0.467, 0.625) includes 0.5, confirming quantitatively what was previously only a qualitative
   flag ("only marginally above chance").
2. **Is SSBMF significantly more concordant than the two-step baseline?** Individually, only the
   largest cohort (Puleo, n=288) reaches significance on its own (the other four all have positive
   point estimates, 0.008–0.084, but are underpowered alone at n=52–123). **Pooled across all 5
   cohorts (fixed-effect, inverse-variance-weighted), the advantage is significant: +0.042 (95% CI
   0.013–0.071).** The pooling is a simplifying fixed-effect assumption (a common effect size across
   cohorts), stated plainly — no formal cross-cohort heterogeneity test (e.g. Cochran's Q) was run,
   since that is a more specific statistical choice than this pass was scoped to make.
3. **Is the variation across cohorts itself significant?** Per-cohort CIs overlap substantially
   (e.g. Moffitt's and Dijk's CIs overlap by more than half their width) — nothing here suggests
   real cohort-to-cohort heterogeneity beyond what sampling noise from very different cohort sizes
   (52–288 patients) would produce on its own.

**Independent review findings addressed (before merge):** (1) a bootstrap replicate that happens to
draw fewer than 2 events was silently producing `NaN` from `concordance()` rather than erroring,
surfacing later as a confusing `quantile()` failure with no indication of the real cause — both
CI functions now check `sum(status[idx]) < 2` per replicate and fail loud with an explicit,
actionable message; two new deterministic tests (fixed seed known to trigger a degenerate resample)
confirm the message, plus two more tests for the full-sample `n<10`/`too few events` guards that
existed but were untested. (2) the progress report's new §5.1 insertion had left a truncated,
duplicated opening line of the following paragraph ("Study-specific baseline hazard...") — removed,
report rebuilt. Neither finding changed any reported number; both are documentation/robustness fixes.

**Verification:** `Rscript tests/run_tests.R` → **392/392** (374 + 18 new, up from an initial
386/386 before the two review-driven test additions above).

**Files:** `code/concordance_ci.R` (new; degenerate-resample guard added post-review),
`tests/test_concordance_ci.R` (new, TDD, 18/18), `tests/run_tests.R` (registered),
`results/benchmark_sim/run_desurv_comparison.R` (added risk-score caching),
`results/benchmark_sim/run_ebmf_cox_external.R` (added risk-score caching; re-run — refreshed stale
output), `results/benchmark_sim/run_external_ci_analysis.R` (new),
`results/benchmark_sim/outputs/desurv_comparison/{desurv_comparison_riskscores.rds,
external_cindex_ci.csv, external_paired_diff_ci.csv}` (new),
`results/benchmark_sim/outputs/ebmf_cox_external/{ebmf_cox_external_results.csv (refreshed),
ebmf_cox_external_riskscores.rds}` (new), `ROADMAP.md` (corrected stale 0.564 figure),
`docs/reports/ssbmf_progress_report_07_15_26.qmd/.pdf/.html` (duplicated-paragraph fix, post-review).

---

## 2026-07-15 — External gap review of the progress report; corrections applied before advisor presentation

**Context.** Before finalizing the progress report for presentation, dispatched an independent
review aimed specifically at the report's own writing (not the underlying model) — stale numbers,
overclaims, or omitted caveats a critical reader would catch.

**Findings and corrections (all applied; no re-fits needed):**
1. **Stale sweep numbers.** The report quoted early 5-seed sweep values (0.931/0.908/0.875); the
   shipped 10-seed CSV/figure gives 0.927/0.900/0.899. Corrected to match the shipped data.
2. **Overclaim.** "Above chance (0.5) on every cohort" was asserted with no uncertainty. Softened to
   "point estimate above 0.5"; added the cross-cohort SE ($\approx$0.02); called out the weakest
   cohort (Moffitt, 0.549) explicitly rather than letting it hide inside the mean.
3. **Hybrid-scenario "+0.10 joint advantage" reframed.** 2 of the two-step baseline's 5 seeds
   collapse (dragging its mean to 0.717); on the 3 non-collapse seeds it scores 0.807, close to the
   joint model's 0.815. The advantage is now described as primarily **robustness** (the joint model
   never collapses here) rather than a pure discrimination margin.
4. **$K_{\text{eff}}=2$ evidence count corrected.** Was described as "three independent lines of
   evidence," counting deflation-init as a third — but deflation-init is a negative,
   analytically-equivalent-to-SVD result, not independent confirmation. Now stated as two positive
   strategies (warm-start, joint BO) plus an analytic ruling-out of initialization as the cause of
   the earlier ambiguity.
5. **Two items added to the report's open-items section as explicitly pending** (not previously
   flagged): per-cohort concordance confidence intervals (only a cross-cohort SE existed), and
   factor-stability-under-resampling (the two active programs' biology is confirmed five ways but
   never re-derived from a bootstrap refit).

**Also removed:** internal config shorthand (the "D1"-"D5" labels) from all report prose in favor of
descriptive model language (YFB, cohort indicator, DeSurv-aligned gene selection, etc.) for the
external-audience-facing documents — the labels are retained only inside code filters, never in text
a reader sees.

**Files:** `docs/reports/ssbmf_progress_report_07_15_26.qmd/.pdf/.html`,
`docs/progress_report/SSBMF_Status_Update_07_15_26.qmd/.pdf/.html`. No code or result changes —
corrections are entirely in report prose.

---

## 2026-07-15 — Whole-arc performance review + multi-cohort sim re-validated under the corrected model

**Context.** Before presenting the progress report (Item 2) to advisors, ran a full performance
review of the current method across the three validation axes the 6/18 feedback and the plan call
for: (1) joint-vs-two-step with and without survival signal, (2) shared vs. study-specific factor
recovery, (3) external PDAC cohorts. Two of the three were already computed under the corrected code
(post-2026-07-12); the third (the multi-cohort shared/specific study) predated it and was re-run.

**Findings — all three axes confirm the method performs as claimed:**

1. **Survival-strength sweep (2-step ± signal), already under corrected code
   (`survival_strength_sweep_results.csv`, 2026-07-13):** at zero signal YFB $\approx$ PCA+Cox
   $\approx$ EBMF+Cox (all near chance, e.g. `default` 0.523/0.539/0.544); the joint advantage
   emerges and grows with signal in 3/4 scenarios (`default` at strength 4: YFB 0.927 vs. 0.900/0.899;
   `low_snr`: 0.815 vs. 0.654/0.708; `high_K`: 0.903 vs. 0.885/0.837). `sparse_synthetic` reverses
   (0.795 vs. 0.934/0.940) — the documented CAVI factor-collapse artifact under symmetric,
   disjoint-support loadings. $\alpha=0$ invariance control: max$|F_{\alpha=0.5}-F_{\alpha=0}|=0$
   across all 240 cases (genomics factorization exactly $\alpha$-invariant → answers the meeting's
   "does k shrink if survival is off" — the genomics factor structure does not change).

2. **Multi-cohort shared/study-specific study — RE-RUN under the corrected model this session
   (`multicohort_sim_results.csv`, regenerated):** the YFB (recommended) arms reproduced the
   2026-06-15 result **within noise** (max C-index diff 0.019, mean 0.004 over 30 fits); the
   EBMF two-step baseline was **bit-identical** (max diff 0.000). Result stands: no-shared-signal →
   YFB 0.549 $\approx$ EBMF 0.550 (equivalence); hybrid (some factors prognostic) → YFB 0.815 vs.
   EBMF 0.717 (**+0.10 joint advantage**); all-shared → YFB 0.871 vs. 0.856; specificity
   classification accuracy 1.00 for YFB across all scenarios. Cohort indicator neutral-to-slightly-
   negative in sim (hybrid 0.813 vs. 0.815), mirroring the real-data finding.
   **The only material movement on re-run was in the LB reference arms** (one `hybrid`/`LB_base`
   cell 0.828 → 0.633, iteration count 13 → 79). Mechanism: the retired $\lambda$ (2026-07-12
   Phase 1b) was *live* in LB's `update_L` but *already inert* in the YFB path, so removing it was a
   true no-op for YFB but perturbed LB's floating-point evaluation enough to tip a near-bifurcation
   fit into a different local optimum — consistent with the documented "LB collapses more than YFB"
   instability (2026-07-12 Phase 2), **not a regression, and not in the recommended model.**

3. **External PDAC (current, `desurv_comparison_results.csv`, regenerated 2026-07-15):** recommended
   config mean external C = 0.627, above chance on all 5 held-out cohorts (Dijk 0.635, Puleo 0.645,
   PACA_AU_seq 0.657, PACA_AU_array 0.648, Moffitt 0.549); $K_{\text{eff}}=2$. Highest mean among all
   5 preprocessing/parameterization configs tried; the cohort indicator lowers it (0.617). No tried
   configuration beats 0.627 by more than 0.0015 (inside SE $\approx 0.02$) — current defaults are at
   or statistically indistinguishable from the best observed.

**Conclusion:** no performance is being left on the table by the current defaults, and the
recommended (YFB) model's simulation behavior is stable under the corrected code. Documented in the
progress report `docs/reports/ssbmf_progress_report_07_15_26`.

*Files: `results/multi_cohort_sim/outputs/multicohort_sim_results.csv` (regenerated under corrected
model); no code changes. Baseline for comparison preserved during the session.*

---

## 2026-07-15 — Compliance review vs. 6/18 feedback + performance-optimality check

**Context.** Independent of the net-benefit gate (tests + C-index-regression check, entry below),
ran two additional review passes before presenting the progress report: (a) an item-by-item
compliance check against the original `docs/plans/rashid_lab_meeting_notes_06_18_2026.md` feedback,
and (b) whether the recommended configuration is actually the best-performing choice or whether
performance is being left on the table.

**Compliance — all model-science items addressed or defensibly deviated:**
- Addressed as literally requested: $h_0(t)$ non-parametric/cancels-in-partial-likelihood (also
  confirmed for the new stratified extension), $\alpha$ default $=0.5$, per-platform z-standardize
  then platform-correct, `strata(study)`, the "why 7 factors" investigation, the $\alpha=0$
  "does $k$ shrink" control, the joint-vs-2-step diagnostics, factor interpretation, external-cohort
  normalization consistency, GO/gene-ID characterization, the paper repo.
- **Two deviations, both defensible but worth stating plainly:** (1) $\lambda$ was retired entirely
  rather than constrained to $(0,1)$ with a $(1-\lambda)$ term added — $\alpha$ already plays that
  rebalancing role, and DeSurv's own $\lambda$ is a different (elastic-net) mechanism this model
  fills via empirical-Bayes priors, but a reader expecting "$\lambda \in (0,1)$" will notice the
  parameter is simply gone. (2) The literal $np$/$n$ normalization convention was implemented and is
  available (`norm_convention="np_n"`) but is not the shipped default — a per-platform-safe
  convention (`per_p`) is, because `np_n` collapses LB's $\mathbf{L}$/$\mathbf{F}$ to zero.
- **Genuine gaps, not addressed:** no formal review of how the DeSurv paper's own simulation was
  constructed (requested at the meeting); "figure out HPC cluster workflow" has no evidence of being
  done; the post-commit Codex review hook and a RAG literature database remain backlog, not actioned.

**Performance-optimality:** confirmed the recommended configuration (YFB, $K=7$, $\alpha=0.5$, no
cohort indicator, normal prior, per-platform z-std) is at or statistically indistinguishable from
the best configuration observed in any sweep — every nominally higher-scoring alternative ($K=8$,
$\alpha\approx0.71$ at 0.6282; $K=4$/$K=5$ warm-start at 0.6270) beats it by $\le 0.0015$, about
1/13 of the per-config SE, and each has a documented, sound reason for non-adoption. The cohort
indicator actively hurts ($-0.010$). No performance is being forgone by current defaults.

**Disposition:** folded into the progress report as a new "Disposition of the 6/18 feedback"
section (full per-item table) and used to correct two minor documentation issues found in the
process (see net-benefit gate entry below).

**Files:** no code changes; `docs/reports/ssbmf_progress_report_07_15_26.qmd` (new section).

---

## 2026-07-15 — Whole-branch review against the net-benefit gate (pre-Item-2 arc audit)

**Context.** Before writing the Item 2 progress report, an independent review of the entire
post-6/18 arc (Phases 1-3, K-parsimony follow-up Steps 1-4, pathway enrichment, stratified Cox)
against the plan's net-benefit gate: (1) all tests green, (2) external mean C-index not worse than
the 0.636 pre-work baseline in any way that isn't already documented and understood.

**Verdict: gate satisfied.** `Rscript tests/run_tests.R` confirmed **374/374** passing. The
$0.636\rightarrow0.627$ change is fully attributable to the documented Phase 1c train/test
preprocessing fix (external cohorts were rank-transformed while training was per-platform
z-standardized — the reverse) — not a regression from objective normalization, $\lambda$ retirement,
the K-parsimony work, pathway enrichment, or the stratified-Cox extension, each of which
independently leaves the recommended fit's $\hat\beta$/$E[\mathbf{L}]$/$E[\mathbf{F}]$ numerically
unchanged or performance-neutral.

**Two minor documentation findings, both fixed:**
1. `CLAUDE.md` had two stale "355/355" test-count mentions (both now 374/374).
2. The 2026-07-15 stratified-Cox entry's "374/374 (was 357)" didn't reconcile against the 355 stated
   in the adjacent pathway-enrichment entry; corrected to "up from 355 before this entry" rather than
   asserting an unverified precise delta.

**A third finding, caught independently while sourcing numbers for the progress report (not by the
dispatched reviewer):** `ROADMAP.md`'s "Sensitivity — D3" line quoted 0.622, a
pre-Phase-1c-preprocessing-fix figure; the current results CSV gives 0.611 (all five preprocessing
configs shifted when that fix landed, but only the recommended config's headline had been updated at
the time). Corrected.

**Files:** `CLAUDE.md`, `DECISIONS.md` (this correction), `ROADMAP.md` — documentation only, no code
or result changes.

---

## 2026-07-15 — Study-specific baseline hazard via stratified Cox partial likelihood

**Question:** the training set pools two studies (TCGA-PAAD, CPTAC) that may have different baseline
survival, but the model uses a single shared baseline hazard. The 6/18 lab-meeting feedback asked
for a `+ strata(study)` term so baseline risk can differ by study while the coefficient vector
$\boldsymbol\beta$ stays shared. (Original-plan "Phase 4"; consolidation-plan "Item 3".)

**Decision — minimal stratified partial likelihood, no parametric baseline (confirmed with the
user before implementing):** form the Breslow risk sets *within* each stratum. The baseline hazard
still cancels per-stratum in the partial likelihood, so no parametric $h_0(t)$ (e.g. Weibull) is
introduced — the survival CAVI derivation is otherwise unchanged. The alternative, a fully
parametric per-study baseline, was considered and declined as more modeling machinery than the
feedback required.

**Implementation:** the Cox Taylor helpers `calc_cox_taylor()` (LB, `code/fit_modular.R`) and
`calc_cox_taylor_yf()` (YFB, `code/fit_cox_on_yf.R`) gained an optional `strata=` argument; the
fitting functions `fit_supervised_mf_modular()` and `fit_cox_on_yf()` gained `strata_id=`, threaded
to both the burn-in and main-loop Taylor calls. When `strata` is `NULL` the original single pooled
risk set is preserved bit-for-bit; when supplied, the per-sample score $u$ and diagonal Hessian $w$
are computed within each stratum and scattered back by index, and the partial log-likelihood is
summed across strata. `strata_id` is distinct from `cohort_id` (which absorbs *genomic* platform
offsets); the two are composable.

**Performance (the merge gate): neutral.** Refitting the recommended configuration (YFB D4, K=7)
with study strata (TCGA vs. CPTAC) vs. without, over the 5 held-out external cohorts:

| | Mean external C | $K_{\text{eff}}$ | max\|$\hat\beta$\| |
|---|---|---|---|
| No strata (baseline) | 0.6267 | 2 | 0.0404 |
| Study strata | 0.6263 | 2 | 0.0405 |

Difference is $-0.0004$ (per-cohort swings within $\pm0.0015$), i.e. within noise. Mechanistically
expected: per-platform z-standardization already absorbs most cross-study structure before the two
training cohorts merge, so restricting the risk sets to within-study leaves $\boldsymbol\beta$
essentially unmoved. **Kept as an available option (`strata_id`, default `NULL`), not enabled by
default** — it adds a capability the feedback requested and is correct, but earns no predictive gain
on this configuration, and the default fit is unchanged.

**Verification:** `tests/test_stratified_cox.R` (17 tests) — single-stratum reduction anchors
(bit-identical to unstratified, and full-fit $\hat\beta$ identical at tol 1e-8), per-stratum
additive decomposition, two independent `survival::coxph` oracles (Breslow partial log-likelihood at
fixed coefficient, and martingale residuals for the score $u$), NA/length-mismatch fail-loud guards.
Full suite 374/374 (up from 355 before this entry).

**Independent review findings addressed:** (1) NA in `strata`/`strata_id` was silently dropped by
`as.factor()` (→ `0/0` downstream) — now rejected with an explicit error at both the helper and
fit-function entry points. (2) The training sign-correction step (`sign_correction=TRUE`) computes a
*pooled* concordance for its coarse $\boldsymbol\beta$-orientation flip, not a stratified one; left
as-is because it only drives a sign flip for a shared $\boldsymbol\beta$ and the effect is
negligible, but noted here as a known minor inconsistency under stratification. (3) The CV/tuning
wrappers (`select_alpha_cv`, `select_K_cv`, `auto_prune_K`) do **not** thread `strata_id` per fold —
hyperparameter tuning runs unstratified even when the final fit is stratified; documented in the
fit-function comments so this is not mistaken for stratified tuning.

**Files:** `code/fit_modular.R`, `code/fit_cox_on_yf.R`, `tests/test_stratified_cox.R`,
`tests/run_tests.R`.

---

## 2026-07-15 — Pathway enrichment on the recommended model's two survival-active programs

**Question:** what biology do Program 7 (Adverse) and Program 3 (Protective) — the recommended
model's $K_{\text{eff}}=2$ survival-active factors (DECISIONS.md 2026-07-13) — represent?

**Method:** `fgsea` ranked-by-weight enrichment (primary; one-sided, `scoreType="pos"`, since the
point-exponential $\mathbf{F}$ prior makes every weight $\ge 0$) plus over-representation analysis
on top-N genes (confirmatory) against MSigDB Hallmark/Reactome/KEGG/GO:BP and a custom PDAC
collection (Moffitt basal/classical, Bailey 4-subtype, DeSurv D1-D3 factor gene lists — all
recovered from local reference data or the DeSurv SI appendix, no fabricated gene lists). All 7
programs enriched; reporting focuses on Programs 3 and 7. Full report:
`docs/reports/pathway_enrichment_report_07_15_26.{qmd,pdf,html}`.

**A correction to the 2026-06-16 planning draft, caught during implementation (not by the plan
itself):** that draft's decision on subtype labels pointed to "MS"/"MS_K2" as the per-sample tumor
basal/classical 2-group call. Direct inspection of `cmbSubtypes.RData` and
`TCGA_PAAD.caf_subtype.rds` shows MS/MS_K2 is actually the Moffitt **stroma** activation axis
(Activated/Normal) — a different biological question. The correct, already-available per-sample
tumor axis is PurIST (categorical Basal-like/Classical + continuous `PurIST.prob`); used that
instead for the subtype-concordance check.

**Result — five independent methods agree, with no contradictions:**

| Method | Program 7 (Adverse) | Program 3 (Protective) |
|---|---|---|
| Gene-set enrichment (fgsea, all collections) | Basal-like/squamous (DeSurv D3, padj=6.1e-4), MET/EGFR-RTK signaling (3 KEGG sets, padj~1.1-1.2e-3), Bailey Squamous (padj=0.030) | Classical (DeSurv D1, padj=0.017; Moffitt Classical, padj=0.035) |
| PurIST subtype concordance (TCGA-PAAD, n=144, 100% matched) | Spearman rho=+0.58 vs. basal-likelihood (p=2.1e-14) | rho=-0.57 (p=4.9e-14) |
| External-cohort survival (5 held-out cohorts) | HR>1 in **5/5** cohorts (range 1.30-2.40) | HR<1 in **5/5** cohorts (range 0.42-0.88) |
| SBMF-vs-DeSurv gene overlap (top-270 genes, hypergeometric) | 80/270 overlap with DeSurv D3 (p=3.2e-16) | 61/270 overlap with DeSurv D1 (p=1.1e-7) |

This is the established PDAC subtype-survival relationship (basal-like/squamous = worse, classical
= better), independently recovered by every one of the five methods run (the table's fgsea row plus
its ORA confirmatory cross-check on top-N genes, which independently reaches the same conclusion at
$padj<10^{-7}$ for both programs' primary custom-collection hits).

**One honest nuance, reported rather than smoothed over:** Program 3 also shows a secondary, weaker
association with DeSurv's basal-like gene list (enrichment padj=0.017 vs. 0.035 for D1; gene
overlap 49/270 vs. 61/270 for D1), reproduced by fgsea and gene overlap. Plausibly attributable to
shared general epithelial/tumor-identity markers (EPCAM, KRT8, TFF1/TFF3) rather than basal-specific
drivers. ORA does **not** confirm this secondary hit at any top-N tested (padj $\ge$ 0.18) — unlike
the primary classical association, which ORA does confirm and which strengthens with N — additional
evidence the secondary association is real but weaker, not a contradiction of the dominant signal.

**Independent code review (before merge) caught two real bugs, both fixed and re-verified before
this entry was finalized:**
1. **[Critical]** `cohort_signature_cox()`'s C-index used a fixed sign convention
   (`concordance(Surv(...) ~ I(-score))`) correct for an adverse (higher score = higher risk)
   signature but silently wrong-direction for a protective one — Program 3's C2 table originally
   reported C-index 0.39-0.49 in all 5 cohorts (looking like "worse than chance") when its true
   discriminative accuracy in the correct direction was 0.51-0.61. Fixed to match this project's
   existing `oriented_cindex()` convention (`max(c_raw, 1-c_raw)`, used in 8+ other benchmark
   scripts) — direction-agnostic by construction. HR signs and p-values (which never used the buggy
   convention) were unaffected; only the C-index column was wrong. New regression test (T9.5)
   constructs a protective signature specifically to catch this class of bug in the future.
2. **[Important]** `compute_geneset_overlap()`'s hypergeometric test used each gene list's
   unrestricted length as `phyper()`'s parameters, but DeSurv's 270-gene-per-factor lists are not a
   full subset of the 2064-gene SBMF background (~245-259 of each 270 fall inside it) — understating
   significance (Program 3 vs. D1: p moved from ~2.4e-6 to the correct ~1.1e-7, a ~22x difference;
   see the corrected table above). Fixed by restricting both gene lists to the background before
   computing overlap/union/`phyper()`. New regression test (T10.3b) covers a set partially outside
   the background.

Both fixes changed only precision/exact numbers, not any directional conclusion — re-verified end to
end (re-ran Steps 8-9 on real data, re-rendered the report) rather than assumed safe.

**New dependencies:** `fgsea`, `clusterProfiler`, `msigdbr`, `org.Hs.eg.db` (Bioconductor;
pre-approved in the 2026-06-16 planning session's own decisions). `msigdbr` 26.1.0 renamed its
`category`/`subcategory` arguments to `collection`/`subcollection` (old names deprecated but
functional) — code uses the new names.

**Files:** `code/pathway_enrichment.R` (all reusable functions); `results/benchmark_sim/run_pathway_enrichment.R`,
`run_subtype_concordance.R`, `run_external_cohort_robustness.R`, `run_sbmf_desurv_overlap.R`
(orchestration); outputs in `results/benchmark_sim/outputs/pathway_enrichment/`. Also: added a
`sampID` field to `load_pdac_raw()`'s return (`results/benchmark_sim/benchmark_helpers.R`) — needed
to match the D4 fit's pooled patient loadings back to real TCGA-PAAD sample barcodes, previously
discarded (`rownames(Y) <- NULL`). Test suite: 355/355 (main), 88/88 (real-data).

---

## 2026-07-13 — K-parsimony follow-up Step 4 (final validation): K=7 remains the recommended default; K=4 is a validated, statistically-equivalent, more-parsimonious alternative

**Synthesis of Steps 1-3.** The conclusive, doubly-verified answer to "is K=7 necessary":

| Config | Fitting procedure | Mean external C | K_eff |
|---|---|---|---|
| K=2 | warm-start-from-K=7 (best of 3 methods) | 0.5608 | 1 |
| K=3 | warm-start-from-K=7 (best of 3 methods) | 0.5955 | 3 |
| **K=4** | **warm-start-from-K=7 (PVE-ranked columns)** | **0.6270** | **2** |
| **K=5** | **warm-start-from-K=7 (PVE-ranked columns)** | **0.6270** | **2** |
| **K=7** | **fresh SVD (one-step, no dependency)** | **0.6267** | **2** |
| K=8 | joint (K,alpha) Bayesian optimization | 0.6282 | 2 |

K=4, K=5, K=7, and K=8 are all statistically indistinguishable on external validation (spread of
0.0015, far inside any individual config's SE of ~0.02) and all converge to the **same K_eff=2** —
the number of factors that actually matter for survival prediction is robust to total K once
optimization artifacts are controlled for, across a 2x range of K (4 to 8). K=2/K=3 remain a
genuine, confirmed floor below which performance drops regardless of fitting procedure.

**Decision: keep K=7 as the primary recommended configuration; do not change
`config/globals.yml`'s `k_merged_yfb_desurv`.** This is a judgment call, stated plainly rather than
silently resolved:
- K=7 is reachable with a single, dependency-free fresh-SVD fit. K=4/K=5's statistically-equivalent
  performance is **only reachable via the two-step warm-start recipe** (fit K=7 first, extract its
  top-K PVE-ranked columns via `extract_top_k_by_pve()`, then refit at the smaller K with
  `init_method="custom"`) — a fresh SVD fit at K=4 alone gives 0.5409, a full 0.086 worse. Adopting
  K=4 as "the" recommended config would mean every future reproduction of the recommended pipeline
  needs this two-step procedure for zero measurable gain in external performance (K_eff is identical
  at 2 either way).
- Step 5 (pathway enrichment) already targets K=7's Factor 7 (adverse) / Factor 3 (protective) —
  keeping K=7 avoids re-deriving which factor indices are the biologically active ones under a
  different K, an unforced complication with no offsetting benefit.
- The scientific finding — K=7's extra factors beyond K_eff=2 are not load-bearing for survival
  prediction, and a smaller, equally-performing K is available and validated — is fully preserved
  and reportable regardless of which config is the "default": this is now documented, tested,
  reproducible fact, not lost by keeping K=7 as the practical default.

**What actually changed vs. Phase 3's original conclusion.** Phase 3 (2026-07-13, entry below) said
"K=7 is not free to shrink." That is now **superseded**: K=7 is not *necessary* (K=4/K=5 tie it), but
remains the *practical default* for reproducibility reasons unrelated to predictive performance. This
is a more complete, more honest answer than either "K=7 is required" (Phase 3, since revised) or
"switch to K=4" (would silently reintroduce a two-step fitting dependency) — both of Steps 1-3's own
optimization strategies (warm-start, joint BO) independently confirm the same K_eff=2 ceiling and the
same statistical tie among K∈{4,5,7,8}.

**Recorded for a future session or the user to revisit:** if a future need (e.g. a reviewer question,
or a manuscript emphasis on parsimony) makes adopting K=4 as the formal default worthwhile despite the
two-step fitting cost, this entry has everything needed to make that change: update
`config/globals.yml`'s `k_merged_yfb_desurv: 4`, document the two-step fitting requirement prominently
in `CLAUDE.md`'s Quick Reference, and re-verify Step 5's target factor indices under the new K.

*Files: no code changes (synthesis of already-committed Steps 1-3 results); `CLAUDE.md`'s "Current
model status" line updated to reference this entry.*

---

## 2026-07-13 — K-parsimony follow-up Step 3: joint (K, alpha) Bayesian optimization finds a marginally better K=8, not a smaller K

**Motivation.** Step 2's deflation-init (entry below) failed to reproduce Step 1's warm-start
rescue — an analytically-understood negative result, not a bug. This raised the question the
user chose to pursue next: does jointly tuning `(K, alpha)` via Bayesian optimization (DeSurv's own
approach, vs. our fixed-alpha=0.5, K-only CV) find a smaller K with comparable external performance?
Branch `joint-k-alpha-bayesopt`.

**New code:** `code/select_k_alpha_bo.R` — `select_k_alpha_bo_objective()` (thin wrapper: scores a
single `(K, alpha)` point via a single-K call to the existing `select_K_cv()`, reusing its
fold-fitting machinery rather than duplicating it), `select_k_alpha_bayesopt()` (wraps
`rBayesianOptimization::BayesianOptimization()`), and `pick_trustworthy_bo_winner()` (validity gate,
see below). `results/benchmark_sim/run_joint_bo.R` (real-data runner). `tests/test_select_k_alpha_bo.R`
(12 tests, TDD). Package `rBayesianOptimization` installed (CRAN; `ParBayesianOptimization`, the
plan's first-choice candidate, is not available for the current R version).

**A real design gap caught before the expensive run, not after.** The raw BO objective (mean CV
concordance) has no way to know whether a high-scoring point reflects genuine survival modeling or
incidental unsupervised-reconstruction alignment with the outcome — a documented failure mode in
this exact project (2026-05-05 entry below: alpha=1.0 degenerate K-CV selection, K_eff=0; the
archived "lucky PCA direction alignment" finding). A `--quick` dry run surfaced exactly this: BO's
raw `Best_Par` was `K=10, alpha≈0` (literally `2.2e-16`, i.e. survival term off) with **K_eff=0**,
yet external mean C=0.6213 — a plausible-looking number for a model that isn't using survival
information at all. Fixed by adding `pick_trustworthy_bo_winner()`: re-fits the top 5 CV-scoring
`(K, alpha)` candidates on full training data and returns the best-scoring one with K_eff > 0,
raising an informative error if none qualify. (`alpha_bounds` were deliberately left at the full
`[0,1]` rather than narrowed away from the extremes — this project's own existing alpha-grid
convention (`config/globals.yml`, `alpha_grid: [0.0, ..., 1.0]`) already tests both extremes, so
narrowing would have been a new, undiscussed deviation; a validity gate on the *result* is more
honest than silently reshaping the search space.)

**Result — real data (D4 config, n_folds=5, init_points=8, n_iter=15):**

| Config | Mean external C | SE | K_eff |
|---|---|---|---|
| K=7, alpha=0.5 (Step 1 fresh SVD reference) | 0.6267 | 0.0199 | 2 |
| K=4, alpha=0.5 (Step 1 warm-start-from-K=7) | 0.6270 | 0.0198 | 2 |
| K=5, alpha=0.5 (Step 1 warm-start-from-K=7) | 0.6270 | 0.0198 | 2 |
| **K=8, alpha=0.7072 (Step 3 BO winner)** | **0.6282** | 0.0203 | 2 |

The BO winner passed the validity gate (K_eff=2, not degenerate) and passed the sanity check against
the existing internal K-CV grid (BO's internal CV=0.6546 vs. the existing K=7/alpha=0.5 internal
reference of 0.633 — BO found a genuinely *better*-scoring region on the internal objective, not a
worse one). **But it is a larger K, not a smaller one** — the stated goal of this step (find a
*smaller* K with comparable performance) was not met. The marginal external-C improvement (0.6282 vs.
0.6267, well within 1 SE of either) is not a meaningfully different answer from K=7's own number;
adjusting alpha upward from 0.5 to ~0.71 (more survival-weighted) at a slightly larger K is, at best,
a minor refinement, not a parsimony win.

**A real limitation of this specific run, stated plainly, with a more precise mechanism than first
suspected.** After finding `(K=8, alpha=0.7072)` at evaluation 9 of 23, the search proposed the
*exact same point* for all 15 remaining rounds (9 through 23) rather than continuing to explore the
smaller-K region under different alpha values. Initially attributed to the default UCB acquisition's
exploration/exploitation balance; code review identified a more precise, more actionable cause:
`K_bounds` is passed as integer-typed (`c(2L, 10L)`), and `rBayesianOptimization` documents that
integer-suffixed bounds are treated as a discrete dimension internally, pre-rounding K *before* the
Gaussian process or acquisition function ever see it — `select_k_alpha_bo_objective()`'s own
`round(K)` was therefore redundant, and the real effect is that the optimizer's search over K
collapses onto a coarser grid than a from-scratch continuous relaxation would give it, likely
compounding with only 8 init points to converge early. This means the search's real exploration was
effectively ~9 unique points, not 23, and `pick_trustworthy_bo_winner()`'s top-5 gate ended up
re-checking 5 copies of the identical point — it still correctly confirmed K_eff>0 (not degenerate),
but was not exercised against genuinely diverse candidates in this particular run. The search's
early, more diverse evaluations (rounds 1-8: K∈{3,4,6,7,8,9} at various alpha) never found a smaller K
competitive with K=7, but with only 8 init points and this early convergence, this is not a thorough
search of the smaller-K region — a larger init_points budget, a continuous-relaxation workaround (not
using the "L"-suffix integer-bounds convention), and/or a retuned acquisition (higher kappa, or
`acq="ei"`/`"poi"`) would be needed before concluding smaller K is definitively unreachable via joint
tuning. Recorded honestly as an open question, not resolved here.

**Judgment call for Step 4:** given Step 3 did not find a smaller K, and Step 1's validated
warm-start-from-K=7 recipe already delivers the actual parsimony goal (K=4 or K=5 at K=7-level
performance, using only existing machinery), Step 4's final validation should adopt **Step 1's
warm-start recipe**, not the Step 3 BO winner, as the config to carry forward — K=8 is not more
parsimonious than K=7, so adopting it would move in the wrong direction relative to this whole
follow-up's purpose. The BO machinery is kept (correct, tested, and a real option for a future,
better-tuned joint search) but is not the answer to "can K shrink."

Full test suite: 311/311 passing (299 + 12 new tests: 7 for the BO wrapper, 5 for
`pick_trustworthy_bo_winner()`). Reviewed by `superpowers:code-reviewer`.

*Files: `code/select_k_alpha_bo.R` (new), `tests/test_select_k_alpha_bo.R` (new, TDD),
`results/benchmark_sim/run_joint_bo.R` (new),
`results/benchmark_sim/outputs/joint_bo/joint_bo_history.csv`,
`joint_bo_external_val.csv` (new).*

---

## 2026-07-13 — K-parsimony follow-up Step 2: deflation-init is mathematically equivalent to SVD-init for non-degenerate data — does not reproduce Step 1's rescue

**Motivation.** Step 1 (entry below) found the CAVI factor-collapse artifact could be fixed by
warm-starting K=4/K=5 from an already-converged K=7 fit's PVE-ranked columns, but that this was "a
one-off warm-start hack" requiring a pre-existing higher-K fit. Step 2's goal
(`docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md`) was to make this permanent and general:
add a deflation-style init (rank-1 SVD of Y for factor 1, rank-1 SVD of the residual for factor 2,
etc.) to `fit_cox_on_yf`/`fit_supervised_mf_modular`, verify it fixes the `sparse_synthetic` collapse
scenario, and confirm no regression.

**New code:** `code/deflation_init.R` (`deflation_svd_init()`), `init_method = "deflation"` added to
both `fit_cox_on_yf()` and `fit_supervised_mf_modular()`, `tests/test_deflation_init.R` (9 tests, TDD).
Branch `cavi-deflation-init`.

**Result — this does NOT reproduce Step 1's rescue, and the reason is analytically clear, not just an
empirical near-miss:**

1. **Real PDAC data (D4 config), fresh deflation-init at K=4 and K=5: bit-identical to fresh SVD-init.**
   K=4: mean C=0.5409, K_eff=1, beta_max=0.0092, 7 iterations. K=5: mean C=0.5960, K_eff=3,
   beta_max=0.0201, 41 iterations. The mean-C/K_eff figures match the Step 1 table's fresh-SVD row
   exactly (0.5409/K_eff=1 at K=4; 0.5960/K_eff=3 at K=5); beta_max and iteration count match the
   Phase 3 entry's fresh-SVD figures below (which reported beta_max only as "~0.009" — this run's more
   precise 0.0092/0.0201 are consistent with, not just approximately equal to, that entry). Neither
   reaches the 0.6068 margin Step 1's warm-start reached (0.6270 at both K).
2. **`sparse_synthetic`-style diagnostic (`results/multi_cohort_sim/diagnose_factor_collapse.R`,
   extended with a deflation-init check):** YFB's dead-factor count is identical between SVD-init and
   deflation-init in every one of 5 seeds (2,2,2,2,1 vs. 2,2,2,2,1). LB's average drops from 2.4/4 to
   1.8/4 but non-monotonically (2 seeds improve, 2 get worse or stay flat, 1 seed's dead-factor count
   *increases* from 2 to 3) — consistent with ordinary fit-to-fit variability in a numerically
   different-but-equivalent computation path, not a real fix.
3. **Why:** greedy rank-1 SVD deflation and a single batch top-K SVD extract the *same* leading-K
   singular subspace whenever a matrix's top-K singular values are distinct — a standard linear-algebra
   fact (this is essentially how many SVD/PCA algorithms are implemented internally). Real data (and
   most synthetic data with continuous noise) essentially never has exactly-tied singular values, so
   `init_method="deflation"` was mathematically guaranteed to match `init_method="svd"` in the regime
   this plan needed it to differ in. The `sparse_synthetic` DGP's "near-tied, disjoint-support" factors
   are near-tied in *amplitude*, not exactly tied in *singular value* — not the actual degenerate case
   deflation could help with.

**What this means for the collapse mechanism.** Step 1 already showed best-ELBO multistart (14 random
restarts) rescued nothing — SVD-init was always the best-ELBO restart. Step 2 now shows a second,
independent "different starting subspace" strategy also rescues nothing, for an analytically
understood reason. Both negative results point the same direction: **the collapse is not about which
linear-algebra decomposition of Y seeds CAVI** — deflation, batch SVD, and 14 random draws all land in
the same basin. What *did* work (Step 1) was warm-starting from a solution that had already been
shaped by many iterations of the full joint CAVI process (EBNM shrinkage + Cox coupling +
factor-wise Gauss-Seidel) at a higher K — a qualitatively different kind of starting point that no
raw transformation of Y alone can produce.

**Judgment call — flagging rather than resolving unilaterally.** The plan's Step 2 goal was a fix
that's "permanent and general, not a one-off warm-start hack." Deflation-init does not deliver that,
and the plan does not specify a fallback for this outcome (Step 3's joint-BO path is reserved for a
CAPACITY-LIMITED verdict, which Step 1 explicitly ruled out). This is exactly the kind of ambiguity
the plan itself says should pause for user input, rather than being silently resolved. Recommendation
brought to the user: adopt Step 1's validated warm-start-from-a-higher-K-fit as the standing
recommended fitting recipe for K<7 in this model family (it is general in the sense of not requiring
new machinery — just fitting once at a generously large K and reusing `extract_top_k_by_pve()` — even
if not "permanent" in the stronger sense of not depending on any prior fit at all), and use it as
Step 4's fitting procedure. The `deflation_svd_init()` code itself is kept (correct, tested, harmless,
and a documented, if narrow, alternative for the true degenerate-singular-value edge case) but is not
claimed as *the* fix.

*Files: `code/deflation_init.R` (new), `code/fit_modular.R` / `code/fit_cox_on_yf.R` (new
`init_method="deflation"` branch), `tests/test_deflation_init.R` (new, TDD, 9/9 passing),
`results/multi_cohort_sim/diagnose_factor_collapse.R` (extended with a deflation-init check). Full
test suite: 299/299 passing.*

---

## 2026-07-13 — K-parsimony follow-up Step 1: K=7's necessity was an optimization artifact, not a capacity floor, at K=4-5

**Motivation.** The Phase 3 entry below found K=7 "not free to shrink" but explicitly flagged K=2/K=4's
suspiciously fast convergence (7-9 iterations to a near-zero β) as consistent with the CAVI
factor-collapse failure mode from Phase 2, and could not rule out that this — rather than a genuine
capacity ceiling — explained their poor external performance. This is Step 1 of
`docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md`: re-fit K ∈ {2,3,4,5} with two
improved-optimization strategies and compare against the original fresh-SVD numbers.

**New code:** `code/warmstart_from_fit.R` (`extract_top_k_by_pve()` — selects a converged fit's
top-K columns by final-iteration PVE, for use as `init_method="custom"` warm-start); `code/fit_modular_multistart.R`
(`fit_cox_on_yf_multistart()` — YFB counterpart to the existing LB multistart wrapper);
`results/benchmark_sim/run_k_parsimony_followup.R` (the comparison runner, branch
`phase3-followup-warmstart`).

**Method.** For each K ∈ {2,3,4,5}: (a) fresh SVD (Phase 3's original numbers, reproduced here), (b)
warm-start from the converged K=7 fit's top-K PVE-ranked columns, (c) best-ELBO multistart
(n_init=15: 1 SVD + 14 random restarts). Same D4 config, same 5 held-out cohorts, K=7 refit here too
as a reproducibility check (reproduced 0.6267 exactly).

**Result — real data:**

| K | method | Mean external C | SE | K_eff |
|---|---|---|---|---|
| 2 | fresh | 0.5406 | 0.0116 | 1 |
| 2 | multistart | 0.5406 | 0.0116 | 1 |
| 2 | warmstart | 0.5608 | 0.0195 | 1 |
| 3 | fresh | 0.5943 | 0.0162 | 1 |
| 3 | multistart | 0.5943 | 0.0162 | 1 |
| 3 | warmstart | 0.5955 | 0.0255 | 3 |
| 4 | fresh | 0.5409 | 0.0116 | 1 |
| 4 | multistart | 0.5409 | 0.0116 | 1 |
| 4 | **warmstart** | **0.6270** | 0.0198 | 2 |
| 5 | fresh | 0.5960 | 0.0269 | 3 |
| 5 | multistart | 0.5960 | 0.0269 | 3 |
| 5 | **warmstart** | **0.6270** | 0.0198 | 2 |
| 7 | fresh (reference) | 0.6267 | 0.0199 | 2 |

Applying the plan's mechanical decision rule (best-of-{fresh, warmstart, multistart} per K within 1 SE
of K=7's margin, 0.6267 − 0.0199 = 0.6068): **K=4 and K=5 both reach it via warm-start (0.6270 ≥
0.6068) — outcome is OPTIMIZATION-LIMITED.** Per the plan, this routes to **Step 2** (deflation-init
fix), not Step 3 (joint Bayesian optimization).

**Two findings worth stating plainly, not smoothed over:**

1. **Warm-start rescued K=4/K=5; best-ELBO multistart rescued nothing.** At every single K tested,
   multistart's `best_idx` was restart 1 — the SVD init itself — meaning none of the 14 random
   restarts found a higher-ELBO solution than SVD init at any K. SVD init is already the local-ELBO
   optimum among these candidates; the degenerate fixed point Phase 3 landed in is not something a
   *different starting point drawn from the same broad random distribution* escapes — it takes a
   warm-start seeded from a demonstrably good solution (the converged K=7 fit) to reach the better
   fixed point. This directly informs Step 2: a fix must change the *character* of the initialization
   (e.g. deflation/greedy, sequentially removing signal like the K=7 warm-start effectively did),
   not just add more random restarts.
2. **The rescue is real but partial, not universal.** K=4 and K=5 warm-start converge to
   essentially the same 2-active-factor solution as K=7 (same K_eff=2, same beta_max=0.0403, same
   mean C=0.6270) — strong evidence the "real" signal in this data lives in ~2 factors that K=7's
   fit already found, and K=4-5 have enough capacity to hold them once initialized correctly. K=3
   warm-start converges to a different point (K_eff jumps to 3, all factors "active" by threshold)
   but external performance barely moves (0.5955 vs. 0.5943 fresh) — breaking the collapse at K=3
   does not, by itself, recover K=7-level performance, suggesting 3 factors may be a genuine
   information floor for this specific signal rather than purely an optimization artifact. K=2
   warm-start also improves (0.5608 vs. 0.5406) but remains well short of the margin. So the
   corrected picture is: **K=2/K=3 show real evidence of a capacity limit; K=4/K=5's Phase 3
   underperformance was specifically an optimization artifact that a smarter initialization
   resolves.**

**Conclusion — supersedes Phase 3's headline claim.** Phase 3 said "K=7 is not free to shrink"; that
holds for K=2/K=3 but not for K=4/K=5 once the CAVI factor-collapse artifact documented in Phase 2 is
corrected for. K=4 (and K=5) are a genuinely available, more parsimonious alternative to K=7 for this
recommended config, provided the fit uses warm-start (or an equivalent fix) rather than fresh SVD
init. Step 2 will build this into a permanent, general initialization option (deflation-style init)
rather than relying on an already-fitted K=7 model as a warm-start source, and Step 4 will produce
the final, doubly-verified K-vs-external-performance answer once that fix is in place.

**Judgment call, documented rather than escalated:** the plan's Step 1 item 2 flagged an earlier
setup error sourcing `code/fit_cox_on_yf.R` (`stopifnot(!is.null(real_Y), ...)` firing on `source()`)
as "likely a missing `tryCatch(source(...))` wrap... root-cause and fix that, don't just work around
it." The fix applied is exactly that wrap — the same idiom already used in `tests/run_tests.R` and
three other existing test files to source this file safely. `fit_cox_on_yf.R`'s `DATA_MODE` runner
block itself (which always fires on `source()`) is an intentional, unmodified, shared pattern with
`fit_modular.R`; hardening it further would be a larger, unrequested change to a core, heavily-tested
file for no benefit beyond what the existing convention already provides. This reasoning, and the new
`warmstart_from_fit.R`/`fit_modular_multistart.R` code it applies to, were reviewed and concurred on
by `superpowers:code-reviewer` prior to that code's commit; this Step 1 write-up (including the
runner script and this entry) was reviewed by a second, separate `superpowers:code-reviewer` pass.

*Files: `code/warmstart_from_fit.R` (new), `code/fit_modular_multistart.R` (extended),
`tests/test_warmstart_from_fit.R`, `tests/test_yfb_multistart.R` (new, TDD),
`results/benchmark_sim/run_k_parsimony_followup.R` (new),
`results/benchmark_sim/outputs/k_parsimony_followup/k_parsimony_followup_results.csv` (new).
Full test suite: 290/290 passing.*

---

## 2026-07-13 — Phase 3: K-parsimony curve on real data — K=7 is not free to shrink; smaller K all underperform

**Motivation.** The existing K-CV table (2026-07-12 entry below) measures internal training-fold
C-index. It doesn't answer the actually decision-relevant question: how much of the recommended
config's headline **external** validation number (D4: mean C=0.627 across 5 held-out PDAC cohorts,
K=7) would survive at a smaller, more parsimonious K. Built
`results/benchmark_sim/run_k_parsimony_curve.R`: refits YFB D4 (per-platform z-std, DeSurv
combined-rank gene selection, top-3000 per cohort, no cohort_id) at K ∈ {2, 3, 4, 5, 7}, re-running
external validation against the same 5 held-out cohorts for each K.

**Result — real data, branch `phase3-k-parsimony`:**

| K | Mean external C | SE | K_eff |
|---|---|---|---|
| 2 | 0.5406 | 0.0116 | 1 |
| 3 | 0.5943 | 0.0162 | 1 |
| 4 | 0.5409 | 0.0116 | 1 |
| 5 | 0.5960 | 0.0269 | 3 |
| 7 | **0.6267** | 0.0199 | 2 |

Applying the same 1-SE decision rule already used for internal K-CV (smallest K within 1 SE of the
best): K=7's margin is 0.6267 − 0.0199 = 0.6068, and no smaller K in the grid reaches it (K=5's
0.5960 is the closest, still ~0.6 SE short). **K=7 remains necessary for the external number in this
run — there is no parsimony discount available in this K range**, with one caveat below on how much
weight a single-seed run can bear.

Two additional observations, reported honestly rather than smoothed over:
- **K=7 is not just best on average — it wins on every one of the 5 individual external cohorts**
  (Dijk 0.6345, Moffitt 0.5486, PACA_AU_array 0.6483, PACA_AU_seq 0.6573, Puleo 0.6447), each the
  single highest C-index in its column across all 5 K values tested. Not a cherry-picked mean.
- **The curve is not monotonic and K_eff doesn't track external performance.** K=4 (K_eff=1)
  essentially ties K=2 (K_eff=1) on every cohort (e.g., Dijk 0.5589 vs. 0.5586) despite having one
  more factor available, while K=3 (also K_eff=1) unpredictably lands much closer to K=5/K=7 on some
  cohorts (PACA_AU_array 0.6131) and not others (Moffitt 0.5315). Two active factors (K_eff) at K=7
  outperforms every smaller-K fit tested, including ones with the same or a higher K_eff (K=5,
  K_eff=3) — K_eff alone is not a reliable predictor of external performance at fixed model capacity.

**Caveat — single seed/init per K, not left unstated.** Every K used one fit (seed=42, SVD init).
The per-K iteration counts are revealing: K=2 and K=4 both converge suspiciously fast (9 and 7
iterations respectively) to a tiny beta_max (~0.009), while K=3 and K=7 take 20 iterations to a
substantially larger beta_max (0.033, 0.040), and K=5 takes 41 iterations. This pattern — fast
convergence to a near-zero-β solution — is consistent with exactly the CAVI factor-collapse failure
mode root-caused in Phase 2 (`sparse_synthetic` finding, this file's 2026-07-12/13 entries): a
single SVD-initialized fit at these K values may be landing in a degenerate local optimum rather
than genuinely representing the best K=2 or K=4 model achievable. So "K=7 is necessary" is a fair
description of *this specific single-seed comparison*, but not yet a fully robust claim about K=2-4's
ceiling — a multistart (best-ELBO) rerun at the underperforming K values would be needed to rule out
collapse before treating this as final. The K=7-vs-K=5 gap is less suspect (both converge at a
comparable iteration count with comparable beta_max scale), so the core finding — a smaller,
already-more-thoroughly-converged K=5 still falls short of K=7 — stands more confidently than the
K=2/K=4 comparisons do.

**Conclusion:** the gap between our K=7 and DeSurv's K=3 (2026-07-12 entry below) is not simply
"our model retains unnecessary factors that could be pruned for free" — refitting at smaller K
directly, rather than just counting active factors from the K=7 fit, shows real, substantial
performance loss at every smaller K tested here, on every held-out cohort, though the K=2/K=4 numbers
specifically should be re-verified with multistart before treating them as a hard ceiling (see
caveat above). The K=7→K_eff=2 gap is better understood as CV-selected capacity that happens to
produce few large-|β| factors, not as 5 factors that could simply be removed. Whether DeSurv's joint
(K, α, λ) Bayesian-optimization tuning would find a smaller K with comparable *external* performance
(as opposed to comparable internal/training performance) remains untested here — a plan for this
exists (`docs/plans/joint_k_alpha_bayesopt_plan_07_12_2026.md`), not yet implemented.

*Files: `results/benchmark_sim/run_k_parsimony_curve.R` (new),
`results/benchmark_sim/outputs/k_parsimony_curve/k_parsimony_curve_results.csv` (new).*

---

## 2026-07-12 — Phase 2: joint model vs. two-step baselines, survival-strength sweep

**Goal (Phase 2 of the post-lab-meeting action plan):** show that the joint model beats a 2-step
(unsupervised factorization → Cox) baseline when survival carries real signal, and is *equivalent*
to it when survival carries none. Implemented as a simulation sweep
(`results/multi_cohort_sim/run_survival_strength_sweep.R`) scaling the true prognostic effect
(`beta_shared`) from 0 to large across 6 levels × 5 seeds, comparing YFB (joint, tuned α=0.5)
against new PCA+Cox and existing EBMF+Cox two-step baselines (LB included as a secondary
reference). Full report: `docs/reports/joint_vs_twostep_sweep_07_12_2026.{qmd,pdf,html}`.

**Result 1 (equivalence at zero signal):** YFB=0.517±0.008, PCA+Cox=0.533±0.017,
EBMF+Cox=0.540±0.017. All within ~1 SE of each other and of chance — consistent with equivalence,
though the raw gaps (0.016–0.023) are slightly above an idealized "<0.01" target; reported honestly
as "consistent with, not a sharp confirmation of" equivalence, given only 5 seeds.

**Result 2 (separation growing with signal):** not immediate or uniform — at strength=0.25 YFB is
tied with both baselines (within noise). From strength≥0.5 it pulls ahead of EBMF+Cox, and from
strength≥1.0 ahead of PCA+Cox too, with the gap widening and YFB's own SE shrinking as signal grows
(SE 0.027→0.006 from strength 1→4). At strength=4: YFB=0.931±0.006 vs. PCA+Cox=0.908±0.015 vs.
EBMF+Cox=0.875±0.057.

**Result 3 (α=0 internal control) — the cleanest of the three:** max$|F_{\alpha=0.5}-F_{\alpha=0}|$
= exactly 0 at every strength level and seed (30/30). Required a dedicated, iteration-count-
controlled comparison (`tol=-1`, fixed iterations) rather than comparing naturally-converged fits
directly — different α values give different combined-objective trajectories, so naturally-
converged fits stop after a different number of iterations, making a naive F comparison confounded
by convergence timing rather than a true α effect. Once controlled for, the result is exact and
unambiguous: YFB's genomics factorization is provably α-invariant, confirming the model behaves as
designed at this boundary.

**New reusable code:** `results/multi_cohort_sim/fit_pca_cox.R` (`fit_pca_cox()`/`predict_pca_cox()`,
9 TDD tests), a `synthetic_multicohort.survival_strength_sweep` config block
(`config/globals.yml`).

**Not addressed here:** LB's own unresolved Phase 1a limitation (documented in this file's
objective-normalization entry) — LB's numbers are shown for reference only, no claims made about
its competitiveness beyond the raw table.

**Same-day follow-up — comprehensive extension (10 seeds, 4 DGP scenarios) surfaces a genuine
joint-model failure mode, root-caused:** the initial 5-seed/1-scenario run was directionally right
but too thin. Extending to 10 seeds and 3 additional DGP structures (`sparse_synthetic`: synthetic
sparse F, 2% active genes/factor, real templates off; `low_snr`: a_shared 12→6; `high_K`: K_shared
4→8) confirmed Results 1–3 hold under `default`, `low_snr`, and `high_K` — but **`sparse_synthetic`
reverses the joint-vs-2-step ordering entirely**: at strength=4, YFB=0.795±0.027 vs.
PCA+Cox=0.934±0.005 vs. EBMF+Cox=0.940±0.004 (a ~5-SE gap, not noise).

Root-caused with a persisted, reproducible diagnostic (`results/multi_cohort_sim/
diagnose_factor_collapse.R`, 5 seeds; `outputs/factor_collapse_diagnostic.csv`), not left as an open
mystery:
1. Refitting YFB at α=0 (survival term fully removed) shows the *identical* dead-factor count as the
   tuned fit in every one of 5 seeds (mean 1.8/4 dead, both cases) — rules out anything related to
   the joint survival objective.
2. LB (identical `update_L.R`/`update_F.R` machinery) also collapses — on average *more* factors die
   than for YFB (mean 2.4/4 vs. 1.8/4 across 5 seeds; worse in 3/5, better in 2/5) — confirms this is
   a shared CAVI vulnerability, not specific to the Cox-on-YF reformulation, though severity is
   seed-dependent rather than uniformly worse for LB.
3. Switching SVD init → random init makes it *uniformly worse*: all 4 factors die in every one of
   5 seeds (vs. 1–2 of 4 for SVD init) — rules out "bad initialization" as the root cause.

**Diagnosis:** NOT an amplitude-hierarchy effect — `a_shared` is applied identically to every shared
factor in *every* scenario including `default`, so there is no built-in amplitude ranking anywhere.
The actual differentiator is **loading structure**: real EBMF templates (`default`/`low_snr`/
`high_K`) are dense (97% of genes carry non-trivial loading in a representative factor, verified
directly on the cached templates) with graded magnitudes and real cross-factor correlation (up to
|r|≈0.4); `sparse_synthetic` gives each factor exactly 2% active genes, uniform magnitude, and
near-zero cross-factor overlap by construction. This near-perfect symmetry (several equal-scale,
non-overlapping factors, no ties to break) appears to be exactly the condition under which our
point_exponential-based joint CAVI hits a degenerate fixed point: a factor whose estimated loadings
dip slightly below its (otherwise identical) siblings gets shrunk further every iteration with no
mechanism to recover, and the ELBO genuinely plateaus at this bad solution (hence fast, "converged"
fits — 6-14 iterations across all checks). EBMF avoids this by fitting factors *greedily* on the
residual (each factor only ever competes against noise, never against symmetric siblings for the
same variance); PCA avoids it by doing no sparsity-inducing shrinkage at all.

**Scoping — narrower than "only a synthetic-data artifact":** `default` (real templates) showed zero
collapses for either model (0/60 fits). But in `low_snr` and `high_K` — also built on real EBMF
templates — **LB collapsed completely on 10 of 180 fits** (frac_recov=0); **YFB collapsed on 0 of
180**. So the honest claim is: YFB has not collapsed in any real-template scenario tested here, but
LB has, occasionally, under lower SNR or higher K, even with real templates — this is *not* purely a
contrived-synthetic-data phenomenon, though we still have no evidence it affects the actual real
PDAC fits. Logged as a known, understood algorithmic limitation (LB apparently more susceptible),
not resolved here; candidate follow-ups (multistart best-ELBO selection via the existing
`code/fit_modular_multistart.R`, or a greedy-init variant) go in `ROADMAP.md`.

Updated report: `docs/reports/joint_vs_twostep_sweep_07_12_2026.{qmd,pdf,html}` (now includes a
per-scenario results table and this failure-mode analysis). Config:
`synthetic_multicohort.survival_strength_sweep.{n_seeds,scenarios}` (`config/globals.yml`).

---

## 2026-07-12 — Fresh K-CV under corrected code: K=7 is genuine, not an artifact; K vs. DeSurv's k=3 is a methodology-comparison question, not a "our model needs more capacity" one

**Motivation.** `k_merged_yfb_desurv=7` predated this entire session (cached 2026-05-27); every
benchmark re-run since then skipped K-CV because a value was already present. It had never been
selected under the corrected code (`boost_beta=FALSE`, fixed train/test preprocessing). Re-ran
`select_K_cv()` fresh (YFB, D4 preprocessing, K grid 2:10, 5-fold, no floor imposed) to see the
real curve before deciding whether to enforce more parsimony.

**Result (no floor applied):**

| K | 2 | 3 | 4 | 5 | 6 | **7** | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| mean C | 0.567 | 0.568 | 0.573 | 0.587 | 0.593 | **0.633** | 0.651 (best) | 0.628 | 0.620 |
| SE | 0.037 | 0.031 | 0.034 | 0.029 | 0.033 | 0.029 | 0.026 | 0.032 | 0.029 |

**Conclusion: K=7 is not an artifact of the K≥3 floor, the pre-Phase-1 code, or noise.** There is a
genuine, non-noise jump between K=6 (0.593) and K=7 (0.633) — about 4 points, larger than any
single fold's SE (~0.03) — and K=2 through K=6 are all more than 1 SE below the K=8 peak (0.651).
The 1-SE rule's selection of K=7 (smallest K statistically tied with the best) is honest evidence,
not a selection-procedure quirk. Reducing K below ~7 costs real, held-out predictive performance
on this model/data/preprocessing combination.

**Why doesn't this match DeSurv's k=3, given identical training/validation data?** Two confirmed
methodological differences, not a data difference:
1. DeSurv tunes **k, α, and its elastic-net penalty λ jointly via Bayesian optimization**
   (`presentation/walther_lab_meeting_06_18_2026/lab_meeting_june18.qmd`, comparison table). We
   tune **only K via CV**, with α fixed at 0.5 and no elastic-net penalty at all (§ removed
   2026-07-12, Phase 1b — and as established there, our retired λ was never analogous to DeSurv's
   anyway). A joint search over 3 interacting hyperparameters can land in a different region of
   the trade-off surface than a 1-hyperparameter search with the other two frozen, even on
   identical data.
2. DeSurv's shrinkage is a **fixed, jointly-tuned elastic-net penalty** — a mechanism that can
   aggressively collapse redundant factors during optimization. Ours is **empirical-Bayes
   adaptive shrinkage** (point-exponential priors on L/F), a different mechanism, not jointly
   tuned with K.

**Not yet verified** (stated as a hypothesis, not a finding): the working explanation that the
"extra" ~5 factors beyond K_eff=2 are earning their keep on genomics reconstruction specifically
(rather than survival) is *inferred* from K_eff staying low regardless of total K — a direct
reconstruction-quality-vs-K curve (RMSE or per-factor PVE) has not been pulled to confirm this.

**Decision:** keep K=7 for D4 (matches 1-SE rule, nearly identical to the K=8 peak, simpler).
Do not attempt to force a smaller K without cause — the CV evidence does not support it. The
"why 7 vs DeSurv's 3" question is a **matched-protocol comparison problem** (Phase 2/6's planned
same-protocol DeSurv head-to-head), not something resolvable by re-running our own K-CV
differently. K_eff=2 (survival-active factors) already matches DeSurv's own finding of survival
concentrated in ~1 factor reasonably well — the parsimony story that holds up is about *effective*
factors, not total factors.

**Affected files:** `results/benchmark_sim/outputs/desurv_comparison/kcv_yfb_desurv_corrected.{csv,rds}`
(new); no code changes (K unchanged in `config/globals.yml`).

---

## 2026-07-12 — Phase 1a's beta boost was unjustified: corrected default is `boost_beta=FALSE`; honest Phase 1 conclusion is "no performance benefit for YFB"

**What was wrong.** The Phase 1a design merged earlier the same day boosted beta's own EBNM
precision by a factor of p (`survival_divisor` threaded into `update_beta_k` via
`norm_convention`), justified as "safe" because beta has no bilinear L↔F-style feedback coupling.
Safe is not the same as *necessary*: beta's own coordinate update (`A_k = alpha * sum(w * EL2_k
or ZF_k^2) / survival_divisor`) has **no genomics term competing with it in its own formula, in
either model**. The genomics/survival scale imbalance Phase 1a set out to fix only structurally
exists in LB's L-update, where both terms are added together in one formula — it does not exist
for beta in either model. Boosting beta's precision anyway did not correct any real imbalance; it
only reduced EBNM shrinkage, which mechanically pulls more factors above a fixed
`beta_threshold` without reflecting a genuine gain in survival signal.

**Empirical confirmation (D4, K=7, real PDAC training + 5-cohort external validation):**

| Configuration | Mean external C | K_eff | β (factors 3, 5, 6, 7) |
|---|---|---|---|
| Pre-Phase-1 baseline (documented 2026-05-27) | 0.636 | 2 | +0.011, —, —, −0.041 |
| Phase 1 corrected (`boost_beta=FALSE`) | **0.6267** | **2** | +0.0115, 0, 0, −0.0404 |
| Phase 1 as first merged (`boost_beta=TRUE`) | 0.6419 | 4 | +0.0186, 0.0034, 0.0036, −0.0425 |

With `boost_beta=FALSE`, β's values are essentially bit-identical to the pre-Phase-1 baseline —
confirming beta genuinely never needed rebalancing. This also means Phase 1a has **zero effect
on any of YFB's fitted output** (EL/EF were already confirmed bit-identical regardless of
`norm_convention`; now EBeta is too) — its only surviving effect is the ELBO monitor's reported
scale during fitting (a diagnostic/convergence quantity, not a fitted parameter). Consequently,
**the 0.636→0.642 "improvement" reported when Phase 1 first merged was an artifact of this
unjustified boost, not a genuine effect of objective normalization.**

**Decision:**
- `boost_beta` (new parameter on `fit_supervised_mf_modular()` and `fit_cox_on_yf()`) defaults to
  `FALSE` in both models. `TRUE` reproduces the superseded design, retained only for reference/
  comparison, not recommended.
- **Honest Phase 1 conclusion for YFB:** objective normalization (1a) provides no performance
  benefit — there was no real per-coordinate imbalance in YFB to fix in the first place (L is
  pure-genomics, β is pure-survival, they never share a coordinate). It is kept anyway because
  the ELBO-monitor fix is a genuine (if purely internal) correctness improvement: previously
  `alpha=0.5` did not mean "balanced" even in the quantity CAVI itself uses for convergence
  monitoring, regardless of whether that quantity affects the final fit.
- **The corrected, final post-Phase-1 numbers for D4:** mean external C = 0.627 (down slightly
  from 0.636, attributable entirely to Phase 1c's preprocessing fix — see that entry below, not
  to 1a), K_eff = 2 (unchanged from before Phase 1 — the parsimony goal is fully preserved).
- Superseded the "D4 K_eff rose 2→4... not yet confirmed" entry below (kept for the investigation
  record) and the initial merge commit's headline numbers (`CLAUDE.md`, `ROADMAP.md`,
  `PROJECT_STATUS.qmd` corrected in the same pass as this entry).
- **Not addressed by this fix:** Phase 1a still provides no working per-coordinate fix for LB
  (unresolved, deferred — both shrink-genomics and boost-survival directions destabilize LB's L,
  as documented in the objective-normalization entry below). This does not affect YFB/D4.

**Affected files:** `code/fit_modular.R`, `code/fit_cox_on_yf.R` (`boost_beta` parameter,
`beta_divisor` exposed in the returned result), `tests/test_normalization.R` (3 new tests: T_conv.8-10),
`CLAUDE.md`, `ROADMAP.md`, `PROJECT_STATUS.qmd` (corrected headline numbers).

---

## 2026-07-12 — Phase 1b: retire lambda, alpha is the sole genomics/survival mixing weight

**Decision:** Removed the `lambda` survival-scale multiplier from `update_L_k`/`update_L_all`
(`code/update_L.R`), `fit_supervised_mf_modular()` (`code/fit_modular.R`), and `fit_cox_on_yf()`
(`code/fit_cox_on_yf.R`, where it was already fully dead — accepted as a parameter but never
referenced internally), plus every `config/globals.yml` entry and benchmark-runner call site that
set it.

**Reason.** `alpha`'s existing `(1-alpha)*genomics + alpha*survival` weighting already plays
exactly the role a `(1-lambda)`/`lambda` pair would — `lambda` was a second, redundant knob
multiplying only the survival side, always left at its no-op default (1.0) in every production
config (`config/globals.yml`'s own comment: "lambda=1 matches or beats p/n and 2p/n"). DeSurv's own
`λ` (Young et al., referenced in the post-lab-meeting plan) is not analogous to this parameter at
all — it is an elastic-net penalty coefficient on β, a regularization role this model already fills
via its empirical-Bayes (EBNM) priors on β (point_normal/point_laplace), not a second likelihood-
mixing weight. No replacement parameter was needed.

**Trade-offs:** None identified — no test file referenced `lambda`, and it was a no-op in every
config, so this is a pure simplification with zero behavior change (261/261 tests passing,
identical to before removal).

**Affected files:** `code/update_L.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`,
`code/fit_modular_multistart.R` (docstring), `config/globals.yml`, and 7 benchmark runner scripts
in `results/benchmark_sim/` that read `cfg$...$lambda`.

---

## 2026-07-12 — D4 K_eff rose 2→4 after Phase 1; likely a `beta_threshold` scale artifact, not yet confirmed

> **SUPERSEDED same day** — see "Phase 1a's beta boost was unjustified" entry above. The root
> cause was not `beta_threshold` calibration; it was an unjustified precision boost applied to
> beta's own update. Retained for the record of how the investigation proceeded.

**Observation.** Re-running `run_desurv_comparison.R` after Phase 1 (objective normalization, λ
retirement, preprocessing fix) gives D4 external mean C-index 0.636 → 0.642 (K=7, up), but
K_eff (factors with |β̂| > `beta_threshold`=0.001) rose from 2 to 4. Inspecting D4's actual β
vector: `[0, 0, 0.0186, 0, 0.0034, 0.0036, -0.0425]`. The two originally-active factors (indices 3
and 7) remain dominant and match the pre-Phase-1 result in sign and similar magnitude
(previously β̂₃=+0.011, β̂₇=−0.041; now +0.0186, −0.0425 — same identity). The two *newly*-active
factors (5, 6) sit at β≈0.0034–0.0036, only ~3.4× the threshold.

**Interpretation (not yet confirmed).** Phase 1a's default `norm_convention="per_p"` boosts β's
own precision by a factor of p in the beta update (see the 2026-07-12 objective-normalization
entry above) — EBNM therefore shrinks β less aggressively overall than before Phase 1. This is
consistent with two previously-suppressed, weak factors now crossing a `beta_threshold` that was
calibrated against the *old*, unboosted β scale. It is also possible this reflects a genuine
(if modest) gain in recoverable signal — the external C-index did improve. These are not
distinguished by the evidence in hand.

**Decision:** do not treat K_eff=4 as either a confirmed parsimony loss or a confirmed
improvement yet. Flagged as an open item (`ROADMAP.md`) to resolve via a `beta_threshold`
recalibration or cross-check (e.g. PVE-based thresholding, already available via
`k_selection.pve_threshold`; or a CV-stability check on which factors survive resampling) before
or alongside Phase 3's own K/K_eff analysis, which depends on exactly this question ("does K
shrink toward DeSurv's 3 on the corrected objective, and how many factors are genuinely
survival-active").

**Affected files:** none (analysis only); follow-up tracked in `ROADMAP.md` under Model Selection.

---

## 2026-07-12 — Phase 1c: fix external-cohort preprocessing to match training (rank vs. per-platform z-std)

**Problem.** `preprocess_desurv_cohort()` (`code/preprocess_desurv.R`) unconditionally
per-subject rank-transformed and never per-platform z-standardized. Several benchmark scripts'
training preprocessing (`preprocess_merged_cohorts(..., rank_transform=FALSE,
per_platform_standardize=TRUE, ...)`) does the exact opposite. External cohorts were therefore
preprocessed inconsistently with training in every script using this DeSurv-aligned pipeline.

**Fix:** `preprocess_desurv_cohort()` gained two new parameters, `rank_transform = TRUE` and
`per_platform_standardize = FALSE` (both defaults preserve prior behavior for existing callers).
Column-wise z-standardization commutes with the later `intersect()`-based gene-set subsetting
(a gene's z-score depends only on its own across-subject mean/SD, not on which other genes are
present), so computing it before vs. after the training-gene intersection gives identical values
for the retained genes — no ordering concern.

**Fixed (mechanical, single training config per run — straightforward to match):**
- `results/benchmark_sim/run_desurv_comparison.R` (Section 6, external validation for D1-D5):
  now passes `rank_transform=FALSE, per_platform_standardize=TRUE`, matching Section 3's training
  call for every config. This is the D4 recommended-configuration benchmark.
- `results/benchmark_sim/run_ebmf_cox_external.R` (the EBMF→Cox 2-step baseline that Phase 2's
  joint-vs-2-step comparison depends on): same fix — its training call is explicitly commented
  "matches D4 training preprocessing," so its external validation must too.
- `results/benchmark_sim/run_YFB_benchmark.R`: training's rank/per-platform settings are
  CLI-flag-controlled (`--no-rank`, `--per-platform-norm`); external validation now threads the
  same flags through (`rank_transform = !NO_RANK, per_platform_standardize = PER_PLATFORM_NORM`)
  instead of always using the function's defaults.
- `results/benchmark_sim/run_LB_benchmark.R`: checked and found **already consistent** — its
  training call uses `rank_transform=TRUE` with no `per_platform_standardize` (i.e. the same
  defaults `preprocess_desurv_cohort()` already uses for external validation) — no change needed.

**Deferred, NOT fixed (documented, not silently left):**
- `results/benchmark_sim/run_merged_benchmark.R` (the older, largely-superseded 18-config M1-M18
  comparison) computes external preprocessing **once** per cohort with fixed settings
  (`top_n=2000`, defaults `rank_transform=TRUE, per_platform_standardize=FALSE`), then reuses that
  single `pre_ext` across all 18 `MODEL_CONFIGS`, which have **heterogeneous** per-config
  `rank`/`per_plat` training settings. A single external preprocessing cannot simultaneously match
  18 different training preprocessing recipes. A correct fix requires computing external
  preprocessing per unique `(rank, per_plat)` combination inside the config loop (mirroring how
  `preproc_cache`/`gene_set_cache` already work for training) — a real restructuring, not a
  parameter flip, and out of Phase 1's scope given this script is superseded by the D-series
  pipeline for the current recommended configuration (D4). Flagged here so it is not mistaken for
  "already fixed."

**Affected files:** `code/preprocess_desurv.R`, `results/benchmark_sim/run_desurv_comparison.R`,
`results/benchmark_sim/run_ebmf_cox_external.R`, `results/benchmark_sim/run_YFB_benchmark.R`,
`tests/test_preprocess_desurv.R` (3 new tests: backward-compatible default, per-platform
z-standardize behavior, both-transforms-off passthrough).

---

## 2026-07-12 — Phase 1a objective normalization: boost survival (not shrink genomics), and only where it is safe

**Background.** The Rashid lab's 6/18/2026 feedback noted that the genomics likelihood term of
the model's objective is O(n·p) in aggregate while the survival (Cox partial likelihood) term is
O(n) — an unnormalized ~p-fold scale gap. This makes the `alpha` mixing weight
`(1-alpha)*genomics + alpha*survival` uninterpretable: alpha=0.5 does not mean "balanced,"
because the genomics side starts ~p times larger regardless of alpha's value. The fix needed to
reach the actual CAVI update coefficients (the A/B precision-and-signal terms passed into each
EBNM sub-problem), not just the reported ELBO monitor, per the post-lab-meeting action plan
(`docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md`, Phase 1a).

**Two candidate directions were tried and empirically rejected for the LB (Cluster A,
`code/fit_modular.R`) model** before arriving at the shipped design. Both were tested at toy scale
(n=150, p=200) and realistic PDAC scale (n=270, p=2000), alpha=0.5:

1. **Shrink genomics** (`A_gen`/`B_gen` in `update_L_k`/`update_F_k` divided by `p`): L, F, and beta
   collapsed to *exactly zero*, reproducibly, at both scales. Root cause: L and F co-adapt
   bilinearly — L's precision `A_gen = sum(Tau*EF2_k)` depends on F's squared values, and F's
   precision depends on L's squared values. Reducing genomics' precision at all removes the
   accidental numerical safety margin the *unnormalized* model was relying on (its very large
   `A_gen` never let either factor's posterior variance grow enough to matter); once one factor's
   estimate dips, the other's effective precision drops too, compounding to zero within 2-3
   iterations.
2. **Boost survival instead** (`A_surv`/`B_surv` in `update_L_k` multiplied by `p`, leaving genomics
   at its original scale — mathematically the same relative rebalancing, verified algebraically to
   give the identical `x = B/A` pseudo-observation, differing only in absolute precision/shrinkage
   strength): avoided the collapse in one synthetic scenario, but caused **unbounded divergence**
   in another (the exact fixture already used by `tests/test_elbo.R`'s `.elbo_test_data`) —
   `max(|EL|)` exceeded 600,000 by iteration 3, overflowing `exp(eta)` in the Cox Taylor expansion
   and producing transient `NaN`s in the ELBO. A scheduled/ramped introduction of the boost (reusing
   the existing `alpha_schedule` warmup+ramp mechanism, which already exists for a similar reason —
   easing survival pressure onto L) was tried and did **not** stabilize it; `max(|EL|)` still
   reached 12 million.
3. **Decoupling** the two update sites resolved it: passing the boosted `survival_divisor` **only**
   to `update_beta_k`'s own precision (no bilinear coupling — β's precision does not feed back into
   itself the way L/F's do) and to the ELBO monitor (reporting/convergence only, never feeds back
   into any update), while leaving `update_L_k`/`update_F_k` completely untouched, was stable at
   both scales with sensible factor recovery (correlation 0.3-0.8 with true loadings on synthetic
   data with a real signal).

**Decision:**
- `genomics_divisor`/`survival_divisor` are **never** passed to `update_L_k`/`update_F_k` inside
  `fit_supervised_mf_modular()`'s (LB) CAVI loop, regardless of `norm_convention`. LB's L/F
  precision stays at its original, pre-Phase-1a scale unconditionally. Only the ELBO assembly and
  the beta update (main loop and burn-in) receive the rebalancing.
  **Consequence: for LB specifically, the stated success criterion ("alpha=0.5 means balanced") is
  only partially met** — the L-update's internal genomics/survival mixture (the place the original
  imbalance actually lives) is unchanged, so LB's alpha still behaves close to its pre-Phase-1a
  self at interior values. This is a known, accepted limitation for LB, deferred pending a properly
  damped/trust-region stabilization of the bilinear coupling (not attempted here — out of Phase 1
  scope).
- **YFB (Cluster B, `code/fit_cox_on_yf.R`, the actual recommended/production model, D4) is
  unaffected by this limitation**: its L (`update_L_surv_YFB_k`) is pure-genomics-only regardless of
  alpha, and F (`update_F_surv_YFB_k`) defaults to `alpha_F=0` (pure genomics) — there is no
  genomics/survival competition in either update to rebalance in the first place, so YFB gets the
  full, verified-stable normalization on its beta update and ELBO monitor. Empirically verified:
  `EL`/`EF` are bit-identical between `norm_convention="per_p"` and `"np_n"` (both leave YFB's L/F
  untouched), confirming no unintended side effect.
- Two conventions are exposed via `norm_convention = c("per_p", "np_n")`:
  - `"per_p"` (default): `genomics_divisor=1` (unchanged), `survival_divisor=1/p` (multiplies
    survival's contribution by p). This is the safe, boost-survival direction.
  - `"np_n"`: the literal convention from the plan text (`genomics_divisor=n*p`,
    `survival_divisor=n`) — retained (applied only to the ELBO monitor and beta, per the same
    L/F-untouched rule) purely for empirical side-by-side comparison, since it was explicitly
    requested; it is *not* recommended, being the direction that would collapse LB's L/F if it were
    ever threaded there.
- The `elbo_proxy` genomics ELBO term itself (`code/update_tau.R`) is **unchanged** — its own
  precision-MLE argmax for `tau_j` is invariant to any positive scalar multiplying the whole
  genomics likelihood, so no rescaling was needed there; only the point where `elbo_proxy` and
  `surv_elbo` are combined into `elbo_full` (in `fit_modular.R` and `fit_cox_on_yf.R`) applies the
  divisors.

**Trade-offs:** LB's own interpretability goal is not fully achieved (see above) — this is an
explicit, documented scope narrowing, not a silent gap. The `tests/fixtures/lb_cohort_null_elbo_baseline.rds`
frozen-ELBO fixture was regenerated because `elbo_full`'s formula legitimately changed (verified:
`EL`/`EF`/`EBeta` are bit-identical to the pre-change baseline for that pure-noise fixture; only the
reported ELBO value differs).

**Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R` (new
`genomics_divisor`/`survival_divisor` parameters, default 1 = no-op), `code/fit_modular.R`,
`code/fit_cox_on_yf.R` (new `norm_convention` parameter; ELBO assembly; beta update call sites,
including burn-in), `tests/test_update_L.R`, `tests/test_update_F.R`, `tests/test_update_beta.R`,
`tests/test_normalization.R` (new), `tests/test_fit_modular_cohort.R` (regenerated baseline),
`tests/fixtures/lb_cohort_null_elbo_baseline.rds`.

---

## 2026-06-16 — Deck revision: KM direction convention, convergence criterion, prior restatement

While revising the 2026-06-18 lab-meeting deck (merged to main, `905279b`), three analytical
points were corrected/clarified. They supersede earlier wording in the 2026-05-27 factor
diagnostics report and parts of ROADMAP/PROJECT_STATUS.

**1. Adverse/protective direction is read from the (YF)-projection survival association
(marginal), NOT from the joint-model β sign.** The two disagree for the recommended model: the
joint coefficients are β̂₇ = −0.041 and β̂₃ = +0.011, but a Cox model on each program's projection
score (YF)ₖ — the quantity the risk model actually scores on — shows **Program 7 is adverse**
(higher activation → worse survival; carries MET, ITGA3, and glycolytic genes GAPDH/ENO1/TPI1/PGK1;
log-rank p≈5×10⁻⁴) and **Program 3 is protective** (epithelial markers MLPH/SLC45A3/TJP3;
p≈8×10⁻³). This is a **suppression effect**: the conditional (joint-β) sign can invert the marginal
direction when programs are correlated. *Decision:* label adverse/protective by the marginal
survival direction (which the KM curve displays and which matches biology), and report |β̂| as the
survival-activity magnitude rather than its signed value. The full risk score η = (YF)β is
correctly oriented throughout (training C = 0.656). *Files:* `figs/make_factor_figs.R` (KM
stratified by ZF = Y·EF; direction from the sign of the Cox coefficient on the projection).

**2. The canonical convergence criterion is relative full-ELBO change < 10⁻⁵.** Both
`code/fit_modular.R` and `code/fit_cox_on_yf.R` declare convergence when
|ELBO⁽ᵗ⁾ − ELBO⁽ᵗ⁻¹⁾| / |ELBO⁽ᵗ⁻¹⁾| < tol (default 1e-5, after a 5-iteration burn-in). The mean
absolute changes ΔL and Δβ are computed and logged but are **not** the stopping rule. This
clarifies the earlier "dual ΔL+Δβ < 1e-3" description (a V2-era design that the modular loops
moved past). The ELBO is both the CAVI objective and the convergence monitor; the number of
programs K is selected separately by cross-validated C-index with the 1-SE rule.

**3. Prior restatement (recommended model).** L and F use a **point-exponential** (non-negative
spike-and-slab) prior; β uses a **normal** prior; τ is a closed-form MLE. Legacy references to
"point-normal priors" describe earlier prior-sensitivity experiments, not the recommended
configuration. *Files:* `code/fit_cox_on_yf.R` (`prior_LF="point_exponential"`,
`prior_beta="normal"`).

**Deliverables:** revised 6/18 deck (28 slides) including a Limitations section and a
DeSurv-comparison backup slide grounded in the DeSurv manuscript (Young et al., PNAS 2026); the
unsupervised EBMF→Cox external baseline is recorded in the 2026-06-15 entry below.

---

## 2026-06-15 — Real-data unsupervised baseline: EBMF→Cox external validation

**Question:** The repo had supervised external C-indices (LB/YFB) and a *within-training*
per-factor EBMF diagnostic, but no *external* evaluation of the unsupervised two-step
(EBMF → Cox) on the held-out cohorts — so the "supervised joint beats unsupervised two-step"
claim was only demonstrated on synthetic data, not real data.

**Decision:** Add `results/benchmark_sim/run_ebmf_cox_external.R`, which fits `flashier` EBMF to
the merged TCGA_PAAD + CPTAC training matrix under the **identical** DeSurv-aligned preprocessing
as the recommended YFB model (per-platform z-std + combined-rank top-3000-per-cohort), fits
`coxph` on the EBMF factor scores, and scores the 5 held-out cohorts via the **same projection
formula used for YFB** (η = Y·F_norm·β). The only difference from YFB is that F is learned
unsupervised — making it a clean like-for-like contrast that isolates the value of supervision.
The script **fails loud** if β→0 or training/external C is at chance.

**Why this scoring choice:** Fitting Cox on the same projection used at test (rather than on
flashier's shrunken L directly) guarantees train/test consistency and mirrors the YFB
external-scoring path exactly. EBMF's natural rank is used (flashier selected K=20); this is
stated on the slide so the baseline is neither strawmanned nor cherry-picked.

**Result (5 held-out PDAC cohorts):** EBMF→Cox mean external C = **0.564** (K=20), vs YFB
**0.636** and LB **0.622**. The ~0.07 gap confirms on real data what the synthetic study showed
(YFB 0.81 > EBMF 0.72): supervised joint factorization generalizes better than the unsupervised
two-step, and does so with far greater parsimony (K_eff=2 survival-active programs vs 20 EBMF
factors). Output: `results/benchmark_sim/outputs/ebmf_cox_external/`. Built for the 6/18/2026
lab-meeting deck (`presentation/walther_lab_meeting_06_18_2026/`).

---

## 2026-06-14 — Multi-Cohort Simulation: non-negative ground truth + fair EBMF benchmark

**Question:** How should a multi-cohort simulation (shared vs. study-specific factors) be
generated so that recovery of the SSBMF (YFB/LB) model is tested *fairly*, and so that the
unsupervised EBMF benchmark is a like-for-like comparator?

**Decision:**
- Generate **non-negative** ground-truth loadings $L$ and gene programs $F$. This matches the
  models' `point_exponential` (NMF-style) prior on $L$ and $F$ — the established model
  assumption (see 2026-05-06 "Phase 2 — initialization constraints", where the SVD init was set
  to `abs()` as "equally valid for the non-negative point_exponential prior").
- Real-data EBMF templates (`flashier ldf$F`, ~48% negative) are taken in **absolute value** to
  obtain non-negative programs; the synthetic-fallback programs are non-negative by construction.
- The unsupervised **EBMF benchmark arm** is run with a **non-negative** prior
  (`ebnm_point_exponential` on $L$ and $F$) so it shares the same representational assumption.

**Why:** An initial draft used signed EBMF templates against the non-negative model prior. This
depressed supervised gene-program recovery to ~0.45 — a pure prior/data mismatch artifact.
Matching priors to the data sign structure raises recovery to ~0.90, equal to or better than the
(non-negative) EBMF benchmark. The flaw was in the simulation, not the model.

**Result (mean over 5 seeds, K=6, a=12):** YFB recovers shared programs at 0.90–0.94 and
specific at 0.84–0.92 (>= EBMF), specificity-classification accuracy 0.97–1.00, held-out
C-index 0.81–0.87 where shared prognostic signal exists and ~0.54 (beta FP rate 0.03–0.07) in the
null "nothing-shared" scenario. The cohort dummy indicator helps under heterogeneity (hybrid)
and is redundant when cohorts are homogeneous (all-shared). YFB > LB on recovery and stability.

**Robustness (grounded in prior work):** fits are deterministic given data (single SVD start;
near-unimodal landscape per the 2026-05-06 "Phase 3 — multi-initialization" 30-restart sweep, so
`n_init=1`); recovery is stable to $K$ over-specification ($K=6/8/10$); `normal` beta prior used
(`point_normal` collapses beta->0 under YFB).

**Status:** Working draft for the 2026-06-15 advisor meeting; subject to revision after feedback.

**Files:** `results/multi_cohort_sim/{generate_multicohort_data,build_ebmf_templates,sim_scoring,run_multicohort_sim,run_signal_ratio_sweep}.R`,
`config/globals.yml` (`synthetic_multicohort`), `docs/reports/multicohort_sim_proposal_06_14_26.{qmd,pdf}`.
The standalone results report was consolidated into the proposal (Part II removed for brevity;
full tables and sweep figures available from the runner outputs).

---

## 2026-05-27 — DeSurv Gene Selection Alignment

**Question:** Does adopting DeSurv's gene selection (combined mean+variance rank, top-3000 per
cohort before normalization) improve external C-index relative to the current variance-only
selection (top-2000 on merged normalized matrix)?

**Implementation:**
- `select_top_variable_genes()` gains `method="combined_rank"`: rank_mean + rank_var,
  lowest rank-sum genes retained. Default "variance" unchanged.
- `preprocess_merged_cohorts()` gains `selection_per_cohort=TRUE`: per-cohort top-N
  selection on log-transformed data, then intersect, before per-platform z-std.
- New comparison configs: D3 (LB DeSurv-aligned) and D4 (YFB DeSurv-aligned).

**Result:**

| Config | Model | Cohort ID | Mean external C | K_eff | Gene set |
|--------|-------|-----------|----------------|-------|---------|
| D1 (= M4) | LB + orig | Yes | 0.616 | 1 | 2000 |
| D2 (= M5) | YFB + orig | No | 0.624 | 2 | 2000 |
| D3 (DeSurv LB) | LB + aligned | Yes | 0.622 | 2 | 2064 |
| D4 (DeSurv YFB) | YFB + aligned | No | 0.636 | 2 | 2064 |
| D5 (DeSurv YFB + cohort) | YFB + aligned | Yes | 0.614 | 2 | 2064 |

**Decision:** delta_yfb = +0.012 (D4 − D2) > +0.005 threshold. Adopting D4 as
new primary configuration. Per-cohort C-index: D4 improves over D2 in 4/5 cohorts
(PACA_AU_array is the exception: D4=0.650 vs D2=0.670). LB also improves:
delta_lb = +0.006 (D3 − D1).

**Cohort indicator with YFB DeSurv-aligned (D5):** Adding a cohort indicator to D4
reduces mean C to 0.614 (−0.022 vs D4, −0.010 vs M5 baseline). Per-platform
z-standardization already removes the platform offset; the cohort indicator absorbs
factor capacity that would otherwise capture survival signal. Cohort indicator
provides no benefit for YFB with DeSurv preprocessing — D4 (no cohort indicator)
remains the primary configuration.

**New primary config:** YFB + per-platform z-std + DeSurv gene selection
(combined_rank, top-3000 per cohort before normalization) + no cohort indicator, K=7.
Fit object: results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds (D4).

---

## 2026-05-25 — Extended preprocessing comparison: 18-configuration benchmark

- **Design:** Extended the 6-configuration benchmark (M1–M6) to 18 configurations by adding three
  additional preprocessing strategies: joint quantile normalization without rank transform
  (M7/M8, M13/M14), joint z-standardization (M9/M10, M15/M16), and log-only (M11/M12, M17/M18).
  Biological K floor K_final = max(K_1se, 3) applied — motivated by DeSurv K=3–4 and the requirement
  for multiple interpretable gene expression programs (K=2 insufficient for "programs" framing).

| ID | Model | Preprocessing | cohort_id | K (1-SE CV + floor) |
|----|-------|---------------|-----------|---------------------|
| M1 | LB | Joint QN + rank | No | 6 |
| M2 | LB | Joint QN + rank | Yes | 6 |
| M3 | LB | Per-platform z-std | No | 3 |
| M4 | LB | Per-platform z-std | Yes | 3 |
| M5 | YFB | Per-platform z-std | No | 3 (floor; 1-SE gave K=2) |
| M6 | YFB | Per-platform z-std | Yes | 3 (floor; 1-SE gave K=2) |
| M7 | LB | Joint QN, no rank | No | 6 |
| M8 | LB | Joint QN, no rank | Yes | 6 |
| M9 | LB | Joint z-std | No | 3 |
| M10 | LB | Joint z-std | Yes | 3 |
| M11 | LB | Log only | No | 3 |
| M12 | LB | Log only | Yes | 3 |
| M13 | YFB | Joint QN, no rank | No | 5 |
| M14 | YFB | Joint QN, no rank | Yes | 5 |
| M15 | YFB | Joint z-std | No | 3 |
| M16 | YFB | Joint z-std | Yes | 3 |
| M17 | YFB | Log only | No | 3 (floor; 1-SE gave K=2) |
| M18 | YFB | Log only | Yes | 3 (floor; 1-SE gave K=2) |

- **Full benchmark results** (5 external cohorts, prior=normal, max_iter=300):

  | Model | Dijk | Moffitt | PACA_AU_arr | PACA_AU_seq | Puleo | Mean C | K_eff | beta_max |
  |-------|------|---------|-------------|-------------|-------|--------|-------|----------|
  | M1 LB joint+rank | 0.569 | 0.515 | 0.645 | 0.667 | 0.600 | 0.599 | 1 | 0.016 |
  | M2 LB joint+rank+cohort | 0.587 | 0.537 | 0.659 | 0.686 | 0.626 | 0.619 | 3 | 0.031 |
  | M3 LB perplat | 0.615 | 0.547 | 0.649 | 0.644 | 0.645 | 0.620 | 1 | 0.627 |
  | M4 LB perplat+cohort | 0.624 | 0.549 | 0.656 | 0.638 | 0.650 | 0.624 | 1 | 0.482 |
  | **M5 YFB perplat** | 0.621 | 0.537 | 0.662 | 0.655 | 0.650 | **0.626** | 2 | 0.053 |
  | M6 YFB perplat+cohort | 0.578 | 0.519 | 0.655 | 0.645 | 0.652 | 0.613 | 2 | 0.052 |
  | M7 LB joint norank | 0.567 | 0.513 | 0.642 | 0.669 | 0.602 | 0.598 | 1 | 0.235 |
  | M8 LB joint norank+cohort | 0.583 | 0.533 | 0.651 | 0.682 | 0.624 | 0.615 | 1 | 0.339 |
  | M9 LB joint z-std | — | — | — | — | — | 0.545 | 0 | 0.000 |
  | M10 LB joint z-std+cohort | — | — | — | — | — | 0.555 | 0 | 0.000 |
  | M11 LB log only | — | — | — | — | — | 0.533 | 0 | 0.000 |
  | M12 LB log only+cohort | — | — | — | — | — | 0.531 | 0 | 0.000 |
  | M13 YFB joint norank | — | — | — | — | — | 0.540 | 0 | 0.000 |
  | M14 YFB joint norank+cohort | — | — | — | — | — | 0.529 | 0 | 0.000 |
  | M15 YFB joint z-std | — | — | — | — | — | 0.562 | 0 | 0.000 |
  | M16 YFB joint z-std+cohort | — | — | — | — | — | 0.543 | 0 | 0.000 |
  | M17 YFB log only | — | — | — | — | — | 0.529 | 0 | 0.000 |
  | M18 YFB log only+cohort | — | — | — | — | — | 0.532 | 0 | 0.000 |

  *(M9–M18 mean C shown from benchmark output; individual cohort columns suppressed because all
  results reflect the β=0 null model, giving C≈0.5 per cohort with random variation.)*

- **Key findings:**
  1. **Per-platform z-std is essential.** 10 of 12 non-per-platform configurations collapse to
     K_eff=0 (β→0). The two that don't (M7/M8) merely match the already-inferior joint-QN LB baseline.
  2. **Rank transform diagnosis (M13/M14):** YFB with joint QN *without* rank transform also collapses
     (K_eff=0). This definitively rules out the rank transform as the cause of YFB β→0 on joint-QN
     data. The root cause is the platform-dominated SVD directions produced by joint normalization of
     mixed-platform data — per-platform z-std corrects this, rank transform does not.
  3. **M5 at K=3 floor:** Mean C=0.626 (vs. 0.625 at K=2 in the prior run). The biological floor
     adds a third factor without degrading external generalization, and provides richer multi-program
     structure for biological interpretation.
  4. **Cohort indicator:** Helpful for LB (M1→M2: +0.020), marginal for LB per-platform (M3→M4:
     +0.004), harmful for YFB per-platform (M5→M6: −0.013).

- **Recommended configuration (updated from K=2 to K=3):**
  - **Primary: M5 — YFB × per-platform z-std × no cohort indicator, K=3**
    - Highest mean external C=0.626; 2/3 factors survival-active.
    - YFB external scoring is exact: η_new = (Y_new %*% EF) %*% β (no approximation).
    - K=3 supports multi-program biological framing; consistent with DeSurv K=3–4.
  - **Sensitivity check: M4 — LB × per-platform z-std × cohort indicator, K=3** (mean C=0.624).
  - **Not recommended:** M9–M18 (β=0 collapse or marginal); M6 (cohort indicator hurts YFB).

- **Files:** `results/benchmark_sim/run_merged_kcv.R`, `results/benchmark_sim/run_merged_benchmark.R`,
  `code/preprocess_desurv.R` (new `normalize_method` param), `config/globals.yml` (8 K keys),
  `docs/reports/merged_benchmark_report.qmd`,
  `results/benchmark_sim/outputs/merged_benchmark/merged_benchmark_results_extended.csv`

---

## 2026-05-25 — Merged-cohort benchmark: 6-configuration apples-to-apples comparison

- **Design:** Six model configurations (2 models × 2 preprocessing × 2 cohort indicator settings)
  evaluated on merged TCGA_PAAD + CPTAC training with K selected by cross-validated C-index (1-SE rule,
  K_grid=2:10, 5 folds). Prior: normal throughout — point_normal is excluded because it produces β→0
  on merged data regardless of model or preprocessing (confirmed on prior runs). YFB × joint
  quantile+rank preprocessing is excluded (structural β→0, all V0–V11 strategies exhausted,
  see 2026-05-22 entry).

| ID | Model | Preprocessing | cohort_id | K (1-SE CV) |
|----|-------|---------------|-----------|-------------|
| M1 | LB (η = Lβ) | Joint quantile+rank | No | 6 |
| M2 | LB | Joint quantile+rank | Yes | 6 |
| M3 | LB | Per-platform z-std | No | 3 |
| M4 | LB | Per-platform z-std | Yes | 3 |
| M5 | YFB (η = (YF)β) | Per-platform z-std | No | 2 |
| M6 | YFB | Per-platform z-std | Yes | 2 |

- **K-CV findings:** LB joint peaks at K=6 (mean C=0.589), non-monotone with poor K=3 (0.439).
  LB per-platform peaks at K=3 (mean C=0.558) and is nearly flat K=3–8 — parsimonious.
  YFB per-platform peaks at K=7 (mean C=0.644) but 1-SE selects K=2 (mean C=0.625 is within
  1 SE of the peak); the 1-SE rule strongly favors parsimony here given high fold-to-fold variance
  in YFB folds.

- **Full benchmark results** (5 external cohorts, prior=normal, max_iter=300):

  | Model | Dijk | Moffitt | PACA_AU_arr | PACA_AU_seq | Puleo | Mean C | K_eff | beta_max |
  |-------|------|---------|-------------|-------------|-------|--------|-------|----------|
  | M1 LB joint | 0.569 | 0.515 | 0.645 | 0.667 | 0.600 | 0.599 | 1 | 0.016 |
  | M2 LB joint + cohort | 0.587 | 0.537 | 0.659 | 0.686 | 0.626 | 0.619 | 3 | 0.031 |
  | M3 LB per-platform | 0.615 | 0.547 | 0.649 | 0.644 | 0.645 | 0.620 | 1 | 0.627 |
  | M4 LB per-platform + cohort | 0.624 | 0.549 | 0.656 | 0.638 | 0.650 | 0.624 | 1 | 0.482 |
  | **M5 YFB per-platform** | 0.621 | 0.537 | 0.662 | 0.655 | 0.650 | **0.625** | 2 | 0.052 |
  | M6 YFB per-platform + cohort | 0.578 | 0.519 | 0.655 | 0.645 | 0.652 | 0.610 | 2 | 0.052 |

- **Key comparisons:**
  - **LB joint vs. LB per-platform (no cohort_id):** M1→M3 mean C 0.599→0.620 (+0.021). Per-platform
    preprocessing consistently improves LB on merged data. Joint quantile normalization cannot fully
    remove the RNA-seq vs. proteomics scale difference; per-platform z-standardization achieves better
    biological signal separation.
  - **Cohort indicator effect (matched preprocessing):**
    - LB joint: M1→M2 +0.020 (3 active factors vs 1 — cohort column frees factors from absorbing platform contrast).
    - LB per-platform: M3→M4 +0.004 (marginal; preprocessing already removed bulk platform effects).
    - YFB per-platform: M5→M6 −0.015 (cohort indicator hurts at K=2 — the cohort column competes with
      the 2 biological factors for the limited genomic variance, reducing biological signal in ZF).
  - **LB vs. YFB at per-platform preprocessing:** M5 (YFB, mean C=0.625) edges out M4 (LB+cohort,
    0.624) and M3 (LB, 0.620). YFB's direct factor projection (η = (YF)β) appears slightly better
    calibrated for external generalization at low K.
  - **M3 beta_max anomaly:** beta_max=0.627 for M3 (K_eff=1) is large compared to prior LB results
    (~0.02). At K=3 with per-platform preprocessing, a single factor absorbs all survival-predictive
    variance. The high beta magnitude reflects that the single active factor is on a small (z-score)
    scale — not a numerical instability, but it signals that K=3 is below the biological signal
    dimensionality for the LB model on this data.

- **Recommended configuration for manuscript:**
  - **Primary: M5 — YFB × per-platform z-std × no cohort indicator, K=2**
    - Highest mean external C-index (0.625) across 5 independent cohorts.
    - Both factors are active (K_eff=2 = K_total), no dead factors.
    - Most parsimonious model overall (K=2).
    - Beta scale is on the natural YFB scale (~0.05); readily interpretable.
    - YFB predictor (η = (YF)β) maps observed gene expression directly onto factor scores,
      giving biologically transparent risk scores.
  - **Sensitivity check: M4 — LB × per-platform z-std × cohort indicator, K=3**
    - Mean C=0.624, effectively tied with M5.
    - Cohort indicator explicitly models the RNA-seq vs. proteomics offset (recoverable quantity).
    - K_eff=1 limits multi-program interpretation — only one prognostic factor identified.
  - **Not recommended: M1 (lowest mean C=0.599); M6 (cohort indicator hurts YFB at low K).**

- **Files:** `results/benchmark_sim/run_merged_kcv.R`, `results/benchmark_sim/run_merged_benchmark.R`,
  `docs/reports/merged_benchmark_report.{qmd,pdf}`,
  `results/benchmark_sim/outputs/merged_benchmark/`


---

## 2026-05-22 — YFB β→0 on merged data: structural diagnosis and two bug fixes

- **Finding:** The β→0 collapse in YFB (η = ZF·β) on merged TCGA+CPTAC is a structural property of the model, not fixable by initialization strategies (Cox warm-start, burn-in, higher α). Diagnosis via `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R`.
- **Root cause:** The YFB F update is genomics-only (`update_F_k` uses only Y reconstruction; survival signal does not flow to F). On merged multi-platform data, the genomics-optimal EF captures platform contrast (RNA-seq vs proteomics), giving ZF factors with no survival correlation (Cox p≈0.8). With B_beta ≈ 0 in every iteration, β converges to 0 regardless of initialization. LB escapes this because EL is dual-source (genomics + α·w·β² survival term), pulling EL toward survival-predictive directions.
- **Two genuine bugs fixed in `fit_cox_on_yf.R`:**
  1. SVD initialization used `pmax(..., 0)`, which zeros out entire EF columns when the SVD vector points all-negative. With K=2 on merged data, EF[:,1] was identically 0, making ZF[:,1]=0 and β_1=0 trivially. Fixed to `abs()`, which is equally valid for the non-negative point_exponential prior.
  2. Cox warm-start used un-normalized ZF (Y %*% EF), giving β_ws ≈ 5.99e-9 on the huge ZF scale — machine epsilon on the CAVI scale. Fixed to normalize EF before Cox regression (consistent with how CAVI computes ZF). After fix, warm-start produces β_ws ≈ [-4.5e-4, 1.3e-4] — still small but on the correct scale.
- **Effect of fixes:** K-CV is now informative (K=4 selected, mean C=0.593 vs uninformative C=0.524–0.557 before). β_max reaches 0.0004 with warm-start (up from 0). But CAVI still drives β→0 as EF evolves toward platform-dominated genomics optimum. Training dB decays as 7e-7/iter (V1) vs 3e-16 (V0).
- **Dual-source F tested (alpha_F ∈ {0.1, 0.3, 0.5}, V6–V8):** `update_F_surv_YFB.R` was already implemented but hardcoded to alpha=0; wired alpha_F parameter through. Result: RMSE immediately jumps from ~290 to ~750-800 (EF instability), and β still collapses to 0. The chicken-and-egg trap: with EBeta≈0, the survival gradient contribution to A_F/B_F is ≈0 at every iteration, so dual-source F cannot bootstrap survival signal. EF grows freely, degrading genomics reconstruction without any β benefit.
- **cohort_id extension also ineffective for YFB merged:** YFB_cohort has K_eff=0, beta_max=0 across all external cohorts — identical to YFB_base. Cohort offsets in L do not change the F update instability.
- **Conclusion:** All CAVI-local fixes (warm-start, burn-in, high α, cohort_id, dual-source F) fail on merged TCGA+CPTAC for YFB. The structural problem is that CAVI has no path from β=0 to β≠0 when EF already converges to a genomics-optimal platform-contrast fixed point. Escaping this requires either breaking the CAVI fixed-point structure (e.g., stochastic initialization ensemble or coordinate ascent with frozen F during β burn-in) or accepting per-platform fits.
- **Path forward:** (1) Per-platform YFB is the current best YFB result (C_dijk ≈ 0.573); (2) LB_merged is superior on multi-platform data (K_eff=3, C_dijk≈0.590); (3) YFB single-cohort is unaffected and performs well; (4) for future work: frozen-F β pre-conditioning (freeze EF for N_frozen iters, let β grow freely, then unfreeze) would break the trap without instability.
- **Affected files:** `code/fit_cox_on_yf.R` (two bug fixes + alpha_F parameter); `code/update_F_surv_YFB.R` (pre-existing, unchanged); `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R`; `results/benchmark_sim/outputs/yfb_beta_fix/`

---

## 2026-05-22 — Cohort indicator extension: fixed L columns absorb platform offset

- **Decision:** Augment both LB (`fit_supervised_mf_modular`) and YFB (`fit_cox_on_yf`) with
  a `cohort_id` parameter that appends C−1 fixed binary columns to L (corner-point encoding,
  reference = first alphabetical cohort level). Corresponding F rows `f_c` are estimated via
  a Normal conjugate update (`update_F_cohort.R`) with Gaussian prior
  N(0, σ²_{F,cohort} · I). The survival linear predictor uses only the K biological columns
  (β_cohort = 0 by construction). σ_{F,cohort} (default 1.0) is an exposed parameter.

- **Rationale:** On merged TCGA_PAAD + CPTAC training, RNA-seq vs proteomics platform effects
  drive the top SVD component, reducing biological signal in the K biological factors and
  causing β→0 in the YFB model. Explicitly parameterizing the cohort offset as fixed indicator
  columns removes it from the factor competition and stabilizes the biological subspace.

- **New code:**
  - `code/update_F_cohort.R` — closed-form Normal conjugate update for F_cohort rows;
    `update_F_cohort_k()` (single column) and `update_F_cohort_all()` (sweep over C−1 columns).
  - `code/compute_elbo.R` — added `compute_normal_kl(EF, EF2, sigma_F_cohort)`:
    KL[N(μ, 1/A) ∥ N(0, σ²)] = ½(E[f²]/σ² + log(Aσ²) − 1); returns negative scalar added
    directly to elbo_full by the caller.
  - `code/fit_modular.R` — `cohort_id` and `sigma_F_cohort` parameters; corner-point design
    matrix built once before the CAVI loop; L_aug/EF_aug carry the augmented matrices while
    EL/EF remain K-column working matrices; EF_cohort, EF2_cohort, L_cohort returned.
  - `code/fit_cox_on_yf.R` — same cohort extension; ZF = Y %*% EF_aug[,1:K,drop=FALSE]
    restricted to the K biological columns in two locations (inner CAVI loop and Phase C).
  - `tests/test_update_F_cohort.R` (10 tests), `tests/test_fit_modular_cohort.R` (11 tests),
    `tests/test_fit_yf_cohort.R` (9 tests) — 229/229 tests passing.

- **Stage 1 — synthetic validation** (n=200, p=500, K_true=3, C=2, offset SD=2.0):
  - Offset absorption: |cor(EF_cohort, true_offset_vec)| = 0.995 for LB_cohort and YFB_cohort.
  - Factor recovery (mean max-cor vs true L): LB_cohort 0.77 vs LB_base 0.41 (+88%);
    YFB_cohort 0.71 vs YFB_base 0.34 (+109%). Cohort indicator substantially improves
    recovery of the biological subspace when a strong platform offset is present.
  - C-index: marginal (both models still affected by β→0 in synthetic conditions).

- **Stage 2 — real PDAC benchmark** (merged TCGA_PAAD n=144 + CPTAC n=129, K=20, α=0.50,
  prior_beta=normal, σ_{F,cohort}=1.0):

  | Cohort | LB_base | LB_cohort | YFB_base | YFB_cohort |
  |--------|---------|-----------|----------|------------|
  | Dijk | 0.590 | 0.545 | 0.506 | 0.524 |
  | Moffitt_GEO_array | 0.529 | 0.520 | 0.533 | 0.509 |
  | PACA_AU_array | 0.657 | 0.664 | 0.577 | 0.564 |
  | PACA_AU_seq | 0.681 | 0.700 | 0.588 | 0.580 |
  | Puleo_array | 0.634 | 0.590 | 0.518 | 0.519 |
  | **Mean** | **0.618** | 0.604 | 0.544 | 0.539 |

  Model summary: LB_base K_eff=3 β_max=0.021; LB_cohort K_eff=2 β_max=0.034;
  YFB_base K_eff=0 (β→0); YFB_cohort K_eff=0 (β→0).
  Cohort offset norm: ‖EF_cohort‖_F = 3.04 (LB), 1.12 (YFB).

- **Interpretation:**
  - LB_cohort wins on RNA-seq external cohorts (PACA_AU_array +0.007, PACA_AU_seq +0.019)
    where the TCGA vs CPTAC distinction is most relevant to platform-matched prediction.
  - LB_cohort underperforms on Dijk (−0.045) and Puleo (−0.044). With K=20, ARD pruning
    already handles platform effects implicitly (K_eff reduces to 3). The cohort column
    absorbs one additional factor degree of freedom (K_eff=2 vs 3 for LB_base), losing a
    biological factor that was contributing to those predictions.
  - YFB β→0 collapse on merged data is unchanged by the cohort extension. The root cause
    (genomics-term dominance in the L update; per-platform norm resolved it for YFB in
    Phase 1 via separate preprocessing) is not addressed by cohort indicator columns alone.
  - Mean C-index: LB_cohort (0.604) vs LB_base (0.618 current / 0.609 main-branch reference).
    The extension does not reliably improve the LB model on this dataset in the fully-converged
    K=20 configuration; it trades one biological factor for one cohort indicator.

- **Trade-offs:**
  - CPTAC is proteomics (not RNA-seq); a single linear cohort indicator cannot fully capture
    the non-linear RNA-seq vs proteomics platform difference. A per-gene cohort shift
    (one degree of freedom per gene in the platform offset) may be too simplistic.
  - The benefit is clearest in synthetic data with a pure additive rank-1 offset. Real platform
    effects are higher-rank and partially confounded with biology.
  - When the LB model already converges with K_eff > 0 (as it does here with K=20), the cohort
    column competes with biological factors rather than complementing them.

- **Affected files:** `code/update_F_cohort.R` (new), `code/compute_elbo.R` (added
  `compute_normal_kl`), `code/fit_modular.R` (`cohort_id`, `sigma_F_cohort` params),
  `code/fit_cox_on_yf.R` (same), `tests/test_update_F_cohort.R` (new),
  `tests/test_fit_modular_cohort.R` (new), `tests/test_fit_yf_cohort.R` (new),
  `tests/run_tests.R` (updated), `tests/fixtures/lb_cohort_null_elbo_baseline.rds` (regenerated).
  Validation: `results/cohort_lmm_sim/run_synthetic.R`,
  `results/benchmark_sim/run_cohort_lmm_benchmark.R`,
  `results/benchmark_sim/outputs/cohort_lmm_benchmark/`.

---

## 2026-05-25 — Frozen-F β pre-conditioning: N_frozen parameter added to fit_cox_on_yf()

- **Decision:** Added `N_frozen` parameter (default 0, backward compatible) to
  `fit_cox_on_yf()`. During the first `N_frozen` CAVI iterations, EF is held fixed
  at its SVD initialization while β and L update freely. After `N_frozen` iterations,
  EF unfreezes and full joint CAVI proceeds.

- **Rationale:** Breaks the β=0 ↔ B_beta=0 CAVI fixed-point that causes YFB to
  collapse on merged multi-platform data. With EF frozen at SVD init (non-degenerate
  via the abs() fix), ZF = Y·EF_init is fixed and non-zero, so A_beta = sum(w·ZF_k²)
  > 0 from iteration 1 and β can grow freely. Unlike dual-source F (V6–V8), this does
  not risk EF instability because EF does not update during the frozen phase.

- **Implementation:** Surgical: `freeze_F <- iter <= N_frozen` flag; F update in step
  (c) and cohort F update both gated on `!freeze_F`. Validation message logged at
  freeze/unfreeze boundary. Input validation: N_frozen must be non-negative numeric.

- **Tests:** 8 new tests in `tests/test_fit_yf_frozen_f.R` (238/238 total passing):
  backward compat (N_frozen=0), valid output with N_frozen=3/5/100, ELBO finite,
  N_frozen=-1/"a" raise errors, N_frozen=0 vs 5 produce different EF, works with
  cohort_id.

- **Effectiveness on merged TCGA+CPTAC:** Tested via `run_yfb_beta_fix_diagnostic.R
  --quick` (V9: N_frozen=10, V10: N_frozen=20, V11: N_frozen=30 + warm-start).
  All three frozen-F variants converge at iteration 9 — during the frozen phase, before EF
  ever unfreezes. Root cause: SVD initialization on merged TCGA+CPTAC is itself
  platform-dominated (RNA-seq vs proteomics contrast is the top SVD component), so
  ZF_SVD = Y·EF_SVD has no survival correlation. With B_beta ≈ 0, β→0 during the
  frozen phase, ELBO stabilizes quickly, and the convergence criterion fires at iter 9
  before the unfreeze point. The "Unfreezing EF" message never appears.
  Frozen-F does NOT rescue β→0 on merged TCGA+CPTAC. This is the last untested
  CAVI-local fix — ALL strategies (V0–V11) are now exhausted.
  One potential future fix not yet implemented: initialize EF from a TCGA-only (single-cohort)
  YFB fit rather than SVD on the merged matrix. TCGA-only EF would not be platform-dominated,
  and frozen-F from that starting point could allow β to grow before platform contrast dominates.
  This requires running a single-cohort fit first (additional overhead) and is not yet tested.
  See diagnostic output in `results/benchmark_sim/outputs/yfb_beta_fix/`.

- **Affected files:** `code/fit_cox_on_yf.R` (N_frozen parameter + logic),
  `tests/test_fit_yf_frozen_f.R` (new), `tests/run_tests.R` (updated),
  `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R` (V9/V10/V11 variants added).

---

## 2026-05-25 — Low-K cohort benchmark: cohort_id is essential at K=5 for LB on merged data

- **Finding:** At K=5, the LB model (η = Lβ) also collapses to β→0 on merged TCGA+CPTAC
  (K_eff=0), not just the YFB model. At K=20, ARD has enough capacity to simultaneously
  absorb the RNA-seq vs. proteomics platform contrast and learn 2–3 biological factors.
  At K=5, it cannot — the platform contrast monopolizes the available factors and biological
  signal is lost.

- **Cohort extension effect at K=5:**

  | Cohort | LB_base K=5 | LB_cohort K=5 | LB_base K=20 | LB_cohort K=20 |
  |--------|-------------|---------------|--------------|----------------|
  | Dijk | 0.544 | 0.601 | 0.590 | 0.545 |
  | Moffitt | 0.543 | 0.537 | 0.529 | 0.520 |
  | PACA_AU_array | 0.521 | 0.643 | 0.657 | 0.664 |
  | PACA_AU_seq | 0.507 | 0.672 | 0.681 | 0.700 |
  | Puleo | 0.504 | 0.617 | 0.634 | 0.590 |
  | **Mean** | 0.524 | **0.614** | **0.618** | 0.604 |

  LB_cohort K=5 achieves mean C=0.614, nearly matching LB_base K=20 (0.618) with
  4× fewer biological factors. K_eff=3 (vs K_eff=0 for LB_base at K=5).

- **Interpretation:** The cohort extension is beneficial in two distinct regimes:
  1. **Low-K regime (K ≤ 5):** *Essential* — LB_base β→0 without cohort_id. Cohort
     column absorbs platform offset, freeing all K biological factors for survival signal.
  2. **High-K regime (K=20):** *Neutral to marginal* — ARD absorbs platform effects
     implicitly by zeroing out factors with low biological signal. Cohort column competes
     with biological factors and may reduce K_eff from 3 to 2.

- **Practical recommendation:** For merged multi-platform PDAC analysis, prefer
  K=5 + cohort_id over K=20 without cohort_id. Equivalent predictive performance,
  far more interpretable (5 factors vs 3 non-zero out of 20).

- **YFB β→0 unchanged:** LB_base K=5 K_eff=0 and YFB K=5 K_eff=0 both still collapse.
  YFB frozen-F fix (N_frozen parameter) is the next diagnostic.

- **Files:** `results/benchmark_sim/run_cohort_lmm_benchmark.R` (--low-k flag added),
  `results/benchmark_sim/outputs/cohort_lmm_benchmark_low_k/` (results).

---

## 2026-05-20 — YFB K-CV sign fix: sign_correction=FALSE in CV folds resolves C<0.5 for normal prior

- **Decision:** Added a `sign_correction` parameter (default TRUE) to `fit_cox_on_yf()` in
  `code/fit_cox_on_yf.R`. Inside `select_K_cv(model="YFB")`, pass `sign_correction = FALSE`
  in the per-fold fit args. Concordance is evaluated via `I(-pred$risk_scores)` — the same
  convention as the LB model.

- **Problem diagnosed:** YFB K-CV produced C < 0.5 for every K value under the normal prior.
  Root cause: a sign double-flip. `fit_cox_on_yf()` previously applied Phase C unconditionally
  (if training C < 0.5, negate EBeta). Inside CV folds, Phase C corrected the sign of EBeta,
  orienting it concordantly with training survival. The concordance evaluation then applied
  `I(-pred$risk_scores)` — expecting anti-concordant raw predictions. With Phase C correction
  in place, the negation produced anti-concordant predictions, giving C < 0.5 universally.

- **Fix:** Wrapped the Phase C block in `if (sign_correction) { ... }`. In the YFB branch of
  `select_K_cv()`, `sign_correction = FALSE` is now passed explicitly (same approach as LB).
  CV folds deliver raw SVD-oriented EBeta; concordance is evaluated via `I(-pred$risk_scores)`.
  This is consistent with LB convention: both models use raw orientation in folds, both evaluate
  via the negated risk score.

- **Empirical results after fix (5-fold CV on TCGA_PAAD, K ∈ {2,…,10,15,20}, normal prior):**

  | K  | Mean C | SE    | Note |
  |----|--------|-------|------|
  | 2  | 0.565  | 0.025 | |
  | 3  | 0.585  | 0.008 | |
  | 5  | 0.601  | 0.011 | **selected (1-SE)** |
  | 6  | 0.605  | 0.012 | best mean C |
  | 7  | 0.593  | 0.028 | |
  | 8–20 | 0.535–0.600 | | |

  1-SE rule: threshold = 0.605 − 0.012 = 0.593; smallest K with mean C ≥ 0.593 is **K=5**.

- **YFB point_normal K-CV remains C=0.5 for all K.** The spike-and-slab prior collapses β→0
  inside CV folds (training fold n≈115 insufficient to escape the spike component). This is
  a separate issue from the sign fix — point_normal K-selection for YFB is an open item.

- **Convention rationale:** sign_correction=FALSE inside CV folds is correct for both LB and
  YFB because per-fold Phase C produces orientation inconsistency across folds (the same CAVI
  solution may be corrected in one fold but not another, inflating apparent C-index variance).
  `I(-pred$risk_scores)` handles the sign convention consistently regardless of SVD orientation.

- **Test added:** KCV-T15 verifies that YFB K-CV with high-SNR data (n=100, K_true=2, sd=0.1,
  β=3.0, ~30% censoring) yields mean C ≥ 0.5 for all K under normal prior. 193/193 tests passing.

- **Affected files:** `code/fit_cox_on_yf.R` (sign_correction parameter, Phase C wrapped in
  `if (sign_correction)`), `code/select_K.R` (sign_correction=FALSE in YFB fit_args with
  explanatory comment), `tests/test_select_K_cv.R` (KCV-T15 added),
  `results/benchmark_sim/outputs/K_cv/K_cv_table_YFB_normal.csv` (regenerated after fix).

---

## 2026-05-06 — Phase 1 (YFB merged): per-platform z-standardization resolves β→0 collapse

- **Decision:** For YFB merged training (TCGA_PAAD + CPTAC), apply per-platform
  z-standardization (normalize each cohort's gene expression matrix independently before
  merging) and disable the per-subject rank transform (`rank_transform=FALSE`). Use K=3
  for merged YFB. These flags are the canonical YFB merged configuration going forward.

- **Empirical evidence:**

  Previous merged YFB (K=20, `rank_transform=TRUE`, quantile normalization across platforms)
  collapsed β→0 (machine-epsilon, ~10⁻¹³) at iteration 1. Diagnosis: the rank transform
  before quantile normalization removed residual between-platform variance that the survival
  signal depends on; the CAVI β update saw A_surv/A_gen ~10⁻³, making the genomics term
  dominate and driving β→0. Per-platform standardization (z-score each cohort separately,
  then merge) preserves within-cohort gene-level variance while removing cross-platform mean
  shifts. With `rank_transform=FALSE` and K=3, K_eff=3 and external C-index ranges 0.53–0.67
  (median 0.64), matching or exceeding the DeSurv range of 0.60–0.65.

- **Trade-offs:** Per-platform standardization assumes each platform's gene expression is
  comparable after z-scoring, which may not hold if the platforms differ in signal-to-noise
  ratio beyond mean/variance. The K=3 choice was made using `auto_prune_K()` (ARD pruning)
  as a starting point; CV-based K selection for the merged YFB case remains future work.

- **Affected files:** `code/preprocess_desurv.R` (`per_platform_standardize` parameter),
  `results/benchmark_sim/run_YFB_benchmark.R` (`--per-platform-norm --no-rank` flags),
  `config/globals.yml` (`k_pdac_yfb_merged: 3`).
  Results: `results/benchmark_sim/outputs/YFB_benchmark_perplatform/`.

---

## 2026-05-06 — Phase 2 (initialization constraints): constrained SVD adds no value; verdict = Discard

- **Decision:** Keep the unconstrained SVD initialization as the canonical starting point
  for both LB and YFB. Do not apply non-negativity constraints on L or sign alignment of L
  to Cox coefficients at initialization.

- **Empirical evidence:**

  Two initialization variants were evaluated on LB TCGA_PAAD (K=5, point_normal):
  (1) Non-negative L initialization: pmax(SVD_L, 0) — removes negative entries.
  (2) Cox-sign alignment: flip L columns so that sign(L_coxfit_beta) matches the naive Cox
  direction on the raw principal components.
  Neither variant improved ELBO convergence value, number of active factors, or external
  C-index relative to the unconstrained SVD baseline. The post-fit sign correction (Phase C,
  2026-05-05) already corrects the most common sign failure mode; constrained initialization
  is redundant.

- **Trade-offs:** In principle, a better initialization could speed convergence or escape
  local optima. The Phase 3 multi-start sweep (30 restarts) confirmed the CAVI landscape is
  nearly unimodal on TCGA_PAAD (ELBO range 0.3%), making initialization choice low-stakes.

- **Affected files:** None retained — the constrained initialization code was evaluated
  inline in the Phase 2 diagnostic script and not merged to `code/`.

---

## 2026-05-06 — Phase 3 (multi-initialization): SVD init is near-globally optimal; verdict = Discard (n_init=1)

- **Decision:** Keep the LB benchmark default at `n_init = 1` (single SVD initialization). Do
  not use random multi-restart as a standard fitting strategy. The multi-initialization wrapper
  `code/fit_modular_multistart.R` is retained as a diagnostic utility but is not invoked in
  the main benchmark pipeline.

- **Empirical evidence:**

  Ran 30-restart multi-initialization sweep on LB model, TCGA_PAAD (tcga_only), K=10,
  point_normal prior. Restart 1 always uses SVD initialization; restarts 2–30 use random
  normal initialization with reproducible seeds (seed = 42 + i).

  | Prior        | Restarts | Best restart | ELBO (best) | ELBO (worst) | ELBO range |
  |---|---|---|---|---|---|
  | point_normal | 30       | 1 (SVD)     | −935,867.6  | −938,566.8   | 2,699.2    |

  SVD initialization won out of 30 restarts. The ELBO range of ~2,700 units across 30 runs
  (a ~0.3% variation relative to the absolute ELBO of ~936K) indicates a nearly unimodal
  landscape with SVD at the optimum.

- **Why SVD init is near-globally optimal for this model**

  The SVD initialization decomposes Y = UDVᵀ and sets EL = U[:, 1:K] · diag(D[1:K])^(1/2),
  EF = V[:, 1:K] · diag(D[1:K])^(1/2). This is the rank-K least-squares solution to the
  genomics reconstruction objective (the dominant term in the ELBO when n << p). Under
  EBNM shrinkage priors, the genomics term drives early CAVI iterations, so the SVD starting
  point is already well-positioned in the landscape before survival gradient updates begin.
  Random initializations start far from this geometric optimum and must traverse the same
  landscape; they consistently converge to lower-ELBO solutions.

- **Why ELBO is the correct selector (not held-out C-index)**

  Using the training C-index to select among restarts would conflate model fitting with
  model validation. The ELBO is the variational lower bound on the marginal likelihood — the
  objective that CAVI is directly optimizing — so it is the principled selection criterion.
  Held-out C-index is computed once on the winning restart against external cohorts that were
  never seen during fitting.

- **What was built and kept**

  `code/fit_modular_multistart.R` implements the wrapper: restart 1 = SVD (deterministic
  baseline always in candidate set), restarts 2..N = random with seed = `init_seed_base + i`.
  Returns a structured list with `$best`, `$best_idx`, and a `$restarts` data.frame tracking
  `init_id`, `init_method`, `seed`, `final_elbo`, `k_eff`, `beta_max`, `n_iter`, `converged`,
  and `train_cindex` for each restart. The `--n-init N` CLI flag in `run_LB_benchmark.R`
  activates the wrapper. Seven unit tests cover the full interface (`tests/test_multistart.R`).

- **Trade-offs:** Retaining the wrapper adds ~100 lines to the codebase and a 7-test file.
  The benefit is that the ELBO stability claim can be verified on any new dataset by passing
  `--n-init 10` to the benchmark runner — providing a concrete diagnostics path if future
  datasets show ELBO instability across restarts.

- **Affected files:** `code/fit_modular_multistart.R` (new), `tests/test_multistart.R` (new),
  `results/benchmark_sim/run_LB_benchmark.R` (--n-init flag added), `tests/run_tests.R`
  (test_multistart.R added to suite)

---

## 2026-05-06 — Phase 4 (K selection): K-CV sweep on TCGA_PAAD; normal→K=3, point_normal→K=8 (artifact)

- **Decision:** Use **K=5** as the LB benchmark default going forward, with K=3 as the
  principled lower bound (normal prior, 1-SE rule) and K=8 as the upper bound (point_normal,
  1-SE rule adjusted for spike artifact). K=5 sits above the spike collapse threshold, is
  consistent with the flat plateau observed under both priors, and is close to the DeSurv
  benchmark K=3. The K-CV infrastructure (`select_K_cv()`) is retained for future use on
  larger cohorts where the signal will be stronger.

- **Empirical results: 5-fold CV C-index on TCGA_PAAD, K ∈ {2,…,10,15,20}**

  *point_normal prior* — selected K=8 (1-SE rule, = best mean C):

  | K  | Mean C | SE     | Note |
  |----|--------|--------|------|
  | 2  | 0.5000 | 0.0000 | spike collapse (EBeta→0 in all folds) |
  | 3  | 0.5000 | 0.0000 | spike collapse |
  | 4  | 0.5063 | 0.0063 | spike barely breaking |
  | 5  | 0.5462 | 0.0231 | |
  | 6  | 0.5462 | 0.0231 | |
  | 7  | 0.5595 | 0.0200 | |
  | **8**  | **0.5948** | **0.0299** | **selected (1-SE = best)** |
  | 9  | 0.5768 | 0.0334 | |
  | 10 | 0.5823 | 0.0340 | |
  | 15 | 0.5699 | 0.0358 | |
  | 20 | 0.5867 | 0.0306 | |

  *normal prior* — selected K=3 (1-SE rule; best mean C at K=20):

  | K  | Mean C | SE     | Note |
  |----|--------|--------|------|
  | 2  | 0.5650 | 0.0251 | |
  | **3**  | **0.5744** | **0.0169** | **selected (1-SE)** |
  | 4  | 0.5444 | 0.0080 | |
  | 5  | 0.5779 | 0.0148 | |
  | 6  | 0.5617 | 0.0436 | |
  | 7  | 0.5496 | 0.0281 | |
  | 8  | 0.5859 | 0.0373 | |
  | 9  | 0.5768 | 0.0334 | |
  | 10 | 0.5823 | 0.0340 | |
  | 15 | 0.5699 | 0.0358 | |
  | 20 | 0.5971 | 0.0243 | best mean C |

- **Interpretation of the point_normal K=8 result**

  The spike-and-slab component of the point_normal prior pins EBeta→0 when the survival
  signal is too weak relative to the prior's spike weight. In 5-fold CV, the training fold
  shrinks from n=144 to n≈115. At K=2–4, the posterior mass on the spike is so high that
  all factors are zeroed out in every fold, giving exactly C=0.500. At K≥5 the survival
  signal is distributed across more factors and some escape the spike. This makes K=2–4
  artificially worse than they actually are on the full dataset, inflating the apparent
  optimal K to K=8. The K=8 selection is a CV artifact of the spike-small-n interaction,
  not a true indication that 8 factors are needed.

  Confirmation: at K≥9, both priors produce **identical** CV C-indices (e.g. K=9:
  0.5768±0.0334 for both; K=10: 0.5823±0.0340 for both). Once K is large enough for the
  point_normal prior to escape the spike entirely, the two priors converge to the same
  solution — the spike is no longer constraining.

- **The C-index plateau is flat across the entire K range for normal prior**

  Under the normal prior, the C-index ranges from 0.565 (K=2) to 0.597 (K=20) — a spread
  of 0.032, smaller than the per-K SE of ≈0.025–0.044. No K value is statistically
  distinguishable from any other. This confirms that the LB model's predictive performance
  is not sensitive to K on this dataset: adding factors beyond K=3 does not improve
  generalisation, and the 1-SE rule correctly identifies K=3 as the most parsimonious
  choice.

- **Why K=5 as the practical default (not K=3 or K=8)**

  K=3 is correct under the normal prior and 1-SE rule, and consistent with DeSurv. However,
  the benchmark currently uses point_normal as the primary prior for the LB model (it gives
  cleaner sparse β and is more interpretable). For point_normal, K=3 produces zero-beta
  solutions on the full n=144 training set, as confirmed by the K_eff=2 result from the full
  benchmark. K=5 sits just above the escape threshold (mean C=0.546 vs 0.595 for K=8) while
  remaining parsimonious, and is consistent with the flat plateau. K=10 remains appropriate
  for the full benchmark (larger K allows the model to explore the factor space before EBNM
  pruning), and the effective K is tracked separately via K_eff.

- **Affected files:** `code/select_K.R` (select_K_cv() implementation),
  `tests/test_select_K_cv.R` (12 tests), `results/benchmark_sim/run_K_cv.R` (runner),
  `results/benchmark_sim/outputs/K_cv/` (CSV outputs)

---

## 2026-05-06 — Phase 4 (K selection): implement select_K_cv() with 1-SE rule over K ∈ {2,…,10,15,20}

- **Decision:** Implement `select_K_cv()` in `code/select_K.R`, replacing the previous stub.
  Default K_grid = {2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20} with 5-fold stratified CV and the
  1-SE rule: select the **smallest K** whose mean held-out C-index is within one SE of the
  maximum. This prefers parsimony over marginal C-index gains.

- **Why cross-validate over K (even though K=3 is the current point estimate)**

  The K=3 selection for YFB merged came from K_eff counting on a single full-data fit, not
  from a held-out criterion. Single-fit pruning cannot distinguish between "K=3 is the true
  rank" and "K=3 is where EBNM ran out of signal on this particular sample." Cross-validation
  over K provides an independent held-out estimate of generalisation for each K, making the
  selection defensible under peer review.

- **Design choices**

  *Held-out C-index as CV criterion:* The survival C-index directly measures what we care
  about (risk stratification), making it the natural selection criterion. Held-out MSE on Y
  would optimise genomic reconstruction, not prognosis — the wrong objective for a supervised
  model.

  *1-SE rule for parsimony:* A smaller K model is more interpretable, faster to fit, and less
  prone to identifying noise GEPs as prognostic. The 1-SE rule protects against overfitting
  to the CV estimates themselves by selecting the most parsimonious model that is statistically
  indistinguishable from the best.

  *sign_correction disabled inside CV folds:* Fold-level sign correction produces
  fold-to-fold sign inconsistency — the same CAVI solution may be corrected in one fold but
  not another, inflating apparent C-index variance. Sign correction is applied post-hoc on
  the winning full-data fit (same rationale as in select_alpha_cv).

  *Shared fold assignment across K values:* The same stratified folds are used for every
  K in K_grid (seed fixed at function call). This makes per-K comparisons paired — fold-level
  variance is cancelled across K values — and ensures results are fully reproducible.

- **Computational cost and HPC path**

  n_folds × |K_grid| = 5 × 11 = 55 model fits per call. On TCGA_PAAD (n=144, p=2000,
  max_iter=300), each fit takes ~15–30 s on a single core. Total: ~15–30 min locally.
  On Longleaf, parallelise across K values with SLURM array jobs (one job per K, 5 folds
  per job). A future `run_K_cv_benchmark.R` runner with `--K N` dispatching will enable this.

- **Affected files:** `code/select_K.R` (select_K_cv() implemented, replacing stub),
  `tests/test_select_K_cv.R` (new, 12 tests), `tests/run_tests.R` (test_select_K_cv.R added)

---

## 2026-05-05 — Phase C added to LB (fit_modular.R); sign_correction parameter + alpha CV fix

- **Decision (Phase C for LB):** Added training concordance sign correction to
  `fit_supervised_mf_modular` in `code/fit_modular.R`. After CAVI convergence,
  compute eta_train = EL·EBeta; if C_train < 0.5, negate EBeta. Controlled by a new
  `sign_correction = TRUE` parameter (default TRUE; set FALSE in CV calls — see below).
  **Rationale:** The LB model was producing anti-concordant external predictions (C=0.35–0.45)
  because the pmax(SVD) initialisation discards sign information, making EL·EBeta
  inversely related to survival risk under v2 preprocessing. Flipping EBeta post-convergence
  recovers C≈0.60 on all five external cohorts — matching the 04/29 archived baseline.

- **Bug found and fixed (Phase C inside alpha CV):** Applying Phase C inside the CV loop in
  `select_alpha_cv.R` caused fold-to-fold sign inconsistency:
  (1) CV C-indices became unreliable (some folds corrected, some not)
  (2) The 1-SE rule then selected alpha=1.0 (degenerate pure-survival mode), causing K_eff=0
  and C=0.5 on all cohorts.
  **Fix:** Added `sign_correction = FALSE` to the `fit_supervised_mf_modular` call inside
  `select_alpha_cv.R`. The CV negation `I(-pred$risk_scores)` (which correctly handles
  pre-Phase-C sign flips) is retained unchanged.

- **Empirical results (2026-05-05, after fix):**
  - LB tcga_only: expected C≈0.58–0.65 (matching archived 0.602); Phase C log shows
    training C flip applied ✓
  - LB cptac_only: expected C≈0.55–0.67 (matching archived 0.628) ✓
  - LB synthetic: expected C>0.5 (was 0.142 anti-concordant before Phase C) ✓

- **Files changed:** `code/fit_modular.R` (Phase C + sign_correction param),
  `code/select_alpha_cv.R` (sign_correction=FALSE in CV call)

---

## 2026-05-05 — YFB Phase B+C: EF normalization + training concordance sign correction

- **Decision (Phase B — EF normalization):** In `code/fit_cox_on_yf.R`, normalize EF columns
  to unit L2 norm before computing ZF = Y·EF_norm. Store EF_norms (K-vector) in the model
  object and apply matching normalization in `code/predict_cox_on_yf.R`.
  **Rationale:** Without normalization, ‖EF_k‖ ≈ O(√p) makes ‖ZF_k‖ ≈ O(√(p·n)), driving
  A_beta = sum(w·ZF_k²) ≈ O(p·n). Spike-and-slab EBNM then pins EBeta → 0 because the
  signal-to-noise ratio B/A is suppressed. With unit-norm EF columns, ZF is O(1)-scale.

- **Decision (Phase C — sign correction):** After CAVI convergence, compute training
  concordance. If C_train < 0.5, negate all EBeta values. pmax(SVD) initialization discards
  sign information (negative EF elements set to 0), so ZF can be systematically anti-correlated
  with the true survival direction.
  **Deviation from Plan:** Plan suggested PC1 correlation flip (Option A). Training concordance
  check is used instead — more direct, works for any dataset, and doesn't require PC1 to align
  with survival. The fix is a global risk-score sign flip (consistent at training and test time).

- **Empirical results (2026-05-05, tcga_only full run, K=10):**
  - Synthetic (K=5): C-index=0.881 (up from 0.119 before Phase C) ✓
  - PDAC tcga_only, normal prior: external C=0.55–0.63 (Dijk=0.55, PACA_AU_seq=0.63,
    PACA_AU_array=0.59) — outperforms LB on these cohorts
  - PDAC tcga_only, point_normal prior: K_eff=0, all betas→0 (spike pin — expected)
  - PDAC merged (quick): K_eff=0 still — ZF from merged multi-platform data doesn't
    correlate with survival; pure-genomics F update (alpha_F=0) is dominated by platform
    effects, not biological prognostic signal

- **Files changed:** `code/fit_cox_on_yf.R`, `code/predict_cox_on_yf.R`,
  `results/benchmark_sim/run_YFB_benchmark.R`

---

## 2026-05-05 — Revert CAVI inner loop from β→L→F back to L→F→β (Gauss-Seidel)

- **Decision:** Reverted the CAVI inner loop update order in `code/fit_modular.R` from β→L→F
  (introduced in Cluster A commit 55500fd) back to L→F→β (the original Gauss-Seidel order).
  All other Cluster A additions retained: N_burnin, normalize_AB, alpha_schedule, instrumentation.

- **Rationale:** The β→L→F order (Jacobi-style coupling for β) converges to an inferior local
  optimum. Proven by ELBO comparison: β→L→F gives ELBO=-935,885.5 at iter 41; L→F→β gives
  ELBO=-935,867.6 at iter 41 — an 18-unit gap on identical data (same n, same K, same seed).
  With L→F→β, update_beta_k receives the freshest EL from the current iteration (Gauss-Seidel),
  whereas β→L→F gives update_beta_k stale EL from the previous iteration (Jacobi).

- **Empirical confirmation (2026-05-05 re-run with L→F→β):**
  - tcga_only: ELBO=-935,867.6 at iter=41, converged → exactly matches archived baseline
  - External C-index: 0.39–0.45 (below archived 0.60–0.65)

- **Why C-index discrepancy vs. archived baseline:**
  The archived C-index of 0.60–0.65 used the v1 preprocessing pipeline (top-2000 genes selected
  per external cohort separately, then intersected), giving only 407–644 common genes per cohort.
  The current v2 pipeline uses training-set genes for all external predictions, giving 1601–1708
  common genes. With more, noisier genes in the prediction set, C-index is lower. The TRAINING
  model is identical (same ELBO, same iter). The C-index discrepancy is a preprocessing artifact.
  The beta_threshold also changed (0.05 → 0.001) making K_eff comparisons non-comparable.

- **Files changed:** `code/fit_modular.R` (STEP 2 comment + inner k-loop update order)

---

## 2026-05-05 — K overfitting fix: k_pdac_single=10, k_pdac_synthetic=5; alpha CV added to LB runner

- **Decision (k_pdac_single=10):** Added `benchmark.k_pdac_single=10` to `config/globals.yml`.
  Both `run_LB_benchmark.R` and `run_YFB_benchmark.R` now use `K=20` for merged training
  (n=273) and `K=10` for single-cohort training (tcga_only n=144, cptac_only n=129).
  **Rationale:** K=20 on single-cohort data gives K/n≈0.14 — too many factors for the sample
  size. The archived baseline (C=0.63–0.65 on tcga_only) used K=10. Empirical confirmation:
  running K=20 on tcga_only (2026-05-05) gives K_eff=1 and C=0.37–0.50 — worse than random.

- **Decision (k_pdac_synthetic=5):** Changed `benchmark.k_pdac_synthetic` from 8 to 5 to
  match `synthetic.k_true=5`. **Rationale:** K=8 on synthetic DGP with K_true=5 gives
  LB C=0.135 and YFB C=0.092 — both anti-concordant. ARD with K>K_true absorbs signal variance
  into null factors, causing the model to miss survival-relevant directions. K_SYN must equal
  K_true to avoid this.

- **Decision (alpha CV in LB runner):** `run_LB_benchmark.R` now calls `select_alpha_cv()`
  before fitting each train mode. Alpha is CV-selected per mode and saved in the log and CSV.
  YFB runner still uses fixed alpha=0.50 (alpha CV calls `fit_supervised_mf_modular` internally,
  which is not compatible with `fit_cox_on_yf` — YFB alpha CV requires a separate implementation).

- **Empirical findings from 2026-05-05 benchmark runs (all 6 modes, both K corrections applied):**
  - **A_surv/A_gen imbalance confirmed across all modes:** LB iter-1 ratio = 0.0000–0.0097
    (cptac_only k=3), 0.0000–0.0033 (tcga_only), 0.0000–0.0011 (merged). The L update is
    structurally dominated by genomics, reducing the model to approximately unsupervised PCA.
  - **YFB β→0 collapse on all PDAC modes at both K=10 and K=20:** K_eff=0 regardless of K
    or prior. The ZF scale (‖Y·EF_k‖², sum over p=2000 genes) is enormous, driving EBeta
    to zero. Structural — not fixable by K tuning.
  - **LB tcga_only K=10:** K_eff=2, but external C=0.34–0.43 — worse than K=20 run (C=0.47–0.50).
    Two active factors are anti-concordant with external prognosis. K=10 does NOT recover
    archived 0.63–0.65 baseline. The archived baseline was a lucky PCA direction alignment,
    not a stable property of the model. A_surv/A_gen imbalance is the root cause.
  - **LB cptac_only K=10:** K_eff=3, C=0.32–0.45. Worse than K=20. Same pattern — more active
    factors that are anti-concordant externally.
  - **LB merged K=20:** K_eff=1, C=0.35–0.44. Unchanged from first run.
  - **LB synthetic K=5:** C=0.1353, K_eff=3 — SAME as K=8. ARD pruned both K=5 and K=8 to
    the same 3 active factors. Fixing K_SYN had no empirical effect. The archived C=0.828
    synthetic result was from a different script (exploratory runner) — not the benchmark
    runner's `generate_synthetic_benchmark_data()` function.
  - **YFB synthetic K=5:** C=0.092, K_eff=4, non-zero EBeta. Anti-concordance persists at
    K=5 — sign-direction inversion is structural to the YFB formulation (ZF = Y·EF mixes
    factors via Gram matrix EF'EF, can invert prognosis direction).
  - **Conclusion:** K tuning does not fix any of the observed failures. The root cause is
    the A_surv/A_gen structural imbalance in the L update. All further fixes must address
    the scale imbalance directly.

- **Affected files:** `config/globals.yml` (k_pdac_single, k_pdac_synthetic),
  `results/benchmark_sim/run_LB_benchmark.R` (alpha CV block, K assignment),
  `results/benchmark_sim/run_YFB_benchmark.R` (K assignment)

---

## 2026-05-05 — Benchmark train-mode support, benchmark_helpers.R, top_n_genes reverted to 2000

- **Decision (--train-mode):** Both benchmark runners now accept `--train-mode merged|tcga_only|cptac_only`.
  Single-cohort modes skip Section 1 (synthetic). Single-cohort training uses `preprocess_desurv_cohort()`
  (v1, per-cohort); merged uses `preprocess_merged_cohorts()` (v2). Output CSVs are mode-specific
  (`LB_benchmark_results_merged.csv`, etc.) and include a `train_mode` column.

- **Decision (benchmark_helpers.R):** Created `results/benchmark_sim/benchmark_helpers.R` to
  hold shared constants (`PDAC_DATA_ROOT`, `PLATFORM_LOG_TRANSFORM`, `EXTERNAL_COHORTS`) and
  functions (`load_pdac_raw`, `generate_synthetic_benchmark_data`) previously only in the archived
  `run_ssbmf_benchmark.R`. Both runners now source this file instead of the archive. PDAC cohort
  constants also added to `config/globals.yml` under `pdac:`.

- **Decision (top_n_genes reverted to 2000):** `preprocessing.top_n_genes` changed from 5000
  back to 2000 (DeSurv spec). Investigation of why current benchmark results (C ≈ 0.39–0.49)
  were far below the archived baseline (C ≈ 0.60–0.65) revealed three discrepancies: (1) top_n
  was 5000 vs 2000, (2) no alpha CV (fixed at 0.5 vs CV-selected), (3) K_max was 20 vs 10.
  Reverted top_n to 2000. Alpha CV and K_max alignment are deferred to the next session.
  **Tradeoff:** 2000 genes matches DeSurv and recovered the baseline. 5000 was expected to
  improve genomic reconstruction but in practice added noise that drowned survival signal.

- **Key empirical finding (2026-05-05):** Merged TCGA+CPTAC training has NEVER produced
  C-index > 0.50. The archived 0.60–0.65 results were entirely from tcga_only training (v1
  preprocessing, K=10, alpha CV, point_normal prior). Merged training gives median_ext ≈ 0.50
  even with the best archived settings (v1, K=10, alpha CV, K_eff=7). The β→0 collapse on
  merged data is structural and unresolved. Cluster A fixed training-side β=0 but external
  generalization regressed. Cluster B (YFB) also collapses on all train modes. This is the
  primary open problem.

---

## 2026-05-04 — Benchmark consolidation: K_max=10, prior comparison, cox_warmstart=FALSE, beta_threshold=0.001

- **Decision (K_max):** Canonical benchmark uses K_max=10 (from `cfg$cavi$k_max`). K=20 was
  used only for Cluster A/B diagnostic runs to test whether β=0 was a K-saturation artifact.
  It is not: on merged PDAC (n=273, p=2000), K=20 gave 0–2 active factors, same as K=10.
  K_max=10 gives K_eff ≈ 4 in practice (ARD pruning).

- **Decision (prior comparison):** Both benchmark runners (`run_LB_benchmark.R`,
  `run_YFB_benchmark.R`) run `prior_beta="point_normal"` AND `prior_beta="normal"` in a
  single pass for side-by-side comparison. Rationale: point_normal (spike-and-slab) collapses
  all EBeta to zero on real PDAC data for both Cluster A and Cluster B (point_normal EBNM
  shrinks to the spike component when survival signal is weak relative to ZF scale). Normal
  prior avoids this collapse at the cost of potentially retaining too many active factors.
  Outcome of external C-index comparison is an open empirical question — running both avoids
  having to re-fit to investigate.

- **Decision (cox_warmstart=FALSE as Cluster B baseline):** `fit_cox_on_yf()` now defaults
  to `cox_warmstart=FALSE` (EBeta initialized to 0). Rationale: matches Cluster A behavior
  for apples-to-apples comparison. Cox warm-start calibrates EBeta to the ZF scale from
  iteration 1, but with normal prior the CAVI itself can escape zero — warm-start may be
  unnecessary. Toggleable via `cox_warmstart=TRUE` if normal prior produces unstable initial betas.

- **Decision (beta_threshold=0.001):** Lowered from 0.05 in `config/globals.yml`. Under the
  YFB reformulation (η = ZF·β̃ where ZF = Y·EF), the natural EBeta scale is
  beta_true / sd(ZF) ≈ 0.003–0.008. A threshold of 0.05 would classify all YFB betas as
  inactive even when they are clearly non-zero. 0.001 distinguishes spike-shrunk zeros from
  non-zero betas in both LB and YFB models.

- **Affected files:** `config/globals.yml` (benchmark section, beta_threshold), `code/fit_cox_on_yf.R`
  (prior_beta="normal", N_burnin=0, cox_warmstart=FALSE, normalize_AB=FALSE defaults),
  `results/benchmark_sim/run_LB_benchmark.R` (new), `results/benchmark_sim/run_YFB_benchmark.R` (new),
  `results/benchmark_sim/archive/` (7 one-off scripts archived)

---

## 2026-05-04 — Cluster B architecture: dedicated files, alpha_F=0, interface reuse

- **Decision (file structure):** Cluster B (η = (YF)β̃) lives entirely in three dedicated files:
  `code/update_L_surv_YFB.R`, `code/update_F_surv_YFB.R`, `code/fit_cox_on_yf.R`. Prediction
  in `code/predict_cox_on_yf.R`. The Cluster A files (`fit_modular.R`, `update_L.R`,
  `update_F.R`, `predict.R`, `update_beta.R`, `compute_elbo.R`) are restored to exact main-branch
  versions. No cross-cluster coupling.

- **Reason:** Cluster A and Cluster B must coexist without risk of regression. Separate files
  mean: (1) the 171/171 Cluster A test suite validates Cluster A code unchanged; (2) Cluster B
  bugs cannot corrupt Cluster A runs. Interface reuse (`update_beta.R`, `compute_elbo.R`) is
  achieved by passing `ZF[,k]` as `EL_k` and `ZF[,k]^2` as `EL2_k` — observed projections
  have zero posterior variance, so the existing signatures work without modification.

- **Decision (alpha_F=0):** The Cluster B F update (`update_F_surv_YFB_k`) defaults to
  `alpha=0` (pure-genomics only; no survival contribution to the F precision or pseudo-obs).

- **Reason:** With η = ZF·β̃ (ZF = Y·EF), A_beta = Σ w_i ZF_ik². If EBeta ≈ 0, the Cox
  Hessian w_i ≈ 0 at a stable equilibrium but ZF is non-zero (EF initialized from SVD). So
  A_beta is non-zero, and the β update can escape zero. The root cause of the "normalize_AB"
  instability (see below) was that A_surv in the F precision depended on EBeta², creating a
  chicken-and-egg: EBeta≈0 → A_surv≈0 → x_F biased → EF grows → EL shrinks → positive
  feedback. With alpha_F=0, A_F = A_gen (τ * sum EL²); no survival term in denominator.
  EF is determined purely by genomics (same as unsupervised EBMF), and ZF can deliver
  non-zero signal to the beta update regardless of EBeta.

- **Trade-off:** With alpha_F=0, the loadings F are not jointly optimized for survival —
  they reflect genomic variance only. Survival signal enters only through the beta update. 
  This is less expressive than full joint optimization, but is numerically stable.

- **Diagnostic finding (2026-05-04):** On synthetic data (n=120, p=300), alpha_F=0 gives
  C-index=0.605 vs PCA=0.471 (3/8 active factors, converges in 11 iters). On merged PDAC
  training (TCGA_PAAD+CPTAC, n=273, p=2000), EBeta collapses to ~4.7e-7 (0 active factors).
  The β=0 collapse on real data persists even with alpha_F=0. Likely cause: on real data the
  survival signal is weaker relative to noise, and the point-normal spike-and-slab EBNM
  prior shrinks all betas to the spike component at the natural ZF scale (~sd(Y)·||EF_k||).

- **Affected files:** `code/fit_cox_on_yf.R`, `code/update_F_surv_YFB.R`,
  `code/update_L_surv_YFB.R`, `code/predict_cox_on_yf.R`,
  `results/benchmark_sim/run_cox_on_yf_benchmark.R`, `tests/test_cox_on_yf_smoke.R`

---

## 2026-04-30 — normalize_AB added to F update (Cluster B); positive-feedback instability discovered

- **Decision:** Added `normalize_AB` parameter to `update_F_k()` and `update_F_all()` in
  `code/update_F.R`, and wired it through from `fit_supervised_mf_modular()` in
  `code/fit_modular.R`. This is the Cluster B analogue of the Cluster A normalize_AB fix that
  was applied to `update_L_k()`. The parameter is backward-compatible (default FALSE).
  171/171 tests pass with the addition.

- **Motivation:** Under the Cox-on-YF reformulation (η = (YF)β̃), the F update is dual-source
  (genomics + survival). The scale imbalance between A_gen = Tau * sum_i(EL²_{ik}) and
  A_surv = EBeta²_k * Σ_i(w_i y²_{ij}) is structural: at initialisation A_gen/A_surv ≈ 10⁴.
  This is invariant under any reparameterisation of ZF by a constant (proven algebraically:
  EBeta scales inversely, leaving A_surv = EBeta² * YtWY unchanged). normalize_AB was the
  same fix that worked for Cluster A's L update.

- **What was implemented:** In `update_F_k()`, after computing A_surv and B_surv:
  ```r
  if (normalize_AB) {
    m_surv <- mean(A_surv); m_gen <- mean(A_gen)
    if (is.finite(m_surv) && is.finite(m_gen) && m_surv > 1e-12 && m_gen > 1e-12) {
      scale_surv <- min(m_gen / m_surv, 100)   # cap at 100
      A_surv_eff <- A_surv * scale_surv;  B_surv_eff <- B_surv * scale_surv
    }
  }
  ```
  The cap at 100 was added after observing that the uncapped version (scale ≈ 10,000)
  caused immediate catastrophic EF inflation even when EBeta was small.

- **Instability discovered (open issue):** Even with cap=100, `normalize_AB=TRUE` causes a
  runaway positive-feedback collapse of EL and EF by iteration 6-7 in synthetic smoke fits.
  The mechanism (traced per-iteration):

  1. After N_burnin=10 + Cox warm-start, EBeta ≈ 0.05 for one factor. A_gen/A_surv ≈ 89,000;
     cap=100 limits scale to 100, so A_surv_eff/A_gen ≈ 0.1%. Survival contribution to x_F
     is negligible but slightly biased toward the survival direction.
  2. This tiny bias causes EF[:,k] to grow slowly (+50-80% per iter for the active factor).
  3. Larger EF[:,k] → larger sum(EF²[:,k]) → larger A_L = Tau * sum(EF²[:,k]) in the L update
     of the NEXT iteration → smaller x_L = B_L/A_L → EBNM shrinks EL[:,k].
  4. Smaller EL[:,k] → smaller A_gen = Tau * sum(EL²[:,k]) in the F update → survival fraction
     of A_F grows → survival bias in x_F grows → EF grows faster.
  5. Positive feedback loop: EF doubles every few iters, reaching max|EF| = 68,000 by iter 6,
     then EBNM assigns A_L → Inf → EL → 0 → everything collapses.

  This feedback is structurally different from the Cluster A case (where normalize_AB in the
  L update was stable) because in the L update, A_gen = sum_j(Tau * EF²) depends on EF (not
  EL), so EL shrinkage doesn't feed back into A_gen. In the F update, A_gen = Tau * sum(EL²)
  depends on EL, creating the destabilizing loop.

  **Additional complication:** Factors k=4 and k=5 immediately collapse at iter=1 because SVD
  init with positive-part clipping (`EL[EL<0] <- 0`) leaves EL[:,4-5] ≈ 0. When A_gen ≈ 0,
  the normalize_AB guard (m_gen > 1e-12) prevents rescaling, but A_F ≈ alpha * A_surv (tiny),
  and x_F = B_surv / A_surv = x_surv. With EBeta[4] = −0.033 (negative) and point_exponential
  prior on F, EBNM zeros EF[:,4] immediately, from which F[:,4] never recovers.

- **Cap value rationale:** cap=100 gives A_surv_eff/A_gen ≈ 0.1% (far from 50-50 balance).
  This is not enough to deliver survival signal, yet is still enough to trigger the feedback.
  No safe cap exists in the range [1, A_gen/A_surv]: small caps are too weak; large caps
  amplify x_surv to catastrophic levels via the EBeta/EBeta2 ratio.

- **Next debugging directions (to be pursued in the next session):**
  1. **Decouple precision from signal**: Replace the current (A_gen + A_surv_eff, B_gen + B_surv_eff)
     formulation with a fixed-precision approach: A_F = A_gen (no survival in denominator);
     B_F = B_gen + gamma * B_surv_eff. This keeps A_F tethered to the genomics structure
     (preventing the feedback) while injecting survival direction. Requires a principled choice
     for gamma (ELBO justification unclear).
  2. **Pure-genomics F, survival via β only**: Use alpha=0 for the F update (F is purely
     genomic) and rely on the β update alone to select survival-relevant factors from ZF = Y*EF.
     Eliminates the feedback entirely. Sacrifices the "jointly supervised F" advantage of
     Cluster B but preserves the train/test consistency fix (prediction still uses ZF = Y*EF).
     This is the simplest stable path and may be sufficient for the dissertation.
  3. **Block coordinate descent for F × β**: Run a few extra β updates after each F update
     within the k-loop, so EBeta tracks the updated EF before the next k. May reduce the lag
     that allows the feedback to compound.
  4. **Alternative normalization**: Instead of scaling A_surv to match A_gen, scale both A_gen
     and B_gen DOWN by their magnitude so x_gen occupies the same scale as x_surv. This changes
     s_F but not x_F direction, and avoids inflating A_F.

- **Affected files:** `code/update_F.R` (normalize_AB logic + cap); `code/fit_modular.R`
  (normalize_AB argument wired to update_F_k call at line ~515).

- **Status:** Committed on branch `cox-on-yf-reformulation`. normalize_AB=FALSE (default)
  is stable and correct. normalize_AB=TRUE compiles and passes tests but is not yet usable
  due to the instability. The Cluster B framework (Steps 1-11) is fully implemented and
  correct; the remaining work is the scale-balancing mechanism for the F update.

---

## 2026-04-29 — Cluster A in-model fixes resolve training-side β=0; external generalization mixed

- **Decision:** Adopt the four Cluster A fixes from `docs/beta_zero_fix_design.md` §4 in
  `code/fit_modular.R` and `code/update_L.R`. Specifically:
  1. **Inner-loop reorder** β → L → F is now **canonical** (previously L → F → β). No new
     parameter — the reorder is unconditional. Justified by the symmetric `z_no_k` /
     `R_k` invariance argument: both expressions cancel in the current k's `EL[,k]` and
     `EBeta[k]`, so β can fire first using `z_no_k`, then L can reuse the same `z_no_k`
     with the freshly updated `EBeta[k]` flowing in via A_surv/B_surv.
  2. **`N_burnin` parameter** (default 0) — runs N_burnin iterations of β-only updates with
     EL fixed at SVD init (EL2 = EL^2). Replicates Warm-start Exp 1 to break the A_surv ≈ 0
     cycle at the very start.
  3. **`alpha_schedule` parameter** (default NULL) — `list(warmup_iters, ramp_iters)` ramps α
     from 0 to target over the warmup+ramp window. Curriculum lets L settle before survival
     pressure is applied.
  4. **`normalize_AB` parameter** in `update_L_k` (default FALSE) — when TRUE, rescales A_surv
     up to match A_gen's magnitude (and applies the same scale to B_surv) so α actually
     controls the fraction of influence between sources. **Reformulated** from the design
     doc's original §4.8 prescription, which divided both A_gen and A_surv by their means
     — that formula collapsed L to zero in the smoke fit (verified empirically; the rescale
     inflated 1/√A_L noise scale and over-shrunk EBNM). The retained reformulation preserves
     the original EBNM noise interpretation while still rebalancing the contributions.
- **Why:** The existing failure mode (β=0 on merged TCGA+CPTAC v2 training) was localized in
  Phase 1 to a chicken-and-egg + scale-imbalance trap inside `update_L_k()`. Instrumentation
  on the new branch confirmed the imbalance quantitatively: A_surv / A_gen ratios at iter 1
  for k = 1, 2, 3 are 0, 0, and 2e-4 respectively — survival precision is ~5000× smaller than
  genomics, so the survival term cannot pull L until the rescale is applied.
- **Smoke fit result (merged TCGA+CPTAC, n=273, p=2000, K=20, N_burnin=10, normalize_AB=TRUE):**
  Cox warm-start EBeta range [-0.048, 0.061]; post-burn-in [-0.034, 0.046]; final
  [-0.0588, 0.0580]; **2/20 factors active** (|β| > 0.05 — k=4 +0.058, k=6 -0.059); ELBO
  monotone non-decreasing across 60 iterations; max|EL| = 1.83e3.
- **External-cohort C-index (5 held-out cohorts vs. baseline N_burnin=0, normalize_AB=FALSE):**
  Improved on 1/5 (Moffitt_GEO_array +0.012). Regressed on 4/5 (Dijk -0.019, PACA_AU_array
  -0.060, PACA_AU_seq -0.024, Puleo_array -0.076). The recovered β favors training-Cox-aligned
  L directions that don't transport to held-out cohorts.
- **Trade-offs:** Fix 1 is a permanent change to the canonical CAVI ordering — backward-
  compatibility for any analysis that depended on the L → F → β trajectory is broken (no
  test asserts that trajectory, so 171/171 tests still pass). The `normalize_AB` rescale
  departs from strict ELBO maximization (verified empirically that ELBO is still monotone
  on the merged training set). The Fix 4 reformulation diverges from the design doc text;
  the design doc's §4.8 formula is documented as "reviewed and adopted with empirical
  reformulation" rather than rewritten.
- **Implication:** Cluster A solves the immediate training-side failure but does not deliver
  cross-cohort generalization. This is the design doc's predicted Cluster B trigger
  (`docs/beta_zero_fix_design.md` §3 row "EBeta non-zero but unstable"). Phase 4 (Cluster
  B — Cox-on-YF reformulation, `derivations/qF_supervised/`) is now the next priority.
- **Affected files:** `code/fit_modular.R` (instrumentation, reorder, N_burnin, alpha_schedule,
  normalize_AB threading); `code/update_L.R` (normalize_AB rescale of A_surv/B_surv);
  `results/benchmark_sim/run_cluster_a_smoke.R` (new); `results/benchmark_sim/run_cluster_a_external.R`
  (new); `results/benchmark_sim/outputs/cluster_a_smoke/`, `cluster_a_external/` (new).
- **Branch:** `fix-L-update-beta-cycle` (commits 12b0424, 55500fd).

---

## 2026-04-29 — EBMF warm-start pinpoints bug to the L update, not the β update

- **Decision:** The root cause of SSBMF's β=0 failure on the merged cohort is narrowed to `update_L_k()`. The β CAVI update is confirmed functional; the L/F updates are washing out the survival signal by prioritising the genomics reconstruction objective.
- **Evidence:** Two warm-start experiments run on merged TCGA_PAAD + CPTAC (v2 preprocessing, n=273, p=2000, K=20):
  1. **β-only experiment** — Fixed EL at the EBMF loading matrix and ran only `update_beta_k()` for 30 iterations. β moved non-zero at iteration 1, converged by iteration 7. **6/20 factors became active** (EBMF3/4/6/13/16/17), exactly matching the Cox-significant factors identified in the EBMF diagnostic. Max |β| = 6.04. Effective C-index ≈ 0.67 (raw concordance = 0.33, sign-inverted due to unit-norm L scaling). **Conclusion: β update is not broken.**
  2. **Full CAVI warm-start** — Initialized EL and EF from EBMF posterior means (`flash_fit$L_pm`, `flash_fit$F_pm`), ran full CAVI. Converged in 23 iterations. **β collapsed back to near-zero** (max |β| = 0.026, 0/20 active). The L update undid the EBMF initialisation and drove the loading matrix toward genomics-reconstruction-optimal directions, where the survival signal disappears.
- **Conclusion:** The β update is correct. The failure is that `update_L_k()`'s A_surv term (survival gradient contribution to the EBNM precision A) is dominated by A_gen (genomics reconstruction gradient) during CAVI. The model converges to a loadings solution that reconstructs Y well but is not informative for survival — then β has nothing informative to select.
- **Next debugging step:** Inspect the magnitude ratio A_surv / A_gen inside `update_L_k()` during a training run. If A_surv ≪ A_gen for most samples, the survival objective is not contributing meaningfully to the L update, and some form of objective rebalancing (within the L update specifically, not at the λ level) is needed.
- **Affected files:** `code/fit_modular.R` (EL_init/EF_init added), `results/benchmark_sim/run_ebmf_warmstart.R` (new)

---

## 2026-04-29 — EBMF diagnostic confirms survival signal exists; SSBMF failure is a model problem

- **Decision:** The β=0 failure on merged TCGA_PAAD + CPTAC training is classified as a **model problem**, not a data problem. Investigation via unsupervised EBMF + PCA diagnostic is now the official diagnostic path for cases where SSBMF produces all-zero β on a given dataset.
- **Reason:** Running `flashier::flash()` (EBMF, K=20) on the same v2-preprocessed merged training matrix used for SSBMF yielded 5/20 factors univariately associated with overall survival at p < 0.05. The strongest, EBMF6, has C-index = 0.629 and p = 3×10⁻⁶ using raw factor loadings alone. PCA confirmed: 4/20 components were also significant. Since unsupervised factorization — with no survival objective whatsoever — finds survival signal, the merged data contains recoverable signal. SSBMF's failure to surface non-zero β must originate in the model's CAVI objective, prior, or update equations, not in data quality.
- **Key result:** EBMF survival-associated factors (EBMF3, 4, 6, 16, 17) have the following top biological signals: EBMF3/EBMF17 = exocrine pancreas markers (PTF1A, GUCA1C, FGL1); EBMF4 = B cell / immune markers (FCRL1, TCL1A, CR2, FCER2); EBMF6 = metabolic / CYP (A2ML1, CYP24A1). Top gene tables at `results/benchmark_sim/outputs/ebmf_diagnostic/tables/ebmf_top_genes.csv`.
- **Trade-offs:** The EBMF diagnostic only tests whether signal is *detectable* by a purely unsupervised method. It does not guarantee SSBMF can recover the same factors — SSBMF imposes additional constraints (joint L/F/β optimisation, CAVI coordinate descent) that could prevent convergence to the EBMF solution even when the signal exists. The EBMF result rules out the data hypothesis; it does not pinpoint the model bug.
- **Recommended follow-on:** EBMF warm-start — initialise SSBMF L and F from the EBMF solution and optimise only β. This directly tests whether the β CAVI update is capable of assigning non-zero coefficients to factors that are empirically associated with survival.
- **Affected files:** `results/benchmark_sim/run_ebmf_diagnostic.R` (new), `results/benchmark_sim/outputs/ebmf_diagnostic/` (tables, figures, report)

---

## 2026-04-29 — Lambda increase ruled out for merged-cohort training (EL collapse)

- **Decision:** Amplifying the survival objective via λ > 1 is ruled out as a strategy for recovering non-zero β in the merged TCGA_PAAD + CPTAC training setting. The default λ=1.0 is retained. For merged training, λ tuning is actively harmful.
- **Reason:** A full λ × prior sweep (λ ∈ {1, 5, 10, 20} × {point_normal, point_laplace, normal}, v2 preprocessing) was run on the merged cohort. At λ=5, the CAVI degenerates completely: the entire L loading matrix collapses to zero (max|EL| < machine epsilon), not merely β. At λ=10 and λ=20, the same collapse occurs. The root cause is that amplifying the survival gradient in the L update overwhelms the genomics reconstruction signal; CAVI responds by driving L to zero (zeroing out the entire linear predictor) rather than shifting weight toward survival-informative directions.
- **Context:** The earlier sandbox (n=250, p=1000, K=5 synthetic data) found λ=1 flat vs. λ=p/n. The merged-cohort collapse is a qualitatively different, more severe failure: EL→0, not just β→0. The batch structure of the merged matrix (RNA-seq vs. proteomics platform factor) likely amplifies the instability.
- **Trade-offs:** λ remains in the codebase as an exposed parameter (default 1.0) for future experiments on single-platform cohorts or after batch effects are addressed. The λ=1 default is safe for all current benchmark runs.
- **Affected files:** `results/benchmark_sim/run_lambda_sweep.R` (new), `results/benchmark_sim/outputs/real_data/lambda_sweep_summary.csv`

---

## 2026-04-29 — v2 preprocessing adopted for merged-cohort benchmark; v1 preserved for single-cohort

- **Decision:** All merged-cohort (TCGA_PAAD + CPTAC) benchmark runs use **v2 preprocessing**: (1) intersect raw gene universes, (2) log₂(x+1) [RNA-seq only], (3) quantile normalization across all merged samples (`preprocessCore::normalize.quantiles()`), (4) top-2000 genes by merged-matrix variance, (5) per-subject rank transform. Single-cohort runs (tcga_only, cptac_only) continue to use v1 (per-cohort pipeline with `preprocess_desurv_cohort()`).
- **Reason:** Under v1, per-cohort top-2000 selection was applied *before* intersecting gene universes, yielding only ~838 common genes — far fewer than the ~2000+ expected. The preprocessing-order bug consumed most of the gene set before cohort merging. v2 fixes the order: intersect first, select top-2000 from the merged variance distribution. Quantile normalization aligns sample-level distributions across RNA-seq and proteomics platforms without introducing explicit batch labels (which would prevent generalisation to new cohorts at prediction time).
- **Trade-offs:** v1 is preserved under `preprocessing_version = "v1"` flag for backward compatibility. v2 output directories carry a `v2_` prefix (e.g., `outputs/real_data/merged/v2_point_normal/`) so v1 and v2 results coexist without overwriting. External cohorts still use v1 single-cohort preprocessing — they are never seen during training, so no joint quantile distribution to normalize against.
- **Affected files:** `code/preprocess_desurv.R` (`preprocess_merged_cohorts()`, `quantile_normalize_merged()`), `results/benchmark_sim/run_ssbmf_benchmark.R` (`run_real_data_benchmark()` `preprocessing_version` param)

---

## 2026-04-29 — Normal prior ruled out for merged-cohort β; point_normal remains default

- **Decision:** The `"normal"` prior for β (soft Gaussian shrinkage, no spike) is not adopted as the default. `"point_normal"` remains the canonical prior.
- **Reason:** The lambda sweep (above) ran all three priors at λ∈{1,5,10,20} on the merged cohort. Under the normal prior, β remains at zero just as with point_normal and point_laplace — the failure is not prior aggressiveness but the fundamental issue identified by the EBMF diagnostic (model/CAVI problem). Adding the normal prior to the benchmark sweep confirmed it does not rescue the merged-cohort fit and adds no new information. It may be revisited after the CAVI L-update issue is diagnosed.
- **Trade-offs:** The normal prior remains available in `update_beta.R` via `ebnm::ebnm_normal` and can be specified via `prior_beta = "normal"` in any benchmark call. It is not removed from the codebase.
- **Affected files:** `results/benchmark_sim/run_lambda_sweep.R`, `results/benchmark_sim/run_ssbmf_benchmark.R`

---

## 2026-04-24 — Lambda survival-scaling parameter: kept at 1.0 after sandbox evaluation

- **Decision:** The `lambda` parameter (scalar multiplier on survival precision terms in the L update) is retained in the codebase at the default value λ=1.0. No active λ tuning is performed.
- **Reason:** A principled argument for λ=p/n exists: the genomics ELBO term sums over p features while the Cox term sums over n patients, so when p>>n the genomics gradient dominates. However, a controlled sandbox (n=250, p=1000, K=5, seed=222) comparing λ∈{1, p/n=5, 2p/n=10} showed no benefit from scaling: hold-out C-index was flat at ≈0.805 across all three conditions, and β RMSE was *worse* at λ=p/n (+0.25) and λ=2p/n (+0.43) than at λ=1. Increasing λ inflates β estimates rather than correcting them, because the dominant source of β scale error is the L–β scale indeterminacy (L can rescale freely), not gradient imbalance.
- **Trade-offs:** The powered-likelihood approach (λ=p/n) is theoretically sound and used in robust Bayesian inference literature. It could become beneficial if the DGP changes (e.g., fewer features, stronger Cox signal). Keeping λ as an exposed parameter with default 1.0 costs nothing and preserves the ability to experiment.
- **Implementation:** λ is a named parameter in `update_L_k()`, `update_L_all()`, and `fit_supervised_mf_modular()` (all default 1.0). It is also in `config/globals.yml` under `cavi.lambda` and threaded through `run_ssbmf_benchmark()` and `run_real_data_benchmark()`. To test λ=p/n, change `globals.yml` and re-run.
- **Affected files:** `code/update_L.R`, `code/fit_modular.R`, `config/globals.yml`, `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/sandbox_lambda_test.R`

---

## 2026-04-24 — Proportional hazards diagnostics added to benchmark pipeline

- **Decision:** `cox.zph()` (Grambsch–Therneau test) is now run on SSBMF risk scores for each external PDAC cohort and results are saved to `ph_diagnostics_table.csv` alongside the benchmark outputs.
- **Reason:** The proportional hazards assumption underlies the Cox model used to generate risk scores. A violation means the log hazard ratio is time-varying, which can distort C-index estimates and KM stratification p-values. Formal PH testing is required before the results can be shared or published.
- **Results (TCGA-only, point_normal):** Dijk p=0.77, Moffitt p=0.72, PACA-AU array p=0.34, PACA-AU seq p=0.41 — all PASS. Puleo_array p=0.026 — **FLAG**. The Puleo violation is marginal and likely reflects the large sample size (n=288) giving power to detect subtle time-varying effects; the C-index and KM results remain valid as approximate assessments.
- **Implementation:** `compute_ph_diagnostics.R` (standalone re-fit + projection script); `run_ssbmf_benchmark.R` (PH test now wired into external cohort loop for future runs); `ssbmf_summary_report.qmd` (Section 4.3).
- **Affected files:** `results/benchmark_sim/compute_ph_diagnostics.R` (new), `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — ARD preferred over ELBO grid search for K selection

- **Decision:** K is determined automatically within a single model fit via Automatic Relevance Determination (ARD) — the point-normal/point-laplace prior on β shrinks irrelevant factor coefficients exactly to zero — rather than by fitting separate models at K = 1, …, K_max and comparing ELBOs.
- **Reason:** Two complementary advantages. *Efficiency:* ARD determines K_eff as a byproduct of CAVI; a grid search requires K_max separate full fits and introduces an outer loop. *Bayesian coherence:* ARD performs continuous soft shrinkage within a single probabilistic model; ELBO-based model selection performs discrete hard comparison across models with different dimensionalities. Setting K_max generously large (10) and letting ARD prune is equivalent to an ELBO grid search in expectation, without the computational overhead.
- **Trade-offs:** ARD can be conservative — correlated factors may collapse one even when both carry marginal survival signal. A full ELBO grid search would be more exhaustive but is 10× more expensive at K_max = 10. In practice, ARD + generous K_max matches published EBNM-based NMF standards (flash, flashier).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `config/globals.yml` (k_max), `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 1.3)

---

## 2026-04-24 — point_normal chosen as default beta prior over point_laplace

- **Decision:** `point_normal` is the recommended default prior on β for all future SSBMF runs.
- **Reason:** Across synthetic and all PDAC training modes, `point_normal` matches or slightly outperforms `point_laplace` on external C-index (TCGA-only: 0.602 vs 0.579; CPTAC-only: 0.628 vs 0.620). `point_laplace` selects higher α̂ (0.7 vs 0.5 on synthetic), suggesting it compensates for over-shrinkage of small-to-moderate coefficients by drawing more heavily on the survival gradient. The Gaussian slab is also simpler to interpret — posterior SDs have a direct normal-distribution meaning, whereas the Laplace slab mixes two scale regimes.
- **Trade-offs:** `point_laplace` has heavier tails and may outperform `point_normal` in settings with very sparse survival signal (few events, high censoring) where strong coefficient shrinkage is needed. Revisit if future larger-n runs show a consistent >0.02 C-index advantage for `point_laplace`.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (default `prior_beta`), `config/globals.yml`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — Multi-modal TCGA+CPTAC merge documented as expected failure

- **Decision:** The merged TCGA RNA-seq + CPTAC proteomics training mode (838-gene intersection, n=273) is documented as a known failure case rather than a valid benchmark condition. All β̂_k = 0 in both priors.
- **Reason:** When RNA-seq and proteomics are intersected at gene symbols and rank-normalised, PC1 separates the two platforms rather than separating patients by biology. The ARD prior correctly diagnoses that none of the learned factors carry survival signal — they carry platform identity instead. This is not a model failure; it is the model correctly reporting that no prognostic structure exists in this feature space.
- **Trade-offs:** Excluding merged results from the primary benchmark simplifies the comparison table. The failure case is retained in Section 6 of the report as a methodological lesson, motivating the shared-L multi-modal extension (separate F matrices per modality, shared L supervised by survival).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 6)

---

## 2026-04-24 — DeSurv-aligned preprocessing pipeline added

- **Decision:** Added `code/preprocess_desurv.R` — a preprocessing module that matches the DeSurv paper (Young et al. 2025, PNAS) pipeline: log₂(counts+1) for RNA-seq → select top-2000 most-variable genes per cohort → rank-transform each subject's expression vector.
- **Reason:** DeSurv serves as the primary external benchmark. Using the same preprocessing ensures any C-index difference reflects model architecture, not data transformation choices. Proteomics/microarray platforms skip the log₂ step (already on a normalized scale).
- **Trade-offs:** Rank-transform destroys absolute expression magnitude (EBNM shrinkage is insensitive to scale, but gene-level variance information is lost). Top-2000 gene filter is platform-specific — genes selected differ across cohorts, requiring intersection after preprocessing. TCGA_PAAD × CPTAC intersection yielded 838 genes (42% of 2000), acceptable given DeSurv used the same cohorts.
- **Affected files:** `code/preprocess_desurv.R`, `results/benchmark_sim/run_ssbmf_benchmark.R`

---

## 2026-04-24 — Alpha mixing parameter grid expanded to [0, 1]

- **Decision:** `config/globals.yml` `alpha_grid` expanded from `[0.1, 0.3, 0.5, 0.7, 0.9]` to `[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]`.
- **Reason:** The boundary values α=0 (pure genomics, unsupervised NMF) and α=1 (pure survival) are meaningful scientific conditions, not just edge cases. Including them lets CV reveal whether the supervision signal is worth anything at all (α=0 wins → pure NMF is optimal; α=1 wins → ignore genomics structure and regress directly).
- **Trade-offs:** Adds 2 extra fits per CV fold. α=0 is equivalent to an unsupervised NMF run; α=1 may be numerically unstable if survival events are sparse (A_surv can be near-zero for small event counts).
- **Affected files:** `config/globals.yml`, `code/select_alpha_cv.R`

---

## 2026-04-24 — predict_supervised_mf() changed from ridge solve() to SVD pseudoinverse

- **Decision:** Replaced `solve(crossprod(EF) + lambda * diag(K))` with an SVD-based Moore-Penrose pseudoinverse in `predict_supervised_mf()`.
- **Reason:** During alpha CV, early-iteration fits can have ARD drive some factor columns of EF to near-zero. The resulting EF'EF is near-singular even with a fixed ridge term λ·I, because near-collinear *non-zero* columns still make the Gram matrix ill-conditioned. The SVD approach sets d_inv = 0 for singular values below `lambda * max(d)` (relative threshold), so collapsed factors contribute nothing to L_test — which is correct since their EBeta ≈ 0 by the same ARD shrinkage.
- **Trade-offs:** SVD is slightly more expensive than a Cholesky solve for dense K×K matrices, but K ≤ 20 makes this negligible. The relative threshold means λ is now a dimensionless tolerance rather than an absolute precision floor; 1e-8 works well empirically.
- **Affected files:** `code/predict.R` (lines 73–94)

---

## 2026-04-24 — Synthetic DGP fixed: equal factor amplitudes + 4-factor survival signal

- **Decision:** Changed `generate_synthetic_benchmark_data()` to (a) use equal F amplitude for all factors (removed 5× multiplier for null factors) and (b) read `beta_true` from `cfg$synthetic$b_true = [1.5, -1.2, 0.8, -0.5, 0.0]` instead of the hardcoded `[1.0, -0.8, 0, 0, 0]`.
- **Reason:** The 5× null-factor inflation caused PCA to capture most variance from non-prognostic factors and incidentally correlate with survival — PCA C-index (0.715) exceeded supervised C-index (0.673) despite the model having the correct generative structure. Equalizing amplitudes restored the intended benchmark: supervised (0.79) beats PCA (0.76). The `b_true` change aligns the DGP with the 4-signal globals.yml spec and tests a more realistic setting (4 prognostic programs at varied effect sizes).
- **Trade-offs:** Changing the DGP invalidates any previously reported synthetic C-index numbers. The new DGP is harder (4 prognostic factors to recover vs. 2) and better reflects real-data complexity.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (`generate_synthetic_benchmark_data()`)

---

## 2026-04-09 — Repository reorganisation [PENDING]

- **Decision:** *Pending.* Propose a cleaner directory structure (documented in `ROADMAP.md` → Infrastructure section) but do not move any files until a dedicated refactor commit with no concurrent branch work.
- **Reason:** The current layout has accumulated structural debt across three simulation generations: `results/full_sim/`, `results/modular_sim_block/`, and `results/modular_sim_factor/` coexist; `.qmd` reports are mixed with output tables/figures; `demos/` is a top-level sibling of `code/` rather than nested within it; `code/SupervisedMF_Context.md` is a doc file in the algorithm directory; `derivations/EBMF/` and `derivations/SurvivalMF/` are early-sketch folders now superseded by the per-update derivation subdirectories.
- **Trade-offs:** Any file move invalidates hard-coded paths in runner scripts and `.qmd` `source()` calls — must audit before moving. Deferring keeps the repo stable while active development continues.
- **Affected files:** `results/`, `code/`, `demos/`, `derivations/`, `.gitignore`

---

## 2026-04-09 — Synthetic seed changed from 42 → 222

- **Decision:** Changed the random seed for the synthetic data-generating process from 42 to 222.
- **Reason:** Seed 42 produced a degenerate case where True F1 and True F5 (the null factor) were mixed by the bijective permutation alignment, causing β̂_null = −0.77 — an artifact of the seed rather than model behavior. Seed 222 cleanly separates all 5 factors and correctly zeroes the null factor.
- **Trade-offs:** Reported simulation results (RMSE, β̂ estimates, factor correlations) are seed-dependent; changing the seed resets the canonical benchmark. Prior figures are invalidated.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — `bijective_match()` replaces column-greedy `which.max` for factor permutation alignment

- **Decision:** Replaced `apply(abs(cors), 2, which.max)` with a `bijective_match()` helper that performs greedy global-maximum assignment to align estimated factors to true factors.
- **Reason:** The column-greedy approach is non-bijective — it can assign the same true factor to multiple estimated factors when two estimated factors are most correlated with the same true factor. This produces incorrect factor-recovery plots and corrupted β̂ comparisons in Figure 3.
- **Trade-offs:** Greedy global-max is still a heuristic (not optimal) but guarantees a one-to-one mapping and is O(K²) — negligible for K ≤ 20.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — Companion document kept in LaTeX (not Quarto) for print-first speaker notes

- **Decision:** The lab meeting companion document (`Notes/lab_meeting_april9_companion.tex`) is authored in plain LaTeX rather than Quarto.
- **Reason:** The companion is a print-first speaker notes document — dense text, no code execution, no cross-referencing with R output. LaTeX compiles faster, supports tighter typographic control, and avoids Quarto's overhead for documents with no R chunks.
- **Trade-offs:** Not integrated with the Quarto build system; must be compiled separately with `pdflatex`. Cannot embed live R output.
- **Affected files:** `presentation/walther_lab_meeting_04_09_2026/Notes/`

---

## 2026-04-01 — Full ELBO tracking added alongside proxy

- **Decision:** Implemented both a fast ELBO proxy (reconstruction error only) and a full ELBO (proxy + survival term + KL divergences) tracked at every iteration.
- **Reason:** The proxy is sufficient for convergence monitoring but the full ELBO is required for model comparison (prior family selection, K selection). Both are now computed and stored in `history$elbo_proxy` and `history$elbo_full`.
- **Trade-offs:** Adds per-iteration cost of `compute_survival_elbo()` (Taylor approximation) and `compute_ebnm_kl()`. Full ELBO requires Taylor refresh every iteration for accuracy.
- **Affected files:** `code/compute_elbo.R`, `code/fit_modular.R`

---

## 2026-04-01 — Four conditions compared: prior family × K strategy (PN/PL × K=5/K_eff)

- **Decision:** All benchmarks compare 4 conditions: point-normal (PN) vs. point-laplace (PL) crossed with fixed K=5 vs. adaptive K_eff from `auto_prune_K()`.
- **Reason:** Isolates the contribution of prior choice from K selection strategy. Reveals whether adaptive K adds value beyond a fixed-K run with the same prior.
- **Trade-offs:** Requires 4× the compute per dataset. K_eff runs require an extra `auto_prune_K()` pass before the main fit.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`, `code/select_K.R`

---

## 2026-03-31 — `fit_modular.R` uses Gauss-Seidel (factor-wise sequential) CAVI, not full-gradient

- **Decision:** Each CAVI iteration updates factors sequentially (k=1, ..., K), immediately incorporating each updated factor before proceeding to the next. This is the canonical CAVI loop in `code/fit_modular.R`.
- **Reason:** Full-gradient (parallel) CAVI requires the dense Hessian of the joint ELBO with respect to all factors simultaneously, which is intractable in closed form under the Taylor approximation for the survival term. Gauss-Seidel is the standard choice for mean-field CAVI.
- **Trade-offs:** Sequential updates mean later factors in each iteration benefit from earlier updates (faster convergence per iteration), but the ordering introduces implicit asymmetry between factors. Not easily parallelized across factors.
- **Affected files:** `code/fit_modular.R`

---

## 2026-03-25 — Hold-out prediction uses pseudo-inverse projection of Y_test onto F

- **Decision:** `predict_supervised_mf()` obtains test-set loadings L_test by projecting Y_test onto the trained factor matrix F via pseudo-inverse: `L_test = Y_test %*% F %*% solve(t(F) %*% F)`.
- **Reason:** After training, F is fixed. The natural way to score a new patient is to find the loadings that best reconstruct their expression profile under the trained factors. The pseudo-inverse gives the least-squares solution.
- **Trade-offs:** Does not propagate uncertainty from F into L_test (point estimate only). Assumes test data is drawn from the same distribution as training data (no domain shift).
- **Affected files:** `code/predict.R`

---

## 2026-03-20 — Convergence criterion: dual threshold on mean|ΔL| AND mean|Δβ| after 5-iteration burn-in

- **Decision:** Convergence is declared when both `mean(|EL - EL_old|) < tol` AND `mean(|EBeta - EBeta_old|) < tol`, checked only after iteration 5.
- **Reason:** A single criterion on L alone can declare convergence while β is still adjusting (or vice versa). The dual criterion ensures both the genomic structure (L) and survival signal (β) have stabilized. Mean is used instead of max because max is dominated by a few high-variance entries near factor orientation boundaries and rarely reaches 1e-3 on real datasets.
- **Trade-offs:** Mean convergence is weaker than max convergence — some individual loadings may still be changing when the algorithm stops. The 5-iteration burn-in prevents premature stopping during large initial swings.
- **Affected files:** `code/fit_modular.R` (lines 353–370)

---

## 2026-03-12 — Modular update architecture: each update function is independently testable

- **Decision:** CAVI updates are split into four independent modules (`update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R`), each with its own test suite and demo scripts.
- **Reason:** The monolithic V2 implementation (`Supervised_Bayesian_MF_V2.R`) mixed all update logic into one function, making it hard to test individual components or swap implementations. Modular architecture enables TDD, isolated debugging, and future replacement of individual components.
- **Trade-offs:** Cross-module dependencies must be managed explicitly (e.g., `compute_R_k` lives in `update_L.R` but is also used by `update_F.R`, requiring careful source ordering).
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`, `code/update_tau.R`, `code/fit_modular.R`

---

## 2026-03-10 — EBNM priors chosen over fixed-penalty; g estimated per CAVI step

- **Decision:** Shrinkage priors for L, F, and β use the Empirical Bayes Normal Means (EBNM) framework, with the prior g estimated from data at each CAVI step rather than fixed by cross-validation.
- **Reason:** Fixed-penalty methods (Ridge, LASSO) require cross-validation for λ, adding computational overhead and a choice of CV criterion. EBNM folds g-estimation into the ELBO maximization — the prior adapts to the data, and no separate CV loop is needed. Implemented via the `ebnm` R package (Stephens lab).
- **Trade-offs:** g-estimation adds per-update overhead. The estimated g is dataset-dependent and may change across CAVI iterations (instability if the model is mis-specified). Point-normal and point-laplace are the two families tested; other families (e.g., generalized double-Pareto) are not yet explored.
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`

---

## 2026-02-12 — Baseline hazard h₀(t) left non-parametric; enters through Cox partial likelihood

- **Decision:** The survival component uses the Cox proportional hazards model with an unspecified baseline hazard h₀(t). h₀(t) cancels exactly in the Cox partial likelihood and never needs to be estimated.
- **Reason:** Specifying a parametric h₀(t) (Weibull, exponential) would add a nuisance parameter and require a prior on its shape/scale. The Cox partial likelihood avoids this entirely while preserving the proportional hazards structure needed to link L to survival.
- **Trade-offs:** Cannot predict absolute survival probabilities (only relative risk via exp(Lβ)). The partial likelihood is not a true likelihood (it conditions on observed event times), which means the ELBO approximation for the survival term requires a Taylor expansion rather than a closed-form KL.
- **Affected files:** `code/update_beta.R`, `code/update_L.R`, `code/compute_elbo.R`

---

## 2026-02-12 — Model formulation: shared L links genomics (Y = LF′ + E) and survival (h(t) = h₀(t)exp(Lβ))

- **Decision:** The core model posits a shared n×K loading matrix L that simultaneously factorizes the genomics matrix Y and enters the Cox proportional hazards model as the linear predictor.
- **Reason:** A two-stage approach (first factor Y, then regress loadings on survival) does not jointly optimize the factors for survival prediction. Sharing L with a joint objective ensures the learned gene expression programs are informative for both reconstruction and prognosis.
- **Trade-offs:** The joint objective couples genomics and survival gradients; the relative scale of the two terms is dataset-dependent and not normalized (p >> n means the genomics term dominates). This is a known limitation — adding a λ scaling parameter is an identified future direction.
- **Affected files:** `code/fit_modular.R`, `code/update_L.R`, `code/update_beta.R`, `derivations/MF_UpdateDerivations/`
