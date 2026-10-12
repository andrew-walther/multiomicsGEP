# Review: `docs/progress_book/chapters/2026-10-02.qmd`

Reviewer: independent check of accuracy and presentation, 2026-10-02.
Version reviewed: the 390-line file, which includes the program table with HRs at lines 231–258. The file changed while this review was running; line numbers refer to that version.

Sources checked:

- **Real-data results:** `real_sweep_summary.csv`, `real_factor_comparison.csv`, `real_fit_cindex_*.csv`, and the four pruned/original fit RDS files in `results/multimodal_real/outputs/fits/`.
- **Simulation sweep:** `pruning_comparison_quick.csv` and `pruning_comparison_summary_quick.csv`.
- **Single-modality external C:** `external_cindex_ci.csv`, `external_paired_diff_ci.csv` and `km_risk_tertile_stats.csv`.
- **Scripts and code:** `make_km_figures.R`, `run_multimodal_real_fit.R`, `run_multimodal_pruning_comparison.R`, `code/fit_multimodal_yfb.R` and `code/multimodal_yfb_updates.R`.
- **Decisions and derivations:** the DECISIONS.md entries for 2026-10-01 and 2026-10-02, and the derivation `.qmd`/`.tex` files.
- **Rendered output:** `_book/chapters/2026-10-02.html`.

Full simulation sweep: `pruning_comparison.csv` is a 72-byte placeholder and no `pruning_comparison_summary.csv` exists yet. The chapter's quick-sweep numbers therefore remain the ones to cite, and nothing needs updating from the full sweep yet. The real-data runs at K_init = 4, 5 and 15 are still in progress, which matches what the chapter says.

---

## Summary verdict

Most numbers match their sources. These were all reproduced: the §4 C-indices and CIs, the K_init table, the variance shares, the KL column, the top genes, the methylation-scale table, the §5 external C, ΔC and KM p-values, and the α_F values.

The chapter is not ready to present as written, for four reasons.

1. **Its headline comparison leaves out the most relevant benchmark.** The project's own recommended single-modality model already reaches C = 0.657 on PACA-AU RNA-seq (n = 52). That is almost certainly the same ICGC donors as the "ICGC primary PDAC" set (n = 50). Against that benchmark, "0.45–0.49 → 0.68" overstates the advance.
2. **Some quick-sweep and real-data numbers or claims are wrong:**
   - the C ranges;
   - the convergence rate of the original fit;
   - which factor is the largest expression factor;
   - the range of loading correlations between the joint and expression-only fits.
3. **Rendering problems:**
   - Five bullet lists have no preceding blank line, so they render as run-on paragraphs with literal " - " (confirmed in the rendered HTML).
   - A planning block and a `[TABLE + FIGURES: …]` placeholder were left in §4.
4. **Several claims are stronger than the evidence:**
   - "separation in all five cohorts" (only 3 of 5 have log-rank p < 0.05);
   - "Every update equation is in the derivation" (the intercept, point-Laplace and ELBO pruning are not in it);
   - "a factor that survives pruning explains variance or survival" (the K_init = 10 fit kept an inactive factor).

---

## Must-fix

**1. The headline gain omits the project's own single-modality result on the same cohort.**
- *Location:* the callout (line 6), "On real data it raises external C on ICGC from about 0.45–0.49 to 0.68". Also line 15 and §4 line 220.
- *Problem:* the comparators at 0.45–0.49 are weak configurations:
  - The single-modality YFB arm (0.491) uses raw log2 inputs, not the recommended per-platform z-standardization (see the header of `run_multimodal_real_fit.R`).
  - It stopped at a 100-iteration cap without converging (`converged = FALSE`, `iterations = 100` in `real_fit_cindex_*_K7.csv`).
  - Meanwhile §5 reports the recommended single-modality model (trained on TCGA + CPTAC) at **C = 0.657 (0.544–0.772) on PACA-AU RNA-seq, n = 52** (`external_cindex_ci.csv`). That is ICGC PACA-AU RNA-seq, the same source as the n = 50 set here. DECISIONS 2026-10-01 says ICGC survival "agrees exactly with `PACA_AU_seq.survival_data.rds`".
  - An advisor will put the two numbers side by side.
- *Fix:*
  - Add the 0.657 to §4 and the executive summary, noting the differences: training cohort (TCGA vs TCGA + CPTAC), preprocessing (raw log2 vs z-std), and n = 50 vs 52.
  - Reword along these lines: "The new framework (0.685) is above the weak in-run baselines (0.42–0.49) and comparable to the recommended single-modality model on PACA-AU RNA-seq (0.657); the intervals overlap."
  - Check the donor overlap between the two sets before the meeting.

**2. The expression-only two-step baseline is missing.**
- *Location:* §4 table, line 214; executive summary line 15 ("the two-step model 0.424").
- *Problem:* the expression-only run also fit a two-step EBMF → Cox on expression alone. It reached **0.525 (0.399–0.651)** on ICGC primary, 0.527 (0.425–0.641) on ICGC all, and 0.670 on TCGA (`real_sweep_summary.csv`, config `pruned_expression_only`). Only the stacked two-step arm (0.424) is reported. Stacking methylation therefore *lowered* the two-step baseline by about 0.10. Leaving this out makes the baselines look uniformly below 0.5.
- *Fix:* add a row "Two-step EBMF → Cox | expression only | 0.525 (0.399–0.651) | 0.527 | 0.670". Revise "Both baselines fall below 0.5 on ICGC" (line 267) to match.

**3. The quick-sweep C ranges are wrong.**
- *Location:* §2 table, lines 150–151.
- *Problem (source: `pruning_comparison_summary_quick.csv`, pruned rows):*
  - "Two prognostic programs … 0.787–0.788": the actual range is **0.785–0.788** (K_init = 15 gives 0.7845).
  - "Adverse + protective … 0.783–0.795": the actual range is **0.766–0.795** (K_init = 3 gives 0.766).
- *Fix:* correct both ranges.

**4. The original fit's convergence rate is misstated.**
- *Location:* line 156, "the original fit converged in 0–50% of runs".
- *Problem:* for adverse + protective at K_init = 10, the original fit converged in 2 of 2 seeds (100%). Overall the original fit converged in **8 of 30 runs**: 0–100% per setting, 0% in the two-prognostic scenario. The statement that every pruned fit converged (30 of 30) is correct.
- *Fix:* "the original fit converged in 8 of 30 runs (none in the two-prognostic scenario)".

**5. Factor 1 is not the largest expression factor.**
- *Location:* line 262, "Its largest expression factor corresponds to Program 3 (|r| = 0.79)".
- *Problem:* the factor that matches P3 is factor 1, which explains 12.8% of the total variance. The chapter's own table (lines 237–240) shows factors 3 (16.7%) and 6 (16.0%) are larger. By expression-only fitted signal (‖L_k‖²‖F_expr,k‖² from the K7 fit), the shares are factor 3 31.5%, factor 6 28.9%, factor 1 23.5%. The sentence also contradicts the table directly above it.
- *Fix:* "Factor 1 (classical/gastric epithelial, 12.8%) corresponds to Program 3 (|r| = 0.79), although here it is not survival-active."

**6. The range of matched loading correlations is overstated.**
- *Location:* line 254, "Matched loadings between the two fits correlate at |r| 0.79–0.91, and the survival-active program matches at 0.88."
- *Problem:* I took the best |cor| of each joint expression factor against the expression-only K7 fit's loadings (3,000 shared genes; both RDS files):

  | Joint factor | Program | Best |r| |
  |---|---|---|
  | 1 | classical | 0.87 |
  | 2 | immune | **0.67** |
  | 3 | mixed | 0.79 |
  | 4 | basal-like | 0.88 (correct) |
  | 5 | stroma | **0.52** |
  | 6 | exocrine | 0.91 |

  The quoted range covers only four of six factors. The two it leaves out are the immune factor (the one shared with methylation) and the stroma factor. Leaving them out weakens the claim that "expression programs are essentially the same with or without methylation".
- *Fix:* "|r| 0.52–0.91 (four of six above 0.79; the immune and stromal programs match less well, 0.67 and 0.52)". Then soften "essentially the same".

**7. A planning block and placeholder were left in §4.**
- *Location:* lines 289–298, the "**Arms.** … **Starting ranks.** … **Metric.** … **Factors.** …" bullets and then `[TABLE + FIGURES: real_sweep_summary.csv, real_sweep_factor_roles.png, real_sweep_cindex.png, real_factor_comparison.csv]`.
- *Problem:* this is pre-results design text. It sits under no heading and has a literal placeholder. It also describes arms ("two-step EBMF → Cox on the stacked modalities") already reported above.
- *Fix:* delete it, or move the Metric bullet under the first §4 table as a methods note. Either insert `real_sweep_cindex.png` or drop the placeholder.

**8. Lists render as run-on paragraphs.**
- *Location:* lines 260–263 ("**The programs match the single-modality ones.**" then "- …"), 265–269 ("**Caveats.**"), 313–315, 319–321 and 364–367 ("**Test design, fixed in advance.**").
- *Problem:* Pandoc needs a blank line between a paragraph and a list. The rendered HTML shows "Caveats. - Small validation set. … - Weak baselines. …" as a single paragraph (`_book/chapters/2026-10-02.html` line 881; also line 1122).
- *Fix:* insert a blank line after each bold lead-in.

**9. The KM claim of separation in all five cohorts is overstated.**
- *Location:* executive summary line 19 ("Held-out Kaplan–Meier curves show survival separation in all five external cohorts"), and the line 313 heading ("The frozen model separates survival in every held-out cohort").
- *Problem:* the tertile log-rank p-values (`km_risk_tertile_stats.csv`) are Dijk 0.00046, Moffitt 0.17, PACA-AU array 0.075, PACA-AU seq 0.021 and Puleo < 0.0001. Only three of five are below 0.05, and §5 itself calls Moffitt the weakest.
- *Fix:* "Risk tertiles are ordered as expected in all five cohorts; log-rank p < 0.05 in three (Puleo, Dijk, PACA-AU seq), p = 0.075 in PACA-AU array and p = 0.17 in Moffitt." First check by eye that the curves are in fact ordered in all five.

**10. The cited derivation does not contain the new framework.**
- *Location:* line 49, "Every update equation is in `derivations/multimodal_YFB/multimodal_YFB_derivation.pdf`."
- *Problem:* the derivation `.qmd` has no intercept μ_m, no point-Laplace prior, no ELBO and no pruning or nullcheck. A search for "intercept", "laplace", "ELBO" and "prun" finds 0 matches. Those pieces exist only in code comments and DECISIONS.md (2026-10-01).
- *Fix:* "The base updates are in …; the intercept, signed priors, ELBO and pruning are documented in DECISIONS.md (2026-10-01) and `code/fit_multimodal_yfb.R`, and will be added to the derivation."

---

## Should-fix

**11. The classification statement is contradicted by the K_init = 10 fit.**
- Line 117 reads: "A factor that survives pruning either explains molecular variance, explains survival … or both."
- The K_init = 10 real fit kept factor 9, which is classified `inactive`: PVE 0.44% (below the 1% cutoff) and prognostic z = 0.
- This also undercuts line 285, "Pruning should start removing factors only at larger K_init, which the sweep will show". At K_init = 10, nullcheck kept even an inactive factor.
- Say this, and present the expectation for larger K_init as a hypothesis.

**12. The pruning rule is described as running after convergence.**
- Line 105 reads "After convergence, remove the factor …".
- In `fit_multimodal_yfb()`, the nullcheck runs when the outer loop ends, whether it converged or hit `max_outer`. The real K7 and K10 fits had not converged when it ran.
- The nullcheck also holds the other factors fixed: their moments are not re-optimized, only μ and τ are re-estimated. The code comment calls this a conservative test.
- State both points. They matter for interpreting "nothing was pruned".

**13. The convergence diagnostic is misdescribed.**
- Line 269 reads: "the largest per-feature reconstruction change stayed around 1%. That maximum over 13,000 features is being held up by a few slowly moving features."
- `multimodal_yfb_relative_change()` returns max|new − old| over all n × p *entries* of LF_mᵀ, divided by the largest absolute entry. The maximum is then taken over the two modalities. It is not a per-feature quantity.
- The final value is 0.0145, and the η change is 3.56 × 10⁻⁴; both are correct as quoted.
- The claim that "a few slowly moving features" drive it is not supported by any saved diagnostic. Either show it, for example with the distribution of per-feature changes, or drop it.
- Two related points:
  - The expression-only framework fit also hit the 500-sweep cap (reconstruction change 2.5 × 10⁻⁴). This is not mentioned.
  - The K_init = 3 fit converged at 398 sweeps.

**14. The quick-sweep summary hides seed-to-seed differences in the low-variance scenario.**
- In the low-variance scenario, held-out C "0.61" is the mean of 0.727 (seed 20261002) and **0.500** (seed 20261003).
- In seed 20261003 the pruned fit has no survival-active factor at all, and the low-variance loading is recovered at |r| = 0.06.
- In seed 20261002 there are two false survival-active factors and the loading is recovered at |r| ≈ 0.44–0.47.
- So "one false survival-active factor" (line 158) is an average of 2 and 0. Report the two seeds separately, or say "0.50–0.73 across seeds".

**15. False survival-active factors also occur in the scenarios where pruning works.**
- In the pruned fits, the mean number of false survival-active factors is 0.5–1 at every K_init in the two-prognostic and adverse + protective scenarios.
  - For example, two-prognostic seed 20261003 has K_survival = 3 at every K_init, so the non-prognostic true factor is flagged survival-active even at K_final = 3.
  - The recovery figure shows this, but the text attributes false survival-active factors only to "duplicates at large K_init" (line 159).
- Executive summary item 2 ("recovered with the correct signs") should read "found survival-active, with the correct sign, in 75–100% of runs (both programs in 50–100% at each K_init)". Adverse + protective at K_init = 7 and 10, seed 20261002, flagged only one of the two programs.

**16. "About the true K = 3 for every starting value" needs a qualifier.**
- Lines 11 and 146 claim this. In the low-variance scenario the pruned K_final is 2–2.5, and adverse + protective at K_init = 3 gives 2.5.
- Add "in the two scenarios where every program carries variance; 2–2.5 when the prognostic program is low-variance".

**17. Minor numeric mismatches.**
- Line 187: median |skewness| on the asin scale is **0.33** (0.3348), not 0.34 (`methylation_distribution_summary.csv`).
- Line 152: "2 of 3" true factors recovered is 1.5 at K_init = 3.
- Line 277: the KL for factor 3 computes to 9,351 from `fit$kl` (rounding; trivial).
- Line 390 says "(21 pp)", but line 27 says "22 pp" and `paper/prelim/project3-ssbmf.pdf` has 22 pages. Use 22.

**18. The source of the orientation choice is ambiguous.**
- Line 329 reads "A single global flip, risk = −η, chosen before seeing any validation outcome".
- The rationale in `make_km_figures.R` and DECISIONS 2026-10-02 cites "Raw concordance with reverse = TRUE is below 0.5 in all five cohorts", which is a validation-set observation. The reported table also already existed.
- The defensible statement is that the flip follows from training-set quantities: the sign of β on the training-defined adverse and protective programs (β₇ = −0.040, β₃ = +0.012). Better still, report the training C of −η.
- Reword to "determined from the training fit's coefficient signs".
- The executive-summary phrase "training-fixed" (line 20) needs the same support.

**19. The program-table values and CpG counts cannot be reproduced from saved outputs.**
- Lines 231–241 and the executive summary (line 17) give HRs (for example ICGC HR 1.99, p = 0.0003) and CpG counts (1,589 and 4,477). No CSV or script in `results/multimodal_real/` produces them.
- With point-Laplace posterior means, almost all CpGs are nonzero: 9,944 of 10,000 for factor 2 at |E f| > 10⁻⁸. The counts therefore depend on an unstated threshold.
- "Uses essentially no CpGs" (line 256) and "Only factor 2 loads on both modalities" (line 287) depend on the same threshold. By fitted methylation signal, factor 6 carries 9.3% of the methylation signal, against 20.9% for factor 2, and factor 4 carries 4.5%.
- Add a script and output (for example `real_program_table.csv`), state the threshold or criterion, and say whether the HR is per SD of the training or the validation projection.
- Note that 7 programs × 2 cohorts makes 14 marginal tests, so the factor 3 ICGC result (p = 0.05) should not be read as a finding.

**20. The scale-imbalance explanation is stated as causal.**
- Line 286 reads "this is also why the simulated low-variance prognostic program is lost".
- It is plausible, but no experiment isolates it, for example by up-weighting survival in the simulation.
- Present it as the leading hypothesis.
- Also, only the 75 events contribute terms to the partial likelihood, not 144.

**21. The single-modality comparator caveat is incomplete.**
- Line 268 lists only "no intercept and no point-Laplace loadings".
- Add: raw log2 rather than the recommended z-standardized inputs, and not converged (100-iteration cap).

**22. The per-run C-index CSVs have stale CIs.**
- In `real_fit_cindex_K7.csv` and the `*_pruned_K3/K7/K10.csv` files, some CIs exclude their own point estimate. For example:
  - two-step, ICGC primary: 0.424 with CI 0.433–0.714;
  - original joint: 0.453 with CI 0.401–0.681.
- These look like per-replicate re-orientation (flip = NULL). `real_sweep_summary.csv`, which the chapter uses, has the correct CIs.
- Regenerate or delete the per-run CSVs so nobody cites them by mistake.

---

## Nice-to-have

**23. Executive summary length.**
- At about 35 lines it is closer to 4–5 minutes than 2–3.
- Suggested cuts:
  - fold "Infrastructure" into one line;
  - drop item 4's orientation sub-bullet (keep it in §5);
  - merge "Open problems" and "For discussion", where the survival-in-selection point appears twice (lines 32 and 43).

**24. Jargon to define for a reader coming in cold.**
- YFB (η = (YF)β), "two-step model" (EBMF → Cox), `ebnm`, "nullcheck" (a flashier term), and α/α_F. These appear in the executive summary on line 32 without definition.
- PVE and "survival-active" (|Eβ|/SDβ ≥ 1.96).

**25. Wording against the no-mannered-prose rule.**
- "Give survival a voice in factor selection" (line 376), "How should survival get a say" (line 43) and "survival still has almost no say" (line 6).
- Suggested literal alternative: "make survival affect which factors are retained".

**26. The β prior acts as ARD.**
- `multimodal_yfb_update_beta_k()` sets ŝ²_β = max(0, x² − 1/A). Whenever the per-factor z < 1, β_k is exactly 0.
- This is why most factors show β = 0.00 and z = 0. One sentence in §1 would explain the many zero rows.

**27. Which programs are survival-active depends on K_init.**
- At K_init = 10, the survival-active factors are:
  - k1 (β = +0.40, best match P3 at |r| = 0.68);
  - k5 (β = +0.24);
  - k10 (β = −0.24, P7 match only 0.38).
- At K_init = 7 there is a single basal-like factor. This is worth one line, because the executive summary presents the P7 match as a finding.

**28. Figure label differs from the script comment.**
- The K_final figure strip says "Moderate noise (factor PVE ~10%)", while the script comment says noise_scale = 0.25 is "~25%" signal-to-noise.
- Reconcile the two if the noise level is quoted.

**29. Trailing whitespace.**
- Line 389 ends with a space.
- The continuation sentence at line 380 attaches to the last nested bullet rather than the parent item. Indent it to the parent's level, or make it its own paragraph.

**30. The internal label is in the linked β derivation, not this chapter.**
- The β derivation this chapter links to uses "Cluster B" (`qBeta_YF_derivation.tex` line 46). That breaks the repo writing rule, but the problem is in the derivation, not here.

---

## Verified claims (correct as written)

**§4 real-data comparison (`real_sweep_summary.csv`)**
- C-indices, ICGC primary:
  - new framework, expression + methylation: 0.685 (0.579–0.779);
  - new framework, expression only: 0.676 (0.569–0.770);
  - original multimodal fit: 0.453 (0.319–0.599);
  - single-modality YFB: 0.491 (0.364–0.633);
  - two-step: 0.424 (0.286–0.567).
- C-indices, ICGC all and TCGA: 0.639 (0.547–0.734), 0.653, 0.442, 0.480 (0.369–0.581) and 0.418 (0.310–0.531); TCGA 0.638, 0.644, 0.597, 0.594 and 0.579.
- Methylation increment: +0.009 ≈ +0.01.
- K_init table: factors kept 3/7/10; survival-active 0/1/3; C 0.557 (0.420–0.699), 0.685 and 0.650 (0.529–0.763); no factors pruned. Only the K_init = 3 fit converged.
- Orientation of the multimodal risk scores is frozen (flip = FALSE in `run_multimodal_real_fit.R`).

**K7 joint fit**
- Variance shares 12.8, 4.3, 16.7, 3.6, 2.6, 16.0 and 4.6% (`factor_pve`).
- Factor 4 is the only survival-active factor: z = 2.56, β = +0.228, adverse direction.
- Top genes for all seven factors match the posterior-mean loadings.
- KL-saved column: 22,602; 16,504; 9,351; 20,323; 14,310; 18,227; 28,895.
- ELBO components: reconstruction 610,423, survival −304, KL 130,212.
- η relative change 3.6 × 10⁻⁴; stopped at the 500-sweep cap.
- 144 × 13,000 ≈ 1.9 M reconstruction terms.

**Factor comparison (`real_factor_comparison.csv`)**
- Factor 4 ↔ P7 at |r| = 0.68 over 1,781 genes; factor 1 ↔ P3 at 0.79.
- K7 factors 1–6 vs `tcga_flash_K14`: best |r| 0.60–0.87. Factor 7 is methylation-only and is not compared, so "every factor" should read "every expression-loading factor".

**Cohorts (DECISIONS 2026-10-01)**
- n and events: 144/75, 50/30, 67/40.
- 6 of 150 TCGA patients dropped.
- 10 + 4 + 2 + 1 = 17 = 67 − 50.
- 305,352 CpGs; 0.07% imputed; screening 3k genes and 10k CpGs on TCGA only.
- DO33168 has two RNA samples.

**Methylation-scale table**
- Matches except the 0.34 noted in item 17: 0.43, 13%, −0.73, 92%; 6.6%, −0.70, 82%.

**Quick sweep**
- 60 fits with no errors.
- Pruned K_final for all three scenarios.
- Prognostic recall 100% (75% at K_init = 15) and 75–100%.
- 0% in the low-variance scenario.
- All 30 pruned fits converged.
- The original fit keeps K_init factors.

**Single-seed table**
- The low-variance |r| = 0.47 and C = 0.728 match quick-sweep seed 20261002 at K_init = 5.
- The other rows are consistent with DECISIONS 2026-10-01. They could not be traced to a saved output.

**§5 single-modality results**
- Mean external C 0.627 vs 0.581.
- YFB is higher in every cohort; only Puleo's paired difference is significant.
- Pooled ΔC = 0.042 (0.013–0.071).
- Cohort Cs: 0.634, 0.549, 0.648, 0.657, 0.645.
- Dijk p = 0.00046 and Moffitt p = 0.17.
- Pooled n = 616; P7 p < 0.0001 and P3 p = 0.0041, with the directions shown in the figure.

**§6 α_F results**
- Mean external C 0.6267, 0.6268, 0.6263 and 0.6066 (DECISIONS 2026-08-20).

**Model and code statements**
- Intercept re-estimated every sweep as colMeans(Y − LFᵀ); a constant η shift leaves the partial likelihood unchanged.
- Survival is absent from the L update.
- Rescaling sets sd(Y F_k) = 1.
- The β update uses A = Σw(E[Z]² + Var Z).
- ELBO = reconstruction + Breslow log PL(E η) − KL.
- The nullcheck re-estimates μ and τ for each candidate set.
- Survival-active cutoff z ≥ 1.96; PVE cutoff 1%.

**Figures**
- All five `../figs/2026-10-02_*.png` paths resolve under `docs/progress_book/figs/`.
- Each figure is byte-identical to its source output.
- The chapter is listed in `_quarto.yml`.
- No "Phase", "Cluster" or "Session N" labels appear in the chapter.
