# Likely Committee Questions for Part III (SSBMF)

Twenty questions a committee member is likely to ask about Chapter 4 at the oral prelim
(2026-12-07). Use them to plan the slide deck: questions the main deck should answer before
they are asked, and backup slides for the rest. They come from the round-5 committee review
of the chapter, part (ii) of
`bios-dissertation/prelim/full-draft/notes/human-writer-review/round5/ch4/committee-proposals.md`
(2026-11-15).

The suggested answers describe the chapter at multiomicsGEP commit 0bb0fc8f, before the
round-5 edits. Check each answer against the current chapter before using it on a slide.

## Where Each Question Belongs in the Deck

Format: **main** (answer on a main-deck slide before it is asked), **backup** (a backup slide
ready to show) or **spoken** (answered aloud, no slide).

| Question | Deck section | Format |
|---|---|---|
| Q-01 What the joint model buys over two steps | 1. Introduction / 4. Discussion | main |
| Q-02 Why not raise $w_S$ | 3. Results, simulation | main |
| Q-03 Partial likelihood in a variational objective | 2. Methods | backup |
| Q-04 $\mathbf{Y}$ used twice | 2. Methods | spoken |
| Q-05 Measurement error and coefficient bias | 1. Introduction | backup |
| Q-06 Positive training concordance difference | 3. Results, PDAC | spoken |
| Q-07 Meaning of survival-active | 2. Methods | backup |
| Q-08 Rank chosen with non-converged fits | 3. Results, PDAC | backup |
| Q-09 Point-Laplace against nonnegative priors | 2. Methods | spoken |
| Q-10 Validation standardization and missing genes | 3. Results, PDAC | spoken |
| Q-11 One baseline hazard | 2. Methods | spoken |
| Q-12 Other comparators | 3. Results, PDAC | spoken |
| Q-13 Dropped Hessian term | 2. Methods | backup |
| Q-14 Program 4 and exocrine content | 3. Results, PDAC | spoken |
| Q-15 Methylation scale and priors | 4. Continuing work | spoken |
| Q-16 What the joint expression and methylation model should show | 4. Continuing work | main |
| Q-17 SafeBayes and Holmes-Walker | 4. Continuing work | backup |
| Q-18 Risk-group cutpoint | 3. Results, PDAC | spoken |
| Q-19 Other Bayesian survival-guided factor models | 1. Introduction | spoken |
| Q-20 What the external concordance claim is | 3. Results, PDAC | main |
| Q-21 Coverage of the coefficient intervals | 3. Results, simulation | backup |

### Q-01. On these data the joint and two-step fits give the same programs and the same predictions. What does the Bayesian joint model buy you?

**Suggested answer.** On PDAC at $w_S=1$, only the coefficients differ: they are shrunk under the
Normal prior, and Program 1's coefficient is set to zero (C-08). The programs are the EBMF
programs (line 289, Appendix C.12). The case for the joint model rests on two things the chapter
shows elsewhere. In simulation near the detection threshold, the joint fit recovered a program
the two-step fit missed (5 of 10 against 0, concordance gain 0.156, Table `tab:strength` panel
B). And the survival-share derivation explains why the PDAC data do not reach that regime: the
share scales with $E[\beta_k^2]$, the real coefficients are at most 0.42, and the programs carry
large variance shares. The proposal is that a data-driven weight will let survival act where the
data support it. The posterior intervals are a second promised benefit and are not yet
calibrated (C-01). Be direct that the current PDAC result is a tie with the two-step fit.

### Q-02. Survival supplies $10^{-6}$ of the loading precision. Why not just raise $w_S$ to $10^4$?

**Suggested answer.** The approximation $\rho_{jk}\approx w_SE[\beta_k^2](D/n)/(\tau_jp^2v_k^2)$
(Appendix C.2) says matching the balance at fixed variance share needs $w_S\propto p^2$. In the
fixed-weight simulation with 2,000 genes, every fit at $w_S\ge40$ stopped at 300 sweeps without
converging and lost the planted program in most cells (line 216). The cause is not identified.
One hypothesis in the derivation (Sec. `sec:strengthen`) is that with $p>n$ the survival-only
estimate can order the training times exactly, and a large weight lets the loadings follow it.
The coefficients in those fits (17 to 26) were near the true values (9 to 26), so the failure was
not runaway coefficients. The continuing work compares calibrated and normalized weights with the
cross-validated one.

### Q-03. Is it legitimate to put a Cox partial likelihood in a Bayesian or variational objective?

**Suggested answer.** Two justifications. Kalbfleisch (1978) and Sinha, Ibrahim and Chen (2003)
show the partial likelihood is the limiting marginal posterior for the coefficients under a
diffuse gamma-process prior on the baseline hazard. Bissiri, Holmes and Walker (2016, Sections 4
and 5.1) treat the negative partial log-likelihood as a loss in a generalized Bayesian update and
apply it to a genetic survival analysis. The variational objective is then the expected loss plus
the KL divergence from the prior, restricted to the mean-field family (their Eq. 7). The fit also
replaces the partial likelihood by a quadratic expansion at each sweep (line 103), so the
objective actually maximized is a surrogate.

### Q-04. You use $\mathbf{Y}$ twice, as the outcome of the factorization and as the covariate in the Cox model. Is that a problem?

**Suggested answer.** The joint density factors as $p(\mathbf Y\mid\mathbf L,\mathbf F,\boldsymbol\tau)\,p(\mathbf t,\boldsymbol\delta\mid\mathbf Y,\mathbf F,\boldsymbol\beta)$
(derivation Eq. `eq:joint`). The survival model conditions on the observed expression, as any
Cox model conditions on its covariates, and $\mathbf F$ appears in both factors. That shared
parameter is how survival informs the loadings. No density for $\mathbf Y$ is counted twice. The
alternative, $\boldsymbol\eta=\mathbf L\boldsymbol\beta$, conditions on the latent scores and gave
nearly the same fit on these data (Appendix C.11).

### Q-05. Your Introduction says two-step pipelines ignore measurement error in the scores. Does your coefficient update correct the bias?

**Suggested answer.** It accounts for loading uncertainty but does not de-attenuate. The update
uses $E[z_{ik}^2]=\bar z_{ik}^2+\sum_jy_{ij}^2\operatorname{Var}(f_{jk})$ in the precision $A_{\beta_k}$,
which adds a ridge-type penalty and shrinks the coefficient toward zero. Classical
measurement-error corrections such as regression calibration move a naive estimate away from
zero. Here the score used for prediction is $\mathbf Y E[\mathbf F]$ and the error is in the
loadings, closer to a Berkson structure. In simulation the joint coefficients were more
attenuated than the two-step coefficients (slopes 1.45 to 3.23 against 0.96 to 1.16), and their
intervals rarely covered the truth (C-01). The honest summary is that the update propagates
uncertainty into the coefficient's posterior and that its calibration is still open.

### Q-06. Why is the training concordance difference positive (0.004) if the programs are identical?

**Suggested answer.** Concordance depends on the relative weights of the program scores. The
joint update shrinks each coefficient by the factor $1-s_k^2/x_k^2$, which differs by program.
Against the two-step Cox coefficients on the same scores, the joint coefficients are 0 against
0.070 for Program 1, 0.101 against 0.147 for Program 3, and 0.209 and 0.231 against 0.247 and
0.272 for Programs 4 and 5, with Program 2 nearly unchanged ($-0.372$ against $-0.388$)
(`k5_working_2026-11-12/joint_vs_twostep_check.txt`). The difference is within one SE and
disappears externally ($-0.001$).

### Q-07. What does "survival-active" mean statistically? Is 1.96 a 95% credible interval?

**Suggested answer.** With the Normal prior refit from one pseudo-observation, the ratio
$E[\beta_k]/\mathrm{SD}(\beta_k)$ equals $\sqrt{x_k^2/s_k^2-1}$, so the rule is a threshold of 2.20
on the working Wald ratio (C-08). It ignores uncertainty in the estimated prior variance and the
mean-field independence of $\beta_k$ from the loadings. It is a screening label. The main-program
rule adds fold stability and contribution to the risk score (Appendix C.6), and only Program 2
passes it.

### Q-08. The rank was chosen with fits that did not converge. Do you trust $K=5$?

**Suggested answer.** At $K=5$ all five fold fits converged in the first partition and all 25 in
the four repeats. Non-convergence grows with $K$ (0 of 5 at $K=12$), and those fits are included
in the curve. The one-SE rule picked $K=5$ on a plateau from 5 to 12 (0.649 to 0.666), and four of
five partitions picked 5 on their own (Appendix C.9). $K=5$ and $K=7$ predict equally well. The
cost of the choice is that Programs 3 and 4 sit near the activity threshold and change label
between ranks (line 664). Rank under a stronger weight is in Continuing Work.

### Q-09. Why point-Laplace priors that give signed contrasts instead of a nonnegative factorization like DeSurv?

**Suggested answer.** The point-exponential prior on nonnegative input lowered training
concordance in every input setup (0.573 to 0.635 against 0.658 at rank 7, 0.519 against 0.655
at rank 5) and at most one of five folds converged. On the divided-by-SD and rank inputs its
programs tracked cohort membership (correlation 0.998 to 1.000, Table `tab:lf-prior-real`). The
cost is interpretation: each program contrasts two poles (line 493).

### Q-10. Each validation cohort is standardized with its own mean and SD, and missing genes are set to zero. How does that affect the scores?

**Suggested answer.** Within-cohort standardization removes platform and batch location and scale
differences, which is necessary across RNA-seq and arrays. Concordance within a series is
unaffected by a shift and nearly unaffected by a common scale. The calibration slope depends on
the scale, which is one reason slopes vary from 0.35 to 1.85. A missing gene contributes zero to
every patient's projection, which lowers the score's variance in that series. Have the per-series
count of missing genes ready (from the external scoring inputs). It is not in the chapter.

### Q-11. Why one baseline hazard for TCGA-PAAD and CPTAC?

**Suggested answer.** Stated as a limitation (line 497). The fit stratified by cohort is the first
open priority in PROJECT_STATUS.qmd. The cohort-strata variant appears only in the
point-exponential comparison (Table `tab:lf-prior-real`). Concordance is evaluated within cohort,
which is unaffected by a cohort-level baseline difference.

### Q-12. Why not compare with supervised principal components, Cox partial least squares, or a two-step fit with a ridge Cox model?

**Suggested answer.** The comparators were chosen as the two-step pipeline the method replaces,
the most closely related supervised factorization (DeSurv), a strong gene-level baseline
(elastic net, which selected the ridge penalty) and a clinical classifier (PurIST). On these data
the joint-minus-two-step difference comes entirely from coefficient shrinkage, and a ridge Cox
two-step would likely match the joint fit. That comparison and supervised principal components
[@Bair2004] are reasonable additions. Neither was run.

### Q-13. You drop part of the Hessian diagonal. Is that accurate in your strong-signal simulations?

**Suggested answer.** The omitted term, $-\hat\theta_i^2\sum_{s\le t_i}d(s)/\hat S(s)^2$, is small
when risk sets are large and the risks are similar. With a true linear predictor SD of 9 to 26,
risk sets late in follow-up are dominated by few patients, and the term may not be small.
Overstated weights $W_{ii}$ strengthen the ridge-type term in the coefficient update and could add
to the attenuation. This has not been checked. A direct test is to refit a few strong-survival
datasets with the exact diagonal.

### Q-14. Program 4's acinar pole is adverse, but Young (2026) found exocrine content contributed little to survival. How do you read Program 4?

**Suggested answer.** Program 4 contrasts acinar enzyme genes with B-cell and interferon genes
(Table `tab:programs`, line 434). A higher score means more acinar and less immune expression.
The adverse direction is likely carried by the immune pole: lower B-cell and interferon
expression is a common poor-prognosis signal in PDAC, and Program 4 correlates negatively with
the Bailey immunogenic signature. Acinar content mostly reflects normal pancreas in the sample.
The label describes genes, not tumor biology (line 434). Program 4 is survival-active but not a
main program (three of five folds).

### Q-15. Methylation was modeled as beta values. Why not M-values? And why does point-exponential expression beat point-Laplace on the 144 patients?

**Suggested answer.** Beta values were filtered to sites with mean between 0.1 and 0.9, where the
variance is less dependent on the mean. M-values (Du et al. 2010, BMC Bioinformatics 11:587) are
the usual choice for Gaussian models and are worth checking in the joint model. On the second
point: on 144 TCGA patients the point-exponential expression fit gave 0.620 against 0.582 for
point-Laplace, opposite to the 273-patient result. The 144-patient comparison has about half the
events, one cohort, and SEs of 0.03 to 0.045. The main prior choice rests on the 273-patient
comparison with both cohorts and the simulation.

### Q-16. What will the joint expression and methylation model show if it works?

**Suggested answer.** Mainly interpretation, not prediction (line 200, D10): which genes and CpG
sites load together, whether methylation adds programs the expression fit lacks, and whether the
joint fit recovers the immune methylation program seen at 100,000 sites ($R^2$ 0.06 on the five
expression programs, line 452). The current fit is dominated by methylation (29% to 97% of
shared-score precision), and the proposed fixes are weighting methylation only in the score
update or summarizing CpG sites by region (line 206).

### Q-17. SafeBayes and Holmes-Walker calibrate for robustness and usually give weights near or below 1. How would they solve a problem that needs $w_S$ near $10^4$?

**Suggested answer.** They would not by themselves. They target the right amount of posterior
uncertainty when the model is misspecified, and SafeBayes never exceeds 1. The point of including
them is the coefficient calibration problem (C-01, C-23), where a calibrated power on the partial
likelihood could correct overconfident or underconfident intervals. Moving the loadings is a
separate problem that the normalized-loss approach (DeSurv-style) and the cross-validated weight
address. The comparison in Continuing Work keeps the two goals separate.

### Q-18. The risk-group cutpoint maximizes the log-rank statistic. Is the high-risk group at the 20% boundary?

**Suggested answer.** The joint cutpoint of 0.618 training SD put 59 of 273 training patients
(21.6%) in the high-risk group, just above the 20% minimum (`prompt7_desurv_style/cutpoints.csv`).
It was chosen on out-of-fold scores and applied unchanged externally, which keeps the external
hazard ratios honest even though maximally selected statistics (Lausen and Schumacher 1992,
Biometrics 48:73) overstate the training contrast. External HRs: 1.69 (1.34 to 2.14) pooled.

### Q-19. Are there really no Bayesian outcome-guided factor models for survival?

**Suggested answer.** The chapter's claim is narrower: no empirical Bayes factorization that fits
several programs with their Cox associations and reports posterior uncertainty (line 55). Bayesian
sparse factor regression with survival outcomes exists in earlier work by West and colleagues on
gene expression, and integrative Bayesian models such as BIP (Chekouo and Safo 2023, Biostatistics
24:124) are framed for general clinical outcomes. Line 53 cites recent supervised Bayesian
factorizations for continuous, binary, ordinal and categorical outcomes. If pressed, concede that
the novelty is the empirical Bayes priors, the projection risk score and the survival-share
analysis, not survival supervision in a Bayesian factor model in general. Andrew should check
whether BIP supports survival outcomes before the defense.

### Q-20. External concordance is 0.630 against 0.614 for DeSurv, with intervals that include zero. What is the claim?

**Suggested answer.** That SSBMF at rank 5 predicts as well as DeSurv, elastic-net Cox and PurIST
trained on the same patients (all paired intervals include zero), and better than larger-rank
EBMF with unpenalized Cox (differences 0.055 and 0.049, intervals excluding zero). The claim is
comparable prediction with interpretable programs and a single fixed score, not superiority.
DeSurv's slopes stay below 1 (0.45 to 0.67) and SSBMF's range from 0.35 to 1.85.

### Q-21. Do the posterior intervals for the coefficients cover the true effects?

**Suggested answer.** Not yet. The Discussion names interval coverage as part of the open problem, and it is not yet planned work. In the
fixed-weight simulation at $w_S=1$ the 95% interval for the planted coefficient contained the
true value in none of the 30 fits with 2,000 genes and in 6 of the 30 with 50 genes, because the
posterior means were shrunk to less than half the true value
(`docs/findings/2026-11-04-survival-weight-p2000-findings.md`, "Per-weight curves"). The cause is
the coefficient update: each coefficient's Normal prior variance is estimated from its single
pseudo-observation, $\hat\sigma_k^2=\max(0,x_k^2-s_k^2)$, and the posterior SD ignores the
uncertainty in that estimate. The strong-survival simulations are far stronger than the PDAC data
(true-predictor SD 9 to 26 against coefficients of at most 0.42), and coverage was measured only
at that strength.
