# Architectural Alignment & Resolution Blueprint: `multiomicsGEP`

**Purpose of this Document:** This specification consolidates all 13 mathematical diagnostics, root-cause analyses, and final architectural decisions for the `multiomicsGEP` project (Supervised Empirical Bayes Matrix Factorization + Cox Proportional Hazards). It serves as the authoritative reference for Claude Code to exit the incremental patch loop (`Task1`–`Task5`), unify the single-modality and multimodal codebases, and finalize the method for publication.

---

## 1. Project North Star & Non-Negotiable Requirements

1. **Joint Bayesian NMF + Projection Cox Model ($\eta = YF\beta$):**
   - Decompose non-negative molecular data $Y_m \approx L F_m^\top + E_m$ into non-negative patient activities ($L \ge 0$) and sparse, non-negative feature loadings ($F_m \ge 0$) using Empirical Bayes Normal Means (`ebnm`), jointly coupled to a Cox proportional hazards model.
   - Parameterize the Cox linear predictor as the feature projection $\eta = \tilde{Z}\beta$ (rather than $L\beta$) so that out-of-sample validation cohorts are scored in closed form ($\hat{\eta}_{\text{new}} = \tilde{Z}_{\text{new}} \bar{\beta}$) without running a test-time optimization loop.
2. **True Survival Supervision of Feature Loadings ($F_m$):**
   - Survival information must actively guide the learning of $F_m$ ($\alpha > 0$) to recover **low-variance prognostic programs** masked by high-variance, outcome-neutral background signals (such as tumor purity or exocrine composition in PDAC).
   - Setting $\alpha = 0$ ($\alpha_F = 0$) is only used as an unsupervised/two-step ablation control, not as the final proposed method.
3. **Principled Rank Selection ($K$) & Role Classification:**
   - Select model complexity using training data only (5-fold CV C-index with the 1-SE parsimony rule).
   - Separate fitted factors into **Survival-Active**, **Genomics-Only (Reconstruction-Active)**, and **Inactive/Dead** roles without destructive post-hoc deletion refits.
4. **Multi-Cohort & Multi-Modal Generalization:**
   - Accommodate pooled training cohorts (TCGA + CPTAC) and matched multi-omics modalities (Gene Expression + DNA Methylation) while remaining numerically stable across both focused feature sets (~2,000–3,000 genes + ~10,000 CpGs) and full feature universes (~17,900 genes + ~305,000 CpGs).

---

## 2. What Must Be Preserved As-Is (Do Not Break)

- **The $\eta = YF\beta$ Projection Contract:** Keep frozen training loading norms $d_{mk}^*$, coefficients $\bar{\beta}$, and orientation signs so external prediction is 100% self-contained.
- **1-Step Unsupervised Feature Pre-Screening (`DeSurv` Rule):** Keep the combined mean + variance rank screening (top 3,000 genes per training cohort before intersection) as the primary pipeline.
- **Three-Way Factor Classification (`classify_factors()`):** Retain the classification of fitted factors into `survival_active`, `genomics_only`, and `dead` on the fitted $K_{\text{init}}$ model.
- **Cohort Handling Strategy:** Keep within-cohort feature scaling, optional cohort-stratified Cox risk sets (`strata_id`), and partially pooled cohort survival coefficients (`beta_cohort_id` as a heterogeneity diagnostic). Do *not* re-enable the $\{0,1\}$ dummy column in $L$ (`cohort_id`) as a default.
- **5-Method Biological Characterization & TDD Suite:** Preserve the `fgsea`, ORA, PurIST subtype concordance, DeSurv gene overlap, external HR pipelines, and existing unit tests.

---

## 3. Complete 13-Point Issue / Diagnosis / Solution Matrix

### Domain A: Survival Supervision of Feature Loadings ($F$)

#### Issue 1: Survival supervision is off in Single-Modality ($\alpha_F = 0$) and swamped ($10^{-7}$) in Multimodal
- **Diagnostic / Root Cause:** In the single-modality benchmark, $\alpha_F = 0$ shuts off survival in the $F$ update completely (making $L$ and $F$ purely unsupervised EBMF). In the multimodal model, reconstruction and survival terms are added with equal 1:1 weight (`w_S = 1`), but because Gaussian reconstruction sums over $N \times P_m$ entries while Cox survival sums over $N$ patients, the median survival fraction of loading precision at terminal state is **$2.03 \times 10^{-7}$ for expression and $9.00 \times 10^{-10}$ for methylation**. Consequently, the model recovers 0% (0/27) of planted low-variance prognostic programs at lower noise.
- **Effective Solution:** Introduce an explicit likelihood-balancing multiplier $\gamma$ on the Cox working precision $w_i$ (analogous to `DeSurv`'s loss scaling):
  $$\gamma = \frac{\alpha}{1 - \alpha} \cdot \frac{\sum_{m=1}^M v_m p_m}{N_{\text{event}}}$$
  where $v_m = \bar{p} / p_m$ is the modality dimension weight and $\alpha \in [0, 1)$ is the supervision strength hyperparameter. When $\alpha \in (0, 1)$, the survival precision in the $F$ update scales to $O(1)$ relative to reconstruction.

#### Issue 2: Uncoupled Jacobi $F$ updates caused RMSE blowup and positive-feedback runaway
- **Diagnostic / Root Cause:** The earlier single-modality survival branch (Appendix A.7, Eq. 18) omitted the within-factor feature coupling term $-\overline{\tilde{\beta}_k^2}\sum_i w_i \frac{y_{ij}}{d_k}\bar{C}_i$ and updated all $p$ genes simultaneously (Jacobi-style). Every correlated gene regressed the full factor residual onto itself at once, overshooting $\eta$ by $O(p)$ and triggering a runaway bilinear feedback loop between $L$ and $F$ (`RMSE` jumping from ~290 to ~790).
- **Effective Solution:** Use the exact coupled coordinate update (sequential Gauss-Seidel sweep over features $j = 1, \dots, p_m$) in both single-modality and multimodal fits, refreshing the partial projection $\bar{\tilde{Z}}_{\cdot k}$ and $\text{Var}(\tilde{Z}_{\cdot k})$ after each feature update, paired with precision damping $\rho \in [0.5, 0.7]$ when $\alpha > 0$.

#### Issue 3: The "Chicken-and-Egg" Gating Trap ($\bar{\beta}_k \approx 0$ prevents $F_{\cdot k}$ from ever finding a low-variance prognostic direction)
- **Diagnostic / Root Cause:** In the $F$ update, the survival gradient is multiplied by $\bar{\beta}_k$ and $\overline{\beta_k^2}$. When $L$ and $F$ are initialized purely from the top $K$ SVD components of $Y$, a low-variance prognostic program lies outside the initial span of $F$. Its initial $\bar{\beta}_k \approx 0$, so its survival gradient is zero and no factor ever rotates toward it.
- **Effective Solution:** Use a **non-negative hybrid initialization**: initialize the first $K-2$ factors from non-negative SVD/NMF of $Y_m$, and seed the last 2 factor columns directly from the positive (adverse) and negative (protective) parts of the univariate Cox score vector:
  $$F_{m, \cdot, K-1}^{(0)} = \max\left(0, +\tilde{Y}_m^\top u^{(0)}\right), \quad F_{m, \cdot, K}^{(0)} = \max\left(0, -\tilde{Y}_m^\top u^{(0)}\right)$$
  followed by a brief $\beta$ warm-up (e.g., 5–10 iterations) so both adverse and protective prognostic directions start with $\bar{\beta}_k^{(0)} \neq 0$ while strictly respecting $F_m \ge 0$.

---

### Domain B: Survival Coefficient ($\beta$) Priors & Collapse to Zero

#### Issue 4: Permanent $\beta \to 0$ collapse on real TCGA multimodal data (`M4` at outer sweep 99) and inside CV folds (`Task2` Fold 1)
- **Diagnostic / Root Cause:** Both $\beta$ prior estimators in the current codebase contain an **absorbing point-mass zero state**:
  1. The per-factor Normal empirical Bayes variance update $\hat{\sigma}_k^2 = \max(0, x_{\beta,k}^2 - s_{\beta,k}^2)$ sets $\hat{\sigma}_k^2 = 0$ (point mass at zero) the instant a factor's survival Wald statistic dips to $|x_{\beta,k} / s_{\beta,k}| \le 1$.
  2. The shared Point-Laplace prior (`shared_point_laplace`), estimated via unconstrained maximum marginal likelihood across only $K = 7$ scalar coordinates, hits slab probability $\hat{\pi}_\beta = 0$ whenever 6 factors are null and the 1 prognostic factor has moderate real-data SNR (especially when $A_{\beta,k}$ is diluted by $\text{Var}(Z_{ik})$ across 13,000 features). Once $\hat{\pi}_\beta = 0$, all $\bar{\beta}_k = 0$ and $\overline{\beta_k^2} = 0$ permanently, shutting off supervision of $F$ for all remaining sweeps and producing constant risk scores (`C = NA`).
- **Effective Solution:** Remove the absorbing zero state from the $\beta$ prior during CAVI iterations by either:
  1. **Constrained Spike-and-Slab (`ebnm_point_normal` / `ebnm_point_laplace`):** Enforce a lower bound on the active slab proportion $\pi_\beta \in [\frac{1}{K}, 0.50]$ (and a minimum slab scale $\sigma_{\text{slab}} \ge 0.05$) so `ebnm` can never set $\hat{\pi}_\beta = 0$ globally across all factors, **OR**
  2. **Inverse-Gamma Regularized Normal Prior:** Update per-factor prior variance via $\hat{\sigma}_k^2 = \frac{\beta_0 + \frac{1}{2}(\bar{\beta}_k^2 + \text{Var}(\beta_k))}{\alpha_0 + 1/2} > 0$ (with small hyper-parameters, e.g., $\alpha_0 = 1, \beta_0 = 0.05$), guaranteeing $\hat{\sigma}_k^2 > 0$ at every sweep while strongly shrinking null coefficients.

---

### Domain C: CAVI Convergence, Scale Canonicalization, & Bookkeeping

#### Issue 5: Multimodal fits hit the 500-sweep cap without converging and oscillate by ~1.9% (`0.01915` vs. `1e-4` tolerance)
- **Diagnostic / Root Cause:** Three interacting bookkeeping conflicts prevent CAVI from ascending a coherent ELBO:
  1. *Data-dependent scale canonicalization:* Rescaling each factor by $c_k = \text{sd}(\bar{Z}_{\cdot k})$ after every factor update depends on the empirical feature covariance $\bar{F}_{mk}^\top \widehat{\text{Cov}}(Y_m) \bar{F}_{mk}$. As features update sequentially, $c_k$ shifts continuously, rescaling $L_{\cdot k}$ and altering the next sweep's reconstruction precision $A_{mjk}^{\text{rec}} = \tau_{mj}\sum_i \overline{l_{ik}^2}$.
  2. *Untransformed $\beta$ priors under rescaling:* Rescaling factor $k$ by $c_k$ maps $\bar{\beta}_k \mapsto c_k \bar{\beta}_k$, but the per-factor Normal prior variance $\sigma_k^2$ is not carried through canonicalization, and a single shared prior $g_\beta$ *cannot* be transformed by $K$ different scale factors $c_1, \dots, c_K$ simultaneously.
  3. *Stale KL caching & normalizer bug:* In both $L$ and $F_m$, the posterior $q$ is updated under the incoming prior, and then the prior $g$ is refitted by marginal likelihood without updating the cached KL divergence. Additionally, the single-modality monitor $J$ counts the Gaussian log-normalizer $\frac{1}{2}\log(A_{L,k}/2\pi)$ once instead of $n$ times (Appendix A.9).
- **Effective Solution:**
  1. **Euclidean Norm Scaling (Completely Decoupled from $\beta$):**
     - After updating $q(L_{\cdot k})$, normalize $L_{\cdot k}$ to unit root-mean-square second moment: let $c_k = \sqrt{\frac{1}{n}\sum_{i=1}^n \overline{l_{ik}^2}}$ and apply the exact pushforward map between $L$ and $F_m$:
       $$\bar{l}_{ik} \leftarrow \bar{l}_{ik} / c_k, \quad \overline{l_{ik}^2} \leftarrow \overline{l_{ik}^2} / c_k^2, \quad \lambda_{L,k} \leftarrow c_k \lambda_{L,k}$$
       $$\bar{f}_{mjk} \leftarrow c_k \bar{f}_{mjk}, \quad \overline{f_{mjk}^2} \leftarrow c_k^2 \overline{f_{mjk}^2}, \quad \lambda_{F,mk} \leftarrow \lambda_{F,mk} / c_k$$
     - In the Cox projection, normalize each modality loading column by its start-of-sweep Euclidean root-mean-square norm $d_{mk} = \sqrt{\frac{1}{p_m}\sum_{j=1}^{p_m} \bar{f}_{mjk}^2 + 10^{-10}}$. Because $F_{mk} / d_{mk}$ is invariant to $c_k$, **$\tilde{Z}_{ik}$, $\beta_k$, and the $\beta$ prior $g_\beta$ never need to be rescaled when $L$ and $F_m$ are canonicalized!**
  2. **Atomic `ebnm` KL Computation:** Compute the exact KL divergence immediately from the `ebnm` return object at each update (`res$log_likelihood` minus the expected Gaussian log-likelihood with the correct $n$ or $p_m$ multiplier), ensuring the posterior $q$ and prior $\hat{g}$ are always evaluated at the exact same state.

---

### Domain D: Rank Selection ($K$), Prior Sparsity, & Factor Pruning

#### Issue 6: Pure Exponential priors never prune factors or features, and Unsupervised EBMF saturates at $K = 20\text{–}40$
- **Diagnostic / Root Cause:** Pure Exponential priors (`ebnm_exponential`) have density $\lambda e^{-\lambda \theta}$ on $(0, \infty)$ with no point mass at zero, so posterior means are strictly positive and never zero out noisy features or redundant factors. Furthermore, in high-dimensional genomics, structured biological/technical variance exceeds i.i.d. Gaussian noise $\tau_j^{-1}$ across dozens of directions, so unsupervised EBMF uses every factor offered up to $K = 40$.
- **Effective Solution:** Use **Point-Exponential priors (`ebnm_point_exponential`)** on both $L$ and $F_m$ to provide an exact spike at zero ($\delta_0$) for sparse gene/CpG program selection and dead-column pruning, paired with a spike-and-slab prior on $\beta$ to distinguish prognostic factors ($K_{\text{surv}}$) from background reconstruction factors ($K_{\text{genomics}}$). Select $(K_{\text{init}}, \alpha)$ via 5-fold training CV C-index with the 1-SE parsimony rule.

#### Issue 7: Post-hoc factor deletion and refitting (`prune_nullcheck`) either deletes zero factors on real multimodal data OR overwrites the planted prognostic program in single-modality simulations
- **Diagnostic / Root Cause:**
  1. *On real multimodal data:* `prune_nullcheck` evaluates candidate deletions holding the other $K-1$ factors frozen (Eq. 36–37). Because reconstruction sums over $144 \times 13,000$ entries ($\approx 610,423$) while survival sums over 75 deaths ($\approx -304$), zeroing any factor while holding others frozen causes a massive drop in the reconstruction score, so zero factors are ever deleted (`Final rank = Starting rank`).
  2. *In single-modality simulations:* A low-variance prognostic factor's $P$-feature KL penalty in $F$ exceeds its unweighted $N$-sample Cox gain, causing unprotected pruning to delete it. When protected pruning deletes a high-variance background factor instead and refits at $K-1$, the protected prognostic factor rotates during the refit to absorb the deleted factor's genomic variance (`|cor|` drops from 0.896 at $K=7$ to 0.391 at $K=6$).
- **Effective Solution:** **Disable post-loop factor deletion and refitting (`prune_nullcheck`) entirely.** Retain the converged $K_{\text{init}}$ fit and use `classify_factors()` to report **Survival-Active**, **Genomics-Only**, and **Dead/Inactive** factors without refitting a reduced rank.

---

### Domain E: Multi-Modal Integration, Feature Scaling, & Pre-Filtering

#### Issue 8: Expanding from 10,000 to 25,000 CpGs caused Cox risk-set underflow at Sweep 6 and drove predictor $\text{SD}(\hat{\eta})$ to $39,086$
- **Diagnostic / Root Cause:** In the multimodal model, $Z_{ik} = \sum_{m,j} \tilde{y}_{mij} f_{mjk}$ and $A_{L,k} = \sum_{m,j} \tau_{mj}\overline{f_{mjk}^2}$ are unnormalized sums over $p_m$ features, and the multimodal Cox working curvature $w_i$ lacked a lower floor. Increasing $p_{\text{meth}}$ from 10,000 to 25,000 inflated $\hat{\eta}_i$ across 5,734 units mid-sweep, underflowing standard $\exp(\hat{\eta}_i)$ risk sets and causing 10,000+ CpGs to dominate 3,000 genes in the shared activity update $L$.
- **Effective Solution:**
  1. Weight each modality $m$'s reconstruction precision by $v_m = \frac{\bar{p}}{p_m}$ and scale each modality's projection by $s_{mk} = \frac{1}{M \sqrt{p_m} d_{mk}}$ (where $d_{mk} = \sqrt{\frac{1}{p_m}\|F_{mk}\|_2^2 + 10^{-10}}$) so that $A_{L,k}$ and $\text{Var}(\tilde{Z}_{ik})$ remain strictly $O(1)$ and invariant to whether $p_{\text{meth}}$ is 10,000, 25,000, or 305,000.
  2. Make log-sum-exp Cox risk-set arithmetic, Breslow tied-event-time handling, and the curvature floor $w_i = \max(a_i, 10^{-6})$ permanent defaults across all models.

#### Issue 9: Architectural split between Single-Modality (z-scored data + non-negative $F$) and Multimodal (unstandardized data + intercepts + signed $F_m$)
- **Diagnostic / Root Cause:** Single-modality z-standardized $Y$ (introducing negative values) without an intercept while constraining $L, F \ge 0$ ($LF^\top \ge 0$ cannot reconstruct negative values). Conversely, multimodal added a feature intercept $\mu_{mj}$ and switched $F_m$ to signed Point-Laplace priors—abandoning parts-based Non-Negative Matrix Factorization on $F_m$.
- **Effective Solution:** Unify into **one codebase** where single-modality is simply the $M = 1$ case of the $M$-modality model:
  - **Non-Negative Reconstruction Input ($Y_m \ge 0$):** Within each cohort, divide each non-negative feature $X_{m, \cdot j} \ge 0$ (where $X_{\text{expr}} = \log_2(\text{TPM}+1)$ and $X_{\text{meth}} = \beta_{\text{CpG}} \in [0, 1]$) by its within-cohort standard deviation $\text{sd}(X_{m,\cdot j})$ *without* subtracting the mean:
    $$Y_{m,ij} = \frac{X_{m,ij}}{\text{sd}_{\text{cohort}}(X_{m,\cdot j})} \ge 0$$
    This harmonizes feature scales across cohorts while keeping $Y_{m,ij} \ge 0$ so **both $L$ and $F_m$ use sparse non-negative Point-Exponential priors (`ebnm_point_exponential`)**.
  - **Centered Cox Projection ($\tilde{Y}_m$):** Define the column-centered matrix $\tilde{Y}_{m,ij} = Y_{m,ij} - \bar{Y}_{m,\cdot j}$ (which has mean 0 and variance 1 in every cohort!) and form the Cox projection $\tilde{Z}$ using $\tilde{Y}_m$. This preserves pure additive NMF ($Y_m \approx L F_m^\top \ge 0$) in reconstruction while giving the Cox model zero-mean, unit-variance-scaled inputs without intercept drift.

#### Issue 10: Dense feature projections dilute survival signal at large $p$, while 2-step supervised pre-filtering undermines the methodological thesis
- **Diagnostic / Root Cause:** As demonstrated in `DeSurv` SI Table S2, projecting $Z = YF$ across thousands of dense background features accumulates noise variance that drowns out the ~250 true marker genes. Conversely, pre-filtering the input matrix $Y$ by univariate Cox survival p-values before running NMF undercuts the core thesis that the joint factorization itself discovers low-variance prognostic programs amidst outcome-neutral background genes.
- **Effective Solution:** Keep **1-step unsupervised feature screening** (~2,000–3,000 genes by combined mean + variance rank and ~10,000 CpGs by variance/promoter mapping) as the primary model, report full-universe sensitivity curves (up to 17,900 genes and 305k CpGs) in the Supplement, and rely on `ebnm_point_exponential` sparsity on $F_m$ (plus optional top-$N$ specificity focusing) to zero out noisy background features during projection. Only evaluate 2-step supervised pre-filtering as an ablation comparator.

---

### Domain F: Multi-Cohort Modeling, Evaluation Orientation, & Repo Hygiene

#### Issue 11: Pooling TCGA + CPTAC without per-cohort normalization collapsed $\beta \to 0$, whereas adding a dummy $\{0,1\}$ cohort indicator to $L$ hurt external C-index
- **Diagnostic / Root Cause:** Without per-cohort scaling, the top SVD component of merged TCGA + CPTAC captures between-cohort technical scale differences rather than biology, driving $\beta \to 0$. Once per-cohort scaling is applied, appending a dummy $\{0,1\}$ cohort column to $L$ (`cohort_id`) wastes a factor degree of freedom and triggers CAVI factor-collapse on fresh SVD (`C` drops from 0.627 to 0.543, and `K_survival_active` drops from 2 to 1).
- **Effective Solution:** Scale features within each cohort prior to pooling (and apply identical within-cohort scaling to validation cohorts), keep `cohort_id = NULL` (no dummy column in $L$) as the default, and handle cohort heterogeneity via **cohort-stratified Cox partial likelihood risk sets (`strata_id`)**, **cohort-specific residual precisions ($\tau_{mjc}$)**, and **partially pooled cohort survival coefficients (`beta_cohort_id`)** as a heterogeneity diagnostic.

#### Issue 12: Inverted `concordance()` sign check in Phase C and circular `max(C, 1-C)` test-set orientation
- **Diagnostic / Root Cause:** `survival::concordance(Surv(time, status) ~ eta)` defaults to assuming a higher predictor indicates *longer* survival (opposite of Cox hazard risk), so omitting `reverse = TRUE` inverted the Phase C training sign-flip check. Furthermore, when Phase C flipped `EBeta`, it failed to flip `EBeta_pooled` and overwrote `EBeta2` with `(-EBeta)^2` (discarding posterior variance). Earlier benchmark scripts masked orientation mismatches by computing `max(C, 1-C)` on test sets, which uses validation outcomes to orient the predictor and biases null C-indices upward to ~0.54.
- **Effective Solution:** Fix risk orientation **once** from training data only using `survival::concordance(..., reverse = TRUE)`, preserve posterior variance `EBeta2 = Var(beta) + EBeta^2` and coherent `EBeta_pooled` signs across flips, and score all external validation cohorts using frozen training orientation (`frozen_reverse_cindex()`).

#### Issue 13: Historical double $\log_2(x+1)$ bug, repeated PACA-AU validation donors, unannotated CpGs, and CSV trace bloat
- **Diagnostic / Root Cause:**
  1. The internal data loader applies $\log_2(x+1)$ to TCGA, Dijk, and PACA-AU RNA-seq, and the benchmark preprocessing wrapper applies $\log_2(x+1)$ a second time under its RNA-seq flag (while CPTAC is logged once).
  2. Summing the 5 external validation cohorts (616 records) double-counts 47 PACA-AU patients present in both the microarray and RNA-seq series plus 1 duplicate RNA-seq donor (568 unique individuals), invalidating naive pooled p-values and concatenated bootstrap CIs.
  3. Multimodal programs lack CpG-to-gene promoter annotation, and 82.3% (144,461 lines) of the 176,342-line branch diff consists of per-iteration CSV trace dumps.
- **Effective Solution:**
  1. Remove the duplicate $\log_2(x+1)$ transformation in the preprocessing wrapper so all RNA-seq cohorts are logged once.
  2. Deduplicate the 568 unique external validation donors for pooled KM/Cox analyses and use cohort-stratified bootstrap resampling (`bootstrap_concordance_diff_ci_stratified()`).
  3. Map top-loading multimodal CpGs (`cgXXXXXXXX`) to gene promoters/enhancers via the Illumina 450K manifest to establish cis-epigenetic regulation of transcriptional programs (especially in Immune Program 2 and Basal-vs-Classical Program 4), and archive/compress the diagnostic CSV trace folders.

---

## 4. Exact Mathematical Equations for the Unified CAVI Engine

To ensure zero ambiguity during implementation, the unified $M$-modality CAVI updates ($M=1$ for single-modality; $M=2$ for Expression + Methylation) are specified below:

### 4.1 Definitions & Scaling Constants
- **Modality Reconstruction Weight:** $v_m = \frac{\bar{p}}{p_m}$ (where $\bar{p} = \frac{1}{M}\sum_{m=1}^M p_m$, or $v_1 = 1$ when $M = 1$).
- **Projection Scaling Constant (Fixed at Start of Outer Sweep):**
  $$d_{mk} = \sqrt{\frac{1}{p_m}\sum_{j=1}^{p_m} \bar{f}_{mjk}^2 + 10^{-10}}, \quad s_{mk} = \frac{1}{M \sqrt{p_m} \, d_{mk}}$$
  Define the scaled, centered design entry $\tilde{a}_{mijk} = s_{mk} \tilde{y}_{mij}$, so the projected score for patient $i$ on factor $k$ is:
  $$\tilde{Z}_{ik} = \sum_{m=1}^M \sum_{j=1}^{p_m} \tilde{a}_{mijk} f_{mjk}, \quad \bar{\tilde{Z}}_{ik} = \sum_{m,j} \tilde{a}_{mijk} \bar{f}_{mjk}, \quad \text{Var}_q(\tilde{Z}_{ik}) = \sum_{m,j} \tilde{a}_{mijk}^2 (\overline{f_{mjk}^2} - \bar{f}_{mjk}^2)$$
- **Tempered Cox Working Quantities (Updated at Start of Outer Sweep):**
  Using log-sum-exp Breslow risk sets (within `strata_id` if supplied), compute score $u_i = \delta_i - a_i$, floored curvature $w_i = \max(a_i, 10^{-6})$, and scaled survival weight:
  $$\gamma = \frac{\alpha}{1 - \alpha} \cdot \frac{\sum_{m=1}^M v_m p_m}{N_{\text{event}}}, \quad \tilde{w}_i = \gamma w_i, \quad \tilde{h}_i = \gamma (u_i + w_i \hat{\eta}_i)$$

### 4.2 Coordinate Updates (For $k = 1, \dots, K$)
1. **Shared Patient Activities $q(L_{\cdot k})$ (`ebnm_point_exponential`):**
   $$A_{L,k} = \sum_{m=1}^M v_m \sum_{j=1}^{p_m} \tau_{mj} \overline{f_{mjk}^2}, \quad B_{L,ik} = \sum_{m=1}^M v_m \sum_{j=1}^{p_m} \tau_{mj} R_{mij}^{(-k)} \bar{f}_{mjk}$$
   Pass $x_{L,ik} = B_{L,ik} / A_{L,k}$ and $s_{L,k} = A_{L,k}^{-1/2}$ to `ebnm_point_exponential`. Then apply the Euclidean norm canonicalization $c_k = \sqrt{\frac{1}{n}\sum_{i=1}^n \overline{l_{ik}^2}}$ to $(L_{\cdot k}, F_{mk})$ and their prior rates $(\lambda_{L,k}, \lambda_{F,mk})$.
2. **Modality Feature Loadings $q(F_{mk})$ (`ebnm_point_exponential`, Sequential Gauss-Seidel over $j = 1, \dots, p_m$):**
   Let $\bar{V}_{i}^{(-k)} = \sum_{l \neq k} \bar{\tilde{Z}}_{il} \bar{\beta}_l$ and $\bar{U}_{ik}^{(-mj)} = \bar{\tilde{Z}}_{ik} - \tilde{a}_{mijk}\bar{f}_{mjk}$. The combined precision and linear signal for feature $j$ of modality $m$ are:
   $$A_{F,mjk} = v_m \tau_{mj} \sum_{i=1}^n \overline{l_{ik}^2} + \overline{\beta_k^2} \sum_{i=1}^n \tilde{w}_i \tilde{a}_{mijk}^2$$
   $$B_{F,mjk} = v_m \tau_{mj} \sum_{i=1}^n R_{mij}^{(-k)} \bar{l}_{ik} + \bar{\beta}_k \sum_{i=1}^n \tilde{a}_{mijk}(\tilde{h}_i - \tilde{w}_i \bar{V}_i^{(-k)}) - \overline{\beta_k^2} \sum_{i=1}^n \tilde{w}_i \tilde{a}_{mijk} \bar{U}_{ik}^{(-mj)}$$
   *If using precision damping $\rho \in (0, 1]$ when $\alpha > 0$:* blend the survival precision/signal by $\rho$ before calling `ebnm_point_exponential` so the returned posterior and `res$log_likelihood` remain in exact KL sync. Immediately refresh $\bar{\tilde{Z}}_{\cdot k}$ and $\text{Var}_q(\tilde{Z}_{\cdot k})$ after each feature $j$.
3. **Survival Coefficients $q(\beta_k)$ (Non-Collapsing Empirical Bayes Prior):**
   $$A_{\beta,k} = \sum_{i=1}^n \tilde{w}_i \left( \bar{\tilde{Z}}_{ik}^2 + \text{Var}_q(\tilde{Z}_{ik}) \right), \quad B_{\beta,k} = \sum_{i=1}^n \bar{\tilde{Z}}_{ik} (\tilde{h}_i - \tilde{w}_i \bar{V}_i^{(-k)})$$
   Update $\beta_k$ via the non-collapsing prior (either `ebnm_point_normal`/`ebnm_point_laplace` with lower-bounded slab weight $\pi_\beta \ge 1/K$, or Normal with Inverse-Gamma variance $\hat{\sigma}_k^2 = \frac{\beta_0 + \frac{1}{2}\overline{\beta_k^2}}{\alpha_0 + 1/2} > 0$).
4. **Feature Residual Precisions $\tau_{mj}$:**
   $$\hat{\tau}_{mj} = \frac{n}{\max(Q_{mj}, \, n \times 10^{-8})}, \quad Q_{mj} = \sum_{i=1}^n \left[ \left(y_{mij} - \sum_{k=1}^K \bar{l}_{ik}\bar{f}_{mjk}\right)^2 + \sum_{k=1}^K \left(\overline{l_{ik}^2}\,\overline{f_{mjk}^2} - \bar{l}_{ik}^2 \bar{f}_{mjk}^2\right) \right]$$

---

## 5. Summary Specification Table

| Model Element | Target Implementation Specification |
| :--- | :--- |
| **Input Data ($Y_m$)** | Strictly non-negative $Y_{m,ij} = X_{m,ij} / \text{sd}_{\text{cohort}}(X_{m,\cdot j}) \ge 0$ (single $\log_2(\text{TPM}+1)$ for RNA-seq, $\beta$-values for methylation, scaled by within-cohort feature SD without subtracting the mean). |
| **Reconstruction Model** | $Y_m = L F_m^\top + E_m$, $(E_m)_{ij} \sim \mathcal{N}(0, \tau_{mj}^{-1})$ with modality weight $v_m = \bar{p} / p_m$. |
| **Priors on $L$ and $F_m$** | Independent `ebnm_point_exponential` priors per factor $k$ on $L_{\cdot k} \ge 0$ and per modality/factor $(m,k)$ on $F_{mk} \ge 0$. |
| **Cox Linear Predictor ($\eta$)** | $\eta_i = \sum_{k=1}^K \tilde{Z}_{ik} \beta_k$, where $\tilde{Z}_{ik} = \sum_{m=1}^M \frac{1}{M\sqrt{p_m}} \left(\tilde{Y}_m \frac{F_{mk}}{d_{mk}}\right)_i$, $\tilde{Y}_m = Y_m - \mathbf{1}\bar{Y}_m^\top$, and $d_{mk} = \sqrt{\frac{1}{p_m}\|F_{mk}\|_2^2 + 10^{-10}}$. |
| **Prior on $\beta$** | Non-collapsing empirical Bayes prior: `ebnm_point_normal` / `ebnm_point_laplace` with lower-bounded slab proportion $\pi_\beta \ge 1/K$, OR Normal with Inverse-Gamma regularized variance $\hat{\sigma}_k^2 > 0$. |
| **Supervision Scaling ($\gamma$)** | Cox working weights $w_i = \max(a_i, 10^{-6})$ (computed via log-sum-exp risk sets with Breslow ties, optionally stratified by `strata_id`) scaled by $\gamma = \frac{\alpha}{1-\alpha}\frac{\sum_m v_m p_m}{N_{\text{event}}}$. |
| **Initialization** | Non-negative hybrid start: $K-2$ factors from non-negative SVD/NMF of $Y_m$, plus 2 factors initialized from positive and negative univariate Cox score directions $\max(0, \pm \tilde{Y}_m^\top u^{(0)})$ with brief $\beta$ warm-up. |
| **CAVI Inner Sweep** | For $k = 1, \dots, K$: (1) Update $q(L_{\cdot k})$ and normalize $\frac{1}{n}\sum_i \overline{l_{ik}^2} = 1$; (2) Update $q(F_{mk})$ sequentially over $j=1\dots p_m$ with exact coupling and precision damping $\rho$; (3) Update $q(\beta_k)$ incorporating $\text{Var}_q(\tilde{Z}_{ik})$. Then update $\tau_{mj}$. Compute KL atomically from `ebnm` output. |
| **Rank Selection & Roles** | Select $(K_{\text{init}}, \alpha)$ via 5-fold training CV C-index (1-SE rule). Disable `prune_nullcheck` deletion refits; classify factors into **Survival-Active**, **Genomics-Only**, and **Dead**. |
| **Evaluation & Hygiene** | Fix predictor orientation once on training data (`reverse = TRUE`), score frozen parameters on the 568 unique external validation donors and 50 matched primary-PDAC ICGC donors, and annotate top CpGs to gene promoters. |