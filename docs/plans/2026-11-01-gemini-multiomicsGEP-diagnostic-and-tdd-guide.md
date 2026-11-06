# Architectural & Mathematical Diagnostic Guide: Resolving the `YFB` Supervision & Convergence Blockers (Plus TDD Gate Suite)

**Target Audience:** Claude Code / Development Team  
**Objective:** Resolve the core mathematical and optimization blockers preventing the `YFB`-parameterized model ($\eta = YF\beta$) from operating as a true, stable, calibrated, joint survival-supervised Bayesian matrix factorization across $M$ modalities.  
**Approach:** Do **not** build a parallel monolithic codebase or add ad-hoc heuristic patches. Instead, understand the exact mathematical root causes of the current failures in `fit_cox_on_yf.R` and `fit_multimodal_yfb.R`, implement the targeted structural solutions below, and verify them against the **7 Red-to-Green TDD Gate Tests** in Section 4 (all of which are specifically designed to **FAIL** if the model is secretly acting like a 2-step unsupervised EBMF $\to$ Cox pipeline or suffering from the October 30 instabilities).

---

## 1. Executive Summary of the Core Problem

The central goal of `multiomicsGEP` is **survival-guided discovery of molecular programs ($F$)**: the Cox proportional hazards likelihood must actively shape the feature loadings $F$ during CAVI so the model can discover:
1. **High/moderate-variance prognostic programs** (that both reconstruct $Y$ and predict survival),
2. **Genomics-only background programs** (that reconstruct outcome-neutral variance in $Y$, such as tumor purity or exocrine tissue composition, with $\beta_k \approx 0$), AND
3. **Low-variance prognostic programs** (that explain very little variance in $Y$—so unsupervised NMF/EBMF completely misses them—but strongly drive patient survival),

while using the projection linear predictor $\eta = YF\beta$ so out-of-sample validation patients are scored in closed form ($\hat{\eta}_{\text{new}} = Y_{\text{new}}\bar{F}\bar{\beta}$).

Currently, the codebase is caught between two unacceptable extremes:
- **Extreme 1 — Unsupervised 2-Step Disguised as a Joint Model ($\alpha_F = 0$ or unweighted 1:1 likelihoods):** When $\alpha_F = 0$ (single-modality benchmark) or unweighted (`w_S = 1` in multimodal), the survival precision in the $F$ update is **$10^{-5}$ to $10^{-10}$ times smaller** than the reconstruction precision (`2.03e-7` in expression, `9.00e-10` in methylation). Supervision has zero practical effect on $F$, reducing the method to a 2-step unsupervised EBMF $\to$ Cox pipeline that recovers **0%** of low-variance prognostic programs.
- **Extreme 2 — Unstable / Uncalibrated Supervision (When survival weight is amplified naively):** When survival influence on $F$ is turned up via naive multipliers or self-normalized updates, four coupled failures occur (documented in the October 30 diagnostics):
  - **Failure A:** Precision imbalance and bilinear $L \leftrightarrow F$ feedback loops.
  - **Failure B:** Loading norm growth ($F_{mjk} \to \infty$) paired with adaptive noise precision collapse ($\tau_{mj} \to 0$), which discounts reconstruction errors and ruins $LF^\top$.
  - **Failure C:** Objective decreases and 500-sweep non-convergence caused by freezing the projection normalizer $d_k(F)$ during a sweep and recomputing it after the sweep (changing the objective between sweeps).
  - **Failure D:** Excessive apparent certainty in $\beta$ (too many false-positive survival-active factors) caused by multiplying the Cox likelihood precision and ignoring correlation among factor projections in mean-field $q(\beta_k)$.

---

## 2. Issue-by-Issue Root-Cause Diagnostics & Proposed Mathematical Solutions

### Issue 1: The Severe $A_{\text{surv}} \ll A_{\text{rec}}$ Precision Imbalance in the $F$ Update

#### Mathematical Diagnosis
In the `YFB` model, $Y_m \approx L F_m^\top + E_m$ ($n \times p_m$) and $\eta_i = \sum_{k=1}^K Z_{ik} \beta_k$, where $Z_{ik} = \sum_{m=1}^M \sum_{j=1}^{p_m} \tilde{a}_{mij} f_{mjk}$ projects centered features $\tilde{y}_{mij}$ onto loading column $f_{m\cdot k}$.
Look at the scaling of the two precision terms in the coordinate update for a single feature loading $f_{mjk}$:
1. **Genomics Reconstruction Precision:**
   $$A_{F,mjk}^{\text{rec}} = \tau_{mj} \sum_{i=1}^n \overline{l_{ik}^2}$$
   When $Y_m$ has unit-variance columns and $L_{\cdot k}$ is normalized to unit root-mean-square ($\frac{1}{n}\sum_{i=1}^n \overline{l_{ik}^2} = 1$), $A_{F,mjk}^{\text{rec}} \sim O(n)$ for every feature $j = 1, \dots, p_m$.
2. **Cox Survival Precision:**
   Because $Z_{ik}$ sums over $p_m$ features, keeping the variance of the projection score $\text{Var}(Z_{\cdot k}) \sim O(1)$ requires scaling each feature's design contribution by $\tilde{a}_{mij} \sim O(1/\sqrt{p_m} \text{ to } 1/p_m)$ (or equivalently, if unscaled, $\beta_k$ shrinks to $O(1/\sqrt{p_m} \text{ to } 1/p_m)$). Consequently, for any parameterization where $\eta_i = \sum_k Z_{ik}\beta_k \sim O(1)$, the survival precision for feature $j$ is:
   $$A_{F,mjk}^{\text{surv}} = \overline{\beta_k^2} \sum_{i=1}^n w_i \tilde{a}_{mij}^2 \sim O\left( \frac{N_{\text{event}}}{p_m} \right)$$
3. **The Resulting Structural Gap:**
   $$\frac{A_{F,mjk}^{\text{surv}}}{A_{F,mjk}^{\text{rec}}} \sim O\left( \frac{N_{\text{event}}}{n \cdot p_m} \right)$$
   When $p_{\text{expr}} = 3,000$ and $p_{\text{meth}} = 10,000\text{–}300,000$, the survival precision for any individual gene or CpG is **$10^4$ to $10^7$ times smaller** than its reconstruction precision! Even though all $p_m$ features collectively determine $Z_{\cdot k}$, the coordinate-wise CAVI update for $f_{mjk}$ evaluates one feature at a time, where $A_{F,mjk}^{\text{surv}}$ is completely drowned out by $A_{F,mjk}^{\text{rec}}$.

#### Why Previous Fixes Failed
- **Why multiplying the Cox likelihood by $\gamma \sim p_m$ failed:** Multiplying the Cox term by a large weight $\gamma \gg 1$ raises $A_{F,mjk}^{\text{surv}}$ to $O(n)$, **but** it also multiplies the precision of $\beta_k$ ($A_{\beta,k}$) by $\gamma \gg 1$. That shrinks the posterior standard error $\text{SD}_q(\beta_k)$ by $1/\sqrt{\gamma}$, artificially inflating Wald statistics $|\bar{\beta}_k|/\text{SD}_q(\beta_k)$ and causing nuisance/noise factors to be flagged as survival-active (Oct 30 Diagnostic #3).
- **Why dividing `A_gen` and `B_gen` by $p$ inside `ebnm` failed in July (Phase 1a):** In Phase 1a, dividing `A_gen` and `B_gen` by $p$ while leaving the `ebnm` prior penalty unweighted meant `ebnm` saw a standard error $s = (\text{A\_gen}/p)^{-1/2} = \sqrt{p} \cdot s_{\text{orig}}$ (as if the sample size were $n / p \approx 0.1$ patients instead of $n = 273$ patients!). Seeing massive noise variance, `ebnm`'s spike-and-slab prior immediately shrank all of $L$ and $F$ to exact zero.

#### Proposed Solution: Generalized-Bayes Fractional Likelihood with Matched Prior Scaling
To balance the $n \times p_m$ genomics likelihood against the $n \times 1$ Cox partial likelihood **without inflating $\beta$'s precision AND without tricking `ebnm` into seeing a tiny sample size $n/p_m$**:
1. **Keep the Cox Partial Likelihood Weight at Exact $1.0$ for $q(\beta)$:**
   In the global Generalized-Bayes objective, weight the Cox partial likelihood by $\gamma_{\text{Cox}} = 1.0$ (so $q(\beta)$ is always updated using the true uninflated event count $N_{\text{event}}$), and weight each genomics modality $m$'s log-likelihood **and its coordinate prior KL penalty** $\text{KL}(q(f_{mjk}) \,\|\, g_{F,mk})$ together by a modality temperature $\kappa_m$:
   $$\mathcal{L}_{\text{ELBO}} = \sum_{m=1}^M \kappa_m \left( \mathbb{E}_q[\ell_{\text{rec}, m}] - \text{KL}(q(F_m) \,\|\, g_{F,m}) \right) - \text{KL}(q(L) \,\|\, g_L) + \left( \mathbb{E}_q[\ell_{\text{Cox}}^{\text{quad}}] - \text{KL}(q(\beta) \,\|\, g_\beta) \right)$$
2. **Look at What Happens When You Divide the $q(f_{mjk})$ Coordinate Objective by $\kappa_m$:**
   Because **both** $\mathbb{E}_q[\ell_{\text{rec}, m}]$ AND $\text{KL}(q(F_m) \,\|\, g_{F,m})$ are multiplied by $\kappa_m$ in $\mathcal{L}_{\text{ELBO}}$, when you optimize $\mathcal{L}_{\text{ELBO}}$ with respect to $q(f_{mjk})$ and divide through by $\kappa_m > 0$:
   - The reconstruction precision and signal inside `ebnm` stay at their **exact, unshrunk natural sample-size scale $O(n)$**:
     $$A_{F,mjk}^{\text{rec}} = \tau_{mj}\sum_{i=1}^n \overline{l_{ik}^2} \sim O(n), \quad B_{F,mjk}^{\text{rec}} = \tau_{mj}\sum_{i=1}^n R_{mij}^{(-k)}\bar{l}_{ik}$$
   - The prior penalty $\text{KL}(q(f_{mjk}) \,\|\, g_{F,mk})$ stays at weight **$1.0$** relative to $A_{F,mjk}^{\text{rec}}$ (so `ebnm` sees the true sample size $n$ and **never** collapses $F_m$ to zero!).
   - The survival precision and signal in $f_{mjk}$ are scaled by $1 / \kappa_m$:
     $$A_{F,mjk}^{\text{surv}} = \frac{1}{\kappa_m} \overline{\beta_k^2} \sum_{i=1}^n w_i \tilde{a}_{mij}^2, \quad B_{F,mjk}^{\text{surv}} = \frac{1}{\kappa_m} \text{SurvGrad}_{mjk}$$
3. **How to Set $\kappa_m$ (Calibrating the Ratio):**
   Choose $\kappa_m$ so that at unit coefficient scale ($\overline{\beta_k^2} = 1$), the median survival precision across features equals $\frac{\alpha}{1-\alpha}$ times the median reconstruction precision:
   $$\kappa_m = \frac{1 - \alpha}{\alpha} \cdot \frac{\text{median}_j\left(\sum_{i=1}^n w_i^{(0)} \tilde{a}_{mij}^2\right)}{\text{median}_j\left(\tau_{mj}^{(0)} n\right)} \quad \text{for } \alpha \in (0, 1)$$
   - **Also include Hybrid Initialization to Break the $\bar{\beta}_k = 0$ Gating Trap:** Because $A_{F,mjk}^{\text{surv}}$ and $B_{F,mjk}^{\text{surv}}$ are multiplied by $\overline{\beta_k^2}$ and $\bar{\beta}_k$, if all $K$ factors are initialized from unsupervised SVD of $Y_m$, a low-variance prognostic program orthogonal to the top $K$ principal components starts with $\bar{\beta}_k^{(0)} \approx 0$ and receives zero survival gradient. Seed $K-2$ factors from non-negative SVD/EBMF and seed 1–2 factors from the positive/negative univariate Cox score directions $\max(0, \pm \tilde{Y}_m^\top u^{(0)})$ so low-variance survival directions start with $\bar{\beta}_k^{(0)} \neq 0$.

---

### Issue 2: Loading Growth & Residual-Precision Collapse ($\tau_{mj} \to 0$ as $\|LF_m^\top\|_F \to \infty$)

#### Mathematical Diagnosis
When survival supervision pulls a subset of genes/CpGs $j$ in $F_{mk}$ away from their purely unsupervised reconstruction values, the squared reconstruction residual $Q_{mj} = \sum_{i=1}^n \mathbb{E}_q[(y_{mij} - l_i^\top f_{mj})^2]$ for those features increases.
Under an unconstrained maximum likelihood precision update $\hat{\tau}_{mj} = n / Q_{mj}$ (or even a standard Gamma prior $\hat{\tau}_{mj} = \frac{a_0 + n/2}{b_0 + Q_{mj}/2}$), increasing $Q_{mj}$ causes $\hat{\tau}_{mj}$ to drop toward zero:
$$Q_{mj} \uparrow \implies \hat{\tau}_{mj} \downarrow 0 \implies A_{F,mjk}^{\text{rec}} = \hat{\tau}_{mj}\sum_{i=1}^n \overline{l_{ik}^2} \downarrow 0$$
Once $A_{F,mjk}^{\text{rec}} \to 0$ for feature $j$, the reconstruction penalty vanishes on feature $j$, so the survival term drives $f_{mjk}$ to grow huge, which makes $Q_{mj}$ even larger and drives $\hat{\tau}_{mj} \to 0$ faster! Normalizing the Cox projection $Z$ hides this growth from the Cox likelihood while the reconstructed matrix $L F_m^\top$ explodes.

#### Proposed Solution: Anchor or Lower-Bound $\tau_{mj}$ Using Total Feature Variance
1. **Why a natural positive lower bound $\tau_{mj}^{\min} > 0$ exists mathematically:**
   In any factor model $Y_{m,\cdot j} = \mu_{mj} + L f_{mj} + e_{mj}$, the residual error variance $\sigma_{mj}^2 = \tau_{mj}^{-1}$ of feature $j$ **cannot exceed the total sample variance of feature $j$** ($\sigma_{mj}^2 \le \text{Var}(Y_{m,\cdot j})$, because a feature that is 100% noise has $\sigma_{mj}^2 = \text{Var}(Y_{m,\cdot j})$). Any state where $\hat{\sigma}_{mj}^2 \gg \text{Var}(Y_{m,\cdot j})$ (i.e., $\hat{\tau}_{mj} \ll 1/\text{Var}(Y_{m,\cdot j})$) means the model's reconstruction is worse than predicting a constant mean!
2. **Actionable Implementation (Two Options to Test in Step 1 of Your Diagnostic Sequence):**
   - **Option 2A (Frozen Unsupervised EBMF Precision — Recommended Primary Fix):** Estimate $\hat{\tau}_{mj}^{(0)}$ from an initial unsupervised EBMF fit ($\alpha = 0$) and **hold $\tau_{mj} = \hat{\tau}_{mj}^{(0)}$ fixed** during the supervised CAVI sweeps. This guarantees that the noise model of each gene/CpG is defined by its intrinsic technical/biological noise floor, and completely prevents the supervised updates from cheating the reconstruction penalty by inflating $\sigma_{mj}^2$.
   - **Option 2B (Variance-Bounded Empirical Bayes $\tau_{mj}$):** If $\tau_{mj}$ is updated during supervised CAVI, enforce the hard lower bound $\tau_{mj} \ge \tau_{mj}^{\min} = \frac{1}{\text{Var}(Y_{m,\cdot j})}$ (along with a feature intercept $\mu_{mj}$ so baseline means do not inflate $Q_{mj}$). With $\tau_{mj} \ge \tau_{mj}^{\min}$, $A_{F,mjk}^{\text{rec}}$ can **never** drop below $n / \text{Var}(Y_{m,\cdot j})$, making loading runaway impossible.

---

### Issue 3: Inconsistent Normalized Objective (Freezing $d_{mk}(F)$ During the Sweep and Recomputing Post-Sweep)

#### Mathematical Diagnosis
In the current single-modality and earlier unified designs, the Cox projection is normalized by the loading column norm:
$$\tilde{Z}_{ik} = \sum_{m=1}^M \tilde{Y}_m \frac{F_{mk}}{d_{mk}(F)}, \quad d_{mk}(F) = \sqrt{\|F_{mk}\|_2^2 + \epsilon}$$
Because $d_{mk}(F)$ is in the denominator, $\eta_i$ is a **non-linear ratio** of $F_{mk}$. When the algorithm holds $d_{mk}$ fixed during the coordinate sweep over $j = 1, \dots, p_m$ and then recomputes $d_{mk}$ at the end of the sweep (or rescales by $c_k = \text{sd}(\bar{Z}_{\cdot k})$ without transforming all priors identically), **the objective function changes between sweeps**. A coordinate step that increases the frozen-$d_{mk}$ surrogate can actually *decrease* the post-normalization objective, leading to non-monotone ELBO trajectories and 500-sweep cap failures.

#### Proposed Solution: Linear Projection Operator with Exact Bijective Scale Canonicalization (Plus Acceptance Backtracking)
You have two mathematically consistent ways to eliminate this discrepancy:

- **Approach 3A (Linear Projection with Fixed Constant Scaling $s_m$ — Simplest & Exact):**
  Once $\tau_{mj}$ is anchored/bounded (Issue 2) and $L_{\cdot k}$ is normalized to $\frac{1}{n}\sum_{i=1}^n \overline{l_{ik}^2} = 1$, the reconstruction term $\|Y_m - L F_m^\top\|_F^2$ **already bounds $\|F_{mk}\|_2$**! You no longer need an $F$-dependent denominator $d_{mk}(F)$ inside the projection $Z_{ik}$.
  Instead, define the projection using a **fixed constant scalar $s_m > 0$** (calibrated once at initialization so initial projections have unit variance, and never changed again):
  $$Z_{ik} = \sum_{m=1}^M s_m (\tilde{Y}_m F_{mk})_i$$
  - **Why this works:** Because $s_m$ is a constant, $Z_{ik}$ is a **strictly linear function of $f_{mjk}$**. There is no denominator $d_{mk}(F)$ to freeze or recompute!
  - Furthermore, the scale transformation $(l_{\cdot k}, f_{mk}, \beta_k) \mapsto (l_{\cdot k}/c_k, \, c_k f_{mk}, \, \beta_k/c_k)$ (accompanied by the exact prior pushforwards $\lambda_{L,k} \mapsto c_k \lambda_{L,k}$, $\lambda_{F,mk} \mapsto \lambda_{F,mk}/c_k$, and $\sigma_{\beta,k}^2 \mapsto \sigma_{\beta,k}^2 / c_k^2$) leaves $L F_m^\top$, $Z\beta$, and every KL divergence **100% invariant to machine precision**.
- **Approach 3B (If Retaining Explicit Normalization $u_k(F) = F_k / \|F_k\|_2$ — Exact Objective + Acceptance Backtracking):**
  If you prefer to keep $u_k(F) = F_{\cdot k} / \|F_{\cdot k}\|_2$ in the definition of the target objective, then:
  1. Evaluate the **actual post-normalization objective** $\mathcal{L}_{\text{ELBO}}^{\text{norm}}$ (using the recomputed norm $d_k(F^{\text{cand}})$, recomputed $\bar{Z}^{\text{cand}}$, and recomputed prior KL) before accepting a column update for $F_{mk}$.
  2. If the full step ($\theta = 1$) decreases $\mathcal{L}_{\text{ELBO}}^{\text{norm}}$ because $d_k(F)$ shifted, backtrack the damping/step-size $\theta \in \{1, 0.5, 0.25, 0.1, 0\}$ (recomputing the posterior moments and KL consistently at the damped precision) until $\mathcal{L}_{\text{ELBO}}^{\text{norm}}(q^{\text{cand}}) \ge \mathcal{L}_{\text{ELBO}}^{\text{norm}}(q^{\text{old}})$.

---

### Issue 4: Too Many Survival-Active Factors (Correlated Projections & Understated $\beta$ Uncertainty)

#### Mathematical Diagnosis
Even after preventing $\beta$ from collapsing to an absorbing all-zero state (by bounding $\pi_\beta \ge 1/K$ or $\hat{\sigma}_{\beta,k}^2 \ge \sigma_{\min}^2 > 0$), two statistical problems cause non-prognostic or mixed factors to be falsely classified as survival-active (`|E[beta_k]| / SD(beta_k) >= 1.96`):
1. **Temperature-inflated precision** (resolved by Issue 1's $\kappa_m$ genomics down-weighting with $\gamma_{\text{Cox}} = 1.0$).
2. **Mean-field factorization $\prod_{k=1}^K q(\beta_k)$ across correlated projections:**
   In real tumor data (and simulations), the projected factor scores $Z_{\cdot 1}, \dots, Z_{\cdot K}$ are often correlated across patients. Under a fully factorized mean-field posterior $\prod_{k=1}^K q(\beta_k)$, each coefficient $\beta_k$'s variance is computed conditionally on the *point means* of the other coefficients:
   $$\text{Var}_{q,\text{MF}}(\beta_k) = \left( \sum_{i=1}^n w_i \overline{Z_{ik}^2} + \sigma_{\beta,k}^{-2} \right)^{-1}$$
   This completely ignores the off-diagonal blocks of the Fisher information matrix $A_{\beta,kl} = \sum_{i=1}^n w_i \bar{Z}_{ik}\bar{Z}_{il}$ ($k \neq l$). It is a classic statistical property of mean-field VI that ignoring posterior correlation between collinear predictors **severely underestimates marginal posterior variances** and allows a single prognostic signal to split across multiple correlated factors.

#### Proposed Solution: Joint $K \times K$ Multivariate Variational Posterior $q(\beta) = \mathcal{N}_K(\bar{\beta}, \mathbf{\Sigma}_\beta)$
Because the number of factors $K$ is very small ($K \in [3, 15]$), there is no computational reason to factorize $q(\beta)$ into $K$ independent scalars!
- Replace $\prod_{k=1}^K q(\beta_k)$ with a **single $K$-dimensional multivariate Gaussian variational factor $q(\beta) = \mathcal{N}_K(\bar{\beta}, \mathbf{\Sigma}_\beta)$**:
  $$(\mathbf{A}_\beta)_{kl} = \sum_{i=1}^n w_i \bar{Z}_{ik}\bar{Z}_{il} + \delta_{kl} \sum_{i=1}^n w_i \text{Var}_q(Z_{ik}), \quad (\mathbf{b}_\beta)_k = \sum_{i=1}^n \bar{Z}_{ik} (u_i + w_i \hat{\eta}_i)$$
  $$\mathbf{\Sigma}_\beta = \left( \mathbf{A}_\beta + \text{diag}(\hat{\sigma}_{\beta,1}^{-2}, \dots, \hat{\sigma}_{\beta,K}^{-2}) \right)^{-1}, \quad \bar{\beta} = \mathbf{\Sigma}_\beta \mathbf{b}_\beta$$
- **Why this directly fixes false-positive survival-active counts:**
  1. Inverting the full $K \times K$ matrix $\mathbf{A}_\beta + \mathbf{\Lambda}_0$ takes $< 0.001$ milliseconds for $K \le 15$.
  2. The marginal posterior standard deviation $\text{SD}_q(\beta_k) = \sqrt{(\mathbf{\Sigma}_\beta)_{kk}}$ now incorporates the **exact variance inflation factor (VIF)** from collinearity with all other factors $l \neq k$! If a nuisance factor is partially correlated with a true prognostic factor, its marginal standard error $\sqrt{(\mathbf{\Sigma}_\beta)_{kk}}$ automatically widens, preventing false-positive survival-active calls.
  3. In the $F_{mk}$ update and the survival ELBO, use the full second-moment matrix $\mathbf{M}_\beta = \mathbb{E}_q[\beta \beta^\top] = \bar{\beta}\bar{\beta}^\top + \mathbf{\Sigma}_\beta$, so cross-factor uncertainty is propagated everywhere.

---

### Issue 5: Separating "Survival-Only Directions" from "Reconstructed Expression Programs" (Rank & Pruning Strategy)

#### Mathematical Diagnosis
1. **Why post-fit factor deletion + refitting (`prune_nullcheck`) is harmful:**
   - On real multimodal data (13,000 features), deleting a factor while holding the other $K-1$ factors frozen always decreases the reconstruction likelihood by thousands of units, so nothing is ever pruned.
   - In single-modality simulations, deleting a high-variance background factor and refitting at rank $K-1$ forces the protected low-variance prognostic factor to rotate to absorb the deleted factor's high-variance genomic signal (`|cor|` dropping from $0.896$ to $0.391$).
2. **Why survival-only directions must not be counted as reconstructed expression programs:**
   - Under strong supervision or hybrid Cox-score initialization, a factor $k$ can align with a low-variance survival gradient ($|\bar{\beta}_k| / \text{SD}_q(\beta_k) \ge 1.96$) while explaining negligible molecular variance ($\text{PVE}_k < 1\%$). Counting such a direction as a major "reconstructed expression program" misleads biological interpretation.

#### Proposed Solution: Disable Destructive Post-Fit Deletion & Use 4-Way Role Classification
1. **Do not delete factors and refit at $K-1$.** Fit at $K_{\text{init}}$ (selected via 5-fold training CV C-index with the 1-SE parsimony rule).
2. **Classify all $K_{\text{init}}$ fitted factors along BOTH axes independently:**
   - **`"prognostic_program"` (Shared Reconstruction + Survival-Active):** $\text{PVE}_k \ge 1\%$ AND $|\bar{\beta}_k|/\text{SD}_q(\beta_k) \ge 1.96$. *(True molecular programs that both reconstruct expression/methylation variance and predict survival.)*
   - **`"survival_only_direction"` (Survival-Active, Low Reconstruction Variance):** $\text{PVE}_k < 1\%$ AND $|\bar{\beta}_k|/\text{SD}_q(\beta_k) \ge 1.96$. *(Reported separately as a pure prognostic projection direction, not counted as a reconstructed molecular program.)*
   - **`"genomics_only_program"` (Reconstruction-Active, Survival-Silent):** $\text{PVE}_k \ge 1\%$ AND $|\bar{\beta}_k|/\text{SD}_q(\beta_k) < 1.96$. *(Outcome-neutral biological variance.)*
   - **`"inactive"` (Dead / Shrunk):** $\text{PVE}_k < 1\%$ AND $|\bar{\beta}_k|/\text{SD}_q(\beta_k) < 1.96$.

---

### Issue 6: Multi-Modal Scaling ($p_1 = 3,000$ Genes vs. $p_2 = 10,000\text{–}300,000$ CpGs) & Data Hygiene

#### Mathematical Diagnosis & Solutions
1. **Equalizing Modality Contributions in Shared Activities $L$ and Projections $Z$:**
   - Without dimension normalization, $10,000\text{–}300,000$ CpGs contribute $3.3\times$ to $100\times$ more terms than $3,000$ genes to $A_{L,k} = \sum_{m,j} \tau_{mj}\overline{f_{mjk}^2}$ and to $Z_{ik}$, drowning out gene expression and causing numerical overflow in $\hat{\eta}_i$ at 25,000+ CpGs.
   - **Fix:** Weight each modality $m$'s reconstruction log-likelihood (and $F_m$ KL penalty) by $v_m = \bar{p}/p_m$, normalize each modality's projection contribution by $p_m$ (or initial projection SD), and use log-sum-exp Breslow risk-set arithmetic with a curvature floor $w_i = \max(a_i, 10^{-6})$.
2. **Reconciling Non-Negative NMF ($L \ge 0, F_m \ge 0$) with Centered Cox Projections:**
   - **Fix:** Keep reconstruction inputs non-negative ($Y_{m,ij} \ge 0$, scaled by within-cohort feature SD with a baseline offset $\mu_{mj}$ so `ebnm_point_exponential` can be used on both $L$ and $F_m$ for parts-based interpretability), and center features **only inside the Cox projection design matrix** ($\tilde{Y}_{m,ij} = Y_{m,ij} - \bar{Y}_{m,\cdot j}$) so every factor projection $Z_{\cdot k}$ has mean zero within each cohort.
3. **Data & Evaluation Hygiene Fixes:**
   - Ensure $\log_2(x + 1)$ is applied **only once** to unlogged RNA-seq cohorts (fixing the historical double-log path on TCGA, Dijk, and PACA-AU RNA-seq).
   - Fix training risk orientation **once** on training data using `survival::concordance(..., reverse = TRUE)` (never `max(C, 1-C)` on validation data).
   - Deduplicate the 568 unique external validation donors (accounting for the 47 shared PACA-AU donors across microarray and RNA-seq plus 1 duplicate RNA-seq record) when computing pooled survival curves and cohort-stratified bootstrap CIs.

---

## 3. Recommended Controlled Experimental Sequence for Claude Code

Follow the exact four-step diagnostic sequence agreed upon in the October 30 review before running any new full grid:

1. **Step 1 — Frozen-vs-Adaptive Precision Test (`Isolating Issue 2`):**
   - Take the saved failing simulation seeds where $\|LF^\top\|_F$ exploded and $\tau_j \to 0$.
   - Freeze $\tau_j$ at its initial unsupervised EBMF estimate $\hat{\tau}_j^{(0)}$ (and compare against variance-bounded $\tau_j \ge 1/\text{Var}(Y_{\cdot j})$).
   - **Verification Gate:** Confirm that freezing/bounding $\tau_j$ completely eliminates reconstruction explosion and loading runaway on those exact seeds.
2. **Step 2 — Single Consistent Objective & Acceptance Check (`Isolating Issue 3`):**
   - Implement the linear constant-scaled projection (or post-normalization objective check with backtracking) and verify that `compute_elbo()` is **monotonically non-decreasing ($\Delta \text{ELBO} \ge -10^{-8}$) on 100% of inner steps** across the previously non-converging seeds.
3. **Step 3 — Generalized-Bayes Temperature ($\kappa_m$) + Joint $K \times K$ Multivariate $q(\beta)$ (`Isolating Issues 1 & 4`):**
   - Implement down-weighted genomics ($\kappa_m$) with uninflated Cox weight ($\gamma_{\text{Cox}} = 1.0$) and replace scalar `update_beta` with the $K \times K$ multivariate Gaussian update $\mathbf{\Sigma}_\beta = (\mathbf{A}_\beta + \mathbf{\Lambda}_0)^{-1}$.
   - **Verification Gate:** In known-loading and planted low-variance simulations, confirm that $\kappa_m$ allows $F$ to recover the planted low-variance prognostic program ($\ge 80\%$ recovery) while $\gamma_{\text{Cox}} = 1.0$ + multivariate $\mathbf{\Sigma}_\beta$ keeps false-positive survival-active null factors within the baseline $+ 0.5$ gate.
4. **Step 4 — Full Simulation & Real-Data Multi-Omics Gate:**
   - Only after Steps 1–3 pass on the controlled fixtures, run the full 20-seed low/high-noise simulation panels and the TCGA $\to$ ICGC multimodal real-data evaluation, reporting `"prognostic_program"`, `"survival_only_direction"`, `"genomics_only_program"`, and `"inactive"` counts separately.

---

## 4. Red-to-Green TDD Gate Suite: 7 Tests That MUST FAIL Until Supervision of $F$ Is Truly Solved

Have Claude Code write and run these **7 specific unit/integration tests** (`tests/test_supervision_of_F_gates.R`). Each test targets a specific failure mode of the current codebase so that **a disguised 2-step unsupervised EBMF $\to$ Cox model or an unstable supervised model will deterministically FAIL**, and only a properly balanced, stable, calibrated supervised-$F$ model will PASS.

### Test 1: The Precision Ratio Gate (`A_surv / A_rec` Order-of-Magnitude Check)
- **What It Tests:** Directly inspects the internal coordinate precision quantities $A_{F,jk}^{\text{rec}}$ and $A_{F,jk}^{\text{surv}}$ passed into the `update_F` EBNM solver on an active prognostic factor ($|\bar{\beta}_k| > 0$) when $\alpha = 0.5$ on realistic dimensions ($n = 150, p = 3,000$ genes or $p = 10,000$ CpGs).
- **Why Current Code FAILS:** In the current benchmark ($\alpha_F = 0$), $A_{\text{surv}} / A_{\text{rec}} = 0$. In the current multimodal code (`w_S = 1`), the audit in Appendix C.1 showed $\text{median}(A_{\text{surv}} / A_{\text{rec}}) = 2.03 \times 10^{-7}$ for expression and $9.00 \times 10^{-10}$ for methylation.
- **Pass Criterion (Assertion):**
  - When $\alpha = 0.5$ and factor $k$ has a unit-scale survival effect, the median effective survival-to-reconstruction precision ratio across active features in `update_F` **must lie in $[0.10, 10.0]$** (i.e., $O(1)$, never $< 10^{-3}$!).
  - Furthermore, the effective sample-size precision seen by `ebnm` relative to its prior KL penalty must remain $O(n)$ ($\text{median}(A_{F,jk}^{\text{rec}}) \ge 0.5 n$), ensuring $F$ does not collapse to zero.

### Test 2: The Single-Step Survival Gradient Perturbation Test (`Does Changing y/status Move F?`)
- **What It Tests:** Holds $Y$, initial $L^{(0)}$, initial $F^{(0)}$, and non-zero $\bar{\beta}^{(0)}$ identical, and runs **a single $F$ update** (`update_F`) under two different survival outcome vectors:
  1. True survival outcomes $(t, \delta)$ strongly aligned with gene set $\mathcal{G}_{\text{prog}}$.
  2. Permuted/inverted survival outcomes $(t_{\text{flip}}, \delta_{\text{flip}})$ aligned with a disjoint gene set $\mathcal{G}_{\text{alt}}$.
- **Why Current Code FAILS:** Under $\alpha_F = 0$, $F^{(1)}$ is bit-for-bit identical (`max(abs(F_1 - F_2)) == 0`). Under unweighted `w_S = 1` at $p = 3,000$, `max(abs(F_1 - F_2)) / max(F_1) < 1e-5` because the survival term is $10^{-7}$ of the reconstruction term.
- **Pass Criterion (Assertion):**
  - At $\alpha = 0.5$, a single sweep of `update_F` must produce a material, directional shift in $F$ toward the survival-aligned genes:
    $$\frac{\|F^{(1)}(t, \delta) - F^{(1)}(t_{\text{flip}}, \delta_{\text{flip}})\|_F}{\|F^{(0)}\|_F} > 0.05$$
    and the mean loading on $\mathcal{G}_{\text{prog}}$ must be strictly larger under $(t, \delta)$ than under $(t_{\text{flip}}, \delta_{\text{flip}})$.

### Test 3: The "Needle-in-a-Haystack" Low-Variance Prognostic Program Recovery Test ($K_{\text{fit}} < K_{\text{true}}$ Bottleneck)
- **What It Tests:** This is the ultimate test separating **joint supervised NMF** from **2-step unsupervised NMF $\to$ Cox**.
  - Simulate a dataset ($n = 200, p = 1,000$) with **$K_{\text{true}} = 5$ true non-negative factors**:
    - Factors 1, 2, 3, 4 are **high-variance outcome-neutral background programs** (each explaining $15\%\text{–}20\%$ of expression variance, with true survival coefficient $\beta_1 = \beta_2 = \beta_3 = \beta_4 = 0$).
    - Factor 5 is a **low-variance prognostic program** (explaining only **$0.8\%$ of expression variance**, with true survival coefficient $\beta_5 = +2.5$).
  - Now fit the model with rank **$K_{\text{fit}} = 3$ (or $K_{\text{fit}} = 4 < 5$)**!
- **Why 2-Step EBMF $\to$ Cox and Current Code Deterministically FAIL:**
  - When $K_{\text{fit}} = 4 < 5$, any unsupervised NMF/EBMF (or $\alpha_F = 0$ / unweighted `w_S = 1`) **must** spend all 4 of its factors capturing the 4 highest-variance programs (Factors 1–4, which explain $75\%$ of the variance in $Y$). Factor 5 ($0.8\%$ variance) is completely ignored as residual noise! Therefore, 2-step EBMF $\to$ Cox achieves **loading correlation $< 0.25$** with Factor 5 and **held-out test C-index $\approx 0.50\text{–}0.55$ (chance)**.
- **Pass Criterion (Assertion):**
  - The supervised `YFB` model ($\alpha = 0.5$) at $K_{\text{fit}} = 4$ must:
    1. Recover the planted low-variance Factor 5 with loading correlation $|\text{cor}(\hat{F}_{\cdot k^*}, F_{\cdot 5}^{\text{true}})| \ge 0.75$.
    2. Achieve held-out test set C-index $C_{\text{test}}^{\text{sup}} \ge 0.75$ and beat 2-step EBMF $\to$ Cox on the exact same test set by **$\Delta C_{\text{test}} = C_{\text{test}}^{\text{sup}} - C_{\text{test}}^{\text{2-step}} \ge +0.15$**!

### Test 4: The Tumor-Purity / Confounded-Marker Decoupling Test (`DeSurv` Primary Scenario)
- **What It Tests:** Simulates a gene set of 50 prognostic marker genes whose expression is **mixed/confounded** with a high-variance outcome-neutral background factor (e.g., 50% marker loading + 50% shared background loading, as in `DeSurv` SI Fig. S2 "Mixed Scenario").
- **Why 2-Step Unsupervised NMF FAILS:** Unsupervised NMF merges the prognostic markers into the dominant outcome-neutral background factor, diluting the projection score $Z = YF$ with 500 non-prognostic background genes and depressing out-of-sample C-index and marker precision.
- **Pass Criterion (Assertion):**
  - Supervised `YFB` ($\alpha > 0$) must isolate the prognostic marker genes into a survival-active factor with **top-50 marker precision $\ge 0.60$** (vs. $< 0.20$ for $\alpha = 0$) and achieve a higher out-of-sample C-index than $\alpha = 0$.

### Test 5: The Reconstruction & Precision Stability Gate (`Anti-Oct-30 Failure #1`)
- **What It Tests:** Runs supervised `YFB` across $\alpha \in \{0.25, 0.50, 0.75\}$ for 200 sweeps on both the low-variance prognostic simulation and real TCGA data, monitoring $\min_j \hat{\tau}_{mj}$ and the relative reconstruction error $\frac{\|Y_m - (\mathbf{1}\hat{\mu}_m^\top + \bar{L}\bar{F}_m^\top)\|_F}{\|Y_m\|_F}$.
- **Why Oct 30 Code FAILS:** Under the Oct 30 run, survival supervision caused loading norms $\|F_{mk}\|_2$ to grow huge while $\hat{\tau}_{mj} \to 0$, causing the reconstructed product $\|LF_m^\top\|_F / \|Y_m\|_F$ to explode to $> 10\times$.
- **Pass Criterion (Assertion):**
  - Across all sweeps $r = 1, \dots, R$:
    1. $\min_{m,j} \hat{\tau}_{mj} \ge \tau_{\min} > 0$ (no precision collapse).
    2. $\frac{\|\mathbf{1}\hat{\mu}_m^\top + \bar{L}\bar{F}_m^\top\|_F}{\|Y_m\|_F} \in [0.1, 2.0]$ (reconstructed matrix norm remains strictly bounded on the scale of $Y_m$).
    3. Reconstructed RMSE at $\alpha = 0.5$ is within $15\%$ of the unsupervised $\alpha = 0$ RMSE.

### Test 6: The Strict Inner-Step Monotone ELBO Gate (`Anti-Oct-30 Failure #2`)
- **What It Tests:** Records the exact unified objective $\mathcal{L}_{\text{ELBO}}$ before and after every individual factor-block update ($q(L_{\cdot k})$, scale canonicalization, $q(F_{mk})$, and $q(\beta)$) across 50 sweeps with fixed Cox working weights $(w, h)$.
- **Why Oct 30 Code FAILS:** Freezing $d_{mk}(F)$ during the $F$ sweep and recomputing it after the sweep caused $\mathcal{L}_{\text{ELBO}}$ to drop on $50/350$ $F$-updates and fail to converge on $153/360$ fits.
- **Pass Criterion (Assertion):**
  - Scale canonicalization changes $\mathcal{L}_{\text{ELBO}}$ by $< 10^{-12}$ (exact invariance).
  - Every single inner update step for $L$, $F_m$, and $\beta$ satisfies $\mathcal{L}_{\text{ELBO}}^{\text{after}} - \mathcal{L}_{\text{ELBO}}^{\text{before}} \ge -10^{-8}$ (100% monotone ascent), and the outer loop converges in $< 150$ sweeps.

### Test 7: The Collinear-Factor False-Positive & Null Calibration Gate (`Anti-Oct-30 Failure #3`)
- **What It Tests:**
  - **Part A (Null Scenario $\beta_{\text{true}} = 0$):** Fit supervised `YFB` ($\alpha = 0.5, K = 7$) on 20 simulated datasets where $Y$ has 4 strong biological factors but survival times are completely independent of $Y$ ($\beta_{\text{true}} = 0$).
  - **Part B (Collinear Nuisance Scenario):** Fit supervised `YFB` ($\alpha = 0.5, K = 5$) where Factor 1 is truly prognostic ($\beta_1 = 2.0$), Factor 2 is a non-prognostic factor ($\beta_2 = 0$) whose patient scores $L_{\cdot 2}$ correlate at $r = 0.65$ with $L_{\cdot 1}$, and Factors 3–5 are orthogonal noise/background.
- **Why Oct 30 Code FAILS:** Temperature-inflated Cox weights ($\gamma \gg 1$) and scalar mean-field $q(\beta_k)$ underestimate $\text{SD}_q(\beta_k)$, flagging Factor 2 (and null factors in Part A) as `"survival_active"` (`|E[beta]| / SD(beta) >= 1.96`).
- **Pass Criterion (Assertion):**
  - In Part A (Null): Mean false-positive survival-active factors per fit is $\le 0.35$, and mean held-out test C-index (using frozen training orientation, never `max(C, 1-C)`) is within $[0.46, 0.54]$.
  - In Part B (Collinear): Multivariate $q(\beta) = \mathcal{N}_K(\bar{\beta}, \mathbf{\Sigma}_\beta)$ with $\gamma_{\text{Cox}} = 1.0$ correctly flags Factor 1 (`|z_1| >= 1.96`) while keeping the collinear nuisance Factor 2 below threshold (`|z_2| < 1.96`) in $\ge 80\%$ of replicates.


---

## 5. Four Final Mathematical Clarifications (Preventing Subtle Implementation Traps)

When implementing the equations in Section 2, verify that the following four mathematical identities are respected:

1. **Why $q(L)$ Must NOT Be Divided by $\kappa_m$ (Preserving the $L$ Sample Size Scale):**
   - In Section 2 Issue 1, the global Generalized-Bayes objective is written as:
     $$\mathcal{L}_{\text{ELBO}} = \sum_{m=1}^M \kappa_m \left( \mathbb{E}_q[\ell_{\text{rec}, m}] - \text{KL}(q(F_m) \,\|\, g_{F,m}) \right) - \text{KL}(q(L) \,\|\, g_L) + \left( \mathbb{E}_q[\ell_{\text{Cox}}^{\text{quad}}] - \text{KL}(q(\beta) \,\|\, g_\beta) \right)$$
   - Look closely at $q(L)$: if $\text{KL}(q(L) \,\|\, g_L)$ is left at weight $1.0$ while $\mathbb{E}_q[\ell_{\text{rec}, m}]$ is multiplied by a small $\kappa_m \sim O(N_{\text{event}} / (n p_m))$, then the precision $A_{L,k} = \sum_m \kappa_m \sum_j \tau_{mj}\overline{f_{mjk}^2}$ would shrink from $O(p)$ down to $O(N_{\text{event}}/n) \approx 0.5$, causing `ebnm` to over-shrink $L$!
   - **Exact Rule:** In the $q(L_{\cdot k})$ update, $L$ only participates in the genomics reconstruction across modalities (never in the Cox predictor $\eta = Z\beta$). Therefore, weight $\text{KL}(q(L) \,\|\, g_L)$ by $\bar{\kappa} = \frac{1}{M}\sum_{m=1}^M \kappa_m p_m / \bar{p}$ (or equivalently, update $q(L_{\cdot k})$ using the standard dimension-balanced modality weights $v_m = \bar{p}/p_m$: $A_{L,k} = \sum_{m=1}^M v_m \sum_{j=1}^{p_m} \tau_{mj}\overline{f_{mjk}^2}$). That way, $q(L)$ always sees the full $O(\bar{p})$ feature precision, $q(F_m)$ sees the full $O(n)$ sample precision, and $q(\beta)$ sees the exact $O(N_{\text{event}})$ event precision!

2. **Exact Cross-Factor Coupling Term in `SurvGrad` When Using Multivariate $q(\beta) = \mathcal{N}_K(\bar{\beta}, \mathbf{\Sigma}_\beta)$:**
   - In Section 2 Issue 1, $\text{SurvGrad}_{mjk}$ is written in shorthand. When $q(\beta) = \mathcal{N}_K(\bar{\beta}, \mathbf{\Sigma}_\beta)$ is a full $K \times K$ multivariate normal (Issue 4) with second-moment matrix $\mathbf{M}_\beta = \mathbb{E}_q[\beta \beta^\top] = \bar{\beta}\bar{\beta}^\top + \mathbf{\Sigma}_\beta$, the off-diagonal expectations $\mathbb{E}_q[\beta_k \beta_l] = (\mathbf{M}_\beta)_{kl}$ are non-zero for $l \neq k$!
   - Therefore, the exact survival precision and linear signal for feature $j$ of modality $m$ in factor $k$ are:
     $$A_{F,mjk}^{\text{surv}} = \frac{1}{\kappa_m} (\mathbf{M}_\beta)_{kk} \sum_{i=1}^n w_i \tilde{a}_{mij}^2$$
     $$B_{F,mjk}^{\text{surv}} = \frac{1}{\kappa_m} \sum_{i=1}^n \tilde{a}_{mij} \left[ \bar{\beta}_k h_i - w_i \sum_{l \neq k} (\mathbf{M}_\beta)_{kl} \bar{Z}_{il} - w_i (\mathbf{M}_\beta)_{kk} \left(\bar{Z}_{ik} - \tilde{a}_{mij}\bar{f}_{mjk}\right) \right]$$
   - Do **not** use $\bar{\beta}_k \bar{\beta}_l$ in place of $(\mathbf{M}_\beta)_{kl}$ in the cross-factor subtraction term; using the full $(\mathbf{M}_\beta)_{kl}$ ensures the coordinate update for $f_{mjk}$ is the exact derivative of $\mathbb{E}_q[\ell_{\text{Cox}}^{\text{quad}}]$!

3. **Guarding the Calibrated Ratio $\kappa_m$ Against Zero Weights:**
   - In Section 2 Issue 1, $\kappa_m$ uses $\text{median}_j(\sum_{i=1}^n w_i^{(0)} \tilde{a}_{mij}^2)$ in the numerator. Because some features in methylation or expression can have near-zero variance or sparse design entries, always compute the reference precision ratio using the **mean across features** (or the median over non-constant features with $\text{Var}(\tilde{y}_{m,\cdot j}) > 10^{-6}$):
     $$\kappa_m = \frac{1 - \alpha}{\alpha} \cdot \frac{\frac{1}{p_m}\sum_{j=1}^{p_m} \left(\sum_{i=1}^n w_i^{(0)} \tilde{a}_{mij}^2\right)}{\frac{1}{p_m}\sum_{j=1}^{p_m} \left(\tau_{mj}^{(0)} \sum_{i=1}^n \overline{l_{ik}^{2,(0)}}\right)}$$
     This is strictly positive and well-defined for any matrix $Y_m$.

4. **Test-Time Projection Invariance to $\kappa_m$:**
   - Confirm that $\kappa_m$ is strictly an inference-time likelihood balancing weight during training; at test time, out-of-sample risk prediction uses only the frozen projection operator $\hat{\eta}_{\text{new}} = \sum_{m=1}^M s_m \tilde{Y}_{m,\text{new}} \bar{F}_m \bar{\beta}$ (with modality availability rescaling if a validation cohort only measures a subset of features/modalities).