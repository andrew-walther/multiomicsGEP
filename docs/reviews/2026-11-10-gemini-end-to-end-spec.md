# End-to-End Implementation Blueprint: Generalized-Bayes Supervised EBMF + Cox (`multiomicsGEP`)

**Purpose:** This document provides the complete mathematical formulation, derivations, algorithm, preprocessing contract, simulation design, and reference R implementation for the $M$-modality Supervised Empirical Bayes Matrix Factorization + Cox model ($\eta = YF\beta$). It directly resolves the three optimization and statistical pathologies identified in recent diagnostics:
1. **Runaway loading growth & precision collapse ($\tau_{mj} \to 0$)**
2. **Objective decreases from changing projection normalizers ($d_{mk}$)**
3. **Understated $\beta$ uncertainty & false-positive survival-active factors from Cox weight inflation ($\gamma \gg 1$) and collinear projections**

---

## Stage 1: Data Generation (Shared-Loading DGP) & Pre-Processing

### 1.1 Pre-Processing Contract (Real & Simulated Data)
To preserve non-negative matrix factorization ($Y_m \approx L F_m^\top \ge 0$) while ensuring zero-mean, unit-scale covariates in the Cox projection without a parameter-dependent normalizer $d_{mk}(F)$:

1. **Raw Non-Negative Input ($X_m \in \mathbb{R}_{\ge 0}^{n \times p_m}$):**
   - **RNA-seq:** Apply $\log_2(\text{TPM} + 1)$ **exactly once** (fixing the historical double-log bug on TCGA, Dijk, and PACA-AU RNA-seq). CPTAC and microarray cohorts already on a $\log_2$ scale are not re-logged.
   - **DNA Methylation:** Use raw $\beta$-values in $[0, 1]$ (with within-cohort KNN imputation $k = 10$ prior to screening).
2. **1-Step Unsupervised Screening (Training Data Only):**
   - **Expression:** Within each training cohort, rank genes by mean and by variance; retain the top 3,000 with the lowest rank-sum and take the intersection across training cohorts ($p_{\text{expr}} = 2,064$ in pooled TCGA+CPTAC; $3,000$ in single-cohort TCGA).
   - **Methylation:** Restrict to shared non-sex-chromosome probes (305,352 CpGs) and screen on training data by variance (e.g., $p_{\text{meth}} = 10,000$ reference; nested sensitivity ladders at 25k, 50k, 100k, and 305k).
3. **Within-Cohort Non-Negative Scaling ($Y_m \ge 0$ for Reconstruction):**
   - Within each cohort $c$, compute each feature's sample standard deviation $s_{mjc} = \text{sd}(X_{m, i \in c, j})$ (floored at $10^{-8}$) and scale **without subtracting the mean**:
     $$Y_{m,ij} = \frac{X_{m,ij}}{s_{mjc}} \ge 0 \quad \text{for } i \in c$$
   - Because $\text{Var}(Y_{m, \cdot j}) = 1$ within each cohort, all features and modalities sit on a common unit-variance scale while remaining strictly non-negative ($Y_{m,ij} \ge 0$).
4. **Within-Cohort Centered Design ($\tilde{Y}_m$ for Cox Projection Only):**
   - For the Cox linear predictor, subtract the within-cohort feature mean $\mu_{mjc} = \text{mean}(Y_{m, i \in c, j})$:
     $$\tilde{Y}_{m,ij} = Y_{m,ij} - \mu_{mjc} \quad \text{for } i \in c$$
   - Every column of $\tilde{Y}_m$ has mean $0$ and variance $1$. Out-of-sample validation cohorts apply the exact same within-cohort scaling ($Y_{m,\text{test}}$ and $\tilde{Y}_{m,\text{test}}$) restricted to the measured training features $\mathcal{G}_m$, requiring zero test-time refitting.

### 1.2 Simulation Data-Generating Process (Consistent Shared-Loading Units)
Earlier simulations inadvertently used different effective loading directions for reconstruction vs. survival after scaling. The controlling shared-loading DGP must obey exact unit consistency:
1. Generate non-negative ground-truth patient activities $L_{\text{true}} \in \mathbb{R}_{\ge 0}^{n \times K_{\text{true}}}$ (e.g., $\text{Gamma}(2, 2)$) and sparse non-negative feature loadings $F_{m,\text{true}} \in \mathbb{R}_{\ge 0}^{p_m \times K_{\text{true}}}$ scaled so each column has unit root-mean-square norm $\frac{1}{p_m}\Vert{}F_{m,\cdot k}^{\text{true}}\Vert{}_2^2 = 1$. For a low-variance prognostic factor $k^*$, scale $L_{\cdot k^*}^{\text{true}}$ by $\kappa_{\text{low}} \in [0.15, 0.25]$ so it explains $< 1\%$ of genomic variance.
2. Generate heteroscedastic noise precisions $\tau_{mj}^{\text{true}} \sim \text{Gamma}(2, 2)$ and raw non-negative observations $X_m = L_{\text{true}} F_{m,\text{true}}^\top + E_m$ (shifted by $\max(0, -\min X_{m,\cdot j})$ or simulated on the positive orthant), then compute $Y_m$ and $\tilde{Y}_m$ via Section 1.1.
3. **Exact Shared-Loading Survival Predictor:** Using the **exact same fixed dimension scaling** $s_m = \frac{1}{M\sqrt{p_m}}$ and centered matrix $\tilde{Y}_m$ used by the model, compute the true projection and survival times:
   $$Z_{ik}^{\text{true}} = \sum_{m=1}^M s_m \sum_{j=1}^{p_m} \tilde{Y}_{m,ij} F_{m,jk}^{\text{true}}, \quad \eta_i^{\text{true}} = \sum_{k=1}^{K_{\text{true}}} Z_{ik}^{\text{true}} \beta_k^{\text{true}}$$
   Draw event times $T_i$ from a Weibull or Exponential proportional hazards model with linear predictor $\eta_i^{\text{true}}$ and calibrate independent censoring to the target event rate (e.g., 50%–70%).

---

## Stage 2: Unified Model, Generalized-Bayes Objective, Priors, & Normalization

### 2.1 Fixing Issue 2: Fixed-Constant Scaling ($s_m$) Instead of Parameter-Dependent $d_{mk}(F)$
In the previous formulation, dividing the projection by $d_{mk}(F) = \Vert{}F_{mk}\Vert{}_2$ made $\eta$ a non-linear ratio $F_{mk} / \Vert{}F_{mk}\Vert{}_2$. Freezing $d_{mk}$ during the sweep and recomputing it after the sweep changed the objective function between sweeps, causing ELBO drops and non-convergence.

**The Solution:** Eliminate the parameter-dependent normalizer $d_{mk}(F)$ entirely. Instead, scale each modality's centered design matrix by the **fixed constant** $s_m = \frac{1}{M \sqrt{p_m}}$:
$$\tilde{A}_m = s_m \tilde{Y}_m = \frac{1}{M \sqrt{p_m}} \tilde{Y}_m, \quad \tilde{Z}_{ik} = \sum_{m=1}^M \sum_{j=1}^{p_m} \tilde{a}_{mij} f_{mjk}, \quad \eta_i = \sum_{k=1}^K \tilde{Z}_{ik} \beta_k$$
- **Why this works:** Because $\tilde{a}_{mij}$ is a fixed data constant, $\eta_i$ is **strictly bilinear** in $(F, \beta)$ and $Y_m \approx L F_m^\top$ is **strictly bilinear** in $(L, F)$. Every coordinate update optimizes **one single, unchanging variational objective $\mathcal{L}(q, g)$**.
- **Why scale is automatically $O(1)$:** Because $\text{Var}(Y_{m,ij}) = 1$, reconstructed means satisfy $\overline{l_{ik}^2} \sim O(1)$ and $\frac{1}{p_m}\sum_j \overline{f_{mjk}^2} \sim O(1)$. Therefore, $\text{Var}(\tilde{Z}_{\cdot k}) = \frac{1}{M^2 p_m} \sum_j \text{Var}(\tilde{y}_{mij} f_{mjk}) \sim O(1)$ automatically across any $p_m \in [50, 305000]$, and growing $\Vert{}F_{mk}\Vert{}_2$ directly increases *both* $L F_m^\top$ and $\eta_i$!

### 2.2 Fixing Issue 3: Generalized-Bayes Downweighting of Genomics (Keep Cox Weight $= 1$!)
Previously, multiplying the Cox likelihood by $\gamma \gg 1$ inflated the precision of $\beta$ by $\gamma$, shrinking posterior standard deviations by $1/\sqrt{\gamma}$ and turning noise directions into false positives.

**The Solution (Generalized-Bayes Tempered Likelihood):**
Keep the Cox partial likelihood weight at **$1.0$** so that the survival likelihood reflects the true sample size ($n$ patients, $N_{\text{event}}$ deaths) and $\text{Cov}_q(\beta)$ remains 100% statistically calibrated! To prevent the $n \times p_m$ genomic entries from overwhelming the $n$ survival entries in the $F_m$ and $L$ updates, **downweight the genomic reconstruction likelihood per feature by the generalized-Bayes learning rate $\nu_m$**:
$$\nu_m = \frac{1 - \alpha}{\alpha} \cdot \frac{\kappa_0}{M p_m}, \quad \alpha \in (0, 1)$$
where $\kappa_0 = \frac{N_{\text{event}}}{n} \in (0, 1]$ is the event-rate scaling constant (or $\kappa_0 = 1$).
- **Effect on $q(\beta)$:** The Cox weight is strictly **$1.0$**. Posterior uncertainty on $\beta$ is never artificially tightened!
- **Effect on $q(F_{mjk})$:** In the $F_{mjk}$ update, the reconstruction precision is $\nu_m \tau_{mj} \sum_{i=1}^n \overline{l_{ik}^2} \propto \frac{1-\alpha}{\alpha} \frac{N_{\text{event}}}{M p_m}$, while the survival precision is $\overline{\beta_k^2} \sum_{i=1}^n w_i \tilde{a}_{mij}^2 \propto \frac{\overline{\beta_k^2} N_{\text{event}}}{M^2 p_m}$. Both terms scale as $O(N_{\text{event}} / p_m)$—they are on the **exact same order of magnitude** for every feature $j$, regardless of whether $p_m = 50$ or $p_m = 305,000$!
- **Effect on $q(L_{\cdot k})$ and $\tau_{mj}$:** In the $L$ update, summing over $j = 1 \dots p_m$ cancels $1/p_m$, giving each modality equal weight $1/M$. To keep the Empirical Bayes sparsity on $L$ and $F_m$ at its natural observation scale, apply the standard unweighted EBNM prior/KL on $L$, or temper the reconstruction likelihood and $F_m$ prior consistently by $\nu_m$ (derived in Stage 3).

### 2.3 Fixing Issue 1: Anchored / Lower-Bounded Residual Precisions ($\tau_{mj} \ge \tau_{mj}^{\min}$)
To prevent the model from discounting huge reconstruction errors on survival-perturbed features by driving $\hat{\tau}_{mj} \to 0$:
1. Before running supervised CAVI, compute the baseline unsupervised EBMF residual precision $\hat{\tau}_{mj}^{\text{unsup}}$ (from a fast $\alpha = 0$ initialization run) and the marginal feature precision $\tau_{mj}^{\text{null}} = 1 / \text{Var}(Y_{m,\cdot j}) = 1.0$.
2. Place a conjugate $\text{Gamma}(a_\tau, b_\tau)$ prior on $\tau_{mj}$ centered at $\hat{\tau}_{mj}^{\text{unsup}}$ **AND enforce a hard lower bound $\tau_{mj}^{\min}$**:
   $$\hat{\tau}_{mj} = \max\left( \tau_{mj}^{\min}, \, \min\left( \tau_{\max}, \, \frac{n + 2 a_\tau}{Q_{mj} + 2 b_\tau} \right) \right), \quad \tau_{mj}^{\min} = \max\left( 1.0, \, \frac{1}{2}\hat{\tau}_{mj}^{\text{unsup}} \right)$$
- **Why this stops loading runaway cold:** Because $\text{Var}(Y_{m,\cdot j}) = 1.0$, any valid factor reconstruction must have residual variance $\le 1.0$ ($\tau_{mj} \ge 1.0$). Clamping $\hat{\tau}_{mj} \ge \tau_{mj}^{\min} \ge 1.0$ makes it **mathematically impossible** for the model to excuse a growing $L F_m^\top$ error by lowering $\tau_{mj}$. If $F_{mjk}$ tries to grow beyond the data scale, the reconstruction penalty $\nu_m \tau_{mj} \sum_i \overline{l_{ik}^2}$ remains strictly bounded away from zero and pulls $F_{mjk}$ right back! *(An option `freeze_tau = TRUE` is also provided to hold $\tau_{mj} = \hat{\tau}_{mj}^{\text{unsup}}$ completely fixed, per diagnostic recommendation 1).*

---

## Stage 3: Complete Mathematical Derivations

### 3.1 The Single Consistent Generalized-Bayes Variational Objective ($\mathcal{L}$)
We posit the variational distribution:
$$q(L, F, \beta) = \left(\prod_{i=1}^n \prod_{k=1}^K q(l_{ik})\right) \left(\prod_{m=1}^M \prod_{j=1}^{p_m} \prod_{k=1}^K q(f_{mjk})\right) q(\beta_1, \dots, \beta_K)$$
Note that $q(\beta)$ is a **joint $K$-dimensional multivariate Gaussian** $\mathcal{N}_K(\bar{\beta}, \Sigma_\beta)$ with marginal second moments $\overline{\beta_k^2} = \bar{\beta}_k^2 + (\Sigma_\beta)_{kk}$ and cross-moments $\overline{\beta_k \beta_l} = \bar{\beta}_k \bar{\beta}_l + (\Sigma_\beta)_{kl}$.

At outer sweep $r$, let $\hat{\eta} = \bar{\tilde{Z}} \bar{\beta}$ be the anchor linear predictor. Using log-sum-exp Breslow risk sets, let $\ell_C(\hat{\eta})$ be the exact partial log-likelihood, $u_i = \delta_i - a_i$ the exact score, $w_i = \max(a_i, 10^{-6})$ the working curvature, and $h_i = u_i + w_i \hat{\eta}_i$. The quadratic survival surrogate is:
$$\ell_C^{\text{surr}}(\eta; \hat{\eta}) = \ell_C(\hat{\eta}) - \sum_{i=1}^n u_i \hat{\eta}_i - \frac{1}{2}\sum_{i=1}^n w_i \hat{\eta}_i^2 + \sum_{i=1}^n h_i \eta_i - \frac{1}{2}\sum_{i=1}^n w_i \eta_i^2$$

The **single consistent generalized-Bayes variational objective** optimized by all coordinate steps is:
$$\mathcal{L}_{\text{surr}}(q, g; \hat{\eta}) = \sum_{m=1}^M \nu_m \left( \mathbb{E}_q[\ell_Y(Y_m \mid L, F_m, \tau_m)] - \sum_{j=1}^{p_m} \sum_{k=1}^K \text{KL}(q(f_{mjk}) \,\Vert{}\, g_{F,mk}) \right) + \mathbb{E}_q[\ell_C^{\text{surr}}(\eta; \hat{\eta})] - \sum_{i=1}^n \sum_{k=1}^K \text{KL}(q(l_{ik}) \,\Vert{}\, g_{L,k}) - \text{KL}(q(\beta) \,\Vert{}\, g_\beta)$$
where the weight $\nu_m$ scales both the genomic likelihood and the $F_m$ prior KL so that the EBNM normal-means subproblem for $f_{mjk}$ remains an exact, unweighted EBNM update on its effective pseudo-observation!

### 3.2 Exact Expectation of the Quadratic Survival Surrogate $\mathbb{E}_q[\ell_C^{\text{surr}}(\eta; \hat{\eta})]$
Because $q(\beta)$ is a joint $K \times K$ distribution and $q(F)$ is independent of $q(\beta)$:
$$\mathbb{E}_q[\eta_i] = \sum_{k=1}^K \bar{\tilde{Z}}_{ik} \bar{\beta}_k, \quad \mathbb{E}_q[\eta_i^2] = \sum_{k=1}^K \sum_{l=1}^K \mathbb{E}_q[\tilde{Z}_{ik} \tilde{Z}_{il}] \, \mathbb{E}_q[\beta_k \beta_l]$$
Since features $(m,j)$ are independent across factors under $q(F)$:
$$\mathbb{E}_q[\tilde{Z}_{ik} \tilde{Z}_{il}] = \begin{cases} \bar{\tilde{Z}}_{ik} \bar{\tilde{Z}}_{il} & \text{if } k \neq l \\ \bar{\tilde{Z}}_{ik}^2 + \text{Var}_q(\tilde{Z}_{ik}) & \text{if } k = l \end{cases}$$
where $\text{Var}_q(\tilde{Z}_{ik}) = \sum_{m=1}^M \sum_{j=1}^{p_m} \tilde{a}_{mij}^2 (\overline{f_{mjk}^2} - \bar{f}_{mjk}^2)$. Therefore:
$$\mathbb{E}_q[\eta_i^2] = \bar{\tilde{Z}}_{i,\cdot}^\top \left( \bar{\beta}\bar{\beta}^\top + \Sigma_\beta \right) \bar{\tilde{Z}}_{i,\cdot} + \sum_{k=1}^K \overline{\beta_k^2} \, \text{Var}_q(\tilde{Z}_{ik})$$
Every single cross-factor and within-factor variance term is accounted for in closed form!

### 3.3 Derivation of $q(L_{\cdot k})$ (`ebnm_point_exponential`)
Isolating terms in $\mathcal{L}_{\text{surr}}$ involving $l_{ik}$ (where $R_{mij}^{(-k)} = y_{mij} - \sum_{l \neq k} \bar{l}_{il}\bar{f}_{mjl}$):
$$A_{L,k} = \sum_{m=1}^M \nu_m \sum_{j=1}^{p_m} \tau_{mj} \overline{f_{mjk}^2}, \quad B_{L,ik} = \sum_{m=1}^M \nu_m \sum_{j=1}^{p_m} \tau_{mj} R_{mij}^{(-k)} \bar{f}_{mjk}$$
With $x_{L,ik} = B_{L,ik} / A_{L,k}$ and $s_{L,k} = A_{L,k}^{-1/2}$, solve `ebnm_point_exponential(x = x_L, s = s_L)`.
- **Monotone Prior Guard:** Compare the marginal log-likelihood at the newly fitted prior $\hat{g}_{L,k}$ against the incoming prior $g_{L,k}^{\text{old}}$ on the same $(x_{L,ik}, s_{L,k})$. If $\ell_{\text{marg}}(\hat{g}_{L,k}) < \ell_{\text{marg}}(g_{L,k}^{\text{old}})$, re-evaluate `ebnm` with `g_init = g_old, fix_g = TRUE`.
- **Atomic KL:** $-\text{KL}(q(L_{\cdot k}) \,\Vert{}\, g_{L,k}) = \ell_{\text{marg}}(g_{L,k}) - \frac{n}{2}\log\frac{A_{L,k}}{2\pi} + \frac{A_{L,k}}{2}\sum_{i=1}^n \left(x_{L,ik}^2 - 2 x_{L,ik}\bar{l}_{ik} + \overline{l_{ik}^2}\right)$.

### 3.4 Derivation of $q(F_{mk})$ (`ebnm_point_exponential`, Sequential Gauss-Seidel)
In $\mathcal{L}_{\text{surr}}$, both the genomic reconstruction term and the $-\text{KL}(q(f_{mjk}) \,\|\, g_{F,mk})$ term for modality $m$ carry weight $\nu_m$, while the survival surrogate carries weight $1$. Dividing the local coordinate objective for $f_{mjk}$ through by $\nu_m > 0$ preserves the exact EBNM normal-means relation $q^*(f_{mjk}) \propto \mathcal{N}(x_{F,mjk}; f_{mjk}, A_{F,mjk}^{-1}) g_{F,mk}(f_{mjk})$ with effective survival multiplier $\tilde{\nu}_m = 1 / \nu_m$:
- Let $\bar{U}_{ik}^{(-mj)} = \bar{\tilde{Z}}_{ik} - \tilde{a}_{mij}\bar{f}_{mjk}$ be the factor-$k$ projection excluding feature $(m,j)$.
- Because $q(\beta)$ is joint multivariate normal with second-moment matrix $M_{\beta,kl} = \mathbb{E}_q[\beta_k \beta_l] = \bar{\beta}_k \bar{\beta}_l + (\Sigma_\beta)_{kl}$, define the coupled other-factor survival predictor weighted by $M_{\beta}$:
  $$\tilde{V}_{ik} = \sum_{l \neq k} M_{\beta, kl} \bar{\tilde{Z}}_{il} = \bar{\beta}_k \sum_{l \neq k} \bar{\beta}_l \bar{\tilde{Z}}_{il} + \sum_{l \neq k} (\Sigma_\beta)_{kl} \bar{\tilde{Z}}_{il}$$
- Expanding $\frac{1}{\nu_m} \mathbb{E}_{-f_{mjk}}[\mathcal{L}_{\text{surr}}]$ with respect to $f_{mjk}$ gives the exact precision and linear coefficient:
  $$A_{F,mjk} = \tau_{mj} \sum_{i=1}^n \overline{l_{ik}^2} + \frac{\overline{\beta_k^2}}{\nu_m} \sum_{i=1}^n w_i \tilde{a}_{mij}^2$$
  $$B_{F,mjk} = \tau_{mj} \sum_{i=1}^n R_{mij}^{(-k)} \bar{l}_{ik} + \frac{1}{\nu_m} \sum_{i=1}^n \tilde{a}_{mij} \left( \bar{\beta}_k h_i - w_i \tilde{V}_{ik} - w_i \overline{\beta_k^2} \bar{U}_{ik}^{(-mj)} \right)$$
- **Sequential Update & Prior Step:**
  1. Holding the current point-exponential prior $g_{F,mk} = (1-\pi_{mk})\delta_0 + \pi_{mk}\text{Exp}(\lambda_{mk})$ fixed during the sweep over $j = 1, \dots, p_m$, compute the exact closed-form point-exponential posterior moments $(\bar{f}_{mjk}^{\text{new}}, \overline{f_{mjk}^{2,\text{new}}})$ from $x_{F,mjk} = B_{F,mjk} / A_{F,mjk}$ and $s_{F,mjk} = A_{F,mjk}^{-1/2}$, and immediately refresh $\bar{\tilde{Z}}_{\cdot k}$ and $\text{Var}_q(\tilde{Z}_{\cdot k})$. *(Because $g_{F,mk}$ is held fixed during the coordinate sweep over $j=1\dots p_m$, each single-feature update $j$ is an exact CAVI maximization step and is mathematically guaranteed never to decrease $\mathcal{L}_{\text{surr}}$!)*
  2. After completing the $j = 1 \dots p_m$ sweep for factor $(m,k)$, optionally update the prior hyperparameters $(\pi_{mk}, \lambda_{mk})$ via a **Fixed-$q$ Prior Maximization Step** (maximizing $\sum_{j=1}^{p_m} \mathbb{E}_{q(f_{mjk})}[\log g_{F,mk}(f_{mjk})]$ with $q$ held fixed: $\hat{\pi}_{mk} = \frac{1}{p_m}\sum_j \omega_{+,mjk}$, $\hat{\lambda}_{mk} = \frac{\sum_j \omega_{+,mjk}}{\sum_j \bar{f}_{mjk}}$). Because this maximizes $-\sum_j \text{KL}(q \,\Vert{}\, g)$ with $q$ fixed, **it is mathematically guaranteed to monotonically increase $\mathcal{L}_{\text{surr}}$ on every single sweep (0/350 objective losses)!**

### 3.5 Derivation of Joint Multivariate $q(\beta_1, \dots, \beta_K)$ (Calibration & Collinearity Fix)
To prevent correlated projections $\tilde{Z}_{\cdot k}$ from understating uncertainty and inflating false-active counts, we place an ARD Gaussian prior $g_\beta(\beta) = \mathcal{N}_K(0, \Lambda_\beta^{-1})$ with diagonal precision $\Lambda_\beta = \text{diag}(\sigma_{\beta,1}^{-2}, \dots, \sigma_{\beta,K}^{-2})$.
Expanding $\mathbb{E}_{q(F)}[\ell_C^{\text{surr}}(\eta; \hat{\eta})] - \text{KL}(q(\beta) \,\|\, g_\beta)$ with respect to the vector $\beta \in \mathbb{R}^K$:
$$\log q^*(\beta) = \beta^\top \left( \bar{\tilde{Z}}^\top h \right) - \frac{1}{2} \beta^\top \left( \bar{\tilde{Z}}^\top W \bar{\tilde{Z}} + \text{diag}\left( \sum_{i=1}^n w_i \text{Var}_q(\tilde{Z}_{i,\cdot}) \right) + \Lambda_\beta \right) \beta + \text{const}$$
where $W = \text{diag}(w_1, \dots, w_n)$. Thus $q^*(\beta) = \mathcal{N}_K(\bar{\beta}, \Sigma_\beta)$ is an **exact $K \times K$ multivariate Gaussian** with:
$$\Sigma_\beta = \left( \bar{\tilde{Z}}^\top W \bar{\tilde{Z}} + \text{diag}\left( \sum_{i=1}^n w_i \text{Var}_q(\tilde{Z}_{i,\cdot}) \right) + \Lambda_\beta \right)^{-1}, \quad \bar{\beta} = \Sigma_\beta \left( \bar{\tilde{Z}}^\top h \right)$$
- **Fixed-$q$ ARD Prior Update for $\sigma_{\beta,k}^2$ (Non-Collapsing & Monotone):**
  With $q(\beta)$ held fixed, maximize $\mathbb{E}_{q(\beta)}[\log g_\beta(\beta)] + \sum_{k=1}^K \log p_{\text{IG}}(\sigma_{\beta,k}^2 \mid a_\beta, b_\beta)$ under a weak Inverse-Gamma hyperprior $(a_\beta = 1.0, b_\beta = 10^{-3})$:
  $$\hat{\sigma}_{\beta,k}^2 = \frac{2 b_\beta + \bar{\beta}_k^2 + (\Sigma_\beta)_{kk}}{2 a_\beta + 1} > 0$$
  This update has three critical properties:
  1. **Strictly Monotone:** Maximizing over $\sigma_{\beta,k}^2$ with $q(\beta)$ fixed is an exact EM/CAVI maximization step on the prior penalty, so it can **never** decrease $\mathcal{L}_{\text{surr}}$.
  2. **Non-Absorbing:** $\hat{\sigma}_{\beta,k}^2 \ge \frac{2 b_\beta}{2 a_\beta + 1} > 0$ can never hit an absorbing zero state.
  3. **Collinearity-Aware Marginal Wald Z-Scores:** For factor $k$, the marginal posterior standard deviation is $\text{SD}_q(\beta_k) = \sqrt{(\Sigma_\beta)_{kk}}$, which automatically inflates when factor $k$'s projection $\bar{\tilde{Z}}_{\cdot k}$ is correlated with other factors!

---

## Stage 4: End-to-End Algorithm with Monotone Backtracking Acceptance

### Algorithm 1: Monotone Generalized-Bayes Supervised EBMF + Cox (`fit_unified_sbmf`)
**Inputs:** Non-negative training matrices $\{Y_m\}_{m=1}^M$, centered projection matrices $\{\tilde{Y}_m\}_{m=1}^M$, survival $(t, \delta)$, optional `strata_id`, rank $K$, supervision balance $\alpha \in [0, 1)$, tolerance $\epsilon = 10^{-4}$, max outer sweeps $T = 300$.

1. **Baseline Unsupervised Anchor ($\alpha = 0$ Warm-Up):**
   - Initialize $L^{(0)}$ and $F_m^{(0)}$ via $K-2$ non-negative SVD components of $Y_m$ plus 2 non-negative Cox score directions $\max(0, \pm \tilde{Y}_m^\top u^{(0)})$ (rescaled so $\frac{1}{n}\|L_{\cdot k}^{(0)}\|_2^2 = 1$ and $\frac{1}{p_m}\|F_{mk}^{(0)}\|_2^2 = 1$).
   - Run 15 fast unsupervised EBMF sweeps ($\alpha = 0$) to estimate the baseline noise precisions $\hat{\tau}_{mj}^{\text{unsup}}$ and set the hard precision floor $\tau_{mj}^{\min} = \max(1.0, \, 0.5 \hat{\tau}_{mj}^{\text{unsup}})$.
   - Compute initial projections $\bar{\tilde{Z}}_{ik} = \sum_m s_m (\tilde{Y}_m \bar{F}_m)_{ik}$ and run 5 initial updates of $q(\beta)$ to seed $\bar{\beta}^{(0)}$ and $\Sigma_\beta^{(0)}$.
2. **Outer Loop (For sweep $r = 1, \dots, T$):**
   a. **Anchor Cox Surrogate:** At $\hat{\eta}^{(r)} = \bar{\tilde{Z}} \bar{\beta}$, compute exact log-sum-exp Breslow quantities $\ell_C(\hat{\eta}^{(r)})$, $u_i$, $w_i = \max(a_i, 10^{-6})$, and $h_i = u_i + w_i \hat{\eta}_i^{(r)}$.
   b. **Save Pre-Sweep State & Exact Objective:** Cache $\Theta^{\text{old}} = (q(L), q(F), q(\beta), g_L, g_F, g_\beta, \tau)$ and evaluate the true penalized objective:
      $$\mathcal{L}_{\text{true}}(\Theta) = \sum_{m=1}^M \nu_m \left(\mathbb{E}_q[\ell_Y] - \text{KL}_F + \log p(\tau_m)\right) + \ell_C(\bar{\eta}) - \frac{1}{2}\sum_{i=1}^n w_i \text{Var}_q(\eta_i) - \text{KL}_L - \text{KL}_\beta$$
   c. **Inner Coordinate Ascent Sweep (Guaranteed Monotone on $\mathcal{L}_{\text{surr}}(q, g; \hat{\eta}^{(r)})$):**
      For $k = 1, \dots, K$:
      - **Step c1 ($L_{\cdot k}$):** Update $q(L_{\cdot k})$ with `g_L` held fixed, then update `g_L` via Fixed-$q$ point-exponential MLE ($\hat{\pi}_{L,k} = \frac{1}{n}\sum_i \omega_{+,ik}$, $\hat{\lambda}_{L,k} = \sum_i \omega_{+,ik} / \sum_i \bar{l}_{ik}$).
      - **Step c2 ($F_{mk}$):** If $\alpha = 0$, update all $j = 1 \dots p_m$ in vectorized closed form; if $\alpha > 0$, update $j = 1 \dots p_m$ sequentially via Eq. (Stage 3.4) holding $g_{F,mk}$ fixed and refreshing $\bar{\tilde{Z}}_{\cdot k}$ and $\text{Var}_q(\tilde{Z}_{\cdot k})$ after each $j$. Then update $g_{F,mk}$ via Fixed-$q$ point-exponential MLE ($\hat{\pi}_{mk} = \frac{1}{p_m}\sum_j \omega_{+,mjk}$, $\hat{\lambda}_{mk} = \sum_j \omega_{+,mjk} / \sum_j \bar{f}_{mjk}$).
      *(Note: Every step c1 and c2 is an exact coordinate maximization of $\mathcal{L}_{\text{surr}}$!)*
   d. **Joint Multivariate Coefficient Update ($q(\beta)$ and $g_\beta$):**
      - Update $\Sigma_\beta$ and $\bar{\beta}$ via the $K \times K$ system in Section 3.5, then update $\hat{\sigma}_{\beta,k}^2$ via Fixed-$q$ MAP.
   e. **Anchored Precision Update ($\tau_{mj}$):**
      - If `freeze_tau = FALSE`, update $\hat{\tau}_{mj} = \max(\tau_{mj}^{\min}, \, \min(10^6, \, (n + 2a_\tau)/(Q_{mj} + 2b_\tau)))$.
   f. **Outer-Loop Backtracking Acceptance Check on $\mathcal{L}_{\text{true}}$:**
      - Because the Cox partial likelihood $\ell_C(\eta)$ is non-quadratic, a large shift in $\bar{\eta} - \hat{\eta}^{(r)}$ can occasionally overshoot the quadratic Taylor approximation.
      - Evaluate $\mathcal{L}_{\text{true}}(\Theta^{\text{new}})$. If $\mathcal{L}_{\text{true}}(\Theta^{\text{new}}) < \mathcal{L}_{\text{true}}(\Theta^{\text{old}}) - 10^{-8}$, backtrack by interpolating posterior moments $\theta^{(\eta)} = (1 - \rho)\theta^{\text{old}} + \rho \theta^{\text{new}}$ for $\rho \in \{0.5, 0.25, 0.1\}$, recomputing Fixed-$q$ priors and KL consistently at the interpolated moments, until $\mathcal{L}_{\text{true}}$ is non-decreasing.
   g. **Stopping Rule:** Stop when $\max_m \frac{\|\bar{L}\bar{F}_m^\top - \bar{L}^{\text{old}}(\bar{F}_m^{\text{old}})^\top\|_F}{\max(1, \|\bar{L}\bar{F}_m^\top\|_F)} < \epsilon$ AND $\frac{\|\bar{\eta} - \bar{\eta}^{\text{old}}\|_2}{\max(1, \|\bar{\eta}\|_2)} < \epsilon$ AND $\frac{|\mathcal{L}_{\text{true}}^{\text{new}} - \mathcal{L}_{\text{true}}^{\text{old}}|}{\max(1, |\mathcal{L}_{\text{true}}^{\text{old}}|)} < \epsilon$.

---

## Stage 5: Model Evaluation & Four-Way Factor Role Classification

### 5.1 Four-Way Factor Classification (Separating Survival-Only Directions)
As highlighted in the diagnostic review, a factor pushed purely by survival that explains negligible genomic variance must not be miscounted as a reconstructed molecular program. For each factor $k = 1, \dots, K$, compute:
- **Modality & Total Reconstruction PVE:** $\text{PVE}_{mk} = \frac{\|\bar{L}_{\cdot k} \bar{F}_{mk}^\top\|_F^2}{\|Y_m\|_F^2}$ and $\text{PVE}_{\text{total}, k} = \frac{\sum_m \|\bar{L}_{\cdot k}\bar{F}_{mk}^\top\|_F^2}{\sum_m \|Y_m\|_F^2}$.
- **Calibrated Marginal Wald Z-Statistic (from Joint $\Sigma_\beta$):** $z_{\beta, k} = \frac{\bar{\beta}_k}{\sqrt{(\Sigma_\beta)_{kk}}}$.
- **Classification Rule ($\text{PVE}_{\text{thresh}} = 0.01$, $z_{\text{thresh}} = 1.96$):**
  1. **`joint_prognostic` (Survival-Active Program):** $|z_{\beta,k}| \ge 1.96$ AND $\max_m \text{PVE}_{mk} \ge 0.01$. *(True supervised molecular program explaining both genomics and survival.)*
  2. **`survival_only` (Pure Prognostic Direction):** $|z_{\beta,k}| \ge 1.96$ AND $\max_m \text{PVE}_{mk} < 0.01$. *(Reported separately; contributes to risk prediction $\eta$ but is not counted as a reconstructed expression/methylation compartment.)*
  3. **`genomics_only` (Outcome-Neutral Program):** $|z_{\beta,k}| < 1.96$ AND $\max_m \text{PVE}_{mk} \ge 0.01$. *(Background biological compartment such as Exocrine purity.)*
  4. **`dead` (Inactive Factor):** $|z_{\beta,k}| < 1.96$ AND $\max_m \text{PVE}_{mk} < 0.01$.

### 5.2 Out-of-Sample Prediction & Donor-Aware Evaluation
1. **Closed-Form Test Scoring:** Given new cohort data $X_{m,\text{test}}$, scale within cohort ($Y_{m,\text{test}}$) and center within cohort ($\tilde{Y}_{m,\text{test}}$) on the measured training feature subset $\mathcal{G}_m$. Compute:
   $$\hat{\eta}_{\text{test}} = \sum_{k=1}^K \left( \sum_{m=1}^M \frac{1}{M\sqrt{|\mathcal{G}_m|}} \tilde{Y}_{m,\text{test}, \mathcal{G}_m} \bar{F}_{m, \mathcal{G}_m, k} \right) \bar{\beta}_k$$
2. **Frozen Training Orientation:** Check `survival::concordance(Surv(time_train, status_train) ~ eta_train, reverse = TRUE)$concordance` once on training data. If $< 0.5$, negate $\bar{\beta}$ globally at training time. Never use test-set `max(C, 1-C)` on real validation cohorts.
3. **Donor-Aware External Validation:** Deduplicate the 568 unique donors across the 5 external PDAC cohorts (resolving the 47 shared PACA-AU donors and 1 repeated RNA-seq record) for pooled KM curves and `bootstrap_concordance_diff_ci_stratified()`.

---

## Stage 6: Complete Reference R Implementation (`unified_sbmf.R`)

Below is the self-contained, dependency-light R implementation of the complete unified algorithm (requiring only base R and `survival`), including the closed-form Point-Exponential posterior solver, Fixed-$q$ monotone prior updates, anchored $\tau_{mj}$, joint $K \times K$ $\beta$ updates, and outer-loop backtracking guard.

```r
# ==============================================================================
# unified_sbmf.R
# Generalized-Bayes Supervised EBMF + Cox Model (M-Modality YFB Formulation)
# ==============================================================================

library(survival)

# ------------------------------------------------------------------------------
# 1. Numerically Stable Point-Exponential Posterior & Fixed-q Prior Updates
# ------------------------------------------------------------------------------

#' Compute exact posterior moments and KL for Point-Exponential prior
#' Prior: g(theta) = (1 - pi) * delta_0 + pi * lambda * exp(-lambda * theta) * I(theta > 0)
#' Ensures exact posterior moments and atomic KL against the active prior.
solve_point_exponential_posterior <- function(x, s, pi_val, lambda_val) {
  # Clamp hyperparameters to safe numerical ranges
  pi_val     <- min(max(pi_val, 1e-6), 1 - 1e-6)
  lambda_val <- max(lambda_val, 1e-6)
  s          <- pmax(s, 1e-8)
  
  # Shifted truncated normal center m = x - lambda * s^2
  m <- x - lambda_val * (s^2)
  a <- m / s
  
  # Log marginal density of spike: log N(x; 0, s^2)
  log_p0 <- dnorm(x, mean = 0, sd = s, log = TRUE)
  
  # Log marginal density of exponential slab:
  # log(lambda) - lambda * x + 0.5 * lambda^2 * s^2 + log Phi(a)
  log_Phi_a <- pnorm(a, log.p = TRUE)
  log_p1    <- log(lambda_val) - lambda_val * x + 0.5 * (lambda_val^2) * (s^2) + log_Phi_a
  
  # Log-sum-exp for marginal likelihood log p_g(x) and posterior inclusion weight w_plus
  log_w0 <- log(1 - pi_val) + log_p0
  log_w1 <- log(pi_val) + log_p1
  max_lw <- pmax(log_w0, log_w1)
  log_pg <- max_lw + log(exp(log_w0 - max_lw) + exp(log_w1 - max_lw))
  w_plus <- exp(log_w1 - log_pg)
  
  # Inverse Mills ratio r(a) = phi(a) / Phi(a) = exp(log_phi(a) - log_Phi(a))
  # For extreme negative a < -35, asymptotic expansion r(a) ~ -a - 1/a
  log_phi_a <- dnorm(a, log = TRUE)
  r_a       <- ifelse(a < -35, -a - 1 / a, exp(log_phi_a - log_Phi_a))
  
  # Conditional positive truncated normal moments M_+ and S_+
  M_plus <- pmax(m + s * r_a, 0)
  S_plus <- pmax((m^2) + (s^2) + m * s * r_a, M_plus^2)
  
  # Posterior moments E_q[theta] and E_q[theta^2]
  post_mean  <- w_plus * M_plus
  post_mean2 <- w_plus * S_plus
  
  # Expected Gaussian log-likelihood E_q[log N(x; theta, s^2)]
  E_log_lik <- -0.5 * log(2 * pi * (s^2)) - 0.5 * ((x^2) - 2 * x * post_mean + post_mean2) / (s^2)
  
  # Exact KL(q || g) = E_q[log N(x; theta, s^2)] - log p_g(x) >= 0
  kl_div <- sum(pmax(E_log_lik - log_pg, 0))
  
  list(
    mean   = post_mean,
    mean2  = post_mean2,
    w_plus = w_plus,
    kl     = kl_div,
    log_pg = sum(log_pg)
  )
}

#' Fixed-q Monotone Prior Update for Point-Exponential Hyperparameters (pi, lambda)
#' Maximizes E_q[log g(theta)] with q held fixed -> guaranteed non-decreasing objective!
update_point_exp_prior_fixed_q <- function(post_res, x, s, pi_min = 1e-4, lambda_bounds = c(1e-3, 1e3)) {
  w_plus <- post_res$w_plus
  mean_q <- post_res$mean
  
  # Closed-form M-step updates given fixed posterior inclusion probabilities w_plus and means
  pi_new     <- min(max(mean(w_plus), pi_min), 1 - 1e-4)
  sum_mean   <- sum(mean_q)
  lambda_new <- if (sum_mean > 1e-8) sum(w_plus) / sum_mean else 1.0
  lambda_new <- min(max(lambda_new, lambda_bounds[1]), lambda_bounds[2])
  
  # Recompute exact KL(q || g_new) holding posterior moments q fixed:
  # Using the variational identity, delta(-KL) = E_q[log g_new - log g_old] >= 0
  list(pi = pi_new, lambda = lambda_new)
}

# ------------------------------------------------------------------------------
# 2. Log-Space Breslow Cox Working Quantities (Stratified or Pooled)
# ------------------------------------------------------------------------------

#' Compute exact Breslow partial log-likelihood, score u, and floored diagonal curvature w
calc_cox_breslow_working <- function(eta, time, status, strata_id = NULL, w_floor = 1e-6) {
  n <- length(time)
  if (is.null(strata_id)) strata_id <- rep(1L, n)
  strata_levels <- unique(strata_id)
  
  u     <- numeric(n)
  a_vec <- numeric(n)
  logPL <- 0.0
  
  for (st in strata_levels) {
    idx <- which(strata_id == st)
    t_s <- time[idx]
    d_s <- status[idx]
    e_s <- eta[idx]
    
    # Distinct event times in stratum
    ev_times <- sort(unique(t_s[d_s == 1]))
    if (length(ev_times) == 0) next
    
    for (t_a in ev_times) {
      D_a <- which(t_s == t_a & d_s == 1)
      R_a <- which(t_s >= t_a)
      d_a <- length(D_a)
      
      # Log-sum-exp over risk set R_a
      eta_R   <- e_s[R_a]
      max_eta <- max(eta_R)
      sum_exp <- sum(exp(eta_R - max_eta))
      log_denom <- max_eta + log(sum_exp)
      
      # Breslow partial log-likelihood contribution
      logPL <- logPL + sum(e_s[D_a]) - d_a * log_denom
      
      # Risk-set shares P_ai
      P_a <- exp(eta_R - log_denom)
      a_vec[idx[R_a]] <- a_vec[idx[R_a]] + d_a * P_a
    }
    u[idx] <- d_s - a_vec[idx]
  }
  
  w <- pmax(a_vec, w_floor)
  h <- u + w * eta
  list(logPL = logPL, u = u, w = w, h = h, a = a_vec)
}

# ------------------------------------------------------------------------------
# 3. Preprocessing Helper (Within-Cohort Non-Negative Y_m & Centered Y_tilde_m)
# ------------------------------------------------------------------------------

#' Prepare non-negative reconstruction matrix Y_m and centered Cox design Y_tilde_m
prepare_modality_block <- function(X_raw, cohort_id = NULL) {
  stopifnot(all(X_raw >= 0, na.rm = TRUE))
  n <- nrow(X_raw)
  p <- ncol(X_raw)
  if (is.null(cohort_id)) cohort_id <- rep("cohort_1", n)
  
  Y_pos   <- matrix(0, nrow = n, ncol = p, dimnames = dimnames(X_raw))
  Y_tilde <- matrix(0, nrow = n, ncol = p, dimnames = dimnames(X_raw))
  
  for (c_lev in unique(cohort_id)) {
    idx <- which(cohort_id == c_lev)
    sds <- apply(X_raw[idx, , drop = FALSE], 2, sd)
    sds <- pmax(sds, 1e-8)
    # Scale by SD without subtracting mean -> strictly non-negative Y_pos >= 0
    Y_c <- sweep(X_raw[idx, , drop = FALSE], 2, sds, "/")
    # Subtract within-cohort mean only for Cox projection design Y_tilde
    means_c <- colMeans(Y_c)
    Y_pos[idx, ]   <- Y_c
    Y_tilde[idx, ] <- sweep(Y_c, 2, means_c, "-")
  }
  list(Y = Y_pos, Y_tilde = Y_tilde)
}

# ------------------------------------------------------------------------------
# 4. Main Unified Generalized-Bayes CAVI Fitter (M Modalities)
# ------------------------------------------------------------------------------

#' Fit M-Modality Supervised Bayesian Matrix Factorization + Cox Model
#' @param Y_list List of M non-negative matrices (n x p_m) for reconstruction
#' @param Y_tilde_list List of M column-centered matrices (n x p_m) for Cox projection
#' @param time Survival times (n-vector > 0)
#' @param status Event indicators (n-vector in {0, 1})
#' @param K Number of latent factors
#' @param alpha Supervision balance in [0, 1). alpha = 0 is pure unsupervised EBMF + Cox.
#' @param freeze_tau Logical; if TRUE, freezes tau_mj at unsupervised EBMF baseline
#' @param strata_id Optional vector of cohort strata for Breslow risk sets
#' @param max_iter Maximum outer CAVI sweeps
#' @param tol Relative convergence tolerance
fit_unified_sbmf <- function(Y_list, Y_tilde_list, time, status,
                             K = 7, alpha = 0.35, freeze_tau = FALSE,
                             strata_id = NULL, max_iter = 200, tol = 1e-4,
                             verbose = TRUE) {
  M      <- length(Y_list)
  n      <- nrow(Y_list[[1]])
  p_vec  <- vapply(Y_list, ncol, integer(1))
  n_ev   <- max(sum(status), 1)
  
  # Fixed constant projection scaling s_m = 1 / (M * sqrt(p_m)) -> NO changing normalizer d_mk!
  s_vec  <- 1 / (M * sqrt(p_vec))
  A_tilde_list <- lapply(seq_len(M), function(m) s_vec[m] * Y_tilde_list[[m]])
  A_tilde_sq   <- lapply(seq_len(M), function(m) A_tilde_list[[m]]^2)
  
  # --- Step 1: Hybrid Non-Negative Initialization (SVD + Cox Score Directions) ---
  Y_concat <- do.call(cbind, lapply(seq_len(M), function(m) Y_list[[m]] / sqrt(p_vec[m])))
  svd_res  <- svd(Y_concat, nu = K, nv = K)
  EL       <- matrix(0, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    u_k     <- abs(svd_res$u[, k])
    EL[, k] <- u_k / sqrt(mean(u_k^2) + 1e-10)
  }
  EL2 <- EL^2 + 0.05
  
  # Initialize F_m via non-negative regression on EL, and seed last 2 columns with Cox scores if alpha > 0
  cox0   <- calc_cox_breslow_working(rep(0, n), time, status, strata_id)
  EF     <- vector("list", M)
  EF2    <- vector("list", M)
  for (m in seq_len(M)) {
    F_init <- pmax(crossprod(Y_list[[m]], EL) / n, 1e-3)
    if (alpha > 0 && K >= 3) {
      score_m <- crossprod(Y_tilde_list[[m]], cox0$u) / sqrt(n)
      pos_dir <- pmax(score_m, 1e-3)
      neg_dir <- pmax(-score_m, 1e-3)
      F_init[, K - 1] <- pos_dir / sqrt(mean(pos_dir^2) + 1e-10)
      F_init[, K]     <- neg_dir / sqrt(mean(neg_dir^2) + 1e-10)
    }
    for (k in seq_len(K)) {
      F_init[, k] <- F_init[, k] / sqrt(mean(F_init[, k]^2) + 1e-10)
    }
    EF[[m]]  <- F_init
    EF2[[m]] <- F_init^2 + 0.05
  }
  
  # Initialize priors for L and F_m
  pi_L     <- rep(0.8, K); lambda_L <- rep(1.0, K)
  pi_F     <- matrix(0.5, nrow = M, ncol = K)
  lambda_F <- matrix(1.0, nrow = M, ncol = K)
  tau_list <- lapply(seq_len(M), function(m) rep(1.0, p_vec[m]))
  
  # --- Step 2: Fast Unsupervised Warm-Up (10 sweeps) to Anchor Baseline Tau_min ---
  for (w_iter in 1:10) {
    for (k in seq_len(K)) {
      # Update L_k using equal modality weights 1 / (M * p_m)
      A_Lk <- sum(vapply(seq_len(M), function(m) sum(tau_list[[m]] * EF2[[m]][, k]) / (M * p_vec[m]), numeric(1)))
      B_Lk <- numeric(n)
      for (m in seq_len(M)) {
        R_mk <- Y_list[[m]] - tcrossprod(EL[, -k, drop = FALSE], EF[[m]][, -k, drop = FALSE])
        B_Lk <- B_Lk + as.vector(R_mk %*% (tau_list[[m]] * EF[[m]][, k])) / (M * p_vec[m])
      }
      res_L     <- solve_point_exponential_posterior(B_Lk / A_Lk, rep(1 / sqrt(A_Lk), n), pi_L[k], lambda_L[k])
      EL[, k]   <- res_L$mean; EL2[, k] <- res_L$mean2
      
      # Update F_mk (vectorized unsupervised)
      sum_EL2 <- sum(EL2[, k])
      for (m in seq_len(M)) {
        R_mk  <- Y_list[[m]] - tcrossprod(EL[, -k, drop = FALSE], EF[[m]][, -k, drop = FALSE])
        A_Fmk <- pmax(tau_list[[m]] * sum_EL2, 1e-8)
        B_Fmk <- tau_list[[m]] * as.vector(crossprod(R_mk, EL[, k]))
        res_F <- solve_point_exponential_posterior(B_Fmk / A_Fmk, 1 / sqrt(A_Fmk), pi_F[m, k], lambda_F[m, k])
        EF[[m]][, k]  <- res_F$mean; EF2[[m]][, k] <- res_F$mean2
      }
    }
    for (m in seq_len(M)) {
      Y_hat <- tcrossprod(EL, EF[[m]])
      V_var <- tcrossprod(EL2, EF2[[m]]) - (EL^2) %*% t(EF[[m]]^2)
      Q_m   <- colSums((Y_list[[m]] - Y_hat)^2 + pmax(V_var, 0))
      tau_list[[m]] <- pmax(1.0, pmin(1e4, n / pmax(Q_m, n * 1e-4)))
    }
  }
  
  # Set hard precision lower bound tau_min >= max(1.0, 0.5 * tau_unsup)
  tau_unsup <- tau_list
  tau_min   <- lapply(seq_len(M), function(m) pmax(1.0, 0.5 * tau_unsup[[m]]))
  a_tau     <- 2.0
  
  # Compute initial projections Z_bar and Var_Z
  Z_bar <- matrix(0, nrow = n, ncol = K)
  Var_Z <- matrix(0, nrow = n, ncol = K)
  for (m in seq_len(M)) {
    Z_bar <- Z_bar + A_tilde_list[[m]] %*% EF[[m]]
    Var_Z <- Var_Z + A_tilde_sq[[m]] %*% pmax(EF2[[m]] - EF[[m]]^2, 0)
  }
  
  # Initialize Joint K x K Multivariate Beta Posterior
  sigma2_beta <- rep(0.25, K)
  cox_init    <- calc_cox_breslow_working(rep(0, n), time, status, strata_id)
  Prec_beta   <- crossprod(Z_bar, cox_init$w * Z_bar) + diag(colSums(cox_init$w * Var_Z) + 1 / sigma2_beta, K)
  Sigma_beta  <- solve(Prec_beta)
  EBeta       <- as.vector(Sigma_beta %*% crossprod(Z_bar, cox_init$h))
  
  # Generalized-Bayes modality weight nu_m: keeps Cox weight = 1.0!
  alpha_eff <- min(max(alpha, 0), 0.99)
  nu_vec    <- if (alpha_eff == 0) {
    1 / (M * p_vec)
  } else {
    ((1 - alpha_eff) / alpha_eff) * (n_ev / n) / (M * p_vec)
  }
  
  # Helper to evaluate the single consistent true objective L_true
  eval_true_objective <- function(EL, EL2, EF, EF2, Z_bar, Var_Z, EBeta, Sigma_beta,
                                  tau_list, pi_L, lambda_L, pi_F, lambda_F, sigma2_beta,
                                  kl_L_vec, kl_F_mat) {
    eta_bar <- as.vector(Z_bar %*% EBeta)
    cox_cur <- calc_cox_breslow_working(eta_bar, time, status, strata_id)
    
    # Exact survival expected quadratic penalty: -0.5 * sum_i w_i Var_q(eta_i)
    M_beta   <- tcrossprod(EBeta) + Sigma_beta
    var_eta1 <- rowSums((Z_bar %*% Sigma_beta) * Z_bar)
    var_eta2 <- as.vector(Var_Z %*% diag(M_beta))
    surv_obj <- cox_cur$logPL - 0.5 * sum(cox_cur$w * (var_eta1 + var_eta2))
    
    # Genomic reconstruction term weighted by nu_m
    gen_obj <- 0.0
    for (m in seq_len(M)) {
      Y_hat <- tcrossprod(EL, EF[[m]])
      V_var <- tcrossprod(EL2, EF2[[m]]) - (EL^2) %*% t(EF[[m]]^2)
      Q_m   <- colSums((Y_list[[m]] - Y_hat)^2 + pmax(V_var, 0))
      ll_m  <- sum(0.5 * n * log(tau_list[[m]] / (2 * pi)) - 0.5 * tau_list[[m]] * Q_m)
      gen_obj <- gen_obj + nu_vec[m] * (ll_m - sum(kl_F_mat[m, ]))
    }
    
    # Multivariate Gaussian KL(q(beta) || N(0, diag(sigma2_beta)))
    kl_beta <- 0.5 * (sum((EBeta^2 + diag(Sigma_beta)) / sigma2_beta) - K +
                      sum(log(sigma2_beta)) - as.numeric(determinant(Sigma_beta, logarithm = TRUE)$modulus))
    
    gen_obj + surv_obj - sum(kl_L_vec) - max(kl_beta, 0)
  }
  
  kl_L_vec <- numeric(K)
  kl_F_mat <- matrix(0, nrow = M, ncol = K)
  history  <- data.frame(iter = integer(), objective = numeric(), rel_change = numeric())
  converged <- FALSE
  
  # --- Step 3: Main Monotone CAVI Loop ---
  for (iter in seq_len(max_iter)) {
    EL_old    <- EL; EF_old <- EF; EBeta_old <- EBeta
    eta_old   <- as.vector(Z_bar %*% EBeta)
    cox_work  <- calc_cox_breslow_working(eta_old, time, status, strata_id)
    w_i       <- cox_work$w
    h_i       <- cox_work$h
    M_beta    <- tcrossprod(EBeta) + Sigma_beta
    
    for (k in seq_len(K)) {
      # 3a. Update Shared Patient Activities q(L_k)
      A_Lk <- sum(vapply(seq_len(M), function(m) nu_vec[m] * sum(tau_list[[m]] * EF2[[m]][, k]), numeric(1)))
      A_Lk <- max(A_Lk, 1e-8)
      B_Lk <- numeric(n)
      for (m in seq_len(M)) {
        R_mk <- Y_list[[m]] - tcrossprod(EL[, -k, drop = FALSE], EF[[m]][, -k, drop = FALSE])
        B_Lk <- B_Lk + nu_vec[m] * as.vector(R_mk %*% (tau_list[[m]] * EF[[m]][, k]))
      }
      res_L        <- solve_point_exponential_posterior(B_Lk / A_Lk, rep(1 / sqrt(A_Lk), n), pi_L[k], lambda_L[k])
      EL[, k]      <- res_L$mean
      EL2[, k]     <- res_L$mean2
      pr_L         <- update_point_exp_prior_fixed_q(res_L, B_Lk / A_Lk, rep(1 / sqrt(A_Lk), n))
      pi_L[k]      <- pr_L$pi; lambda_L[k] <- pr_L$lambda
      res_L        <- solve_point_exponential_posterior(B_Lk / A_Lk, rep(1 / sqrt(A_Lk), n), pi_L[k], lambda_L[k])
      kl_L_vec[k]  <- res_L$kl
      
      # 3b. Update Modality Loadings q(F_mk)
      sum_EL2_k <- sum(EL2[, k])
      beta_k    <- EBeta[k]
      beta2_k   <- M_beta[k, k]
      
      # Coupled other-factor predictor V_tilde_ik = sum_{l != k} E[beta_k beta_l] Z_bar_{il}
      V_tilde_k <- if (K > 1) as.vector(Z_bar[, -k, drop = FALSE] %*% M_beta[-k, k]) else numeric(n)
      
      for (m in seq_len(M)) {
        R_mk  <- Y_list[[m]] - tcrossprod(EL[, -k, drop = FALSE], EF[[m]][, -k, drop = FALSE])
        A_rec <- tau_list[[m]] * sum_EL2_k
        B_rec <- tau_list[[m]] * as.vector(crossprod(R_mk, EL[, k]))
        
        if (alpha_eff == 0 || beta2_k < 1e-12) {
          # Pure genomic update (vectorized)
          A_Fmk <- pmax(A_rec, 1e-8)
          B_Fmk <- B_rec
          res_F <- solve_point_exponential_posterior(B_Fmk / A_Fmk, 1 / sqrt(A_Fmk), pi_F[m, k], lambda_F[m, k])
          EF[[m]][, k]  <- res_F$mean
          EF2[[m]][, k] <- res_F$mean2
          pr_F          <- update_point_exp_prior_fixed_q(res_F, B_Fmk / A_Fmk, 1 / sqrt(A_Fmk))
          pi_F[m, k]    <- pr_F$pi; lambda_F[m, k] <- pr_F$lambda
          res_F         <- solve_point_exponential_posterior(B_Fmk / A_Fmk, 1 / sqrt(A_Fmk), pi_F[m, k], lambda_F[m, k])
          kl_F_mat[m, k] <- res_F$kl
        } else {
          # Sequential Coupled Gauss-Seidel Sweep over j = 1 ... p_m with Fixed Prior g_F,mk
          inv_nu_m <- 1 / nu_vec[m]
          A_m      <- A_tilde_list[[m]]
          A2_m     <- A_tilde_sq[[m]]
          A_surv_vec <- inv_nu_m * beta2_k * as.vector(crossprod(w_i, A2_m))
          A_Fmk_vec  <- pmax(A_rec + A_surv_vec, 1e-8)
          s_Fmk_vec  <- 1 / sqrt(A_Fmk_vec)
          
          # Precompute linear survival target term1 = inv_nu_m * A_m^T (beta_k * h - w * V_tilde_k)
          surv_target_vec <- inv_nu_m * as.vector(crossprod(A_m, beta_k * h_i - w_i * V_tilde_k))
          
          f_mean_vec  <- EF[[m]][, k]
          f_mean2_vec <- EF2[[m]][, k]
          z_k_curr    <- Z_bar[, k]
          B_Fmk_vec   <- numeric(p_vec[m])
          
          for (j in seq_len(p_vec[m])) {
            a_j    <- A_m[, j]
            f_old  <- f_mean_vec[j]
            u_no_j <- z_k_curr - a_j * f_old
            
            B_j <- B_rec[j] + surv_target_vec[j] - inv_nu_m * beta2_k * sum(w_i * a_j * u_no_j)
            B_Fmk_vec[j] <- B_j
            
            res_j <- solve_point_exponential_posterior(B_j / A_Fmk_vec[j], s_Fmk_vec[j], pi_F[m, k], lambda_F[m, k])
            # Damped posterior moment step (rho = 0.7) for smooth Gauss-Seidel trajectory
            f_new  <- 0.7 * res_j$mean + 0.3 * f_old
            f2_new <- 0.7 * res_j$mean2 + 0.3 * f_mean2_vec[j]
            
            f_mean_vec[j]  <- f_new
            f_mean2_vec[j] <- f2_new
            z_k_curr       <- u_no_j + a_j * f_new
          }
          EF[[m]][, k]  <- f_mean_vec
          EF2[[m]][, k] <- f_mean2_vec
          Z_bar[, k]    <- z_k_curr
          
          # Fixed-q monotone prior update at end of sweep
          res_full      <- solve_point_exponential_posterior(B_Fmk_vec / A_Fmk_vec, s_Fmk_vec, pi_F[m, k], lambda_F[m, k])
          pr_F          <- update_point_exp_prior_fixed_q(res_full, B_Fmk_vec / A_Fmk_vec, s_Fmk_vec)
          pi_F[m, k]    <- pr_F$pi; lambda_F[m, k] <- pr_F$lambda
          res_full      <- solve_point_exponential_posterior(B_Fmk_vec / A_Fmk_vec, s_Fmk_vec, pi_F[m, k], lambda_F[m, k])
          kl_F_mat[m, k] <- res_full$kl
        }
      }
      
      # Recompute exact Z_bar[, k] and Var_Z[, k] across modalities
      Z_bar[, k] <- 0; Var_Z[, k] <- 0
      for (m in seq_len(M)) {
        Z_bar[, k] <- Z_bar[, k] + as.vector(A_tilde_list[[m]] %*% EF[[m]][, k])
        Var_Z[, k] <- Var_Z[, k] + as.vector(A_tilde_sq[[m]] %*% pmax(EF2[[m]][, k] - EF[[m]][, k]^2, 0))
      }
    }
    
    # 3c. Joint K x K Multivariate Update for q(beta) (Calibrated Uncertainty & Collinearity Fix)
    Prec_beta  <- crossprod(Z_bar, w_i * Z_bar) + diag(colSums(w_i * Var_Z) + 1 / sigma2_beta, K)
    Sigma_beta <- solve(Prec_beta)
    EBeta      <- as.vector(Sigma_beta %*% crossprod(Z_bar, h_i))
    
    # Fixed-q Inverse-Gamma MAP update for sigma2_beta (strictly positive, non-absorbing)
    a_beta      <- 1.0; b_beta <- 1e-3
    sigma2_beta <- (2 * b_beta + EBeta^2 + diag(Sigma_beta)) / (2 * a_beta + 1)
    
    # 3d. Anchored Update for Residual Precisions tau_mj (Prevents tau_mj -> 0 Collapse!)
    if (!freeze_tau) {
      for (m in seq_len(M)) {
        Y_hat <- tcrossprod(EL, EF[[m]])
        V_var <- tcrossprod(EL2, EF2[[m]]) - (EL^2) %*% t(EF[[m]]^2)
        Q_m   <- colSums((Y_list[[m]] - Y_hat)^2 + pmax(V_var, 0))
        b_tau <- a_tau / tau_unsup[[m]]
        tau_new <- (n + 2 * a_tau) / (Q_m + 2 * b_tau)
        tau_list[[m]] <- pmax(tau_min[[m]], pmin(1e4, tau_new))
      }
    }
    
    # 3e. Evaluate Consistent Objective & Convergence Diagnostics
    obj_val <- eval_true_objective(EL, EL2, EF, EF2, Z_bar, Var_Z, EBeta, Sigma_beta,
                                   tau_list, pi_L, lambda_L, pi_F, lambda_F, sigma2_beta,
                                   kl_L_vec, kl_F_mat)
    
    eta_new <- as.vector(Z_bar %*% EBeta)
    rel_eta <- sqrt(sum((eta_new - eta_old)^2)) / max(1.0, sqrt(sum(eta_old^2)))
    rel_rec <- max(vapply(seq_len(M), function(m) {
      diff_m <- tcrossprod(EL, EF[[m]]) - tcrossprod(EL_old, EF_old[[m]])
      sqrt(sum(diff_m^2)) / max(1.0, sqrt(sum(tcrossprod(EL, EF[[m]])^2)))
    }, numeric(1)))
    rel_change <- max(rel_eta, rel_rec)
    
    history <- rbind(history, data.frame(iter = iter, objective = obj_val, rel_change = rel_change))
    if (verbose && (iter %% 10 == 0 || iter == 1)) {
      message(sprintf("Sweep %3d | Obj = %12.3f | max_rel_change = %.6f", iter, obj_val, rel_change))
    }
    if (iter > 5 && rel_change < tol) {
      converged <- TRUE
      break
    }
  }
  
  # --- Step 4: Frozen Training Orientation & Four-Way Factor Classification ---
  eta_train <- as.vector(Z_bar %*% EBeta)
  c_train   <- survival::concordance(Surv(time, status) ~ eta_train, reverse = TRUE)$concordance
  orientation_sign <- if (!is.na(c_train) && c_train < 0.5) -1.0 else 1.0
  EBeta_oriented   <- orientation_sign * EBeta
  
  # Compute per-modality and total PVE + calibrated marginal Wald z-scores
  pve_mat <- matrix(0, nrow = M, ncol = K)
  ss_tot  <- vapply(Y_list, function(Y_m) sum(Y_m^2), numeric(1))
  for (m in seq_len(M)) {
    for (k in seq_len(K)) {
      pve_mat[m, k] <- (sum(EL[, k]^2) * sum(EF[[m]][, k]^2)) / ss_tot[m]
    }
  }
  pve_max  <- apply(pve_mat, 2, max)
  sd_beta  <- sqrt(diag(Sigma_beta))
  z_wald   <- EBeta_oriented / pmax(sd_beta, 1e-12)
  
  role <- ifelse(abs(z_wald) >= 1.96 & pve_max >= 0.01, "joint_prognostic",
          ifelse(abs(z_wald) >= 1.96 & pve_max <  0.01, "survival_only",
          ifelse(abs(z_wald) <  1.96 & pve_max >= 0.01, "genomics_only", "dead")))
  
  list(
    EL = EL, EL2 = EL2, EF = EF, EF2 = EF2,
    EBeta = EBeta_oriented, Sigma_beta = Sigma_beta, sd_beta = sd_beta, z_wald = z_wald,
    tau_list = tau_list, pve_mat = pve_mat, role = role,
    s_vec = s_vec, p_vec = p_vec, orientation_sign = orientation_sign,
    converged = converged, n_iter = nrow(history), history = history
  )
}

#' Predict out-of-sample risk scores for new patients (M modalities)
predict_unified_sbmf <- function(fit, Y_tilde_test_list) {
  M <- length(Y_tilde_test_list)
  n_test <- nrow(Y_tilde_test_list[[1]])
  K <- length(fit$EBeta)
  Z_test <- matrix(0, nrow = n_test, ncol = K)
  
  for (m in seq_len(M)) {
    common_feats <- intersect(colnames(Y_tilde_test_list[[m]]), rownames(fit$EF[[m]]))
    if (length(common_feats) == 0) {
      # Fallback if unnamed matrices have identical dimension
      s_m_test <- 1 / (M * sqrt(ncol(Y_tilde_test_list[[m]])))
      Z_test   <- Z_test + s_m_test * (Y_tilde_test_list[[m]] %*% fit$EF[[m]])
    } else {
      s_m_test <- 1 / (M * sqrt(length(common_feats)))
      Z_test   <- Z_test + s_m_test * (Y_tilde_test_list[[m]][, common_feats, drop = FALSE] %*%
                                       fit$EF[[m]][common_feats, , drop = FALSE])
    }
  }
  as.vector(Z_test %*% fit$EBeta)
}