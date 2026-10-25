# Gemini review thread — Supervised EBMF Cox Model (2026-10-17)

- **Source:** shared Gemini conversation, <https://share.gemini.google/gsEy6kv6Z6wZ> (created 2026-10-17, published 2026-10-18).
- **Transcribed:** 2026-10-18 from the rendered share page. Math is the LaTeX source Gemini stored for each equation; tables are converted to Markdown. Gemini's wording, emphasis and equations are unchanged, including any errors in them, except for three display corrections: two absolute values that Gemini stored as norms (`\Vert`; the Laplace density in Turn 1 and the survival-active rule in Turn 5's final table) now use `\vert`, and one garbled expression in Turn 4 (§3, $Z_{ik}$, Eq. 24) is restored to valid LaTeX.
- **Not included:** the attached PDFs themselves, and the interactive prior-comparison widget in Turn 1, which never finished generating in the share.
- **Status:** external advice, not project decisions. Gemini read repository documents as of 2026-10-17; several of its claims about the code and results have not been checked against the repository. Decisions belong in `DECISIONS.md`.

---

## Turn 1

### User

I'm working on building a joint NMF (EBMF) + Cox survival model to do (supervised) survival informed matrix factorization of expression and/or methylation data.

the matrix factorization component has matrices Y = LF + E where LF requires a non-negative matrix factorization with a bayesian formulation and the survival component is cox model.

the priors placed on L & F are exponential and the possible priors placed on B in the cox model are normal, point normal, Laplace, or point laplace

illustrate comparisons between these distributions and provide a recommendation for the priors to use

also provide a brief derivation of the implementation of the updates for each parameter (L/F/B/Tau) (consider the EBMF "Empirical Bayes Matrix Factorization" paper in literature), and illustrate how the model can be setup so that we have supervision that informs the fitting of the "F" parameter. also important to note is that "YFB" should be the linear predictor in the Cox model so that we can project YF to get exact risk scores (instead of having to re-estimate with using LB)

### Gemini

Building a supervised Empirical Bayes Matrix Factorization (EBMF) model where the Cox linear predictor is parameterized as $\eta = YF\beta$ rather than $L\beta$ creates a self-contained projection operator $W = F\beta \in \mathbb{R}^{P \times 1}$. At test time, a new patient's molecular profile $y_{\text{new}} \in \mathbb{R}^{1 \times P}$ yields an exact out-of-sample risk score $\hat{\eta}_{\text{new}} = y_{\text{new}} F \beta$ via simple matrix multiplication, avoiding an iterative inner optimization step to estimate $l_{\text{new}}$.

### 1. Comparison of Prior Distributions for $\beta$ and $L, F$

#### Exponential Priors on $L$ and $F$ (Non-Negative Matrix Factorization)

To enforce non-negativity ($L_{ik} \ge 0, F_{jk} \ge 0$) within the EBMF framework (Wang & Stephens, 2021), we assign independent exponential priors (or point-exponential priors for sparse NMF):

- **Formulation:** $g_L(l_{ik}) = \text{Exp}(\lambda_k^L)$ and $g_F(f_{jk}) = \text{Exp}(\lambda_k^F)$ on $[0, \infty)$, where rates $\lambda_k^L, \lambda_k^F > 0$ are estimated via empirical Bayes.
- **Role:** Constrains the latent space to an additive parts-based representation (e.g., co-expression modules or hypomethylation/hypermethylation signatures) and acts as an $L_1$-like shrinkage on the positive orthant.

#### Candidate Priors for Survival Coefficients ($\beta_k$)

In the Cox model, $\beta \in \mathbb{R}^K$ weights the $K$ projected latent factors $(YF)_{ik}$. Because $K$ is typically moderate ($K \in [5, 50]$), while some factors capture technical batch effects or non-prognostic biological variance, the prior $g_\beta(\beta_k)$ controls factor-level survival selection and shrinkage.

| Prior Family | Density $g(\beta_k)$ | Shrinkage Profile & Posterior Mean Behavior | Posterior Mode (MAP) Sparsity | Variational EBMF Tractability |
|---|---|---|---|---|
| **Normal** | $\mathcal{N}(0, \sigma_\beta^2)$ | **Linear ($L_2$ / Ridge):** Shrinks all coefficients proportionally by $\frac{s_k^2}{s_k^2 + \sigma_\beta^2}$. Over-shrinks strong true survival signals when noise is present. | **Dense:** Never sets $\beta_k = 0$. Every factor contributes to the risk score. | **Closed-form:** Conjugate to the Gaussian working likelihood; trivial moment updates. |
| **Laplace** | $\frac{\lambda_\beta}{2} \exp(-\lambda_\beta \vert{}\beta_k\vert{})$ | **Constant ($L_1$ / Lasso):** Applies constant shrinkage bias $\pm \lambda_\beta s_k^2$ to large signals while strongly shrinking small signals toward zero. | **Sparse (MAP only):** Posterior mode can be zero, but the **posterior mean $\mathbb{E}_q[\beta_k]$ is strictly non-zero**. | **Piecewise Gaussian:** Requires truncated normal moments on $(-\infty, 0)$ and $(0, \infty)$, or a scale-mixture of normals augmentation. |
| **Point-Normal** | $(1 - \pi_\beta)\delta_0 + \pi_\beta \mathcal{N}(0, \sigma_\beta^2)$ | **Spike-and-Slab ($L_0 + L_2$):** Adaptive thresholding via posterior inclusion probability (PIP $\alpha_k$). Weak factors are zeroed out; strong factors receive mild linear shrinkage. | **True Bayesian Sparsity:** $\alpha_k \to 0$ drives $\mathbb{E}_q[\beta_k] \to 0$, cleanly decoupling non-prognostic factors. | **Closed-form:** Standard EBNM analytic updates for $\alpha_k$, posterior mean, and second moment. |
| **Point-Laplace** | $(1 - \pi_\beta)\delta_0 + \pi_\beta \text{Laplace}(0, \lambda_\beta)$ | **Heavy-Tailed Spike-and-Slab:** Combines exact zero-inflation ($\delta_0$) with constant tail bias. Preserves large prognostic effects without the quadratic attenuation of Point-Normal. | **True Bayesian Sparsity + Tail Robustness:** Eliminates noise factors while leaving dominant survival factors virtually unshrunk. | **Closed-form (Piecewise):** Solved analytically via the `ebnm` package using normal-means convolution with Laplace slabs. |

#### Prior Recommendation

- **Primary Recommendation: Point-Normal (or Point-Laplace if $K$ is large).**

    - **Why:** In genomic datasets (expression/methylation), the top variance components in $Y$ are frequently driven by cell-type composition, sex, age, or batch effects that may be completely orthogonal to survival. Because your Cox predictor is $\eta = YF\beta$, any non-zero $\beta_k$ actively pulls $F_{\cdot k}$ toward survival supervision during coordinate ascent. A **spike-and-slab prior (Point-Normal or Point-Laplace)** allows the model to set $\mathbb{E}_q[\beta_k] \approx 0$ (via PIP $\alpha_k \to 0$) for purely unsupervised variance factors, preventing the survival likelihood from distorting biological factors that have no true prognostic value.
    - **Point-Normal vs. Point-Laplace:** Use **Point-Normal** as your default—its second moments $\mathbb{E}_q[\beta_k^2]$ are smoother and less prone to numerical instability during the Taylor-expanded Cox updates. Upgrade to **Point-Laplace** only if you set $K$ very large ($K > 30$) and observe that the quadratic penalty of the Gaussian slab is under-estimating the hazard ratios of 1–2 dominant prognostic factors.

### 2. Joint Model Setup & How Supervision Informs $F$

#### Generative Specification

Let $Y \in \mathbb{R}^{N \times P}$ be the observed molecular matrix ($N$ subjects, $P$ features), decomposed via $K$ latent factors:

1. **Matrix Factorization Likelihood:**

   $$
   Y_{ij} = \sum_{k=1}^K L_{ik} F_{jk} + E_{ij}, \quad E_{ij} \sim \mathcal{N}(0, \tau_j^{-1})
   $$

   where $\tau_j$ is the feature-specific precision (or scalar $\tau$ if homoscedastic), $L \in \mathbb{R}_+^{N \times K}$, and $F \in \mathbb{R}_+^{P \times K}$.
2. **Projected Survival Linear Predictor:** Define the projected patient score matrix $Z = YF \in \mathbb{R}^{N \times K}$, with elements $Z_{ik} = \sum_{j=1}^P Y_{ij} F_{jk}$. The Cox proportional hazards linear predictor for patient $i$ is:

   $$
   \eta_i = (YF\beta)_i = \sum_{k=1}^K \beta_k \left( \sum_{j=1}^P Y_{ij} F_{jk} \right)
   $$

   Given observed survival times $T_i$ and event indicators $\delta_i \in \{0, 1\}$ with risk sets $R_i = \{m : T_m \ge T_i\}$, the Breslow partial log-likelihood is:

   $$
   \ell_{\text{Cox}}(\eta) = \sum_{i=1}^N \delta_i \left[ \eta_i - \log \left( \sum_{m \in R_i} \exp(\eta_m) \right) \right]
   $$

#### How Supervision Flows into $F$ (Dual Likelihood Coupling)

In standard unsupervised EBMF, $F$ only appears in the reconstruction term $\Vert{}Y - LF\Vert{}_F^2$. When we set $\eta = YF\beta$, **the matrix $F$ appears simultaneously in both the Gaussian reconstruction likelihood and the Cox partial likelihood**:

- **Reconstruction Pull:** $F_{\cdot k}$ is pulled toward matching the sample loadings $L_{\cdot k}$ so that $L_{\cdot k} F_{\cdot k}^T \approx R^{(k)}$ (the $k$-th partial residual of $Y$).
- **Survival Pull:** Because $\eta = YF\beta$, each feature weight $F_{jk}$ directly scales feature $j$'s contribution to the hazard ratio via $\beta_k Y_{ij} F_{jk}$. When $\mathbb{E}_q[\beta_k] \neq 0$, the gradient of the Cox likelihood with respect to $F_{jk}$ pushes $F_{jk}$ to up-weight genes/CpGs $j$ whose profiles $Y_{\cdot j}$ covary with survival residuals.
- **Identifiability & Scale Note:** Because $\eta = YF\beta$ depends on the product $F_{\cdot k} \beta_k$ while $Y \approx LF^T$ depends on $L_{\cdot k} F_{\cdot k}^T$, you should **fix the scale of $L_{\cdot k}$ or $F_{\cdot k}$** (e.g., rescale $\mathbb{E}_q[L_{\cdot k}]$ to unit $L_2$ norm at the end of each iteration and absorb the inverse scale into $F_{\cdot k}$ and $\lambda_k$) so that $F$ cannot shrink arbitrarily small while $L$ and $\beta$ inflate to compensate.

### 3. Derivation of Variational EBMF + Cox Coordinate Ascent Updates

Following Wang & Stephens (2021), we posit a mean-field variational posterior that factorizes across all latent variables and parameters:

$$
q(L, F, \beta) = \prod_{k=1}^K \left[ \prod_{i=1}^N q(l_{ik}) \prod_{j=1}^P q(f_{jk}) \right] q(\beta_k)
$$

For any parameter $\theta$, we denote its variational first moment (posterior mean) as $\bar{\theta} = \mathbb{E}_q[\theta]$ and its second moment as $\overline{\theta^2} = \mathbb{E}_q[\theta^2]$, with variance $V_q(\theta) = \overline{\theta^2} - \bar{\theta}^2$.

#### Step 0: Quadratic Approximation to the Cox Partial Likelihood

Because $\ell_{\text{Cox}}(\eta)$ is non-conjugate to Gaussian/Exponential variational updates, we apply a second-order Taylor expansion (Laplace/IRLS working response, as in `glmnet` and Bayesian survival variational inference) around the current expected linear predictor $\hat{\eta}_i = \sum_{k=1}^K (Y \bar{F})_{ik} \bar{\beta}_k$:

- **Gradient ($u_i$) and Negative Diagonal Hessian ($w_i$):**

  $$
  u_i = \frac{\partial \ell_{\text{Cox}}}{\partial \eta_i}\Big\vert{}_{\hat{\eta}} = \delta_i - \sum_{m : i \in R_m} \delta_m \frac{\exp(\hat{\eta}_i)}{\sum_{r \in R_m} \exp(\hat{\eta}_r)}
  $$

  $$
  w_i = -\frac{\partial^2 \ell_{\text{Cox}}}{\partial \eta_i^2}\Big\vert{}_{\hat{\eta}} = \sum_{m : i \in R_m} \delta_m \left[ \frac{\exp(\hat{\eta}_i)}{\sum_{r \in R_m} \exp(\hat{\eta}_r)} - \left( \frac{\exp(\hat{\eta}_i)}{\sum_{r \in R_m} \exp(\hat{\eta}_r)} \right)^2 \right]
  $$
- **Working Survival Response ($z_i$):**

  $$
  z_i = \hat{\eta}_i + \frac{u_i}{w_i}
  $$

  This replaces the Cox log-likelihood in the ELBO with a weighted Gaussian pseudo-likelihood:

  $$
  \ell_{\text{Cox}}(\eta) \approx -\frac{1}{2} \sum_{i=1}^N w_i (z_i - \eta_i)^2 + \text{const}
  $$

#### Step 1: Updating $q(L_{\cdot k})$ (Pure EBNM Exponential Update)

Because the survival predictor $\eta = YF\beta$ does **not** contain $L$, the variational update for $L_{\cdot k}$ is identical to standard EBMF. Define the partial residual excluding factor $k$:

$$
R_{ij}^{(-k)} = Y_{ij} - \sum_{k' \neq k} \bar{l}_{ik'} \bar{f}_{jk'}
$$

Expanding the expected log-likelihood with respect to $l_{ik}$ yields a quadratic form $-\frac{1}{2} A_{ik}^L (l_{ik} - x_{ik}^L)^2$, which reduces to an **Empirical Bayes Normal Means (EBNM)** problem with exponential prior $g_L \in \text{Exp}(\lambda_k^L)$:

- **Precision ($A_{ik}^L$) and Pseudo-Observation ($x_{ik}^L$):**

  $$
  A_{ik}^L = \sum_{j=1}^P \tau_j \overline{f_{jk}^2}, \quad s_{ik}^L = (A_{ik}^L)^{-1/2}
  $$

  $$
  x_{ik}^L = \frac{\sum_{j=1}^P \tau_j R_{ij}^{(-k)} \bar{f}_{jk}}{A_{ik}^L}
  $$
- **Empirical Bayes Hyperparameter & Posterior Moments:** Maximize the marginal log-likelihood to update the prior rate $\hat{\lambda}_k^L = \arg\max_{\lambda > 0} \sum_{i=1}^N \log \int_0^\infty \mathcal{N}(x_{ik}^L \mid l, (s_{ik}^L)^2) \lambda e^{-\lambda l} dl$. Given $\hat{\lambda}_k^L$, the posterior $q(l_{ik}) = \mathcal{N}_+(l_{ik} \mid \mu_{ik}^L, (s_{ik}^L)^2)$ is a **Truncated Normal on $[0, \infty)$** with location $\mu_{ik}^L = x_{ik}^L - \hat{\lambda}_k^L (s_{ik}^L)^2$. Let $\gamma_{ik} = \frac{\mu_{ik}^L}{s_{ik}^L}$ and $r(\gamma) = \frac{\phi(\gamma)}{\Phi(\gamma)}$ be the inverse Mills ratio:

  $$
  \bar{l}_{ik} = \mu_{ik}^L + s_{ik}^L r(\gamma_{ik}), \quad \overline{l_{ik}^2} = (\mu_{ik}^L)^2 + (s_{ik}^L)^2 + s_{ik}^L \mu_{ik}^L r(\gamma_{ik})
  $$

#### Step 2: Updating $q(F_{\cdot k})$ (Supervised EBNM Update)

Because $F_{jk}$ appears in both $Y \approx LF^T$ and $\eta = YF\beta$, we must expand $\mathbb{E}_q[(z_i - \eta_i)^2]$ with respect to $f_{jk}$. Let $\eta_i^{(-k)} = \sum_{k' \neq k} \beta_{k'} \sum_{j=1}^P Y_{ij} f_{jk'}$ and define the expected survival residual excluding factor $k$:

$$
\bar{z}_i^{(-k)} = z_i - \sum_{k' \neq k} \bar{\beta}_{k'} \sum_{j=1}^P Y_{ij} \bar{f}_{jk'}
$$

Furthermore, within factor $k$, when isolating feature $j$, define the factor-$k$ partial projection excluding feature $j$:

$$
\tilde{\eta}_{i,k}^{(-j)} = \sum_{j' \neq j} Y_{ij'} \bar{f}_{j'k}
$$

Collecting all quadratic and linear terms in $f_{jk}$ from both the Gaussian reconstruction ELBO and the Cox quadratic ELBO yields:

- **Combined Precision ($A_{jk}^F$):**

  $$
  A_{jk}^F = \underbrace{\tau_j \sum_{i=1}^N \overline{l_{ik}^2}}_{\text{Expression/Methylation Info}} + \underbrace{\overline{\beta_k^2} \sum_{i=1}^N w_i Y_{ij}^2}_{\text{Cox Survival Info}}
  $$
- **Combined Pseudo-Observation ($x_{jk}^F$):**

  $$
  x_{jk}^F = \frac{1}{A_{jk}^F} \left[ \underbrace{\tau_j \sum_{i=1}^N R_{ij}^{(-k)} \bar{l}_{ik}}_{\text{Reconstruction Alignment}} + \underbrace{\sum_{i=1}^N w_i Y_{ij} \left( \bar{\beta}_k \bar{z}_i^{(-k)} - \overline{\beta_k^2} \tilde{\eta}_{i,k}^{(-j)} \right)}_{\text{Survival Gradient Alignment}} \right]
  $$

> **Computational Vectorization Trick ($O(NP)$):** Computing $\tilde{\eta}_{i,k}^{(-j)}$ naively takes $O(NP^2)$. Instead, maintain the full projected vector $\bar{Z}_{\cdot k} = Y \bar{F}_{\cdot k} \in \mathbb{R}^N$. Then for each feature $j$, $\tilde{\eta}_{i,k}^{(-j)} = \bar{Z}_{ik} - Y_{ij} \bar{f}_{jk}$. Because $f_{jk}$ must be updated coordinate-wise (or via a damped block-jacobi step if $P$ is massive) whenever $\overline{\beta_k^2} > 0$, updating $\bar{Z}_{\cdot k} \leftarrow \bar{Z}_{\cdot k} + Y_{\cdot j}(\bar{f}_{jk}^{\text{new}} - \bar{f}_{jk}^{\text{old}})$ keeps the entire sweep over $j=1\dots P$ strictly $O(NP)$.

- **EBNM Step for $F_{\cdot k}$:** Pass observations $x_{\cdot k}^F$ and standard errors $s_{jk}^F = (A_{jk}^F)^{-1/2}$ into the Exponential EBNM solver to estimate $\hat{\lambda}_k^F$ and compute Truncated Normal moments $\bar{f}_{jk}$ and $\overline{f_{jk}^2}$ exactly as in Step 1.

#### Step 3: Updating $q(\beta_k)$ (Cox Projection EBNM Update)

To update $\beta_k$, we expand $\mathbb{E}_q\left[\sum_{i=1}^N w_i (z_i - \sum_{k=1}^K \beta_k Z_{ik})^2\right]$ where $Z_{ik} = \sum_{j=1}^P Y_{ij} f_{jk}$. Note that because $f_{jk}$ is a random variable under $q$, the second moment $\mathbb{E}_q[Z_{ik}^2]$ includes the posterior variance of $F_{\cdot k}$:

$$
\bar{Z}_{ik} = \sum_{j=1}^P Y_{ij} \bar{f}_{jk}, \quad \overline{Z_{ik}^2} = \bar{Z}_{ik}^2 + \sum_{j=1}^P Y_{ij}^2 V_q(f_{jk})
$$

Isolating $\beta_k$ yields a scalar (or vector across $k$) EBNM problem with precision $A_k^\beta$ and pseudo-observation $x_k^\beta$:

- **Precision ($A_k^\beta$) and Pseudo-Observation ($x_k^\beta$):**

  $$
  A_k^\beta = \sum_{i=1}^N w_i \overline{Z_{ik}^2}, \quad s_k^\beta = (A_k^\beta)^{-1/2}
  $$

  $$
  x_k^\beta = \frac{\sum_{i=1}^N w_i \bar{Z}_{ik} \bar{z}_i^{(-k)}}{A_k^\beta}
  $$
- **Posterior Moments under Your Chosen Prior $g_\beta$:**

    - **Under Point-Normal Prior $g_\beta \sim (1-\pi_\beta)\delta_0 + \pi_\beta \mathcal{N}(0, \sigma_\beta^2)$:** Define the posterior variance given inclusion $v_k = \left( A_k^\beta + \sigma_\beta^{-2} \right)^{-1}$ and posterior mean given inclusion $\mu_k = v_k A_k^\beta x_k^\beta$. The log Bayes Factor for factor $k$ is:

      $$
      \text{BF}_k = \sqrt{\frac{v_k}{\sigma_\beta^2}} \exp\left( \frac{\mu_k^2}{2 v_k} \right)
      $$

      The Posterior Inclusion Probability (PIP) $\alpha_k$ and resulting moments are:

      $$
      \alpha_k = \frac{\pi_\beta \text{BF}_k}{(1 - \pi_\beta) + \pi_\beta \text{BF}_k}, \quad \bar{\beta}_k = \alpha_k \mu_k, \quad \overline{\beta_k^2} = \alpha_k (\mu_k^2 + v_k)
      $$
    - **Under Point-Laplace Prior:** Pass $(x_{1:K}^\beta, s_{1:K}^\beta)$ directly to `ebnm::ebnm_point_laplace()`, which computes the exact piecewise-Gaussian posterior moments $\bar{\beta}_k$ and $\overline{\beta_k^2}$.

#### Step 4: Updating Residual Precision $\tau_j$

The noise precision $\tau_j$ only appears in the Gaussian matrix factorization likelihood. Using the standard EBMF second-moment identity to account for posterior uncertainty in $L$ and $F$:

$$
\bar{R}_{ij}^2 = \mathbb{E}_q\left[ \left( Y_{ij} - \sum_{k=1}^K l_{ik} f_{jk} \right)^2 \right] = \left( Y_{ij} - \sum_{k=1}^K \bar{l}_{ik} \bar{f}_{jk} \right)^2 + \sum_{k=1}^K \left( \overline{l_{ik}^2} \, \overline{f_{jk}^2} - \bar{l}_{ik}^2 \bar{f}_{jk}^2 \right)
$$

The closed-form maximum likelihood / variational update for column-specific precision $\tau_j$ is:

$$
\tau_j = \frac{N}{\sum_{i=1}^N \bar{R}_{ij}^2}
$$

*(If using a Gamma prior $\text{Gamma}(a_0, b_0)$ on $\tau_j$, update as $\tau_j = \frac{a_0 + N/2}{b_0 + \frac{1}{2}\sum_{i=1}^N \bar{R}_{ij}^2}$.)*

### 4. Algorithmic Summary & Practical Implementation Tips

1. **Initialization:** Run standard unsupervised NMF (e.g., via `NNLM::nnmf` or `flashier` with `ebnm_point_exponential`) for 20–30 iterations on $Y$ to initialize $\bar{L}, \bar{F}, \tau_j$. Initialize $\bar{\beta}$ via a ridge or lasso Cox regression (`glmnet`) of $(T, \delta)$ on $Y\bar{F}$.
2. **Relative Scaling Weight ($\gamma$):** In high-dimensional genomics ($P \gg N$), the Gaussian likelihood $\ell_{\text{Gauss}}(Y)$ sums over $N \times P$ terms, whereas $\ell_{\text{Cox}}(\eta)$ sums over $N$ individuals. Although $x_{jk}^F$ naturally balances feature $j$'s column sum (over $N$ samples) against its survival contribution (over $N$ samples), the variance of $Y_{ij}$ affects the relative pull. Standardizing columns of $Y$ (or scaling $Y$ by median non-zero expression for NMF) ensures $w_i Y_{ij}^2$ and $\tau_j \overline{l_{ik}^2}$ operate on commensurate scales.
3. **Out-of-Sample Risk Prediction:** Once converged, compute the composite feature weight vector $\hat{w} = \bar{F}\bar{\beta} \in \mathbb{R}^P$. For any validation cohort $Y_{\text{test}} \in \mathbb{R}^{M \times P}$, the exact linear predictor is $\hat{\eta}_{\text{test}} = Y_{\text{test}} \hat{w}$, and factor-specific risk contributions are $Z_{\text{test}, k} \bar{\beta}_k = (Y_{\text{test}} \bar{F}_{\cdot k}) \bar{\beta}_k$.

### Interactive Exploration: Comparing EBNM Shrinkage & Prior Profiles

Use the interactive tool below to inspect how the four candidate priors for $\beta$ (Normal, Laplace, Point-Normal, and Point-Laplace) transform the pseudo-observations $x_k^\beta$ into posterior expectations $\mathbb{E}_q[\beta_k \mid x_k^\beta, s_k^\beta]$, and how the spike probability ($1 - \pi_\beta$) suppresses noise factors while preserving prognostic signals.

*[Interactive widget not rendered in the share — see note at top.]*

Want R or Python code implementing this CAVI loop with flashier/ebnm?

---

## Turn 2

### User

*Attachments:* `paper.pdf` (DeSurv, Young et al. PNAS 2026), `si_appendix.pdf` (DeSurv SI)

some other notes regarding software and logistical issues:

we use the "flashier" package from the EBMF paper
I have alot of derivations/code already written that I could potentially share so you can understand what I've done, what works, and where issues are showing up
the project is essentially a bayesian implementation with EBMF+Cox of this project/R package that was completed with a frequentist approach: https://github.com/rashidlab/DeSurv-paper/tree/public

additional questions to consider and address in plain easy to understand language so that I am confident I understand the model that I've built. For all of these items, please give a thorough explanation and then wrap up by giving a brief summary where each question is answered in a brief line/bullet point (a couple if necessary):

in a bayesian framework how should we select "K", the number of factors? my current setup involved doing a grid search for "K_init" (the number of initial factors to fit, then selected "K_init" based on best cross-validated C-index on external data. Then we obtained a "K_eff" or the actual retained factors (that explain the survival signal & reconstruction) based on factors that either explained variance above a threshold in the reconstruction and/or factors with a survival Beta above a specified threshold. is this the best way? is there a better way to select K? can we select K via a prior?
I've had issues with involving supervision in the fitting of "F", where it led to convergence issues or the Betas collapsing to 0. after some feature preprocessing I eventually avoided this issue, but I found that including supervision of the fitting of F didn't lead to a performance boost in the cross-validated c-index. however, a potential aim of the project is to identify low-variance prognostic factors (which requires supervision instead of just finding the factors that explain the variance in the data matrix) and to outperform a "2-step" method where we first to NMF/EBMF then "plug" those factors into a cox model (essentially, we want to show that our joint bayesian model performs better at the prediction task & identifies better explanatory factors than doing a 2-step MF then cox method)
if we are interested in identifying a subset of factors that explain the data reconstruction and survival/prognostic signal (different factors can be focused on either task), what prior(s) make most sense to apply to L/F/B/Tau?
my initial attempts at this model involved fitting a "single modality" model where we train the model on 2 pooled cohorts of expression data and then do validation by predicting on 5 cohorts of expression data. there are two extensions I'm interested in applying for the sake of a "value add" of this method over other recent existing work (and we place more emphasis on the 2nd of these extensions)
(1) extend the method to "multiple cohorts" where we specify cohort membership in training data (instead of simply pooling two expression cohorts together). how would we do this and is it feasible? in initial attempts, we added a "dummy indicator" to the loadings matrix where "1" indicated that the sample was included in a cohort and "0" indicated that the sample isn't in the cohort (one column needed to distinguish a pair of cohorts). From my recollection, this attempt didn't lead to any increase in prediction performance in the survival task. If anything, we encountered issues due the need to dedicate an additional degree of freedom to the cohort membership indicator. How should we approach this extension task? does the issue we faced make sense? should we be able to expect an increase in prediction performance by specifying cohort membership vs. pooled data?
(2) the second extension that we are currently working on is extending the model to accommodate multiple data modalities (expression - genes & methylation - CpGs). we initially setup this model by allowing for a separate prior on "F" for each modality. We then plan to evaluate this model by training on one cohort with matched samples of expression & methylation data and then validating the model an another cohort with matched samples of expression & methylation data. In brief that means: train on Cohort A w/ matched expression & methylation samples THEN validate on Cohort B w/ matched expression & methylation samples. At present, I've had convergence issues with the model and also had issues where the fitted model for any specified initial "K" value just retains that K value (and doesn't prune off any uninformative factors - but essentially we've found that maybe every factor allowed has a large enough "contribution" that it is retained). How would you go about approaching this implementation? are the issues we've faced valid & something we can address? I wonder if the K selection procedure could also impact the efficacy of this method? Also, initial implementations of this extension were done with only 10,000 CpGs out of a set of about 300k available CpGs (would like to increase this if it helps), and I found that in the retained factor there was 1 or 2 factors that were related to the methylation data but we didn't see a boost in survival prediction performance
another suggestion for the project was to increase the number of features considered from the current set (about 2000 genes & 10,000 CpGs) to utilize the full set (or at least consider how the model performs across the full range of available features (where we have about 18000 total genes & 300k total CpGs). However, advice was also given that suggests it would be better to settle on a model that uses fewer features if it performs as well in validation as a model its more features & if the factors selected present similar biological signal/interpretation
Additionally, we did initially preprocessing on the gene set by considering the intersection of the top 3000 genes in terms of both variance and mean expression (and I think ended up with 2064 genes). then we z-standardized the expression values (they were originally log transformed as they are RNA-seq data) which helped with some of the original fitting issues (like Betas collapsing to 0). However, it was suggested that I could compare the outcomes of (1) the one step unsupervised preprocessing considering variance & mean expression vs. (2) a two step preprocessing (feature selection) where first filter in an unsupervised manner by variance and mean expression then by a supervised filtering via association to the supervised survival outcome. what do you think of this? essentially this is just filtering down to the features that explain the prognostic signal earlier on right?

I'm curious what your thoughts/ideas/suggestions/concerns/recommendations are for all of the items I discussed above. I've been working on implementation for awhile with Claude & Codex, but there are some larger big-picture problems that require consideration and discussion as well for the sake of strong implementation, personal understanding, and building a project that is targeted in the goals it seeks to accomplish instead of just becoming a mess of every attempt possible at trying to make something happen

In short, the ideas I shared are:

issues with supervision of the matrix factorization parameters ("turning off" supervision essentially leads to a 2-step MF then Cox implementation?)
method for how to select "K", the number of latent factors that explain the data reconstruction and prognostic survival signal
extensions to "multiple cohorts" and "multiple modalities"
considering more/less features for fitting the models
considering adding supervision (with survival outcome) to initial feature selection rather than filtering features based on unsupervised mean/variance

Other thoughts:

while I've implemented this method with Claude/Codex and have substantial results with simulated & real data, I'm also concerned with the amount of bloat thats built up in the repository for this project
I'm considering trying a "fresh start" to build up the project again from the lessons we've learned so far for the sake of eliminating technical debt
the overall goal is to build a model that accommodates multiple cohorts & multiple modalities where the input is expression and/or methylation data with survival outcomes and then the model is able to identify a set of informative latent factors that explain the survival signal as well as the data reconstruction and the factors can be analyzed to understand their relevant biological character (ie. adverse vs. protective for survival, tissue types, etc.)
presumably, we want the method to perform better than simply doing unsupervised NMF then plugging those factors into a cox model for the "survival" component
we also want the model to use bayesian implementation so that it is flexible and requires a minimal amount of tuning/parameter specification.
presumably we could also build results similar to or more interesting than done in the "DeSurv" method (paper attached w/ supplementary info) like a comparison of factors identified by our model vs EBMF, prediction performance vs. 2-step EBMF then Cox, etc.

### Gemini

Having the full `DeSurv` manuscript and SI Appendix provides the exact context needed to diagnose every single issue you have encountered. Translating a frequentist, penalized optimization method (`DeSurv`) into an Empirical Bayes Variational Inference framework (`flashier` / EBMF) introduces subtle mathematical differences in how likelihoods scale, how sparsity operates, and how non-negativity interacts with preprocessing.

Below is a comprehensive breakdown of each question in plain, intuitive language, connecting what happened in `DeSurv` to what is happening inside your Bayesian `flashier` model, followed by a concise bulleted summary at the end.

#### 1. Selecting the Number of Factors ($K$) in a Bayesian EBMF Framework

##### Why your current approach ($K_{\text{init}}$ grid search + post-hoc thresholding for $K_{\text{eff}}$) feels clunky

- In frequentist `DeSurv`, NMF has no built-in mechanism to zero out an entire column of $W$ or $H$. Thus, the authors *had* to run Bayesian Optimization (BO) over a grid of $K$ using cross-validated C-index and the 1-SE rule.
- In `flashier` (Wang & Stephens, 2021), the defining selling point of Empirical Bayes Matrix Factorization is **automatic rank determination via spike-and-slab priors**.
- Why did your factors fail to prune automatically (forcing you to apply post-hoc variance/$\beta$ thresholds)? Because you placed a **pure Exponential prior** (`ebnm_exponential`) on $L$ and $F$. A pure exponential prior has density $\lambda e^{-\lambda x}$ on $[0, \infty)$; it has **no point mass (spike) at zero**. Consequently, the variational posterior mean $\bar{l}_{ik}$ and $\bar{f}_{jk}$ is *strictly positive* for every factor, and `flashier` can never set a factor's loading/factor column to exact zero!

##### How to select $K$ via a prior (The Native EBMF Way)

- **Switch from `ebnm_exponential` to `ebnm_point_exponential`:** The Point-Exponential prior is $(1 - \pi_k)\delta_0 + \pi_k \text{Exp}(\lambda_k)$. When `flashier` updates a factor $k$ that only captures noise, the Empirical Bayes step estimates the active proportion $\hat{\pi}_k \to 0$. Once $\hat{\pi}_k = 0$ on either $L_{\cdot k}$ or $F_{\cdot k}$, the entire factor zeros out automatically during the CAVI loop.
- **Greedy Addition + Backfitting (`flash_greedy` $\to$ `flash_backfit`):** Instead of a grid search over $K_{\text{init}}$, you can either:

    1. **Over-specified $K_{\max}$ with automatic pruning:** Initialize at a generous upper bound (e.g., $K_{\max} = 15$ or $20$) with `ebnm_point_exponential` on $L$ and $F$ and `ebnm_point_normal` (or `ebnm_point_laplace`) on $\beta$. Let the ELBO automatically zero out redundant factors.
    2. **Sequential Greedy + Supervised Backfit:** Add factors one by one until the joint ELBO stops improving (or the newly added factor's $\hat{\pi}_k \to 0$), then run joint supervised backfitting (`flash_backfit`) where $F$ and $\beta$ refine simultaneously.
- **Why is a factor still retained in genomics even with Point-Exponential?** In bulk transcriptomics and methylation data, structured biological noise (e.g., cell cycle, batch, minor cell types, GC content) is *not* i.i.d. Gaussian noise. Because EBMF keeps adding factors as long as a rank-1 component explains more covariance than i.i.d. column noise $\tau_j^{-1}$, an unsupervised likelihood on 10,000+ features will happily retain 20+ biological variance factors. To keep $K_{\text{eff}}$ compact:

    - **Distinguish $K_{\text{total}}$ from $K_{\text{prog}}$ automatically via the $\beta$ prior:** If you use a **Point-Normal or Point-Laplace prior on $\beta$**, each factor $k$ gets a Posterior Inclusion Probability ($\text{PIP}_k = \alpha_k \in [0, 1]$) for survival. You don't need arbitrary ad-hoc thresholds on $\beta$: factors with $\text{PIP}_k > 0.5$ (or $\bar{\beta}_k \neq 0$) are your **prognostic factors** ($K_{\text{prog}}$), while factors with $\bar{L}_{\cdot k} \neq 0$ but $\text{PIP}_k \approx 0$ are your **outcome-neutral background/reconstruction factors** (like the Exocrine factor N2 in `DeSurv`).

#### 2. Supervision of $F$: Convergence Issues, $\beta \to 0$ Collapse, and Beating the 2-Step Baseline

Understanding why your supervised $F$ updates either collapsed $\beta \to 0$, failed to converge, or behaved identically to a 2-step (unsupervised EBMF $\to$ Cox) model is the single most important technical diagnostic in your project. Look closely at how `DeSurv` solved three problems in its Supplementary Methods (SI Sections 1, 2.3, and 16) that must be translated into your Bayesian model:

##### Diagnosis A: The Likelihood Scale Imbalance (Why supervision on $F$ didn't move $F$ away from the 2-step solution)

- Look at the `DeSurv` objective in SI Eq. 1:

  $$
  \mathcal{L}(W, H, \beta) = (1 - \alpha)\frac{1}{2NP}\Vert{}X - WH\Vert{}_F^2 - \alpha \frac{2}{N_{\text{event}}}\ell_{\text{Cox}}(W, \beta)
  $$

  And look at SI Eq. 4: even *after* dividing the NMF loss by $NP$ and the Cox loss by $N_{\text{event}}$, the authors still had to multiply the Cox gradient by a dynamic **gradient-balancing factor** $\zeta^{(t)} = \frac{\Vert{}(1-\alpha)\nabla_W \mathcal{L}_{\text{NMF}}\Vert{}_F^2}{\Vert{}\alpha \nabla_W \mathcal{L}_{\text{Cox}}\Vert{}_F^2}$ because the NMF reconstruction gradient naturally dwarfs the Cox gradient!
- Now look at what happens in a naive Bayesian EBMF + Cox model without tempering:

    - Your expression matrix $Y$ has $N \times P$ data points (e.g., $273 \times 2,000 \approx 546,000$ observations).
    - Your survival outcome $(T, \delta)$ has only $N = 273$ observations ($139$ events).
    - In our Step 2 update for $F_{jk}$ from earlier, the precision is $A_{jk}^F = \tau_j \sum_i \overline{l_{ik}^2} + \overline{\beta_k^2} \sum_i w_i Y_{ij}^2$. Because $Z = YF$ sums over $P = 2,000$ genes, the variance of $Z_{ik}$ is roughly $O(P)$ times larger than a single gene $Y_{ij}$. Consequently, to keep the linear predictor $\eta_i = \sum_k Z_{ik}\beta_k$ on a reasonable hazard ratio scale ($O(1)$), each $\beta_k$ shrinks to $O(1/P)$ or $O(1/\sqrt{P})$!
    - When $\beta_k$ is $10^{-3}$, $\overline{\beta_k^2}$ is $10^{-6}$. Therefore, in the $F_{jk}$ update, the survival term $\overline{\beta_k^2} \sum_i w_i Y_{ij}^2$ is **3 to 6 orders of magnitude smaller** than the reconstruction term $\tau_j \sum_i \overline{l_{ik}^2}$!
    - **The result:** The survival supervision has virtually zero pull on $F$. Your "joint" model is effectively doing 99.99% unsupervised EBMF on $F$, which explains why its cross-validated C-index is identical to the 2-step EBMF $\to$ Cox method!

##### Diagnosis B: Why $\beta$ Collapsed to 0 (and Why Optimization Oscillated)

- **Scale Non-Identifiability:** In $Y \approx LF^T$ and $\eta = YF\beta$, if you don't normalize the columns of $F$ (or $L$) at every single iteration (as `DeSurv` does in Algorithm S2, lines 6–7 via $D = \text{diag}(\Vert{}W_{\cdot k}\Vert{}_2)$), $F$ and $\beta$ can play a tug-of-war. Worse, if $Z_{ik} = (YF)_{ik}$ is not standardized or scaled by $1/P$ (or unit-norm columns of $F$), the raw variance of $Z_{\cdot k}$ is huge, which forces the EBNM prior on $\beta$ to see a tiny signal-to-noise ratio or forces the Taylor expansion weights $w_i$ to explode ($\exp(\hat{\eta}_i) \to \infty$), destabilizing the CAVI loop.
- **The `z-standardization` Band-Aid Warning:** You noted that z-standardizing the genes helped stop $\beta$ from collapsing to 0. Why? Because centering each gene to mean 0 made the projection $Z_{ik} = \sum_j Y_{ij} F_{jk}$ have mean 0 and much smaller variance, stabilizing the Cox linear predictor $\eta_i$! **However, there is a catch:** z-standardizing $Y$ introduces negative numbers into $Y$. If $Y$ has negative values and you constrain $L \ge 0$ and $F \ge 0$, then $LF^T \ge 0$ can *only* reconstruct the positive entries of $Y$ and systematically misses all negative entries!

    - **How `DeSurv` avoided negative numbers while fixing scale:** In `DeSurv` (Materials & Methods Section E / SI Section 7), they kept $Y \ge 0$ by using **within-sample rank transformation** (scaled to $[0, 1]$ or divided by $P$) for the NMF matrix, AND in SI Section 8 they **standardized the projected scores $Z = YF$ using training mean and SD** before multiplying by $\beta$ in the Cox model!
    - **How you should do this cleanly in Bayesian EBMF+Cox:** Keep $Y \ge 0$ (e.g., $\log_2(\text{TPM}+1)$ scaled by column standard deviation without subtracting the mean, OR within-sample rank in $[0,1]$) so that $Y \approx LF^T$ is a valid non-negative factorization. Then, in the Cox linear predictor, **center the columns of $Y$ inside the projection** (i.e., define $\tilde{Y}_{ij} = Y_{ij} - \bar{Y}_{\cdot j}$ so $\eta = \tilde{Y}F\beta$, which is mathematically equivalent to centering $Z_{\cdot k}$ in the Cox model since baseline hazard absorbs constants!) and **enforce column normalization $\Vert{}F_{\cdot k}\Vert{}_2 = 1$** after each $F$ update.

##### Diagnosis C: The Secret Weapon in `DeSurv` — Gene Focusing ($n_{\text{top}} = 270$) and Tempered Likelihood

- Look at **SI Section 16 (Tables S2–S4)** of the attached paper: When `DeSurv` used all 1,970 genes ($n_{\text{top}} = \text{NULL}$) to project $Z = YF$, **it failed adjusted validation in 30 out of 35 $(K, \alpha)$ settings** (and $K=3$ was non-significant at $P = 0.097$)! Only when they truncated each column of $F$ (`W` in their notation) to its top $n_{\text{top}} = 270$ genes did the model validate across cohorts.
- Why? Because in dense NMF, thousands of background genes have small positive weights $F_{jk} > 0$. When you project $Z = YF$, summing across 1,700 noisy background genes drowns out the 250 true marker genes!
- **How your Bayesian model solves this naturally (and beats 2-step NMF):**

    1. **Use `ebnm_point_exponential` on $F$:** Unlike frequentist `DeSurv` which had to do an ad-hoc hard truncation to $n_{\text{top}} = 270$ genes *after* fitting, a **Point-Exponential prior on $F$** forces exact zeros ($F_{jk} = 0$) on non-marker genes *during* variational updates!
    2. **Introduce a Power/Tempering Weight $\gamma$ on the Cox Likelihood:** Just as `DeSurv` uses $\alpha$ to balance $\mathcal{L}_{\text{NMF}}$ and $\mathcal{L}_{\text{Cox}}$, in a joint variational ELBO where $P \gg N$, you should weight the Cox pseudo-likelihood by a balancing scalar $\gamma$ (typically $\gamma \propto \frac{P}{\bar{n}_{\text{event}}}$ or tuned across a small grid of 3–4 values like $\gamma \in \{1, 10, 50, 100\}$). Without $\gamma > 1$, the $N \times P$ Gaussian reconstruction likelihood overwhelms the $N \times 1$ Cox likelihood, preventing supervision from uncovering low-variance prognostic factors.

#### 3. Recommended Priors for $L, F, \beta, \tau$ When Factors Have Dual Roles (Reconstruction vs. Survival)

You want the model to discover a compact set of factors where **some factors explain pure biological variance/reconstruction** (e.g., Exocrine purity, Basal/Classical bulk variance) and **specific factors isolate low-variance prognostic signals** (e.g., iCAF/tumor-stroma coupling). Here is the exact prior configuration that enables this division of labor:

| Parameter | Role in Model | Recommended Prior (`ebnm` family) | Why This Prior Achieves Your Goal |
|---|---|---|---|
| **$L \in \mathbb{R}_+^{N \times K}$** | Sample mixture loadings | **Point-Exponential** (`ebnm_point_exponential`) | Enforces non-negativity ($L_{ik} \ge 0$) so sample loadings represent additive mixture proportions, while the point mass at $0$ allows samples that completely lack a tissue compartment (e.g., zero exocrine contamination) to have exact $L_{ik} = 0$. Also enables automatic factor pruning when $\hat{\pi}_k^L \to 0$. |
| **$F \in \mathbb{R}_+^{P \times K}$** | Feature (Gene / CpG) programs | **Point-Exponential** (`ebnm_point_exponential`) | **Critical for out-of-sample projection ($Z = YF$):** Forces each gene program $F_{\cdot k}$ to be sparse (only true pathway/marker genes are non-zero). Replaces `DeSurv`'s ad-hoc $n_{\text{top}} = 270$ post-hoc truncation with principled Bayesian variable selection during fitting. |
| **$\beta \in \mathbb{R}^K$** | Cox survival weights on $Z = YF$ | **Point-Normal** (`ebnm_point_normal`) or **Point-Laplace** (`ebnm_point_laplace`) | **Separates prognostic factors from reconstruction-only factors:** High-variance non-prognostic factors (like exocrine purity) get $\text{PIP}_k \to 0$ ($\bar{\beta}_k = 0$), so they can freely reconstruct $Y$ without corrupting the Cox risk score $\eta = YF\beta$. Prognostic factors get $\text{PIP}_k \to 1$ ($\bar{\beta}_k \neq 0$) and actively pull $F_{\cdot k}$ toward survival supervision. |
| **$\tau_j \in \mathbb{R}_+^P$** | Feature-specific residual precision | **Column-specific Empirical Bayes / ML** (`var_type = 2` in `flashier`) | Different genes and CpGs have vastly different residual noise variances. Estimating a precision $\tau_j$ per feature prevents noisy genes from dominating the reconstruction loss. |

#### 4. Extension 1: Accommodating Multiple Training Cohorts (TCGA + CPTAC)

##### Why adding a dummy $\{0, 1\}$ cohort indicator column to $L$ failed

- When you fix a column of $L$ to be $1$ for Cohort A and $0$ for Cohort B, its corresponding factor column $F_{\cdot, \text{cohort}}$ simply learns the **mean shift (differential expression) between Cohort A and Cohort B** ($\Delta \mu_j = \mu_{j, \text{Cohort A}} - \mu_{j, \text{Cohort B}}$).
- Why didn't this help out-of-sample survival prediction?

    1. **Within-cohort harmonization already removed global platform shifts:** If your data are rank-transformed within samples (like `DeSurv` did) or column-standardized within each cohort prior to pooling, there is zero mean shift left for the dummy factor to absorb!
    2. **Projection leakage if not zeroed out:** If that dummy factor was accidentally included in $YF\beta$ or interacted with shared factors, it learns cohort-specific technical noise rather than transferable biology.
    3. **Cohort effects in NMF are multiplicative/variance-based, not just additive:** In RNA-seq vs. microarray or TCGA vs. CPTAC, batch differences affect gene-wise noise variance ($\tau_j$) and baseline hazard rates, not just a single rank-1 additive offset on $Y$.

##### How to properly handle multiple cohorts in your Bayesian EBMF+Cox model

Instead of wasting a latent factor column on a dummy indicator, incorporate cohort structure at the **likelihood and precision levels** (which takes zero extra latent degrees of freedom):

1. **Stratified Cox Partial Likelihood (Essential):** Different cohorts have different patient recruitment criteria, follow-up lengths, and baseline survival rates (just as `DeSurv` uses `strata(dataset)` in all pooled Cox analyses). In your Step 0 Cox Taylor expansion, compute the risk sets $R_i$, gradients $u_i$, and Hessian weights $w_i$ **separately within each training cohort $c \in \{1, \dots, C\}$**:

   $$
   \ell_{\text{Cox}}^{\text{strat}}(\eta) = \sum_{c=1}^C \sum_{i \in \text{Cohort } c} \delta_i \left[ \eta_i - \log \left( \sum_{m \in R_i^{(c)}} \exp(\eta_m) \right) \right]
   $$

   This prevents survival supervision on $F$ from being confounded by survival rate differences between TCGA and CPTAC!
2. **Cohort-Specific Residual Precisions ($\tau_{jc}$):** A gene $j$ might be cleanly measured in TCGA but noisy in CPTAC. Allowing precision $\tau_{jc}$ to be estimated per gene *per cohort* ($\tau_{jc} = N_c / \sum_{i \in c} \bar{R}_{ij}^2$) automatically down-weights noisy features in a specific cohort during the $L$ and $F$ updates without adding any extra columns to $L$ or $F$.
3. **Should you expect a massive C-index boost over pooled data?** Be realistic here: if TCGA and CPTAC are both RNA-seq PDAC cohorts that have already been harmonized, specifying cohort membership primarily acts as a **guardrail against batch confounding** rather than a magic booster of C-index. That is why your intuition to **prioritize Extension 2 (Multi-Modality)** as the centerpiece methodological contribution is spot on.

#### 5. Extension 2: Multi-Modal Integration (Expression + Methylation)

Integrating matched gene expression ($Y^{(1)} \in \mathbb{R}_+^{N \times P_1}$) and DNA methylation ($Y^{(2)} \in \mathbb{R}_+^{N \times P_2}$) with survival supervision is a fantastic, high-impact extension over `DeSurv` (the `DeSurv` discussion explicitly notes multi-omics integration as an important open next step!). Let's diagnose the exact three issues you hit: (1) convergence issues, (2) retaining all $K$ factors without pruning, and (3) methylation only showing up in 1–2 factors without improving C-index.

##### Issue A: Why $Y^{(2)}$ (10,000 CpGs) drowned out $Y^{(1)}$ (2,000 Genes) and caused convergence issues

- In a multi-modal shared-$L$ model:

  $$
  Y^{(1)} = L (F^{(1)})^T + E^{(1)}, \quad Y^{(2)} = L (F^{(2)})^T + E^{(2)}
  $$

  Notice what happens when you update the shared patient loadings $L_{ik}$:

  $$
  A_{ik}^L = \underbrace{\sum_{j=1}^{P_1} \tau_j^{(1)} \overline{(f_{jk}^{(1)})^2}}_{\text{2,000 Genes}} + \underbrace{\sum_{m=1}^{P_2} \tau_m^{(2)} \overline{(f_{mk}^{(2)})^2}}_{\text{10,000 (or 300k) CpGs}}
  $$
- Because $P_2$ (10,000 CpGs) is $5\times$ larger than $P_1$ (2,000 genes)—and would be $150\times$ larger if you used 300k CpGs!—the sum over $P_2$ completely dominates $A_{ik}^L$ and $x_{ik}^L$.
- **The Fix (Modality Likelihood Balancing):** Just as `MOFA` (Multi-Omics Factor Analysis) and `DeSurv` normalize Frobenius losses by dimension $(N P_m)$, you **must weight each modality's log-likelihood by $w_{\text{mod}}^{(m)} = \frac{\bar{P}}{P_m}$** (where $\bar{P}$ is a reference dimension like $2,000$, or $1/P_m$). Practically, this simply scales the effective precision $\tilde{\tau}_j^{(m)} = \left(\frac{\bar{P}}{P_m}\right) \tau_j^{(m)}$ so that 2,000 genes and 10,000 CpGs have equal voice (50% / 50%) in determining the shared patient loadings $L$!

##### Issue B: How to structure the Cox Linear Predictor $\eta$ across Modalities

- In your multi-modal model, you have two feature factor matrices: $F^{(1)} \in \mathbb{R}_+^{P_1 \times K}$ (genes) and $F^{(2)} \in \mathbb{R}_+^{P_2 \times K}$ (CpGs). How is $\eta$ formed?

    - **Option 1 (Concatenated Projection per Factor $k$):** Define the $k$-th factor's multi-modal projected score as the weighted average of its normalized projections:

      $$
      Z_{ik} = \frac{1}{2} \left( \tilde{Y}^{(1)} F_{\cdot k}^{(1)} + \tilde{Y}^{(2)} F_{\cdot k}^{(2)} \right), \quad \eta_i = \sum_{k=1}^K Z_{ik} \beta_k
      $$

      *(Or allow separate survival weights $\beta_k^{(1)}$ and $\beta_k^{(2)}$, or let the relative norms of $F_{\cdot k}^{(1)}$ and $F_{\cdot k}^{(2)}$ determine whether factor $k$ is driven by expression, methylation, or both!)*
    - **Column Normalization in Multi-Modal EBMF:** To prevent scale oscillation between $L$, $F^{(1)}$, $F^{(2)}$, and $\beta$, **normalize the rows/columns of $L$ instead of $F$** (set $\Vert{}L_{\cdot k}\Vert{}_2 = \sqrt{N}$ after Step 1, and absorb the scale into $F_{\cdot k}^{(1)}$ and $F_{\cdot k}^{(2)}$). Why? Because if you fix the scale of $L_{\cdot k}$, then $F_{\cdot k}^{(1)}$ and $F_{\cdot k}^{(2)}$ are naturally allowed to have *different* norms! If factor $k$ is a pure methylation factor, $\Vert{}F_{\cdot k}^{(1)}\Vert{}_2 \to 0$ (genes turn off) while $\Vert{}F_{\cdot k}^{(2)}\Vert{}_2 > 0$ (CpGs stay active)!

##### Issue C: Why every factor was retained ($K_{\text{eff}} = K_{\text{init}}$) and why Methylation didn't boost C-index

1. **Why no factors pruned:** As explained in Question 1, using a pure Exponential prior (`ebnm_exponential`) instead of `ebnm_point_exponential` prevents any factor from ever reaching exact zero. Furthermore, in 10,000 CpGs + 2,000 genes, there are easily 20+ orthogonal axes of real biological/epigenetic covariation (cell type methylation, age-related CpG drift, island vs. open-sea methylation). Every factor you initialize *will* find a real covariance pattern in 10,000 CpGs that beats random Gaussian noise!

     - **Solution:** Use `ebnm_point_exponential` on $F^{(1)}$ and $F^{(2)}$ AND use `ebnm_point_normal` on $\beta$. Even if all $K = 8$ factors are retained for *reconstruction* (because there are 8 real biological patterns across expression and methylation), the **Point-Normal prior on $\beta$ will prune the survival coefficients $\bar{\beta}_k \to 0$** on the 5–6 factors that don't predict survival!
2. **Why methylation only appeared in 1–2 factors and didn't boost C-index yet:**

     - DNA methylation beta-values lie in $[0, 1]$ and are strongly bimodal (most CpGs are either permanently unmethylated $\approx 0.05$ or permanently methylated $\approx 0.95$). If you don't filter CpGs properly (e.g., mapping promoter/enhancer CpGs or top variable CpGs) or if the methylation scale is much smaller than log-expression, a single factor will just soak up global methylation level while the rest focus on expression.
     - More importantly, in many cancers (including PDAC), transcriptomic subtypes (Classical/Basal/iCAF) already capture a large share of the prognostic variance. **Do not judge multi-modal success solely by whether C-index jumps from 0.65 to 0.72!** Even if the multi-modal C-index is similar (or slightly better, e.g., $+0.01$ to $+0.02$), showing that **joint EBMF+Cox links specific epigenetic CpG hyper/hypo-methylation modules directly to transcriptional tumor/stroma programs within the same prognostic factor $k$** is a massive biological value-add that single-modality `DeSurv` cannot do!

#### 6. Feature Count (2,000 Genes & 10,000 CpGs vs. 18,000 Genes & 300k CpGs)

Your advisor/collaborator's advice to **settle on a focused feature set (e.g., ~2,000–3,000 genes and ~5,000–10,000 CpGs) if validation performance and biological interpretation match the larger set** is 100% the right call, for three concrete reasons:

1. **Signal-to-Noise in Projection ($Z = YF$):** Remember that your linear predictor projects all features via $Z = YF$. If $P = 300,000$ CpGs, even tiny non-zero posterior noise on 290,000 irrelevant CpGs will accumulate massive variance in the sum $\sum_{j=1}^{300,000} Y_{ij} F_{jk}$, diluting the true prognostic signal (the exact phenomenon proven in `DeSurv` SI Table S2 where 1,970 genes failed without focusing to 270 genes).
2. **Reconstruction Swamping:** With 300,000 CpGs, roughly 250,000 are housekeeping/constitutive sites or sex-chromosome/age-related drift. Including them forces the NMF reconstruction term to waste factors modeling non-prognostic epigenetic background.
3. **Computational Speed & Iteration Agility:** Fitting CAVI across 2,000 genes + 5,000–10,000 CpGs runs in seconds to minutes, allowing thorough cross-validation and simulation benchmarking.

- **How to frame this in your paper:** Keep ~2,000–3,000 genes and ~5,000–10,000 CpGs (ideally promoter/enhancer-mapped or most variable CpGs) as your primary architecture, and include **one small sensitivity figure/table in the Supplement** showing that expanding from 2k to 5k/10k genes (or 10k to 25k/50k CpGs) yields highly correlated factor loadings ($r > 0.85$) and no meaningful gain in validation C-index, while increasing runtime $10\times$. Reviewers love that justification.

#### 7. 1-Step Unsupervised vs. 2-Step Supervised Feature Pre-Filtering

You asked about comparing:

- **(1) 1-Step Unsupervised Pre-filtering:** Filter to the top ~2,000–3,000 features by mean & variance (exact `DeSurv` recipe: rank by mean + rank by variance, take intersection across training cohorts).
- **(2) 2-Step Supervised Pre-filtering:** First filter by unsupervised mean & variance, then filter further by univariate Cox association with survival (essentially the Bair & Tibshirani 2004 / 2006 "Supervised PCA" pre-screening idea, cited as Refs 2 & 16 in `DeSurv`).

##### My strong recommendation: Keep (1) Unsupervised Mean/Variance as your primary pipeline, and only test (2) as an ablation/comparison!

Here is why:

1. **It undermines your core methodological thesis if used as the default:** The entire thesis of `DeSurv` (and your Bayesian EBMF+Cox model) is that **the joint factorization model itself discovers low-variance prognostic programs amidst high-variance background genes without needing pre-filtering by survival**! If you pre-filter the input matrix $Y$ using survival p-values before running your model, a skeptical reviewer will immediately ask: *"If you already threw away all the outcome-neutral genes using univariate Cox p-values in preprocessing, why do you need a joint supervised NMF model? Wouldn't standard unsupervised EBMF on those pre-selected survival genes do the exact same thing?"*
2. **Univariate Cox screening misses cooperative/masked signals:** In PDAC, a gene expressed in both tumor and stroma might have a weak marginal univariate Cox p-value because its expression is confounded by tissue purity, yet when deconvolved into an iCAF factor by NMF, it becomes strongly prognostic. Pre-filtering by marginal survival association risks throwing away genes needed to cleanly deconvolve the background compartments (like Exocrine or Basal/Classical) that allow the prognostic factor to separate!
3. **Where 2-step pre-filtering *is* useful:** You can include it as a **benchmark comparison** (or use a very gentle semi-supervised screen *only* for reducing 300k CpGs down to 10k CpGs if unsupervised variance alone picks up too many noisy intergenic probes—though even for CpGs, filtering to gene-promoter/enhancer islands with high variance is cleaner and strictly unsupervised).

#### 8. Recommendation on a "Fresh Start" to Eliminate Repository Bloat

When working iteratively with AI coding assistants over weeks of debugging, repositories almost always accumulate "scar tissue": 15 flags for different normalizations, dead fallback branches, ad-hoc thresholds, and conflicting helper scripts.

**Do a clean, modular "v2" rebuild—but don't throw away your old repo.** Keep the old repository as an archive/reference, and create a clean, minimal R package structure with just **four core files**:

1. `R/preprocess.R`: Clean, reproducible functions for non-negative feature filtering and scaling (keeping $Y \ge 0$ for NMF while centering $\tilde{Y}$ for the Cox projection).
2. `R/cavi_updates.R`: The core CAVI loop (`update_L`, `update_F_supervised`, `update_beta_cox`, `update_tau`) with explicit likelihood tempering weight $\gamma$ and modality balancing weights $w_{\text{mod}}^{(m)}$.
3. `R/predict_and_eval.R`: Out-of-sample projection $\hat{\eta}_{\text{new}} = \tilde{Y}_{\text{new}} \bar{F} \bar{\beta}$, C-index calculation, and per-factor variance ($R^2$) vs. survival ($\Delta \ell_{\text{Cox}}$ and $\text{PIP}_k$) metrics (matching Fig. 3C of `DeSurv`).
4. `documentation.md`: A living mathematical and architectural specification documenting every equation and design choice.

#### Summary of Key Takeaways

- **Selecting $K$ in a Bayesian Framework:**

    - Replace pure Exponential (`ebnm_exponential`) with **Point-Exponential (`ebnm_point_exponential`)** on $L$ and $F$ and **Point-Normal (`ebnm_point_normal`)** on $\beta$ so the model has a spike at zero ($0$) to prune dead factors and zero out noisy genes automatically.
    - Let the **Posterior Inclusion Probability ($\text{PIP}_k$) of $\beta_k$** naturally separate your retained factors into **prognostic factors** ($\text{PIP}_k > 0.5$) and **reconstruction-only background factors** ($\text{PIP}_k \approx 0$), eliminating arbitrary post-hoc cutoffs.
- **Why Supervision of $F$ Struggled & How to Beat the 2-Step Baseline:**

    - The Gaussian likelihood sums over $N \times P$ terms while the Cox likelihood sums over $N$ terms, making the survival gradient $10^3\text{–}10^6\times$ smaller than the reconstruction gradient in the $F$ update.
    - Fix this by: (1) adding a **likelihood tempering weight $\gamma$** (analogous to `DeSurv`'s $\alpha$ and $\zeta^{(t)}$) to scale up the Cox information in the $F$ update, (2) enforcing **unit-norm scaling** on factor columns at every iteration, and (3) using **Point-Exponential priors on $F$** so background genes are zeroed out before projection $Z = YF$ (performing `DeSurv`'s $n_{\text{top}} = 270$ gene focusing automatically inside the EBNM step).
    - Avoid z-standardizing $Y$ directly (which creates negative values that violate NMF $LF^T \ge 0$); keep $Y \ge 0$ for the NMF reconstruction and only center $Y$ inside the Cox linear predictor ($\eta = (Y - \bar{Y})F\beta$).
- **Best Priors for $L, F, \beta, \tau$:**

    - **$L$:** Point-Exponential (non-negative sparse sample mixtures + factor pruning).
    - **$F$:** Point-Exponential (non-negative sparse gene/CpG programs $\to$ clean out-of-sample projection $YF$).
    - **$\beta$:** Point-Normal (or Point-Laplace) (zeroes out $\beta_k$ on non-prognostic variance factors so they don't corrupt survival prediction).
    - **$\tau_j$:** Feature-specific Empirical Bayes/ML precision (handles heteroscedastic noise across genes/CpGs).
- **Extension 1 (Multiple Training Cohorts):**

    - Drop the dummy $\{0,1\}$ column in $L$ (which only captures an additive mean shift and wastes a factor).
    - Instead, use a **cohort-stratified Cox partial likelihood** (separate risk sets $R_i^{(c)}$ per cohort) and **cohort-specific residual precisions $\tau_{jc}$**. Expect this to prevent batch confounding rather than dramatically inflating C-index.
- **Extension 2 (Multiple Modalities: Expression + Methylation):**

    - Weight each modality's reconstruction precision by $w_{\text{mod}}^{(m)} = \bar{P} / P_m$ so 10,000 CpGs do not drown out 2,000 genes when updating shared loadings $L$.
    - Normalize $L_{\cdot k}$ rather than $F_{\cdot k}$ so that modality-specific factor loadings $F_{\cdot k}^{(1)}$ (genes) and $F_{\cdot k}^{(2)}$ (CpGs) can freely take different magnitudes (allowing pure-expression, pure-methylation, or joint factors).
    - Emphasize **biological integration** (linking CpG modules to transcriptional programs within prognostic factors) alongside C-index parity/improvement.
- **Feature Count (Focused Set vs. All 18k Genes / 300k CpGs):**

    - Stick with the **focused feature set (~2,000–3,000 genes and ~5,000–10,000 promoter/enhancer CpGs)** as your main model to prevent projection variance dilution in $Z = YF$, and include a small supplementary sensitivity check across larger feature counts.
- **1-Step Unsupervised vs. 2-Step Supervised Feature Pre-Filtering:**

    - Use **1-step unsupervised filtering (mean + variance)** as your primary pipeline. Pre-filtering genes by survival p-values before running your model undercuts the core claim that your joint factorization discovers prognostic factors on its own without pre-screening.

Would you like to share the core update functions from your current codebase next so we can inspect how the Cox weights and `flashier` EBNM calls are wired and sketch out the clean v2 architecture?

---

## Turn 3

### User

*Attachments:* `writeup.pdf` (`paper/writeup/writeup.pdf` as of 2026-10-17)

oh also - here is a recent writeup (build by codex) that explains the current implementation of my method with method/derivations/results/biology (may not be entirely accurate) hopefully this can help get an understand of the "gap" between how the method should be implemented versus how it currently is implemented!

### Gemini

Sharing this Codex writeup is a game-changer. It acts like an X-ray of your entire codebase: it exposes the exact lines of math and code where your current implementation diverges from the model you *intended* to build, and explains every single anomalous result you have seen.

Here is a plain-language diagnostic of the "gap" between your current implementation and the ideal architecture, followed by a blueprint for your clean v2 rebuild and a concise summary.

#### 1. The 5 "Smoking Guns" in the Writeup (Why Your Code Is Behaving Strangely)

##### Smoking Gun #1: Supervision is literally turned off in your Single-Modality benchmark ($\alpha_F = 0$)—and effectively $10^{-7}$ in your Multimodal model!

- **What the writeup reveals:**

    - In your single-modality model (Section 1.2 & Appendix A.7, Eq. 18), the weight mixing the reconstruction and survival terms in the $F$ update is set to **$\alpha_F = 0$**. That means survival information **never enters $F$ at all**! Your single-modality model is currently doing 100% unsupervised EBMF on $L$ and $F$, and just running a Cox update on $\beta$ alongside it.
    - Even worse, look at **Appendix A.7 (Eq. 18)** for when you *did* try turning $\alpha_F > 0$ earlier:

      $$
      B_{F,jk} = (1 - \alpha_F)B_{F,jk}^{\text{rec}} + \alpha_F \bar{\tilde{\beta}}_k \sum_i w_i y_{ij}(z_i - \bar{\eta}_{i,-k})
      $$

      As noted in the writeup right below Eq. 18, **that formula dropped the coupling term $-\overline{\tilde{\beta}_k^2}\sum_i w_i \frac{y_{ij}}{d_k}\bar{C}_i$ AND updated all $P$ features simultaneously (Jacobi-style) instead of sequentially**! Updating all 2,064 features simultaneously along the survival gradient without accounting for how changing gene $1$ changes the linear predictor $\eta_i$ for gene $2$ causes the linear predictor $\eta$ to overshoot by $2,064\times$ on the next iteration, which immediately causes oscillation or forces $\beta \to 0$!
    - Now look at **Appendix C.1 (page 47)** for the Multimodal model (where the coupling term *was* included via sequential Gauss-Seidel updates in Eq. 30–31):

      > *"A separate read-only audit of the frozen joint fit's sole credible program finds median survival fractions of loading precision of **$2.03 \times 10^{-7}$ for expression** and **$9.00 \times 10^{-10}$ for methylation**; fractions of the absolute numerator terms are $7.27 \times 10^{-5}$ and $4.40 \times 10^{-6}$."*
      >
      > Because the multimodal update (Eq. 30–31) adds the reconstruction and survival terms with **equal weight (1:1)**, the survival precision is **ten million times smaller** ($2 \times 10^{-7}$) than the reconstruction precision! That is why the multimodal simulation found **0.00 (0 out of 27)** of the planted low-variance prognostic programs at lower noise (Table 6)!

##### Smoking Gun #2: You Have Two Completely Different Models for Single-Modality vs. Multimodal

Look at **Table 1 (page 3)**: instead of the multimodal model being a multi-block extension of the single-modality model, Codex built two incompatible architectures:

- **Data & Priors:** Single-modality uses z-standardized data (which has negative numbers!) with **non-negative Point-Exponential priors** on $F$ and no intercept ($Y = LF^T + E$). Multimodal uses unstandardized data with a feature intercept $\mu_m$ and **signed Point-Laplace priors** on $F_m$ ($Y_m = 1\mu_m^T + LF_m^T + E_m$). **Wait:** your multimodal model isn't even doing Non-Negative Matrix Factorization on $F_m$ right now—it's doing semi-NMF ($L \ge 0$, $F_m \in \mathbb{R}$ signed)!
- **Scale Normalization:** Single-modality divides by Euclidean column norm $d_k = \Vert{}F_{\cdot k}\Vert{}_2$. Multimodal rescales by $c_k = \text{sd}(\bar{Z}_{\cdot k})$.
- **Factor Pruning:** Single-modality never deletes factors (just labels them post-hoc). Multimodal uses an ad-hoc post-loop deletion score $S$ (Eq. 36).

##### Smoking Gun #3: Why the Multimodal Model Fails to Converge (500-Sweep Cap) and Never Removes Any Factors on Real Data

Look at **Table 8 (page 20)** and **Appendix B.11 / C.2**:

1. Why it never converges (hits the 500-sweep cap on 7 of 9 ranks):

     - In Algorithm 2 (page 39), step 2b runs $S=5$ inner sweeps where for every feature $j$, it updates $q(f_{mjk})$, and **at the end of each column it refits the prior $g_{mk}$ AND rescales the factor by $c_k = \text{sd}(\bar{Z}_{\cdot k})$**.
     - Look at **Appendix B.9 (page 44)** and **Appendix C.2 (page 47)**: when it rescales $(L_{\cdot k}, F_{mk}, \beta_k)$ by $c_k$, it **forgets to rescale the Normal prior variance $\sigma_k^2$ on $\beta_k$** AND when it refits the prior $g_{mk}$, it **caches the KL divergence against the *old* prior instead of the newly refitted prior**! Because the prior parameters and scale are shifting out of sync with the variational objective every inner sweep, the coordinate updates are not ascending a valid ELBO and enter a small limit cycle (oscillating by ~0.005–0.019 forever instead of settling below $10^{-4}$)!
2. Why it never removes a single factor on real data (Final rank = Starting rank in Table 8):

     - Look at Eq. 36 (page 45): your factor removal check computes $\Delta_k = S(\text{factors} \neq k) - S(\text{all factors})$ **while holding all other $K-1$ factors completely frozen**!
     - In a matrix of $144 \times 13,000$ features ($1.87$ million entries), if you suddenly delete factor $k$ and hold the other $K-1$ factors frozen without letting them adapt, the sum of squared errors $\sum_{m,i,j} \tau_{mj}(\tilde{y}_{mij} - l_i^T f_{mj})^2$ increases massive amounts (because any factor explaining even 0.5% of variance across 13,000 features contributes thousands of log-likelihood units, as shown by the 610,423 reconstruction score on page 20)! Thus $\Delta_k$ is *always* strongly negative when evaluated with frozen remaining factors, and zero factors are ever removed!

##### Smoking Gun #4: Why the 25,000-CpG Run Exploded ($\text{SD}(\hat{\eta}) = 39,086$)

Look at **Section 3.6 (page 22)** and **Appendix C.1 (page 47)**:

- In the multimodal model, methylation is unstandardized $\beta$-values in $[0, 1]$, and $Z_{ik} = \sum_{m,j} \tilde{y}_{mij} f_{mjk}$ is an **unnormalized sum over $P_m$ features**.
- When you increased from 10,000 CpGs to 25,000 CpGs, the raw sum $Z_{ik}$ grew $2.5\times$ larger during the inner feature loop *before* the end-of-factor rescaling step $c_k = \text{sd}(\bar{Z}_{\cdot k})$ could catch it. That caused $\hat{\eta}_i$ to jump to $5,734$ at outer sweep 6, which caused $\exp(\hat{\eta}_i)$ to overflow/underflow in the Cox risk set denominator! Even after fixing the log-sum-exp arithmetic, the risk score SD blew up to $39,086$ because **there is no step-size damping or column-norm constraint during the inner $F$ sweep**, and **no floor on the curvature $w_i$** in the multimodal Cox approximation (Appendix B.2 explicitly notes: *"There is no floor on $w_i$"*, whereas Single-Modality floored $w_i$ at $10^{-6}$!).

##### Smoking Gun #5: The "Double Log-Transformation" Bug in Preprocessing

Look at **Section 1.4 (page 5)** and **Appendix C.3 (page 49)**:

- The writeup discovered that your data loader applied $\log_2(x + 1)$ to TCGA, Dijk, and PACA-AU RNA-seq, and then your benchmark preprocessing script applied $\log_2(x + 1)$ **a second time** ($\log_2(\log_2(\text{TPM}+1)+1)$), while CPTAC was only log-transformed once!
- Fixing that in your clean v2 rebuild is trivial and will immediately make your TCGA and CPTAC training distributions match properly.

#### 2. How to Unify & Fix Everything in a Clean "v2" Rebuild

Instead of maintaining two contradictory codebases (Single-Modality vs. Multimodal) with broken KL bookkeeping and ad-hoc removal scores, **you only need ONE unified model** where Single-Modality is simply the $M = 1$ special case of the $M$-modality model!

Here is how every piece of your unified v2 model should work so that it is mathematically airtight, converges monotonically, and achieves your exact goals:

##### A. Unified Preprocessing & Non-Negativity (Keeping NMF Interpretable!)

1. **Fix the Double-Log Bug:** Ensure all RNA-seq cohorts are $\log_2(\text{TPM}+1)$ (or within-sample rank in $[0,1]$ as in `DeSurv`). For methylation, convert CpG $\beta$-values (or M-values shifted to be $\ge 0$, or $\beta$-values scaled by their within-cohort SD without subtracting the mean) so that **all input blocks $Y_m \ge 0$ are strictly non-negative**.
2. **Why $Y_m \ge 0$ without centering in the reconstruction matters:** If you keep $Y_m \ge 0$ (divide each feature $j$ by its training standard deviation $s_{mj}$ so features have unit variance, *without* subtracting the mean), then you can use **Point-Exponential (`ebnm_point_exponential`) priors on BOTH $L$ and $F_m$**! Every factor loading $F_{mjk} \ge 0$ is then a true non-negative sparse gene/CpG program (just like `DeSurv`), rather than a signed PCA-like factor that is hard to interpret biologically!
3. **Where centering happens (Inside the Cox Projection Only!):** In the Cox model, define the centered feature matrix $\tilde{Y}_{m,ij} = Y_{m,ij} - \bar{Y}_{m,\cdot j}$. The linear predictor is:

   $$
   \eta_i = \sum_{k=1}^K Z_{ik} \beta_k, \quad Z_{ik} = \sum_{m=1}^M w_{\text{mod}, m} \left( \tilde{Y}_m \frac{F_{mk}}{\Vert{}F_{mk}\Vert{}_2} \right)_i
   $$

   Because $\tilde{Y}_m$ is centered, $\bar{Z}_{\cdot k}$ has mean 0 in every cohort, completely eliminating intercept drift in the Cox model while preserving 100% pure NMF ($Y_m \approx L F_m^T$) in the reconstruction!

##### B. Fixing the Likelihood Imbalance (So Supervision Actually Guides $F$ and Finds Low-Variance Factors!)

Why did the multimodal model have a $10^{-7}$ survival fraction in the $F$ update? Because the reconstruction precision sums over $N$ samples with $\tau_{mj} \sim O(1)$, while the survival precision is $\overline{\beta_k^2} \sum_i w_i \tilde{y}_{mij}^2$ unscaled by feature count! To make the survival gradient and reconstruction gradient commensurate (just like `DeSurv` uses $(1-\alpha)\frac{1}{NP}$ vs. $\alpha\frac{2}{N_{\text{event}}}$):

1. **Modality Balancing Weight ($v_m$):** Weight each modality $m$'s reconstruction likelihood by $v_m = \frac{\bar{P}}{P_m}$ (e.g., $\bar{P} = 2,000$). That way, 10,000 CpGs don't drown out 2,000 genes by $5:1$ when updating $L_{ik}$.
2. **Survival Supervision Weight ($\gamma$):** Multiply the Cox working precision $w_i$ by a supervision scaling weight:

   $$
   \gamma = \frac{\alpha}{1 - \alpha} \cdot \frac{\sum_m v_m P_m}{N_{\text{event}}}
   $$

   where $\alpha \in [0, 1)$ is your single intuitive supervision strength knob (identical to `DeSurv`!).

     - When $\alpha = 0 \implies \gamma = 0$: Pure unsupervised EBMF + 2-step Cox projection.
     - When $\alpha \in (0, 1)$ (e.g., $\alpha = 0.35$ as in `DeSurv`): The survival precision in the $F_{mjk}$ update is on the **exact same order of magnitude** ($O(1)$) as the reconstruction precision, rather than $10^{-7}$!
3. **Why planted low-variance programs weren't found when $\bar{\beta}_k = 0$:** In Section 5.3 (Aim 1, page 29), the writeup makes a brilliant observation:

   > *"the survival contribution is multiplied by the current coefficient, so a direction no program already points along receives none"*
   >
   > If you initialize $L$ and $F$ purely from the top $K$ SVD components of $Y$ (which only capture the highest-variance directions!), a low-variance prognostic program isn't in the initial span of $F$, so its initial $\bar{\beta}_k \approx 0$, and therefore its survival gradient $\bar{\beta}_k (\dots)$ is zero!

     - **The Fix (Supervised Warm-Start Initialization for 1–2 Factors):** Initialize the first $K-1$ factors from NMF/SVD of $Y$, and initialize at least one factor's loading vector $F_{\cdot K}^{(0)}$ using the positive/negative univariate Cox score vector $\max(0, \pm Y^T u^{(0)})$ (the initial Cox martingale residuals)! That guarantees at least one factor starts pointed along the survival gradient with $\bar{\beta}_K^{(0)} \neq 0$, so CAVI immediately locks onto low-variance prognostic programs!

##### C. Fixing Convergence & Scale Canonicalization (No More 500-Sweep Caps!)

To eliminate the oscillation documented in Appendix B.9 and C.2:

1. **Use exact `ebnm` objects for KL (Never cache stale priors!):** When you call `res <- ebnm(x, s, prior_family = "point_exponential")`, `ebnm` returns `res$log_likelihood` (the exact marginal log-likelihood $\sum_h \log p_{\hat{g}}(x_h)$ under the **newly fitted** prior $\hat{g}$!). By the EBMF identity (Appendix A.3, Eq. 11):

   $$
   -\text{KL}(q \,\Vert{}\, \hat{g}) = \log p_{\hat{g}}(x) - \mathbb{E}_q[\log \mathcal{N}(x; \theta, s^2)]
   $$

   If you compute the KL immediately from the `ebnm` return object at the moment of the update, there is **zero mismatch** between the posterior and the refitted prior!
2. **Fix the Scale of $L_{\cdot k}$ (Unit Root-Mean-Square $\frac{1}{N}\Vert{}L_{\cdot k}\Vert{}_2^2 = 1$) instead of rescaling by $\text{sd}(Z_{\cdot k})$ mid-loop:**

     - In standard `flashier`, scale non-identifiability is handled cleanly between $L$ and $F$. If instead you normalize the projection covariate passed to Cox by dividing $Z_{\cdot k}$ by its fixed start-of-sweep norm $d_k^{(r)} = \sqrt{\frac{1}{N}\sum_{i=1}^N \bar{Z}_{ik}^2 + 10^{-8}}$ (and floor the Cox weights $w_i \ge 10^{-6}$ using log-sum-exp arithmetic), then $\eta_i$ **can never explode** to $39,086$, even with 25,000 or 100,000 CpGs!

##### D. Fixing Rank Selection ($K$) Without Ad-Hoc Post-Loop Deletion Scores

- Why did you need the complicated removal score $S$ (Eq. 36) in the multimodal model? Because single-modality used a fixed $K$ and multimodal tried to delete factors *after* the loop by zeroing a factor while holding the other $K-1$ factors frozen.
- **How to do it naturally in v2:**

    1. During the CAVI loop, because both $L_{\cdot k}$ and $F_{mk}$ have **Point-Exponential priors** (`ebnm_point_exponential`), `ebnm` automatically estimates the active proportion $\pi_{L,k}$ and $\pi_{F,mk}$. If $\bar{L}_{\cdot k}$ or all $\bar{F}_{mk}$ shrink to zero (or explain $< 0.1\%$ of variance with $\text{PIP}(\beta_k) < 0.05$), drop that factor *inside* the loop and let the remaining factors continue adapting!
    2. For selecting the starting rank $K_{\text{init}}$ (or tuning $\alpha$ across $\{0, 0.25, 0.5, 0.75\}$): use 5-fold cross-validated C-index on the **training set only** with the 1-SE rule (just like `DeSurv` did, without peeking at the external validation cohorts!).
    3. Use a **Point-Normal prior (`ebnm_point_normal`) on $\beta$** (where all $K$ coefficients share an empirical Bayes prior or have a spike at zero) so every retained factor has a true **Posterior Inclusion Probability ($\text{PIP}_k \in [0, 1]$)** for survival, rather than an arbitrary $\vert{}\bar{\beta}_k\vert{} > 0.001$ or $\vert{}\bar{\beta}_k\vert{}/\text{sd}(\beta_k) \ge 1.96$ rule.

#### 3. Addressing the Specific Audit Findings in the Writeup (Donor Overlap & Biology)

The writeup also flagged two really important data-hygiene details that are great to clean up now so your paper is bulletproof:

1. Repeated PACA-AU Donors in External Validation (Section 2.1, Table 2): The writeup notes that PACA-AU microarray ($n=63$) and PACA-AU RNA-seq ($n=52$) share 47 overlapping patients (and 1 patient has 2 RNA-seq records).

     - **Easy Fix:** Deduplicate that 1 repeated donor in PACA-AU RNA-seq ($n=51$), and when you report a pooled external validation HR or C-index across cohorts, either report the 5 cohorts separately (as in your Table 4) or restrict the pooled Kaplan-Meier / Cox model to the **568 unique donors** (e.g., using RNA-seq for the 47 shared PACA-AU donors and microarray for the 16 microarray-only PACA-AU donors). That completely eliminates the "donor dependence" caveat in the writeup!
2. Multimodal CpG Annotation & Biological Story (Section 4.2, Table 11):

     - Look at **Table 11 & Figure 7 (pages 26–27)**: even in your unconverged multimodal fit, **Program 2 (Immune: *CCL19, PTPRC, CCL21, LTB*) has 5.1% expression share and 20.9% methylation share**, and **Program 4 (Basal-like vs. Classical: *S100A2, KRT5, KRT17* vs. *CLDN18, REG4*) is strongly survival-active in both TCGA ($p=0.0055$) and ICGC ($p=0.00033$, HR $= 1.99$) with balanced expression (6.2%) and methylation (4.5%) shares**!
     - Once you map those top-loading CpG probes (`cgXXXXXXXX`) to their corresponding genes/promoters using the Illumina 450K annotation manifest (`IlluminaHumanMethylation450kanno.ilmn12.hg19`), you will be able to show which genes in Program 4 and Program 2 are epigenetically regulated (e.g., promoter hypomethylation of Basal keratins or *S100A2* coupled with expression activation)! That is a killer biological result for your paper.

#### Summary of the "Gap" (Current Implementation vs. Target v2 Method)

- **Gap 1 — Supervision of $F$ ($\alpha_F$):**

    - *Current:* Set to $\alpha_F = 0$ in single-modality (completely off), and unweighted in multimodal (making the survival precision $10^{-7}\times$ smaller than reconstruction precision, so it fails to find low-variance prognostic factors).
    - *Target v2:* Weight the Cox pseudo-likelihood by $\gamma = \frac{\alpha}{1-\alpha}\frac{\sum_m v_m P_m}{N_{\text{event}}}$, use sequential coordinate updates with the exact coupling term, and initialize at least one factor along the Cox score residual so low-variance prognostic programs are actively recovered.
- **Gap 2 — Single-Modality vs. Multimodal Discrepancy:**

    - *Current:* Two incompatible codebases—single-modality uses z-scored data (with negative values) + non-negative Point-Exponential $F$, while multimodal uses unstandardized data + intercepts + signed Point-Laplace $F_m$ (no longer NMF on $F_m$).
    - *Target v2:* Unify into **one codebase** where $Y_m \ge 0$ is kept non-negative (scaled by feature SD without centering) with **Point-Exponential priors on both $L$ and $F_m$** (true sparse NMF across all modalities), while centering $\tilde{Y}_m = Y_m - \bar{Y}_m$ *only* inside the Cox projection $\eta = \sum_m \tilde{Y}_m F_m \beta$.
- **Gap 3 — Multimodal Non-Convergence (500-Sweep Cap) & 25k-CpG Explosion:**

    - *Current:* Mid-loop rescaling by $\text{sd}(\bar{Z}_{\cdot k})$ fails to transform the $\beta_k$ prior and caches KL against stale priors (causing limit-cycle oscillations), while unnormalized sums in $Z_{ik}$ with unfloored Cox weights $w_i$ cause numerical overflow at 25,000 CpGs.
    - *Target v2:* Floor $w_i \ge 10^{-6}$ with log-sum-exp risk sets, normalize projection columns by start-of-sweep norms $d_{mk}$, weight modalities by $v_m = \bar{P}/P_m$, and compute KL directly from the `ebnm` return object at each update.
- **Gap 4 — Rank Selection ($K$) & Factor Removal:**

    - *Current:* Single-modality uses post-hoc $\vert{}\bar{\beta}_k\vert{} > 0.001$ cutoffs; multimodal evaluates factor deletion $\Delta_k$ holding the other $K-1$ factors frozen (which never deletes any factor on real data because of the 13,000-feature reconstruction penalty).
    - *Target v2:* Let `ebnm_point_exponential` prune dead factors *during* CAVI when $\hat{\pi}_k \to 0$, use `ebnm_point_normal` on $\beta$ to give exact Posterior Inclusion Probabilities ($\text{PIP}_k$) for survival, and select $(K_{\text{init}}, \alpha)$ strictly via training CV (1-SE rule) without touching external validation cohorts.
- **Gap 5 — Data Hygiene & Preprocessing Bugs:**

    - *Current:* TCGA, Dijk, and PACA-AU RNA-seq are accidentally $\log_2(x+1)$-transformed twice, and the 616 external validation records include 48 repeated PACA-AU donor records.
    - *Target v2:* Fix the loader to apply a single $\log_2(x+1)$ transform across all RNA-seq cohorts, deduplicate the 568 unique external validation donors for pooled tests, and annotate the multimodal CpG probes to gene promoters.

Would you like me to write out the clean, self-contained v2 R implementation (`cavi_ebmf_cox.R`) and its companion `documentation.md` so you have a rock-solid, bug-free foundation to test on your simulated and real PDAC datasets?

---

## Turn 4

### User

*Attachments:* `ROADMAP.md`, `PROJECT_STATUS.pdf`, `DECISIONS.md` (as of 2026-10-17)

These attached files should also help give context to the progression of the progress from beginning to present, the current state of the projects, and some decisions that have been made along the way (not claiming by any stretch that those decisions are correct, effective, or should stay as the final decision)

### Gemini

Reading `ROADMAP.md`, `PROJECT_STATUS.pdf`, and `DECISIONS.md` from March through October 16, 2026 puts the entire puzzle together. You can see the exact evolutionary arc of the project:

- **March–April (Cluster A / `LB`):** Started with $\eta = L\beta$, discovered the $5000\times$ scale imbalance between $A_{\text{gen}}$ and $A_{\text{surv}}$ in $L$, and found that test-time pseudo-inverse projection didn't match training.
- **May–August (Cluster B / `YFB` with $\alpha_F = 0$):** Switched to $\eta = (YF)\beta$, found that z-standardizing per cohort (`D4`) + normalizing $F$ columns to unit $L_2$ norm made $\alpha_F = 0$ super stable (`C = 0.627`, finding Program 7 Basal-like adverse and Program 3 Classical protective).
- **September–October 16 (Multimodal + Workstream C `supervised-F`):** Advisor meetings on 9/18 and 10/2 established that $\alpha_F = 0$ **cannot be the final method** because it never uses survival to shape $F$, and introduced the multimodal expression + methylation model (`fit_multimodal_yfb.R`) with intercepts, signed Point-Laplace priors, scale canonicalization $\text{sd}(\bar{Z}_{\cdot k}) = 1$, and ELBO nullcheck pruning.

Most importantly, the October 15–16 entries in `DECISIONS.md` show why the repository feels like it has accumulated massive complexity (710 tests, 176,000+ lines of CSV traces, `Task2-1` through `Task2-9` diagnostic patches) while **three core blockers remain unresolved**:

#### 1. The Exact Root Cause of Your Current October 15–16 Blockers

##### Blocker 1: Why the "Shared $\beta$ Prior" (`shared_point_laplace`) Collapses All $\beta \to 0$ at Outer 99 on Real Multimodal Data (`Task2` / M4)

Look at the `2026-10-15` entries in `DECISIONS.md`:

- In simulation (`n = 400`, strong planted $\beta_3 = 3.3$), `M4` (`init_supervised = TRUE`, `beta_prior = "shared_point_laplace"`, `beta_warmup = 10`) worked great (recovered 18/20 and 20/20).
- But when you ran that exact `M4` configuration on **real TCGA multimodal data** (`n = 144`, 75 deaths, $K = 7$), **at outer sweep 99 / inner sweep 1, the shared Point-Laplace slab probability $\hat{\pi}_\beta$ hit exact $0$, permanently zeroing out every single $\beta_k = 0$ for the rest of the 500 sweeps** (producing constant risk scores and `C = NA` on ICGC)!
- **Why did that happen mathematically?**

    1. **Empirical Bayes prior estimation on $K = 7$ observations has a phase transition at small SNR:** When you fit a shared `ebnm_point_laplace` (or `ebnm_point_normal`) prior by maximum marginal likelihood across only $K = 7$ scalar coordinates $(\beta_1, \dots, \beta_7)$, if 6 of the 7 factors are non-prognostic ($\beta_k = 0$) and the 1 prognostic factor (Program 4) has a modest real-data z-score (e.g., $|x_{\beta, 4} / s_{\beta, 4}| \approx 2.1$ in $n = 144$ with 75 deaths, especially when $A_{\beta, k} = \sum_i w_i (\bar{Z}_{ik}^2 + \text{Var}(Z_{ik}))$ is inflated by $\text{Var}(Z_{ik})$ across 13,000 features!), the marginal log-likelihood $\sum_{k=1}^7 \log \left( (1-\pi)\mathcal{N}(x_k \mid 0, s_k^2) + \pi \text{Slab}(x_k) \right)$ is **globally maximized at $\hat{\pi} = 0$** whenever $\sum_{k=1}^7 (\text{BF}_k - 1) \le 0$!
    2. **The Absorbing State Trap:** Once `ebnm` returns $\hat{\pi}_\beta = 0$ on a single inner sweep, the posterior moments for all $k=1\dots 7$ become $\bar{\beta}_k = 0$ and $\overline{\beta_k^2} = 0$. On the very next sweep, $F$ receives **zero** survival gradient ($\bar{\beta}_k = 0$), so $F$ drifts purely toward unsupervised reconstruction (which is dominated by the 10,000 CpGs!), which reduces $x_{\beta, k}$ even further, locking $\hat{\pi}_\beta = 0$ forever!
    3. **How to fix this permanently (No Absorbing Zero State):**

        - **Never let $\pi_\beta$ be estimated unconstrained by MLE on $K=7$ noisy points!** Either:

            - Fix the prior inclusion probability $\pi_\beta \in [0.25, 0.50]$ (e.g., expecting ~2 out of 7 factors to be prognostic) or enforce a lower bound $\pi_\beta \ge \frac{1}{K}$ in `ebnm` (`g_init` / constrained mixture proportion), **OR**
            - Use a **zero-mean Normal prior with a weak Inverse-Gamma / half-Cauchy prior (or lower floor $\sigma_{\beta, k}^2 \ge \sigma_{\min}^2 > 0$) per factor** instead of the unregularized $\hat{\sigma}_k^2 = \max(0, x_k^2 - s_k^2)$ formula! Notice that $\hat{\sigma}_k^2 = \max(0, x_k^2 - s_k^2)$ (Appendix A.5 / B.7) sets $\hat{\sigma}_k^2 = 0$ (a point mass at zero!) the instant $|x_k| \le s_k$ (i.e. whenever a factor's survival Wald $|z| \le 1$)! In real cancer survival data with $n=144$ (or CV fold $n=115$ as in `Task2-3` Fold 1), a truly prognostic factor can easily dip below $|z| = 1$ during early CAVI sweeps while $F$ is still rotating—and $\max(0, x^2 - s^2)$ instantly kills it!

##### Blocker 2: Why Canonicalization ($\text{sd}(\bar{Z}_{\cdot k}) = 1$) + `w_S = 1` Still Fails to Converge (500-Sweep Cap) and Oscillates on Both Single & Multimodal Fits

Look at `Task2-6`, `Task2-7`, and `Task2-9` (October 16 entries in `DECISIONS.md` and `PROJECT_STATUS.pdf`):

- Codex built an exact complete-cycle observer (`Task2-7`) and found:

  > *"No late fixed-prior posterior step loses the reference objective; **F-prior refits lose it 50/350 times and working refreshes 18/50. Joint L/beta comparisons have unresolved losses; transformed reference priors do not certify a common implemented objective.**"*
- And in `Task2-9`, even after fixing the point-Laplace prior step (`fixed-q`), **the relative output movement stayed stuck at $0.01915$ against the $10^{-4}$ tolerance**!
- **Why are your outputs oscillating by ~1–2% every sweep instead of converging?** There are **three interacting reasons** in `fit_multimodal_yfb.R` / `fit_cox_on_yf(canonicalize = TRUE)`:

    1. **Rescaling by $c_k = \text{sd}(\bar{Z}_{\cdot k})$ is NOT a likelihood-invariant transformation of the variational problem when $Z_{ik} = \tilde{Y}_i^T F_{\cdot k}$ is a data projection!**

        - Wait—look closely at Appendix B.9 (Eq. 35): the code claims that $(L_{\cdot k}, F_{mk}, \beta_k) \mapsto (c_k L_{\cdot k}, F_{mk} / c_k, c_k \beta_k)$ leaves $L F_m^T$ and $\eta = \tilde{Y} F \beta$ unchanged.
        - That is true for the *point estimates*, **but look at what it does to the second moments and the next iteration's precision**: In Step 1 (`update_L`), $A_{L,k} = \sum_{m,j} \tau_{mj} \overline{f_{mjk}^2}$. When you divide $F_{mk}$ by $c_k = \text{sd}(\bar{Z}_{\cdot k})$, $\bar{Z}_{\cdot k} = \sum_m \tilde{Y}_m \bar{F}_{mk}$ depends on the **cross-feature covariance of $\tilde{Y}_m$ along $\bar{F}_{mk}$** ($\bar{F}_{mk}^T \widehat{\text{Cov}}(Y_m) \bar{F}_{mk}$)! Every time $F_{mk}$ rotates slightly toward correlated genes/CpGs during the Gauss-Seidel sweep, $\text{sd}(\bar{Z}_{\cdot k})$ jumps up, which divides $F_{mk}$ down, which inflates $L_{\cdot k}$ up, which changes $A_{mjk}^{\text{rec}} = \tau_{mj} \sum_i \overline{l_{ik}^2}$ on the next sweep!
        - Even worse, in `M4` (`shared_point_laplace` on $\beta$), when factor $k$ is rescaled by $c_k$ (which is *different* for each factor $k=1\dots K$), **a single shared prior $g_\beta$ across all $K$ factors cannot be pushforward-transformed by 7 different scale constants $c_1, \dots, c_7$ simultaneously** (the exact bug flagged in `DECISIONS.md` 2026-10-15: *"how factor-specific canonicalization interacts with an untransformed common beta prior"*)!
    2. **The Undamped Gauss-Seidel Sweep Over 13,000 Correlated Features:**

        - In each inner sweep, updating 3,000 genes and 10,000 CpGs sequentially while holding $L$ and $\beta$ fixed causes the projection $\bar{Z}_{\cdot k}$ to drift substantially before $\beta_k$ and $L_{\cdot k}$ get to react.
    3. **Modality Imbalance ($10,000$ CpGs vs. $3,000$ Genes) in `update_L` and `prune_nullcheck`:**

        - Look at Eq. 29 (Appendix B.5): $A_{L,k} = \sum_{m,j} \tau_{mj} \overline{f_{mjk}^2}$ and $B_{L,ik} = \sum_{m,j} \tau_{mj} R_{mij}^{(-k)} \bar{f}_{mjk}$ sum over 3,000 genes and 10,000 CpGs **without dividing by $P_m$**!
        - And in `Open item 2 (M5 decision)` (`2026-10-15` in `DECISIONS.md`): when you ran `prune_yf_nullcheck()` on single-modality, **it deleted the planted low-variance prognostic program** because *"the planted program's gain in the Cox term is smaller than its KL and reconstruction cost"*! Why? Because in the ELBO removal score $S$ (Eq. 36), the KL divergence $\sum_{j=1}^P \text{KL}(q(f_{jk}) \,\|\, g_F)$ sums over $P = 2,064$ (or $13,000$) features, while $\log \text{PL}(\bar{\eta})$ is a single survival log-likelihood over $75$ events ($\approx -304$ vs. $610,423$ reconstruction and $130,212$ KL, as shown on page 20 of your writeup!)! Any low-variance prognostic factor pays a KL penalty across $P$ features to have non-zero loadings $F_{\cdot k}$, which easily exceeds its $+15$ improvement in unweighted Cox partial log-likelihood!

#### 2. Strategic Assessment: How to Get Out of the "Patch Spiral"

Right now, your workflow with Claude and Codex is trapped in what software engineers call a **patch spiral**:

- Every time a simulation cell or CV fold fails one of the rigid gate criteria (e.g., Fold 1 hits $\beta=0$, or seed 20261026 oscillates at `19/20` converged, or `prune_nullcheck` deletes the planted factor), a new opt-in patch is bolted onto the existing code (`init_supervised`, `beta_warmup`, `protected_pruning`, `full_state_refit`, `profiled_curvature`, `fixed_q_loading_prior`, `log_space_cox`).
- Because those patches are layered on top of a CAVI loop whose **scale canonicalization ($c_k = \text{sd}(Z_k)$), unweighted likelihoods ($O(NP)$ vs $O(N)$), and point-mass variance estimates ($\max(0, x^2 - s^2)$)** are fundamentally fighting each other, fixing one symptom breaks another gate!

Here is the clearest, most honest way to resolve your big-picture questions and bring this project to a clean, publishable finish line:

##### 1. How to Make Supervision of $F$ Work Reliably (Without 10 Fragile Knobs)

Instead of relying on `canonicalize = TRUE` ($\text{sd}(\bar{Z}_{\cdot k}) = 1$) + `w_S = 1` + `shared_point_laplace`, which has proven unstable on real TCGA multimodal data:

- **Fix the scale of $L_{\cdot k}$ (or $F_{\cdot k}$) using Euclidean $L_2$ norm, NOT data-dependent $\text{sd}(Y F_{\cdot k})$!**

    - Why was your original `D4` normalization ($d_k = \|F_{\cdot k}\|_2$, or in multimodal, setting $\frac{1}{N}\sum_{i=1}^N \overline{l_{ik}^2} = 1$) so much more stable than $c_k = \text{sd}(\bar{Z}_{\cdot k})$? Because $\|L_{\cdot k}\|_2$ and $\|F_{\cdot k}\|_2$ are **pure coordinate norms**—they do not depend on the empirical covariance matrix $Y^T Y$!
    - If you fix $\frac{1}{N} \sum_{i=1}^N \overline{l_{ik}^2} = 1$ at the end of each $L_{\cdot k}$ update (transferring the scalar multiplier $c_k$ into $F_{mk}$ and its prior rate $\lambda_{mk}$, while leaving $Z_{ik} / d_k$ normalized by a fixed sweep-level scale $d_k$), then **$\beta_k$ is completely decoupled from the $L \leftrightarrow F$ rescaling**! You never have to rescale $\beta_k$ or its prior variance $\sigma_k^2$!
- **Replace $\hat{\sigma}_k^2 = \max(0, x_{\beta, k}^2 - s_{\beta, k}^2)$ with a Non-Collapsing Prior Floor (or Fixed-Slab Point-Normal with $\pi_\beta \ge 1/K$):**

    - As your `2026-10-15` diagnosis proved, both $\max(0, x^2 - s^2)$ and unconstrained `shared_point_laplace` collapse $\beta \to 0$ on real data whenever the initial SNR is modest.
    - If you use a per-factor Normal prior $\beta_k \sim \mathcal{N}(0, \sigma_k^2)$ with an Inverse-Gamma $(\alpha_0, \beta_0)$ empirical Bayes update:

      $$
      \hat{\sigma}_k^2 = \frac{\beta_0 + \frac{1}{2}(\bar{\beta}_k^2 + \text{Var}(\beta_k))}{\alpha_0 + \frac{1}{2}}
      $$

      (with small prior pseudo-counts, e.g., $\alpha_0 = 1, \beta_0 = 0.05$), **$\hat{\sigma}_k^2$ is strictly positive ($\ge \frac{\beta_0}{\alpha_0 + 0.5} > 0$) at every iteration**! It can **never** hit an absorbing zero state at outer sweep 99 or in CV Fold 1, while still strongly shrinking null factors toward zero!

##### 2. How to Handle $K$ Selection & Factor Classification (Stop Deleting Factors Mid-Refit!)

Look at what your October 15 experiments proved (`Task1-4`, `Task1-5`, `DECISIONS.md` 2026-10-15):

- Whenever you delete a factor using `prune_nullcheck` and refit the $K-1$ model, either:

    1. Unprotected pruning **deletes the low-variance prognostic factor** because its $P$-dimensional KL cost in $F$ outweighs its $N$-dimensional Cox log-likelihood gain, OR
    2. Protected pruning deletes a background factor, and during the $K-1$ refit, **the protected prognostic factor gets hijacked to explain the deleted factor's genomic variance** (`loading |cor| 0.896 at K=7 becomes 0.391 at K=6`)!
- **Think about what that means scientifically:** That is not a bug—it is a fundamental mathematical property of low-rank matrix factorization! If the data truly has 6 high-variance biological compartments (e.g., Classical, Basal, Exocrine, Immune, Stroma, Methylation-drift) plus 1 low-variance prognostic program ($K_{\text{true}} = 7$), and you force the model down to $K = 5$ or $6$ factors, the low-variance factor **must** get overwritten by the larger variance components!
- **Your advisor's comment on 2026-10-15 in `DECISIONS.md` already gave you the exact right answer**:

  > *"Andrew also clarifies that retaining seven joint programs can be appropriate if their expression, methylation and prognostic contributions justify them. Audit roles and redundancy rather than treating deletion as the objective. The Single write-up's four retained programs are classifications of a K=7 fit, not four factors after deletion/refitting."*
- **Recommendation:** **Drop post-hoc factor deletion + refitting (`prune_nullcheck`) entirely!**

    - Fit at $K_{\text{init}}$ (selected via training 5-fold CV 1-SE rule, e.g., $K = 7$).
    - Let the coordinate-wise Point-Exponential / Point-Laplace priors naturally zero out dead columns if $K_{\text{init}}$ is truly over-specified.
    - Classify the fitted factors into **Survival-Active**, **Genomics-Only (Reconstruction-Active)**, and **Inactive/Dead** without deleting and refitting! That immediately eliminates the entire `Task1` protected-pruning refit loss problem!

##### 3. Balancing Modalities in the Multimodal Model (Expression + Methylation)

In `fit_multimodal_yfb.R`, why does the model struggle when you expand from 10,000 CpGs to 25,000+ CpGs?

- Because in $Z_{ik} = \sum_{m,j} \tilde{y}_{mij} f_{mjk}$ (Eq. 24) and $A_{L,k} = \sum_{m,j} \tau_{mj} \overline{f_{mjk}^2}$ (Eq. 29), **every CpG probe is given the same weight as a gene**! When $p_{\text{meth}} = 25,000$ (or $300,000$) and $p_{\text{expr}} = 3,000$, methylation has $8\times$ to $100\times$ more terms in the sum than expression!
- **The Fix:** Normalize each modality block's contribution to both the reconstruction likelihood and the linear predictor $Z_{ik}$ by $1 / p_m$ (or $\bar{p} / p_m$). Once you scale by $1 / p_m$:

    - Increasing CpGs from 10,000 to 25,000, 50,000, or 300,000 **keeps the variance of $Z_{ik}$ and the magnitude of $A_{L,k}$ completely invariant to $p_m$**!
    - Your 25k/50k/100k CpG sensitivity curve (`Task5`) will never blow up its predictor SD to $39,086$!

#### Summary of Where Your Project Stands & What to Do Next

- **Where the Project Actually Stands (The Big Picture):**

    - **Single-Modality (`D4`, $\alpha_F = 0$, $K=7$):** Fully working, converged, and biologically validated across 5 external cohorts (`mean C = 0.627`, beating the $K=40$ LASSO two-step comparator `0.601`, and cleanly recovering Basal-like Adverse Program 7 and Classical Protective Program 3).
    - **Supervised $F$ (`Workstream C`) & Multimodal (`M4`):** Currently stuck in a loop of small diagnostic patches (`Task1`–`Task5`) because of three interacting mathematical choices in the October codebase: (1) unconstrained shared-prior MLE on $K=7$ betas (`shared_point_laplace`) and $\max(0, x^2 - s^2)$ which hit an absorbing $\beta = 0$ state when real-data SNR is moderate; (2) data-covariance-dependent rescaling $c_k = \text{sd}(\bar{Z}_{\cdot k})$ every sweep which fights the priors and prevents $10^{-4}$ convergence; and (3) `prune_nullcheck` factor deletion + refitting which allows deleted high-variance factors to overwrite low-variance prognostic factors.
- **How to Resolve Each Open Item Cleanly:**

    - **Selecting $K$ & Factor Pruning:** Stop deleting factors and refitting (`prune_nullcheck`), which destroys low-variance prognostic programs when high-variance compartments are removed. Keep your three-way classification (**Survival-Active**, **Genomics-Only**, **Dead/Inactive**) on the single $K_{\text{init}}$ fit (chosen by training CV 1-SE rule), as already endorsed in the 2026-10-15 `DECISIONS.md` entry.
    - **Fixing $\beta \to 0$ Collapse in Real Multimodal & CV Folds (`Task2`):** Replace the absorbing point-mass variance estimator $\hat{\sigma}_k^2 = \max(0, x_k^2 - s_k^2)$ and unconstrained `shared_point_laplace` ($\hat{\pi}_\beta \to 0$) with an **Inverse-Gamma regularized Normal prior variance** $\hat{\sigma}_k^2 = \frac{\beta_0 + \frac{1}{2}\overline{\beta_k^2}}{\alpha_0 + \frac{1}{2}} > 0$ (or a lower floor on $\pi_\beta \ge 1/K$) so $\beta_k$ can never permanently die during early CAVI sweeps.
    - **Fixing Convergence Caps (`500-sweep cap` & `0.019` oscillation):** Stop rescaling $\bar{Z}_{\cdot k}$ and $\beta_k$ by the data-dependent $\text{sd}(\bar{Z}_{\cdot k})$ inside the inner loop. Instead, fix the Euclidean scale of $L_{\cdot k}$ ($\frac{1}{N}\Vert{}L_{\cdot k}\Vert{}_2^2 = 1$), normalize the projection by a fixed sweep-level feature scale, and apply a small damping step ($\rho \in [0.5, 0.7]$) on the sequential $F$ sweep when supervision is active.
    - **Fixing the 25k–300k CpG Expansion (`Task5`):** Scale each modality $m$'s precision and projection contribution by $\bar{p} / p_m$ so that expanding $p_{\text{meth}}$ from $10,000$ to $25,000$ or $300,000$ does not multiply the linear predictor variance or $L$-update precision by $2.5\times$ to $30\times$.
    - **Eliminating Repo Bloat:** As your October 16 audit found, **82.3% of the 176,342 added lines are CSV trace dumps** and only **0.9% (~1,500 lines) are production R code**. You don't need to throw away your test suite or reports—just archive the diagnostic trace CSVs and consolidate `fit_cox_on_yf.R` and `fit_multimodal_yfb.R` around a single, stable mathematical core.

Would you like to walk through the exact 3-line mathematical change to `fit_multimodal_yfb.R` / `fit_cox_on_yf.R` that stops the `Task2` $\beta \to 0$ collapse and the `0.019` oscillation, or would you prefer to do the plain-English "teach-back" walkthrough of the 5 biological characterization methods flagged in your `ROADMAP.md`?

---

## Turn 5

### User

I haven’t actually read your responses yet. However, based on all of the information I’ve provided to you your task now is to succinctly summarize the following:

the goals of the project
the current implementation and pinpointed strengths and weaknesses where issues present themselves
recommendations of how to address each issue identified and which components of the project that should remain as is
a final summary that presents the key recommended characteristics of the “best case” implementation that addresses the requirements/features that I specified are necessary
this will end up being a story of “what do we need”, “where are we now”, “what needs to change”, “what are the details of the final best solution”

### Gemini

Here is the complete synthesis of your project, structured around the four-part story you requested: **What We Need**, **Where We Are Now**, **What Needs to Change (vs. Stay)**, and **The Final "Best-Case" Blueprint**.

### 1. What Do We Need? (Project Goals & Core Requirements)

- **Joint Bayesian NMF + Survival Architecture ($\eta = YF\beta$):**

    - Decompose molecular data $Y \approx LF^T + E$ into non-negative patient activities ($L$) and non-negative feature programs ($F$) using Empirical Bayes Matrix Factorization (EBMF), jointly coupled to a Cox proportional hazards model.
    - Parameterize the Cox linear predictor as the feature projection $\eta = YF\beta$ (rather than $L\beta$) so that out-of-sample validation cohorts can be scored directly via $\hat{\eta}_{\text{new}} = Y_{\text{new}}\bar{F}\bar{\beta}$ without running an inner optimization loop on test patients.
- **True Survival Supervision of Feature Programs ($F$):**

    - Allow the Cox survival likelihood to actively guide the learning of $F$ ($\alpha_F > 0$) so the model isolates **low-variance prognostic gene/CpG programs** that are masked by high-variance, outcome-neutral signals (like tumor purity or exocrine tissue composition in PDAC).
    - Demonstrate superior prognostic discrimination (C-index) and cleaner biological program recovery compared to a two-step unsupervised pipeline (unsupervised EBMF followed by Cox regression).
- **Principled Bayesian Rank Selection ($K$) & Division of Factor Roles:**

    - Select model complexity using training data only (avoiding validation leakage) and leverage Empirical Bayes shrinkage priors to automatically separate factors into **Survival-Active programs** (prognostic), **Genomics-Only programs** (outcome-neutral background reconstruction), and **Inactive/Dead factors**.
- **Multi-Cohort & Multi-Modal Generalization:**

    - **Multiple Training Cohorts:** Train across pooled cohorts (TCGA + CPTAC, $n = 273$) without batch offsets corrupting the biological programs or collapsing survival coefficients.
    - **Multiple Modalities (Expression + Methylation):** Jointly factorize matched gene expression and DNA methylation matrices ($Y_{\text{expr}}, Y_{\text{meth}}$) using shared patient loadings $L$ and modality-specific feature loadings $F_m$ to uncover coordinated epigenetic and transcriptional survival programs.
- **Feature Scalability & Biological Interpretability:**

    - Scale reliably across focused feature sets (~2,000–3,000 genes and ~10,000 CpGs) and full feature universes (~17,900 genes and ~305,000 CpGs) without numerical overflow.
    - Produce biologically interpretable signatures validated via pathway enrichment (`fgsea`/ORA), PDAC subtype classifiers (`PurIST`/`DeCAF`), and out-of-sample survival stratification.

### 2. Where Are We Now? (Current Implementation: Strengths & Pinpointed Weaknesses)

#### What Is Working Well (Strengths)

- **Strong Single-Modality Predictive Benchmark (`D4`):**

    - Training on pooled TCGA + CPTAC ($n = 273$, 2,064 DeSurv-screened genes, per-cohort z-standardization, $K = 7$) converges in 20 sweeps and achieves a mean external C-index of **0.627** across five independent PDAC cohorts, outperforming the 40-factor EBMF + LASSO Cox two-step comparator (**0.601**).
- **Validated Single-Modality PDAC Biology:**

    - Five independent biological checks confirm that the retained factors capture established PDAC biology: **Program 7** is an adverse Basal-like/squamous and glycolytic/RTK-signaling program ($\text{HR} > 1$ in 5/5 external cohorts), **Program 3** is a protective Classical epithelial program ($\text{HR} < 1$ in 5/5 cohorts), and **Programs 5 & 6** capture genomics-only stromal (EMT/ECM) and immune/vascular variation.
- **Rigorous Testing & Audit Infrastructure:**

    - The repository has a modular CAVI design backed by 710 passing unit and integration tests, exact reproducibility checks, and honest diagnostic logging.

#### Where Issues Present Themselves (Pinpointed Weaknesses)

- **Issue 1 — Supervision of $F$ Is Turned Off in the Benchmark ($\alpha_F = 0$) and Swamped ($10^{-7}$) in Multimodal:**

    - In the single-modality benchmark, $\alpha_F = 0$, meaning survival never enters the $F$ update and the loadings are purely unsupervised EBMF.
    - Earlier attempts to turn on $\alpha_F > 0$ in single-modality dropped the within-factor feature coupling term and updated all genes simultaneously (Eq. 18), causing immediate overshooting and RMSE blowup.
    - In the multimodal model, the exact sequential Gauss-Seidel update was implemented with equal 1:1 weighting (Eq. 30–31), but because the reconstruction precision sums over $N \times P$ while survival sums over $N$, the **median survival fraction of loading precision is $2.03 \times 10^{-7}$ for expression and $9.00 \times 10^{-10}$ for methylation**. Consequently, the model recovers **0%** of planted low-variance prognostic programs in lower-noise multimodal simulations.
- **Issue 2 — $\beta \to 0$ Absorbing Collapse on Real Multimodal Data & CV Folds:**

    - Both the empirical Bayes Normal variance estimator $\hat{\sigma}_k^2 = \max(0, x_k^2 - s_k^2)$ (which sets $\hat{\sigma}_k^2 = 0$ whenever $\vert{}x_k\vert{} \le s_k$) and the shared Point-Laplace prior (`shared_point_laplace` in `M4`) have an **absorbing point-mass zero state**.
    - When `M4` is fit on real TCGA multimodal data ($n = 144$, $K = 7$), the shared slab probability $\hat{\pi}_\beta$ drops to exact $0$ at outer sweep 99 / inner sweep 1, permanently locking all $\beta_k = 0$ and producing constant risk scores (`C = NA`) on ICGC validation.
- **Issue 3 — Multimodal Non-Convergence (500-Sweep Cap) & Prior/Scale Bookkeeping Conflicts:**

    - Almost all real-data multimodal fits fail to converge and hit the 500-sweep cap.
    - **Root cause:** Rescaling each factor by data-dependent $c_k = \text{sd}(\bar{Z}_{\cdot k})$ after every factor update changes the coordinate scale as features rotate, fails to transform a shared $\beta$ prior across different $c_k$'s, and caches KL divergences against stale priors before refitting $g_{mk}$. This prevents the algorithm from ascending a single coherent ELBO and traps it in a ~1.9% oscillation (`0.01915` vs. `1e-4` tolerance).
- **Issue 4 — Counterproductive Post-Hoc Factor Deletion (`prune_nullcheck`):**

    - On real multimodal data, `prune_nullcheck` evaluates deleting factor $k$ while holding the other $K-1$ factors frozen (Eq. 36–37); because 13,000 features contribute ~610,423 to the reconstruction score versus ~$-304$ for survival, deleting any factor heavily penalizes reconstruction, so **zero factors are ever removed** (`Final rank = Starting rank` for all $K$).
    - Conversely, in single-modality simulations, `prune_nullcheck` either deletes the planted low-variance prognostic program (whose $P$-feature KL cost exceeds its $N$-sample Cox gain) or, under protected pruning, loses the planted program during the $K-1$ refit (`|cor|` drops from 0.896 to 0.391) because the remaining factors rotate to absorb the deleted high-variance factor's signal.
- **Issue 5 — Unweighted Modality Sums & 25,000-CpG Numerical Explosion:**

    - In the multimodal model, $Z_{ik} = \sum_{m,j} \tilde{y}_{mij} f_{mjk}$ and $A_{L,k} = \sum_{m,j} \tau_{mj}\overline{f_{mjk}^2}$ sum directly over all features without normalizing by $p_m$.
    - Increasing methylation from 10,000 to 25,000 CpGs inflates $\hat{\eta}_i$ to a span of 5,734 units by sweep 6 (underflowing standard exponential risk sets) and, even with log-space arithmetic, drives the predictor SD to **39,086**.
- **Issue 6 — Split Codebases, Non-Negativity Mismatch, & Data Hygiene Bugs:**

    - Single-modality and multimodal use incompatible models (Table 1): single-modality puts non-negative Point-Exponential priors on z-standardized data (which contains negative values and no intercept), whereas multimodal uses intercepts and signed Point-Laplace priors on $F_m$ (abandoning non-negative parts-based NMF on $F_m$).
    - The internal data loader and preprocessing flags accidentally apply $\log_2(x+1)$ **twice** to TCGA, Dijk, and PACA-AU RNA-seq, and the 616 external validation records include 47 shared PACA-AU donors across microarray/RNA-seq plus 1 duplicate RNA-seq donor (568 unique people).

### 3. What Needs to Change vs. What Should Remain As Is

#### Components That Should Remain As Is (Keep!)

- **The $\eta = YF\beta$ Projection Formulation:** Keep the projection predictor as the core survival link so test patients are scored in closed form.
- **1-Step Unsupervised Feature Pre-Screening (Mean + Variance Rank):** Keep the `DeSurv`-aligned unsupervised gene screening (top 3,000 per training cohort by combined mean and variance rank) as your primary filter. Do *not* switch to 2-step supervised pre-filtering by survival p-values as your default, as pre-screening by survival undercuts the core thesis that your joint factorization discovers prognostic programs without pre-filtering.
- **Three-Way Post-Fit Factor Classification (Without Deletion Refits):** Keep `classify_factors()` (**Survival-Active**, **Genomics-Only**, **Inactive/Dead**) on the fitted $K_{\text{init}}$ model. As noted in your October 15 advisor decision, retaining ~7 joint programs and classifying their distinct reconstruction and survival roles is scientifically preferable to forced deletion.
- **Cohort Handling via Per-Cohort Scaling & Optional Stratification:** Keep per-cohort standardization and optional cohort-stratified Cox risk sets (`strata_id`) / partially pooled cohort betas (`beta_cohort_id` for diagnostics) rather than adding dummy $\{0,1\}$ cohort indicator columns to $L$, which waste factor capacity and hurt external C-index.
- **Biological Characterization Suite & Unit Tests:** Keep the 5-method biological validation pipeline (`fgsea`, ORA, PurIST concordance, DeSurv gene overlap, external HR checks) and the modular test suite.

#### Targeted Changes Needed (Fix!)

- **Fix 1 — Unify Single & Multimodal into One Non-Negative Architecture:**

    - Keep the reconstruction matrices $Y_m \ge 0$ strictly non-negative (e.g., $\log_2(\text{TPM}+1)$ and methylation $\beta$-values scaled by their within-cohort standard deviation *without* subtracting the mean, or using non-negative baseline shifts) so that **both $L$ and $F_m$ use sparse non-negative Point-Exponential priors (`ebnm_point_exponential`)**.
    - Center the features **only inside the Cox projection** ($\tilde{Y}_{m,ij} = Y_{m,ij} - \bar{Y}_{m,\cdot j}$, so $Z_{ik} = \sum_m v_m \tilde{Y}_m \tilde{F}_{mk}$ has mean zero), giving you true non-negative sparse gene/CpG programs ($F_m \ge 0$) *and* a centered, intercept-free Cox linear predictor simultaneously.
- **Fix 2 — Balance Modality Dimensions ($v_m$) and Cox Supervision Strength ($\gamma$):**

    - Weight each modality $m$'s reconstruction precision and projection sum by $v_m = \frac{\bar{p}}{p_m}$ (or $\frac{1}{p_m}$) so that 10,000–300,000 CpGs and 2,000–3,000 genes contribute equally to $L$ and $Z$, preventing predictor explosion when expanding CpG counts.
    - Weight the Cox working likelihood by $\gamma = \frac{\alpha}{1-\alpha} \frac{\sum_m v_m p_m}{N_{\text{event}}}$ (where $\alpha \in [0, 1)$ is tuned via training CV, matching `DeSurv`'s loss scaling) so the survival precision in the $F$ update is $O(1)$ relative to reconstruction rather than $10^{-7}$.
    - Initialize 1–2 factor columns of $F$ using the positive/negative univariate Cox score directions $\max(0, \pm \tilde{Y}^T u^{(0)})$ alongside the leading SVD/NMF components so low-variance prognostic directions start with $\bar{\beta}_k^{(0)} \neq 0$ and receive immediate survival supervision.
- **Fix 3 — Prevent $\beta \to 0$ Collapse via Non-Absorbing $\beta$ Priors:**

    - Replace the absorbing variance estimate $\hat{\sigma}_k^2 = \max(0, x_k^2 - s_k^2)$ and unconstrained 7-point MLE `shared_point_laplace` with either an **Inverse-Gamma regularized Normal prior** $\hat{\sigma}_k^2 = \frac{\beta_0 + \frac{1}{2}\overline{\beta_k^2}}{\alpha_0 + 1/2} > 0$ or a **Point-Normal / Point-Laplace prior with a lower-bounded active proportion $\pi_\beta \ge \frac{1}{K}$**. This guarantees $\overline{\beta_k^2} > 0$ during early CAVI sweeps so survival supervision never shuts off mid-fit.
- **Fix 4 — Stabilize CAVI Convergence (Euclidean Norm Scaling + Log-Space Cox + Atomic KL):**

    - Replace data-covariance-dependent rescaling $c_k = \text{sd}(\bar{Z}_{\cdot k})$ with **Euclidean column-norm scaling** ($d_{mk} = \sqrt{\frac{1}{p_m}\Vert{}F_{mk}\Vert{}_2^2 + 10^{-8}}$, or fixing $\frac{1}{n}\Vert{}L_{\cdot k}\Vert{}_2^2 = 1$), which decouples $\beta_k$ and its prior from inner-loop feature rotations.
    - Use log-sum-exp arithmetic with a curvature floor $w_i \ge 10^{-6}$ in the Cox working approximation, apply mild damping ($\rho = 0.5$) on the sequential $F$ update when $\alpha > 0$, and compute KL directly from the returned `ebnm` object so the ELBO is strictly consistent.
- **Fix 5 — Clean Up Data Preprocessing, Donor Deduplication, & Repo Bloat:**

    - Remove the duplicate $\log_2(x+1)$ call in the RNA-seq loader/preprocessing path.
    - Deduplicate the 568 unique external validation donors (resolving the 47 shared PACA-AU donors and 1 repeated RNA-seq record) for all pooled confidence intervals and log-rank tests.
    - Archive or compress the 144,000+ lines of diagnostic CSV traces that currently make up 82.3% of the branch diff.

### 4. Final Summary: Key Characteristics of the "Best-Case" Implementation

| Component | Recommended "Best-Case" Specification | Why This Solves Your Project Requirements |
|---|---|---|
| **1. Generative & Survival Model** | • **Reconstruction:** $Y_m = L F_m^T + E_m$, $(E_m)_{ij} \sim \mathcal{N}(0, \tau_{mj}^{-1})$ for modalities $m \in \{1, \dots, M\}$ ($M=1$ for single-modality; $M=2$ for Expr + Meth). • **Cox Linear Predictor:** $\eta_i = \sum_{k=1}^K \tilde{Z}_{ik} \beta_k$, where $\tilde{Z}_{ik} = \sum_{m=1}^M v_m \left(\tilde{Y}_m \frac{F_{mk}}{d_{mk}}\right)_i$ uses centered features $\tilde{Y}_m = Y_m - 1\bar{Y}_m^T$ and Euclidean loading norms $d_{mk}$. | Unifies single-modality and multimodal into **one codebase**; preserves non-negative NMF on $Y_m \ge 0$ while ensuring zero-mean projections $\tilde{Z}_{\cdot k}$ in the Cox model without intercept drift. |
| **2. Priors ($L, F_m, \beta, \tau_{mj}$)** | • **$L \in \mathbb{R}_+^{n \times K}$:** `ebnm_point_exponential` per factor. • **$F_m \in \mathbb{R}_+^{p_m \times K}$:** `ebnm_point_exponential` per modality and factor. • **$\beta \in \mathbb{R}^K$:** `ebnm_point_normal` (or `ebnm_point_laplace`) with lower-bounded slab weight $\pi_\beta \ge 1/K$, OR Normal with Inverse-Gamma regularized variance $\hat{\sigma}_k^2 > 0$. • **$\tau_{mj} \in \mathbb{R}_+^{p_m}$:** Feature-specific ML precision with floor $Q_{mj} \ge n \times 10^{-8}$. | Point-Exponential on $F_m$ enforces sparse, non-negative gene/CpG programs (automatically focusing on top marker features during projection), while the non-collapsing $\beta$ prior separates prognostic factors from reconstruction-only factors without ever hitting the $\beta = 0$ absorbing trap. |
| **3. Likelihood Balancing & Initialization** | • **Modality Weight:** $v_m = \bar{p} / p_m$ applied to modality $m$'s reconstruction precision and projection. • **Supervision Weight:** Cox working weights scaled by $\gamma = \frac{\alpha}{1-\alpha}\frac{\sum_m v_m p_m}{N_{\text{event}}}$. • **Hybrid Initialization:** $K-1$ factors initialized from non-negative SVD/NMF of $Y_m$; 1 factor initialized from univariate Cox score residuals $\max(0, \pm \tilde{Y}^T u^{(0)})$. | Prevents 10k–300k CpGs from drowning out 2k–3k genes; brings the $F$ survival precision from $10^{-7}$ up to $O(1)$ so $\alpha \in (0,1)$ genuinely shapes $F$; ensures low-variance prognostic factors are captured from iteration 1. |
| **4. CAVI Loop & Convergence** | • **Update Order:** Outer sweep updates stratified/pooled log-sum-exp Cox quantities $(w_i \ge 10^{-6}, z_i)$; inner loop cycles $k=1\dots K$ updating $q(L_{\cdot k}) \to$ sequential damped Gauss-Seidel $q(F_{mk}) \to q(\beta_k)$, followed by $\tau_{mj}$. • **KL & Stopping:** Compute KL immediately from `ebnm` output; stop when relative change in ELBO and $\bar{\eta}$ is $< 10^{-4}$. | Eliminates the 500-sweep oscillation caused by $\text{sd}(\bar{Z}_{\cdot k})$ rescaling and stale prior KL caching, guaranteeing smooth convergence across all feature counts. |
| **5. Selecting $K$ & Factor Classification** | • **Selecting $(K_{\text{init}}, \alpha)$:** 5-fold cross-validated C-index on training data only using the 1-SE parsimony rule (e.g., $K_{\text{init}} = 7$). • **No Post-Hoc Deletion Refits:** Disable `prune_nullcheck` refitting; classify fitted factors into **Survival-Active** ($\text{PIP}_k > 0.5$ or $\vert{}\bar{\beta}_k\vert{}/\text{SD}(\beta_k) \ge 1.96$), **Genomics-Only** ($\text{PVE}_k \ge 1\%$), and **Dead**. | Avoids validation leakage, prevents low-variance prognostic factors from being overwritten during $K-1$ refits, and cleanly reports both prognostic and background biological compartments. |
| **6. Preprocessing & Evaluation Hygiene** | • **Single Log Transform:** Apply $\log_2(x+1)$ exactly once to unlogged RNA-seq cohorts. • **Primary Feature Set:** ~2,000–3,000 DeSurv-ranked genes + ~10,000 variable/promoter CpGs (with full-universe sensitivity curves in Supplement). • **Donor-Aware Validation:** Score all 568 unique external PDAC donors and 50 matched primary-PDAC ICGC donors via frozen training parameters. | Fixes the historical double-log bug, removes repeated-donor bias in pooled p-values/CIs, and delivers a clean, reproducible story from simulation to multi-omics PDAC biology. |
