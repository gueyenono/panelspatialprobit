# 1. Model specification

This file defines the statistical models the package will estimate, the identification restrictions, default priors, and the quantities reported to users (probabilities and impacts). Estimation algorithms are in `02_estimation_algorithms.md`.

## 1.1 Notation and data layout

- Units $i = 1,\dots,N$, periods $t = 1,\dots,T$, $n = NT$ latent observations.
- **Stacking: time is the slow index and unit the fast index**, so $z = (z_1', \dots, z_T')'$ with $z_t \in \mathbb{R}^N$. This is the convention in Kapoor–Kelejian–Prucha, Baltagi–Egger–Kesina (BEK) and `spmixW`, and it makes every period an $N$-block.
- $W_t$ ($N \times N$): spatial weights for period $t$, with a zero diagonal and usually row-standardised. In the time-invariant case $W_t = W$, so the $n \times n$ operator is $I_T \otimes W$. That Kronecker product is **never formed**: the code always works per period with the $N \times N$ matrix.
- $M_t$: weights for the spatial error process (default $M_t = W_t$).
- $\Delta = \iota_T \otimes I_N$ maps unit effects $\alpha \in \mathbb{R}^N$ to the $n$-vector.

## 1.2 The general latent model (v1.0 family)

$$
z_t = \rho W_t z_t + X_t\beta + W_t X^{(d)}_t\theta + \gamma\, y_{t-1} + \eta\, W_t y_{t-1} + \alpha + \tau_t \iota_N + u_t,
\qquad u_t = \lambda M_t u_t + \varepsilon_t,
\qquad \varepsilon_{it} \sim N(0, v_{it}).
$$

Observation rules:

- **Binary:** $y_{it} = 1\{z_{it} > 0\}$.
- **Ordered** ($J \ge 3$ categories): $y_{it} = j \iff \kappa_{j-1} < z_{it} \le \kappa_j$, with $\kappa_0 = -\infty$, $\kappa_1 = 0$ (because $X$ has an intercept) and $\kappa_J = +\infty$. The free cut-points are $\kappa_2 < \dots < \kappa_{J-1}$. The summary also reports the "no intercept, free first cut-point" parameterisation used by `MASS::polr` / `ordinal::clm`.

Writing $S_t(\rho) = I - \rho W_t$ and $R_t(\lambda) = I - \lambda M_t$, the conditional latent density is

$$
R_t\big(S_t z_t - X_t\beta - W_tX^{(d)}_t\theta - \dots - \alpha\big) = \varepsilon_t .
$$

**Every variant in this file is a special case of $A z - B = \varepsilon$**, with sparse $A = \operatorname{blockdiag}_t(R_t S_t)$. The estimation engine is built on this observation (see `02`).

### Named spatial specifications (`model =` argument)

| `model` | $\rho$ | $\lambda$ | $\theta$ (Durbin $WX$) | Spillovers in impacts |
|---|---|---|---|---|
| `"sar"` | free | 0 | 0 | global, via $S^{-1}$ |
| `"sdm"` | free | 0 | free | global, via $S^{-1}(\beta_r I + \theta_r W)$ |
| `"sem"` | 0 | free | 0 | none (direct only) |
| `"sarar"` | free | free | 0 | global |
| `"sdem"` | 0 | free | free | local, via $\theta_r W$ |
| `"slx"` | 0 | 0 | free | local |
| `"none"` | 0 | 0 | 0 | none (non-spatial panel probit; useful as a baseline and for tests) |

### Unobserved heterogeneity (`effects =` argument)

| `effects` | Specification | Notes |
|---|---|---|
| `"pooled"` | $\alpha = 0$ | Equivalent to `spatialprobit` on the block-diagonal $I_T \otimes W$; used for cross-checks. |
| `"random"` (default) | $\alpha_i \overset{iid}{\sim} N(0, \sigma_\alpha^2)$ | Kakamu & Wago (2005), BEK (2016). |
| `"random"` + `cre = ~ x1 + x2` | $\alpha_i = \bar x_i'\pi + a_i$ (Mundlak), or all $x_{i1..iT}$ (Chamberlain) | Implemented by **appending columns to $X$**, with $\pi$ inside $\beta$. This is the standard nonlinear-panel answer to "I need fixed effects" (Wooldridge 2010, ch. 15). |
| `"spatial_random"` (Phase 4, optional) | $\alpha = \phi W\alpha + a$ | Spatially autocorrelated unit effects, as in the Wang & Kockelman (2009) regional effects. |

**Fixed effects are deliberately not offered.** With fixed $T$, unit intercepts in a probit cause incidental-parameter bias (Neyman–Scott), and flat priors on $\alpha_i$ do not remove it. The documentation will say this and point to correlated random effects (CRE). A frequentist large-$T$ fixed-effects estimator exists (concentrated partial ML, Arbia, Billé & Leorato 2023), but it is out of scope.

**Time effects** (`time_effects = TRUE`): $\tau_t$ enter as period dummies in $X$. With $N$ observations per period this causes no incidental-parameter problem.

### Dynamics (`dynamic =` argument, Phase 5)

- **Observed state dependence:** $\gamma y_{i,t-1}$, plus optionally a space-time lag $\eta (W y_{t-1})_i$. Both are data columns, so the sampler is unchanged. The first period is used only for conditioning, leaving $T-1$ estimation periods.
- **Initial conditions:** Wooldridge (2005), i.e. $\alpha_i = \psi y_{i0} + \bar x_i'\pi + a_i$. It is implemented, like the CRE, by appending columns, and is the default whenever `dynamic` is set.
- **Latent AR(1)** (stretch goal): $\gamma z_{i,t-1}$, as in Wang & Kockelman (2009). Here $A$ becomes block lower-bidiagonal. It stays sparse with the same log-determinant ($\log|A| = \sum_t \log|S_t|$, since $A$ is block-triangular), so the generic engine still applies. The initial latent $z_0$ needs a separate treatment.

### Error distribution

- Default: $v_{it} = 1$ (probit).
- Phase 4 option: Student-$t$ errors through the scale mixture $v_{it} \sim IG(r/2, r/2)$ (Geweke 1993; also used in LeSage's `sar_g`). This gives a robit-type link that is robust to outliers in the latent space. $r$ is user-fixed (default 4) to keep the scale identified.

## 1.3 Identification and parameter space

- **Scale:** $\operatorname{Var}(\varepsilon_{it}) = 1$ (probit normalisation). For ordered models, location is fixed by $\kappa_1 = 0$.
- **$\rho, \lambda$ range:** for row-standardised $W$ the admissible interval is $(1/\omega_{\min}, 1)$, where $\omega$ are the eigenvalues of $W$. The default is $(-1, 1)$ like `spatialprobit`, and users can set the range. For non-row-standardised $W$ the bounds come from the extreme eigenvalues, with a check at fit time.
- **SARAR with $M = W$:** $\rho$ and $\lambda$ are separately identified only through the $X$ structure. The model warns when $M$ is identical to $W$ and there are no Durbin terms.
- **Weak information:** binary data are not very informative about $\rho$. The panel dimension helps, but extreme prevalence (mostly 0s or mostly 1s) needs a warning and a prior-sensitivity note in the docs.
- **Time-invariant regressors** are fine under RE/CRE; the docs explain this.
- **Isolates** (zero rows in $W$): allowed, but they generate a warning, and row-standardisation leaves those rows at zero.

## 1.4 Default priors

All defaults are **proper**, which simulation-based calibration (SBC) requires (see `04`). Every prior can be changed through `psp_priors()`.

| Parameter | Default prior | Rationale |
|---|---|---|
| $\beta$, $\theta$, $\gamma$, $\eta$ | $N(0, 10^2 I)$ | Vague on the probit scale, but proper (`spatialprobit` uses an effectively improper $10^{12}$). |
| $\rho$, $\lambda$, $\phi$ | Beta$(1.01, 1.01)$ rescaled to $(r_{\min}, r_{\max})$ | The LeSage–Pace convention. Nearly flat, but keeps mass away from the boundary. |
| $\sigma_\alpha$ | Half-Cauchy$(0, 1)$ | Gelman (2006). Implemented conjugately via the inverse-gamma mixture of Wand et al. (2011); an IG$(a, b)$ option on $\sigma_\alpha^2$ is kept for comparability with Kakamu & Wago and BEK. |
| $\kappa_2 < \dots < \kappa_{J-1}$ | Independent $N(0, 10^2)$ truncated to the ordering constraint | Proper, like the normal threshold priors of Wang & Kockelman. A flat (improper) option is kept for comparability with Albert & Chib and `spatialprobit`. |

## 1.5 Reported quantities

### Probabilities

Write $m_t = S_t^{-1}(X_t\beta + \dots)$ for the systematic part of the reduced form, excluding $\alpha$. The general reduced-form variance per observation is

$$
s_{it}^2 = \sigma_\alpha^2\,[S_t^{-1}S_t^{-\top}]_{ii} + [A_t^{-1} V_t A_t^{-\top}]_{ii},
$$

which simplifies to $(1+\sigma_\alpha^2)[(S'S)^{-1}]_{ii}$ for SAR with $V = I$.

- **Population-averaged probability (default):** $P(y_{it} = 1) = \Phi(m_{it}/s_{it})$, with $\alpha$ integrated out.
- **Unit-specific probability:** condition on $\alpha$, adding $S_t^{-1}\alpha$ to the mean and dropping the $\sigma_\alpha^2$ term.
- **Ordered:** $P(y_{it} = j) = \Phi((\kappa_j - m_{it})/s_{it}) - \Phi((\kappa_{j-1} - m_{it})/s_{it})$.

### Impacts (marginal effects on probabilities)

For covariate $r$ in period $t$:

$$
\frac{\partial P(y_{it} = 1)}{\partial x_{jt,r}} = \frac{\phi(m_{it}/s_{it})}{s_{it}}\,\big[S_t^{-1}(\beta_r I + \theta_r W_t)\big]_{ij}.
$$

- **Average direct impact:** the mean over $(i, t)$ of the $i = j$ terms. **Average total impact:** the mean over $(i, t)$ of the row sums. **Indirect = total − direct.** Each is reported with its posterior mean, standard deviation and 95% credible interval, both averaged over periods and (optionally) per period.
- **Ordered:** category-specific impacts, replacing $\phi(\cdot)/s$ with $[\phi((\kappa_{j-1} - m)/s) - \phi((\kappa_j - m)/s)]/s$. The category impacts sum to zero across $j$, and the tests check this.
- **SEM:** no spillovers, so only direct effects are reported (with the SEM-specific $s_{it}$).
- **Dynamic models:** short-run (within-period) impacts in closed form. Longer-run responses come from posterior predictive simulation of paths (Phase 5), because nonlinear dynamic spatial models have no closed-form long-run multiplier on probabilities.
- **Difference from `spatialprobit`:** it uses $\phi(m_i)$ without dividing by $s_i$, i.e. it ignores the reduced-form heteroskedasticity created by $S^{-1}$. We implement the formula above and offer the LeSage–Pace convention as an option (`impacts(..., scale = "lesage_pace")`) for comparison with published results (see Lacombe & LeSage 2013/2018 on interpretation).
- **Dummy regressors:** discrete-change impacts are a Phase 2 stretch goal.

## 1.6 Explicitly out of scope for v1.0

The following are left for v2 or for other packages:

- Frequentist estimators: Mendell–Elston approximate ML (`ProbitSpatial`), GHK simulated ML, partial ML (`pmlsp`), GMM (`spldv`, `McSpatial`).
- Fixed-effects probit.
- Bivariate, multivariate and multinomial outcomes (the BEK system, Gao's MSOP).
- Spatial Tobit.
- Convex combinations of several $W$ (Debarsy & LeSage; `spmixW`).
- Panels with missing covariates. Missing **outcomes** are in scope; see Phase 6.
