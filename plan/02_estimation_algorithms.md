# 2. Estimation algorithms

Estimation is Bayesian MCMC with data augmentation (Albert & Chib 1993; LeSage 2000), extended to panels as in Kakamu & Wago (2005) and Baltagi, Egger & Kesina (2016). This file specifies each Gibbs block, the computational strategy, and the changes relative to `spatialprobit`.

## 2.1 Sampler state and block order

State: $(z, \rho, \beta, \lambda, \alpha, \sigma_\alpha^2, \xi, \kappa, v)$, where $\xi$ is the auxiliary variable of the half-Cauchy prior.

One iteration, in this order:

1. $z \mid \cdot$ — one or more **warm-started** Gibbs sweeps over the latent vector (§2.3).
2. $\rho \mid z, \alpha, \lambda$ **with $\beta$ integrated out** — griddy Gibbs (§2.4).
3. $\beta \mid z, \rho, \lambda, \alpha$ — Gaussian.
4. $\lambda \mid z, \rho, \beta, \alpha$ — griddy Gibbs (SEM/SARAR only).
5. $\alpha \mid z, \rho, \lambda, \beta, \sigma_\alpha^2$ — Gaussian (independent across $i$, or sparse joint under SEM/SARAR).
6. $\sigma_\alpha^2, \xi \mid \alpha$ — inverse gamma.
7. $\kappa, z \mid \cdot$ — cut-point move (ordered only, §2.7).
8. $v \mid \cdot$ — inverse gamma ($t$ errors only).

**Why this order matters.** Step 2 marginalises $\beta$, so $\beta$ must be redrawn (step 3) before any block conditions on it. This is the partially collapsed Gibbs rule of van Dyk & Park (2008). Breaking the rule can change the stationary distribution (see §2.10).

## 2.2 The core abstraction

Conditional on $(\rho, \lambda)$, every model in `01` has the latent kernel

$$
p(z \mid \cdot) \propto \exp\!\Big(-\tfrac12 (Az - b)' V^{-1} (Az - b)\Big)\ \text{restricted to the truncation box defined by } y,
$$

where $A = \operatorname{blockdiag}_t(R_t S_t)$ and $b = R(X\beta + \Delta\alpha)$. The Jacobian of the latent density is $\log|A| = \sum_t \big(\log|I - \rho W_t| + \log|I - \lambda M_t|\big)$.

The engine therefore needs only four operations:

- a sparse $A$, rebuilt when $\rho$ or $\lambda$ changes (only the values change, not the sparsity pattern);
- a residual vector $r = Az - b$;
- per-coordinate access to the columns of $A$;
- the log-determinant grids.

Adding a model variant means saying how to build $A$ and $b$, not writing a new sampler.

## 2.3 Latent block: precision-free, residual-update Gibbs sweep

For coordinate $k$ with column $a_k$ of $A$ (non-zeros: itself plus the units that list $k$ as a neighbour):

```
h_k   = a_k' V^{-1} a_k                    # conditional precision (cache; changes only with rho/lambda/v)
m_k   = z_k - (a_k' V^{-1} r) / h_k        # conditional mean
z_k'  ~ TN(m_k, 1/h_k; lower_k, upper_k)   # bounds from y (and kappa for ordered)
r    += a_k * (z_k' - z_k)                 # O(nnz(a_k)) update
```

- **Cost:** $O(\mathrm{nnz}(A))$ per sweep. The precision matrix $H = A'V^{-1}A$ and the mean $\mu = A^{-1}b$ are **never formed or solved for**. This removes the per-iteration sparse QR that `spatialprobit` performs.
- **Exactness:** each coordinate update is an exact draw from its full conditional, starting from the current state. One sweep per iteration is therefore a valid Gibbs kernel; `control$z_sweeps > 1` only improves mixing per iteration.
- **Missing outcomes** (Phase 6): set bounds to $(-\infty, \infty)$. This is the Bayesian imputation of $z$, and it keeps the full spatial operator intact for unbalanced panels.
- **SARAR:** $A = RS$ has second-order neighbours. Form the product once per $(\rho, \lambda)$ update with Eigen's sparse–sparse product; the cost is about $\mathrm{nnz}(W) \times$ average degree.
- **Univariate truncated normal:** an exact algorithm that is stable in the tails, using Robert (1995) exponential rejection for one-sided tails and Chopin (2011) or Botev & L'Ecuyer (2017) for two-sided intervals. Means 40 sd from the bound must not under- or overflow. It uses R's RNG (`unif_rand`, `norm_rand`) so `set.seed()` reproduces runs.
- **Sweep order:** sequential in v1. Checkerboard or colouring orders are possible later for parallelism.

## 2.4 $\rho$: collapsed griddy Gibbs

For SAR, SDM and SARAR (given $\lambda$), integrate $\beta \sim N(c, \Sigma_\beta)$ out. With $R$-transformed data $\tilde X = RX$, set $a_0 = R(z - \Delta\alpha) - \tilde X c$ and $a_1 = RWz$, stacked over $t$. Then

$$
\log p(\rho \mid z, \alpha, \lambda) = \sum_t \ell_t(\rho) - \tfrac12\big(q_{00} - 2\rho\, q_{01} + \rho^2 q_{11}\big) + \log\pi(\rho) + \text{const},
$$

where

$$
q_{ab} = a_a'V^{-1}a_b - (\tilde X'V^{-1}a_a)'\,B\,(\tilde X'V^{-1}a_b), \qquad B = (\tilde X'V^{-1}\tilde X + \Sigma_\beta^{-1})^{-1}
$$

and $\ell_t(\rho) = \log|I - \rho W_t|$.

- This is exact for $\sigma_\varepsilon = 1$ (Woodbury identity; $|I + \tilde X\Sigma_\beta\tilde X'|$ does not depend on $\rho$).
- **Cost:** the $q_{ab}$ are three scalars from $O(nk)$ work, and evaluating the grid is $O(\text{grid size})$.
- **Two-stage grid.** With $n = 10^5$ the posterior sd of $\rho$ can be about 0.003, which is finer than a fixed 0.001 grid can resolve well. So: evaluate a coarse grid (about 200 points); refine about 500 points over mode ± 8 sd; get $\ell(\rho)$ by cubic-spline interpolation of the precomputed log-determinant grid; draw by inverse CDF with linear interpolation inside cells, giving a continuous draw rather than a grid point.
- **Fallback:** an adaptive random-walk Metropolis step (Kakamu & Wago; BEK), for cases where a log-determinant grid is impractical (huge $N$ with many distinct $W_t$). Tuning targets 0.4–0.6 acceptance and **adapts only during warm-up**; adaptation is frozen afterwards so the chain is a valid MCMC.

## 2.5 $\lambda$: griddy Gibbs, not collapsed

Given $\beta$, let $e = Sz - X\beta - \Delta\alpha$. Then $R e = e - \lambda Me$, so $Q(\lambda) = e'V^{-1}e - 2\lambda\, e'V^{-1}Me + \lambda^2 (Me)'V^{-1}(Me)$. This is quadratic and uses the same grid machinery with $\sum_t \log|I - \lambda M_t|$.

It is not collapsed over $\beta$, because the design $RX$ depends on $\lambda$. After each $\lambda$ draw, refresh $\tilde X = RX$ and $\tilde X'V^{-1}\tilde X$ ($k \times k$).

## 2.6 Random-effects blocks

- **$\alpha$ under SAR, SDM or `none`** (no $R$): with $e_{it} = (S_t z_t - X_t\beta)_i$,
  $\alpha_i \mid \cdot \sim N\!\left(\frac{\sum_t e_{it}/v_{it}}{\sum_t 1/v_{it} + \sigma_\alpha^{-2}},\ \frac{1}{\sum_t 1/v_{it} + \sigma_\alpha^{-2}}\right)$. The $\alpha_i$ are independent across $i$, so the cost is $O(n)$ (Kakamu & Wago 2005, §2.2.2).
- **$\alpha$ under SEM or SARAR:** the precision is $P = \sum_t R_t'V_t^{-1}R_t + \sigma_\alpha^{-2}I_N$, a sparse $N \times N$ matrix. Use a sparse Cholesky (Eigen `SimplicialLLT`): the symbolic factorisation is done once because the pattern is fixed, and the numeric factorisation is redone per $\lambda$.
- **`spatial_random`:** the prior precision becomes $\sigma_\alpha^{-2}(I - \phi W)'(I - \phi W)$, which is also sparse, with $\phi$ drawn by griddy Gibbs using $\log|I - \phi W|$.
- **$\sigma_\alpha^2$ (half-Cauchy$(0, A)$):** $\sigma_\alpha^2 \mid \alpha, \xi \sim IG\big(\tfrac{N+1}{2}, \tfrac12\sum_i\alpha_i^2 + 1/\xi\big)$ and $\xi \mid \sigma_\alpha^2 \sim IG(1, 1/\sigma_\alpha^2 + 1/A^2)$ (Wand et al. 2011).
- **Mixing risk.** $\alpha$, $z$ and the intercept are strongly coupled when $T$ is small or $\sigma_\alpha$ is large. Phase 1 ships the simple blocks and **measures ESS per second**. Phase 6 adds whichever of the following the benchmarks justify:
  1. a joint $(\beta, \alpha)$ block, whose precision has an arrowhead structure;
  2. centred versus non-centred parameterisations of $\alpha$ (Papaspiliopoulos, Roberts & Sköld 2007);
  3. marginal or parameter-expanded data augmentation for the probit scale (Liu & Wu 1999; Imai & van Dyk 2005).

### Mixing experiments (2026-10-28)

Item 3 shipped in Phase 1 as the scale move. Three further moves were then implemented and measured. Each left the posterior unchanged: the C++ vs reference-sampler comparisons passed with every move. **All three were reverted** because none improved effective sample size. The test cases had 4,000 draws × 2 chains, using bulk ESS:
- static SAR with random effects, N = 100, T = 5, balanced outcomes;
- the same with rare outcomes (about 5% ones), binary and ordered;
- a latent-AR(1) case on a 4 × 4 grid with T = 4 and $\gamma = 1.5$, with near-separated outcomes.

| Move | What it did | Result |
|---|---|---|
| Latent-AR time-profile move | Generalised Gibbs on $z_t \to g^{t-1}z_t$, $\gamma \to g\gamma$ (sign-preserving, binary only) | ESS of $\gamma$ unchanged within noise (e.g. 41–64 vs 58–137 at $\gamma = 0.8$; 9 vs 6 in the worst case); $\gamma$–$\rho$ correlation unchanged |
| Scale move including $\sigma_\alpha$ | Rescale $(z, \beta, \alpha, \kappa, \sigma_\alpha)$ jointly; $\alpha/\sigma_\alpha$ invariant | No change (rare case: $\sigma_\alpha$ ESS 73 → 56, intercept 38 → 35). The unit error variance pins the overall latent scale, so $g$ stays within about $1 \pm 1/\sqrt{2n}$ |
| Interweaving for $\sigma_\alpha$ (ASIS; Yu & Meng 2011) | Redraw $\sigma_\alpha$ given $\alpha/\sigma_\alpha$ and $z$ (slice sampler) | No change (rare case: $\sigma_\alpha$ ESS 73 → 33, intercept 38 → 34; latent worst case 5 → 5) |

Four sweeps of the latent block per iteration did not help the latent-AR case either.

**Not tried, and not worth trying: a joint update of each unit's latent path** $z_{i,1..T}$ (2026-10-29). Any such move, however exact (for example exact HMC for the truncated normal path), draws from the same conditional $p(z \mid \theta, y)$ that repeated sweeps converge to. Its best case is therefore an exact draw of $z$ given the parameters, and four sweeps per iteration already approximate that without improving ESS. The slow mode is the dependence between $z$ and the parameters (diagnosis below), which no move on $z$ alone can remove.

**Diagnosis.** In the rare-outcome case (4.6% ones), the posterior correlations are intercept–$\rho$ 0.84, intercept–$\sigma_\alpha$ −0.69 and $\sigma_\alpha$–$x_2$ −0.57. The autocorrelation times are 160 for the intercept and 174 for $\sigma_\alpha$, against 39 for $x_1$. Even the best-identified combination, the mean latent level intercept / $(1-\rho)$, has ESS 135 out of 16,000. The slow mode is therefore not one parameter pair that a reparameterisation could separate. It is the known slowness of data augmentation with imbalanced outcomes (Johndrow, Smith, Pillai & Dunson 2019): given the latent values of the many zeros, the intercept and spatial level are tightly pinned, while given $y$ they are not. Strong latent AR dynamics and near-separation amplify it.

The remedies in that literature are calibrated data augmentation, or Metropolis/HMC on the marginal posterior. Both need the likelihood $p(y \mid \theta)$, an $NT$-dimensional rectangle probability that is available here only as a noisy GHK estimate (`psp_loglik()`). There is no cheap exact fix within this sampler. The package therefore documents that rare outcomes, near-separation and strong latent dynamics need longer chains; the R-hat / ESS warning flags these fits.

## 2.7 Cut-points (ordered model)

- **Baseline: Albert–Chib Gibbs.** $\kappa_j \sim U\big(\max\{z: y = j\} \vee \kappa_{j-1},\ \min\{z: y = j+1\} \wedge \kappa_{j+1}\big)$, times the prior density if it is proper. This is simple, but it is known to mix badly when $n$ is large, because the allowed interval shrinks like $O(1/n)$. That is exactly the panel setting.
- **Cowles (1996) collapsed MH is not directly available.** It needs $P(y \mid \kappa, \beta, \dots)$ with $z$ integrated out, which in a spatial model is an $n$-dimensional normal rectangle probability. A GHK or Mendell–Elston approximation would break exactness.
- **Proposed default: a joint transformation move on $(\kappa_j, z)$** (generalised Gibbs; Liu & Sabatti 2000). Propose a new $\kappa_j'$ and map the latent axis piecewise-linearly so that every $z$ in categories $j$ and $j+1$ stays in its category:
  - **Interior $j$:** random walk on $\operatorname{logit}$ of the relative position of $\kappa_j$ in $(\kappa_{j-1}, \kappa_{j+1})$. Category $j$ is rescaled onto $(\kappa_{j-1}, \kappa_j']$ and category $j+1$ onto $(\kappa_j', \kappa_{j+1}]$.
  - **Top cut-point:** random walk on $\log(\kappa_{J-1} - \kappa_{J-2})$. Category $J-1$ is rescaled and category $J$ is translated.
  - **MH ratio:** $\dfrac{p(z' \mid \cdot)}{p(z \mid \cdot)} \times |\text{Jacobian}| \times \dfrac{\pi(\kappa')}{\pi(\kappa)}$. The Jacobian is the product of per-category scale factors raised to the category counts, times the Jacobian of the $\kappa$ reparameterisation. $p(z \mid \cdot)$ needs only $\|r\|^2_{V^{-1}}$, because the log-determinant does not change, so one evaluation costs $O(\mathrm{nnz}(A))$.
- **Decision rule.** Implement both. Run SBC on both, which must pass. Make the move with the higher ESS per second the default. Keep Albert–Chib behind `control$cutpoint_sampler = "albert_chib"`.

## 2.8 Other blocks

- **$\beta$:** $N\big(B(\tilde X'V^{-1}\tilde y + \Sigma_\beta^{-1}c),\ B\big)$ with $\tilde y = R(Sz - \Delta\alpha)$. Under SAR with $V = I$, $X'X$ is constant, so its Cholesky factor is precomputed.
- **$v_{it}$** ($t$ errors): $IG\big(\tfrac{r+1}{2}, \tfrac{r + \varepsilon_{it}^2}{2}\big)$. After each update, refresh the $h_k$ cache and the weighted cross-products.

## 2.9 Log-determinants

- **Time-invariant $W$:** compute $\ell(\rho) = \log|I_N - \rho W|$ once on a grid and multiply by $T$. Method by size (`control$lndet = "auto"`):
  - $N \le$ about 3,000: eigenvalues, exact. $\sum_k \log(1 - \rho\,\omega_k)$, using the real part of the complex log so asymmetric $W$ is handled, and symmetrising when $W$ is similar to a symmetric matrix.
  - Larger $N$: sparse LU or Cholesky on grid points (Pace & Barry 1997), then a spline.
  - Very large $N$: Chebyshev (Pace & LeSage 2004) or Monte Carlo (Barry & Pace 1999) approximations.
- **Time-varying $W_t$:** $\sum_t \ell_t(\rho)$. Grids are cached by hashing the matrices, so identical periods are computed once; approximate methods are used by default when $N \times T$ is large.
- **Reuse.** The `lndet` argument accepts a precomputed grid (as `spmixW` does). `spatialreg::do_ldet` is used only as a test oracle (Suggests), not as a runtime dependency, because its interface is semi-internal.

## 2.10 Lessons from the `spatialprobit` source (1.0.4)

These are observations from reading the code. Each will be **confirmed empirically in Phase 1 by running SBC on `spatialprobit::sarprobit`** before we rely on it in documentation or a paper.

| # | Observation | Where | Our design |
|---|---|---|---|
| 1 | The $\rho$ draw uses the linear-SAR concentrated form, $-\tfrac{n-k}{2}\log Q(\rho)$, which integrates out $\sigma^2$. In a probit $\sigma = 1$ is fixed, so the exact kernel is $\exp(-\tfrac12 Q(\rho))$. | `draw_rho`, `R/SpatialProbit-MCMC.R` (lines 63–90, 355) | Exact collapsed conditional (§2.4). |
| 2 | Block order z → β → ρ(β integrated) → z(using the stale β) breaks the partially collapsed Gibbs rule. | `R/SpatialProbit-MCMC.R` (lines 302–368) | Order z → ρ → β (§2.1). |
| 3 | With `m > 1` the truncated-normal inner Gibbs chain is restarted each iteration (no `start.value`), so the z-update is only approximately a full-conditional draw. | `rtmvnorm.sparseMatrix` calls (lines 314–320) | Warm-started exact sweeps (§2.3). |
| 4 | A sparse QR of $S$ and a solve every iteration. | lines 359–368 | Precision-free residual updates. |
| 5 | Ordered probit: the prior-mean variable `c` is overwritten by the posterior mean inside the loop (harmless only under a diffuse prior); `phidraw` indexes `B` without a comma; the `2:(J-1)` loop runs backwards when `J = 2`; there are no impacts. | `R/sarorderedprobit.R` (lines 217, 242, 323) | Fixed by design and tested. |
| 6 | Impacts use $\phi(m_i)$ without the reduced-form sd $s_i$. | lines 398–415 | Formula in `01` §1.5, with a LeSage–Pace compatibility option. |
| 7 | `logLik` uses `diag(S %*% t(S))`, but $\operatorname{Var}(z) = S^{-1}S^{-\top}$; it is also a product of marginal probabilities (a pseudo-likelihood). | lines 809–829 | Approximate observed-data likelihood (§2.12), clearly labelled. |
| 8 | No panel structure, no random effects, no multiple chains or convergence diagnostics. | — | Core scope. |

**SBC result** (2026-10-27; `validation/sbc/sbc_spatialprobit.R`, results in `validation/sbc/results/`). The design was 8×8 rook grid, T = 1, pooled SAR probit, β ~ N(0, I), ρ ~ Beta(2, 2), 500 replications, start values β = 0 and ρ = 0.

- **Failures:** `spatialprobit::sar_probit_mcmc` returned NaN draws in **264 of 500** replications.
- **ρ:** among the 236 successful replications, ρ is clearly miscalibrated (χ² p < 0.0001, mean rank 0.61, a pile-up in the top bin: the posterior understates ρ). The intercept is borderline (p = 0.037).
- **Caveat:** keeping only the successful runs could itself distort the ranks, so this is strong evidence of a problem but not a clean measurement of its size.
- **Control:** `panelspatialprobit` on the identical design had 0 failures and all p > 0.08.

These results support observations 1–3 above. Which of them causes the bias has not been isolated.

### SBC results after the phase 2-6 work (2026-10-27)

The engine was rewritten around the general operator $A = RS$ (§2.2), so every variant was re-run: 21 variants, 500 replications each, $N = 16$, $T = 4$. The variants are:
- {sar, sdm, slx, sem, sdem, sarar, none} with random effects;
- pooled sar and sem;
- Mundlak and Chamberlain CRE;
- ordered $J = 3, 4$ with both cut-point samplers, and sem with $J = 3$;
- dynamic with Wooldridge initial conditions;
- 20% missing outcomes (binary and ordered).

There were no sampler failures. Of the 109 rank-uniformity tests, 6 have $p < 0.05$ (about 5 expected under uniformity); the smallest is $p = 0.011$. Results are in `validation/sbc/results/`. The dynamic variant comes from a separate re-run (`sbc_2026-10-27_sar_dyn_bin.txt`): in the main run its data generator errored on panels with a constant outcome, which voided most of its replications.

**Later variants** (same design, 500 replications each):
- **Time-varying weights.** A different fifth of the rook links is removed in each period. `sar_tv_re_bin` passes. In the first run of `sarar_tv_re_ord3`, 4 of 7 tests had $p < 0.05$, with no systematic histogram shape. An independent replication (`SBC_SEED_OFFSET=1`) passed (min $p = 0.064$), and so did the pooled 1,000 replications (min $p = 0.064$). The non-time-varying `sarar_re_ord3` was run alongside and passed, apart from one $p = 0.019$.
- **Spatial random effects.** `sar_sre_bin` and `sem_sre_ord3` pass (min $p = 0.089$).
- **Student-t errors** (df = 4). `sar_t_re_bin` and `sem_t_re_ord3` pass. `sarar_t_re_ord3` gave $\lambda$ $p = 0.024$ and $0.030$ in two runs (pooled $p = 0.002$). With chains twice as long, $\lambda$ passes ($p = 0.775$), and $\rho$ is at $p = 0.027$, as in normal-error SARAR ($p = 0.019$). This is slow $\rho/\lambda$ mixing in SARAR with $M = W$ at the default length, not a t-specific error.
- **Latent AR(1)** ($\gamma \sim N(0, 1)$). `sem_lat_re_ord3` passes. `sar_lat_re_bin` fails at the default length ($\gamma$ $p = 0.001$). The failure is concentrated in the replications with max $\hat R \ge 1.02$, with U-shaped ranks. With chains three times as long it passes (min $p = 0.090$), so $\gamma$ mixes slowly against $\rho$ when drawn given $z$. A joint $(\gamma, z)$ move is a candidate improvement.

## 2.11 Impacts computation

Impacts are computed **after sampling** on a thinned subset of draws (`impacts_draws`, default 1,000). Population-averaged impacts need only $(\beta, \theta, \rho, \lambda, \sigma_\alpha)$. Unit-specific impacts need $\alpha$, so the sampler stores $\alpha$ draws for that subset only (for example $N = 5{,}000$ × 1,000 draws = 40 MB).

Per draw (one sparse LU of $S$ per draw for time-invariant $W$, reused across $t$):

- $m_t = S^{-1}X_t\beta$: one triangular solve pair per $t$.
- $s_{it}$ requires $\operatorname{diag}(S^{-1}S^{-\top})$ (and $A^{-1}VA^{-\top}$ for SEM/SARAR):
  - exact (dense) when $N \le$ `impacts_exact_max` (default 500; lowered from 2,000 after the Phase 1 benchmark, where the stochastic path was 8x faster at N = 1,000 with error about 5% of the posterior sd);
  - otherwise a stochastic diagonal estimator (Bekas, Kokiopoulou & Saad 2007) with Rademacher probes, where each probe costs two sparse solves. The number of probes is configurable, and the reported MC error is tested against the exact value on medium $N$.
- **Direct:** needs $\operatorname{diag}(S^{-1})$ (and $\operatorname{diag}(S^{-1}W)$ for SDM), computed the same way (exact or stochastic; a Neumann/power series is an alternative for $|\rho| < 1$).
- **Total:** $\frac1N\, g' S^{-1}(\beta_r\iota + \theta_r W\iota)$, one solve, where $g_i = \phi(m_i/s_i)/s_i$.

The output is an object of class `psp_impacts` holding draws of the average direct, indirect and total impacts per covariate (per category for ordered models, optionally per period), with summary and plot methods.

## 2.12 Diagnostics and model assessment

- **Convergence:** 4 chains by default. Rank-normalised split-$\hat R$ and bulk/tail ESS (Vehtari et al. 2021) via `posterior`. The fit warns when $\hat R > 1.01$ or ESS < 400 for any parameter. The fit object stores acceptance rates for MH and transformation moves, and timing.
- **Posterior predictive checks** (`pp_check`-style, Phase 6): simulate $y^{rep}$ from the draws and compare:
  - prevalence by period;
  - the transition matrix $y_{t-1} \to y_t$;
  - per-period join-count or Moran statistics of $y$;
  - category shares for ordered models.
- **Approximate likelihood for model comparison** (Phase 6, experimental). The observed-data likelihood is an $n$-dimensional orthant probability. Given $\alpha$ the periods are independent, so it factorises into $T$ $N$-dimensional problems. Approximate these with Mendell–Elston (Martinetti & Geniaux 2017; `ProbitSpatial`'s C++ code is GPL ≥ 2) or GHK (Pace & LeSage 2016). This is used for a DIC/WAIC-type comparison and documented as approximate.
- **In-sample classification:** hit rates and AUC.

## 2.13 Performance targets

These are to be re-baselined after the Phase 1 benchmark.

- **Per-iteration cost:** $O(\mathrm{nnz}(W)\,T + nk)$; no $n \times n$ object is ever formed.
- **Memory:** draws are stored for parameters only; $z$ is kept as a running posterior mean unless `save_latent = TRUE`.
- **Targets on a laptop core:** SAR-RE binary with $N = 1{,}000$, $T = 10$, $k = 5$ runs 10,000 iterations in about a minute or less, and $N = 10{,}000$, $T = 20$ runs in about 30 minutes or less.
- **Parallelism:** chains run in parallel on Rust threads (`psp_control(cores = )`, every OS), each with its own L'Ecuyer-CMRG stream from `.psp_streams()`, so the draws do not depend on `cores`. Within a chain the sweep is sequential (a period-parallel sweep is R6 in `06`). Post-processing (impacts, predict, PPC, response paths, GHK) spreads blocks of draws over the user's `future::plan()`, each Rust call single-threaded (`06`, R-D11). *Until 2026-10-31 (C++ engine): chains ran in parallel at the R level with `parallel::mcmapply` (forking, so not on Windows).*
