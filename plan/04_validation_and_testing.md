# 4. Validation and testing

Validation has three layers:

1. **Is the code doing what we think?** Unit tests, which ship with the package.
2. **Does the sampler target the right posterior?** Simulation-based calibration (SBC), the release gate.
3. **Is the estimator useful?** Monte Carlo studies, cross-package agreement, empirical replications.

Layer 1 runs on CRAN and CI. Layers 2 and 3 live in `validation/` (excluded from the build) and run on demand or nightly.

## 4.1 Unit tests (`tests/testthat/`, CRAN-safe, under 2 minutes in total)

| Area | Test |
|---|---|
| Truncated normal | Moments versus analytic values for one-sided, two-sided and extreme-tail cases (bounds up to 40 sd from the mean). KS test with a fixed seed. No NaN or Inf. |
| Latent sweep | (a) **Invariant:** after a sweep, `r == A z - b` to 1e-10. (b) **Stationarity:** for $n = 5$, a long run of sweeps matches the TMVN mean and covariance from `tmvtnorm` or importance sampling, within MC error. (c) SAR, SEM, SARAR and $V \neq I$ variants of (b). |
| Log-determinants | Each method versus `determinant()` on small $W$ (ring, rook grid, kNN; symmetric and asymmetric); approximations within stated tolerance; spline-interpolation error bound. |
| Collapsed $\rho$ kernel | $q_{ab}$ and the kernel versus brute-force dense marginalisation of $\beta$ on small problems. |
| Griddy draw | Draws from a known 1-D density (Beta, truncated normal) pass a KS test; continuous within-cell interpolation. |
| RE blocks | $\alpha$ conditional mean and variance versus dense Gaussian conditioning; the SEM sparse-Cholesky path versus dense. |
| Cut-point moves | A detailed-balance check on a toy ordered model, with the target computed by enumeration and quadrature for $n = 3$. Every proposal keeps $z$ inside its category. |
| Impacts | **Finite differences:** for small $N$, numerically differentiate the exact $P(y_{it} = 1)$ (dense) with respect to each $x_{jt,r}$ and compare with our direct, indirect and total impacts. The stochastic diagonal estimator is within tolerance of the exact value. Ordered category impacts sum to 0. |
| Data pipeline | Permuting the rows of the input data gives identical internal $(y, X)$ and identical results under the same seed. Mismatched $W$ IDs, duplicate (id, time) pairs, unbalanced panels, empty ordered categories and $J = 2$ in `psoprobit()` all raise informative errors. |
| Reproducibility | Same seed gives identical draws; different `cores` give the same per-chain results (L'Ecuyer streams). |
| Methods | Snapshot tests for `print`/`summary`; `coef`/`vcov`/`confint` shapes; `as_draws` round-trip. |
| Simulator | `sim_psprobit()` hits target prevalences and cut-point shares; parameters are stored as attributes. |

## 4.2 Sampler correctness: SBC (`validation/sbc/`)

- **Method:** simulation-based calibration (Talts et al. 2018; Modrák et al. 2025, *Bayesian Analysis*, doi:10.1214/23-BA1404). Draw $\theta \sim$ prior, simulate $y \mid \theta$ with `sim_psprobit()`, fit, and compute rank statistics.
  - **Test quantities:** every scalar parameter plus data-dependent quantities, such as the joint log-density of $(z, y)$ and the average impacts, which Modrák et al. show catch errors that parameter ranks miss.
  - **Displays:** ECDF-difference plots and rank histograms.
- **Design:** $N = 25$ (5×5 rook grid, row-standardised), $T = 4$, 500 replications per variant (1,000 for the core), and about 200 roughly independent draws per fit after thinning.
- **Variant matrix**, run as each phase lands:
  - specification: {`none`, `sar`, `sdm`, `sem`, `sarar`};
  - effects: {`pooled`, `random`, `cre`};
  - outcome: {binary, ordered $J = 3, 5$};
  - dynamics: {static, `lag_y` + Wooldridge initial conditions};
  - samplers: both cut-point samplers.
- **Complement during development:** the Geweke (2004) "getting it right" joint-distribution test, which is cheaper and finds bugs earlier.
- **Experiment on existing code:** run the same SBC on `spatialprobit::sarprobit` and `sarorderedprobit` with $T = 1$, using proper priors set through their `prior` argument. This tests the hypotheses in `02` §2.10 and either motivates or retracts the "why a new sampler" argument.
- **Gate:** no variant is documented as supported until its SBC passes, meaning no systematic ECDF departure beyond the 95% simultaneous band.
- **Tooling:** the `SBC` package from GitHub, or our own rank computation (about 50 lines) to avoid a non-CRAN dependency.

## 4.3 Frequentist Monte Carlo (`validation/montecarlo/`; paper-grade)

- **Grid:**
  - $N \in \{100, 400, 1600\}$, $T \in \{5, 10\}$;
  - $\rho \in \{0.2, 0.5, 0.8\}$ (plus a negative case), $\sigma_\alpha \in \{0, 0.5, 1\}$;
  - $k = 2$ (one continuous and one binary regressor), prevalence ≈ 0.3–0.5;
  - $W$ from rook grid, kNN ($k = 6$) and a real county map;
  - 200 replications.
- **Metrics:** bias, RMSE and 95% credible-interval coverage for $\beta$, $\rho$, $\sigma_\alpha$, $\kappa$ and the average direct and indirect impacts; ESS per second.
- **Comparators** (misspecified, to show the value of the model):
  - non-spatial pooled probit (`glm`);
  - non-spatial RE probit (`lme4::glmer(..., family = binomial("probit"))` or `MCMCglmm`);
  - pooled SAR probit ignoring the panel structure (`spatialprobit::sarprobit` on $I_T \otimes W$);
  - for ordered models, `MCMCpack::MCMCoprobit` and `spatialprobit::sarorderedprobit`.
- **Replications of published designs:**
  - Kakamu & Wago (2005, §3): $N = 50$, $\beta = (1, 1)$, $\alpha_i \sim N(2, 2)$;
  - the univariate restriction of the Baltagi–Egger–Kesina (2016) basic design;
  - the Wang & Kockelman (2009) simulated dynamic spatial ordered probit design ($S = 3$, $\gamma_2 = 2.1$). *Done 2026-10-29: `validation/wang_kockelman/run_wk.R`; see roadmap 5.3.*

## 4.4 Cross-package agreement (`validation/crosscheck/`)

- **$T = 1$, pooled, SAR binary:** our posterior versus `spatialprobit::sarprobit`, on the `Katrina` data shipped in `spatialprobit` and `ProbitSpatial`. ProbitSpatial's approximate-ML point estimates serve as a second reference.
  - Expect agreement within MC error **plus** the systematic differences predicted in `02` §2.10. If SBC shows that `spatialprobit` is miscalibrated, document the gap rather than forcing agreement.
- **$T = 1$, SAR ordered:** our model versus `sarorderedprobit`. *Done 2026-10-29 (`spatialprobit_ordered.R`): means within 0.07 sd; see roadmap 3.6.*
- **$\rho$ fixed at 0:** our pooled binary and ordered models versus `MCMCpack::MCMCprobit` and `MCMCoprobit` (same priors).
- **Log-determinants:** versus `spatialreg::do_ldet` / `jacobianSetup` methods.

## 4.5 Performance (`validation/bench/`)

- Use `bench::mark` and ESS per second for each block. Time per iteration as a function of $N$, $T$ and nnz($W$).
- Record baselines per release in `validation/bench/history.csv` so performance regressions are caught.
- Profile with `perf` or `gprof` only if a target is missed.

## 4.6 Empirical examples (vignettes, paper)

These are candidates to settle in Phase 7.

- **Binary:** the Elhorst, Heijnen, Samarina & Jacobs (2017) inflation-targeting adoption panel (58 countries, 1985–2008). The replication data are on the ZBW Journal Data Archive under CC-BY 4.0, so redistribution as package data is allowed with attribution.
- **Ordered:** US Drought Monitor county categories (D0–D4) by week or month. These are freely available (the producers request credit; check the redistribution terms before shipping them as package data), give a natural ordered spatial panel, and can be built reproducibly from a script in `data-raw/`.
- **Cross-sectional sanity example:** the Katrina business re-opening data, also as an ordered outcome at 3/6/12 months (LeSage et al. 2011).
