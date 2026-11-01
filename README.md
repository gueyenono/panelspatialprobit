# panelspatialprobit

Bayesian estimation of **panel spatial probit** and **panel spatial ordered
probit** models in R.

For unit $i$ and period $t$ the latent outcome follows

$$z_t = \rho W z_t + X_t\beta + W X_t\theta + \alpha + u_t,\qquad u_t = \lambda M u_t + \varepsilon_t,\qquad \varepsilon_t \sim N(0, I_N),$$

with $y_{it} = 1\{z_{it} > 0\}$ (binary) or $y_{it} = j$ when
$\kappa_{j-1} < z_{it} \le \kappa_j$ (ordered). Supported options:

- **Spatial specification:** `model = "sar"`, `"sdm"`, `"slx"`, `"sem"`,
  `"sdem"`, `"sarar"` or `"none"`, with separate error weights `M` if needed.
- **Spatial weights:** one matrix, or one per period (time-varying weights).
- **Heterogeneity:** pooled, random or spatially autocorrelated random effects
  (`effects =`). Correlated random effects use Mundlak means (`cre = ~ x1`) or
  Chamberlain period values (`cre = list(type = "chamberlain", vars = ~ x1)`),
  and `time_effects = TRUE` adds period dummies.
- **Impacts:** average direct, indirect and total effects on probabilities
  (per category for ordered models), population-averaged or conditional on the
  random effects, overall or per period.
- **Dynamics:** lagged outcomes and their spatial lag, with Wooldridge
  initial conditions (`dynamic = psp_dynamic()`), or latent AR(1) dynamics
  (`psp_dynamic(latent = TRUE)`), and dynamic response paths
  (`psp_response_path()`).
- **Errors:** normal (probit) or Student-t (`errors = psp_t(df)`, robit).
- **Missing outcomes** are imputed during sampling.
- **Fit and checking:** predicted and counterfactual probabilities
  (`predict(newdata = )`), hit rate and AUC (`psp_classification()`),
  posterior predictive checks (`psp_ppc()`), and an experimental simulated
  likelihood with DIC (`psp_loglik()`, `psp_dic()`).
- **Sampler:** written in Rust. Chains run in parallel threads (`cores`) and
  are reproducible from one seed, whatever the number of cores; impacts,
  predictions and the other post-processing steps run in parallel over
  `future::plan()`. Convergence diagnostics come from the `posterior` package.

The package is in early development (see `plan/` for the roadmap). Rows
must form a balanced panel (outcomes may be missing).

## Installation

The package is not on CRAN yet. Installing it from source needs a Rust
toolchain (Cargo and rustc >= 1.84; see <https://rustup.rs>):

```r
# install.packages("remotes")
remotes::install_github("gueyenono/panelspatialprobit")
```

## Example

```r
library(panelspatialprobit)

W <- psp_weights(rook_weights(10, 10), style = "W")
d <- sim_psprobit(W, T = 5, beta = c(0, 1, -1), rho = 0.6, sigma_alpha = 0.8, seed = 42)

fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                control = psp_control(seed = 1))
summary(fit)          # posterior summary, R-hat, ESS, impacts
psp_impacts(fit)      # direct / indirect / total effects on P(y = 1)
predict(fit)          # posterior mean probabilities, in the row order of `d`

# ordered outcome
d2 <- sim_psprobit(W, T = 5, beta = c(0, 1, -1), rho = 0.5, cutpoints = c(0.8, 1.6), seed = 7)
fit2 <- psoprobit(y ~ x1 + x2, data = d2, index = c("id", "time"), W = W)
```

## How the sampler works

Each model reduces, given its spatial parameters, to a Gaussian kernel
$\lVert Az - b\rVert^2$ for the latent outcomes with a sparse operator
$A = (I - \lambda M)(I - \rho W)$, so one sampler handles every variant.
Each iteration draws, in order:

1. the latent outcomes, with exact, warm-started truncated-normal Gibbs sweeps
   that update residuals in place (no matrix solves);
2. $\rho$, by griddy Gibbs with $\beta$ integrated out;
3. $\beta$, then $\lambda$ (griddy Gibbs), then the random effects and their
   variance (half-Cauchy prior);
4. for ordered models, the cut-points: an Albert–Chib step plus a joint move of
   cut-points and latent values;
5. a scale move that rescales all of the above jointly (Liu & Sabatti, 2000),
   which improves mixing.

`plan/02_estimation_algorithms.md` gives the details.

## Validation

- **Reference sampler:** the Rust engine is tested against a dense pure-R sampler
  that uses different algorithms (`tests/testthat/test-reference.R`).
- **Simulation-based calibration:** every supported model variant passes SBC
  (`validation/sbc/run_sbc.R`).
- **Monte Carlo study:** bias, RMSE and coverage at fixed parameters, including
  state dependence versus heterogeneity in dynamic models
  (`validation/montecarlo/run_mc.R`).
- **Impacts:** checked against finite differences of the exact probabilities.
- **Cut-point moves:** a stationarity test on a toy ordered model
  (`tests/testthat/test-cutmove.R`).
- **Cross-checks:** agreement with `spatialprobit` on the Katrina data and with
  `MCMCpack::MCMCoprobit` for non-spatial ordered probit
  (`validation/crosscheck/`).

The vignettes `vignette("getting-started", package = "panelspatialprobit")` and
`vignette("ordered-outcomes", package = "panelspatialprobit")` walk through
complete analyses.
