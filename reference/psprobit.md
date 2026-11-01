# Bayesian panel spatial probit and ordered probit models

Fits spatial probit (`psprobit()`) and spatial ordered probit
(`psoprobit()`) models to balanced panel data by Markov chain Monte
Carlo. The latent outcome of unit \\i\\ in period \\t\\ follows \$\$z_t
= \rho W z_t + X_t\beta + W X_t\theta + \alpha + u_t, \quad u_t =
\lambda M u_t + \varepsilon_t, \quad \varepsilon_t \sim N(0, I_N),\$\$
with \\y\_{it} = 1\\z\_{it} \> 0\\\\ for binary outcomes and \\y\_{it} =
j\\ if \\\kappa\_{j-1} \< z\_{it} \le \kappa_j\\ for ordered outcomes
(\\\kappa_1 = 0\\).

## Usage

``` r
psprobit(
  formula,
  data,
  index,
  W,
  model = c("sar", "sdm", "slx", "sem", "sdem", "sarar", "none"),
  effects = c("random", "pooled", "spatial_random"),
  cre = NULL,
  durbin = NULL,
  time_effects = FALSE,
  M = NULL,
  dynamic = NULL,
  errors = "normal",
  priors = psp_priors(),
  control = psp_control(),
  impacts = TRUE
)

psoprobit(
  formula,
  data,
  index,
  W,
  model = c("sar", "sdm", "slx", "sem", "sdem", "sarar", "none"),
  effects = c("random", "pooled", "spatial_random"),
  cre = NULL,
  durbin = NULL,
  time_effects = FALSE,
  M = NULL,
  dynamic = NULL,
  errors = "normal",
  priors = psp_priors(),
  control = psp_control(),
  impacts = TRUE
)
```

## Arguments

- formula:

  Model formula with an intercept. For `psprobit()` the response is 0/1,
  logical or a two-level factor; for `psoprobit()` an (ordered) factor
  or integer codes with at least three categories.

- data:

  A data frame with one row per unit and period (balanced panel).
  Outcomes may be missing (`NA`): their latent values are imputed during
  sampling, which keeps the spatial structure intact. Regressors may
  not.

- index:

  Character vector naming the unit and time columns of `data`.

- W:

  Spatial weights for the units (see
  [`psp_weights()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_weights.md));
  its row names identify the units. It is used as given, so
  row-standardise it beforehand (e.g. `psp_weights(W, style = "W")`) if
  desired. For weights that change over time, pass a list with one
  matrix per period, named by period (or in period order); \\W_t\\ then
  replaces \\W\\ in period \\t\\, including in the Durbin terms and the
  spatially lagged outcome of
  [`psp_dynamic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_dynamic.md).

- model:

  Spatial specification:

  - `"sar"`: spatial lag of the latent outcome;

  - `"sdm"`: spatial Durbin, SAR plus \\WX\\;

  - `"slx"`: \\WX\\ only;

  - `"sem"`: spatial error, \\(I - \lambda M)(z_t - X_t\beta - \alpha) =
    \varepsilon_t\\;

  - `"sdem"`: spatial Durbin error, SEM plus \\WX\\;

  - `"sarar"`: spatial lag and spatial error;

  - `"none"`.

- effects:

  Unit heterogeneity: `"random"` (\\\alpha_i \sim N(0,
  \sigma\_\alpha^2)\\), `"pooled"` (\\\alpha = 0\\) or
  `"spatial_random"` (spatially autocorrelated effects \\\alpha = \phi
  W\alpha + a\\, \\a_i \sim N(0, \sigma\_\alpha^2)\\, as in Wang and
  Kockelman 2009; \\W\\ must then be time-invariant and \\\phi\\ uses
  `rho_range`).

- cre:

  Correlated random effects: a one-sided formula (or `TRUE` for all
  regressors) giving the regressors whose unit means are added to the
  model (Mundlak device, columns `mean_<x>`), or
  `list(type = "chamberlain", vars = ~ x1)` to add the values of the
  regressors in every period instead (Chamberlain device, columns
  `cre_<x>_<period>`). Chamberlain needs \\T\\ extra columns per
  regressor.

- durbin:

  For `"sdm"`/`"slx"`: a one-sided formula selecting the regressors
  whose spatial lags enter the model (default: all).

- time_effects:

  Logical; add period dummies.

- M:

  Spatial weights for the error term of `"sem"`, `"sdem"` and `"sarar"`
  (default: `W`). Its row names must match those of `W`. Like `W`, it
  may be a list with one matrix per period.

- dynamic:

  Optional dynamic specification from
  [`psp_dynamic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_dynamic.md)
  (lagged outcomes and initial conditions). The first period is then
  used only as the initial condition, and
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns values for
  the other periods.

- errors:

  Distribution of the latent errors: `"normal"` (probit) or Student-t
  from
  [`psp_t()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_t.md).

- priors:

  Priors, see
  [`psp_priors()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_priors.md).

- control:

  MCMC settings, see
  [`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md).

- impacts:

  Logical; compute impacts (average direct, indirect and total marginal
  effects on probabilities) after sampling.

## Value

An object of class `psprobit` (and `psoprobit` for ordered models).

## Convergence

The sampler augments the data with the latent outcomes (Albert and Chib
1993), which mixes well for balanced outcomes but slowly in some
situations:

- rare outcomes (roughly under 10% in one category);

- near-separation (outcomes that barely vary within units or periods);

- strong latent dynamics (`psp_dynamic(latent = TRUE)` with
  \\\|\gamma\|\\ near or above 1);

- `"sarar"` with `M = W`.

The intercept, \\\rho\\ and \\\sigma\_\alpha\\ are then strongly
correlated a posteriori and their chains are autocorrelated. This
slowness is inherent to data augmentation with imbalanced outcomes
(Johndrow et al. 2019), and reparameterisations do not remove it. Run
longer chains (`psp_control(iter = , warmup = )`, with `thin` to save
memory) until the warning about R-hat or bulk ESS disappears.

## References

Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *JASA* 88, 669-679.

Johndrow, J. E., Smith, A., Pillai, N. and Dunson, D. B. (2019). MCMC
for imbalanced categorical data. *JASA* 114, 1394-1403.

Kakamu, K. and Wago, H. (2005). Bayesian spatial panel probit model with
an application to business cycle in Japan. *MODSIM 2005*, 856-862.

LeSage, J. P. and Pace, R. K. (2009). *Introduction to Spatial
Econometrics*. CRC Press.

Wang, X. and Kockelman, K. M. (2009). Bayesian inference for ordered
response data with a dynamic spatial-ordered probit model. *Journal of
Regional Science* 49, 877-913.

## Examples

``` r
W <- psp_weights(rook_weights(5, 5), style = "W")
d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.5, sigma_alpha = 0.5,
                  seed = 1)
fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                control = psp_control(iter = 300, warmup = 200, chains = 2, seed = 1))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 200) for: (Intercept), x1, x2, rho, sigma_alpha. Consider more iterations.
summary(fit)
#> Bayesian panel spatial probit (model = sar, effects = random) 
#> Call: psprobit(formula = y ~ x1 + x2, data = d, index = c("id", "time"), 
#>     W = W, control = psp_control(iter = 300, warmup = 200, chains = 2, 
#>         seed = 1))
#> Units: 25, periods: 4, observations: 100
#> Chains: 2, posterior draws: 600
#> Outcome distribution:
#> 
#>  0  1 
#> 42 58 
#> 
#> Posterior summary:
#>               mean    sd   q2.5  q97.5  rhat ess_bulk ess_tail
#> (Intercept)  0.247 0.217 -0.147  0.686 1.027       94      147
#> x1           1.632 0.410  0.954  2.522 1.210        8      153
#> x2          -1.119 0.276 -1.716 -0.628 1.098       20       47
#> rho          0.080 0.156 -0.225  0.387 1.013      108      122
#> sigma_alpha  0.667 0.406  0.072  1.471 1.477        4       14
#> 
#> In-sample classification: hit rate 0.800, AUC 0.895
#> 
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 500 draws)
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1  0.280 0.031  0.219  0.335
#> x2 -0.194 0.031 -0.246 -0.129
#> 
#> indirect:
#>      mean    sd   q2.5 q97.5
#> x1  0.031 0.052 -0.053 0.150
#> x2 -0.021 0.038 -0.115 0.038
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1  0.311 0.063  0.204  0.441
#> x2 -0.215 0.050 -0.315 -0.124
```
