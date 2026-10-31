# Simulate panel spatial probit and ordered probit data

Draws data from the latent model \$\$z_t = \rho W z_t + X_t\beta + W
X_t\theta + \alpha + u_t, \quad u_t = \lambda M u_t + \varepsilon_t\$\$
with standard normal regressors. The outcome is binary
(`cutpoints = NULL`) or ordered with cut-points \\(0, \kappa_2, \dots,
\kappa\_{J-1})\\.

## Usage

``` r
sim_psprobit(
  W,
  T,
  beta = c(0, 1, -1),
  rho = 0.5,
  sigma_alpha = 1,
  theta = NULL,
  cutpoints = NULL,
  seed = NULL,
  lambda = 0,
  M = NULL,
  gamma = NULL,
  eta = NULL,
  phi = 0,
  df = NULL,
  gamma_z = 0
)
```

## Arguments

- W:

  Spatial weights (`N x N`), used as given, or a list of `T` such
  matrices (time-varying weights, period by period).

- T:

  Number of periods.

- beta:

  Coefficients; the first is the intercept, the others belong to
  regressors `x1`, `x2`, ...

- rho:

  Spatial autoregressive parameter (0 for `"slx"`/`"none"`).

- sigma_alpha:

  Standard deviation of the random effects (0 for pooled).

- theta:

  Optional Durbin coefficients (one per regressor), giving the effects
  of \\WX\\.

- cutpoints:

  Optional free cut-points \\\kappa_2 \< \dots \< \kappa\_{J-1}\\, all
  positive (the first cut-point is fixed at 0). For example
  `cutpoints = 1` gives \\J = 3\\ categories.

- seed:

  Optional seed.

- lambda:

  Spatial error parameter (0 for no spatial error term).

- M:

  Error weights (default `W`); a matrix or a list of `T` matrices.

- gamma:

  Optional coefficient(s) of the lagged outcome \\y\_{i,t-1}\\ (one per
  category above the first for ordered outcomes), as in
  [`psp_dynamic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_dynamic.md).
  The first period is generated without lags.

- eta:

  Optional coefficient(s) of the spatially lagged outcome \\(W
  y\_{t-1})\_i\\.

- phi:

  Spatial parameter of the random effects, \\\alpha = \phi W \alpha +
  a\\ with \\a_i \sim N(0, \sigma\_\alpha^2)\\ (as for
  `effects = "spatial_random"`; needs a single `W`).

- df:

  Optional degrees of freedom of Student-t latent errors (as for
  `errors = psp_t(df)`); normal errors by default.

- gamma_z:

  Coefficient of latent AR(1) dynamics, \\\gamma z\_{t-1}\\ (as for
  `psp_dynamic(latent = TRUE)`; no lag in the first period).

## Value

A data frame with columns `id`, `time`, `y` and the regressors, with
attributes `truth` (parameters, random effects and latent values) and
`W` (the weights as passed, after
[`psp_weights()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_weights.md)).

## Examples

``` r
W <- psp_weights(rook_weights(3, 3), style = "W")
d <- sim_psprobit(W, T = 2, beta = c(0, 1), rho = 0.3, seed = 1)
head(d)
#>   id time y         x1
#> 1  1    1 0 -0.6264538
#> 2  2    1 1  0.1836433
#> 3  3    1 1 -0.8356286
#> 4  4    1 1  1.5952808
#> 5  5    1 1  0.3295078
#> 6  6    1 0 -0.8204684
```
