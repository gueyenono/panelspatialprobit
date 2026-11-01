# Dynamic specification

Adds observed state dependence to
[`psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
/
[`psoprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md):
the previous period's outcome \\y\_{i,t-1}\\ and, optionally, its
spatial lag \\(W y\_{t-1})\_i\\ enter the latent equation as regressors.
The first period is used only as the initial condition, so estimation
uses periods \\2, \dots, T\\. For ordered outcomes the lag enters as one
indicator per category above the first.

## Usage

``` r
psp_dynamic(
  lag_y = !latent,
  lag_Wy = FALSE,
  initial = c("wooldridge", "none"),
  latent = FALSE
)
```

## Arguments

- lag_y:

  Logical; include \\y\_{i,t-1}\\ (columns `lag_*`).

- lag_Wy:

  Logical; include \\(W y\_{t-1})\_i\\ (columns `W_lag_*`).

- initial:

  `"wooldridge"` or `"none"` (treat the first-period outcome as
  exogenous). Ignored with `latent = TRUE`.

- latent:

  Logical; latent AR(1) dynamics instead of lagged outcomes.

## Value

An object of class `psp_dynamic`.

## Details

With `initial = "wooldridge"` (the default) the random effects are
modelled as \\\alpha_i = \psi' y\_{i1} + \bar x_i'\pi + a_i\\
(Wooldridge 2005): the first-period outcome (columns `init_*`) and the
unit means of all regressors over the estimation periods (columns
`mean_*`) are added to the model. This handles the initial-conditions
problem of dynamic random-effects models.

Impacts of the regressors are short-run (within-period) impacts.

With `latent = TRUE` the dynamics act on the latent outcome instead,
\\z_t = \rho W z_t + \gamma z\_{t-1} + X_t\beta + \alpha + u_t\\ (Wang
and Kockelman 2009). The coefficient \\\gamma\\ is reported as `lag_z`.
Every period is modelled; the first has no lagged term, so its latent
outcome serves as the initial condition and shares the random effect
with the later periods (no separate initial-conditions device is
needed). Probabilities, impacts and predictions integrate over the
latent history, and
[`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md)
gives the propagation of a shock. Latent dynamics need time-invariant
weights and cannot be combined with lagged outcomes. \\\gamma\\ mixes
slowly against \\\rho\\ (it is drawn given the latent values), so plan
on longer chains than for static models and check the convergence
diagnostics.

## References

Wang, X. and Kockelman, K. M. (2009). Bayesian inference for ordered
response data with a dynamic spatial-ordered probit model. *Journal of
Regional Science* 49, 877-913.

Wooldridge, J. M. (2005). Simple solutions to the initial conditions
problem in dynamic, nonlinear panel data models with unobserved
heterogeneity. *Journal of Applied Econometrics* 20, 39-54.

## Examples

``` r
W <- psp_weights(rook_weights(5, 5), style = "W")
d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1, -1), rho = 0.3, gamma = 0.8, seed = 1)
fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, dynamic = psp_dynamic(),
                impacts = FALSE,
                control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, x2, lag_y, init_y, mean_x1, mean_x2, rho, sigma_alpha. Consider more iterations.
coef(fit)
#> (Intercept)          x1          x2       lag_y      init_y     mean_x1 
#>  -0.8280242   1.3557468  -1.6341862   1.9007703   0.3654488   1.4835307 
#>     mean_x2         rho sigma_alpha 
#>   1.0649738   0.1102372   1.9517472 
```
