# Dynamic response paths

Traces how a change in a regressor propagates over time in a dynamic
model
([`psp_dynamic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_dynamic.md)):
the shock changes the outcome probabilities in the period it occurs (the
short-run impact), which changes the outcomes, which enter the next
periods through the lagged outcomes (`lag_y`, `W_lag_y`), and so on.

## Usage

``` r
psp_response_path(
  fit,
  variable,
  shock = 1,
  start = NULL,
  units = NULL,
  type = c("temporary", "permanent"),
  draws = 200,
  nsim = 20,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted `psprobit` object.

- variable:

  Name of the regressor to shock.

- shock:

  Size of the change in `variable`.

- start:

  Period (a value of the time index) in which the shock starts; default
  the first modelled period.

- units:

  Identifiers of the shocked units (default: all units).

- type:

  `"temporary"` (the shock lasts one period) or `"permanent"` (it lasts
  from `start` to the last period).

- draws:

  Number of posterior draws.

- nsim:

  Number of simulated outcome histories per draw.

- seed:

  Optional seed for the simulation.

## Value

An object of class `psp_response`: a list with `summary`, a data frame
with the posterior mean, sd and 95% interval of the response by
`horizon` (0 = start period), `period`, `category` (ordered models) and
`units` (`"all"` averages the change over every unit; when only some
units are shocked, `"shocked"` averages over them and `"others"` over
the remaining ones), and `draws`, the responses for each posterior draw.
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods are
available.

## Details

For each posterior draw, outcome histories are simulated forward from
the start period under the baseline (the observed regressors) and under
the shock, with the same random numbers in both scenarios. Outcomes
before the start period are the observed ones. In every period the
outcome probabilities given the simulated history are computed exactly,
and the response is their average difference between the two scenarios.
The simulation conditions on the posterior draws of the random effects
kept by the sampler (as
[`psp_ppc()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_ppc.md)
does), so it describes the sample units. Correlated-random-effects
columns (`mean_x`, `cre_x_*`) and the initial conditions are held fixed.

With latent dynamics (`psp_dynamic(latent = TRUE)`) the shock propagates
through the latent mean, \\\Delta m\_{t+h} = (\gamma S^{-1})^h \Delta
m_t\\, and the responses are computed exactly from the reduced form (no
simulation; `nsim` is unused). Static models are accepted too; their
responses vanish after the shock ends.

The posterior draws are processed in parallel under the user's
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(see "Parallel post-processing" in
[`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md)).

## Examples

``` r
W <- psp_weights(rook_weights(5, 5), style = "W")
d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1, -1), rho = 0.3, gamma = 0.8, seed = 1)
fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, dynamic = psp_dynamic(),
                impacts = FALSE,
                control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, x2, lag_y, init_y, mean_x1, mean_x2, rho, sigma_alpha. Consider more iterations.
psp_response_path(fit, "x1", draws = 50, nsim = 5, seed = 1)
#> Response of P(y = 1) to a temporary shock of 1 in x1 (50 draws x 5 histories)
#> 
#> Average change over all units:
#>  horizon period   mean     sd   q2.5  q97.5
#>        0      2 0.2352 0.0671 0.1195 0.3668
#>        1      3 0.1068 0.0622 0.0217 0.2519
#>        2      4 0.0463 0.0382 0.0029 0.1422
#>        3      5 0.0248 0.0238 0.0000 0.0734
#>        4      6 0.0123 0.0151 0.0000 0.0476
```
