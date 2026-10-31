# Posterior predictive checks

Simulates replicated outcomes from the posterior and compares summary
statistics of the replicates with those of the observed outcomes
(Gelman, Meng and Stern 1996). Replicates are drawn for the same units,
periods and regressors, conditional on the posterior draws of the random
effects kept by the sampler (pooled models have none). Missing outcomes
are left out of every statistic.

## Usage

``` r
psp_ppc(fit, draws = NULL, seed = NULL)
```

## Arguments

- fit:

  A fitted `psprobit` object.

- draws:

  Number of posterior draws to replicate (default
  `control$impacts_draws`).

- seed:

  Optional seed for the replicates.

## Value

An object of class `psp_ppc`: a list with, for each statistic, the
observed values (`observed`), the replicated values (`replicated`, draws
x components) and posterior predictive p-values \\P(T(y^{rep}) \ge
T(y))\\ (`ppp`). Its [`print()`](https://rdrr.io/r/base/print.html)
method shows a summary.

## Details

The statistics are:

- `share`: the share of each outcome category (of \\y = 1\\ for binary
  models) in each period;

- `moran`: Moran's I of the outcome (category codes for ordered models)
  in each period, which checks the spatial dependence;

- `persistence`: the share of units whose outcome equals that of the
  previous period, which checks the dependence over time (\\T \ge 2\\).

## References

Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica* 6, 733-807.

## Examples

``` r
W <- psp_weights(rook_weights(5, 5), style = "W")
d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.5, seed = 1)
fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, impacts = FALSE,
                control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, x2, rho, sigma_alpha. Consider more iterations.
psp_ppc(fit, draws = 50, seed = 1)
#> Posterior predictive checks (50 replicates)
#> ppp = P(T(y_rep) >= T(y)); values near 0 or 1 signal misfit.
#> 
#> share:
#>            observed rep_mean rep_q2.5 rep_q97.5  ppp
#> P(y=1) t=1     0.48    0.552    0.320     0.680 0.84
#> P(y=1) t=2     0.64    0.630    0.449     0.791 0.60
#> P(y=1) t=3     0.52    0.502    0.360     0.680 0.50
#> P(y=1) t=4     0.60    0.529    0.258     0.680 0.34
#> 
#> moran:
#>           observed rep_mean rep_q2.5 rep_q97.5  ppp
#> Moran t=1    0.145    0.084   -0.213     0.332 0.32
#> Moran t=2    0.282    0.149   -0.090     0.394 0.18
#> Moran t=3   -0.100   -0.009   -0.290     0.276 0.70
#> Moran t=4    0.283    0.103   -0.146     0.386 0.08
#> 
#> persistence:
#>             observed rep_mean rep_q2.5 rep_q97.5  ppp
#> persist t=2     0.68    0.568    0.400     0.782 0.26
#> persist t=3     0.64    0.621    0.480     0.760 0.50
#> persist t=4     0.60    0.609    0.449     0.760 0.66
```
