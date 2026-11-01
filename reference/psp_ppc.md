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

The posterior draws are processed in parallel under the user's
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(see "Parallel post-processing" in
[`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md)).

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
#> P(y=1) t=1     0.48    0.538    0.236     0.720 0.80
#> P(y=1) t=2     0.64    0.653    0.449     0.800 0.70
#> P(y=1) t=3     0.52    0.493    0.320     0.640 0.54
#> P(y=1) t=4     0.60    0.523    0.360     0.671 0.28
#> 
#> moran:
#>           observed rep_mean rep_q2.5 rep_q97.5  ppp
#> Moran t=1    0.145    0.098   -0.138     0.307 0.36
#> Moran t=2    0.282    0.122   -0.189     0.384 0.14
#> Moran t=3   -0.100    0.016   -0.233     0.343 0.78
#> Moran t=4    0.283    0.095   -0.208     0.314 0.04
#> 
#> persistence:
#>             observed rep_mean rep_q2.5 rep_q97.5  ppp
#> persist t=2     0.68    0.579    0.440     0.813 0.24
#> persist t=3     0.64    0.600    0.378     0.751 0.50
#> persist t=4     0.60    0.614    0.409     0.760 0.60
```
