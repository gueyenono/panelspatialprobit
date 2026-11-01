# Impacts (marginal effects) on outcome probabilities

Computes average direct, indirect (spillover) and total impacts of each
regressor on the outcome probabilities, following LeSage and Pace (2009,
ch. 10) but accounting for the reduced-form standard deviation of the
latent outcome. Impacts are population-averaged: the random effects are
integrated out. For ordered models impacts are computed per category
(they sum to zero over categories).

## Usage

``` r
psp_impacts(
  fit,
  draws = NULL,
  scale = c("reduced_form", "lesage_pace"),
  type = c("population", "unit"),
  by_period = FALSE
)

# S3 method for class 'psprobit'
impacts(obj, ...)
```

## Arguments

- fit:

  A fitted `psprobit` object.

- draws:

  Number of posterior draws to use (default `control$impacts_draws` of
  the fit).

- scale:

  `"reduced_form"` (default) divides by \\s_i\\ as above;
  `"lesage_pace"` uses \\\phi(m\_{it})\\ without it, as in the
  `spatialprobit` package, for comparison with published results.

- type:

  `"population"` (default, random effects integrated out) or `"unit"`
  (conditional on the random effects; needs `effects = "random"`).

- by_period:

  Logical; also report impacts averaged within each period (the default
  averages over all observations).

- obj:

  A fitted `psprobit` object (method for
  [`spatialreg::impacts`](https://r-spatial.github.io/spatialreg/reference/impacts.html)).

- ...:

  Unused.

## Value

An object of class `psp_impacts` with draws and a summary table.

## Details

With `type = "unit"` the impacts are instead conditional on the
posterior draws of the random effects: \\m_t = S^{-1}(X_t\beta +
\alpha)\\ and \\s_i^2 = \[S^{-1}S^{-\top}\]\_{ii}\\. They describe the
sample units rather than the population and use the random-effect draws
kept by the sampler (about `control$impacts_draws` of them).

For regressor \\r\\ and period \\t\\, \$\$\partial P(y\_{it} = 1) /
\partial x\_{jt,r} = \phi(m\_{it}/s_i)/s_i \\ \[S^{-1}(\beta_r I +
\theta_r W)\]\_{ij},\$\$ where \\m_t = S^{-1}X_t\beta\\ and \\s_i^2 =
(1 + \sigma\_\alpha^2) \[S^{-1}S^{-\top}\]\_{ii}\\. Direct impacts
average the \\i = j\\ terms, total impacts average row sums, and
indirect = total - direct.

The posterior draws are processed in parallel under the user's
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(see "Parallel post-processing" in
[`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md)).

## Examples

``` r
W <- psp_weights(rook_weights(4, 4), style = "W")
d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
                control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, rho, sigma_alpha. Consider more iterations.
psp_impacts(fit, draws = 50)
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 50 draws)
#> 
#> direct:
#>      mean     sd  q2.5  q97.5
#> x1 0.1517 0.0314 0.107 0.2208
#> 
#> indirect:
#>       mean     sd    q2.5  q97.5
#> x1 -0.0022 0.0309 -0.0344 0.0815
#> 
#> total:
#>      mean     sd   q2.5  q97.5
#> x1 0.1496 0.0474 0.0914 0.2843
```
