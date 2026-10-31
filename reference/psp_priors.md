# Priors for panel spatial probit models

Priors for panel spatial probit models

## Usage

``` r
psp_priors(
  beta = prior_normal(0, 10),
  rho = prior_beta(1.01, 1.01),
  lambda = prior_beta(1.01, 1.01),
  sigma_alpha = prior_half_cauchy(1),
  cutpoints = prior_normal(0, 10),
  phi = prior_beta(1.01, 1.01)
)
```

## Arguments

- beta:

  Normal prior for all regression coefficients (including Durbin,
  correlated-random-effects and time-effect columns). `mean` and `sd`
  may be scalars or vectors with one entry per column of the design
  matrix.

- rho:

  Beta prior for the spatial autoregressive parameter, rescaled to the
  interval given by `rho_range` in
  [`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md).

- lambda:

  Beta prior for the spatial error parameter (models `"sem"`, `"sdem"`
  and `"sarar"`), rescaled to `lambda_range` in
  [`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md).

- sigma_alpha:

  Prior for the standard deviation of the random effects:
  [`prior_half_cauchy()`](https://gueyenono.github.io/panelspatialprobit/reference/priors.md)
  (default) or
  [`prior_inv_gamma()`](https://gueyenono.github.io/panelspatialprobit/reference/priors.md)
  on the variance.

- cutpoints:

  Prior for the free cut-points of ordered models:
  [`prior_normal()`](https://gueyenono.github.io/panelspatialprobit/reference/priors.md)
  (truncated to the ordering constraint) or
  [`prior_flat()`](https://gueyenono.github.io/panelspatialprobit/reference/priors.md).

- phi:

  Beta prior for the spatial parameter of spatial random effects
  (`effects = "spatial_random"`), rescaled to `rho_range`.

## Value

An object of class `psp_priors`.

## Examples

``` r
psp_priors(beta = prior_normal(0, 5))
#> $beta
#> normal(mean = 0, sd = 5)
#> 
#> $rho
#> beta(a = 1.01, b = 1.01)
#> 
#> $lambda
#> beta(a = 1.01, b = 1.01)
#> 
#> $sigma_alpha
#> half_cauchy(scale = 1)
#> 
#> $cutpoints
#> normal(mean = 0, sd = 10)
#> 
#> $phi
#> beta(a = 1.01, b = 1.01)
#> 
#> attr(,"class")
#> [1] "psp_priors"
```
