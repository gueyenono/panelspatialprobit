# Prior distributions

Constructors for the prior distributions used by
[`psp_priors()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_priors.md).

## Usage

``` r
prior_normal(mean = 0, sd = 10)

prior_beta(a = 1.01, b = 1.01)

prior_half_cauchy(scale = 1)

prior_inv_gamma(shape = 1, scale = 1)

prior_flat()
```

## Arguments

- mean, sd:

  Mean and standard deviation of a normal prior.

- a, b:

  Shape parameters of a beta prior, rescaled to the admissible interval
  of the spatial parameter.

- scale:

  Scale of a half-Cauchy prior on a standard deviation.

- shape:

  Shape of an inverse-gamma prior on a variance; `scale` is then its
  scale.

## Value

An object of class `psp_prior`.

## Examples

``` r
prior_normal(0, 5)
#> normal(mean = 0, sd = 5)
prior_half_cauchy(1)
#> half_cauchy(scale = 1)
```
