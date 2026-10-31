# Methods for fitted panel spatial probit models

Methods for fitted panel spatial probit models

## Usage

``` r
# S3 method for class 'psprobit'
print(x, digits = 3, ...)

# S3 method for class 'psprobit'
summary(object, classification = TRUE, ...)

# S3 method for class 'psprobit'
coef(object, ...)

# S3 method for class 'psprobit'
vcov(object, ...)

# S3 method for class 'psprobit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'psprobit'
nobs(object, ...)

# S3 method for class 'psprobit'
predict(
  object,
  newdata = NULL,
  type = c("prob", "link", "class"),
  draws = NULL,
  ...
)

# S3 method for class 'psprobit'
fitted(object, ...)

# S3 method for class 'psprobit'
plot(x, pars = NULL, ...)

# S3 method for class 'psprobit'
as_draws(x, ...)

# S3 method for class 'psprobit'
tidy(x, conf.level = 0.95, ...)

# S3 method for class 'psprobit'
glance(x, ...)
```

## Arguments

- x, object:

  A fitted `psprobit` object.

- digits:

  Number of digits to print.

- ...:

  Unused.

- classification:

  Logical; report in-sample classification statistics (see
  [`psp_classification()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_classification.md))
  in [`summary()`](https://rdrr.io/r/base/summary.html).

- parm:

  Parameters for which intervals are returned (default: all).

- level:

  Credible level.

- newdata:

  Optional data frame for counterfactual predictions: the same units and
  periods as the fitted data, with new regressor values (the response
  column may be left out). Correlated-random-effects columns are kept at
  their fitted values, since they stand for fixed unit characteristics.
  Not available for models with lagged outcomes.

- type:

  For [`predict()`](https://rdrr.io/r/stats/predict.html): `"prob"`
  (posterior mean probabilities; a matrix with one column per category
  for ordered models), `"link"` (posterior mean of the systematic part
  of the latent outcome) or `"class"` (most probable outcome).

- draws:

  Number of posterior draws used by
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- pars:

  Parameters to plot (default: all).

- conf.level:

  Credible level for
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html).

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix,
[`confint()`](https://rdrr.io/r/stats/confint.html) a matrix of interval
bounds, [`summary()`](https://rdrr.io/r/base/summary.html) an object of
class `summary.psprobit`,
[`predict()`](https://rdrr.io/r/stats/predict.html)/[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
vectors or matrices in the row order of the data.
