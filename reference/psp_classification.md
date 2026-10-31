# In-sample classification statistics

Compares posterior mean outcome probabilities (`predict(type = "prob")`)
with the observed outcomes (missing outcomes are left out).

## Usage

``` r
psp_classification(fit, draws = NULL)
```

## Arguments

- fit:

  A fitted `psprobit` object.

- draws:

  Number of posterior draws used for the probabilities.

## Value

A list with `hit_rate` (share of observations whose most probable
outcome is the observed one), `auc` (binary models: area under the ROC
curve; `NA` for ordered models), `brier` (mean squared error of the
probability vector) and `confusion` (observed by predicted table).

## Examples

``` r
W <- psp_weights(rook_weights(4, 4), style = "W")
d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
                control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, rho, sigma_alpha. Consider more iterations.
psp_classification(fit, draws = 50)
#> $hit_rate
#> [1] 0.6041667
#> 
#> $auc
#> [1] 0.6381119
#> 
#> $brier
#> [1] 0.2349346
#> 
#> $confusion
#>         predicted
#> observed  0  1
#>        0 20  6
#>        1 13  9
#> 
```
