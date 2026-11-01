# Student-t latent errors

Replaces the normal errors \\\varepsilon\_{it} \sim N(0, 1)\\ of
[`psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
/
[`psoprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
by Student-t errors with `df` degrees of freedom, which makes the model
robust to outlying observations (for binary outcomes this is the robit
model of Liu 2004; `df = 7` approximates the logit link). The errors are
sampled as the scale mixture \\\varepsilon\_{it} \sim N(0, v\_{it})\\,
\\v\_{it} \sim IG(\nu/2, \nu/2)\\; the posterior means of \\v\_{it}\\
are returned in `fit$v_mean` (large values flag observations the model
fits poorly).

## Usage

``` r
psp_t(df = 7)
```

## Arguments

- df:

  Degrees of freedom, larger than 2.

## Value

An object of class `psp_errors`.

## Details

The degrees of freedom are fixed rather than estimated: with discrete
outcomes the data carry very little information about them.

With t errors the latent outcome is no longer normal, because the
spatial filter mixes the errors of several units. Probabilities and
impacts are therefore averaged over `control$mixture_draws` simulated
sets of \\v\\ (fixed per fit, so results are reproducible): given \\v\\
the latent outcome is normal and the usual formulas apply. Above
`impacts_exact_max` units the reduced-form variance treats the own error
term exactly and the neighbours' terms at their mean.

## References

Liu, C. (2004). Robit regression: a simple robust alternative to
logistic and probit regression. In Gelman, A. and Meng, X.-L. (eds.),
*Applied Bayesian Modeling and Causal Inference from Incomplete-Data
Perspectives*, 227-238. Wiley.

## Examples

``` r
W <- psp_weights(rook_weights(5, 5), style = "W")
d <- sim_psprobit(W, T = 4, beta = c(0, 1, -1), rho = 0.4, df = 5, seed = 1)
fit <- psprobit(y ~ x1 + x2, d, c("id", "time"), W, errors = psp_t(df = 5),
                impacts = FALSE,
                control = psp_control(iter = 300, warmup = 200, chains = 1, seed = 1))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, x2, rho, sigma_alpha. Consider more iterations.
coef(fit)
#> (Intercept)          x1          x2         rho sigma_alpha 
#>   0.1290682   1.3537639  -1.1033764   0.3442226   1.5044671 
```
