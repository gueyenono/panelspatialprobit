# Interpreting impacts

The coefficients of a spatial probit model are not marginal effects.
Three things stand between a coefficient and the change in a
probability:

1.  **the spatial multiplier:** a change in unit $`j`$’s regressor moves
    $`j`$’s latent outcome, which moves its neighbours’ through
    $`\rho W`$, which moves their neighbours’, and so on;
2.  **the probit link:** a change in the latent outcome moves the
    probability by an amount that depends on where the unit sits on the
    normal curve;
3.  **the scale of the latent outcome:** the spatial process and the
    random effects add variance to the latent outcome, so a unit change
    in it moves the probability less than in a plain probit.

This vignette explains what
[`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md)
computes, how the answer depends on the model, and how to work with the
draws.
[`vignette("getting-started")`](https://gueyenono.github.io/panelspatialprobit/articles/getting-started.md)
has the short version.

``` r

library(panelspatialprobit)
W <- psp_weights(rook_weights(7, 7), style = "W")
ctrl <- psp_control(iter = 3000, warmup = 1000, chains = 2, seed = 1, progress = FALSE)
```

## What is computed

Write $`S = I - \rho W`$. Integrating out the random effects, the latent
outcome of period $`t`$ has mean $`m_t = S^{-1}(X_t\beta + WX_t\theta)`$
and the standard deviation of unit $`i`$’s latent outcome is
$`s_i = \sqrt{(1 + \sigma_\alpha^2)\,[S^{-1}S^{-\top}]_{ii}}`$ (with
spatial errors, $`S^{-1}R^{-1}`$ takes the place of $`S^{-1}`$ in the
error part). So $`P(y_{it} = 1) = \Phi(m_{it}/s_i)`$ and

``` math
\frac{\partial P(y_{it} = 1)}{\partial x_{jt,r}} =
\frac{\phi(m_{it}/s_i)}{s_i}\,\big[S^{-1}(\beta_r I + \theta_r W)\big]_{ij}.
```

This is an $`N \times N`$ matrix of effects for every regressor, period
and posterior draw. Following LeSage and Pace (2009), it is summarised
by three averages:

- **direct:** the average diagonal element, the effect of a unit’s own
  regressor on its own probability (including the feedback that returns
  to it through its neighbours);
- **total:** the average row sum, the effect on a unit’s probability of
  the same change in the regressor of every unit;
- **indirect:** total minus direct, the spillover.

The averages run over units and periods; `by_period = TRUE` also reports
them period by period. Every quantity is computed for each posterior
draw, so the impacts come with posterior intervals.

## A spatial Durbin example

We simulate data in which $`x_1`$ also acts through the neighbours’
values ($`\theta_1 = 0.5`$) and $`x_2`$ does not ($`\theta_2 = 0`$), and
fit a spatial Durbin model:

``` r

d <- sim_psprobit(W, T = 5, beta = c(-0.2, 1, -0.8), theta = c(0.5, 0), rho = 0.4,
                  sigma_alpha = 0.7, seed = 1)
fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                model = "sdm", control = ctrl)
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 200) for:
#> rho. Consider more iterations.
fit$impacts
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 500 draws)
#> 
#> direct:
#>       mean     sd    q2.5   q97.5
#> x1  0.1883 0.0201  0.1480  0.2265
#> x2 -0.1413 0.0202 -0.1814 -0.1029
#> 
#> indirect:
#>       mean     sd    q2.5  q97.5
#> x1  0.2874 0.0581  0.1738 0.4018
#> x2 -0.0725 0.0538 -0.1741 0.0329
#> 
#> total:
#>       mean     sd    q2.5   q97.5
#> x1  0.4757 0.0588  0.3562  0.5913
#> x2 -0.2138 0.0611 -0.3295 -0.0939
```

The spillover of $`x_1`$ is larger than its direct effect: its own
Durbin term adds to the multiplier. In a SAR model (no $`\theta`$) the
matrix in brackets is $`\beta_r S^{-1}`$ for every regressor, so the
ratio of indirect to direct impacts is the same for all regressors. That
is a strong restriction; the Durbin terms relax it.

## Which models have spillovers

The model fixes which channels exist:

| model | spillovers through | indirect impacts |
|----|----|----|
| `"sar"` | $`\rho`$: global (all units, decaying with distance) | yes, same ratio for all regressors |
| `"sdm"` | $`\rho`$ and $`\theta`$ | yes, free ratio |
| `"slx"` | $`\theta`$ only: neighbours of the changed unit | local only |
| `"sdem"` | $`\theta`$ only; $`\lambda`$ is in the errors | local only |
| `"sem"` | none; $`\lambda`$ is in the errors | zero |
| `"sarar"` | $`\rho`$ | as `"sar"` |

A spatial error term ($`\lambda`$) does not carry a change in a
regressor from one unit to another; it only changes $`s_i`$. The SEM’s
indirect impacts are therefore exactly zero:

``` r

fit_sem <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                    model = "sem", control = ctrl)
fit_sem$impacts
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 500 draws)
#> 
#> direct:
#>       mean     sd    q2.5   q97.5
#> x1  0.1888 0.0272  0.1305  0.2365
#> x2 -0.1497 0.0227 -0.1955 -0.1055
#> 
#> indirect:
#>    mean sd q2.5 q97.5
#> x1    0  0    0     0
#> x2    0  0    0     0
#> 
#> total:
#>       mean     sd    q2.5   q97.5
#> x1  0.1888 0.0272  0.1305  0.2365
#> x2 -0.1497 0.0227 -0.1955 -0.1055
```

Here the SEM is the wrong model: it attributes the spillover of $`x_1`$
to nothing, and its direct impacts are close to the SDM’s. With real
data, compare the specifications with
[`psp_ppc()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_ppc.md)
and
[`psp_dic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_loglik.md)
(see
[`vignette("getting-started")`](https://gueyenono.github.io/panelspatialprobit/articles/getting-started.md))
rather than choosing by the size of the impacts.

## The scale

The division by $`s_i`$ matters. In a plain probit the latent error has
variance one, so $`\beta`$ is measured in standard deviations of the
latent outcome. Here the variance of the latent outcome is larger, by
the factor $`(1 + \sigma_\alpha^2)[S^{-1}S^{-\top}]_{ii}`$, and a
probability changes less for the same latent change. The `spatialprobit`
package (Wilhelm and de Matos 2013) leaves this factor out
($`P = \Phi(m_{it})`$). `scale = "lesage_pace"` reproduces its
convention, for comparison with published results:

``` r

psp_impacts(fit, scale = "lesage_pace")
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 500 draws)
#> Scale: LeSage-Pace (no reduced-form sd)
#> 
#> direct:
#>       mean     sd    q2.5   q97.5
#> x1  0.2160 0.0223  0.1709  0.2584
#> x2 -0.1621 0.0219 -0.2074 -0.1221
#> 
#> indirect:
#>       mean     sd    q2.5  q97.5
#> x1  0.3305 0.0696  0.2022 0.4781
#> x2 -0.0828 0.0616 -0.1979 0.0388
#> 
#> total:
#>       mean     sd    q2.5   q97.5
#> x1  0.5464 0.0713  0.4161  0.6909
#> x2 -0.2449 0.0685 -0.3750 -0.1120
```

The impacts are larger, because they ignore the extra variance. The
default (`"reduced_form"`) gives the derivatives of the probabilities
the model actually implies; we recommend it for interpretation.

## Population and unit impacts

The population impacts integrate the random effects out: they describe a
population of units like those in the sample. `type = "unit"` instead
conditions on the posterior draws of the random effects, so that
$`m_t = S^{-1}(X_t\beta + WX_t\theta + \alpha)`$ and
$`s_i^2 = [S^{-1}S^{-\top}]_{ii}`$. These are the impacts for the units
in the sample:

``` r

psp_impacts(fit, type = "unit", draws = 300)
#> Impacts on P(y = 1) (posterior mean, sd, 95% interval; 300 draws)
#> Unit-specific (conditional on the random effects)
#> 
#> direct:
#>       mean     sd    q2.5   q97.5
#> x1  0.1968 0.0206  0.1574  0.2376
#> x2 -0.1463 0.0222 -0.1881 -0.1003
#> 
#> indirect:
#>       mean     sd    q2.5  q97.5
#> x1  0.3108 0.0650  0.1865 0.4418
#> x2 -0.0750 0.0545 -0.1702 0.0451
#> 
#> total:
#>       mean     sd    q2.5   q97.5
#> x1  0.5076 0.0663  0.3632  0.6382
#> x2 -0.2214 0.0646 -0.3334 -0.0827
```

Unit impacts are larger than population impacts. A unit with a known
random effect has less latent uncertainty left, so its probability
responds more sharply. Averaging over the random effects flattens the
response. Both are correct answers to different questions.

## Working with the draws

`imp$draws` is an array with dimensions draw × variable × category ×
effect, so any derived quantity gets a posterior distribution. For
example, the share of the total effect of $`x_1`$ that is spillover:

``` r

dr <- fit$impacts$draws
share <- dr[, "x1", 1, "indirect"] / dr[, "x1", 1, "total"]
quantile(share, c(0.05, 0.5, 0.95))
#>        5%       50%       95% 
#> 0.4976683 0.6036779 0.6879632
```

`summary(imp)` returns the table as a data frame, one row per variable,
category and effect (and period with `by_period = TRUE`).

## Ordered outcomes

For an ordered model the impacts are computed for every category:
$`\partial P(y_{it} = k)/\partial x_{jt,r}`$. Raising a regressor with a
positive coefficient lowers the probabilities of the low categories and
raises those of the high ones, and the impacts sum to zero over
categories. See
[`vignette("ordered-outcomes")`](https://gueyenono.github.io/panelspatialprobit/articles/ordered-outcomes.md).

## Derivatives and discrete changes

Impacts are derivatives: they describe small changes. For a discrete
change, such as a policy that raises $`x_1`$ by one unit everywhere,
compare predictions directly. `predict(newdata = )` takes the same units
and periods with new regressor values:

``` r

d_up <- d
d_up$x1 <- d_up$x1 + 1
mean(predict(fit, newdata = d_up) - predict(fit))
#> [1] 0.4231255
```

The average change in probability is close to the total impact of
$`x_1`$ above, but smaller: probabilities are bounded, and the normal
curve flattens as units move into the tails. For changes of more than a
fraction of a standard deviation of the regressor, prefer the discrete
comparison.

Changing the regressor for some units only works the same way: the
predicted change for the other units is then their spillover.

## Dynamic models

With a lagged outcome, the impacts above are short-run: they describe
the period in which the regressor changes. The change also works through
later periods via the outcome history.
[`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md)
traces that; see
[`vignette("dynamics")`](https://gueyenono.github.io/panelspatialprobit/articles/dynamics.md).

## References

LeSage, J. P. and Pace, R. K. (2009). *Introduction to Spatial
Econometrics*. CRC Press.

Wilhelm, S. and de Matos, M. G. (2013). Estimating spatial probit models
in R. *The R Journal* 5(1), 130–143.
