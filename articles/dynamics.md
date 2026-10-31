# Dynamic panel spatial probit

Outcomes observed over time are usually persistent: a region in
recession this year is likely to be in recession next year, a firm that
exported last year is likely to export again. Three different mechanisms
produce this pattern, and they have different implications:

- **state dependence:** the outcome itself changes the next period’s
  probability (exporting builds contacts that make the next export
  cheaper);
- **heterogeneity:** some units are persistently more likely to have the
  outcome, for reasons the regressors do not capture (the random effects
  $`\alpha_i`$);
- **persistent shocks:** the unobserved latent factors are themselves
  autocorrelated, so the propensity persists even when the outcome does
  not feed back.

In a spatial panel these mix with a fourth: neighbours’ outcomes. This
vignette shows how `panelspatialprobit` separates them and how to read
the dynamic effects.

``` r

library(panelspatialprobit)
W <- psp_weights(rook_weights(7, 7), style = "W")
ctrl <- psp_control(iter = 3000, warmup = 1000, chains = 2, seed = 1, progress = FALSE)
```

## State dependence: lagged outcomes

`dynamic = psp_dynamic()` adds the previous period’s outcome to the
model:

``` math
z_{it} = \rho (Wz_t)_i + \gamma y_{i,t-1} + \eta (Wy_{t-1})_i + x_{it}'\beta + \alpha_i +
\varepsilon_{it},
```

with $`\gamma`$ reported as `lag_y` and, with `lag_Wy = TRUE`, $`\eta`$
reported as `W_lag_y` (the neighbours’ previous outcomes, a lagged
spatial diffusion channel). For ordered outcomes the lag enters as one
indicator per category above the first.

The first period has no previous outcome. It serves as the **initial
condition** and is not modelled itself, so a panel of $`T`$ periods
gives $`T - 1`$ modelled periods.

### The initial conditions problem

The random effect $`\alpha_i`$ influenced the first outcome $`y_{i1}`$
too, so the lagged outcome is correlated with $`\alpha_i`$. Treating
$`y_{i1}`$ as fixed and independent of $`\alpha_i`$ attributes part of
the heterogeneity to state dependence and overstates $`\gamma`$. The
default, `psp_dynamic(initial = "wooldridge")`, follows Wooldridge
(2005): it models $`\alpha_i`$ given the initial outcome and the
regressors,
``` math
\alpha_i = \xi_0 y_{i1} + \bar x_i'\xi + a_i,
```
by adding `init_y` and the unit means of the regressors (`mean_x1`, …)
to the model. `a_i` is then the random effect, independent of the
initial condition.

``` r

d <- sim_psprobit(W, T = 8, beta = c(-0.4, 1, -0.8), rho = 0.4, sigma_alpha = 0.6,
                  gamma = 0.8, seed = 3)
fit <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                dynamic = psp_dynamic(), impacts = FALSE, control = ctrl)
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 200) for:
#> sigma_alpha. Consider more iterations.
round(coef(fit), 2)
#> (Intercept)          x1          x2       lag_y      init_y     mean_x1 
#>       -0.59        0.74       -0.84        0.80        0.45        0.35 
#>     mean_x2         rho sigma_alpha 
#>        0.48        0.37        0.35
```

The warning concerns $`\sigma_\alpha`$, which the Wooldridge terms leave
little to explain; in an application, run longer chains until it
disappears. The true $`\gamma`$ is 0.8. A pooled model that ignores the
heterogeneity puts the persistence due to $`\alpha_i`$ into `lag_y`:

``` r

fit_pooled <- psprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                       effects = "pooled", dynamic = psp_dynamic(initial = "none"),
                       impacts = FALSE, control = ctrl)
round(coef(fit_pooled), 2)
#> (Intercept)          x1          x2       lag_y         rho 
#>       -0.52        0.71       -0.71        0.98        0.38
```

One data set proves little, but the pattern is systematic. In a Monte
Carlo study with 200 replications
([`vignette("validation")`](https://gueyenono.github.io/panelspatialprobit/articles/validation.md)),
the Wooldridge model recovered $`\gamma = 0.8`$ with a bias of 0.01 and
93% coverage of 90% intervals, while the pooled model overstated it by
0.31, with 11% coverage.

The coefficients of `init_y` and `mean_*` describe the random effects,
not causal effects, and they receive no impacts.

## Effects over time

Because today’s outcome is in tomorrow’s model, a change in a regressor
keeps acting after it ends. The impacts of
[`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md)
(and `fit$impacts`) describe the period of the change only, the
**short-run** effect.
[`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md)
traces the rest. It simulates the outcome histories with and without the
change, using the same random numbers in both, and averages the
difference in the probabilities period by period over histories and
posterior draws.

A **temporary** change lasts one period (here, one unit of $`x_1`$ in
the first modelled period, for all units):

``` r

psp_response_path(fit, "x1", draws = 100, nsim = 10, seed = 1)
#> Response of P(y = 1) to a temporary shock of 1 in x1 (100 draws x 10 histories)
#> 
#> Average change over all units:
#>  horizon period   mean     sd   q2.5  q97.5
#>        0      2 0.2443 0.0313 0.1896 0.3052
#>        1      3 0.0757 0.0232 0.0348 0.1247
#>        2      4 0.0233 0.0128 0.0053 0.0569
#>        3      5 0.0081 0.0062 0.0009 0.0263
#>        4      6 0.0030 0.0030 0.0000 0.0092
#>        5      7 0.0010 0.0015 0.0000 0.0050
#>        6      8 0.0004 0.0006 0.0000 0.0022
```

Horizon 0 is the short-run total impact. The later horizons come only
from the lagged outcome: the change made the outcome more likely in
period 2, which makes it more likely in period 3, and so on, decaying at
a rate set by `lag_y` and the curvature of the probit.

A **permanent** change stays in place, so its effect builds up as the
outcome history adjusts. With `units =` only some units are changed, and
the effects on the changed units and on the others (the spillover) are
reported separately:

``` r

psp_response_path(fit, "x1", type = "permanent", units = c("1", "2", "3"),
                  draws = 100, nsim = 10, seed = 1)
#> Response of P(y = 1) to a permanent shock of 1 in x1 (100 draws x 10 histories)
#> 
#> Average change over all units:
#>  horizon period   mean     sd   q2.5  q97.5
#>        0      2 0.0116 0.0026 0.0071 0.0170
#>        1      3 0.0156 0.0031 0.0097 0.0218
#>        2      4 0.0155 0.0032 0.0097 0.0210
#>        3      5 0.0225 0.0043 0.0141 0.0298
#>        4      6 0.0232 0.0049 0.0139 0.0326
#>        5      7 0.0216 0.0054 0.0121 0.0325
#>        6      8 0.0143 0.0040 0.0067 0.0217
#> 
#> Average change over the shocked units:
#>  horizon period   mean     sd   q2.5  q97.5
#>        0      2 0.1539 0.0398 0.0866 0.2363
#>        1      3 0.1885 0.0372 0.1253 0.2697
#>        2      4 0.1928 0.0430 0.1164 0.2743
#>        3      5 0.3034 0.0600 0.1810 0.4120
#>        4      6 0.2947 0.0600 0.1679 0.4113
#>        5      7 0.2566 0.0657 0.1284 0.3749
#>        6      8 0.1653 0.0526 0.0591 0.2590
#> 
#> Average change over the other units:
#>  horizon period   mean     sd   q2.5  q97.5
#>        0      2 0.0024 0.0008 0.0012 0.0040
#>        1      3 0.0043 0.0017 0.0016 0.0082
#>        2      4 0.0039 0.0017 0.0015 0.0077
#>        3      5 0.0042 0.0018 0.0016 0.0090
#>        4      6 0.0055 0.0023 0.0024 0.0106
#>        5      7 0.0063 0.0028 0.0026 0.0113
#>        6      8 0.0044 0.0022 0.0015 0.0106
```

For the changed units the effect grows over the first periods as the
outcome history adjusts. It then varies from period to period, because
the slope of the probit depends on where the units’ probabilities sit,
which their regressors in each period determine. The effect on the other
units is small: only three units changed, and only their neighbours
respond.

`start =` sets the period of the change, and `$draws` holds the
posterior draws of every horizon, for derived quantities such as
cumulative effects.

## Persistent shocks: latent dynamics

`psp_dynamic(latent = TRUE)` puts the persistence in the latent outcome
instead (Wang and Kockelman 2009):

``` math
z_t = \rho W z_t + \gamma z_{t-1} + X_t\beta + \alpha + \varepsilon_t,
```

with $`\gamma`$ reported as `lag_z`. The outcome does not feed back;
what persists is the unobserved propensity, including the part of last
period’s shock that the outcome did not reveal. Every period is modelled
(the first one has no lag and shares $`\alpha`$ with the others), and
probabilities and impacts integrate over the latent history: the
variance of the latent outcome grows with $`t`$ as the lagged shocks
accumulate.

``` r

d_lat <- sim_psprobit(W, T = 6, beta = c(-0.2, 1, -0.8), rho = 0.3, sigma_alpha = 0.5,
                      gamma_z = 0.5, seed = 4)
fit_lat <- psprobit(y ~ x1 + x2, data = d_lat, index = c("id", "time"), W = W,
                    dynamic = psp_dynamic(latent = TRUE),
                    control = psp_control(iter = 8000, warmup = 2000, thin = 4,
                                          chains = 2, seed = 4, progress = FALSE))
round(coef(fit_lat), 2)
#> (Intercept)          x1          x2       lag_z         rho sigma_alpha 
#>       -0.34        1.12       -0.94        0.51        0.27        0.87
```

The sampler updates the latent outcomes given $`\gamma`$ and $`\gamma`$
given the latent outcomes, and the two are strongly dependent, so latent
dynamics need longer chains than the other models (above, four times the
usual length, thinned). Watch the convergence warning, and see the
“Convergence” section of
[`?psprobit`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md).

[`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md)
works here too; it computes the paths exactly from the reduced form,
without simulating outcome histories. Latent dynamics need
time-invariant weights and cannot be combined with lagged outcomes.

## Which kind of persistence?

The three mechanisms produce similar outcome patterns, so the data
separate them only through the model’s structure. Some guidance:

- If theory says the outcome itself changes the conditions for the next
  period (adoption of a technology, participation, a sunk entry cost),
  use `lag_y`, with random effects and the Wooldridge initial condition.
- If the persistence is in slowly moving unobserved conditions, use
  latent dynamics.
- Always include random effects: without them, both kinds of dynamics
  absorb the heterogeneity.

[`psp_ppc()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_ppc.md)
reports, for each period, the share of units whose outcome repeats from
the previous period. A model whose replicated persistence falls short of
the observed one (posterior predictive p-values near 0) is missing a
source of persistence. Here we compare the dynamic model with two static
models fitted to the same periods, one with random effects and one
pooled:

``` r

persistence <- function(f) {
  p <- psp_ppc(f, draws = 200, seed = 1)$persistence
  round(cbind(observed = p$observed, replicated = colMeans(p$replicated), ppp = p$ppp), 2)
}
d_est <- d[d$time > 1, ]
fit_static <- psprobit(y ~ x1 + x2, data = d_est, index = c("id", "time"), W = W,
                       impacts = FALSE, control = ctrl)
fit_static_pooled <- psprobit(y ~ x1 + x2, data = d_est, index = c("id", "time"), W = W,
                              effects = "pooled", impacts = FALSE, control = ctrl)
persistence(fit)
#>             observed replicated  ppp
#> persist t=3     0.57       0.60 0.72
#> persist t=4     0.59       0.54 0.24
#> persist t=5     0.65       0.59 0.18
#> persist t=6     0.63       0.58 0.28
#> persist t=7     0.63       0.58 0.28
#> persist t=8     0.69       0.63 0.22
persistence(fit_static)
#>             observed replicated  ppp
#> persist t=3     0.57       0.59 0.64
#> persist t=4     0.59       0.55 0.28
#> persist t=5     0.65       0.57 0.16
#> persist t=6     0.63       0.58 0.25
#> persist t=7     0.63       0.58 0.26
#> persist t=8     0.69       0.60 0.15
persistence(fit_static_pooled)
#>             observed replicated  ppp
#> persist t=3     0.57       0.52 0.25
#> persist t=4     0.59       0.50 0.12
#> persist t=5     0.65       0.53 0.04
#> persist t=6     0.63       0.50 0.06
#> persist t=7     0.63       0.50 0.06
#> persist t=8     0.69       0.53 0.01
```

The pooled static model has no source of persistence and falls short in
most periods. The static random-effects model reproduces the persistence
about as well as the dynamic model, through the random effects. This is
the point made at the start: the outcome patterns of state dependence
and heterogeneity are similar, and a check like this one cannot choose
between them. The choice rests on the substantive question and on
handling the initial condition properly.

## Restrictions

- Dynamic models need a balanced panel with complete outcomes (missing
  outcomes would make the lagged regressors missing).
- `predict(newdata = )` is not available with lagged outcomes, because
  new regressors would change the outcome history. Use
  [`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md).
- The impacts of
  [`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md)
  are short-run for lagged outcomes. For latent dynamics they are the
  effects of a change in one period on that period’s probabilities,
  integrating over the latent history.

## References

Wang, X. and Kockelman, K. M. (2009). Bayesian inference for ordered
response data with a dynamic spatial-ordered probit model. *Journal of
Regional Science* 49, 877–913.

Wooldridge, J. M. (2005). Simple solutions to the initial conditions
problem in dynamic, nonlinear panel data models with unobserved
heterogeneity. *Journal of Applied Econometrics* 20, 39–54.
