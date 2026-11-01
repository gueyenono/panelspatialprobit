# Ordered outcomes

Many spatial panel outcomes are ordered categories: credit ratings,
survey scales, damage levels, land-use intensity.
[`psoprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
fits the spatial ordered probit model. It uses the same arguments as
[`psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
(see
[`vignette("getting-started")`](https://gueyenono.github.io/panelspatialprobit/articles/getting-started.md)),
so this vignette concentrates on what is specific to ordered outcomes.

``` r

library(panelspatialprobit)
```

## The model

The latent model is the same as in the binary case,

``` math
z_t = \rho W z_t + X_t\beta + W X_t\theta + \alpha + \varepsilon_t,
\qquad \varepsilon_t \sim N(0, I_N),
```

and the observed category is $`y_{it} = j`$ when
$`\kappa_{j-1} < z_{it} \le \kappa_j`$, for $`j = 1, \dots, J`$, with
$`\kappa_0 = -\infty`$ and $`\kappa_J = +\infty`$.

**Identification.** The model includes an intercept, so the first
cut-point is fixed at $`\kappa_1 = 0`$, and the free cut-points are
$`0 < \kappa_2 < \dots < \kappa_{J-1}`$ (named `cut[2|3]`, `cut[3|4]`,
…). This is the parameterisation of `MCMCpack::MCMCoprobit`. Packages
such as [`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html) and
`ordinal::clm` drop the intercept and estimate all $`J - 1`$ thresholds
instead; [`summary()`](https://rdrr.io/r/base/summary.html) reports the
thresholds in that parameterisation too,
$`\zeta_j = \kappa_j - \beta_0`$.

The response can be an ordered factor, a factor (its level order is
used) or integer codes. Every category must occur in the data.

## Example

We simulate $`J = 4`$ categories on a 10 × 10 grid observed over 5
periods:

``` r

W <- psp_weights(rook_weights(10, 10), style = "W")
d <- sim_psprobit(W, T = 5, beta = c(0.2, 1, -0.8), rho = 0.5, sigma_alpha = 0.6,
                  cutpoints = c(0.8, 1.6), seed = 7)
d$y <- factor(d$y, labels = c("none", "low", "medium", "high"), ordered = TRUE)
table(d$y)
#> 
#>   none    low medium   high 
#>    216     80     68    136
```

``` r

fit <- psoprobit(y ~ x1 + x2, data = d, index = c("id", "time"), W = W,
                 control = psp_control(iter = 3000, warmup = 1000, chains = 2,
                                       seed = 1, progress = FALSE))
summary(fit)
#> Bayesian panel spatial ordered probit (model = sar, effects = random) 
#> Call: psoprobit(formula = y ~ x1 + x2, data = d, index = c("id", "time"), 
#>     W = W, control = psp_control(iter = 3000, warmup = 1000, 
#>         chains = 2, seed = 1, progress = FALSE))
#> Units: 100, periods: 5, observations: 500
#> Chains: 2, posterior draws: 6000
#> Outcome distribution:
#> 
#>   none    low medium   high 
#>    216     80     68    136 
#> 
#> Posterior summary:
#>                    mean    sd   q2.5  q97.5  rhat ess_bulk ess_tail
#> (Intercept)       0.122 0.095 -0.062  0.315 1.001      841     1219
#> x1                1.076 0.087  0.909  1.248 1.002     1209     2415
#> x2               -0.932 0.086 -1.102 -0.766 1.005     1144     2634
#> rho               0.511 0.048  0.414  0.601 1.000     1398     2839
#> sigma_alpha       0.744 0.111  0.538  0.968 1.003      586      915
#> cut[low|medium]   0.861 0.089  0.698  1.044 1.001      821     1347
#> cut[medium|high]  1.711 0.127  1.467  1.963 1.003      817     1540
#> 
#> Thresholds in the MASS::polr parameterisation (no intercept, zeta_j = cut_j - (Intercept)):
#>               mean    sd   q2.5 q97.5  rhat ess_bulk ess_tail
#> none|low    -0.122 0.095 -0.315 0.062 1.001      841     1219
#> low|medium   0.739 0.118  0.517 0.976 1.000      985     1809
#> medium|high  1.589 0.144  1.311 1.871 1.001      930     1710
#> 
#> In-sample classification: hit rate 0.614
#> 
#> Impacts on category probabilities (posterior mean, sd, 95% interval; 500 draws)
#> 
#> Category none
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1 -0.212 0.011 -0.230 -0.190
#> x2  0.185 0.011  0.163  0.206
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1 -0.189 0.032 -0.252 -0.136
#> x2  0.165 0.028  0.116  0.222
#> 
#> total:
#>     mean    sd   q2.5  q97.5
#> x1 -0.40 0.034 -0.466 -0.339
#> x2  0.35 0.032  0.292  0.413
#> 
#> Category low
#> 
#> direct:
#>      mean    sd   q2.5 q97.5
#> x1  0.001 0.007 -0.011 0.015
#> x2 -0.001 0.006 -0.014 0.011
#> 
#> indirect:
#>      mean    sd   q2.5 q97.5
#> x1  0.001 0.006 -0.011 0.014
#> x2 -0.001 0.006 -0.012 0.010
#> 
#> total:
#>      mean    sd   q2.5 q97.5
#> x1  0.002 0.013 -0.023 0.030
#> x2 -0.002 0.012 -0.025 0.021
#> 
#> Category medium
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1  0.031 0.007  0.018  0.043
#> x2 -0.027 0.006 -0.038 -0.016
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1  0.027 0.007  0.015  0.043
#> x2 -0.024 0.006 -0.037 -0.013
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1  0.058 0.013  0.033  0.083
#> x2 -0.050 0.011 -0.073 -0.030
#> 
#> Category high
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1  0.180 0.013  0.154  0.204
#> x2 -0.157 0.012 -0.180 -0.134
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1  0.161 0.029  0.113  0.223
#> x2 -0.140 0.026 -0.198 -0.097
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1  0.341 0.036  0.277  0.416
#> x2 -0.298 0.032 -0.367 -0.240
```

The truth is $`\beta = (0.2, 1, -0.8)`$, $`\rho = 0.5`$,
$`\sigma_\alpha = 0.6`$ and cut-points $`(0, 0.8, 1.6)`$.

## Impacts per category

A regressor moves probability mass between categories: if it raises the
latent outcome, the top categories become more likely and the bottom
ones less likely. Impacts are therefore reported **per category** and,
for each regressor and effect type, sum to zero over the categories.

``` r

imp <- fit$impacts
s <- summary(imp)
s[s$variable == "x1" & s$effect == "total", c("category", "mean", "q2.5", "q97.5")]
#>    category         mean        q2.5       q97.5
#> 17     none -0.400464702 -0.46609936 -0.33897627
#> 19      low  0.002032479 -0.02273775  0.02978470
#> 21   medium  0.057716139  0.03259027  0.08328227
#> 23     high  0.340716084  0.27657645  0.41639009
sum(s$mean[s$variable == "x1" & s$effect == "total"])
#> [1] -4.336809e-19
```

The options of
[`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md)
(`type = "unit"`, `by_period = TRUE`, `scale = "lesage_pace"`) work as
for binary models.

## Predicted probabilities

`predict(type = "prob")` returns an $`NT \times J`$ matrix of posterior
mean category probabilities, in the row order of the data;
`type = "class"` returns the most probable category.

``` r

head(round(predict(fit), 3))
#>       none   low medium  high
#> [1,] 0.015 0.038  0.092 0.855
#> [2,] 0.432 0.226  0.177 0.165
#> [3,] 0.449 0.229  0.173 0.149
#> [4,] 0.761 0.142  0.067 0.031
#> [5,] 0.327 0.229  0.208 0.236
#> [6,] 0.096 0.140  0.208 0.556
psp_classification(fit)$confusion
#>         predicted
#> observed none low medium high
#>   none    194   0      0   22
#>   low      56   0      0   24
#>   medium   32   0      0   36
#>   high     23   0      0  113
```

## The cut-point sampler

Updating the cut-points one at a time given the latent values (Albert
and Chib 1993) mixes slowly when there are many observations, because
each cut-point is squeezed between the latent values of the neighbouring
categories. The default sampler (`cutpoint_sampler = "transform"`) adds
a joint Metropolis–Hastings move that shifts a cut-point and rescales
the latent values of the two adjacent categories with it (Liu and
Sabatti 2000), so the cut-point can move far in one step. Its step size
is tuned during warm-up; the acceptance rates after warm-up are stored
in the fit:

``` r

fit$cut_accept
#>       [,1]      [,2]
#> [1,] 0.233 0.3223333
#> [2,] 0.354 0.4400000
```

`cutpoint_sampler = "albert_chib"` uses the Gibbs step alone. Both
samplers target the same posterior and pass simulation-based
calibration; the transformation move gives a higher effective sample
size per second for the cut-points, by up to two orders of magnitude for
interior cut-points (see
[`vignette("validation")`](https://gueyenono.github.io/panelspatialprobit/articles/validation.md)).

## Priors on the cut-points

The default prior is $`N(0, 10^2)`$ on each free cut-point, restricted
to the ordering constraint. `psp_priors(cutpoints = prior_flat())` uses
a flat prior, as `MCMCpack::MCMCoprobit` does.

## References

Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
polychotomous response data. *Journal of the American Statistical
Association* 88, 669–679.

LeSage, J. P., Pace, R. K., Lam, N., Campanella, R. and Liu, X. (2011).
New Orleans business recovery in the aftermath of Hurricane Katrina.
*Journal of the Royal Statistical Society A* 174, 1007–1027.

Liu, J. S. and Sabatti, C. (2000). Generalised Gibbs sampler and
multigrid Monte Carlo for Bayesian computation. *Biometrika* 87,
353–369.
