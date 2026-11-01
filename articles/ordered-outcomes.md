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
#> (Intercept)       0.126 0.097 -0.060  0.319 1.002      737     1444
#> x1                1.072 0.088  0.903  1.249 1.003     1030     2184
#> x2               -0.933 0.086 -1.102 -0.763 1.001     1110     2569
#> rho               0.511 0.049  0.413  0.605 1.000     1202     2434
#> sigma_alpha       0.751 0.109  0.547  0.974 1.004      498      873
#> cut[low|medium]   0.868 0.091  0.701  1.053 1.002      638     1198
#> cut[medium|high]  1.720 0.126  1.481  1.973 1.001      705     1473
#> 
#> Thresholds in the MASS::polr parameterisation (no intercept, zeta_j = cut_j - (Intercept)):
#>               mean    sd   q2.5 q97.5  rhat ess_bulk ess_tail
#> none|low    -0.126 0.097 -0.319 0.060 1.002      737     1444
#> low|medium   0.742 0.121  0.501 0.985 1.002      882     1758
#> medium|high  1.594 0.147  1.303 1.892 1.001      980     1849
#> 
#> In-sample classification: hit rate 0.614
#> 
#> Impacts on category probabilities (posterior mean, sd, 95% interval; 500 draws)
#> 
#> Category none
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1 -0.211 0.011 -0.234 -0.188
#> x2  0.183 0.012  0.159  0.206
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1 -0.187 0.033 -0.254 -0.131
#> x2  0.162 0.027  0.110  0.217
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1 -0.397 0.036 -0.467 -0.336
#> x2  0.345 0.031  0.288  0.407
#> 
#> Category low
#> 
#> direct:
#>    mean    sd   q2.5 q97.5
#> x1    0 0.007 -0.014 0.015
#> x2    0 0.007 -0.013 0.013
#> 
#> indirect:
#>    mean    sd   q2.5 q97.5
#> x1    0 0.007 -0.014 0.013
#> x2    0 0.006 -0.011 0.012
#> 
#> total:
#>    mean    sd   q2.5 q97.5
#> x1    0 0.014 -0.028 0.027
#> x2    0 0.012 -0.025 0.025
#> 
#> Category medium
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1  0.030 0.007  0.017  0.043
#> x2 -0.026 0.006 -0.039 -0.014
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1  0.026 0.007  0.014  0.041
#> x2 -0.023 0.006 -0.035 -0.012
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1  0.057 0.013  0.031  0.081
#> x2 -0.049 0.012 -0.071 -0.027
#> 
#> Category high
#> 
#> direct:
#>      mean    sd   q2.5  q97.5
#> x1  0.180 0.013  0.154  0.205
#> x2 -0.157 0.013 -0.181 -0.132
#> 
#> indirect:
#>      mean    sd   q2.5  q97.5
#> x1  0.160 0.031  0.108  0.224
#> x2 -0.139 0.025 -0.191 -0.094
#> 
#> total:
#>      mean    sd   q2.5  q97.5
#> x1  0.340 0.038  0.279  0.420
#> x2 -0.296 0.032 -0.358 -0.235
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
#>    category          mean        q2.5       q97.5
#> 17     none -0.3970782465 -0.46653413 -0.33591238
#> 19      low  0.0004226489 -0.02845663  0.02705857
#> 21   medium  0.0565979051  0.03058806  0.08141688
#> 23     high  0.3400576925  0.27881770  0.41981368
sum(s$mean[s$variable == "x1" & s$effect == "total"])
#> [1] -3.122502e-17
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
#> [1,] 0.015 0.039  0.092 0.854
#> [2,] 0.426 0.227  0.179 0.168
#> [3,] 0.442 0.230  0.175 0.153
#> [4,] 0.752 0.146  0.070 0.033
#> [5,] 0.324 0.229  0.209 0.239
#> [6,] 0.097 0.141  0.208 0.555
psp_classification(fit)$confusion
#>         predicted
#> observed none low medium high
#>   none    194   0      0   22
#>   low      56   0      0   24
#>   medium   31   0      0   37
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
#>           [,1]      [,2]
#> [1,] 0.2343333 0.3120000
#> [2,] 0.3496667 0.4353333
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
