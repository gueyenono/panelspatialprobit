# Approximate log-likelihood and DIC (experimental)

The likelihood of a spatial probit model is an \\NT\\-dimensional normal
rectangle probability with no closed form. `psp_loglik()` estimates it
with the GHK simulator (Geweke, Hajivassiliou and Keane), applied to the
sparse precision matrix of the latent outcomes so that it scales to
large \\N\\. For random-effects models the random effects are added to
the latent vector and integrated out with it (`type = "marginal"`), or
held at their posterior draws (`type = "conditional"`).

## Usage

``` r
psp_loglik(
  fit,
  draws = 100,
  reps = 200,
  type = c("marginal", "conditional"),
  seed = 1
)

# S3 method for class 'psprobit'
logLik(object, reps = 500, seed = 1, ...)

psp_dic(fit, draws = 100, reps = 200, seed = 1)
```

## Arguments

- fit, object:

  A fitted `psprobit` object.

- draws:

  Number of posterior draws to evaluate.

- reps:

  Number of GHK replications per evaluation.

- type:

  `"marginal"` (random effects integrated out) or `"conditional"` (given
  the random-effect draws kept by the sampler). The two coincide for
  pooled models.

- seed:

  Seed for the GHK replications (the same random numbers are used for
  every draw).

- ...:

  Unused.

## Value

`psp_loglik()`: a vector of log-likelihood values, one per draw, with
attribute `mcse`. [`logLik()`](https://rdrr.io/r/stats/logLik.html): a
`logLik` object. `psp_dic()`: a list with `dic`, `p_d` (effective number
of parameters), `d_bar` (posterior mean deviance) and `d_hat` (deviance
at the posterior mean).

## Details

The GHK estimate of the likelihood is unbiased, but its logarithm is
biased downwards when the replications vary a lot; increase `reps` until
the estimates settle, and compare models with the same `reps` and `seed`
(common random numbers). `attr(, "mcse")` gives a delta-method Monte
Carlo standard error for each value, and a warning is given when it
exceeds 0.5. Without random effects to integrate out (pooled models, or
`type = "conditional"`) the periods are independent and are simulated
separately, which keeps the error small even for large panels. The
marginal likelihood of a random-effects model is a single
\\NT\\-dimensional problem and becomes unreliable beyond a few hundred
observations.

[`logLik()`](https://rdrr.io/r/stats/logLik.html) evaluates the marginal
log-likelihood at the posterior mean of the parameters, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work. `psp_dic()` computes
the deviance information criterion (Spiegelhalter et al. 2002) from the
marginal likelihood.

The posterior draws are processed in parallel under the user's
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(see "Parallel post-processing" in
[`psp_control()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_control.md)).

## References

Geweke, J. (1991). Efficient simulation from the multivariate normal and
Student-t distributions subject to linear constraints. *Computer Science
and Statistics: Proceedings of the 23rd Symposium on the Interface*,
571-578.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
(2002). Bayesian measures of model complexity and fit. *JRSS B* 64,
583-639.

## Examples

``` r
W <- psp_weights(rook_weights(4, 4), style = "W")
d <- sim_psprobit(W, T = 3, beta = c(0, 1), rho = 0.4, seed = 2)
fit <- psprobit(y ~ x1, data = d, index = c("id", "time"), W = W, impacts = FALSE,
                control = psp_control(iter = 200, warmup = 100, chains = 1, seed = 2))
#> Warning: Possible convergence problems (R-hat > 1.01 or bulk ESS < 100) for: (Intercept), x1, rho, sigma_alpha. Consider more iterations.
logLik(fit, reps = 200)
#> 'log Lik.' -35.86657 (df=4)
psp_dic(fit, draws = 20, reps = 100)
#> Warning: The GHK estimates are unreliable (Monte Carlo se up to 1.00 on the log scale): increase 'reps'. GHK needs many replications for large N x T.
#> Warning: Negative p_D: the GHK error dominates; increase 'reps'.
#> $dic
#> [1] 61.75925
#> 
#> $p_d
#> [1] -4.844922
#> 
#> $d_bar
#> [1] 66.60417
#> 
#> $d_hat
#> [1] 71.44909
#> 
```
