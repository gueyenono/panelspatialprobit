# MCMC and computation settings

MCMC and computation settings

## Usage

``` r
psp_control(
  iter = 2000,
  warmup = 1000,
  thin = 1,
  chains = 4,
  cores = 1,
  seed = NULL,
  z_sweeps = 1,
  rho_range = c(-1, 1),
  lambda_range = c(-1, 1),
  lndet = c("auto", "eigen", "sparse"),
  cutpoint_sampler = c("transform", "albert_chib"),
  scale_move = TRUE,
  impacts_draws = 500,
  impacts_exact_max = 500,
  impacts_probes = 50,
  mixture_draws = 50,
  progress = interactive()
)
```

## Arguments

- iter:

  Number of post-warm-up iterations per chain.

- warmup:

  Number of warm-up iterations per chain (discarded).

- thin:

  Thinning interval for stored draws.

- chains:

  Number of chains.

- cores:

  Number of chains run at once, on separate threads (on every operating
  system, Windows included). Results do not depend on `cores`. Under
  `R CMD check` (`_R_CHECK_LIMIT_CORES_`) at most 2 are used.

- seed:

  Optional integer seed. Given `seed`, results are reproducible on a
  given platform, whatever `cores` and the user's
  [`RNGkind()`](https://rdrr.io/r/base/Random.html): each chain draws
  from its own L'Ecuyer-CMRG stream derived from `seed`, and the
  session's random number state is left unchanged. Without `seed`, one
  is drawn from the session's generator (`fit$seed`).

- z_sweeps:

  Number of Gibbs sweeps over the latent variables per iteration.

- rho_range:

  Admissible interval for the spatial autoregressive parameter.

- lambda_range:

  Admissible interval for the spatial error parameter.

- lndet:

  Method for log-determinants: `"auto"`, `"eigen"` or `"sparse"`.

- cutpoint_sampler:

  `"transform"` (Albert-Chib Gibbs step plus a joint transformation move
  of cut-points and latent values) or `"albert_chib"`.

- scale_move:

  Logical; add a generalised Gibbs move that rescales the latent values,
  coefficients, random effects and cut-points jointly (Liu and Sabatti,
  2000). It leaves the posterior unchanged and greatly improves the
  mixing of the coefficients.

- impacts_draws:

  Number of posterior draws used to compute impacts and predictions.

- impacts_exact_max:

  Largest number of units for which diagonals of matrix inverses are
  computed exactly; above it, a stochastic estimator is used.

- impacts_probes:

  Number of probe vectors of the stochastic estimator.

- mixture_draws:

  Number of simulated sets of mixing variances over which probabilities
  and impacts are averaged for Student-t errors
  ([`psp_t()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_t.md)).

- progress:

  Logical; print sampling progress, for each chain.

## Value

An object of class `psp_control`.

## Parallel post-processing

[`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md),
[predict()](https://gueyenono.github.io/panelspatialprobit/reference/psprobit-methods.md),
[`psp_ppc()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_ppc.md),
[`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md),
[`psp_loglik()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_loglik.md),
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`psp_dic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_loglik.md)
work through the posterior draws in blocks and follow the
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
set by the user: for example
`future::plan("multisession", workers = 4)`, or
`future::plan("multicore", workers = 4)` on Unix-alikes. The default
sequential plan runs in the current R session. Results do not depend on
the plan or the number of workers.

## Examples

``` r
psp_control(iter = 500, warmup = 250, chains = 2)
#> $iter
#> [1] 500
#> 
#> $warmup
#> [1] 250
#> 
#> $thin
#> [1] 1
#> 
#> $chains
#> [1] 2
#> 
#> $cores
#> [1] 1
#> 
#> $seed
#> NULL
#> 
#> $z_sweeps
#> [1] 1
#> 
#> $rho_range
#> [1] -1  1
#> 
#> $lambda_range
#> [1] -1  1
#> 
#> $lndet
#> [1] "auto"
#> 
#> $cutpoint_sampler
#> [1] "transform"
#> 
#> $scale_move
#> [1] TRUE
#> 
#> $impacts_draws
#> [1] 500
#> 
#> $impacts_exact_max
#> [1] 500
#> 
#> $impacts_probes
#> [1] 50
#> 
#> $mixture_draws
#> [1] 50
#> 
#> $progress
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "psp_control"
```
