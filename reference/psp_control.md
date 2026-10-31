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

  Number of cores used to run chains in parallel (forking, so only
  effective on Unix-alikes). Results do not depend on `cores`.

- seed:

  Optional integer seed. Results are reproducible given `seed`.

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

  Logical; print sampling progress (only when chains run sequentially,
  i.e. `cores = 1`).

## Value

An object of class `psp_control`.

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
