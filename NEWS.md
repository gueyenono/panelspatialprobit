# panelspatialprobit 0.0.0.9000

* First implementation (plan phases 0-3):
  * Model fitting:
    * `psprobit()` and `psoprobit()` for balanced panels.
    * SAR, SDM, SLX and non-spatial specifications.
    * Pooled or random effects, Mundlak correlated random effects and time
      effects.
  * Methods and helpers:
    * `psp_impacts()`, with the `spatialreg::impacts()` method, for direct,
      indirect and total effects on probabilities.
    * Methods: `summary()`, `predict()`, `coef()`, `vcov()`, `confint()`,
      `plot()`, `tidy()`, `glance()` and `as_draws()`.
    * `sim_psprobit()`, `psp_weights()` and `rook_weights()`.
  * Engine: a sampler with residual-update latent sweeps, a collapsed
    griddy Gibbs step for rho, a half-Cauchy random-effects prior, a
    transformation move for the cut-points and a generalised-Gibbs scale move.
* Completion of phases 2-3:
  * Chamberlain correlated random effects:
    `cre = list(type = "chamberlain", vars = ~ x1)`.
  * `psp_impacts()` gains `type = "unit"` (impacts conditional on the random
    effects, using random-effect draws kept by the sampler) and
    `by_period = TRUE` (impacts within each period).
  * `psp_classification()`: in-sample hit rate, AUC, Brier score and confusion
    table; reported by `summary()`.
  * `summary()` of ordered models also shows the thresholds in the
    `MASS::polr` parameterisation.
  * `psp_control(progress = )` prints sampling progress.
  * Vignettes "Getting started with panel spatial probit" and "Ordered
    outcomes".
* Spatial error models (phase 4): `model = "sem"`, `"sdem"` and `"sarar"`,
  with optional error weights `M` and a `lambda` prior (`psp_priors(lambda =)`,
  `psp_control(lambda_range =)`). The engine now works with the general
  operator `A = (I - lambda M)(I - rho W)`, and impacts and predictions account
  for the error term in the latent variance.
* Dynamic models (phase 5): `dynamic = psp_dynamic()` adds the lagged outcome
  (and optionally its spatial lag), with Wooldridge (2005) initial conditions
  by default. `sim_psprobit()` gains `gamma` and `eta` for lagged outcomes and
  `lambda` / `M` for spatial errors.
* Missing outcomes are imputed through their latent values (phase 6).
* `predict(newdata = )` gives counterfactual predictions for the fitted units
  and periods (phase 6).
* `psp_ppc()`: posterior predictive checks of category shares, Moran's I and
  persistence per period (phase 6).
* The griddy Gibbs step for `rho` and `lambda` locates the posterior on a coarse
  pass before refining, about 5x faster for small panels.
* Time-varying spatial weights (phase 6): `W` (and `M`) may be a list with one
  matrix per period. Durbin terms, the spatial lag of the lagged outcome,
  impacts, predictions and posterior predictive checks use the weights of
  each period.
* Spatially autocorrelated random effects (phase 4):
  `effects = "spatial_random"` gives `alpha = phi W alpha + a`, with a `phi`
  prior in `psp_priors(phi = )`. `sim_psprobit()` gains `phi`.
* `psp_response_path()` (phase 5): how a temporary or permanent change in a
  regressor propagates over time in a dynamic model, by simulating the outcome
  histories under the baseline and the shock with common random numbers.
* Experimental approximate likelihood (phase 6): `psp_loglik()`, `logLik()`
  (and so `AIC()` / `BIC()`) and `psp_dic()`. The likelihood is estimated with
  the GHK simulator on the sparse precision of the latent outcomes, with the
  random effects integrated out.
* Student-t latent errors (phase 4): `errors = psp_t(df)` (robit for binary
  outcomes), sampled as a scale mixture of normals with fixed degrees of
  freedom. Probabilities and impacts average over simulated mixing variances
  (`psp_control(mixture_draws = )`); `fit$v_mean` flags outlying observations.
  `sim_psprobit()` gains `df`.
* Latent AR(1) dynamics (phase 5): `psp_dynamic(latent = TRUE)` gives
  `z_t = rho W z_t + gamma z_{t-1} + X_t beta + alpha + u_t` (coefficient
  `lag_z`), with probabilities, impacts, predictions, response paths and the
  approximate likelihood integrating over the latent history.
  `sim_psprobit()` gains `gamma_z`.
* New vignettes (phase 7): "Interpreting impacts", "Dynamic panel spatial
  probit" and "How we validated the samplers" (reference sampler, SBC, Monte
  Carlo study, a replication of the Wang and Kockelman (2009) simulation, and
  cross-checks; drawn from summaries written by
  `validation/vignette_data.R`).
* `summary()` prints effective sample sizes as whole numbers.
* The engine is written in Rust (it was C++ with RcppEigen). Installing from
  source needs Cargo and rustc >= 1.84 (<https://rustup.rs>). Changes:
  * Chains run in parallel threads on every operating system, Windows
    included (`psp_control(cores = )`); the draws do not depend on `cores`.
  * New `psp_control(threads_per_chain = )`: the sweep over the latent
    variables updates the periods of a chain in parallel, each period with its
    own random number substream, so the draws do not depend on it either.
    With 4 threads, a binary SAR fit with N = 10,000 and T = 20 runs 3 times
    as fast as with the C++ engine.
  * The sampler is faster on one core: per iteration, 1.3 to 1.5 times faster
    than the C++ engine at N = 1,000 and 1.4 to 2 times at N = 10,000 (SEM
    models gain most, from a supernodal Cholesky factor of the random-effects
    precision).
  * `psp_impacts()`, `predict()`, `psp_ppc()`, `psp_response_path()` and
    `psp_loglik()` run in Rust (1.3 to 18 times faster on one core) and in parallel over the
    workers of `future::plan()` (see `?psp_control`); their results do not
    depend on the plan.
  * The sampler always uses inversion normals, whatever
    `RNGkind(normal.kind = )` says.
  * For the same seed, fits, `psp_ppc()`, `psp_response_path()` and
    `psp_loglik()` give different (still reproducible) numbers: the latent
    sweep draws from one substream per period, and the post-processing from
    package streams rather than R's generator.
