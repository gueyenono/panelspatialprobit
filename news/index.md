# Changelog

## panelspatialprobit 0.0.0.9000

- First implementation (plan phases 0-3):
  - Model fitting:
    - [`psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
      and
      [`psoprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/psprobit.md)
      for balanced panels.
    - SAR, SDM, SLX and non-spatial specifications.
    - Pooled or random effects, Mundlak correlated random effects and
      time effects.
  - Methods and helpers:
    - [`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md),
      with the
      [`spatialreg::impacts()`](https://r-spatial.github.io/spatialreg/reference/impacts.html)
      method, for direct, indirect and total effects on probabilities.
    - Methods: [`summary()`](https://rdrr.io/r/base/summary.html),
      [`predict()`](https://rdrr.io/r/stats/predict.html),
      [`coef()`](https://rdrr.io/r/stats/coef.html),
      [`vcov()`](https://rdrr.io/r/stats/vcov.html),
      [`confint()`](https://rdrr.io/r/stats/confint.html),
      [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
      [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
      [`glance()`](https://generics.r-lib.org/reference/glance.html) and
      `as_draws()`.
    - [`sim_psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/sim_psprobit.md),
      [`psp_weights()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_weights.md)
      and
      [`rook_weights()`](https://gueyenono.github.io/panelspatialprobit/reference/rook_weights.md).
  - Engine: a C++ sampler with residual-update latent sweeps, a
    collapsed griddy Gibbs step for rho, a half-Cauchy random-effects
    prior, a transformation move for the cut-points and a
    generalised-Gibbs scale move.
- Completion of phases 2-3:
  - Chamberlain correlated random effects:
    `cre = list(type = "chamberlain", vars = ~ x1)`.
  - [`psp_impacts()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_impacts.md)
    gains `type = "unit"` (impacts conditional on the random effects,
    using random-effect draws kept by the sampler) and
    `by_period = TRUE` (impacts within each period).
  - [`psp_classification()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_classification.md):
    in-sample hit rate, AUC, Brier score and confusion table; reported
    by [`summary()`](https://rdrr.io/r/base/summary.html).
  - [`summary()`](https://rdrr.io/r/base/summary.html) of ordered models
    also shows the thresholds in the
    [`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html)
    parameterisation.
  - `psp_control(progress = )` prints sampling progress.
  - Vignettes “Getting started with panel spatial probit” and “Ordered
    outcomes”.
- Spatial error models (phase 4): `model = "sem"`, `"sdem"` and
  `"sarar"`, with optional error weights `M` and a `lambda` prior
  (`psp_priors(lambda =)`, `psp_control(lambda_range =)`). The engine
  now works with the general operator `A = (I - lambda M)(I - rho W)`,
  and impacts and predictions account for the error term in the latent
  variance.
- Dynamic models (phase 5): `dynamic = psp_dynamic()` adds the lagged
  outcome (and optionally its spatial lag), with Wooldridge (2005)
  initial conditions by default.
  [`sim_psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/sim_psprobit.md)
  gains `gamma` and `eta` for lagged outcomes and `lambda` / `M` for
  spatial errors.
- Missing outcomes are imputed through their latent values (phase 6).
- `predict(newdata = )` gives counterfactual predictions for the fitted
  units and periods (phase 6).
- [`psp_ppc()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_ppc.md):
  posterior predictive checks of category shares, Moran’s I and
  persistence per period (phase 6).
- The griddy Gibbs step for `rho` and `lambda` locates the posterior on
  a coarse pass before refining, about 5x faster for small panels.
- Time-varying spatial weights (phase 6): `W` (and `M`) may be a list
  with one matrix per period. Durbin terms, the spatial lag of the
  lagged outcome, impacts, predictions and posterior predictive checks
  use the weights of each period.
- Spatially autocorrelated random effects (phase 4):
  `effects = "spatial_random"` gives `alpha = phi W alpha + a`, with a
  `phi` prior in `psp_priors(phi = )`.
  [`sim_psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/sim_psprobit.md)
  gains `phi`.
- [`psp_response_path()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_response_path.md)
  (phase 5): how a temporary or permanent change in a regressor
  propagates over time in a dynamic model, by simulating the outcome
  histories under the baseline and the shock with common random numbers.
- Experimental approximate likelihood (phase 6):
  [`psp_loglik()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_loglik.md),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (and so
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) /
  [`BIC()`](https://rdrr.io/r/stats/AIC.html)) and
  [`psp_dic()`](https://gueyenono.github.io/panelspatialprobit/reference/psp_loglik.md).
  The likelihood is estimated with the GHK simulator on the sparse
  precision of the latent outcomes, with the random effects integrated
  out.
- Student-t latent errors (phase 4): `errors = psp_t(df)` (robit for
  binary outcomes), sampled as a scale mixture of normals with fixed
  degrees of freedom. Probabilities and impacts average over simulated
  mixing variances (`psp_control(mixture_draws = )`); `fit$v_mean` flags
  outlying observations.
  [`sim_psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/sim_psprobit.md)
  gains `df`.
- Latent AR(1) dynamics (phase 5): `psp_dynamic(latent = TRUE)` gives
  `z_t = rho W z_t + gamma z_{t-1} + X_t beta + alpha + u_t`
  (coefficient `lag_z`), with probabilities, impacts, predictions,
  response paths and the approximate likelihood integrating over the
  latent history.
  [`sim_psprobit()`](https://gueyenono.github.io/panelspatialprobit/reference/sim_psprobit.md)
  gains `gamma_z`.
- New vignettes (phase 7): “Interpreting impacts”, “Dynamic panel
  spatial probit” and “How we validated the samplers” (reference
  sampler, SBC, Monte Carlo study, a replication of the Wang and
  Kockelman (2009) simulation, and cross-checks; drawn from summaries
  written by `validation/vignette_data.R`).
- [`summary()`](https://rdrr.io/r/base/summary.html) prints effective
  sample sizes as whole numbers.
