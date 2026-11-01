# 3. Package architecture

## 3.1 User-facing API

```r
fit <- psprobit(                         # binary;  psoprobit() for ordered (same signature)
  formula      = y ~ x1 + x2 + x3,
  data         = df,                     # data.frame, tibble, or plm::pdata.frame
  index        = c("id", "year"),        # unit and time identifiers (as in plm/splm)
  W            = W,                      # N x N (matrix, dgCMatrix, spdep listw/nb) or list of T matrices
  M            = NULL,                   # error-process weights; defaults to W
  model        = "sar",                  # "sar","sdm","sem","sarar","sdem","slx","none"
  effects      = "random",               # "pooled","random","spatial_random"
  cre          = NULL,                   # ~ x1 + x2 (Mundlak means) or list(type = "chamberlain", ...)
  durbin       = NULL,                   # ~ x1 + x2 : WX terms (implied for "sdm"/"sdem"/"slx" if NULL = all)
  time_effects = FALSE,
  dynamic      = NULL,                   # psp_dynamic(lag_y = TRUE, lag_Wy = FALSE, initial = "wooldridge")
  priors       = psp_priors(),           # beta = prior_normal(0, 10), rho = prior_beta(1.01, 1.01), ...
  control      = psp_control(),          # iter, warmup, thin, chains, cores, seed, z_sweeps, lndet, ...
  impacts      = TRUE                    # compute impacts after sampling (can be redone with impacts())
)
```

- **Two front-ends, one engine.** `psprobit()` and `psoprobit()` validate and build the design, then call the internal engine `psp_fit()`. A low-level matrix interface, `psp_fit_matrix(y, X, W, N, T, ...)`, is exported for power users and simulations, in the same spirit as `sar_probit_mcmc()`.
- `psoprobit()` accepts an ordered factor or integer codes $1..J$. It errors on empty categories and when $J = 2$, pointing the user to `psprobit()`.
- **Constructors** keep the call readable and validate early:
  - `psp_priors()` with helpers `prior_normal()`, `prior_beta()`, `prior_half_cauchy()`, `prior_inv_gamma()`;
  - `psp_control()`;
  - `psp_dynamic()`.

### Methods

| Method | Notes |
|---|---|
| `print`, `summary` (+ `print.summary`) | Coefficient table with mean, sd, 95% CrI, $\hat R$, bulk/tail ESS; model, effects and data dimensions; prevalence or category table. |
| `coef`, `vcov`, `confint`, `nobs`, `formula`, `model.matrix` | Posterior means, posterior covariance, credible intervals. |
| `impacts()` | S3 method registered on `spatialreg::impacts` with **delayed registration** (`S3method(spatialreg::impacts, psprobit)` in NAMESPACE, R ≥ 3.6), so it works whether or not `spatialreg` is attached. Also exported as `psp_impacts()` to avoid masking. |
| `fitted`, `predict(type = c("link", "prob", "class"), newdata, marginal = TRUE)` | In-sample first (Phase 2); `newdata` for the same units/periods with new X (Phase 6). |
| `simulate` | Posterior predictive $y^{rep}$ (used by PPCs). |
| `plot` | Trace and density plots; uses `bayesplot` if installed, otherwise base graphics. |
| `as_draws()` (`posterior`), `as.mcmc.list()` (`coda`) | Interoperability. |
| `tidy()`, `glance()` (`generics`) | broom-style output. |
| `logLik` | Phase 6, approximate, with an attribute recording the approximation used. |
| `update` | Standard. |

### Simulation and weights utilities (exported)

- `sim_psprobit(N, T, W, X = NULL, beta, rho, lambda, theta, sigma_alpha, cutpoints, model, effects, dynamic, seed)`: returns a data.frame with the true parameters as an attribute. It is the backbone of tests, SBC, Monte Carlo and vignettes.
- `psp_weights(W, ids = NULL, style = c("asis", "W"))`:
  - coerces matrix / `dgCMatrix` / `listw` / `nb` to a sparse `dgCMatrix`;
  - checks that $W$ is square, has a zero diagonal and contains no NA;
  - optionally row-standardises;
  - detects isolates;
  - computes the admissible range of $\rho$;
  - attaches unit IDs.

## 3.2 Data pipeline (`model_frame.R`)

1. Evaluate `formula` on `data`; extract `index` into unit and time keys. Error on duplicated (id, time) pairs.
2. Check balance. Phase 1–5 requires a balanced panel with a clear error message; Phase 6 adds missing-$y$ support.
3. **Match units to $W$.** Units are matched by `rownames(W)` or `psp_weights(ids = )`. If $W$ has no names, require the data's unit IDs to be exactly the sorted unique IDs and warn. Silent mis-alignment of $W$ and data is the most likely user error, so this step is strict.
4. Sort internally by (time, unit-in-$W$-order) and keep the permutation, so fitted values and predictions come back in the user's row order.
5. Build the regressor columns:
   - $X$ via `model.matrix`;
   - Durbin columns $W_tX_t$ per period;
   - CRE means or Chamberlain columns;
   - time dummies;
   - dynamic lags, dropping the first period.
6. Check the rank of $X$ and drop or flag aliased columns, as `lm` does. Typical causes are intercept + time dummies or time-invariant CRE means.
7. Code $y$: binary accepts 0/1, logical or a two-level factor; ordered maps to $1..J$ and keeps the labels.

## 3.3 Internal layout

```
DESCRIPTION / NAMESPACE / NEWS.md / README.Rmd / LICENSE
R/
  psprobit.R          # front-ends psprobit(), psoprobit()
  psp_fit.R           # engine: chains, seeds, parallel, assembling the fit object
  model_frame.R       # §3.2
  weights.R           # psp_weights(), W validation/coercion
  logdet.R            # grids: eigen, sparse LU, Chebyshev, MC; spline interpolation
  priors.R  control.R dynamic.R
  impacts.R           # §2.11 + psp_impacts class/methods
  methods.R summary.R predict.R simulate.R plot.R tidy.R
  diagnostics.R ppc.R loglik_approx.R
  reference_sampler.R # pure-R, dense, slow ORACLE used only in tests (internal)
  utils.R  zzz.R      # .onLoad delayed S3 registration
  000-wrappers.R      # generated by savvy-cli
  rust.R              # conversions for Rust, .psp_post_map()
src/                  # the Rust engine since 2026-10-31; layout in 06 §6.4
  Makevars.in, Makevars.win.in, init.c   # build (configure) and savvy registration
  rust/pspcore/       # numerics: rng, dist, sparse, sampler, rf, post, simulate, ghk
  rust/pspr/          # savvy bindings (psp_rs_* functions)
tests/testthat/       # fast CRAN-safe tests (see 04)
vignettes/            # getting-started, ordered, impacts, validation summary
inst/CITATION
validation/           # NOT in the build: SBC, Monte Carlo, benchmarks, replication scripts + results
plan/  literature/  related_packages/   # NOT in the build
```

**Key design choices**

- **The full MCMC loop runs in compiled code** (one call runs every chain, on threads). This avoids per-iteration R overhead when $n$ is large, and each chain keeps a single RNG stream. Model variants are runtime flags. *The engine was C++ (RcppEigen) until 2026-10-31 and is now Rust (`06`); the last C++ version is tagged `cpp-final`.*
- **The pure-R reference sampler** (`reference_sampler.R`) is written first, in Phase 1. It is a direct, readable transcription of `02` using dense algebra, and it is the oracle for the compiled engine. Seeds cannot match across implementations, so the two are compared **distributionally**: both pass SBC, and long-run posterior summaries agree on fixed datasets.
- Sparse matrices cross into Rust as CSC arrays built from `dgCMatrix` (`.psp_csc()`); sparse LU and Cholesky come from `faer` (`06` §6.8). *Decision D3 (RcppEigen) is superseded by `06` R-D1 (savvy) and the Rust port.*
- **Toolchain:** Rust ≥ 1.84 (`SystemRequirements`), built by `configure` and `src/Makevars.in`; crates vendored for CRAN (`06` §6.9). No OpenMP.

## 3.4 Dependencies

| Field | Packages | Why |
|---|---|---|
| Depends | R (≥ 4.1) | Native pipe in docs/examples; delayed S3 registration. |
| Imports | `Matrix`, `future`, `stats`, `methods`, `parallel`, `posterior`*, `generics` | Sparse matrices, parallel post-processing, RNG streams, diagnostics, tidy generics. |
| SystemRequirements | Cargo, rustc ≥ 1.84 | The Rust engine (`06`); no LinkingTo since the C++ engine was removed. |
| Suggests | `testthat` (≥ 3.0.0), `spdep`, `spatialreg`, `sf`, `coda`, `bayesplot`, `knitr`, `rmarkdown`, `spatialprobit`, `MCMCpack`, `tmvtnorm`, `withr` | Weights helpers, test oracles, cross-package checks, plots, vignettes. |

\* `posterior` in Imports or Suggests is decision D4. Without it, $\hat R$ and ESS would need about 100 lines of our own code.

**`SBC`** (the Simulation-Based Calibration package) is **not on CRAN**, so it is used only in `validation/`, never in DESCRIPTION.

## 3.5 Repository and tooling conventions

- `usethis` scaffolding, `roxygen2` with markdown, testthat 3rd edition, snapshot tests for print/summary.
- **Style:** tidyverse style (`styler`), `lintr` config; snake_case. Internal R helpers are prefixed with `.psp_`; Rust functions exposed to R are prefixed `psp_rs_`; `cargo fmt` and `cargo clippy -D warnings`.
- **GitHub Actions:**
  - `R-CMD-check` on ubuntu / macOS / windows × (release, devel, oldrel-1);
  - `test-coverage` (covr);
  - `pkgdown`;
  - a manual or nightly `validation` workflow that runs the SBC smoke subset.
- **`.Rbuildignore`:** `plan/`, `literature/`, `related_packages/`, `validation/`, `CLAUDE.md`, `.github/`, `README.Rmd`.
- **`.gitignore`:**
  - compiled objects (`src/*.o`, `src/*.so`, `src/*.dll`);
  - `literature/*.pdf`, because publisher and working-paper PDFs should not be redistributed in a public repo; keep `literature/README.md`;
  - `related_packages/*/` sources and tarballs, which can be re-downloaded and are third-party; keep `related_packages/README.md`.
- **Synced-folder caveat:** the project currently sits in a Google Drive folder. Compiling there makes the sync client churn on `.o`/`.so` files and can cause file-lock errors during `R CMD INSTALL`. Recommended: keep the git working copy on a local disk with a GitHub remote (decision D5).
- **Lifecycle:** functions added after Phase 3 start as `experimental` (`lifecycle` badges in docs only; no runtime dependency).
- **License:** GPL (≥ 3) (decision D2). This is compatible with porting ideas or code from `spatialprobit` (GPL ≥ 2), `ProbitSpatial` (GPL ≥ 2) and `spmixW` (GPL ≥ 3). Any ported code keeps its original copyright notice in `inst/COPYRIGHTS` and in DESCRIPTION `Authors@R` (`role = "cph"`).
