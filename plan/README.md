# panelspatialprobit: development plan

**Goal:** an R package for **Bayesian estimation of panel spatial probit and panel spatial ordered probit models**.

The survey in `../related_packages/README.md` and `../literature/README.md` found methods in the literature but **no released implementation** for panels:

- `spatialprobit`, `ProbitSpatial`, `pmlsp`, `spldv` and `McSpatial` are cross-sectional;
- `spmixW` handles panels, but only with continuous outcomes.

## Files

| File | Contents |
|---|---|
| [`01_model_specification.md`](01_model_specification.md) | The model family (SAR/SDM/SEM/SARAR/SDEM/SLX × pooled/RE/CRE × binary/ordered × static/dynamic); identification; default priors; probability and impact formulas; out-of-scope list. |
| [`02_estimation_algorithms.md`](02_estimation_algorithms.md) | Gibbs blocks and their order; the generic $Az - b$ engine; the precision-free latent sweep; the collapsed griddy $\rho$ draw; the RE blocks; cut-point samplers; log-determinants; impact computation; diagnostics; performance targets; **problems found in the `spatialprobit` source**. |
| [`03_package_architecture.md`](03_package_architecture.md) | User API and methods, data pipeline, file layout, R/C++ split, dependencies, repo and tooling conventions. |
| [`04_validation_and_testing.md`](04_validation_and_testing.md) | Unit tests; **simulation-based calibration as the release gate**; Monte Carlo design; cross-package checks; benchmarks; empirical examples. |
| [`05_roadmap.md`](05_roadmap.md) | Phases 0–7 with task checklists and exit criteria; risk register; reading list per phase. |
| [`06_rust_integration.md`](06_rust_integration.md) | Decided 2026-10-30 (savvy, Rust 1.84, Rust before the first CRAN submission); R0 under way on branch `rust` (C++ reference tagged `cpp-final`): moving the sampler and post-processing from C++ and R loops to Rust; measured baselines, decisions R-D1 to R-D10, phases R0–R7. |

## The plan in one paragraph

Every supported model is reduced, conditional on its spatial parameters, to a latent Gaussian kernel $\exp(-\frac12\|Az - b\|^2)$ with a sparse, period-block-diagonal $A$. A single C++ engine then serves all variants:

- **latent block:** an exact, warm-started, residual-update truncated-normal Gibbs sweep with $O(\mathrm{nnz})$ cost and no matrix solves;
- **$\rho$:** an exact collapsed griddy-Gibbs draw with $\beta$ integrated out;
- **random effects:** closed-form blocks with a half-Cauchy prior on $\sigma_\alpha$;
- **ordered models:** a joint $(\kappa, z)$ transformation move against the slow Albert–Chib cut-point sampler.

Impacts on probabilities are computed after sampling, accounting for the reduced-form variance. A pure-R reference sampler acts as the oracle. **No model variant is documented as supported until it passes simulation-based calibration.**

## Phases at a glance

| Phase | Scope | Size | Milestone |
|---|---|---|---|
| 0 | Setup: git, skeleton, CI, Rcpp/RcppEigen, CLAUDE.md | S | Clean `R CMD check` |
| 1 | Core: binary SAR probit, pooled and RE, C++ engine, R oracle, SBC | L | SBC passes |
| 2 | Impacts, CRE, time effects, SDM/SLX, `predict` | M | v0.1.0 (GitHub) |
| 3 | Ordered probit, cut-point samplers | M | **v0.5.0: first CRAN candidate** |
| 4 | SEM, SARAR, SDEM, $t$ errors, spatial RE | M | SBC per variant |
| 5 | Dynamics: $y_{t-1}$, $Wy_{t-1}$, Wooldridge initial conditions, (latent AR) | M–L | SBC passes |
| 6 | Time-varying $W$, missing $y$, parallel chains, mixing upgrades, PPC, approximate logLik | L | Performance targets |
| 7 | Vignettes, pkgdown, paper, CRAN 1.0 | M | Release |

## Decisions needed from you

Each decision has a recommendation; the plan assumes the recommended option unless you say otherwise.

| # | Decision | Recommendation | Alternatives |
|---|---|---|---|
| D1 | Function names | `psprobit()` (binary) + `psoprobit()` (ordered), sharing one engine | A single `spprobit_panel(family = )`; or `panel_sprobit()` / `panel_soprobit()` |
| D2 | License | GPL (≥ 3): compatible with porting from `spatialprobit`, `ProbitSpatial` and `spmixW` | GPL (≥ 2) (no porting from `spmixW`); MIT (clean-room only) |
| D3 | C++ backend | `RcppEigen` (zero-copy `dgCMatrix`, sparse LU and Cholesky; also used by `ProbitSpatial`). *Superseded 2026-10-31: the engine is Rust (`06`, R-D1).* | `RcppArmadillo` (weaker sparse solvers) |
| D4 | Diagnostics dependency | `posterior` in Imports (rank-normalised $\hat R$, ESS) | Suggests + own minimal $\hat R$/ESS |
| D5 | Repository location | Local disk + GitHub remote; keep Drive for papers only | Stay in Google Drive (risk of sync and lock issues during compilation) |
| D6 | Fixed effects | Not offered; CRE (Mundlak/Chamberlain) instead, with documentation on why | Offer with an incidental-parameters warning for large $T$ |
| D7 | Priority after the core | Ordered probit (Phase 3) before SEM/SARAR (Phase 4) and dynamics (Phase 5) | Dynamics first, if your application needs state dependence |
| D8 | Publication target | CRAN + R Journal (or JSS) paper | CRAN only |

## Additions to the earlier survey

- **`spldv`** (CRAN 0.1.3; Piras & Sarrias, JSS 2023, doi:10.18637/jss.v107.i08): GMM spatial probit and logit, cross-sectional. It was not in `related_packages/` and is worth adding as a comparator.
- `SBC` (simulation-based calibration) is GitHub-only, not on CRAN, so it is used only in `validation/`.
- The name `panelspatialprobit` is free on CRAN (checked 2026-10-27).
