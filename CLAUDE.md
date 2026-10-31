# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`panelspatialprobit` is an R package (R + C++ via RcppEigen) for Bayesian MCMC estimation of panel spatial probit and panel spatial ordered probit models.

The design and roadmap are in `plan/`:

| File | Contents |
|---|---|
| `plan/01_model_specification.md` | Models |
| `plan/02_estimation_algorithms.md` | Sampler math |
| `plan/03_package_architecture.md` | API and layout |
| `plan/04_validation_and_testing.md` | Testing |
| `plan/05_roadmap.md` | Phases and checklists |
| `plan/06_rust_integration.md` | Rust migration (decided 2026-10-30; R0 under way on branch `rust`, C++ reference tagged `cpp-final`) |

Read the relevant plan file before changing a sampler block. `literature/` and `related_packages/` hold reference papers and third-party package sources; both are excluded from the build.

## Commands

The project directory is on an NTFS mount where every file appears root-owned, so:

- `sed -i` fails with "preserving permissions"; edit files with python or the Edit tool;
- `git` complains about "dubious ownership" unless the directory is added to `safe.directory`.

Install into a scratch library rather than the user library.

Remove `src/*.o` and `src/*.so` before installing. Plain `roxygen2::roxygenise()` (and `pkgload::load_all()`) compile `src/` at `-O0` and leave the object files there. `R CMD INSTALL` then reuses them, and the installed sampler runs about 10x slower.

```sh
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'   # after changing C++ exports or roxygen docs (load_code = "source" breaks S3 registration)
rm -f src/*.o src/*.so && R CMD INSTALL --no-test-load -l <lib> .   # build + install
R_LIBS=<lib> Rscript -e 'library(panelspatialprobit); testthat::test_dir("tests/testthat", load_package = "none")'   # the package must be attached (helpers use its functions)
NOT_CRAN=true ...                                                  # also runs the C++ vs pure-R reference comparison (~2 min)
R_LIBS=<lib> Rscript -e 'library(panelspatialprobit); testthat::test_dir("tests/testthat", load_package = "none", filter = "impacts")'   # one file (test_file() skips helper-fits.R)
R CMD build . && R CMD check --as-cran --no-manual panelspatialprobit_*.tar.gz              # expected: 1 NOTE (new submission)
Rscript validation/sbc/run_sbc.R <reps> <cores> <lib> [variants]   # simulation-based calibration; writes validation/sbc/results/ (SBC_SEED_OFFSET=k for a replication)
Rscript validation/montecarlo/run_mc.R <reps> <cores> <lib>        # Monte Carlo study (bias, RMSE, coverage); writes validation/montecarlo/results/
Rscript validation/crosscheck/spatialprobit_katrina.R <lib>        # agreement with spatialprobit::sarprobit (T = 1)
Rscript validation/crosscheck/spatialprobit_ordered.R <lib> [cores] # agreement with spatialprobit::sarorderedprobit (T = 1, ~5 min)
Rscript validation/bench/bench_models.R <lib>                    # paper Table 5: SAR, ordered SAR and SEM timings and ESS/s (~6 min)
Rscript validation/wang_kockelman/run_wk.R <reps> <cores> <lib>  # Wang & Kockelman (2009) simulation (~1 h for 50 reps on 7 cores)
Rscript validation/vignette_data.R                                 # refresh vignettes/validation-*.csv after new SBC / Monte Carlo / benchmark results
Rscript validation/bench/bench_postproc.R <lib>                  # impacts, predict, PPC, loglik and log-det timings (~20 s)
```

**Rust (branch `rust`, plan/06).** The package builds the Rust workspace in `src/rust` through `configure` and `src/Makevars.in`. Until phase R5 it still holds the C++ engine too, and `src/registration.c` registers both (it compiles savvy's generated `init.c` itself, so new `#[savvy]` functions are picked up automatically). `R/rust.R` converts R objects for Rust (`.psp_csc()`, `.psp_rs_weights()`); the `psp_rs_*` functions are internal. The frozen C++ reference (tag `cpp-final`) is installed in `~/.cache/panelspatialprobit/lib_cpp`.

```sh
export PSP_CARGO_TARGET_DIR=$HOME/.cache/panelspatialprobit/target   # always: keeps cargo's ~1 GB output off the Drive mount and builds incrementally
(cd src/rust && cargo test -p pspcore && cargo clippy --workspace --all-targets -- -D warnings && cargo fmt --all --check)   # cargo is in ~/.cargo/bin
~/.cache/panelspatialprobit/tools/bin/savvy-cli update . src/rust/pspr   # after changing a #[savvy] function; then roxygenise
sh tools/vendor.sh   # writes src/rust/vendor.tar.xz (not committed) before R CMD build: R CMD check warns when an install downloads crates
Rscript tools/rust_authors.R   # rewrite inst/AUTHORS after any change to src/rust/Cargo.lock
(cd src/rust && MIRIFLAGS=-Zmiri-deterministic-floats cargo +nightly miri test -p pspcore -- --skip par::)   # undefined-behaviour check (~1 min)
Rscript validation/rust_port/check_interrupt.R <lib>   # Ctrl-C during a multi-threaded Rust job (unix)
Rscript validation/rust_port/compare_sampler.R <lib_cpp> <lib_rust> [iter warmup chains] [case ...]   # Rust vs frozen C++ sampler: same draws, or posterior means (~2 min; 3000 1000 2 for the sparse-RE cases)
Rscript validation/rust_port/check_fit_threads.R <lib>   # Ctrl-C, progress, mclapply and cores with the Rust sampler (unix)
```

## Architecture

- **Data pipeline.** `R/model_frame.R` (`.psp_design`):
  - matches data units to the row names of `W` (strict; errors on mismatch);
  - requires a balanced panel;
  - **sorts observations time-slow / unit-fast**, the stacking convention everywhere (C++, impacts, predictions);
  - keeps `row_order` so `predict()` returns results in the user's row order (for dynamic models, the order of the non-initial rows);
  - aligns the error weights `M` (sem/sdem/sarar) with `W`;
  - accepts time-varying weights: `W`/`M` may be a list with one matrix per period (`.psp_align_weights`); a list of identical matrices collapses to one matrix. Code reaches the weights of period `t` through `.psp_Wt(W, t)`, never `W` directly.
  
  It builds the extra regressor columns, tagged in `roles` (`intercept`, `x`, `durbin`, `cre`, `time`, `dynamic`, `initial`):
  - Durbin (`W_x`);
  - Mundlak (`mean_x`) or Chamberlain (`cre_x_<period>`);
  - lagged outcomes (`lag_y`, `W_lag_y`) and Wooldridge initial conditions (`init_y`) for `dynamic = psp_dynamic()`, after which the first period is dropped;
  - time dummies (`time_<level>`).
  
  Only `x` columns get impacts; `durbin_of` maps `W_x` columns back to their base variables.
- **Engine.** `R/psprobit.R` (`.psp_fit`):
  - computes the log-determinant grids for `W` and `M` (`.psp_lndet` in `R/weights.R`);
  - runs one C++ call per chain with L'Ecuyer-CMRG streams, so results are independent of `cores` and the user's RNG state is restored;
  - assembles a `posterior` draws_array.
- **Sampler.** On branch `rust` the sampler is the Rust port (`src/rust/pspcore/src/sampler/`, called through `psp_rs_sample()`; plan/06), block for block the C++ one below; `options(panelspatialprobit.engine = "cpp")` still selects `src/sampler.cpp` until R5. Chains run on Rust threads (`cores`) on every OS.
- **C++ sampler.** `src/sampler.cpp` (`psp_sampler_cpp`) runs one chain. Every model is the kernel `|A z - b|^2` with block-diagonal `A = R S`, `S_t = I - rho W_t`, `R_t = I - lambda M_t` (plan/02 s2.2). Weights arrive as lists (`Ops`, one matrix or one per period; `at(ops, t)`), and the log-determinant grids are averages over periods, so `T * lndet` is the sum. Block order matters: z sweep → rho (beta integrated out, given lambda) → beta → lambda → alpha → sigma_alpha → cut-points → scale move. Rho is collapsed over beta, so beta must be redrawn immediately after rho (partially collapsed Gibbs rule).
  - **Latent sweep:** maintains the residual `r = A z - b` and updates it in place using the columns of the sparse `A`, rebuilt once per iteration (no matrix solves).
  - **Random effects:** independent normal draws without an error term; a sparse Cholesky of `sum_t R_t'R_t + B'B / sigma^2` (pattern analysed once) with one or with spatial random effects (`effects = "spatial_random"`, `B = I - phi W`, phi drawn on the `W` grid after alpha).
  - **Scale move:** a generalised Gibbs move on (z, beta, alpha, kappa).
  - **Truncated normal sampler:** `src/truncnorm.h`.
  - **Student-t errors** (`errors = psp_t(df)`, fixed df): precision weights `wt = 1/v` enter every quadratic form (sweep, rho/beta/lambda, alpha, cut move, scale move); `v` is drawn after alpha. With t errors `B` is recomputed every iteration.
  - **Latent AR(1)** (`psp_dynamic(latent = TRUE)`): `gamma` is the coefficient of the design column `lag_z`, which the sampler fills with `z_{t-1}` (`refresh_lag`, before every residual computation). The full `A` gains `-gamma R_t` below the diagonal: the sweep and the cut move include the period-`t+1` coupling, and the scale move leaves `gamma` unscaled.
- **Impacts and prediction.** `R/impacts.R`: `.psp_reduced_form` gives `m = S^{-1} X beta`, the reduced-form sd (latent covariance `sigma_alpha^2 S^{-1}S^{-T} + K K'`, `K = S^{-1} R^{-1}`) and the diagonals needed for impacts. `type = "unit"` impacts condition on the random-effect draws that the sampler keeps (`alpha_draws`, `alpha_index`). It is exact (dense inverse) when `N <= impacts_exact_max`, otherwise it uses a stochastic Rademacher-probe estimator. `predict()` in `R/methods.R` reuses it.
- **Reduced form variants.** `.psp_rf_block` takes `mix` (t errors: `s` becomes an `n x Q` matrix over simulated mixing variances, averaged by `.psp_mix_mean`); `.psp_rf_latent_ar` handles latent AR (mean recursion, lagged-error and random-effect variance sums; probes then have `P (T + 1)` columns).
- **Response paths.** `R/response.R` (`psp_response_path`) simulates outcome histories forward (baseline vs shock, common random numbers), rebuilding the `lag_*`/`W_lag_*` columns, and computes probabilities exactly given each history via `.psp_rf_block`.
- **Approximate likelihood.** `R/loglik.R` builds the sparse precision of `(z, alpha)` and runs the GHK simulator in `src/ghk.cpp`. Unrestricted coordinates (alpha, missing outcomes) are ordered first so GHK integrates them out exactly; independent periods are simulated as separate blocks.
- **Test oracle.** `R/reference_sampler.R` is a deliberately different pure-R sampler (dense precision conditionals, Metropolis for rho and lambda), used only as an oracle in `test-reference.R`.
- **Cut-point kernels.** `psp_cut_kernel_cpp` exposes the cut-point updates on a toy model for the stationarity test in `test-cutmove.R`.

## Conventions

- **Mixing changes.** A new sampler move must show an ESS gain on the cases in plan/02 §2.6 ("Mixing experiments") and pass the reference-sampler tests. That section lists three moves that were tried and reverted; read it before trying another.
- **Model support gate.** New model variants must pass SBC (`validation/sbc/run_sbc.R`, add the variant) before they are documented as supported.
- **Parameter names.** Draw columns are named: regressor names, `rho`, `lambda`, `sigma_alpha`, `cut[a|b]`. Impacts and prediction look parameters up by these names.

## Still needed (as of 2026-10-30)

Phases 1–6 and vignettes (7.1) are done; `R CMD check --as-cran` gives 1 NOTE (new submission). The repository is on GitHub (`gueyenono/panelspatialprobit`, remote `origin`, branch `master`). The items below are open; details are in `plan/05_roadmap.md`.

**Needs the maintainer (do not do these without being asked):**

1. **Git setup** (Phase 0). Done 2026-10-30: `safe.directory` is set, and the repository is on GitHub with remote `origin`. Do not change the global git config on the user's behalf. Commit or push only when asked.
2. **Author details** (Phase 0, 7.3). Replace the placeholder `person("Package", "Author", ...)` in `DESCRIPTION` (`Authors@R`; add an ORCID if wanted) and in `inst/CITATION` (both `author` and `textVersion`). Once the paper exists, add it to `inst/CITATION`.
3. **Repository URLs** (7.3). Done 2026-10-30: `url:` in `_pkgdown.yml` and `URL:`/`BugReports:` in `DESCRIPTION` point to `gueyenono/panelspatialprobit` and its GitHub Pages site. The three C++-era workflows in `.github/workflows/` pass (commit `968ea92`); the pkgdown workflow deletes `CLAUDE.md` before building, because pkgdown renders every top-level `.md` file. GitHub Pages serves the `gh-pages` branch at that URL (enabled 2026-10-30).
4. **Empirical example data** (7.2). Choose the dataset(s) from `plan/04` §4.6:
   - binary: the Elhorst et al. (2017) inflation-targeting panel, CC-BY 4.0;
   - ordered: US Drought Monitor, after checking its redistribution terms.

   Build them in `data/` from scripts in `data-raw/`, with the licence noted, and document them. Then use them in the vignettes, which so far rely on simulated data only.
5. **Platform checks** (7.4). Run the `rhub` and win-builder checks. They email the maintainer, so the user must start them, and only after item 2. Deferred with item 8; they run on the Rust build.
6. **Paper** (7.5). A draft is in `paper/paper.md`: the abstract, outline and full text of Sections 1–6, 8 and 9 (all but the Section 7 application, which needs item 4). Choose the venue, then convert it to that venue's format. Section 7 needs writing, the abstract and conclusion each have a placeholder sentence for it, and the `[TODO]` marks need checking. The paper is for the R Journal or JSS, covering:
   - the models;
   - the `Az - b` engine and the cut-point move;
   - SBC and Monte Carlo evidence, as summarised in `vignette("validation")`;
   - an empirical application, which needs item 4.
7. **Reference papers** (Phase 0). Obtain the paywalled papers in the "Reading to obtain" table of `plan/05_roadmap.md`. Then check the model specifications against the published versions, especially Baltagi, Egger & Kesina (2018) and Wang & Kockelman (2009). The Wang & Kockelman replication (`validation/wang_kockelman/`) follows the 2008 TRB working-paper version in `literature/`; check its design and Table 1 against the published paper.
8. **CRAN submission.** Deferred by the maintainer (2026-10-30), with no date. When it happens, it comes after items 1–5 and the Rust migration (`plan/06_rust_integration.md`: phases R0–R5, then its "Before the CRAN submission" list). The first CRAN release is the Rust version (decision R-D8).

**Optional work that Claude can do on request:**

- Latent-AR mixing: a joint update of each unit's latent path is ruled out (plan/02 §2.6: four sweeps per iteration already approximate its best case without an ESS gain). A real gain needs a move that changes the parameters and the latent values together, or the marginal likelihood; see the diagnosis there before trying one.
- Convert `README.md` to `README.Rmd` (`.Rbuildignore` already lists it) once real example data exist.
- After any new SBC, Monte Carlo, Wang & Kockelman, `bench_models.R` or `bench_cutpoints.R` run, rerun `validation/vignette_data.R` so that `vignette("validation")` stays current.
