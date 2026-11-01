# 6. Rust integration

*Written 2026-10-30. Status: decisions R-D1, R-D7 and R-D8 taken by the maintainer on 2026-10-30 (savvy; Rust 1.84; Rust before the first CRAN submission); the CRAN submission itself is deferred, with no date. R0-R4 done 2026-10-31 (the last C++ commit is tagged `cpp-final`; work continues on branch `rust`); R5 under way (C++ removed, docs updated; GitHub Actions and the merge into `master` pending). Measurements were taken on the development machine (8 cores, R 4.6.0, rustc 1.97.1) on 2026-10-30. The spike crates are not in the repository.*

**Goal.** Move the heavy numerical work from C++ (RcppEigen) and from R-level loops into Rust, use Rust's safe threading where it pays off, and then delete the C++ code. The statistics stay fixed:

- the same models (`01`);
- the same sampler blocks, in the same order (`02` §2.1);
- the same validation gates (`04`, and the SBC gate in CLAUDE.md).

## 6.1 Summary

- **Two crates** in `src/rust/`:
  - `pspcore`: pure Rust with no R dependency. It holds the RNG, distributions, sparse operators, sampler, reduced form, impacts, prediction, simulation kernels and GHK, and is tested with `cargo test`.
  - `pspr`: the R binding, built with [savvy](https://github.com/yutannihilation/savvy). It converts R objects, runs threads, checks for interrupts and prints progress.
- **Moves to Rust:**
  - the MCMC sampler (every variant), the truncated normal sampler and the cut-point test harness;
  - the sparse log-determinant grid;
  - the reduced form, impacts and `predict()`;
  - the simulation loops of `psp_ppc()` and `psp_response_path()`;
  - the whole GHK likelihood.
- **Stays in R:**
  - formula handling and the design;
  - weights preparation;
  - `eigen()`;
  - seeds and streams, probes and t-mixing variances;
  - `sim_psprobit()`;
  - the reference sampler;
  - summaries and methods.

  The reasons are in §6.5.
- **Random numbers.** The RNG is an exact port of R's L'Ecuyer-CMRG generator and of R's normal (inversion), exponential and gamma generators. Seeds keep their meaning, chains stay independent of `cores`, and during the port the Rust sampler can be compared draw by draw with the C++ one.
- **Threads.** Scoped `std` threads only, with no global pool and no rayon. They run:
  - chains in parallel on every OS (Windows runs them sequentially today);
  - posterior draws in parallel during post-processing;
  - later (R6), periods in parallel inside the latent sweep, which gives the same Markov kernel.
- **Linear algebra.** [faer](https://crates.io/crates/faer) 0.24: sparse LU and Cholesky with reusable symbolic analysis, AMD ordering and dense factorisations, in pure Rust.
- **Migration.** Work on a feature branch. The last C++ version is frozen and installed in its own library as the reference. C++ is deleted only after the full SBC, Monte Carlo and cross-check gate passes on Rust (R4).
- **Phases R0 to R7** are in §6.11; R0 to R5 are required, R6 and R7 come after.

## 6.2 What Rust can and cannot gain here

### Baseline (current C++ and R code)

| Workload | Now | Where the time goes |
|---|---|---|
| SAR-RE binary, N = 1,000, T = 10, k = 5 | 1.31 ms/iter (`bench_models.R`); 1.71 ms in a 2,000-iteration run | One latent sweep costs 0.65 ms (measured as the cost of an extra `z_sweeps`), so 38-50% of an iteration. The rest is residual passes, `build_A` (a new sparse matrix every iteration) and the rho grid. |
| SEM-RE binary, same size | 2.69-2.86 ms/iter | Sweep 0.75 ms. The rest: the sparse product `R S` every iteration, the sparse Cholesky of the RE precision, and `M X` with a k x k inverse after every lambda draw. |
| SAR binary, N = 10,000, T = 20 | 32.5 ms/iter | Same split at larger n. |
| `psp_impacts()`, N = 10,000, 50 draws | 19 s | Two sparse LUs per draw (`solve(t(S), .)` and `solve(S, .)`), each with 50-222 right-hand sides, in an R loop. |
| `psp_impacts()`, N = 1,000, 100 draws | 2.0 s | Same, at the stochastic-diagonal size (N > `impacts_exact_max`). |
| `predict()`, N = 1,000, 100 draws | 1.8 s | Reduced form per draw. |
| `psp_ppc()`, N = 1,000, T = 10, 100 draws | 3.7 s | T sparse LU factorisations per draw (one `Matrix::solve` per period). |
| `psp_loglik(type = "conditional")`, n = 10,000, 10 draws x 100 reps | 2.7 s | Precision and ordering in R/Matrix, then the C++ GHK loop. |
| Sparse log-det grid, N = 1,000 | 0.4 s | About 124 `Matrix::determinant` calls. |
| 4 chains, `cores = 4`, on Windows | sequential | `mcmapply` is unix-only. |

### Spike results

The spike was faer 0.24.4 with savvy 0.11.0 in a scratch crate.

- **Solve workload.** The impacts solve was timed per draw: assemble `S = I - rho W`, factorise it numerically with the symbolic analysis reused, and solve 222 right-hand sides (rook grid, row-standardised W):

  | N | faer | `Matrix::solve` |
  |---|---|---|
  | 10,000 | 0.079 s | 0.242 s |
  | 1,024 | 0.0106 s | 0.0113 s |

  The symbolic analysis took 6 ms once.
- **Build.** An offline release build (`cargo build --release --offline -j 2`) took 82 s and peaked at about 1 GB of memory.
- **Vendored size.** The full `cargo vendor` tarball is 5.5 MB with xz. It drops to 3.6 MB once `windows-sys` is stubbed out: it is reachable only through `cfg(apple)`, then `cfg(windows)`, dependencies of faer's `gemm-common`, so it is never compiled on any target.
- **Shared object.** A shared object that uses only the sparse LU and Cholesky is 9.9 MB, or 5.0 MB stripped.

### What follows from the measurements

- **The sampler port is not where the speed comes from.** The sampler is already compiled. Rust running the same algorithm runs at about the same speed. The gains come from changes made during the port:
  - value-only updates of `A` and of the RE precision on precomputed sparsity patterns, so nothing is allocated per iteration;
  - precomputed truncation bounds;
  - fewer residual passes.

  After that, R6 runs the periods of the latent sweep in parallel; the sweep is 38-50% of an iteration and is sequential today.
- **Post-processing gains most.** It gets:
  - one LU per draw, serving both `S` and `S'`;
  - the symbolic analysis done once per fit;
  - faer's faster factorisation at large N (3x at N = 10,000, about the same at N = 1,000);
  - posterior draws spread over threads.

  At N = 1,000 the single-thread gain is small (about 1.5x); most of it comes from threads.
- **Windows users get parallel chains.**

### Costs to accept

- Installing from source needs Rust ≥ 1.84. CRAN's Windows and macOS binaries do not.
- Installation takes about 1.5 min longer.
- The source tarball is about 3.6 MB larger.
- The installed library is several MB.
- `Depends: R (>= 4.5)`, because savvy requires it.
- Draws differ from the C++ version for the same seed.
- CRAN reviews Rust packages more closely (§6.9).

### Targets (to confirm in R0, enforced in R2/R3/R6)

| Workload | Port target (R2/R3) | After R6 |
|---|---|---|
| Sampler, any variant, 1 thread | within ±10% of C++ per iteration; never slower by more than 10% | SAR ≥ 1.3x, SEM/SARAR ≥ 1.5x faster than C++ |
| Sampler, N = 10,000, T = 20, 4 threads per chain | n/a | ≥ 2.5x faster than C++ |
| Impacts, N = 10,000, 50 draws | ≤ 6 s (1 thread), ≤ 2 s (4 threads) | |
| Impacts / predict, N = 1,000, 100 draws | ≤ 1.3 s (1 thread), ≤ 0.5 s (4 threads) | |
| `psp_ppc()`, N = 1,000, T = 10, 100 draws | ≤ 0.4 s (1 thread) | |
| GHK conditional, n = 10,000, 10 draws x 100 reps | ≤ 2.5 s (1 thread), ≥ 3x faster with 4 threads | |

## 6.3 Decisions

Each decision has a recommendation; this plan assumes it unless the maintainer says otherwise. **R-D1, R-D7 and R-D8 were decided by the maintainer on 2026-10-30.**

| # | Decision | Recommendation | Alternatives |
|---|---|---|---|
| R-D1 | R binding layer | **Decided: savvy 0.11** (MIT, 2026-09; MSRV 1.81; needs R ≥ 4.5). It generates the C and R wrappers (`savvy-cli update`). Its R API calls run inside `R_UnwindProtect`. A panic is caught and becomes an R error (needs `panic = "unwind"`), and errors are raised by the C wrapper after Rust returns, so no R `longjmp` crosses Rust frames (checked in the vendored source). Its README lists polars, prqlr, string2path and SedonaDB as users. | **extendr** (rextendr 0.5.0 is installed here; `extendr-api` 0.9 has MSRV 1.71 and supports older R; heavier macro layer, faster-moving API). **Hand-written C ABI** (no binding dependency, any R version; about 500 lines of C glue and an `unsafe` Rust boundary). Pick the C ABI if R 4.1-4.4 support matters. |
| R-D2 | Linear algebra | **faer 0.24.4, pinned exactly** (`=0.24.4`), `default-features = false`, features `std`, `linalg`, `sparse`, `sparse-linalg` (no `rayon`, `rand` or `npy`). It covers sparse LU and Cholesky with a symbolic/numeric split, AMD and custom orderings (`SymmetricOrdering::Custom`), transpose solves, and dense LU and Cholesky. SIMD kernels are chosen at run time, so `-C target-cpu=native` is never needed. | `nalgebra` plus hand-written sparse code (no sparse LU: we would write and maintain one); a Rust port of CSparse (LGPL-2.1+; smaller, but ours to maintain). |
| R-D3 | Random numbers | **Exact port of R's generators** (§6.6). | `rand`/`rand_distr`: faster normals (ziggurat), but seeds lose their meaning and the port can only be checked distributionally. Can be revisited after R5 (R7) if profiling shows the RNG dominates. |
| R-D4 | Threads | **`std::thread::scope`**, with a fixed worker count and work queues. Threads never outlive a call. rayon's global pool keeps threads alive, and a later `fork()` (e.g. `parallel::mclapply()` calling `psprobit()`) then deadlocks in the child. | rayon with a local `ThreadPool` built and dropped per call (adds a dependency for no gain here). |
| R-D5 | Migration strategy | **Feature branch `rust` plus a frozen C++ reference build.** Tag the last C++ commit (e.g. `cpp-final`), install it into its own library (`lib_cpp`), and compare across R processes (§6.10). This means two builds, not two backends: each R caller switches to Rust as its port lands, and there is no run-time switch between engines. Until R5, though, the branch package still holds the C++ code next to the Rust code, and Rcpp and savvy would each define `R_init_panelspatialprobit`. *Found in R0 (2026-10-30):* `src/registration.c` registers both sets of routines in one table. `Rcpp::compileAttributes()` then leaves registration to it; it compiles savvy's generated `init.c` itself, with that file's `R_init` renamed, and copies its routine table at load time (R1). Both it and the `OBJECTS` line of the Makevars are removed in R5. | A dual backend inside one package, selectable at run time (useful only if the port must ship piecemeal). |
| R-D6 | Scope | §6.5. | Port everything numeric (rejected: loses the reference sampler as an oracle and changes `sim_psprobit()` data). |
| R-D7 | Minimum Rust version | **Decided: 1.84**, set by faer (released 2025-01-09). CRAN asks for testing with a cargo at least two years old, and preferably four; 1.84 turns two in January 2027. | An older faer line (older API, no fixes); dropping faer (R-D2 alternatives). |
| R-D8 | Order relative to the first CRAN submission | **Decided: the first CRAN submission is the Rust version**, after R5. The submission itself is deferred, with no date (maintainer, 2026-10-30). The C++ version is never released on CRAN. Consequences: the first-submission review also covers the Rust build. The build follows CRAN's rules from R0 on (cargo home inside the build tree, `-j 2`, rustc version in the log, Rust version check), but the CRAN-only work waits for §6.11 "Before the CRAN submission": vendoring, the CRAN flavour check, win-builder/mac-builder/rhub, and `cran-comments.md`. | Submit the C++ version first and ship Rust in the next release (the original recommendation: earlier submission, and the build system is not reviewed with a first submission). |
| R-D9 | Reproducibility contract | **Same seed and same machine and build give identical draws, whatever `cores`/threads**, as tested now (`test-fit.R`). Not promised: identity with the C++ version, or across CPUs. faer picks SIMD kernels at run time, and the eigen log-det path already depends on the BLAS. Documented in NEWS and `?psp_control`. | Bitwise identity across CPUs (would rule out faer in the sampler and the log-det; not worth it). |
| R-D10 | Data across the boundary | **R unpacks `dgCMatrix` into plain vectors** (`p`, `i`, `x`, `n`), and Rust copies its inputs into owned structures once per call. That costs O(nnz + nk), which is negligible against the work. Outputs are built in Rust and copied into R objects on the main thread. **No persistent external pointers**, since they do not survive `saveRDS()`. | Zero-copy borrowing of R memory (lifetimes tied to R objects; R memory touched from worker threads). |
| R-D11 | Post-processing parallelism | **Decided by the maintainer (2026-10-31): the {future} package.** `psp_impacts()`, `predict()`, `psp_ppc()`, `psp_response_path()`, `psp_loglik()` (and `logLik()`, `psp_dic()`) follow the user's `future::plan()`; the default `sequential` plan runs in the calling process. R splits the draws into fixed chunks (whatever the plan, so results do not depend on it), and each future calls a single-threaded Rust function on its chunk; chunks that draw random numbers get their own L'Ecuyer-CMRG stream from R (`.psp_streams()`), so `future.seed` is not needed. `future` goes in `Imports`. The sampler's chains stay on Rust threads (`cores`, R-D4). | `fit$control$cores` with an `options(panelspatialprobit.threads = )` override, on Rust threads (the earlier recommendation; lower overhead, no new dependency); a `cores =` argument per function. |

## 6.4 Target layout

```
DESCRIPTION                 # SystemRequirements: Cargo (Rust's package manager), rustc >= 1.84
configure, configure.win    # one line each: run tools/config.sh
cleanup, cleanup.win        # remove the generated Makevars and src/.cargo
tools/
  config.sh                 # find cargo, log versions, check MSRV, write src/Makevars(.win)
  vendor.sh                 # cargo vendor -> stub unreachable crates -> src/rust/vendor.tar.xz
  rust_authors.R            # inst/AUTHORS from `cargo metadata` (crate, version, authors, licence)
inst/AUTHORS, inst/COPYRIGHTS
R/
  000-wrappers.R            # generated by savvy-cli (replaces RcppExports.R)
  rust.R                    # .psp_csc(), .psp_rs_weights(), .psp_rs_model(), .psp_threads()
  reference_postproc.R      # the current R reduced form / impacts / predict, kept as a test oracle
src/
  Makevars.in, Makevars.win.in
  init.c                    # generated by savvy-cli
  registration.c            # R0-R4 only: registers the Rcpp and savvy routines together (R-D5)
  rust/
    Cargo.toml              # workspace: pspcore, pspr; resolver 3; [profile.release]
    Cargo.lock              # committed
    deny.toml               # cargo deny licence allow-list (CI only; .Rbuildignore'd)
    vendor-config.toml      # [source] replacement used only for vendored (CRAN) builds
    vendor.tar.xz           # NOT in git; produced by tools/vendor.sh before R CMD build
    api.h                   # generated by savvy-cli
    pspcore/
      Cargo.toml            # depends on faer only
      src/
        lib.rs              # #![forbid(unsafe_code)]
        error.rs            # PspError (validation, numerical failure, cancelled)
        rng/                # mrg32k3a.rs (streams, substreams), normal.rs, exp.rs, gamma.rs, rt.rs
        dist/               # truncnorm.rs, pnorm.rs (R's Cody algorithm incl. log), dnorm.rs
        sparse/             # csc.rs, weights.rs (unique matrices + period map), pattern.rs (value-only updates)
        dense.rs            # small k x k Cholesky/inverse helpers (or faer::Mat)
        sampler/            # mod.rs (run_chain), state.rs, input.rs, latent.rs, griddy.rs, beta.rs,
                            # lambda.rs, re.rs, tmix.rs, cut.rs, scale.rs, store.rs
        lndet.rs
        rf/                 # block.rs (.psp_rf_block), latent_ar.rs, mix.rs, solver.rs (LU cache)
        impacts.rs, predict.rs, ppc.rs, response.rs
        ghk/                # precision.rs, ordering.rs, simulate.rs
        par.rs              # scoped worker pool, fixed chunking, cancel flag, progress counters
      tests/                # integration tests; fixtures generated from R (tests/data/*.txt)
    pspr/
      Cargo.toml            # crate-type = ["staticlib", "lib"]; depends on savvy, pspcore
      api.h                 # second copy written by savvy-cli; .gitignore'd
      src/
        lib.rs              # #[savvy] functions only
        convert.rs          # R list -> pspcore input structs, with field-named errors
        run.rs              # threads + interrupt polling + progress printing (main thread)
        interrupt.rs        # R_CheckUserInterrupt inside savvy::unwind_protect
rust-dev/                   # optional, .Rbuildignore'd: criterion benches, proptest (dev-deps kept
                            # out of the vendored workspace, since cargo vendor includes dev-deps)
```

**Release profile** (`src/rust/Cargo.toml`):

```toml
[profile.release]
opt-level = 3
codegen-units = 1     # R0: measure against 16 (install time vs speed)
lto = "thin"          # R0: measure off / thin / fat
panic = "unwind"      # required: savvy turns panics into R errors; "abort" would kill the R session
debug = false
```

**Lints:**

- `pspcore` forbids `unsafe`.
- `pspr` denies:
  - `clippy::print_stdout` and `clippy::print_stderr`, because CRAN forbids writing to the console other than through R's own print functions;
  - `clippy::unwrap_used` and `clippy::expect_used`, because errors must reach R as messages.

**Core types (sketch):**

```rust
pub struct Csc { pub n: usize, pub colptr: Vec<usize>, pub rowind: Vec<usize>, pub val: Vec<f64> }
pub struct Weights { pub mats: Vec<Csc>, pub period: Vec<usize> } // period t uses mats[period[t]]

pub trait Rng {                       // implemented by RMrg32k3a (R-compatible)
    fn unif(&mut self) -> f64;        // unif_rand()
    fn norm(&mut self) -> f64;        // norm_rand(), inversion
    fn exp(&mut self) -> f64;         // exp_rand()
    fn gamma(&mut self, shape: f64, scale: f64) -> f64; // rgamma()
}

pub fn run_chain(input: &SamplerInput, rng: &mut RMrg32k3a, ctl: &ChainControl)
    -> Result<ChainOutput, PspError>; // ChainControl: cancel flag, progress counter
```

## 6.5 What moves to Rust, what stays in R

### Moves to Rust

| # | Component | Today | Rust | Notes |
|---|---|---|---|---|
| 1 | MCMC sampler, all variants | `src/sampler.cpp` | `pspcore::sampler` | §6.7 |
| 2 | Truncated normal | `src/truncnorm.h` | `pspcore::dist::truncnorm` | Same algorithms (normal, Robert 1995 exponential and uniform rejection). |
| 3 | Cut-point kernel harness | `psp_cut_kernel_cpp` | `psp_rs_cut_kernel` | For `test-cutmove.R`. |
| 4 | Sparse log-det grid | `Matrix::determinant` at about 124 points (`.psp_lndet`) | `pspcore::lndet` | One symbolic LU, then a numeric LU per point, points spread over threads; log\|det\| from the diagonal of U. The spline stays in R. |
| 5 | Eigen log-det grid evaluation | `vapply` over about 4,000 points | `pspcore::lndet` | Eigenvalues still come from R's `eigen()`. |
| 6 | Reduced form | `.psp_reduced_form`, `.psp_rf_block`, `.psp_rf_latent_ar` | `pspcore::rf` | §6.8 |
| 7 | Impacts | draw loop in `psp_impacts()`, `.psp_density_weights` | `pspcore::impacts` | Summary table stays in R. |
| 8 | Prediction | draw loop in `predict.psprobit()` | `pspcore::predict` | `psp_classification()` benefits too. |
| 9 | PPC replicates | `vapply` in `psp_ppc()` | `pspcore::ppc` | Returns the `n x draws` outcome matrix; statistics stay in R. |
| 10 | Response paths | simulation and probability loop in `psp_response_path()` | `pspcore::response` | Summary stays in R. |
| 11 | GHK likelihood | `.psp_ghk_loglik` and `psp_ghk_cpp` | `pspcore::ghk` | Precision, ordering, Cholesky and simulation; the final log-mean-exp can stay in R. |

### Stays in R

| Component | Why |
|---|---|
| Formula, model frame, design (`.psp_design`), weights alignment, `psp_weights()`, `rook_weights()` | R's job: `model.matrix`, factors, strict ID matching, `Matrix` coercions. Runs once per fit. |
| `eigen()` for the eigenvalue log-det | LAPACK; runs once; Rust adds nothing. |
| Spline interpolation of the sparse log-det grid (`splinefun`) | Cheap, and already validated. |
| `.psp_streams()` (chain seeds), `.psp_probes()`, `.psp_t_mix()` | Cheap; keeping them in R keeps the current seed semantics. |
| `sim_psprobit()` | Must keep producing the same data for a given seed. The tests, the SBC and Monte Carlo results in `validation/` and the vignettes depend on it. |
| **The reference sampler** (`R/reference_sampler.R`) | Its value as an oracle is that it shares no code with the engine. **Never port it.** |
| The current R reduced form, impacts and predict code | Moved to `R/reference_postproc.R` and kept as a deterministic test oracle for items 6-8 (§6.10). |
| posterior summaries and diagnostics; `print`, `summary`, `plot`, `tidy`, `glance`; PPC statistics; impacts and response summaries | Cheap, and they belong with R's object system. |

### The R-Rust interface

These are internal functions; savvy generates the wrappers. `model` is built once per call by `.psp_rs_model(fit_or_design)`. It is a list holding:

- `y` codes, and `X` as a column-major double matrix;
- `W` and `M` as `list(mats = <CSC lists>, period = <0-based index per period>)`, with identical matrices collapsed as `.psp_lndet_avg` already does;
- `N`, `T`, `J`;
- the model flags;
- the column roles: indices of the `x` columns, the `durbin_of` map, and the `lag_z`, `lag_*` and `W_lag_*` columns.

| R wrapper | Replaces | Returns |
|---|---|---|
| `psp_rs_sample(model, spec, seeds, threads, progress)` | `psp_sampler_cpp` per chain plus `mcmapply` | A list of chains with exactly today's fields (`draws`, `z_mean`, `alpha_mean`, `alpha_sd`, `cut_accept`, `cut_step`, `alpha_draws`, `v_mean`), so the assembly in `.psp_fit` is unchanged. `seeds` is a 6 x chains integer matrix from `.psp_streams()`. |
| `psp_rs_rng(kind, n, par, seed)` | new, for tests | Draws from `unif`, `norm`, `exp`, `gamma` or `t`. |
| `psp_rs_rtnorm(n, mean, sd, lo, hi, seed)` | `psp_rtnorm_cpp` | |
| `psp_rs_cut_kernel(...)` | `psp_cut_kernel_cpp` | Same arguments. |
| `psp_rs_lndet_sparse(W, points, threads)`, `psp_rs_lndet_eigen(re, im, grid)` | the two branches of `.psp_lndet` | |
| `psp_rs_reduced_form(model, par, alpha, probes, mix, need_impacts)` | `.psp_reduced_form` (one draw) | For tests and for any remaining R caller. |
| `psp_rs_impacts(model, par, alpha, probes, mix, opts, threads)` | the loop in `psp_impacts()` | Today's 5-d array. |
| `psp_rs_predict(model, par, probes, mix, threads)` | the loop in `predict()` | Mean link and probabilities. |
| `psp_rs_ppc(model, par, alpha, seed, threads)` | the `vapply` in `psp_ppc()` | Integer `n x draws`. |
| `psp_rs_response(model, par, alpha, shock, seed, threads)` | the loop in `psp_response_path()` | Today's 4-d array. |
| `psp_rs_ghk(model, par, alpha, reps, seed, threads)` | `.psp_ghk_loglik` and `psp_ghk_cpp` | Per draw: `ll`, `mcse`, or the log weights. |

`spec` is a named list: priors, grids (`rho_grid`, `lndet_grid` and the lambda equivalents), `nu`, `lag_col`, MCMC settings and starting values. In effect it is the argument list of `psp_sampler_cpp` today, minus the data. Rust parses it field by field and reports a missing or ill-typed field by name.

## 6.6 Random numbers

**Generator.** `RMrg32k3a` is an exact port of the `LECUYER_CMRG` branch of R's `RNG.c`, with state of six 32-bit integers. Its draws:

| Method | Port of | Notes |
|---|---|---|
| `unif()` | `unif_rand()` | |
| `norm()` | `norm_rand()` with `INVERSION` | `u = (int)(2^27 unif) + unif` then `qnorm5(u / 2^27)`; port AS241 as in `qnorm.c`. The extreme-tail branch added in R 4.3 is unreachable here, because u ≥ 2^-59. |
| `exp()` | `exp_rand()` | Ahrens and Dieter 1972, with R's `q[]` table. |
| `gamma(a, scale)` | `rgamma()` in `rgamma.c` | GD for a ≥ 1, GS for a < 1. R's static cache is per-call state and does not change results. |
| `t(df)` | `rt()` | `norm / sqrt(rchisq(df) / df)`, with `rchisq` = `rgamma(df/2, 2)`; used by PPC and response paths for t errors. |

The R code involved is GPL-2+. Ported files carry R's copyright notice, and the R Core Team is listed as copyright holder (`inst/COPYRIGHTS`, `Authors@R`).

**Streams.**

- **Chains.** As now, `.psp_streams(seed, chains)` runs `parallel::nextRNGStream()` in R and passes the six integers per chain. R stores them as signed `int`, and Rust reinterprets them as `u32`.
- **Inside Rust.** `next_stream()` and `next_substream()` are ports of `RngStream`'s jump matrices (A^(2^127) and A^(2^76) mod m), checked against `parallel::nextRNGStream()` and `nextRNGSubStream()`.
- **Post-processing.** `psp_ppc()`, `psp_response_path()` and `psp_loglik()` take one base stream from `.psp_streams(seed, 1)`. Rust gives posterior draw d (or GHK replication block b) the d-th stream after it, so results do not depend on the thread count. With `seed = NULL`, R draws a seed with `sample.int()` as `.psp_fit` does.
- **Periods (R6).** In the parallel latent sweep, period t of chain c uses substream t of chain c's stream.

**Behaviour changes** (go in NEWS):

- The sampler always uses inversion normals. It ignores a non-default `RNGkind(normal.kind = )`, which the C++ code honoured.
- `psp_ppc()`, `psp_response_path()` and `psp_loglik()` give different (but still reproducible) numbers for the same seed, because they move from R's default generator to package streams.
- Sampler draws for a given seed differ from the C++ version, because floating-point summation order and the AMD ordering of the RE Cholesky differ.

**Tests.** `test-rng.R` runs under `RNGkind("L'Ecuyer-CMRG", "Inversion")` with a fixed seed. It compares `psp_rs_rng()` with `runif()`, `rnorm()`, `rexp()`, `rgamma()` (shapes 0.3, 1, 2.5, 500) and `rt()` over 10^5 draws each, with `expect_identical`; if a platform's `libm` differs, allow 1 ulp. It also compares `next_stream` and `next_substream` with R's.

## 6.7 Sampler port (`pspcore::sampler`)

The block order and every formula follow `src/sampler.cpp` and `plan/02` exactly. Comments cite the plan sections, as the C++ does.

| Block (C++) | Rust | Changes during the port | Check against the frozen C++ build (§6.10) |
|---|---|---|---|
| `build_A`, `colsq` | `sparse::pattern` | For each distinct pair (W_t, M_t), precompute the merged pattern of I, W, M and MW and the four value vectors on it. An update is then `A.x = e_I - rho e_W - lambda e_M + rho lambda e_MW`, with column sums of squares in the same pass: no allocation, O(nnz). | Exact for `A = I - rho W` (same arithmetic as Eigen); 1e-14 for `R S`. |
| `bounds` | `latent.rs` | Precomputed `lo`/`hi` arrays. Binary bounds are fixed; ordered bounds are refreshed after each cut-point update. Branches leave the inner loop. | Same draws. |
| `latent_sweep` | `latent.rs` | Same loop, including the t weights and the latent-AR coupling to t + 1. | Same draws (1e-10 over the first iterations). |
| `compute_resid`, `apply_op` | `state.rs` | Shared buffers. For latent AR, the `lag_z` column is a separate vector filled from z, so `X` is not copied per chain. | 1e-12 |
| `draw_rho`, `draw_from_grid`, `lndet_interp` | `griddy.rs` | none | Same draws. |
| `refresh_design`, beta draw | `beta.rs` | Same explicit inverse for now. Replacing it with Cholesky solves is an R6 item, gated by SBC. | 1e-12 |
| lambda block | `lambda.rs` | none | |
| RE block, diagonal | `re.rs` | none | Same draws. |
| RE block, sparse (error term or spatial RE) | `re.rs` | faer `SymbolicLlt` once per chain. P values come from precomputed contribution plans for `sum_t R_t' D_t R_t` and `B'B / sigma^2`, so the t-weight case needs no sparse products. Numeric factorisation each iteration; the draw is the mean plus the solve of L' u = e, unpermuted. | Mean and precision 1e-12. Draws only distributionally, because faer's AMD ordering differs from Eigen's. |
| phi, sigma_alpha (half-Cauchy, inverse gamma) | `re.rs` | none | Same draws. |
| t mixing variances | `tmix.rs` | none | Same draws. |
| `cut_ac_step`, `cut_transform_move` with step adaptation | `cut.rs` | none | Same draws. |
| scale move | `scale.rs` | none | Same draws. |
| storage, thinning, `alpha_every`, acceptance rates | `store.rs` | Output layout identical to today's list. | identical |

**Threads, interrupts and progress** (`pspr::run`):

```
main thread (R)                      worker threads (std::thread::scope)
---------------                      -----------------------------------
copy inputs -> SamplerInput (Sync)   take chain c from the queue; own state + RMrg32k3a(seed_c)
spawn min(threads, chains) workers   every 100 iterations: check cancel flag, bump progress[c]
loop every 100 ms:                   finish -> send ChainOutput (or PspError) on a channel
  R_CheckUserInterrupt() inside
    savvy::unwind_protect            (workers never call R)
  -> on interrupt: set cancel, join, return the abort so R unwinds after Rust frames drop
  print "Chain c: x% (warm-up|sampling)" when a chain crosses 10% (also for parallel chains)
join; convert outputs to R lists
```

- **Thread count.** `cores` keeps its meaning as the number of chains run at once, but now on every OS. `parallel::mcmapply` and the `.Platform$OS.type == "unix"` test are removed.
- **CRAN limit.** When `_R_CHECK_LIMIT_CORES_` is set, `.psp_threads()` caps the count at 2, as `parallel` does.
- **Determinism.** Each chain is sequential with its own stream, so draws do not depend on the thread count. `test-fit.R` already tests this.

## 6.8 Post-processing port

**Reduced-form engine (`pspcore::rf`).** A `RfContext` is built once per call. It holds:

- the CSC matrices;
- the merged patterns of `S = I - rho W_u`, `R = I - lambda M_u` and `B = I - phi W`, for each distinct matrix u;
- one faer `SymbolicLu` per pattern.

Each worker owns its numeric LUs and solve workspaces (faer `MemBuffer`), allocated once. Per draw:

- **Factorise once.** Factorise S numerically once per distinct W_u. Solve `S x = b` and `S' x = b` with the same factors (faer transpose solve). Today's code factorises `t(S)` separately.
- **Exact path** (N ≤ `impacts_exact_max`): a dense inverse of S from faer's dense LU, then today's `rowSums` formulas.
- **Probe path:** today's single multi-RHS solve, with the same column layout (`xb`, `1`, `W1`, probes, `W probes`, `UG`, `UK`, `RV`), so it can be checked column by column.
- **Latent AR:** the recursions of `.psp_rf_latent_ar` (mean, lagged-error and random-effect variances; probe recursions with `P (T + 1)` columns).
- **t errors:** `s` becomes `n x Q` over the mixing variances passed in from `.psp_t_mix()`; averages follow `.psp_mix_mean`.

`pnorm`/`dnorm` are ports of R's `pnorm_both` (Cody), including the log versions, so probabilities agree with R to about 1e-15.

**Determinism.** Draws are split into fixed chunks of 8, whatever the thread count. Per-chunk partial sums (`predict()` accumulators) are added in chunk order, so results do not depend on the thread count.

**Per function:**

| Function | Rust part | Notes |
|---|---|---|
| `psp_impacts()` | Per draw: reduced form, then the density weights (`n x ncat`), then direct, total and indirect for each `x` variable (with theta from `durbin_of`), optionally by period | Returns today's array; `.psp_impacts_summary` stays in R. Unit type: `alpha` rows passed in. |
| `predict()` | Per draw: reduced form and category probabilities, accumulated | `back()` reordering stays in R. |
| `psp_ppc()` | Per draw: errors (normal or t), `R_t^{-1}`, the latent-AR term, `S_t^{-1}` with one LU per distinct matrix (today: one per period), then cut | Statistics stay in R. |
| `psp_response_path()` | Per draw: forward simulation with common random numbers (the same stream for baseline and shock), `.psp_lag_effect`, exact probabilities from `rf::block`; the latent-AR branch reuses the reduced form | Summary and `plot()` stay in R. |
| `psp_loglik()`, `logLik()`, `psp_dic()` | Per draw: block lower-bidiagonal A, `m = A^{-1}(R X beta)`, `Q = A'A` with the alpha block; ordering with the free coordinates first (AMD of their block), then AMD of `abs(Q_rr) + abs(C) abs(C)'` for the rest, as today; simplicial or supernodal Cholesky with that ordering (`SymmetricOrdering::Custom`); then the GHK loop | Parallel over draws when draws ≥ threads, otherwise over replications (each block of replications has its own stream). `.psp_ghk_check` stays in R. |

## 6.9 Build system and CRAN compliance

**DESCRIPTION**

- `SystemRequirements: Cargo (Rust's package manager), rustc >= 1.84`: CRAN's wording, plus the version.
- `Depends: R (>= 4.5)` (savvy).
- Remove `Rcpp` from Imports and `Rcpp` and `RcppEigen` from LinkingTo (R5); drop `importFrom(Rcpp, sourceCpp)`.
- `Authors@R`:
  - add `person("R Core Team", role = "cph", comment = "RNG and nmath code ported to Rust")`;
  - add `person("Authors of the dependency Rust crates", role = "cph", comment = "see inst/AUTHORS")`.

  `inst/AUTHORS` is generated by `tools/rust_authors.R`.
- Licence check: all 74 crates in the spike's vendor tree can be used under MIT, alone or as one of the options. `unicode-ident`, a build-time dependency, also carries the permissive Unicode-3.0 licence. All of these are compatible with GPL (≥ 3). A CI step (`cargo deny check licenses`, not vendored) fails on anything else.

**`configure` / `configure.win`** (portable `/bin/sh`):

1. Look for `cargo` and `rustc` on `PATH`, then in `~/.cargo/bin`. CRAN asks for this, and it is needed on this machine, where `~/.cargo/bin` is not on the non-interactive `PATH` and R comes from nix.
2. Print `rustc --version` and `cargo --version` before compiling (CRAN asks for the rustc version in the install log).
3. Stop with a clear message if `rustc` < 1.84 or cargo is missing. Suggest `rustup` or the system package manager; never install Rust from the script (CRAN forbids it).
4. If `src/rust/vendor.tar.xz` exists (a CRAN tarball), build offline from the vendored sources (`--offline`). Otherwise (a git checkout) use the registry with `--locked`.
5. Write `src/Makevars` from `Makevars.in`, substituting:
   - the cargo path;
   - the target directory: `PSP_CARGO_TARGET_DIR` if set, else `src/rust/target`;
   - the vendoring step;
   - the Rust target (Windows: `x86_64-pc-windows-gnu`).

*As built in R0* (`tools/config.sh`, which `configure` and `configure.win` both run): all five steps. With `src/rust/vendor.tar.xz` the build unpacks it, copies `src/rust/vendor-config.toml` to `src/.cargo/config.toml`, runs cargo with `--offline`, and removes both after linking; a git checkout uses the registry with `--locked`. `CARGO_HOME` is `src/.cargo` (relative to `src/`, where make runs) unless `PSP_CARGO_TARGET_DIR` is set and the build is not vendored: a developer build then keeps the user's own cargo home and registry cache. All paths in the Makevars are relative or come from `configure`, so the Makevars use no GNU make extensions (no `$(CURDIR)`, no `$(subst)`).

**`src/Makevars.in`** (sketch; portable make, no GNU extensions):

```make
LIBDIR = @TARGET_DIR@/@RUST_TARGET_SUBDIR@release
STATLIB = $(LIBDIR)/libpspr.a
PKG_LIBS = -L$(LIBDIR) -lpspr @EXTRA_LIBS@

all: $(SHLIB) @CLEAN_TARGET@

$(SHLIB): $(STATLIB)

$(STATLIB):
	@VENDOR_STEP@
	CARGO_HOME="$(CURDIR)/.cargo" "@CARGO@" build --lib --release -j 2 @OFFLINE@ \
	  @RUST_TARGET_FLAG@ --manifest-path=rust/Cargo.toml --target-dir="@TARGET_DIR@"

clean_target: $(SHLIB)
	rm -rf .cargo rust/vendor rust/target
```

- **`CARGO_HOME` inside the build tree.** Packages must not write to the home directory during installation.
- **`-j 2`** (CRAN).
- **`CLEAN_TARGET`** is set only for vendored builds with the in-tree target, so development builds stay incremental.
- **cargo runs on every build.** As in savvy's template, `libpspr.a` is deleted after linking, so make always calls cargo, which decides what to rebuild.
- **Registration (R0-R4).** `src/registration.c` `#include`s savvy's `init.c` with its `R_init` renamed and registers savvy's and Rcpp's routines together; the Makevars set `OBJECTS` without `init.o` (R-D5).
- **Windows.** `@EXTRA_LIBS@` is `-lws2_32 -ladvapi32 -luserenv -lbcrypt -lntdll`. Follow the current savvy Windows template for the linker setting and test on win-builder in R0. If CRAN adds Windows on ARM, add `aarch64-pc-windows-gnullvm`.
- **macOS.** CRAN's builder has Rust in `~/.cargo/bin`, which step 1 covers. Build for both macOS architectures CRAN uses, and pass R's `MACOSX_DEPLOYMENT_TARGET` through.
- **Portable code only.** No `-C target-cpu=native`: CRAN binaries run on many CPUs, and faer and pulp dispatch AVX2/AVX-512 at run time anyway.

**Vendoring (`tools/vendor.sh`)**, run before `R CMD build` for a release:

1. `cargo vendor` from `src/rust/Cargo.toml`.
2. Stub out crates that no CRAN target compiles. Today that is `windows-sys` and `windows-link`, which cut the tarball from 5.5 to 3.6 MB. `cargo vendor-filterer` with `--platform` for the CRAN targets can do it, or a script that empties the crates and rewrites their checksums.
3. `tar -cJf src/rust/vendor.tar.xz`.

The tarball is not committed. The whole source tarball should stay under 5 MB, and `R CMD build` must be followed by a check that it does.

**Checks that commonly give NOTEs for Rust packages**, to verify in R0 and explain in `cran-comments.md` if they appear:

- `checking compiled code`: `abort`/`exit` symbols from Rust's standard library;
- `installed package size`: the library is several MB, about 5 MB stripped in the spike.

**Threads policy.** The default stays at 1 (`psp_control(cores = 1)`), and tests and examples use at most 2.

**Development on this machine.**

- **Target directory.** `export PSP_CARGO_TARGET_DIR=$HOME/.cache/panelspatialprobit/target` keeps cargo's output (about 1 GB, many small files) off the Google Drive/`fuseblk` mount. That avoids sync churn and file-lock problems, and builds stay incremental: a change to `pspcore` alone rebuilds in seconds.
- **No more -O0 pitfall.** The Rust code is always built with the release profile from the Makevars, so the `-O0` problem of `roxygenise()`/`load_all()` no longer applies to it. Only the tiny `init.c` is affected.
- **Wrapper generation.** `savvy-cli update . src/rust/pspr` regenerates `R/000-wrappers.R`, `src/init.c` and `src/rust/api.h` after any change to a `#[savvy]` function. The second argument is needed because `src/rust/Cargo.toml` is a workspace, not the binding crate. It replaces `Rcpp::compileAttributes()`, and its version must match the savvy crate (0.11.0). Until R5, `src/registration.c` picks up the new routines from `init.c` by itself. Here `savvy-cli` is built from crates.io into `~/.cache/panelspatialprobit/tools/bin` (`cargo install savvy-cli --version 0.11.0 --locked --root ~/.cache/panelspatialprobit/tools`).
- **Lock file.** The workspace uses resolver 3, so `cargo generate-lockfile` and `cargo update` pick dependency versions whose `rust-version` is at most 1.84 (for example `generativity` 1.1.0, not 1.2.1, which needs 1.85).

**GitHub Actions.**

- Cache cargo builds (`Swatinem/rust-cache`) in the `rust` workflow only. R-CMD-check runs without a cache, so it builds as CRAN does (target and `CARGO_HOME` inside the build tree).
- On Windows, `rustup target add x86_64-pc-windows-gnu`.
- Add a `rust` job: `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, `cargo test -p pspcore`, `cargo deny check licenses`.
- Add an MSRV job: `cargo +1.84 build -p pspr`.
- Keep the three existing workflows; they first ran on 2026-10-30 (commit `05391bb`, C++ version). R-CMD-check and test-coverage also run on pushes to `rust`.

*As built in R0:* `.github/workflows/rust.yaml` holds the `rust` and `msrv` jobs.

## 6.10 Testing and validation

| Layer | What | Where | When |
|---|---|---|---|
| 1. Rust unit tests | RNG known-answer vectors (fixtures generated once from R); truncnorm moments and tail cases; griddy draws against exact quantiles; pattern updates against dense formulas; LU/Cholesky residuals; reduced form against dense algebra | `pspcore/tests`, `cargo test` | every commit (CI `rust` job) |
| 2. Memory and UB | `cargo miri test` on `rng`, `dist`, `sparse` (miri is installed here); `R -d valgrind` on a small fit, impacts and GHK call (valgrind is installed here) for the FFI layer | local | R1, R5 |
| 3. Cross-build equivalence (port only) | Two libraries, `lib_cpp` (frozen C++) and `lib_rust`, used in separate R processes with the same seeds. **Same draws:** `rtnorm` over a grid of cases; the cut kernel; the first 20-50 saved draws of short chains (`warmup = 0`) for none, sar, sdm and slx (pooled and diagonal RE, binary and ordered, t errors), to 1e-8. **Distributional:** sparse-RE variants (sem, sdem, sarar, spatial RE), posterior means within 4 MCSE (as `compare_to_reference`). **Post-processing:** a fit saved with `saveRDS` in one build is reloaded in the other; impacts, predict and logLik are compared (deterministic parts to 1e-10, GHK within MC error) | `validation/rust_port/` (not in the build; like `validation/crosscheck/`) | R2, R3 |
| 4. Existing suite, unchanged | reference sampler (`test-reference.R`), `test-cutmove.R`, `test-truncnorm.R`, finite-difference impacts, loglik against `mvtnorm` and quadrature, reproducibility and `cores` independence, restoration of the user's RNG state | `tests/testthat` | every commit |
| 5. New package tests | `test-rng.R` (§6.6); `test-threads.R` (sampler, impacts, predict, PPC, response and GHK identical with 1 and 2 threads); `test-postproc-oracle.R` (Rust reduced form, impacts and predict equal to `R/reference_postproc.R` to 1e-10 on small fits, including time-varying W, latent AR, t errors, spatial RE and both probe and exact paths); an input-validation error from Rust reaches R as an error, not a crash | `tests/testthat` | from R1 |
| 6. Manual | Ctrl-C during a long fit (serial and parallel chains) returns to the prompt promptly and leaves no threads; progress output with parallel chains; `psprobit()` inside `parallel::mclapply()` | checklist in the R2 PR | R2, R5 |
| 7. Statistical gate | full SBC for every variant (`validation/sbc/run_sbc.R`); Monte Carlo (`run_mc.R`); `spatialprobit` cross-checks (Katrina, ordered); a Wang and Kockelman subset; `bench_models.R` and `bench_large.R`; then `vignette_data.R` | `validation/` | R4, and again after any R6 change |
| 8. Performance gate | the targets in §6.2, recorded as `validation/bench/rust_<date>.txt` next to the C++ baseline | `validation/bench/` | R0 (baseline), R2, R3, R6 |

`R CMD check --as-cran` runs layers 4 and 5 (skipping the slow ones on CRAN, as now). Layers 1 and 2 run outside R. Layer 3 exists only while the C++ build is the reference, and is archived in R5.

## 6.11 Phases

Sizes are relative (S < M < L), as in `05`. R0 to R5 run in order, and each phase starts when the previous exit criterion is met. After R5, R6, R7 and the CRAN list can be done in any order.

### Action sequence

One step per row, in the order to do them.

- **After** lists the steps that must be finished first. Steps whose prerequisites are met can run alongside each other.
- **Who:** "You" is the maintainer. "Claude, when asked" marks git operations, which CLAUDE.md reserves for an explicit request; work on the `rust` branch is committed the same way, as steps finish.
- The phase checklists below give the details of each step.

| # | Phase | Step | After | Who |
|---|---|---|---|---|
| 1 | R0 | Check that the first CI run on `master` (C++ version, commit `05391bb`) passes, and fix it if not, so later failures can be blamed on Rust alone. *Done 2026-10-30: all three workflows pass on `968ea92`, after the pkgdown fixes.* | — | Claude |
| 2 | R0 | Tag the result `cpp-final` and create the branch `rust` from it. *Done 2026-10-30: `cpp-final` is the commit that records steps 1 and 2 (code identical to `968ea92`).* | 1 | You, or Claude when asked |
| 3 | R0 | Install `cpp-final` into a scratch library `lib_cpp`: the frozen reference. *Done 2026-10-30: `~/.cache/panelspatialprobit/lib_cpp`, from `git archive cpp-final` (source copy in `~/.cache/panelspatialprobit/src_cpp`).* | 2 | Claude |
| 4 | R0 | Record the C++ baseline: `bench_models.R`, `bench_large.R`, a callgrind profile of the sampler blocks (`validation/bench/rust_baseline_<date>.txt`). *Done 2026-10-30: `validation/bench/rust_baseline_2026-10-30.txt`, plus the new `bench_postproc.R`. callgrind stalls under this nix-built R, so the block split is by ablation.* | 3 | Claude |
| 5 | R0 | Workspace skeleton: `pspcore`, `pspr`, `configure`/`configure.win`/`cleanup`, `Makevars.in`/`Makevars.win.in`, savvy wrappers, one exported function; `.gitignore` and `.Rbuildignore` entries; `SystemRequirements` and `Depends: R (>= 4.5)`. *Done 2026-10-30: the function is `psp_rs_rng("unif", ...)` (MRG32k3a, equal to `runif()`; `test-rng.R`); `src/registration.c` (R-D5); vendoring moved here (below).* | 2 | Claude |
| 6 | R0 | Local `R CMD check --as-cran` on Linux; measure `codegen-units` and `lto`, install time and installed size. *Done 2026-10-30: see "R0 results" below.* | 5 | Claude |
| 7 | R0 | Add Rust to the GitHub workflows: cargo cache, Windows target, a `rust` job (fmt, clippy, tests, licences), an MSRV job; add `rust` to the `push` branches of R-CMD-check and test-coverage (they run on `main`/`master` and pull requests only). *Done 2026-10-30 (`rust.yaml`; R-CMD-check vendors the crates first); not yet run.* | 5 | Claude |
| 8 | R0 | Commit and push `rust`, so GitHub Actions checks the skeleton on Linux, macOS and Windows. *Done 2026-10-31 (`95ebe48`).* | 6, 7 | You, or Claude when asked |
| 9 | R0 | Fix platform failures until the R0 exit is met (each fix pushed as in step 8). *Done 2026-10-31: nothing to fix; every workflow passed on the first run, and R-CMD-check gives 1 NOTE (new submission) on all five configurations. R0 exit met.* | 8 | Claude |
| 10 | R1 | `rng`, with `test-rng.R` and Rust unit tests. *Done 2026-10-31: bit-identical to R on x86_64 for `runif`, `rnorm`, `rexp`, `rgamma` (shapes 0.3-500), `rt` and both stream jumps.* | 9 | Claude |
| 11 | R1 | `dist` (truncnorm, `pnorm`, `dnorm`); `psp_rs_rtnorm()` replaces `psp_rtnorm_cpp` in `test-truncnorm.R`; draws checked against `lib_cpp`. *Done 2026-10-31: also `qnorm`; `pnorm`/`dnorm`/`qnorm` bit-identical to R over both tails and log scales; truncated normal draws identical to the C++ sampler's (in-package test until R5).* | 10 | Claude |
| 12 | R1 | `sparse` (CSC, weights, pattern updates) and the R helpers `.psp_csc()` and `.psp_rs_weights()`. *Done 2026-10-31 (`test-rust-sparse.R`): `I - rho W` identical to Matrix, `(I - lambda M)(I - rho W)` within 1e-14.* | 9 | Claude |
| 13 | R1 | `par` and `pspr::run`: worker pool, cancel flag, interrupt polling, progress; interrupt check with `SIGINT`. *Done 2026-10-31 (`test-rust-par.R`, `validation/rust_port/check_interrupt.R`): a 4-thread job stops 2 ms after SIGINT with an R `interrupt` condition and no threads left.* | 9 | Claude |
| 14 | R1 | miri on `rng`, `dist` and `sparse`; R1 exit. *Done 2026-10-31: no findings (nightly toolchain with miri, `-Zmiri-deterministic-floats`); R1 exit met.* | 11, 12, 13 | Claude |
| 15 | R2 | Cross-build harness `validation/rust_port/`: `lib_cpp` against the branch build, in separate R processes with the same seeds. *Done 2026-10-31: `validation/rust_port/compare_sampler.R` (with `cases.R`, `fit_cases.R`): 17 same-draw cases and 8 sparse random-effects cases.* | 3, 14 | Claude |
| 16 | R2 | Switch `.psp_fit` to `psp_rs_sample()` and port binary SAR/SLX/none with pooled and diagonal RE; unported variants stop with "not yet ported" on the branch. *Done 2026-10-31. Unported variants ran on the C++ engine instead of stopping (`.psp_engine()`), so the branch stayed green; every variant runs on Rust since step 21. Same draws as C++ to 1e-11.* | 15 | Claude |
| 17 | R2 | Ordered outcomes, both cut-point samplers, scale move; `psp_rs_cut_kernel()` replaces `psp_cut_kernel_cpp` in `test-cutmove.R`. *Done 2026-10-31 (same draws as C++, both cut-point samplers; `test-cutmove.R` on `psp_rs_cut_kernel()`).* | 16 | Claude |
| 18 | R2 | t errors. *Done 2026-10-31 (same draws as C++).* | 17 | Claude |
| 19 | R2 | SDM (design path only). *Done 2026-10-31 (same draws as C++).* | 18 | Claude |
| 20 | R2 | SEM, SDEM and SARAR (lambda, sparse RE block). *Done 2026-10-31: pooled SEM/SDEM/SARAR give the same draws as C++; with random effects (sparse precision, faer simplicial Cholesky with AMD) posterior means agree within 2.1 MCSE (4,000 iterations).* | 19 | Claude |
| 21 | R2 | Spatial random effects (phi). *Done 2026-10-31 (posterior means within 2.4 MCSE, also with t errors).* | 20 | Claude |
| 22 | R2 | Time-varying weights. *Done 2026-10-31 (same draws as C++; with sparse random effects within 0.5 MCSE).* | 21 | Claude |
| 23 | R2 | Latent AR. *Done 2026-10-31 (same draws as C++ for SAR random effects and pooled SEM; SEM random effects within 1.5 MCSE).* | 22 | Claude |
| 24 | R2 | Full testthat suite with `NOT_CRAN=true`; `bench_models.R` against step 4; the SBC smoke subset; interrupt, progress and `mclapply` checks; R2 exit. *Done 2026-10-31: see "R2 results" below; R2 exit met.* | 23 | Claude |
| 25 | R3 | Choose how users set the post-processing threads. *Decided 2026-10-31: {future} (R-D11).* | — | You |
| 26 | R3 | Move the current R reduced form, impacts and predict code to `R/reference_postproc.R`; start `test-postproc-oracle.R` *Done 2026-10-31: `.psp_ref_*` (log-det, reduced form, latent AR, density weights, impacts and predict loops), verbatim.* | 24 | Claude |
| 27 | R3 | Log-determinant grids (sparse; eigen grid evaluation) *Done 2026-10-31: `pspcore::lndet` (simplicial LU per point; eigen grid sum).* | 24 | Claude |
| 28 | R3 | Reduced form (`rf`: block, t mixture, latent AR), with oracle tests *Done 2026-10-31: `pspcore::rf`; oracle cases: SAR/SEM-RE, SDM, SARAR, spatial RE, t errors, time-varying W, latent AR (also with SEM and t), exact and probe paths, with and without alpha.* | 26 | Claude |
| 29 | R3 | Impacts *Done 2026-10-31 (`pspcore::post`).* | 25, 28 | Claude |
| 30 | R3 | `predict()` *Done 2026-10-31: block sums of 8 draws, added in R in block order.* | 25, 28 | Claude |
| 31 | R3 | `psp_ppc()` replicates *Done 2026-10-31 (`pspcore::simulate::ppc_draw`; one stream per draw).* | 25, 28 | Claude |
| 32 | R3 | Response paths *Done 2026-10-31 (`pspcore::simulate::response_draw`; latent-AR branch via `psp_rs_probs`).* | 25, 28 | Claude |
| 33 | R3 | GHK likelihood *Done 2026-10-31 (`pspcore::ghk`; faer AMD on the free and restricted blocks, unordered simplicial Cholesky).* | 24, 25 | Claude |
| 34 | R3 | `test-threads.R`; post-processing benchmarks against §6.2; R3 exit *Done 2026-10-31: see "R3 results" below; R3 exit met except GHK with 4 workers (2.7-2.9x for ≥ 3x).* | 27, 29-33 | Claude |
| 35 | R4 | Full SBC for every variant, with replication runs for marginal results *Done 2026-10-31: 30 variants x 500, replications for 17; see "R4 results" below.* | 34 | Claude |
| 36 | R4 | Monte Carlo, `spatialprobit` cross-checks, Wang and Kockelman subset, benchmarks (alongside 35) *Done 2026-10-31 (the full Wang and Kockelman grid).* | 34 | Claude |
| 37 | R4 | `vignette_data.R`; compare the vignette numbers with the C++ ones and explain any change; R4 exit *Done 2026-10-31: R4 exit met.* | 35, 36 | Claude |
| 38 | R5 | Delete the C++ files; drop Rcpp and RcppEigen from DESCRIPTION and NAMESPACE; archive the `validation/rust_port/` results *Done 2026-10-31 (`validation/rust_port/README.md` marks `compare_sampler.R` historical).* | 37 | Claude |
| 39 | R5 | Update the documentation listed under R5 *Done 2026-10-31.* | 38 | Claude |
| 40 | R5 | Push; GitHub Actions clean on the three OSes (R5 exit) | 39 | You, or Claude when asked |
| 41 | R5 | Merge `rust` into `master` | 40 | You, or Claude when asked |
| 42 | after R5 | R6 optimisations, one at a time, each gated by ESS per second, the reference tests and SBC | 41 | Claude |
| 43 | after R5 | R7 extensions, as wanted | 41 | Claude, when asked |
| 44 | after R5 | "Before the CRAN submission" list, once you decide to submit | 41; CLAUDE.md items 1-5 | Claude (vendoring, `cran-comments.md`); you (win-builder, mac-builder, rhub, the submission) |

### R0: Decisions, skeleton and platform proof (M)

Steps 1-9.

- [x] Maintainer answers R-D1, R-D7 and R-D8. *Decided 2026-10-30: savvy, Rust 1.84, Rust before the first CRAN submission.*
- [x] The first CI run on `master` (C++ version) passes, or its failures are fixed first. *Done 2026-10-30: R-CMD-check and test-coverage passed on the first run; pkgdown failed (`psp_t` missing from the reference index) and passes since `968ea92`, which is the commit to tag.*
- [x] Tag the passing commit `cpp-final`; branch `rust` from it; install it into `lib_cpp` (scratch library). Tags, commits and pushes are made by the maintainer, or by Claude on request. *Done 2026-10-30.*
- [x] Profile the C++ sampler (`valgrind --tool=callgrind` on `bench_models.R` cases) and record the per-block breakdown and the baseline benchmark file. This can run alongside the skeleton work. *Done 2026-10-30, by ablation instead of callgrind (see step 4).*
- [x] Workspace skeleton: `pspcore`, `pspr`, `configure`/`configure.win`/`cleanup`, `Makevars.in`/`Makevars.win.in`, savvy wrappers, one exported function (e.g. `psp_rs_rng("unif", ...)`).
- [x] Measure the build options: `codegen-units` 1 vs 16 and `lto` off/thin/fat, against install time and the speed of a toy kernel. Record the installed size.
- [x] Rust in the GitHub workflows (§6.9).
- [x] Cross-platform proof: `R CMD check --as-cran` with the skeleton on Linux locally, and on GitHub Actions for macOS and Windows. The repository has existed since 2026-10-30 (`gueyenono/panelspatialprobit`).

**Exit:** the skeleton package passes `R CMD check --as-cran` on Linux, macOS and Windows (GitHub Actions), with only NOTEs that are understood; it installs in under 3 min with `-j 2`; decisions are recorded in this file.

**R0 results (2026-10-30, Linux, steps 5-6):**

- `R CMD check --as-cran` on the vendored tarball: **1 NOTE** (new submission), as for the C++ version. "checking compiled code" is OK (no `abort`/`exit` NOTE). The tests take 20 s.
- **Vendoring moved into R0.** R 4.6's check has a "checking Rust compilation" step that gives a WARNING ("Downloads Rust crates") when the install fetches crates, so every checked tarball must carry them. `tools/vendor.sh` (with the `windows-sys`/`windows-link` stubs) and the offline path in `tools/config.sh` are done; the R-CMD-check workflow runs `tools/vendor.sh` before building. Source tarball 3.7 MB (vendored crates 3.4 MB).
- **Install:** the Rust part builds in 68 s from cold at `-j 2`, peak memory about 1 GB; a full `R CMD INSTALL` takes about 1.5 min (39 s with a warm `PSP_CARGO_TARGET_DIR`).
- **Installed size:** reported as INFO, 27.7 MB, almost all debug information from this machine's `-g` (the C++-only library is already 22.6 MB); stripped, the library is 1.7 MB.
- **Release profile:** `codegen-units` 1 vs 16 and `lto` off/thin/fat make no difference to the toy kernel (sparse LU and 222-column solve at N = 10,000: 0.077-0.083 s; uniform draw 4.6-4.8 ns) and little to the build time (67-69 s; 76 s for 16 units with LTO). Kept: `codegen-units = 1`, `lto = "thin"`. Recheck once the sampler exists (R2 benchmarks).
- **Toolchain on this machine:** configure finds nix's cargo 1.95 on `PATH` before rustup's 1.97 in `~/.cargo/bin`; both are above 1.84.

### R1: Foundations (M)

Steps 10-14.

- [x] `rng`: MRG32k3a, inversion normal (AS241), `exp_rand`, `rgamma` (GD and GS), `rt`, `next_stream`, `next_substream`; `test-rng.R`.
- [x] `dist`: truncnorm, R's `pnorm` (including log) and `dnorm`. `psp_rs_rtnorm()` replaces `psp_rtnorm_cpp` in `test-truncnorm.R`.
- [x] `sparse`: `Csc`, `Weights`, conversion from `.psp_csc()`, mat-vec per period, merged patterns with value updates.
- [x] `par`: scoped worker pool, fixed chunking, cancel flag, progress counters. `pspr::run`: interrupt polling and progress printing.
- [x] Rust unit tests alongside each module; miri on `rng`/`dist`/`sparse`.

**Exit:** R's RNG streams reproduced (`expect_identical`, or ≤ 1 ulp where `libm` differs); `rtnorm` draws identical to `lib_cpp`; interrupting a dummy multi-thread job returns promptly on Linux.

**R1 results (2026-10-31):**

- **Exactness.** Every generator and normal function matches R bit for bit on x86_64 (`test-rng.R`; Rust known-answer tests in `pspcore/tests/r_reference.rs`). On other architectures `test-rng.R` compares to 1e-13, because a C compiler may fuse multiply-adds in R's own code (Apple clang does by default).
- **Cross-build exactness.** On GitHub's Ubuntu runners the Rust and C++ truncated normal draws differ in the last bit in the uniform-rejection branch: R's compiler flags there let g++ fuse `a + (b - a) * u` into a multiply-add, which Rust never does (R's own nmath was not affected: `test-rng.R` is exact there). The in-package Rust/C++ comparison therefore uses a 1e-14 tolerance, and the R2 cross-build checks compare with tolerances (§6.10), never bit for bit.
- **Registration.** `src/registration.c` now compiles savvy's `init.c` itself and copies its `CallEntries` at load time, so new `#[savvy]` functions need no manual edit; `src/Makevars.in` lists `OBJECTS` without `init.o` (this replaces the R0 `PKG_CPPFLAGS` rename).
- **Interrupt test harness.** A child started with `&` from a non-interactive shell inherits SIGINT as ignored, and R keeps that; `check_interrupt.R` forks with `parallel::mcparallel` instead, which also exercises Rust after `fork()`.
- **Authorship.** `tools/rust_authors.R` writes `inst/AUTHORS` (72 crates) from `cargo metadata`; `inst/COPYRIGHTS` lists the R Core code ported to Rust; `Authors@R` has both copyright-holder entries.
- **miri** needs the nightly toolchain with the `miri` and `rust-src` components (installed here 2026-10-31): `MIRIFLAGS=-Zmiri-deterministic-floats cargo +nightly miri test -p pspcore -- --skip par::` (about 1 min; the `par` tests are too slow under miri).

### R2: Sampler (L)

Steps 15-24. The first variant starts by switching `.psp_fit` to `psp_rs_sample()`:

- remove `mcmapply`, and `.psp_with_stream()` around the sampler;
- remove the `!parallel_run` restriction on `progress`;
- variants not yet ported stop with a "not yet ported" error on the branch.

Port in this order, each with its cross-build check (§6.7) before the next:

1. binary SAR/SLX/none, pooled and RE (diagonal);
2. ordered, with both cut-point samplers and the scale move; `psp_rs_cut_kernel()` replaces `psp_cut_kernel_cpp` in `test-cutmove.R`;
3. t errors;
4. SDM (no sampler change; check the design path);
5. SEM, SDEM and SARAR (lambda, sparse RE block);
6. spatial RE (phi);
7. time-varying weights;
8. latent AR.

**Exit:**

- cross-build checks pass for every variant;
- the full testthat suite passes with `NOT_CRAN=true`, including `test-reference.R`;
- `bench_models.R` per-iteration time within ±10% of C++ (or faster) for SAR, ordered SAR and SEM;
- the SBC smoke subset passes;
- the interrupt, progress and `mclapply` checks (§6.10 layer 6) pass.

**R2 results (2026-10-31):**

- **Engine switch.** `.psp_fit` builds `.psp_rs_model()` and a settings list and calls `psp_rs_sample()`; the C++ chain runner is kept as `.psp_run_cpp()`, used only when `options(panelspatialprobit.engine = "cpp")` (comparisons; removed in R5). Unported variants ran on C++ while the port went on, instead of stopping, so the branch stayed green throughout. Chains now run in parallel on every OS (`cores`), capped at 2 under `_R_CHECK_LIMIT_CORES_` (`.psp_threads()`).
- **Same draws.** `validation/rust_port/compare_sampler.R`: 17 cases give the C++ draws to at most 9e-12 (relative): SAR, SLX, none and SDM (pooled, random), missing outcomes, ordered with both cut-point samplers, t errors, pooled SEM, SDEM and SARAR, time-varying W, latent AR (with and without an error term), lagged outcomes.
- **Sparse random effects** (SEM, SDEM, SARAR with random effects; spatial random effects; with t errors, time-varying M and latent AR): faer's simplicial Cholesky with AMD ordering, symbolic analysis once per chain. Draws differ from C++ (another ordering), and posterior means over 4,000 iterations agree within 2.2 MCSE in all 8 cases. For t errors, P is assembled from precomputed contributions d_i R_ia R_ib, with no sparse products.
- **Tests.** The full suite with `NOT_CRAN=true` passes on the Rust sampler, including the comparison with the pure-R reference sampler for every variant. `test-cutmove.R` uses `psp_rs_cut_kernel()`.
- **Speed** (`validation/bench/rust_r2_2026-10-31.txt`): per iteration, C++ -> Rust: SAR 1.30 -> 1.37 ms (+5%), ordered 1.63 -> 1.77 (+9%), SEM 2.68 -> 2.43 (-9%); at N = 10,000, T = 20: SAR 32.5 -> 33.6 (+3%), ordered 41.0 -> 44.3 (+8%), SEM 72.1 -> 74.3 (+3%). Within the ±10% port target. Two changes made SEM faster than the C++: M X is computed once (only a latent AR lag column is refreshed), and without t weights X̃'X̃ comes from the Gram matrices X'X, X'MX and (MX)'MX; the random-effect right-hand side applies R'R once to the summed residuals when M is time-invariant. Dot products use four partial sums (vectorised), as Eigen's do.
- **SBC smoke subset** (100 replications: `sar_re_ord4`, `sarar_re_ord3`, `sar_sre_bin`, `sarar_t_re_ord3`, `sem_lat_re_ord3`): no rank test below the multiple-testing threshold (smallest p = 0.011 of 33); the three variants with p < 0.05 were rerun with `SBC_SEED_OFFSET=1` and 300 replications, all p > 0.05.
- **Threads, interrupts, fork** (`validation/rust_port/check_fit_threads.R`): Ctrl-C during a 4-chain, 4-thread fit stops within 0.1 s with an R interrupt and no threads left; progress is printed for every parallel chain; `psprobit()` inside `mclapply()` gives the parent's draws; `cores = 4` equals `cores = 1`. On this machine a forked child that calls `eigen()` after the parent did deadlocks in the BLAS (also in plain R, with no package code), so the check uses the sparse log-determinant; the C++ engine's `mcmapply` path had the same exposure.

**R3 results (2026-10-31):**

- **Parallelism (R-D11).** Rust post-processing calls are single-threaded; `.psp_post_map()` (`R/rust.R`) cuts the draws into blocks (1 draw, or 8 for `predict()`'s sums) and runs one `future::future()` per group of consecutive blocks, or calls Rust directly under the sequential plan. Per-draw results and block sums do not depend on the plan: `test-threads.R` checks impacts, unit impacts, predict, PPC, response paths and the GHK likelihood for identity between `sequential` and `multisession` with 2 workers (also checked by hand with `multicore`, 3 workers). Random numbers: one L'Ecuyer-CMRG stream per draw (PPC, response paths), or the same stream for every draw (GHK, common random numbers), from `.psp_rs_seeds()`. With `seed = NULL`, PPC and response paths draw their base seed from the current RNG state, which is left unchanged, as before; the replicates differ from the R version for the same seed.
- **Oracle** (`test-postproc-oracle.R`, 375 checks): Rust log-det grids, reduced form (every quantity, both paths, with and without alpha), impacts (population, unit, LeSage-Pace, by period) and predict match `R/reference_postproc.R` to 1e-10 (observed differences about 3e-16), for binary and ordered models.
- **Engine.** `pspcore::sparse::SparseLu`: faer LU with the column ordering computed once; `new()` lets faer choose a supernodal factor (2x faster than simplicial for the 202-column probe solve at N = 10,000), `new_simplicial()` exposes `U` for the log-determinant. `RfContext` holds one filter per distinct W and M and refactorises only when the coefficient changes. *Found:* faer 0.24.4's `SymmetricOrdering::Identity` misreads a lower-triangle input (panic in the symbolic Cholesky); `SparseChol::new_unordered()` passes the identity as a custom ordering.
- **Speed** (`validation/bench/rust_r3_2026-10-31.txt`): impacts at N = 1,000 0.61 s (C++/R 2.11 s), 0.24-0.33 s with 4 workers; predict 0.39 s (1.85); PPC replicates 0.25 s, whole call 0.47 s (3.71; the statistics stay in R); response paths 3.7 s (66 s), 1.0 s with 4 workers; GHK 2.2 s (2.8), 0.74-0.80 s with 4 workers (2.7-2.9x: ten draws over four workers is at most 3.3x); impacts at N = 10,000 4.1 s (19 s), 2.0 s with 4 multicore workers. Port changes behind the gains, all with identical results: no solves for probe columns that are not needed (impact columns outside impacts; UK when it equals UG), the response-path sd once per draw and weight matrix, one tail for one-sided GHK intervals.
- **On this machine** the nix `.libPaths()` is 55,536 characters, and parallelly launches all PSOCK workers in one command, which fails beyond 2 workers (the 128 KiB argument limit). `R_PARALLELLY_MAKENODEPSOCK_SETUP_STRATEGY=sequential` avoids it.

### R3: Post-processing (L)

Steps 25-34.

- [x] Decide how users set post-processing threads (maintainer). Decided 2026-10-31: {future} (R-D11); document `future::plan()` in `?psp_control` and in each post-processing function.
- [x] Move the current R reduced form, impacts and predict code to `R/reference_postproc.R` (oracle).
- [x] Port in this order, switching the R caller over at each step: `lndet` (sparse, and the eigen grid evaluation), `rf` (block, latent AR, mixture), then `impacts`, `predict`, `ppc`, `response` and `ghk`. Each port adds its cases to `test-postproc-oracle.R` and `test-threads.R` as it lands.

**Exit:** oracle tests pass (1e-10); the existing impacts, loglik, PPC and response tests pass; the §6.2 post-processing targets are met.

### R4: Validation gate (M; mostly compute time)

Steps 35-37.

- [x] Full SBC for every variant, with the same replication counts as the recorded runs, and a replication run (`SBC_SEED_OFFSET`) for any marginal result.
- [x] Monte Carlo, the `spatialprobit` cross-checks, the Wang and Kockelman subset (e.g. the rho = 0.6 and 0.9 cells), `bench_models.R`, `bench_large.R`, `bench_cutpoints.R`.
- [x] `validation/vignette_data.R`; compare the vignette numbers with the C++ ones and explain any change.

**Exit:** every documented variant passes SBC (CLAUDE.md model support gate), and the Monte Carlo coverage and bias stay within the ranges recorded in `05`.

**R4 results (2026-10-31):** all on the Rust engine (`validation/bench/rust_r4_2026-10-31.txt`); R4 exit met.

- **SBC.** All 30 variants, 500 replications each (`sbc_2026-10-31`): no failed fit, median max R-hat at most 1.017; 17 of 185 rank tests have p < 0.05, the smallest is `sar_lat_re_bin` `lag_z` (p = 0.0008, as with C++ at 0.001: slow mixing at the default chain length; `sar_lat_re_bin_long` passes), and none is below the Bonferroni threshold (2.7e-4). Replications (`SBC_SEED_OFFSET=1`, 500 each) for the 15 variants with any p < 0.05 and for `sem_lat_re_ord3`, and offset 2 for `sarar_t_re_ord3`: smallest p = 0.012 of 121, 4 below 0.05; `lag_z` gives p = 0.081. Many p-values in the first run repeat the C++ ones, since seeds, true values and data are the same and most variants give the C++ draws (R2).
- **Monte Carlo** (200 replications per design, `mc_2026-10-31`): every bias, RMSE and coverage within 0.01-0.02 of the C++ run (`mc_2026-10-27`); static designs biases at most 0.017, 90% coverage 0.88-0.92, within the `05` ranges; the dynamic designs as before (Wooldridge recovers `lag_y`, the pooled model overstates it by 0.31).
- **Wang and Kockelman** (all 12 cells, 50 data sets): average posterior means within 0.05 and coverages within 0.06 of the C++ run (`wk_2026-10-29`; overall coverage 0.89 in both); at rho = 0.9 64-80% of fits have max R-hat >= 1.05 (C++ 64-70%), Monte Carlo variation in a regime that needs longer chains.
- **Cross-checks.** `sarorderedprobit`: identical to the C++ record (same draws). Katrina: means differ from the C++ run by at most 0.03 posterior sd (Monte Carlo error; ESS about 1,000-2,600), and from `spatialprobit` by at most 0.59 sd, as before.
- **Benchmarks** (`models_`, `large_`, `cutpoints_2026-10-31`): per iteration against the C++ record, SAR 1.31 -> 1.35 ms, ordered 1.66 -> 1.75, SEM 2.69 -> 2.41; at N = 10,000, 32.5 -> 33.3, 40.7 -> 43.4, 74.4 -> 74.9; cut-point ESS ratios unchanged. A first `bench_models.R` run right after the 2-hour SBC and Monte Carlo jobs was 8-15% slower; C++ and Rust run alternately in one session gave the R2 gaps (SAR +4%, ordered +8%, SEM -8%), so it was rerun on an idle machine. `bench_large.R`: impacts for 50 draws 19 s -> 4 s.
- **Vignette** (`vignette_data.R` now keeps the newest run per variant and seed offset): every SBC row is from Rust; the text changed only where a number did (one rank test below 0.01 instead of three; Wooldridge coverage 0.92; Katrina 8x faster; the rho = 0.9 R-hat shares and RMSEs in the Wang and Kockelman section; SEM about twice the SAR cost).

### R5: Remove C++, documentation (S)

Steps 38-41.

- [x] Delete `src/sampler.cpp`, `src/ghk.cpp`, `src/truncnorm.h`, `src/RcppExports.cpp`, `R/RcppExports.R` and `src/registration.c`, and the `OBJECTS` line in both Makevars templates; update DESCRIPTION and NAMESPACE.
- [x] Archive `validation/rust_port/` results (keep the scripts, mark them historical).
- [x] Update the docs:
  - NEWS (behaviour changes in §6.6);
  - README (installing from source needs Rust);
  - `?psp_control` (`cores` on Windows, threads, the reproducibility contract);
  - `vignettes/validation.Rmd`;
  - `paper/paper.md`: "C++" at about lines 30, 181, 700 and 870, the software versions at about line 1762, and the Rcpp/RcppEigen references;
  - `plan/02` §2.13 (parallelism), `plan/03` §3.3-3.4 and D3, `plan/05`;
  - CLAUDE.md (commands, architecture).
- [ ] Merge `rust` into `master` once GitHub Actions is clean.

**Exit:** `R CMD check --as-cran` on the three OSes (GitHub Actions) gives only the expected NOTEs; the docs no longer mention the C++ engine.

### R6: Optimisations that Rust makes easy (M)

Step 42. Each item must show an ESS-per-second gain on the `02` §2.6 cases and pass the reference-sampler tests and SBC (CLAUDE.md "Mixing changes" convention). Read §2.6 "Mixing experiments" first.

- [ ] **Period-parallel latent sweep.** Without latent AR, the z_t are conditionally independent given everything else, and each period's sweep touches only its own residual segment. Sweeping periods in parallel is therefore the same Markov kernel. With latent AR, z_t depends only on z_{t-1} and z_{t+1}: sweep all even periods in parallel, then all odd ones (chromatic Gibbs). Each period writes only r_t and r_{t+1}, so there are no races. Substream per (chain, period) (§6.6). New `psp_control(threads_per_chain = 1)`.
- [ ] Specialise `latent_sweep` on (t errors, latent AR) with const generics.
- [ ] Cholesky solves instead of the explicit k x k inverse in the beta block.
- [ ] Fuse residual passes (`compute_resid` runs up to five times per iteration today).
- [ ] Rerun R4's gate on the final configuration.

**Exit:** the §6.2 "after R6" targets met, or re-baselined with a written reason; SBC passes.

### R7: Optional extensions

Step 43.

- Chromatic (graph-coloured) sweep within a period for T = 1 or small T. Units j and k can be updated at the same time when `(A'A)_jk = 0`, which keeps the kernel exact.
- Chebyshev or Barry-Pace Monte Carlo log-determinants for N beyond about 50,000 (`02` §2.9).
- A GHK variant with inverse-CDF draws, whose estimate is continuous in the parameters under common random numbers (smoother DIC differences).
- A faster RNG (R-D3 alternative), if profiling after R6 shows the RNG dominates the sweep.
- `pspcore` as a standalone crate (e.g. for Python bindings later).

### Before the CRAN submission (deferred, no date)

Step 44. These matter only for CRAN, so they wait for the submission (R-D8). They need R5 (not R6 or R7) and the maintainer items 1-5 in CLAUDE.md.

- [x] `tools/vendor.sh` with stubbing (§6.9); an offline install from the resulting tarball; the source tarball ≤ 5 MB. *Done in R0 (2026-10-30): 3.7 MB; check again before submitting.*
- [ ] Find the rustc version on each CRAN check flavour (install logs of existing Rust packages on the CRAN check pages) and confirm ≥ 1.84 everywhere, or revisit R-D7.
- [ ] win-builder (release and devel), mac-builder and rhub on the vendored tarball. They email the maintainer, so the maintainer starts them (CLAUDE.md item 5, after item 2).
- [ ] `cran-comments.md`: Rust, vendoring, installed size, any `compiled code` NOTE.

## 6.12 Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| The first CRAN submission is reviewed together with the Rust build (R-D8) | Medium | Prove the build on GitHub Actions for all three OSes in R0, before porting any real code; R5 exit needs clean checks there; the CRAN-only checks are in §6.11 "Before the CRAN submission". |
| A CRAN flavour has rustc < 1.84 | Low-Medium | R0 checks flavour logs; `configure` fails with a clear message; fallback is R-D7's alternative. |
| faer API changes (0.x) | High over time | Pin `=0.24.4`, commit `Cargo.lock`, upgrade deliberately with the full test and benchmark run. |
| Vendored size grows past the 5 MB tarball guidance | Medium | Stub unreachable crates; keep faer features minimal; CRAN allows hosting the bundle elsewhere with a checksum (last resort). |
| Install time or memory on CRAN machines (82 s and about 1 GB at `-j 2`) | Low | Within normal limits; measure `codegen-units`/`lto` in R0. |
| The port changes results subtly | Medium | Cross-build draw-by-draw checks, the deterministic post-processing oracle, the pure-R reference sampler, full SBC before C++ is removed. |
| Windows toolchain problems (gnu target, linker, system libraries) | Medium | Prove it in R0 with the skeleton, before any real code; follow savvy's template. |
| Hanging or slow interrupts with threads | Low | Workers poll a cancel flag every 100 iterations; the main thread polls every 100 ms; manual checklist (§6.10 layer 6). |
| Panics or UB at the FFI boundary | Low | `panic = "unwind"` with savvy's `catch_unwind`; `unsafe` forbidden in `pspcore` and confined to `pspr::interrupt`; miri and valgrind (§6.10 layer 2). |
| `fork()` after threads (mclapply) | Low | Scoped threads only, no global pool (R-D4); a manual test. |
| Cargo on the Google Drive mount (lock files, sync churn) | High here | `PSP_CARGO_TARGET_DIR` off the mount; `src/rust/target` in `.gitignore` and `.Rbuildignore`; or decision D5 (local working copy). |
| Maintainer unfamiliar with Rust | Medium | `pspcore` mirrors `sampler.cpp` block by block, with plan references in comments; rustdoc on every public item; this file kept current. |
| R ≥ 4.5 requirement excludes users | Low (no CRAN users yet) | Stated in R-D1; the C ABI alternative removes it. |

## 6.13 File checklist

| File | Change | Phase |
|---|---|---|
| `DESCRIPTION` | SystemRequirements, Depends R ≥ 4.5, Authors@R cph entries, drop Rcpp/RcppEigen | R0 (partly), R5 |
| `NAMESPACE` | Drop `importFrom(Rcpp, sourceCpp)`; keep `useDynLib(panelspatialprobit, .registration = TRUE)` | R5 |
| `configure`, `configure.win`, `cleanup` | New. Git mode 100755: on the NTFS mount every file looks executable, so set it explicitly with `git update-index --chmod=+x` | R0 |
| `tools/config.sh`, `cleanup.win` | New (the shared configure script; Windows cleanup) | R0 |
| `src/Makevars.in`, `src/Makevars.win.in`, `src/init.c`, `src/rust/**` | New | R0+ |
| `src/registration.c` | New in R0 (both engines' routines, R-D5); deleted in R5 | R0, R5 |
| `tools/rust_authors.R`, `inst/AUTHORS`, `inst/COPYRIGHTS` | New (done in R1) | R1 |
| `tools/vendor.sh`, `src/rust/vendor-config.toml` | New | R0 (R CMD check needs vendored crates) |
| `cran-comments.md` | New | before the CRAN submission |
| `.Rbuildignore` | Add `^src/rust/target$`, `^src/rust/vendor$`, `^src/\.cargo$`, `^src/Makevars$`, `^src/Makevars\.win$`, `^rust-dev$`, `^tools/vendor\.sh$`, `^cran-comments\.md$` | R0 |
| `.gitignore` | Add `src/rust/target/`, `src/rust/vendor/`, `src/rust/vendor.tar.xz`, `src/.cargo/`, `src/Makevars`, `src/Makevars.win`, `src/*.a` | R0 |
| `R/psprobit.R` | `.psp_fit` calls `psp_rs_sample()`; `mcmapply` and `.psp_with_stream()` for the sampler removed | R2 |
| `R/weights.R`, `R/impacts.R`, `R/methods.R`, `R/ppc.R`, `R/response.R`, `R/loglik.R` | Call the Rust functions; old code moved to `R/reference_postproc.R` | R3 |
| `R/rust.R`, `R/000-wrappers.R`, `R/reference_postproc.R` | New | R1-R3 |
| `src/*.cpp`, `src/truncnorm.h`, `R/RcppExports.R` | Delete | R5 |
| `tests/testthat/test-rng.R`, `test-threads.R`, `test-postproc-oracle.R` | New; `test-truncnorm.R` and `test-cutmove.R` switched to the Rust hooks | R1-R3 |
| `validation/rust_port/` | New (cross-build checks) | R1-R3 |
| `.github/workflows/*` | Rust cache, Windows target, `rust` job, MSRV job | R0 |
| `NEWS.md`, `README.md`, `vignettes/validation.Rmd`, `paper/paper.md`, `plan/02`, `plan/03`, `plan/05`, `CLAUDE.md` | §6.11 R5 | R5 |

**Commands after the migration** (for CLAUDE.md, R5):

```sh
export PATH="$HOME/.cargo/bin:$PATH"                                   # cargo is not on the default PATH here
export PSP_CARGO_TARGET_DIR="$HOME/.cache/panelspatialprobit/target"   # keep cargo output off the Drive
(cd src/rust && cargo test -p pspcore)                                 # Rust unit tests
(cd src/rust && cargo clippy --all-targets -- -D warnings && cargo fmt --check)
savvy-cli update . src/rust/pspr          # after changing #[savvy] functions (replaces Rcpp::compileAttributes)
R CMD INSTALL --no-test-load -l <lib> .   # build + install (cargo runs from the Makevars)
R CMD build . && R CMD check --as-cran --no-manual panelspatialprobit_*.tar.gz
tools/vendor.sh                           # CRAN tarballs only: writes src/rust/vendor.tar.xz before R CMD build
```
