# 6. Rust integration

*Written 2026-10-30. Status: decisions R-D1, R-D7 and R-D8 taken by the maintainer on 2026-10-30 (savvy; Rust 1.84; Rust before the first CRAN submission); implementation not started. Measurements were taken on the development machine (8 cores, R 4.6.0, rustc 1.97.1) on 2026-10-30. The spike crates are not in the repository.*

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
| R-D5 | Migration strategy | **Feature branch `rust` plus a frozen C++ reference build.** Tag the last C++ commit (e.g. `cpp-final`), install it into its own library (`lib_cpp`), and compare across R processes (§6.10). This means two builds, not two backends. Rcpp's `RcppExports.cpp` and savvy's `init.c` would each define `R_init_panelspatialprobit`, and keeping both in one package would add work that is thrown away in R5. | A dual backend inside one package behind an internal option (needs merged native registration; useful only if the port must ship piecemeal). |
| R-D6 | Scope | §6.5. | Port everything numeric (rejected: loses the reference sampler as an oracle and changes `sim_psprobit()` data). |
| R-D7 | Minimum Rust version | **Decided: 1.84**, set by faer (released 2025-01-09). CRAN asks for testing with a cargo at least two years old, and preferably four; 1.84 turns two in January 2027. | An older faer line (older API, no fixes); dropping faer (R-D2 alternatives). |
| R-D8 | Order relative to the first CRAN submission | **Decided: the first CRAN submission is the Rust version**, after R5. The C++ version is never released on CRAN. Consequences: the submission waits for R0 to R5; the first-submission review also covers the Rust build, so R0's platform proof and R5's `cran-comments.md` matter more; the rhub and win-builder checks (CLAUDE.md item 5) are run on the Rust build. | Submit the C++ version first and ship Rust in the next release (the original recommendation: earlier submission, and the build system is not reviewed with a first submission). |
| R-D9 | Reproducibility contract | **Same seed and same machine and build give identical draws, whatever `cores`/threads**, as tested now (`test-fit.R`). Not promised: identity with the C++ version, or across CPUs. faer picks SIMD kernels at run time, and the eigen log-det path already depends on the BLAS. Documented in NEWS and `?psp_control`. | Bitwise identity across CPUs (would rule out faer in the sampler and the log-det; not worth it). |
| R-D10 | Data across the boundary | **R unpacks `dgCMatrix` into plain vectors** (`p`, `i`, `x`, `n`), and Rust copies its inputs into owned structures once per call. That costs O(nnz + nk), which is negligible against the work. Outputs are built in Rust and copied into R objects on the main thread. **No persistent external pointers**, since they do not survive `saveRDS()`. | Zero-copy borrowing of R memory (lifetimes tied to R objects; R memory touched from worker threads). |

## 6.4 Target layout

```
DESCRIPTION                 # SystemRequirements: Cargo (Rust's package manager), rustc >= 1.84
configure, configure.win    # find cargo, log versions, check MSRV, write src/Makevars(.win)
cleanup                     # remove generated Makevars, src/rust/target, src/.cargo, vendor dir
tools/
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
  rust/
    Cargo.toml              # workspace: pspcore, pspr; [profile.release]
    Cargo.lock              # committed
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
      Cargo.toml            # crate-type = ["staticlib"]; depends on savvy, pspcore
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
- **Wrapper generation.** `savvy-cli update .` regenerates `R/000-wrappers.R`, `src/init.c` and `src/rust/api.h` after any change to a `#[savvy]` function. It replaces `Rcpp::compileAttributes()`, and its version must match the savvy crate.

**GitHub Actions.**

- Cache cargo builds (`Swatinem/rust-cache`).
- On Windows, `rustup target add x86_64-pc-windows-gnu`.
- Add a `rust` job: `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings`, `cargo test -p pspcore`, `cargo deny check licenses`.
- Add an MSRV job: `cargo +1.84 build -p pspr`.
- Keep the three existing workflows. They have never run yet; that needs the GitHub repository (CLAUDE.md item 1).

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

Sizes are relative (S < M < L), as in `05`. Each phase starts when the previous exit criterion is met.

### R0: Decisions, skeleton and platform proof (M)

- [x] Maintainer answers R-D1, R-D7 and R-D8. *Decided 2026-10-30: savvy, Rust 1.84, Rust before the first CRAN submission.*
- [ ] Branch `rust` from the current commit; tag that commit `cpp-final`; install it into `lib_cpp` (scratch library). Commits and tags are made by the maintainer, or on request.
- [ ] Profile the C++ sampler (`valgrind --tool=callgrind` on `bench_models.R` cases) and record the per-block breakdown. Record the baseline benchmark file.
- [ ] Workspace skeleton: `pspcore`, `pspr`, `configure`/`configure.win`/`cleanup`, `Makevars.in`/`Makevars.win.in`, savvy wrappers, one exported function (e.g. `psp_rs_rng("unif", ...)`).
- [ ] Measure the build options: `codegen-units` 1 vs 16 and `lto` off/thin/fat, against install time and the speed of a toy kernel. Record the installed size.
- [ ] `tools/vendor.sh` with stubbing; confirm a source tarball ≤ 5 MB.
- [ ] Find the rustc version on each CRAN check flavour (install logs of existing Rust packages on the CRAN check pages) and confirm ≥ 1.84 everywhere, or revisit R-D7.
- [ ] Cross-platform proof: `R CMD check --as-cran` with the skeleton on Linux locally, and on GitHub Actions for macOS and Windows. **Needs the GitHub repository** (CLAUDE.md items 1 and 3). win-builder and mac-builder email the maintainer, so the maintainer starts them (item 5).

**Exit:** the skeleton package passes `R CMD check --as-cran` on Linux, macOS and Windows, with only NOTEs that are understood; it installs in under 3 min with `-j 2`; the tarball is ≤ 5 MB; decisions are recorded in this file.

### R1: Foundations (M)

- [ ] `rng`: MRG32k3a, inversion normal (AS241), `exp_rand`, `rgamma` (GD and GS), `rt`, `next_stream`, `next_substream`.
- [ ] `dist`: truncnorm, R's `pnorm` (including log) and `dnorm`.
- [ ] `sparse`: `Csc`, `Weights`, conversion from `.psp_csc()`, mat-vec per period, merged patterns with value updates.
- [ ] `par`: scoped worker pool, fixed chunking, cancel flag, progress counters. `pspr::run`: interrupt polling and progress printing.
- [ ] `test-rng.R`, Rust unit tests, miri on `rng`/`dist`/`sparse`.

**Exit:** R's RNG streams reproduced (`expect_identical`, or ≤ 1 ulp where `libm` differs); `rtnorm` draws identical to `lib_cpp`; interrupting a dummy multi-thread job returns promptly on Linux.

### R2: Sampler (L)

Port in this order, each with its cross-build check (§6.7) before the next:

1. binary SAR/SLX/none, pooled and RE (diagonal);
2. ordered, with both cut-point samplers and the scale move;
3. t errors;
4. SDM (no sampler change; check the design path);
5. SEM, SDEM and SARAR (lambda, sparse RE block);
6. spatial RE (phi);
7. time-varying weights;
8. latent AR.

- [ ] `.psp_fit` calls `psp_rs_sample()`; remove `.psp_with_stream()` around the sampler and `mcmapply`.
- [ ] Remove the `!parallel_run` restriction on `progress`.
- [ ] Replace `psp_rtnorm_cpp` and `psp_cut_kernel_cpp` in the tests with the Rust versions.

**Exit:**

- cross-build checks pass for every variant;
- the full testthat suite passes with `NOT_CRAN=true`, including `test-reference.R`;
- `bench_models.R` per-iteration time within ±10% of C++ (or faster) for SAR, ordered SAR and SEM;
- the SBC smoke subset passes.

### R3: Post-processing (L)

- [ ] Move the current R reduced form, impacts and predict code to `R/reference_postproc.R` (oracle).
- [ ] `lndet` (sparse and eigen grid evaluation), `rf` (block, latent AR, mixture), `impacts`, `predict`, `ppc`, `response`, `ghk`; R callers switched over one at a time.
- [ ] Decide how users set post-processing threads. Recommended: `fit$control$cores`, overridable with `options(panelspatialprobit.threads = )`. Document in `?psp_control`.
- [ ] `test-postproc-oracle.R`, `test-threads.R`.

**Exit:** oracle tests pass (1e-10); the existing impacts, loglik, PPC and response tests pass; the §6.2 post-processing targets are met.

### R4: Validation gate (M; mostly compute time)

- [ ] Full SBC for every variant, with the same replication counts as the recorded runs, and a replication run (`SBC_SEED_OFFSET`) for any marginal result.
- [ ] Monte Carlo, the `spatialprobit` cross-checks, the Wang and Kockelman subset (e.g. the rho = 0.6 and 0.9 cells), `bench_models.R`, `bench_large.R`, `bench_cutpoints.R`.
- [ ] `validation/vignette_data.R`; compare the vignette numbers with the C++ ones and explain any change.

**Exit:** every documented variant passes SBC (CLAUDE.md model support gate), and the Monte Carlo coverage and bias stay within the ranges recorded in `05`.

### R5: Remove C++, documentation (S)

- [ ] Delete `src/sampler.cpp`, `src/ghk.cpp`, `src/truncnorm.h`, `src/RcppExports.cpp` and `R/RcppExports.R`; update DESCRIPTION and NAMESPACE.
- [ ] Archive `validation/rust_port/` results (keep the scripts, mark them historical).
- [ ] Update the docs:
  - NEWS (behaviour changes in §6.6);
  - README (installing from source needs Rust);
  - `?psp_control` (`cores` on Windows, threads, the reproducibility contract);
  - `vignettes/validation.Rmd`;
  - `paper/paper.md`: "C++" at about lines 30, 181, 700 and 870, the software versions at about line 1762, and the Rcpp/RcppEigen references;
  - `plan/02` §2.13 (parallelism), `plan/03` §3.3-3.4 and D3, `plan/05`;
  - CLAUDE.md (commands, architecture).
- [ ] `cran-comments.md`: Rust, vendoring, installed size, any `compiled code` NOTE.

**Exit:** `R CMD check --as-cran` on the three OSes gives only the expected NOTEs; the docs no longer mention the C++ engine.

### R6: Optimisations that Rust makes easy (M)

Each item must show an ESS-per-second gain on the `02` §2.6 cases and pass the reference-sampler tests and SBC (CLAUDE.md "Mixing changes" convention). Read §2.6 "Mixing experiments" first.

- [ ] **Period-parallel latent sweep.** Without latent AR, the z_t are conditionally independent given everything else, and each period's sweep touches only its own residual segment. Sweeping periods in parallel is therefore the same Markov kernel. With latent AR, z_t depends only on z_{t-1} and z_{t+1}: sweep all even periods in parallel, then all odd ones (chromatic Gibbs). Each period writes only r_t and r_{t+1}, so there are no races. Substream per (chain, period) (§6.6). New `psp_control(threads_per_chain = 1)`.
- [ ] Specialise `latent_sweep` on (t errors, latent AR) with const generics.
- [ ] Cholesky solves instead of the explicit k x k inverse in the beta block.
- [ ] Fuse residual passes (`compute_resid` runs up to five times per iteration today).
- [ ] Rerun R4's gate on the final configuration.

**Exit:** the §6.2 "after R6" targets met, or re-baselined with a written reason; SBC passes.

### R7: Optional extensions

- Chromatic (graph-coloured) sweep within a period for T = 1 or small T. Units j and k can be updated at the same time when `(A'A)_jk = 0`, which keeps the kernel exact.
- Chebyshev or Barry-Pace Monte Carlo log-determinants for N beyond about 50,000 (`02` §2.9).
- A GHK variant with inverse-CDF draws, whose estimate is continuous in the parameters under common random numbers (smoother DIC differences).
- A faster RNG (R-D3 alternative), if profiling after R6 shows the RNG dominates the sweep.
- `pspcore` as a standalone crate (e.g. for Python bindings later).

## 6.12 Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| The first CRAN submission is reviewed together with the Rust build (R-D8) | Medium | Prove the build on CRAN-like platforms in R0 (GitHub Actions, then win-builder, mac-builder and rhub) before porting any real code; R5 exit needs clean checks on all three OSes; `cran-comments.md` explains the vendoring, the sizes and any NOTEs. |
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
| `src/Makevars.in`, `src/Makevars.win.in`, `src/init.c`, `src/rust/**` | New | R0+ |
| `tools/vendor.sh`, `tools/rust_authors.R`, `inst/AUTHORS`, `inst/COPYRIGHTS` | New | R0, R1 |
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
savvy-cli update .                        # after changing #[savvy] functions (replaces Rcpp::compileAttributes)
R CMD INSTALL --no-test-load -l <lib> .   # build + install (cargo runs from the Makevars)
tools/vendor.sh && R CMD build . && R CMD check --as-cran --no-manual panelspatialprobit_*.tar.gz
```
