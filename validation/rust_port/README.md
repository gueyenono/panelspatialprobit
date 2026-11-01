# Rust port checks (plan/06)

Scripts written while the sampler moved from C++ to Rust (phases R1-R4,
2026-10-30/31). Their results are recorded in `plan/06_rust_integration.md`
("R2 results" to "R4 results") and in `validation/bench/rust_*.txt`.

**Historical** (needs the frozen C++ build, tag `cpp-final`, installed in its
own library, e.g. `~/.cache/panelspatialprobit/lib_cpp`):

- `compare_sampler.R` (with `cases.R` and `fit_cases.R`): fits each case with
  the C++ and the Rust build and compares the draws (exact cases, to 1e-8) or
  the posterior means (sparse random-effect cases, within 4 MCSE).

**Still current** (the Rust build alone; unix):

- `check_interrupt.R`: Ctrl-C during a multi-threaded Rust job;
- `check_fit_threads.R`: Ctrl-C, progress, `mclapply()` and `cores` with the
  sampler.
