//! Numerical core of panelspatialprobit (plan/06 §6.4).
//!
//! Pure Rust with no R dependency: everything here is testable with
//! `cargo test`. The R binding lives in the `pspr` crate.

#![forbid(unsafe_code)]
// Constants ported from R's C sources, and test values printed from R with
// 17 digits, are kept digit for digit (including R's spelled-out log(2)).
// Index loops mirror the C/C++ code being ported and CSC column pointers.
#![allow(
    clippy::excessive_precision,
    clippy::approx_constant,
    clippy::needless_range_loop
)]

pub mod dist;
pub mod error;
pub mod par;
pub mod rng;
pub mod sparse;

pub use error::PspError;
