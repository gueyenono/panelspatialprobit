//! Numerical core of panelspatialprobit (plan/06 §6.4).
//!
//! Pure Rust with no R dependency: everything here is testable with
//! `cargo test`. The R binding lives in the `pspr` crate.

#![forbid(unsafe_code)]
// Constants ported from R's C sources, and test values printed from R with
// 17 digits, are kept digit for digit.
#![allow(clippy::excessive_precision)]

pub mod error;
pub mod rng;

pub use error::PspError;
