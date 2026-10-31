//! Exact ports of R's random number generators (plan/06 §6.6), so that seeds
//! keep their meaning and draws can be compared with R's own.

mod mrg32k3a;

pub use mrg32k3a::RMrg32k3a;
