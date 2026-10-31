//! Exact ports of R's random number generators (plan/06 §6.6), so that seeds
//! keep their meaning and draws can be compared with R's own.
//!
//! [`RMrg32k3a`] is R's L'Ecuyer-CMRG uniform generator. The [`Rng`] trait
//! builds R's other generators on any uniform source, as R's nmath does on
//! `unif_rand()`: normals by inversion (`norm_rand()` with
//! `normal.kind = "Inversion"`), `exp_rand()`, `rgamma()`, `rchisq()` and
//! `rt()`.

mod exp;
mod gamma;
mod mrg32k3a;
mod normal;

pub use mrg32k3a::RMrg32k3a;

/// R's generators on top of a uniform source.
pub trait Rng {
    /// One uniform draw on (0, 1): `unif_rand()`.
    fn unif(&mut self) -> f64;

    /// A standard normal draw by inversion: `norm_rand()`.
    fn norm(&mut self) -> f64
    where
        Self: Sized,
    {
        normal::norm_rand(self)
    }

    /// A standard exponential draw: `exp_rand()`.
    fn exp(&mut self) -> f64
    where
        Self: Sized,
    {
        exp::exp_rand(self)
    }

    /// A gamma draw with the given shape and scale: `rgamma(shape, scale)`.
    fn gamma(&mut self, shape: f64, scale: f64) -> f64
    where
        Self: Sized,
    {
        gamma::rgamma(self, shape, scale)
    }

    /// A chi-squared draw: `rchisq(df)`.
    fn chisq(&mut self, df: f64) -> f64
    where
        Self: Sized,
    {
        if !df.is_finite() || df < 0.0 {
            return f64::NAN;
        }
        self.gamma(df / 2.0, 2.0)
    }

    /// A Student-t draw: `rt(df)`.
    fn t(&mut self, df: f64) -> f64
    where
        Self: Sized,
    {
        if df.is_nan() || df <= 0.0 {
            return f64::NAN;
        }
        if !df.is_finite() {
            return self.norm();
        }
        // The numerator is drawn first, as in R.
        let num = self.norm();
        num / (self.chisq(df) / df).sqrt()
    }
}

impl Rng for RMrg32k3a {
    fn unif(&mut self) -> f64 {
        RMrg32k3a::unif(self)
    }
}
