//! Exact univariate truncated normal sampling, ported from the C++ engine's
//! `src/truncnorm.h` with the same algorithms and the same order of draws:
//! naive normal rejection when the support holds enough mass, Robert (1995)
//! exponential rejection for one-sided tails, and Robert's uniform rejection
//! for short two-sided intervals.

use crate::rng::Rng;

const SQRT_2PI: f64 = 2.5066282746310002;

/// Standard normal truncated to `[a, inf)`, `a` finite.
pub fn rtnorm_std_lower<R: Rng>(rng: &mut R, a: f64) -> f64 {
    if a <= 0.0 {
        // Acceptance probability >= 0.5.
        loop {
            let x = rng.norm();
            if x >= a {
                return x;
            }
        }
    }
    // Translated exponential proposal with the optimal rate.
    let lambda = 0.5 * (a + (a * a + 4.0).sqrt());
    loop {
        let x = a + rng.exp() / lambda;
        let d = x - lambda;
        if rng.unif() <= (-0.5 * d * d).exp() {
            return x;
        }
    }
}

/// Standard normal truncated to `[a, b]`, both finite, `a < b`.
pub fn rtnorm_std_two_sided<R: Rng>(rng: &mut R, mut a: f64, mut b: f64) -> f64 {
    if a < 0.0 && b > 0.0 {
        if b - a >= SQRT_2PI {
            // Normal rejection is efficient.
            loop {
                let x = rng.norm();
                if x >= a && x <= b {
                    return x;
                }
            }
        }
        // Uniform rejection, mode inside the interval.
        loop {
            let x = a + (b - a) * rng.unif();
            if rng.unif() <= (-0.5 * x * x).exp() {
                return x;
            }
        }
    }
    // Interval on one side of zero: work with a >= 0 and mirror if needed.
    let flip = b <= 0.0;
    if flip {
        let tmp = a;
        a = -b;
        b = -tmp;
    }
    // Robert's rule for choosing the uniform or the exponential proposal.
    let sq = (a * a + 4.0).sqrt();
    let bound = 2.0 * std::f64::consts::E.sqrt() / (a + sq) * (0.25 * (a * a - a * sq)).exp();
    let x = if b - a <= bound {
        loop {
            let x = a + (b - a) * rng.unif();
            if rng.unif() <= (0.5 * (a * a - x * x)).exp() {
                break x;
            }
        }
    } else {
        loop {
            let x = rtnorm_std_lower(rng, a);
            if x <= b {
                break x;
            }
        }
    };
    if flip {
        -x
    } else {
        x
    }
}

/// N(mean, sd^2) truncated to `[lo, hi]`; `lo` and `hi` may be infinite.
pub fn rtnorm<R: Rng>(rng: &mut R, mean: f64, sd: f64, lo: f64, hi: f64) -> f64 {
    let a = (lo - mean) / sd;
    let b = (hi - mean) / sd;
    let x = match (a.is_finite(), b.is_finite()) {
        (false, false) => rng.norm(),
        (true, false) => rtnorm_std_lower(rng, a),
        (false, true) => -rtnorm_std_lower(rng, -b),
        (true, true) => {
            if b <= a {
                return mean + sd * a; // degenerate interval
            }
            rtnorm_std_two_sided(rng, a, b)
        }
    };
    mean + sd * x
}
