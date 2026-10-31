// Port of dnorm4() in R's src/nmath/dnorm.c (the accurate version, with
// Welinder's split of x for |x| >= 5).
// Copyright (C) 2000-2024 The R Core Team; GPL (>= 2).

const M_1_SQRT_2PI: f64 = 0.398942280401432677939946059934;
const M_LN_SQRT_2PI: f64 = 0.918938533204672741780329736406;

/// `dnorm(x, mu, sigma, give_log)`.
pub fn dnorm(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return x + mu + sigma;
    }
    let d0 = if give_log { f64::NEG_INFINITY } else { 0.0 };
    if sigma < 0.0 {
        return f64::NAN;
    }
    if !sigma.is_finite() {
        return d0;
    }
    if !x.is_finite() && mu == x {
        return f64::NAN;
    }
    if sigma == 0.0 {
        return if x == mu { f64::INFINITY } else { d0 };
    }
    let x = (x - mu) / sigma;
    if !x.is_finite() {
        return d0;
    }
    let x = x.abs();
    if x >= 2.0 * f64::MAX.sqrt() {
        return d0;
    }
    if give_log {
        return -(M_LN_SQRT_2PI + 0.5 * x * x + sigma.ln());
    }
    if x < 5.0 {
        return M_1_SQRT_2PI * (-0.5 * x * x).exp() / sigma;
    }
    // sqrt(-2 log(2) (DBL_MIN_EXP + 1 - DBL_MANT_DIG)), DBL_MIN_EXP = -1021.
    if x > (-2.0 * std::f64::consts::LN_2 * (-1021.0 + 1.0 - 53.0)).sqrt() {
        return 0.0;
    }
    let x1 = (x * 65536.0).round_ties_even() / 65536.0;
    let x2 = x - x1;
    M_1_SQRT_2PI / sigma * ((-0.5 * x1 * x1).exp() * ((-0.5 * x2 - x1) * x2).exp())
}
