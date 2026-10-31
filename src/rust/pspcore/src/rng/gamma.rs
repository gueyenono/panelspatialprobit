// Port of rgamma() in R's src/nmath/rgamma.c: GD (Ahrens and Dieter, 1982)
// for shape >= 1 and GS (Ahrens and Dieter, 1974) for shape < 1.
// Copyright (C) 2000-2024 The R Core Team; GPL (>= 2).
//
// R caches the shape-dependent constants in static variables between calls;
// here they are recomputed per call, which gives the same values.

use super::Rng;

const SQRT32: f64 = 5.656854;
const EXP_M1: f64 = 0.36787944117144232159; // exp(-1)

// q0 = sum(q[k] a^(-k)); q = q0 + (t^2 / 2) sum(a[k] v^k).
const Q1: f64 = 0.04166669;
const Q2: f64 = 0.02083148;
const Q3: f64 = 0.00801191;
const Q4: f64 = 0.00144121;
const Q5: f64 = -7.388e-5;
const Q6: f64 = 2.4511e-4;
const Q7: f64 = 2.424e-4;

const A1: f64 = 0.3333333;
const A2: f64 = -0.250003;
const A3: f64 = 0.2000062;
const A4: f64 = -0.1662921;
const A5: f64 = 0.1423657;
const A6: f64 = -0.1367177;
const A7: f64 = 0.1233795;

/// `rgamma(a, scale)`: a gamma draw with shape `a` and scale `scale`.
pub(super) fn rgamma<R: Rng>(rng: &mut R, a: f64, scale: f64) -> f64 {
    if a.is_nan() || scale.is_nan() {
        return f64::NAN;
    }
    if a <= 0.0 || scale <= 0.0 {
        if scale == 0.0 || a == 0.0 {
            return 0.0;
        }
        return f64::NAN;
    }
    if !a.is_finite() || !scale.is_finite() {
        return f64::INFINITY;
    }

    if a < 1.0 {
        // GS algorithm for a < 1.
        let e = 1.0 + EXP_M1 * a;
        let x = loop {
            let p = e * rng.unif();
            if p >= 1.0 {
                let x = -((e - p) / a).ln();
                if rng.exp() >= (1.0 - a) * x.ln() {
                    break x;
                }
            } else {
                let x = (p.ln() / a).exp();
                if rng.exp() >= x {
                    break x;
                }
            }
        };
        return scale * x;
    }

    // GD algorithm for a >= 1.
    // Step 1.
    let s2 = a - 0.5;
    let s = s2.sqrt();
    let d = SQRT32 - s * 12.0;

    // Step 2: t is a standard normal deviate, x a (s, 1/2)-normal deviate;
    // immediate acceptance.
    let mut t = rng.norm();
    let mut x = s + 0.5 * t;
    let ret_val = x * x;
    if t >= 0.0 {
        return scale * ret_val;
    }

    // Step 3: squeeze acceptance.
    let mut u = rng.unif();
    if d * u <= t * t * t {
        return scale * ret_val;
    }

    // Step 4.
    let r = 1.0 / a;
    let q0 = ((((((Q7 * r + Q6) * r + Q5) * r + Q4) * r + Q3) * r + Q2) * r + Q1) * r;
    let (b, si, c) = if a <= 3.686 {
        (0.463 + s + 0.178 * s2, 1.235, 0.195 / s - 0.079 + 0.16 * s)
    } else if a <= 13.022 {
        (1.654 + 0.0076 * s2, 1.68 / s + 0.275, 0.062 / s + 0.024)
    } else {
        (1.77, 0.75, 0.1515 / s)
    };

    let quotient = |t: f64| -> f64 {
        let v = t / (s + s);
        if v.abs() <= 0.25 {
            q0 + 0.5
                * t
                * t
                * ((((((A7 * v + A6) * v + A5) * v + A4) * v + A3) * v + A2) * v + A1)
                * v
        } else {
            q0 - s * t + 0.25 * t * t + (s2 + s2) * (1.0 + v).ln()
        }
    };

    // Steps 5-7: quotient acceptance, only if x is positive.
    if x > 0.0 {
        let q = quotient(t);
        if (1.0 - u).ln() <= q {
            return scale * ret_val;
        }
    }

    // Steps 8-11: double-exponential rejection.
    loop {
        let e = rng.exp();
        u = rng.unif();
        u = u + u - 1.0;
        t = if u < 0.0 { b - si * e } else { b + si * e };
        if t >= -0.71874483771719 {
            let q = quotient(t);
            if q > 0.0 {
                let w = q.exp_m1();
                if c * u.abs() <= w * (e - 0.5 * t * t).exp() {
                    break;
                }
            }
        }
    }
    x = s + 0.5 * t;
    scale * x * x
}
