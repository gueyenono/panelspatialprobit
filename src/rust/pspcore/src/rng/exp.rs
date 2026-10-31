// Port of exp_rand() in R's src/nmath/sexp.c (Ahrens and Dieter, 1972).
// Copyright (C) 2000-2024 The R Core Team; GPL (>= 2).

use super::Rng;

// q[k-1] = sum(log(2)^k / k!), k = 1, ..., n; the highest n (16) is set by
// q[n-1] = 1.0 within double precision.
const Q: [f64; 16] = [
    0.6931471805599453,
    0.9333736875190459,
    0.9888777961838675,
    0.9984959252914960,
    0.9998292811061389,
    0.9999833164100727,
    0.9999985691438767,
    0.9999998906925558,
    0.9999999924734159,
    0.9999999995283275,
    0.9999999999728814,
    0.9999999999985598,
    0.9999999999999289,
    0.9999999999999968,
    0.9999999999999999,
    1.0000000000000000,
];

/// `exp_rand()`: a standard exponential draw.
pub(super) fn exp_rand<R: Rng>(rng: &mut R) -> f64 {
    let mut a = 0.0;
    let mut u = rng.unif();
    while u <= 0.0 || u >= 1.0 {
        u = rng.unif();
    }
    loop {
        u += u;
        if u > 1.0 {
            break;
        }
        a += Q[0];
    }
    u -= 1.0;

    if u <= Q[0] {
        return a + u;
    }

    let mut i = 0;
    let mut ustar = rng.unif();
    let mut umin = ustar;
    loop {
        ustar = rng.unif();
        if umin > ustar {
            umin = ustar;
        }
        i += 1;
        if u <= Q[i] {
            break;
        }
    }
    a + umin * Q[0]
}
