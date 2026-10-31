// Port of pnorm5() and pnorm_both() in R's src/nmath/pnorm.c (Cody's
// rational Chebyshev approximations, with R's log and tail extensions).
// Copyright (C) 2000-2024 The R Core Team; GPL (>= 2).

const A: [f64; 5] = [
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113,
];
const B: [f64; 4] = [
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956,
];
const C: [f64; 9] = [
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8,
];
const D: [f64; 8] = [
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727,
];
const P: [f64; 6] = [
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303,
];
const Q: [f64; 5] = [
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5,
];

const M_SQRT_32: f64 = 5.656854249492380195206754896838;
const M_1_SQRT_2PI: f64 = 0.398942280401432677939946059934;

/// Which tail(s) `pnorm_both` computes, as R's `i_tail` (0, 1, 2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tail {
    Lower,
    Upper,
    Both,
}

/// `pnorm_both(x, &cum, &ccum, i_tail, log_p)`: returns `(cum, ccum)`, the
/// lower and upper tail probabilities of the standard normal at `x` (only the
/// requested ones are meaningful).
pub fn pnorm_both(x: f64, tail: Tail, log_p: bool) -> (f64, f64) {
    if x.is_nan() {
        return (x, x);
    }
    let lower = tail != Tail::Upper;
    let upper = tail != Tail::Lower;
    let eps = f64::EPSILON * 0.5;
    let (mut cum, mut ccum) = (0.0, 0.0);
    let y = x.abs();

    // R's do_del(X) followed by swap_tail.
    let del_swap = |xx: f64, temp: f64, cum: &mut f64, ccum: &mut f64| {
        let xsq = (xx * 16.0).trunc() / 16.0;
        let del = (xx - xsq) * (xx + xsq);
        if log_p {
            *cum = (-xsq * (xsq / 2.0)) - del / 2.0 + temp.ln();
            if (lower && x > 0.) || (upper && x <= 0.) {
                *ccum = (-(-xsq * (xsq / 2.0)).exp() * (-del / 2.0).exp() * temp).ln_1p();
            }
        } else {
            *cum = (-xsq * (xsq / 2.0)).exp() * (-del / 2.0).exp() * temp;
            *ccum = 1.0 - *cum;
        }
        if x > 0. {
            let t = *cum;
            if lower {
                *cum = *ccum;
            }
            *ccum = t;
        }
    };

    if y <= 0.67448975 {
        let (xnum, xden) = if y > eps {
            let xsq = x * x;
            let mut xnum = A[4] * xsq;
            let mut xden = xsq;
            for i in 0..3 {
                xnum = (xnum + A[i]) * xsq;
                xden = (xden + B[i]) * xsq;
            }
            (xnum, xden)
        } else {
            (0.0, 0.0)
        };
        let temp = x * (xnum + A[3]) / (xden + B[3]);
        if lower {
            cum = 0.5 + temp;
        }
        if upper {
            ccum = 0.5 - temp;
        }
        if log_p {
            if lower {
                cum = cum.ln();
            }
            if upper {
                ccum = ccum.ln();
            }
        }
    } else if y <= M_SQRT_32 {
        let mut xnum = C[8] * y;
        let mut xden = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);
        del_swap(y, temp, &mut cum, &mut ccum);
    } else if (log_p && y < 1e170)
        || (lower && -38.4674 < x && x < 8.2924)
        || (upper && -8.2924 < x && x < 38.4674)
    {
        let xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }
        let mut temp = xsq * (xnum + P[4]) / (xden + Q[4]);
        temp = (M_1_SQRT_2PI - temp) / y;
        del_swap(x, temp, &mut cum, &mut ccum);
    } else {
        let (d0, d1) = if log_p {
            (f64::NEG_INFINITY, 0.0)
        } else {
            (0.0, 1.0)
        };
        if x > 0. {
            cum = d1;
            ccum = d0;
        } else {
            cum = d0;
            ccum = d1;
        }
    }
    (cum, ccum)
}

/// `pnorm(x, mu, sigma, lower_tail, log_p)`.
pub fn pnorm(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        return x + mu + sigma;
    }
    let (d0, d1) = if log_p {
        (f64::NEG_INFINITY, 0.0)
    } else {
        (0.0, 1.0)
    };
    let (dt0, dt1) = if lower_tail { (d0, d1) } else { (d1, d0) };
    if !x.is_finite() && mu == x {
        return f64::NAN;
    }
    if sigma <= 0.0 {
        if sigma < 0.0 {
            return f64::NAN;
        }
        return if x < mu { dt0 } else { dt1 };
    }
    let p = (x - mu) / sigma;
    if !p.is_finite() {
        return if x < mu { dt0 } else { dt1 };
    }
    let (cum, ccum) = pnorm_both(p, if lower_tail { Tail::Lower } else { Tail::Upper }, log_p);
    if lower_tail {
        cum
    } else {
        ccum
    }
}
