// Port of the INVERSION branch of norm_rand() in R's src/nmath/snorm.c.
// Copyright (C) 1998-2024 The R Core Team; GPL (>= 2).

use super::Rng;
use crate::dist::qnorm_std;

const BIG: f64 = 134_217_728.0; // 2^27

/// `norm_rand()` with `normal.kind = "Inversion"`. One uniform alone does not
/// have enough bits, so two are combined.
pub(super) fn norm_rand<R: Rng>(rng: &mut R) -> f64 {
    let u1 = rng.unif();
    let u1 = (BIG * u1) as i32 as f64 + rng.unif();
    qnorm_std(u1 / BIG)
}
