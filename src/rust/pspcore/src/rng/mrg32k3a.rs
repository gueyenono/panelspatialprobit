// Port of the LECUYER_CMRG branch of unif_rand() in R's src/main/RNG.c, and
// of nextStream()/nextSubStream() in R's src/library/parallel/src/rngstream.c.
// Copyright (C) 1997-2024 The R Core Team; GPL (>= 2). Ported to Rust for
// panelspatialprobit (plan/06 §6.6).

use crate::PspError;

const M1: i64 = 4_294_967_087;
const M2: i64 = 4_294_944_443;
const NORMC: f64 = 2.328306549295727688e-10;
const A12: i64 = 1_403_580;
const A13N: i64 = 810_728;
const A21: i64 = 527_612;
const A23N: i64 = 1_370_589;

// Jump matrices A1^(2^76), A2^(2^76) (substreams) and A1^(2^127), A2^(2^127)
// (streams) mod m1 and m2, from rngstream.c.
const A1P76: [[u64; 3]; 3] = [
    [82758667, 1871391091, 4127413238],
    [3672831523, 69195019, 1871391091],
    [3672091415, 3528743235, 69195019],
];
const A2P76: [[u64; 3]; 3] = [
    [1511326704, 3759209742, 1610795712],
    [4292754251, 1511326704, 3889917532],
    [3859662829, 4292754251, 3708466080],
];
const A1P127: [[u64; 3]; 3] = [
    [2427906178, 3580155704, 949770784],
    [226153695, 1230515664, 3580155704],
    [1988835001, 986791581, 1230515664],
];
const A2P127: [[u64; 3]; 3] = [
    [1464411153, 277697599, 1610723613],
    [32183930, 1464411153, 1022607788],
    [2824425944, 32183930, 2093834863],
];

fn jump(s: &[u32; 6], a1: &[[u64; 3]; 3], a2: &[[u64; 3]; 3]) -> [u32; 6] {
    let mut out = [0u32; 6];
    for i in 0..3 {
        let mut tmp: u64 = 0;
        for j in 0..3 {
            // tmp < m1 < 2^32 and both factors are below 2^32, so the sum
            // stays below 2^64, as in R's uint_least64_t arithmetic.
            tmp += a1[i][j] * s[j] as u64;
            tmp %= M1 as u64;
        }
        out[i] = tmp as u32;
    }
    for i in 0..3 {
        let mut tmp: u64 = 0;
        for j in 0..3 {
            tmp += a2[i][j] * s[j + 3] as u64;
            tmp %= M2 as u64;
        }
        out[i + 3] = tmp as u32;
    }
    out
}

/// L'Ecuyer's MRG32k3a, exactly as R implements it for
/// `RNGkind("L'Ecuyer-CMRG")`.
///
/// The state is the six integers that follow the kind code in `.Random.seed`
/// (or in the output of `parallel::nextRNGStream()`). R stores them as signed
/// `int`; they are reinterpreted here as `u32`, as R's code does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RMrg32k3a {
    s: [u32; 6],
}

impl RMrg32k3a {
    /// Builds the generator from R's six seed integers, with R's validity
    /// check (`FixupSeeds`): the first three below m1, the last three below
    /// m2, and neither triple all zero.
    pub fn from_seed(seed: [i32; 6]) -> Result<Self, PspError> {
        let s = seed.map(|x| x as u32);
        let bad1 = s[..3].iter().any(|&x| x as i64 >= M1) || s[..3].iter().all(|&x| x == 0);
        let bad2 = s[3..].iter().any(|&x| x as i64 >= M2) || s[3..].iter().all(|&x| x == 0);
        if bad1 || bad2 {
            return Err(PspError::InvalidInput(
                "invalid L'Ecuyer-CMRG seed".to_string(),
            ));
        }
        Ok(Self { s })
    }

    /// The current state, in R's signed representation.
    pub fn seed(&self) -> [i32; 6] {
        self.s.map(|x| x as i32)
    }

    /// The generator at the start of the next stream, 2^127 steps on:
    /// `parallel::nextRNGStream()`.
    pub fn next_stream(&self) -> Self {
        Self {
            s: jump(&self.s, &A1P127, &A2P127),
        }
    }

    /// The generator at the start of the next substream, 2^76 steps on:
    /// `parallel::nextRNGSubStream()`.
    pub fn next_substream(&self) -> Self {
        Self {
            s: jump(&self.s, &A1P76, &A2P76),
        }
    }

    /// One uniform draw on (0, 1): `unif_rand()`.
    pub fn unif(&mut self) -> f64 {
        let s = &mut self.s;

        // `%` truncates like C's `k = p / m; p -= k * m`.
        let mut p1 = (A12 * s[1] as i64 - A13N * s[0] as i64) % M1;
        if p1 < 0 {
            p1 += M1;
        }
        s[0] = s[1];
        s[1] = s[2];
        s[2] = p1 as u32;

        let mut p2 = (A21 * s[5] as i64 - A23N * s[3] as i64) % M2;
        if p2 < 0 {
            p2 += M2;
        }
        s[3] = s[4];
        s[4] = s[5];
        s[5] = p2 as u32;

        (if p1 > p2 { p1 - p2 } else { p1 - p2 + M1 }) as f64 * NORMC
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // From R 4.6.0:
    //   RNGkind("L'Ecuyer-CMRG"); set.seed(42); s <- .Random.seed[-1]
    //   sprintf("%.17g", runif(5))
    const SEED_42: [i32; 6] = [
        -2_133_391_687,
        507_561_766,
        1_260_545_903,
        1_362_917_092,
        -1_772_566_379,
        -1_344_458_670,
    ];
    const UNIF_42: [f64; 5] = [
        0.17384558454153168,
        0.55474009676509084,
        0.48337712221370116,
        0.73748307381674638,
        0.79656476776243001,
    ];

    #[test]
    fn matches_r_runif() {
        let mut rng = RMrg32k3a::from_seed(SEED_42).unwrap();
        for &u in &UNIF_42 {
            assert_eq!(rng.unif().to_bits(), u.to_bits());
        }
    }

    #[test]
    fn rejects_invalid_seeds() {
        assert!(RMrg32k3a::from_seed([0, 0, 0, 1, 1, 1]).is_err());
        assert!(RMrg32k3a::from_seed([1, 1, 1, 0, 0, 0]).is_err());
        // 4294967087 = m1 is not a valid first-triple value.
        assert!(RMrg32k3a::from_seed([-209, 1, 1, 1, 1, 1]).is_err());
        assert!(RMrg32k3a::from_seed([1, 1, 1, 1, 1, 1]).is_ok());
    }
}
