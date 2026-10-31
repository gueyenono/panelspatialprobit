//! Known-answer tests against R 4.6.0 (plan/06 §6.10, layer 1). The values
//! were printed from R with `sprintf("%.17g", ...)` after
//! `RNGkind("L'Ecuyer-CMRG", "Inversion"); set.seed(42)`; each generator
//! starts from that seed. The R-level tests (`test-rng.R`) compare 10^5 draws.

// R values are printed with 17 digits and kept as printed.
#![allow(clippy::excessive_precision, clippy::approx_constant)]

use pspcore::dist::{dnorm, pnorm, qnorm_std};
use pspcore::rng::{RMrg32k3a, Rng};

const SEED_42: [i32; 6] = [
    -2_133_391_687,
    507_561_766,
    1_260_545_903,
    1_362_917_092,
    -1_772_566_379,
    -1_344_458_670,
];

fn rng() -> RMrg32k3a {
    RMrg32k3a::from_seed(SEED_42).unwrap()
}

fn assert_bits(got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert_eq!(g.to_bits(), w.to_bits(), "got {g:e}, want {w:e}");
    }
}

#[test]
fn rnorm_inversion() {
    let mut r = rng();
    let got: Vec<f64> = (0..4).map(|_| r.norm()).collect();
    assert_bits(
        &got,
        &[
            -0.93907707864928736,
            -0.041679428023489604,
            0.82941349301613676,
            -0.43935819856743974,
        ],
    );
}

#[test]
fn rexp() {
    let mut r = rng();
    let got: Vec<f64> = (0..4).map(|_| r.exp()).collect();
    assert_bits(
        &got,
        &[
            1.777059037452144,
            0.10948019353018168,
            0.71696346396214339,
            1.0139511987908763,
        ],
    );
}

#[test]
fn rgamma_gs_and_gd() {
    let cases: [(f64, [f64; 3]); 4] = [
        (
            0.3,
            [
                0.0083136786750979703,
                0.25130511346442225,
                1.3283667551796301,
            ],
        ),
        (
            1.0,
            [0.11287733908271214, 2.1008782586437995, 3.8765483573808277],
        ),
        (
            2.5,
            [1.7848217983431849, 5.9997596005712666, 0.50819950213925025],
        ),
        (
            500.0,
            [957.46513285445428, 1027.6129046670587, 919.29284232364546],
        ),
    ];
    for (shape, want) in cases {
        let mut r = rng();
        let got: Vec<f64> = (0..3).map(|_| r.gamma(shape, 2.0)).collect();
        assert_bits(&got, &want);
    }
}

#[test]
fn rt() {
    let mut r = rng();
    let got: Vec<f64> = (0..3).map(|_| r.t(5.0)).collect();
    assert_bits(
        &got,
        &[-1.0656230093887549, -1.7247047180526254, 1.1281898976923819],
    );
}

#[test]
fn streams_and_substreams() {
    let r = rng();
    assert_eq!(
        r.next_stream().seed(),
        [
            -94575036,
            -24861725,
            357812871,
            346211554,
            -1301022835,
            440649733
        ]
    );
    assert_eq!(
        r.next_substream().seed(),
        [
            -1871565002,
            1968597656,
            -505666782,
            1055078269,
            658163990,
            1090264187
        ]
    );
}

const X: [f64; 10] = [-40.0, -9.0, -3.0, -0.5, 0.0, 1e-17, 0.3, 2.0, 7.0, 10.0];

#[test]
fn pnorm_lower_and_upper_log() {
    let got: Vec<f64> = X.iter().map(|&x| pnorm(x, 0.0, 1.0, true, false)).collect();
    assert_bits(
        &got,
        &[
            0.0,
            1.1285884059538408e-19,
            0.0013498980316300946,
            0.30853753872598694,
            0.5,
            0.5,
            0.61791142218895267,
            0.97724986805182079,
            0.99999999999872013,
            1.0,
        ],
    );
    let got: Vec<f64> = X.iter().map(|&x| pnorm(x, 0.0, 1.0, false, true)).collect();
    assert_bits(
        &got,
        &[
            -0.0,
            -1.1285884059538408e-19,
            -0.0013508099647481938,
            -0.36894641528865652,
            -0.69314718055994529,
            -0.69314718055994529,
            -0.96210281816885079,
            -3.7831843336820317,
            -27.384307498811076,
            -53.23128515051247,
        ],
    );
}

#[test]
fn dnorm_and_log() {
    let got: Vec<f64> = X.iter().map(|&x| dnorm(x, 0.0, 1.0, false)).collect();
    assert_bits(
        &got,
        &[
            0.0,
            1.0279773571668917e-18,
            0.0044318484119380075,
            0.35206532676429952,
            0.3989422804014327,
            0.3989422804014327,
            0.38138781546052414,
            0.053990966513188063,
            9.1347204083645936e-12,
            7.6945986267064199e-23,
        ],
    );
    let got: Vec<f64> = X.iter().map(|&x| dnorm(x, 0.0, 1.0, true)).collect();
    assert_bits(
        &got,
        &[
            -800.91893853320471,
            -41.418938533204674,
            -5.4189385332046731,
            -1.0439385332046727,
            -0.91893853320467278,
            -0.91893853320467278,
            -0.96393853320467282,
            -2.9189385332046727,
            -25.418938533204674,
            -50.918938533204674,
        ],
    );
}

#[test]
fn qnorm_all_branches() {
    let p = [1e-300, 1e-20, 1e-5, 0.07, 0.3, 0.5, 0.9, 0.99999];
    let got: Vec<f64> = p.iter().map(|&p| qnorm_std(p)).collect();
    assert_bits(
        &got,
        &[
            -37.047096299361201,
            -9.2623400897984052,
            -4.2648907939228256,
            -1.4757910281791702,
            -0.52440051270804067,
            0.0,
            1.2815515655446008,
            4.2648907939238399,
        ],
    );
}
