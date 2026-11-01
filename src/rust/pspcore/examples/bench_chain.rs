//! Sampler timings without R (plan/06 R6): the three models of
//! `validation/bench/bench_models.R` (SAR binary, SAR ordered with J = 4, SEM
//! binary; random effects, k = 5) on a rook grid, one chain on one thread.
//! Prints ms per iteration and the posterior means, so that a change to the
//! sampler can be timed and checked against the previous build.
//!
//!   PSP_THREADS=1 cargo run --release --example bench_chain [side_r side_c T iter [model ...]]
//!
//! Defaults: 25 x 40 grid, T = 10, 2,000 + 2,000 iterations, all three
//! models; `PSP_THREADS` threads for the latent sweep.

use pspcore::lndet::lndet_sparse;
use pspcore::par;
use pspcore::rng::{RMrg32k3a, Rng};
use pspcore::sampler::{run_chain, CutSampler, GridPrior, KappaPrior, SamplerInput, SigmaPrior};
use pspcore::sparse::{Csc, Weights};
use std::time::{Duration, Instant};

fn rook(nr: usize, nc: usize) -> Csc {
    let mut trip = Vec::new();
    for r in 0..nr {
        for c in 0..nc {
            let i = r * nc + c;
            let mut nb = Vec::new();
            if r > 0 {
                nb.push(i - nc);
            }
            if r + 1 < nr {
                nb.push(i + nc);
            }
            if c > 0 {
                nb.push(i - 1);
            }
            if c + 1 < nc {
                nb.push(i + 1);
            }
            let w = 1.0 / nb.len() as f64;
            for j in nb {
                trip.push((i, j, w));
            }
        }
    }
    Csc::from_triplets(nr * nc, nr * nc, &trip).unwrap()
}

/// z = (I - rho W)^{-1} (X beta + alpha + u), u = (I - lambda W)^{-1} eps, by
/// fixed-point iteration (|rho|, |lambda| < 1, W row-standardised).
fn solve_filter(w: &Csc, c: f64, b: &[f64]) -> Vec<f64> {
    let mut z = b.to_vec();
    let mut wz = vec![0.0; b.len()];
    for _ in 0..200 {
        w.mul_vec(&z, &mut wz);
        for i in 0..b.len() {
            z[i] = b[i] + c * wz[i];
        }
    }
    z
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: usize| args.get(i).map_or(d, |s| s.parse().unwrap());
    let (nr, nc, tt, iter) = (num(0, 25), num(1, 40), num(2, 10), num(3, 2000));
    let threads = std::env::var("PSP_THREADS").map_or(1, |s| s.parse().unwrap());
    let models: Vec<String> = if args.len() > 4 {
        args[4..].to_vec()
    } else {
        vec!["sar".into(), "ord".into(), "sem".into()]
    };
    let w = rook(nr, nc);
    let nn = nr * nc;
    let n = nn * tt;
    let k = 5;
    let grid: Vec<f64> = (0..=198).map(|i| -0.99 + 0.01 * i as f64).collect();
    let t0 = Instant::now();
    let lndet = lndet_sparse(&w, &grid).unwrap();
    eprintln!("lndet grid: {:.2} s", t0.elapsed().as_secs_f64());
    let gp = GridPrior {
        a: 1.0,
        b: 1.0,
        min: -1.0,
        max: 1.0,
        grid,
        lndet,
    };
    let weights = Weights::new(vec![w.clone()], vec![0; tt]).unwrap();

    for model in &models {
        let (sem, ordered) = (model == "sem", model == "ord");
        let mut rng = RMrg32k3a::from_seed([12345; 6]).unwrap();
        let beta = [0.0, 1.0, -1.0, 0.5, 0.5];
        let mut x = vec![1.0; n * k];
        for v in x[n..].iter_mut() {
            *v = rng.norm();
        }
        let alpha: Vec<f64> = (0..nn).map(|_| 0.7 * rng.norm()).collect();
        let mut zs = vec![0.0; n];
        for t in 0..tt {
            let r = t * nn..(t + 1) * nn;
            let mut b: Vec<f64> = (0..nn)
                .map(|i| {
                    let idx = t * nn + i;
                    (0..k).map(|j| x[j * n + idx] * beta[j]).sum::<f64>() + alpha[i]
                })
                .collect();
            let eps: Vec<f64> = (0..nn).map(|_| rng.norm()).collect();
            let u = if sem {
                solve_filter(&w, 0.5, &eps)
            } else {
                eps
            };
            for i in 0..nn {
                b[i] += u[i];
            }
            let z = if sem { b } else { solve_filter(&w, 0.5, &b) };
            zs[r].copy_from_slice(&z);
        }
        let y: Vec<i32> = zs
            .iter()
            .map(|&z| {
                if ordered {
                    1 + (z > 0.0) as i32 + (z > 0.8) as i32 + (z > 1.6) as i32
                } else {
                    (z > 0.0) as i32
                }
            })
            .collect();
        let inp = SamplerInput {
            y,
            x: std::mem::take(&mut x),
            k,
            n_units: nn,
            periods: tt,
            n_cat: if ordered { 4 } else { 2 },
            ordered,
            has_rho: !sem,
            has_lambda: sem,
            has_re: true,
            spatial_re: false,
            w: weights.clone(),
            m: weights.clone(),
            beta_mean: vec![0.0; k],
            beta_prec: (0..k * k)
                .map(|i| if i % (k + 1) == 0 { 1e-4 } else { 0.0 })
                .collect(),
            rho: gp.clone(),
            lambda: gp.clone(),
            phi_a: 1.0,
            phi_b: 1.0,
            nu: 0.0,
            lag_col: None,
            sigma_prior: SigmaPrior::HalfCauchy { scale: 1.0 },
            kappa_prior: KappaPrior::Flat,
            iter,
            warmup: iter,
            thin: 1,
            z_sweeps: 1,
            threads,
            cut_sampler: CutSampler::Transform,
            scale_move: true,
            beta_start: vec![0.0; k],
            rho_start: 0.0,
            lambda_start: 0.0,
            sigma2_start: 1.0,
            kappa_start: if ordered { vec![1.0, 2.0] } else { Vec::new() },
            alpha_every: 0,
        };
        let seed = [1, 2, 3, 4, 5, 6];
        let t0 = Instant::now();
        let out = par::run(
            1,
            1,
            Duration::from_millis(100),
            |_, ctx| {
                let mut rng = RMrg32k3a::from_seed(seed)?;
                run_chain(&inp, &mut rng, ctx)
            },
            |_| Ok::<(), ()>(()),
        )
        .unwrap()
        .pop()
        .unwrap();
        let el = t0.elapsed().as_secs_f64();
        let means: Vec<String> = (0..out.n_par)
            .map(|p| {
                let col = &out.draws[p * out.nsave..(p + 1) * out.nsave];
                format!("{:.4}", col.iter().sum::<f64>() / out.nsave as f64)
            })
            .collect();
        println!(
            "{model} N={nn} T={tt} threads={threads}: {:.3} ms/iter; means {}",
            1000.0 * el / (2 * iter) as f64,
            means.join(" ")
        );
    }
}
