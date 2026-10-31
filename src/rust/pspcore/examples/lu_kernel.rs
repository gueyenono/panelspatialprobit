//! Toy kernel for comparing release-profile settings (plan/06 R0): sparse LU
//! factorisations and solves of S = I - rho W, with W the row-standardised
//! rook matrix of a 100 x 100 grid (the impacts workload of §6.2), and a loop
//! of uniform draws.
//!
//!   cargo run --release --example lu_kernel

use faer::linalg::solvers::SolveCore;
use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
use faer::sparse::{SparseColMat, Triplet};
use pspcore::rng::RMrg32k3a;
use std::time::Instant;

fn rook_s(side: usize, rho: f64) -> SparseColMat<usize, f64> {
    let n = side * side;
    let mut trip = Vec::with_capacity(5 * n);
    for r in 0..side {
        for c in 0..side {
            let i = r * side + c;
            let mut nb = Vec::with_capacity(4);
            if r > 0 {
                nb.push(i - side);
            }
            if r + 1 < side {
                nb.push(i + side);
            }
            if c > 0 {
                nb.push(i - 1);
            }
            if c + 1 < side {
                nb.push(i + 1);
            }
            trip.push(Triplet::new(i, i, 1.0));
            let w = rho / nb.len() as f64;
            for j in nb {
                trip.push(Triplet::new(i, j, -w));
            }
        }
    }
    match SparseColMat::try_new_from_triplets(n, n, &trip) {
        Ok(s) => s,
        Err(e) => panic!("triplets: {e:?}"),
    }
}

fn main() {
    let side = 100;
    let reps = 20;
    let rhs = 222;
    let s0 = rook_s(side, 0.5);
    let symbolic = match SymbolicLu::try_new(s0.symbolic()) {
        Ok(s) => s,
        Err(e) => panic!("symbolic LU: {e:?}"),
    };

    let t0 = Instant::now();
    let mut check = 0.0;
    for k in 0..reps {
        let s = rook_s(side, 0.3 + 0.02 * k as f64);
        let lu = match Lu::try_new_with_symbolic(symbolic.clone(), s.as_ref()) {
            Ok(lu) => lu,
            Err(e) => panic!("numeric LU: {e:?}"),
        };
        let mut b = faer::Mat::<f64>::ones(side * side, rhs);
        lu.solve_in_place_with_conj(faer::Conj::No, b.as_mut());
        check += b[(0, 0)];
    }
    let lu_s = t0.elapsed().as_secs_f64() / reps as f64;

    let seed = [1, 2, 3, 4, 5, 6];
    let mut rng = match RMrg32k3a::from_seed(seed) {
        Ok(r) => r,
        Err(e) => panic!("{e}"),
    };
    let t1 = Instant::now();
    let mut sum = 0.0;
    for _ in 0..50_000_000 {
        sum += rng.unif();
    }
    let unif_ns = t1.elapsed().as_secs_f64() * 1e9 / 5e7;

    println!(
        "lu+solve {:.4} s/draw (N = {}, {} rhs); unif {:.2} ns/draw; check {:.6} {:.3}",
        lu_s,
        side * side,
        rhs,
        unif_ns,
        check,
        sum / 5e7
    );
}
