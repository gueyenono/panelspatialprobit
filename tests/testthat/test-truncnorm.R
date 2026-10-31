# Rust truncated normal sampler (pspcore::dist::truncnorm, plan/06 step 11).
rtn <- function(n, mean, sd, lo, hi, seed = 1) {
  panelspatialprobit:::psp_rs_rtnorm(as.integer(n), mean, sd, lo, hi,
                                     panelspatialprobit:::.psp_streams(seed, 1)[[1]][-1])
}

tn_mean <- function(mu, s, a, b) {
  al <- (a - mu) / s
  be <- (b - mu) / s
  mu + s * (dnorm(al) - dnorm(be)) / (pnorm(be) - pnorm(al))
}

test_that("truncated normal draws respect bounds and match analytic means", {
  cases <- list(c(0, 1, 0, Inf), c(0, 1, -Inf, -1), c(0, 1, -0.5, 0.7), c(1, 2, 3, 3.2),
                c(0, 1, 6, Inf), c(0, 1, -Inf, -8), c(0, 1, 2, 9), c(-3, 0.5, -10, -2.9),
                c(0, 1, -5, 5))
  for (cs in cases) {
    x <- rtn(20000, cs[1], cs[2], cs[3], cs[4])
    expect_true(all(x >= cs[3] & x <= cs[4]))
    expect_true(all(is.finite(x)))
    expect_lt(abs(mean(x) - tn_mean(cs[1], cs[2], cs[3], cs[4])), 0.02 * cs[2])
  }
})

test_that("extreme tails are finite", {
  x <- rtn(1000, 0, 1, 40, Inf)
  expect_true(all(is.finite(x) & x >= 40))
  x <- rtn(1000, 50, 1, -Inf, 0)
  expect_true(all(is.finite(x) & x <= 0))
})

# Until the C++ engine is removed (plan/06 R5): the Rust sampler gives the
# C++ one's draws, given the same L'Ecuyer-CMRG state and inversion normals.
test_that("Rust and C++ truncated normal draws are identical", {
  skip_if_not(R.version$arch == "x86_64", "bitwise identity is checked on x86_64")
  cases <- list(c(0, 1, -Inf, Inf), c(0, 1, 0, Inf), c(0, 1, -Inf, -1), c(0, 1, -0.5, 0.7),
                c(1, 2, 3, 3.2), c(0, 1, 6, Inf), c(0, 1, -Inf, -8), c(0, 1, 2, 9),
                c(-3, 0.5, -10, -2.9), c(0, 1, -5, 5), c(0, 1, 1, 1))
  for (cs in cases) {
    s <- panelspatialprobit:::.psp_streams(7, 1)[[1]]
    cpp <- panelspatialprobit:::.psp_with_rng_restore({
      RNGkind("L'Ecuyer-CMRG", "Inversion")
      assign(".Random.seed", s, envir = globalenv())
      panelspatialprobit:::psp_rtnorm_cpp(5000, cs[1], cs[2], cs[3], cs[4])
    })
    rs <- panelspatialprobit:::psp_rs_rtnorm(5000L, cs[1], cs[2], cs[3], cs[4], s[-1])
    expect_identical(rs, cpp)
  }
})
