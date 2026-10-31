# The Rust ports of R's generators and normal distribution functions
# reproduce R's values exactly (plan/06 §6.6). Bitwise identity is checked on
# x86_64; other platforms may differ in the last bits where the C compiler
# fuses multiply-adds or libm differs, so there the check is to 1e-13.

rs_rng <- function(kind, n, par, seed) {
  panelspatialprobit:::psp_rs_rng(kind, as.integer(n), as.double(par), seed)
}

# R's draws under L'Ecuyer-CMRG and inversion normals, and the seed they
# started from (without the kind code). `expr` is a promise, evaluated after
# set.seed().
r_draws <- function(seed, expr) {
  panelspatialprobit:::.psp_with_rng_restore({
    RNGkind("L'Ecuyer-CMRG", "Inversion")
    set.seed(seed)
    s <- .Random.seed[-1]
    list(seed = s, x = expr)
  })
}

expect_r_identical <- function(object, expected) {
  if (R.version$arch == "x86_64") expect_identical(object, expected)
  else expect_equal(object, expected, tolerance = 1e-13)
}

test_that("uniform, normal and exponential draws match R", {
  for (seed in c(1L, 42L, 2026L)) {
    r <- r_draws(seed, runif(1e5))
    expect_identical(rs_rng("unif", 1e5, numeric(0), r$seed), r$x)
    r <- r_draws(seed, rnorm(1e5))
    expect_r_identical(rs_rng("norm", 1e5, numeric(0), r$seed), r$x)
    r <- r_draws(seed, rexp(1e5))
    expect_r_identical(rs_rng("exp", 1e5, numeric(0), r$seed), r$x)
  }
  expect_identical(rs_rng("unif", 0, numeric(0), 1:6), numeric(0))
})

test_that("gamma and t draws match R", {
  for (shape in c(0.3, 1, 2.5, 500)) {
    r <- r_draws(3L, rgamma(1e5, shape, scale = 0.7))
    expect_r_identical(rs_rng("gamma", 1e5, c(shape, 0.7), r$seed), r$x)
  }
  for (df in c(1, 4, 30)) {
    r <- r_draws(4L, rt(1e5, df))
    expect_r_identical(rs_rng("t", 1e5, df, r$seed), r$x)
  }
})

test_that("stream and substream jumps match parallel's", {
  s <- r_draws(5L, NULL)$seed
  for (i in 1:5) {
    expect_identical(panelspatialprobit:::psp_rs_next_stream(s, FALSE),
                     parallel::nextRNGStream(c(10407L, s))[-1])
    expect_identical(panelspatialprobit:::psp_rs_next_stream(s, TRUE),
                     parallel::nextRNGSubStream(c(10407L, s))[-1])
    s <- parallel::nextRNGStream(c(10407L, s))[-1]
  }
})

test_that("pnorm, dnorm and qnorm match R", {
  fun <- function(f, x, lower = TRUE, log = FALSE) {
    panelspatialprobit:::psp_rs_norm_fun(f, x, lower, log)
  }
  x <- c(-Inf, -1e200, -40, -38.5, -20, seq(-9, 9, by = 0.0137), 20, 38.5, 40, 1e200, Inf,
         0, 1e-300, -1e-17, NaN)
  for (lower in c(TRUE, FALSE)) for (lg in c(FALSE, TRUE)) {
    expect_r_identical(fun("pnorm", x, lower, lg), pnorm(x, lower.tail = lower, log.p = lg))
  }
  for (lg in c(FALSE, TRUE)) expect_r_identical(fun("dnorm", x, log = lg), dnorm(x, log = lg))
  p <- c(0, 1e-320, 1e-300, 1e-100, 1e-20, seq(1e-6, 1 - 1e-6, length.out = 2001),
         1 - 1e-12, 1, NaN)
  expect_r_identical(fun("qnorm", p), qnorm(p))
})

test_that("errors from Rust reach R as errors", {
  s <- 1:6
  expect_error(rs_rng("unif", 1, numeric(0), s[1:5]), "6 elements")
  expect_error(rs_rng("unif", -1, numeric(0), s), "non-negative")
  expect_error(rs_rng("none", 1, numeric(0), s), "unknown `kind`")
  expect_error(rs_rng("gamma", 1, 1, s), "2 element")
  expect_error(rs_rng("unif", 1, numeric(0), c(0L, 0L, 0L, 1L, 1L, 1L)),
               "invalid L'Ecuyer-CMRG seed")
})
