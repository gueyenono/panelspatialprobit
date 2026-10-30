tn_mean <- function(mu, s, a, b) {
  al <- (a - mu) / s
  be <- (b - mu) / s
  mu + s * (dnorm(al) - dnorm(be)) / (pnorm(be) - pnorm(al))
}

test_that("truncated normal draws respect bounds and match analytic means", {
  set.seed(1)
  cases <- list(c(0, 1, 0, Inf), c(0, 1, -Inf, -1), c(0, 1, -0.5, 0.7), c(1, 2, 3, 3.2),
                c(0, 1, 6, Inf), c(0, 1, -Inf, -8), c(0, 1, 2, 9), c(-3, 0.5, -10, -2.9),
                c(0, 1, -5, 5))
  for (cs in cases) {
    x <- panelspatialprobit:::psp_rtnorm_cpp(20000, cs[1], cs[2], cs[3], cs[4])
    expect_true(all(x >= cs[3] & x <= cs[4]))
    expect_true(all(is.finite(x)))
    expect_lt(abs(mean(x) - tn_mean(cs[1], cs[2], cs[3], cs[4])), 0.02 * cs[2])
  }
})

test_that("extreme tails are finite", {
  x <- panelspatialprobit:::psp_rtnorm_cpp(1000, 0, 1, 40, Inf)
  expect_true(all(is.finite(x) & x >= 40))
  x <- panelspatialprobit:::psp_rtnorm_cpp(1000, 50, 1, -Inf, 0)
  expect_true(all(is.finite(x) & x <= 0))
})
