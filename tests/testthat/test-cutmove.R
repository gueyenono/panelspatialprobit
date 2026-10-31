# Stationarity of the cut-point updates on a toy ordered model (plan/04 s4.1).
# If (kappa, z) is drawn from the joint prior p(kappa) prod N(z_i; mu_i, 1) and
# y is the implied category vector, then (kappa, z) | y follows the posterior.
# A kernel that leaves every posterior p(kappa, z | y) invariant must therefore
# leave the joint distribution, and so the prior marginal of kappa, unchanged.

# The Rust cut-point kernel (pspcore::sampler::cut_kernel), seeded from R's
# RNG so that set.seed() fixes the whole test.
cut_kernel <- function(y, mu, z, kap, reps, ac, move) {
  panelspatialprobit:::psp_rs_cut_kernel(y, mu, z, kap, 4L, 1, 1, 0.7, as.integer(reps), ac,
                                         move, sample.int(.Machine$integer.max, 6))
}

draw_prior_kappa <- function(kmean, ksd) {
  repeat {
    k <- rnorm(2, kmean, ksd)
    if (k[1] > 0 && k[2] > k[1]) return(k)
  }
}

run_kernel <- function(M, ac, move, reps = 3, n = 12) {
  mu <- seq(-1, 2.5, length.out = n)
  out <- matrix(NA_real_, M, 2)
  for (m in seq_len(M)) {
    kap <- draw_prior_kappa(1, 1)
    z <- rnorm(n, mu)
    y <- findInterval(z, c(0, kap)) + 1L
    res <- cut_kernel(y, mu, z, kap, reps, ac, move)
    y_new <- findInterval(res$z, c(0, res$kappa)) + 1L
    if (!identical(y_new, y)) stop("the kernel changed a category")
    out[m, ] <- res$kappa
  }
  out
}

test_that("cut-point kernels leave the joint distribution invariant", {
  set.seed(11)
  M <- 3000
  ref <- t(replicate(M, draw_prior_kappa(1, 1)))
  for (cfg in list(c(ac = FALSE, move = TRUE), c(ac = TRUE, move = FALSE),
                   c(ac = TRUE, move = TRUE))) {
    out <- run_kernel(M, cfg[["ac"]], cfg[["move"]])
    for (j in 1:2) {
      p <- suppressWarnings(ks.test(out[, j], ref[, j])$p.value)
      expect_gt(p, 1e-3)
    }
  }
})

test_that("the transformation move actually moves the cut-points", {
  set.seed(12)
  n <- 12
  mu <- seq(-1, 2.5, length.out = n)
  kap <- c(0.8, 1.6)
  z <- rnorm(n, mu)
  y <- findInterval(z, c(0, kap)) + 1L
  res <- cut_kernel(y, mu, z, kap, 20, FALSE, TRUE)
  expect_false(isTRUE(all.equal(res$kappa, kap)))
})
