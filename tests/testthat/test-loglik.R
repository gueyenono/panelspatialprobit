test_that("GHK is exact without spatial dependence or random effects", {
  s <- small_data(seed = 7, rho = 0, sigma_alpha = 0)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "none",
                                   effects = "pooled", impacts = FALSE,
                                   control = fast_control()))
  par <- fit$coefficients
  xb <- as.vector(fit$X %*% par[colnames(fit$X)])
  exact <- sum(pnorm(ifelse(fit$y == 1, xb, -xb), log.p = TRUE))
  expect_equal(as.numeric(logLik(fit, reps = 5)), exact, tolerance = 1e-10)
  expect_equal(attr(logLik(fit, reps = 5), "df"), 3)
  s2 <- small_data(seed = 7, rho = 0, sigma_alpha = 0, cutpoints = 0.8)
  f2 <- suppressWarnings(psoprobit(y ~ x1 + x2, s2$d, c("id", "time"), s2$W, model = "none",
                                   effects = "pooled", impacts = FALSE,
                                   control = fast_control()))
  par <- f2$coefficients
  xb <- as.vector(f2$X %*% par[colnames(f2$X)])
  kap <- c(-Inf, 0, par[["cut[2|3]"]], Inf)
  exact <- sum(log(pnorm(kap[f2$y + 1] - xb) - pnorm(kap[f2$y] - xb)))
  expect_equal(as.numeric(logLik(f2, reps = 5)), exact, tolerance = 1e-10)
})

test_that("the marginal likelihood integrates out the random effects", {
  s <- small_data(seed = 8, rho = 0, sigma_alpha = 0.8, T = 4)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "none",
                                   impacts = FALSE, control = fast_control()))
  par <- fit$coefficients
  xb <- matrix(fit$X %*% par[colnames(fit$X)], fit$N)
  sgn <- matrix(2 * fit$y - 1, fit$N)
  sa <- par[["sigma_alpha"]]
  exact <- sum(vapply(seq_len(fit$N), function(i) {
    log(integrate(function(a) {
      vapply(a, function(ai) prod(pnorm(sgn[i, ] * (xb[i, ] + ai))), 0) * dnorm(a, 0, sa)
    }, -Inf, Inf)$value)
  }, 0))
  ll <- logLik(fit, reps = 4000)
  expect_lt(abs(as.numeric(ll) - exact), 4 * attr(ll, "mcse") + 0.01)
  # conditional on alpha: exact given the stored draws
  llc <- psp_loglik(fit, draws = 2, reps = 3, type = "conditional")
  d <- panelspatialprobit:::.psp_draw_subset_from(seq_along(fit$alpha_index), 2)
  pm <- panelspatialprobit:::.psp_par_matrix(fit)[fit$alpha_index[d], ]
  ex <- vapply(1:2, function(j) {
    m <- matrix(fit$X %*% pm[j, colnames(fit$X)], fit$N) + fit$alpha_draws[d[j], ]
    sum(pnorm(sgn * m, log.p = TRUE))
  }, 0)
  expect_equal(as.numeric(llc), ex, tolerance = 1e-10)
})

test_that("GHK agrees with mvtnorm for a spatial model", {
  skip_if_not_installed("mvtnorm")
  s <- small_data(seed = 9, rho = 0.6, sigma_alpha = 0, T = 1)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, effects = "pooled",
                                   impacts = FALSE, control = fast_control()))
  par <- fit$coefficients
  Si <- solve(diag(16) - par[["rho"]] * as.matrix(fit$W))
  m <- as.vector(Si %*% fit$X %*% par[colnames(fit$X)])
  set.seed(1)
  ex <- log(mvtnorm::pmvnorm(lower = ifelse(fit$y == 1, 0, -Inf),
                             upper = ifelse(fit$y == 1, Inf, 0), mean = m,
                             sigma = tcrossprod(Si),
                             algorithm = mvtnorm::GenzBretz(maxpts = 2e5, abseps = 1e-7)))
  ll <- logLik(fit, reps = 4000)
  expect_lt(abs(as.numeric(ll) - as.numeric(ex)), 4 * attr(ll, "mcse") + 0.01)
  dic <- psp_dic(fit, draws = 10, reps = 100)
  expect_named(dic, c("dic", "p_d", "d_bar", "d_hat"))
  expect_true(is.finite(dic$dic))
  expect_true(is.finite(AIC(fit)))
})
