test_that("t-error probabilities match the t distribution without spatial terms", {
  s <- small_data(seed = 20, rho = 0, sigma_alpha = 0)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "none",
                                   effects = "pooled", errors = psp_t(5), impacts = FALSE,
                                   control = fast_control(mixture_draws = 20000)))
  expect_equal(fit$errors$df, 5)
  expect_length(fit$v_mean, 48)
  par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
  fit$draws <- posterior::subset_draws(fit$draws, chain = 1, iteration = 1)
  p <- predict(fit)
  xb <- as.vector(fit$X %*% par[colnames(fit$X)])
  exact <- pt(xb, 5)
  exact[fit$row_order] <- exact
  expect_equal(p, exact, tolerance = 0.01)
})

test_that("t-error impacts are derivatives of the predicted probabilities", {
  s <- small_data(seed = 21)
  for (mod in c("sdm", "sarar")) {
    fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = mod,
                                     errors = psp_t(4), impacts = FALSE,
                                     control = fast_control(mixture_draws = 10)))
    fit$draws <- posterior::subset_draws(fit$draws, chain = 1, iteration = 1)
    imp <- psp_impacts(fit)
    h <- 1e-5
    shift <- function(eps) {
      f <- fit
      f$X[, "x1"] <- f$X[, "x1"] + eps
      if ("W_x1" %in% colnames(f$X)) {
        f$X[, "W_x1"] <- f$X[, "W_x1"] + eps * rep(Matrix::rowSums(f$W), f$T)
      }
      mean(predict(f))
    }
    expect_equal(imp$draws[1, "x1", 1, "total"], (shift(h) - shift(-h)) / (2 * h),
                 tolerance = 1e-5)
  }
  expect_error(logLik(fit), "normal errors")
  expect_error(psprobit(y ~ x1, s$d, c("id", "time"), s$W, errors = "cauchy"), "psp_t")
  expect_error(psp_t(2), "larger than 2")
})

test_that("the stochastic reduced form approximates the exact one for t errors", {
  s <- small_data(seed = 22)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "sarar",
                                   errors = psp_t(4), impacts = FALSE,
                                   control = fast_control(mixture_draws = 200)))
  par <- panelspatialprobit:::.psp_par_matrix(fit)[1, ]
  mix <- panelspatialprobit:::.psp_t_mix(fit)
  ex <- panelspatialprobit:::.psp_reduced_form(fit, par, NULL, mix = mix)
  fit$control$impacts_exact_max <- 1
  fit$control$impacts_probes <- 3000L
  st <- panelspatialprobit:::.psp_reduced_form(fit, par, panelspatialprobit:::.psp_probes(fit),
                                               mix = mix)
  p_ex <- rowMeans(pnorm(ex$m / ex$s))
  p_st <- rowMeans(pnorm(st$m / st$s))
  expect_lt(max(abs(p_ex - p_st)), 0.02)
})
