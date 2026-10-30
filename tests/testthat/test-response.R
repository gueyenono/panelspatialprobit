test_that("the first response equals the unit-specific impact", {
  s <- small_data(seed = 6, T = 3, cutpoints = 0.8)
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                    control = fast_control()))
  r <- psp_response_path(fit, "x1", shock = 1e-5, draws = 5, nsim = 2, seed = 1)
  imp <- summary(psp_impacts(fit, draws = 5, type = "unit", by_period = TRUE))
  imp <- imp[imp$variable == "x1" & imp$effect == "total" & imp$period == "1", ]
  r0 <- r$summary[r$summary$horizon == 0, ]
  expect_equal(r0$mean / 1e-5, imp$mean, tolerance = 1e-3)
  # static model: a temporary shock has no lasting effect
  expect_true(all(r$summary$mean[r$summary$horizon > 0] == 0))
  expect_equal(unique(r$summary$units), "all")
})

test_that("dynamic responses persist and are reproducible", {
  W <- psp_weights(rook_weights(5, 5), style = "W")
  d <- sim_psprobit(W, T = 6, beta = c(-0.3, 1.5, -1), rho = 0.2, sigma_alpha = 0.3,
                    gamma = 1.5, seed = 4)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), W,
                                   dynamic = psp_dynamic(lag_Wy = TRUE), impacts = FALSE,
                                   control = fast_control()))
  r <- psp_response_path(fit, "x1", draws = 30, nsim = 10, seed = 2)
  expect_equal(r$summary$horizon, 0:4)
  expect_gt(r$summary$mean[1], 0)
  expect_gt(r$summary$mean[2], 0)
  expect_identical(r, psp_response_path(fit, "x1", draws = 30, nsim = 10, seed = 2))
  rp <- psp_response_path(fit, "x1", units = c("1", "7"), start = 4, type = "permanent",
                          draws = 10, nsim = 5, seed = 3)
  expect_setequal(unique(rp$summary$units), c("all", "shocked", "others"))
  expect_equal(unique(rp$summary$period), c("4", "5", "6"))
  expect_error(psp_response_path(fit, "lag_y"), "variable")
  expect_error(psp_response_path(fit, "x1", start = 1), "modelled periods")
  expect_output(print(r), "temporary shock")
})
