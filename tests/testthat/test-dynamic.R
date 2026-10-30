test_that("dynamic design adds lag and initial-condition columns", {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = 4, beta = c(0, 1, -0.5), rho = 0.3, gamma = 0.7, seed = 1)
  des <- panelspatialprobit:::.psp_design(y ~ x1 + x2, d, c("id", "time"), W, "sar", NULL,
                                          NULL, FALSE, "binary", NULL,
                                          psp_dynamic(lag_Wy = TRUE))
  expect_equal(des$T, 3)
  expect_equal(nrow(des$X), 48)
  expect_true(all(c("lag_y", "W_lag_y", "init_y", "mean_x1", "mean_x2") %in% colnames(des$X)))
  # lag of period 2 (first estimation period) is the period-1 outcome
  expect_equal(unname(des$X[1:16, "lag_y"]), d$y[1:16])
  expect_equal(unname(des$X[33:48, "lag_y"]), d$y[33:48])
  expect_equal(unname(des$X[33:48, "W_lag_y"]), as.vector(as.matrix(W %*% d$y[33:48])))
  expect_equal(unname(des$X[17:32, "init_y"]), d$y[1:16])
  expect_equal(des$y, d$y[17:64])
  expect_equal(unname(des$roles[c("lag_y", "init_y")]), c("dynamic", "initial"))

  d$y <- d$y + 1L
  d$y[d$y == 2 & d$x1 > 0.5] <- 3L
  des_o <- panelspatialprobit:::.psp_design(y ~ x1, d, c("id", "time"), W, "sar", NULL, NULL,
                                            FALSE, "ordered", NULL, psp_dynamic())
  expect_true(all(c("lag_y2", "lag_y3", "init_y2", "init_y3") %in% colnames(des_o$X)))
  expect_error(panelspatialprobit:::.psp_design(y ~ x1, d[d$time <= 2, ], c("id", "time"), W,
                                                "sar", NULL, NULL, FALSE, "ordered", NULL,
                                                psp_dynamic()), "three periods")
})

test_that("dynamic fits recover state dependence and predict in data order", {
  W <- psp_weights(rook_weights(8, 8), style = "W")
  d <- sim_psprobit(W, T = 8, beta = c(-0.3, 1, -0.8), rho = 0.3, sigma_alpha = 0.5,
                    gamma = 1, seed = 2)
  perm <- sample(nrow(d))
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, d[perm, ], c("id", "time"), W,
                                   dynamic = psp_dynamic(), impacts = FALSE,
                                   control = psp_control(iter = 1000, warmup = 500, chains = 2,
                                                         seed = 1, progress = FALSE)))
  m <- posterior::as_draws_matrix(fit$draws)
  expect_lt(abs(mean(m[, "lag_y"]) - 1) / sd(m[, "lag_y"]), 3.5)
  p <- predict(fit, draws = 10)
  expect_length(p, 64 * 7)
  # aligned with the non-initial rows in data order: compare with a fit on the
  # unpermuted data (same sorted design, so the same draws)
  fit0 <- suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), W,
                                    dynamic = psp_dynamic(), impacts = FALSE,
                                    control = psp_control(iter = 1000, warmup = 500, chains = 2,
                                                          seed = 1, progress = FALSE)))
  dd <- d[perm, ]
  kept <- dd[dd$time > 1, ]
  d_kept <- d[d$time > 1, ]
  idx <- match(paste(kept$id, kept$time), paste(d_kept$id, d_kept$time))
  expect_equal(predict(fit, draws = 10), predict(fit0, draws = 10)[idx])
  cl <- psp_classification(fit, draws = 10)
  expect_equal(sum(cl$confusion), 64 * 7)
  expect_output(print(fit), "dynamic")
  expect_error(psprobit(y ~ x1, d, c("id", "time"), W, dynamic = "yes"), "psp_dynamic")
})
