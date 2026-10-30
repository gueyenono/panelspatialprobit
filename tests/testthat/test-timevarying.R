test_that("time-varying weights are aligned by period name and unit", {
  s <- tv_data(seed = 2, T = 3)
  Wn <- stats::setNames(rev(s$W), c("3", "2", "1"))
  perm <- c(16:1)
  Wn <- lapply(Wn, function(w) w[perm, perm])
  des <- suppressWarnings(panelspatialprobit:::.psp_design(
    y ~ x1, s$d, c("id", "time"), Wn, "sdm", NULL, NULL, FALSE, "binary"))
  expect_true(is.list(des$W))
  expect_equal(names(des$W), c("1", "2", "3"))
  for (t in 1:3) {
    expect_equal(as.matrix(des$W[[t]])[rownames(s$W[[t]]), rownames(s$W[[t]])],
                 as.matrix(s$W[[t]]))
    rows <- (t - 1) * 16 + 1:16
    expect_equal(unname(des$X[rows, "W_x1"]),
                 unname(as.vector(des$W[[t]] %*% des$X[rows, "x1"])))
  }
  # identical matrices collapse to a single one
  same <- rep(list(s$W[[1]]), 3)
  des <- suppressWarnings(panelspatialprobit:::.psp_design(
    y ~ x1, s$d, c("id", "time"), same, "sar", NULL, NULL, FALSE, "binary"))
  expect_false(is.list(des$W))
  expect_error(psprobit(y ~ x1, s$d, c("id", "time"), s$W[1:2]), "list of 2")
  expect_error(psprobit(y ~ x1, s$d, c("id", "time"), stats::setNames(s$W, c("a", "b", "c"))),
               "periods")
})

test_that("time-varying fits give the same results as a list of identical matrices", {
  s <- small_data(seed = 3)
  f1 <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                  control = fast_control()))
  f2 <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), rep(list(s$W), 3),
                                  impacts = FALSE, control = fast_control()))
  expect_equal(unclass(f1$draws), unclass(f2$draws))
})

test_that("the log-determinant grid averages the periods", {
  s <- tv_data(seed = 4, T = 3)
  avg <- panelspatialprobit:::.psp_lndet_avg(s$W, -1, 1, "eigen")
  one <- lapply(s$W, panelspatialprobit:::.psp_lndet, rmin = -1, rmax = 1, method = "eigen")
  expect_equal(avg$lndet, Reduce(`+`, lapply(one, `[[`, "lndet")) / 3)
})

test_that("models with time-varying weights fit, predict and check", {
  s <- tv_data(seed = 5, T = 3, cutpoints = 0.8)
  fit <- suppressWarnings(psoprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "sarar",
                                    control = fast_control()))
  expect_s3_class(fit$impacts, "psp_impacts")
  p <- predict(fit, draws = 20)
  expect_equal(dim(p), c(48L, 3L))
  expect_equal(unname(rowSums(p)), rep(1, 48))
  pp <- psp_ppc(fit, draws = 10, seed = 1)
  expect_length(pp$moran$observed, 3)
  d_dyn <- sim_psprobit(tv_weights(5, 6), T = 5, beta = c(-0.3, 1, -1), rho = 0.3,
                        gamma = 0.8, seed = 6)
  fit_dyn <- suppressWarnings(psprobit(y ~ x1 + x2, d_dyn, c("id", "time"), tv_weights(5, 6),
                                       dynamic = psp_dynamic(lag_Wy = TRUE), impacts = FALSE,
                                       control = fast_control()))
  expect_length(fit_dyn$W, 4)
  expect_equal(names(fit_dyn$W), as.character(2:5))
})

test_that("spatial random effects need time-invariant weights", {
  s <- tv_data(seed = 2, T = 3)
  expect_error(psprobit(y ~ x1, s$d, c("id", "time"), s$W, effects = "spatial_random"),
               "time-invariant")
})
