# Rust sparse operators (pspcore::sparse, plan/06 step 12) against Matrix.

dense <- function(a, n) {
  as.matrix(Matrix::sparseMatrix(i = a$i + 1L, p = a$p, x = a$x, dims = c(n, n)))
}
filter <- function(W, M, rho, lambda, t = 1L, T = 1L) {
  w <- panelspatialprobit:::.psp_rs_weights(W, T)
  m <- if (is.null(M)) w else panelspatialprobit:::.psp_rs_weights(M, T)
  panelspatialprobit:::psp_rs_filter(w, m, !is.null(M), rho, lambda, as.integer(t))
}

test_that("I - rho W matches Matrix", {
  W <- psp_weights(rook_weights(6, 7), style = "W")
  n <- nrow(W)
  a <- filter(W, NULL, 0.45, 0)
  want <- as.matrix(Matrix::Diagonal(n) - 0.45 * W)
  dimnames(want) <- NULL
  expect_identical(dense(a, n), want)
  expect_equal(a$colsq, colSums(want^2), tolerance = 1e-15)
})

test_that("(I - lambda M)(I - rho W) matches Matrix", {
  W <- psp_weights(rook_weights(6, 7), style = "W")
  M <- psp_weights(Matrix::t(rook_weights(7, 6)), style = "W")
  n <- nrow(W)
  for (par in list(c(0.4, -0.3), c(-0.9, 0.7), c(0, 0))) {
    a <- filter(W, M, par[1], par[2])
    want <- as.matrix((Matrix::Diagonal(n) - par[2] * M) %*% (Matrix::Diagonal(n) - par[1] * W))
    dimnames(want) <- NULL
    expect_equal(dense(a, n), want, tolerance = 1e-14)
    expect_equal(a$colsq, colSums(want^2), tolerance = 1e-14)
  }
})

test_that("time-varying weights are stored once per distinct matrix and applied per period", {
  W1 <- psp_weights(rook_weights(4, 5), style = "W")
  W2 <- psp_weights(rook_weights(5, 4), style = "W")
  Wl <- list(W1, W2, W1, W2, W2)
  w <- panelspatialprobit:::.psp_rs_weights(Wl, 5L)
  expect_length(w$mats, 2L)
  expect_identical(w$period, c(0L, 1L, 0L, 1L, 1L))
  expect_identical(panelspatialprobit:::.psp_rs_weights(W1, 3L)$period, rep(0L, 3))
  v <- sin(seq_len(20 * 5))
  want <- unlist(lapply(1:5, function(t) as.vector(Wl[[t]] %*% v[(t - 1) * 20 + 1:20])))
  expect_equal(panelspatialprobit:::psp_rs_weights_apply(w, v), want, tolerance = 1e-15)
  a <- filter(Wl, NULL, 0.3, 0, t = 2L, T = 5L)
  want <- as.matrix(Matrix::Diagonal(20) - 0.3 * W2)
  dimnames(want) <- NULL
  expect_identical(dense(a, 20), want)
})

test_that("malformed inputs are rejected with the field named", {
  W <- psp_weights(rook_weights(3, 3), style = "W")
  w <- panelspatialprobit:::.psp_rs_weights(W, 2L)
  bad <- w
  bad$mats[[1]]$i <- rev(bad$mats[[1]]$i)
  expect_error(panelspatialprobit:::psp_rs_weights_apply(bad, numeric(18)), "W\\$mats\\[\\[1\\]\\]")
  bad <- w
  bad$period <- c(0L, 1L)
  expect_error(panelspatialprobit:::psp_rs_weights_apply(bad, numeric(18)), "period map")
  bad <- w
  bad$mats[[1]]$x <- NULL
  expect_error(panelspatialprobit:::psp_rs_weights_apply(bad, numeric(18)), "missing field `x`")
  expect_error(panelspatialprobit:::psp_rs_weights_apply(w, numeric(5)), "length N T")
})
