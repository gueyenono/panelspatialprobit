test_that("psp_weights coerces, validates and row-standardises", {
  W <- rook_weights(3, 3)
  expect_s4_class(W, "dgCMatrix")
  expect_equal(rownames(W), as.character(1:9))
  Ws <- psp_weights(W, style = "W")
  expect_equal(unname(Matrix::rowSums(Ws)), rep(1, 9))
  expect_equal(rownames(Ws), as.character(1:9))
  Wd <- psp_weights(as.matrix(W), ids = letters[1:9])
  expect_equal(rownames(Wd), letters[1:9])
  bad <- as.matrix(W)
  diag(bad) <- 1
  expect_error(psp_weights(bad), "zero diagonal")
  expect_error(psp_weights(matrix(0, 2, 3)), "square")
  expect_error(psp_weights(W, ids = 1:3), "one entry per row")
})

test_that("isolates keep zero rows", {
  W <- as.matrix(rook_weights(2, 2))
  W[1, ] <- 0
  W[, 1] <- 0
  Ws <- psp_weights(W, style = "W")
  expect_equal(unname(Matrix::rowSums(Ws)), c(0, 1, 1, 1))
})

test_that("log-determinant grids match determinant()", {
  W <- psp_weights(rook_weights(4, 5), style = "W")
  ld_e <- panelspatialprobit:::.psp_lndet(W, -1, 1, "eigen")
  ld_s <- panelspatialprobit:::.psp_lndet(W, -1, 1, "sparse")
  for (r in c(-0.9, -0.3, 0.2, 0.7, 0.95)) {
    i <- which.min(abs(ld_e$grid - r))
    exact_i <- as.numeric(determinant(diag(20) - ld_e$grid[i] * as.matrix(W))$modulus)
    expect_equal(ld_e$lndet[i], exact_i, tolerance = 1e-8)
    expect_equal(ld_s$lndet[i], exact_i, tolerance = 1e-4)
  }
})

test_that("log-determinant grid agrees with spatialreg", {
  skip_if_not_installed("spatialreg")
  skip_if_not_installed("spdep")
  W <- psp_weights(rook_weights(5, 6), style = "W")
  nb <- spdep::mat2listw(as.matrix(rook_weights(5, 6)), style = "B")$neighbours
  lw <- spdep::nb2listw(nb, style = "W")
  ld <- panelspatialprobit:::.psp_lndet(W, -1, 1, "eigen")
  env <- new.env()
  assign("listw", lw, envir = env)
  assign("can.sim", TRUE, envir = env)
  assign("similar", FALSE, envir = env)
  assign("family", "SAR", envir = env)
  assign("n", 30, envir = env)
  assign("verbose", FALSE, envir = env)
  assign("method", "eigen", envir = env)
  spatialreg::eigen_setup(env)  # modifies env in place
  for (r in c(-0.8, -0.2, 0.3, 0.9)) {
    i <- which.min(abs(ld$grid - r))
    expect_equal(ld$lndet[i], spatialreg::do_ldet(ld$grid[i], env), tolerance = 1e-8)
  }
})
