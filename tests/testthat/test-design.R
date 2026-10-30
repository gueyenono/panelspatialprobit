test_that("row order of the data does not change the internal design", {
  s <- small_data()
  d <- s$d
  perm <- sample(nrow(d))
  a <- panelspatialprobit:::.psp_design(y ~ x1 + x2, d, c("id", "time"), s$W, "sar",
                                        NULL, NULL, FALSE, "binary")
  b <- panelspatialprobit:::.psp_design(y ~ x1 + x2, d[perm, ], c("id", "time"), s$W, "sar",
                                        NULL, NULL, FALSE, "binary")
  expect_identical(a$X, b$X)
  expect_identical(a$y, b$y)
  expect_identical(perm[b$row_order], a$row_order)
})

test_that("design adds Durbin, Mundlak and time columns", {
  s <- small_data()
  des <- panelspatialprobit:::.psp_design(y ~ x1 + x2, s$d, c("id", "time"), s$W, "sdm",
                                          ~ x1, ~ x2, TRUE, "binary")
  expect_true(all(c("W_x2", "mean_x1", "time_2", "time_3") %in% colnames(des$X)))
  expect_equal(unname(des$roles[c("W_x2", "mean_x1", "time_2")]), c("durbin", "cre", "time"))
  # W_x2 is W times x2 within each period
  rows <- 1:16
  expect_equal(unname(des$X[rows, "W_x2"]), as.vector(s$W %*% des$X[rows, "x2"]))
  # Mundlak mean is constant within units
  m <- matrix(des$X[, "mean_x1"], 16)
  expect_equal(m[, 1], m[, 3])
})

test_that("informative errors for bad panels", {
  s <- small_data()
  d <- s$d
  expect_error(psprobit(y ~ x1, d[-1, ], c("id", "time"), s$W), "balanced")
  d2 <- d
  d2$time[2] <- d2$time[1]
  d2$id[2] <- d2$id[1]
  expect_error(psprobit(y ~ x1, d2, c("id", "time"), s$W), "Duplicated")
  d3 <- d
  d3$id[d3$id == "1"] <- "zz"
  expect_error(psprobit(y ~ x1, d3, c("id", "time"), s$W), "not found")
  expect_error(psprobit(y ~ x1 - 1, d, c("id", "time"), s$W), "intercept")
  expect_error(psprobit(y ~ x1, d, c("id", "time"), s$W, model = "sem",
                        M = rook_weights(3, 3)), "rows")
  expect_error(psoprobit(y ~ x1, d, c("id", "time"), s$W), "at least 3")
  d4 <- d
  d4$y <- factor(d4$y + 1, levels = 1:3)
  expect_error(psoprobit(y ~ x1, d4, c("id", "time"), s$W), "Empty")
})

test_that("Chamberlain CRE adds one column per regressor and period", {
  s <- small_data(seed = 2)
  des <- panelspatialprobit:::.psp_design(y ~ x1 + x2, s$d, c("id", "time"), s$W, "sar",
                                          list(type = "chamberlain", vars = ~ x1), NULL,
                                          FALSE, "binary")
  cc <- paste0("cre_x1_", 1:3)
  expect_true(all(cc %in% colnames(des$X)))
  expect_true(all(des$roles[cc] == "cre"))
  # column for period 2 holds each unit's period-2 value in every period
  x1_t2 <- des$X[17:32, "x1"]  # rows are sorted time-slow, unit-fast
  expect_equal(unname(des$X[1:16, "cre_x1_2"]), unname(x1_t2))
  expect_equal(unname(des$X[33:48, "cre_x1_2"]), unname(x1_t2))
  mund <- panelspatialprobit:::.psp_design(y ~ x1 + x2, s$d, c("id", "time"), s$W, "sar",
                                           list(type = "mundlak"), NULL, FALSE, "binary")
  expect_true(all(c("mean_x1", "mean_x2") %in% colnames(mund$X)))
})

test_that("error weights are aligned with W and validated", {
  s <- small_data()
  M <- s$W[16:1, 16:1]
  des <- panelspatialprobit:::.psp_design(y ~ x1, s$d, c("id", "time"), s$W, "sem", NULL,
                                          NULL, FALSE, "binary", M)
  expect_equal(rownames(des$M), rownames(s$W))
  expect_equal(as.matrix(des$M), as.matrix(s$W))
  expect_null(suppressWarnings(panelspatialprobit:::.psp_design(
    y ~ x1, s$d, c("id", "time"), s$W, "sar", NULL, NULL, FALSE, "binary"))$M)
  expect_warning(panelspatialprobit:::.psp_design(y ~ x1, s$d, c("id", "time"), s$W, "sar",
                                                  NULL, NULL, FALSE, "binary", M), "ignored")
  expect_warning(panelspatialprobit:::.psp_design(y ~ x1, s$d, c("id", "time"), s$W, "sarar",
                                                  NULL, NULL, FALSE, "binary"), "weakly")
  M2 <- M
  dimnames(M2) <- list(letters[1:16], letters[1:16])
  expect_error(panelspatialprobit:::.psp_design(y ~ x1, s$d, c("id", "time"), s$W, "sem",
                                                NULL, NULL, FALSE, "binary", M2), "match")
})
