test_that("predict(newdata) reproduces in-sample predictions and responds to X", {
  s <- small_data()
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, model = "sdm",
                                   cre = ~ x1, impacts = FALSE, control = fast_control()))
  nd <- s$d[sample(nrow(s$d)), ]
  nd$y <- NULL
  idx <- match(paste(nd$id, nd$time), paste(s$d$id, s$d$time))
  expect_equal(predict(fit, newdata = nd, draws = 10), predict(fit, draws = 10)[idx])
  nd2 <- s$d
  nd2$x1 <- nd2$x1 + 1  # beta_x1 > 0 in the simulation
  up <- predict(fit, newdata = nd2, draws = 10)
  expect_gt(mean(up), mean(predict(fit, draws = 10)))
  # CRE columns stay at their fitted values
  expect_equal(mean(predict(fit, newdata = s$d, draws = 10)), mean(predict(fit, draws = 10)))
  expect_error(predict(fit, newdata = s$d[s$d$time < 3, ]), "periods")

  so <- small_data(cutpoints = c(0.7, 1.5))
  fo <- suppressWarnings(psoprobit(y ~ x1 + x2, so$d, c("id", "time"), so$W, model = "sem",
                                   impacts = FALSE, control = fast_control()))
  ndo <- so$d
  ndo$y <- NULL
  expect_equal(predict(fo, newdata = ndo, draws = 5), predict(fo, draws = 5))
})

test_that("posterior predictive checks return calibrated statistics", {
  s <- small_data(T = 4)
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, s$d, c("id", "time"), s$W, impacts = FALSE,
                                   control = fast_control()))
  pp <- psp_ppc(fit, draws = 50, seed = 1)
  expect_s3_class(pp, "psp_ppc")
  expect_equal(names(pp$share$observed), paste0("P(y=1) t=", 1:4))
  expect_equal(dim(pp$moran$replicated), c(50, 4))
  expect_length(pp$persistence$observed, 3)
  expect_equal(unname(pp$share$observed), as.vector(tapply(s$d$y, s$d$time, mean)))
  expect_true(all(pp$share$ppp > 0 & pp$share$ppp < 1))
  expect_output(print(pp), "ppp")
  expect_identical(psp_ppc(fit, draws = 20, seed = 2), psp_ppc(fit, draws = 20, seed = 2))

  so <- small_data(cutpoints = c(0.7, 1.5))
  so$d$y[c(2, 9)] <- NA
  fo <- suppressWarnings(psoprobit(y ~ x1 + x2, so$d, c("id", "time"), so$W, effects = "pooled",
                                   model = "sarar", impacts = FALSE, control = fast_control()))
  ppo <- psp_ppc(fo, draws = 20, seed = 1)
  expect_length(ppo$share$observed, 4 * 3)
})
