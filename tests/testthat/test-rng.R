# The Rust ports of R's generators reproduce R's draws exactly (plan/06 §6.6).

with_lecuyer <- function(seed, code) {
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old_seed, envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion")
  set.seed(seed)
  code(.Random.seed[-1])
}

test_that("uniform draws match runif() under L'Ecuyer-CMRG", {
  for (seed in c(1L, 42L, 2026L)) {
    with_lecuyer(seed, function(s) {
      expect_identical(panelspatialprobit:::psp_rs_rng("unif", 1e5L, s), runif(1e5))
    })
  }
  s <- parallel::nextRNGStream(c(10407L, 1L, 2L, 3L, 4L, 5L, 6L))[-1]
  expect_identical(panelspatialprobit:::psp_rs_rng("unif", 0L, s), numeric(0))
})

test_that("errors from Rust reach R as errors", {
  s <- c(1L, 2L, 3L, 4L, 5L, 6L)
  expect_error(panelspatialprobit:::psp_rs_rng("unif", 1L, s[1:5]), "6 elements")
  expect_error(panelspatialprobit:::psp_rs_rng("unif", -1L, s), "non-negative")
  expect_error(panelspatialprobit:::psp_rs_rng("none", 1L, s), "unknown `kind`")
  expect_error(panelspatialprobit:::psp_rs_rng("unif", 1L, c(0L, 0L, 0L, 1L, 1L, 1L)),
               "invalid L'Ecuyer-CMRG seed")
})
