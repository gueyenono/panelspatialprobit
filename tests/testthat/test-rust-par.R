# Rust worker threads (pspcore::par, pspr::run; plan/06 step 13). Interrupts
# are checked by validation/rust_port/check_interrupt.R, which needs a signal.

job <- function(tasks, threads, iters, progress = FALSE) {
  panelspatialprobit:::psp_rs_dummy_job(as.integer(tasks), as.integer(threads),
                                        as.integer(iters), 1:6, progress)
}

test_that("results do not depend on the number of threads", {
  one <- job(5, 1, 20)
  expect_length(one, 5)
  expect_identical(job(5, 2, 20), one)
  expect_identical(job(5, 8, 20), one)
  expect_identical(job(0, 2, 20), numeric(0))
})

test_that("task i uses the i-th stream after the seed", {
  s <- 1:6
  for (i in 1:2) s <- panelspatialprobit:::psp_rs_next_stream(s, FALSE)
  want <- sum(panelspatialprobit:::psp_rs_rng("norm", 3000L, numeric(0), s))
  expect_equal(job(2, 2, 3)[2], want, tolerance = 1e-12)
})

test_that("progress is printed from the main thread", {
  out <- capture.output(job(2, 2, 10, progress = TRUE))
  expect_true(all(c("Task 1: 100%", "Task 2: 100%") %in% out))
  expect_false(any(duplicated(out)))
})
