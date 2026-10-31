# Cross-build check of the Rust sampler against the frozen C++ build
# (plan/06 steps 15-24, §6.10 layer 3): each case is fitted by both libraries
# in separate R processes with the same seed. Exact cases must give the same
# draws to 1e-8 (relative); the others are compared by posterior means within
# 4 Monte Carlo standard errors (needs longer chains: iter >= 2000).
# Usage: Rscript validation/rust_port/compare_sampler.R <lib_cpp> <lib_rust>
#          [iter] [warmup] [chains] [case ...]
a <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
iter <- if (length(a) >= 3) a[3] else "40"
warmup <- if (length(a) >= 4) a[4] else "20"
chains <- if (length(a) >= 5) a[5] else "2"
which <- if (length(a) > 5) a[-(1:5)] else character(0)
run <- function(lib) {
  out <- tempfile(fileext = ".rds")
  st <- system2("Rscript", c(file.path(here, "fit_cases.R"), lib, out, iter, warmup, chains,
                             which))
  if (st != 0) stop("fitting failed with ", lib)
  readRDS(out)
}
cpp <- run(a[1])
rs <- run(a[2])
source(file.path(here, "cases.R"))
library(panelspatialprobit, lib.loc = a[2])
exact <- vapply(psp_cases(), `[[`, logical(1), "exact")
rel <- function(x, y) max(abs(x - y) / pmax(1, abs(y)), 0)
rows <- lapply(names(cpp), function(nm) {
  x <- rs[[nm]]
  y <- cpp[[nm]]
  if (!is.null(x$error) || !is.null(y$error)) {
    return(data.frame(case = nm, exact = exact[[nm]], check = "error", value = NA, ok = FALSE,
                      cpp_s = NA, rust_s = NA, note = paste(x$error, y$error)))
  }
  if (exact[[nm]]) {
    v <- max(rel(unclass(x$draws), unclass(y$draws)), rel(x$z_mean, y$z_mean))
    data.frame(case = nm, exact = TRUE, check = "max rel diff", value = signif(v, 3),
               ok = v < 1e-8, cpp_s = y$seconds, rust_s = x$seconds, note = "")
  } else {
    mx <- posterior::as_draws_array(x$draws)
    my <- posterior::as_draws_array(y$draws)
    m <- posterior::summarise_draws(mx, "mean", "mcse_mean")
    n <- posterior::summarise_draws(my, "mean", "mcse_mean")
    z <- abs(m$mean - n$mean) / sqrt(m$mcse_mean^2 + n$mcse_mean^2)
    data.frame(case = nm, exact = FALSE, check = "max |z| of means", value = round(max(z), 2),
               ok = max(z) < 4, cpp_s = y$seconds, rust_s = x$seconds, note = "")
  }
})
res <- do.call(rbind, rows)
print(res, row.names = FALSE)
cat(if (all(res$ok)) "ALL OK\n" else "FAILURES\n")
