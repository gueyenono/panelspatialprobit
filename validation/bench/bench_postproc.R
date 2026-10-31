# Post-processing timings for the Rust migration (plan/06 §6.2): impacts,
# predict(), psp_ppc() and psp_loglik() on a SAR random-effects binary fit at
# N = 1,000 (25 x 40 rook grid), T = 10, k = 5, plus the log-determinant grids.
# Impacts use the stochastic diagonals (N > impacts_exact_max). Each task is
# timed under future::plan("sequential") and, with [workers] > 1, under
# multicore and multisession plans with that many workers (R-D11; the
# log-determinant grids run in the session). Writes
# validation/bench/postproc_<date>.txt.
# Usage: Rscript validation/bench/bench_postproc.R [lib] [workers]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args) >= 1) args[1])
W <- psp_weights(rook_weights(25, 40), style = "W")
d <- sim_psprobit(W, T = 10, beta = c(0, 1, -1, 0.5, 0.5), rho = 0.5, sigma_alpha = 0.7, seed = 1)
fit <- suppressWarnings(psprobit(
  y ~ x1 + x2 + x3 + x4, d, c("id", "time"), W, impacts = FALSE,
  control = psp_control(iter = 1000, warmup = 1000, chains = 1, seed = 1, progress = FALSE)))
time <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  proc.time()[["elapsed"]] - t0
}
workers <- if (length(args) >= 2) as.integer(args[2]) else 4L
tasks <- function() c(
  "psp_impacts(), 100 draws" = time(psp_impacts(fit, draws = 100)),
  "predict(), 100 draws" = time(predict(fit, draws = 100)),
  "psp_ppc(), 100 draws" = time(psp_ppc(fit, draws = 100, seed = 1)),
  "psp_response_path(), 100 draws x 20 sims" =
    time(psp_response_path(fit, "x1", draws = 100, nsim = 20, seed = 1)),
  "psp_loglik(type = \"conditional\"), 10 draws x 100 reps" =
    time(suppressWarnings(psp_loglik(fit, draws = 10, reps = 100, type = "conditional"))))
future::plan("sequential")
res <- data.frame(task = names(tasks()), sequential = round(unname(tasks()), 2))
if (workers > 1) {
  for (pl in c("multicore", "multisession")) {
    if (pl == "multicore" && !future::supportsMulticore()) next
    future::plan(pl, workers = workers)
    invisible(tasks())  # start the workers
    res[[sprintf("%s_%d", pl, workers)]] <- round(unname(tasks()), 2)
  }
  future::plan("sequential")
}
lnd <- c("log-determinant grid, sparse" = time(panelspatialprobit:::.psp_lndet(W, -1, 1, "sparse")),
         "log-determinant grid, eigen" = time(panelspatialprobit:::.psp_lndet(W, -1, 1, "eigen")))
out <- sprintf("validation/bench/postproc_%s.txt", format(Sys.Date()))
sink(out)
cat(sprintf("panelspatialprobit %s, R %s; SAR-RE binary fit, N = 1000, T = 10, k = 5; seconds\n\n",
            packageVersion("panelspatialprobit"), getRversion()))
print(res, row.names = FALSE)
cat("\n")
print(data.frame(task = names(lnd), seconds = round(unname(lnd), 2)), row.names = FALSE)
sink()
cat(readLines(out), sep = "\n")
