# Roadmap 6.7: large benchmark, SAR random-effects binary, N = 10,000
# (100 x 100 rook grid), T = 20, k = 5 (200,000 observations). Runs a shorter
# chain and reports the time per iteration; the plan/02 s2.13 target is
# 10,000 iterations in about 30 minutes.
# Usage: Rscript validation/bench/bench_large.R [lib] [iterations]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args) >= 1) args[1])
iters <- if (length(args) >= 2) as.integer(args[2]) else 1000L
W <- psp_weights(rook_weights(100, 100), style = "W")
d <- sim_psprobit(W, T = 20, beta = c(0, 1, -1, 0.5, 0.5), rho = 0.5, sigma_alpha = 0.7,
                  seed = 1)
t0 <- proc.time()
fit <- suppressWarnings(psprobit(
  y ~ x1 + x2 + x3 + x4, d, c("id", "time"), W, impacts = FALSE,
  control = psp_control(iter = iters, warmup = iters, chains = 1, seed = 1, progress = FALSE)))
el <- (proc.time() - t0)[["elapsed"]]
t1 <- proc.time()
imp <- psp_impacts(fit, draws = 50)
el_imp <- (proc.time() - t1)[["elapsed"]]
per_iter <- el / (2 * iters)
out <- sprintf("validation/bench/large_%s.txt", format(Sys.Date()))
sink(out)
cat(sprintf("N = 10000, T = 20, k = 5 (n = 200000): %d iterations in %.0f s (%.1f ms/iter,",
            2L * iters, el, 1000 * per_iter),
    sprintf("including the log-determinant grid, %s);", fit$lndet$method),
    sprintf("10,000 iterations would take about %.1f min.", 10000 * per_iter / 60), "\n")
cat(sprintf("Impacts (50 draws, stochastic diagonals): %.0f s\n\n", el_imp))
s <- fit$diagnostics
print(data.frame(variable = s$variable, mean = round(colMeans(posterior::as_draws_matrix(fit$draws)), 3),
                 truth = c(0, 1, -1, 0.5, 0.5, 0.5, 0.7), ess_bulk = round(s$ess_bulk)),
      row.names = FALSE)
sink()
cat(readLines(out), sep = "\n")
