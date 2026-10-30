# Paper Table 5: time per iteration and ESS per second for SAR binary, SAR
# ordered (J = 4) and SEM binary random-effects models, k = 5, at N = 1,000
# (25 x 40 rook grid, T = 10, 5,000 + 5,000 iterations) and N = 10,000
# (100 x 100 rook grid, T = 20, 1,000 + 1,000 iterations). One chain on one
# core; times include warm-up and the log-determinant grid, and ESS per second
# divides the bulk ESS by the whole run time (as in bench_core.R).
# Writes validation/bench/models_<date>.txt (report) and .csv (per parameter).
# Usage: Rscript validation/bench/bench_models.R [lib]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args) >= 1) args[1])
sizes <- list(list(N = 1000, grid = c(25, 40), T = 10, iter = 5000),
              list(N = 10000, grid = c(100, 100), T = 20, iter = 1000))
beta <- c(0, 1, -1, 0.5, 0.5)
f <- y ~ x1 + x2 + x3 + x4
rows <- list()
for (sz in sizes) {
  W <- psp_weights(rook_weights(sz$grid[1], sz$grid[2]), style = "W")
  ctl <- psp_control(iter = sz$iter, warmup = sz$iter, chains = 1, seed = 1, progress = FALSE)
  designs <- list(
    sar_binary = list(fun = psprobit, model = "sar",
                      d = sim_psprobit(W, T = sz$T, beta = beta, rho = 0.5, sigma_alpha = 0.7,
                                       seed = 1)),
    sar_ordered = list(fun = psoprobit, model = "sar",
                       d = sim_psprobit(W, T = sz$T, beta = beta, rho = 0.5, sigma_alpha = 0.7,
                                        cutpoints = c(0.8, 1.6), seed = 1)),
    sem_binary = list(fun = psprobit, model = "sem",
                      d = sim_psprobit(W, T = sz$T, beta = beta, rho = 0, lambda = 0.5,
                                       sigma_alpha = 0.7, seed = 1)))
  for (dn in names(designs)) {
    ds <- designs[[dn]]
    t0 <- proc.time()
    fit <- suppressWarnings(ds$fun(f, ds$d, c("id", "time"), W, model = ds$model,
                                   impacts = FALSE, control = ctl))
    el <- (proc.time() - t0)[["elapsed"]]
    s <- fit$diagnostics
    rows[[length(rows) + 1]] <- data.frame(
      design = dn, N = sz$N, T = sz$T, n = sz$N * sz$T, iterations = 2 * sz$iter,
      seconds = round(el, 1), ms_per_iter = round(1000 * el / (2 * sz$iter), 2),
      variable = s$variable, ess_bulk = round(s$ess_bulk), ess_per_sec = round(s$ess_bulk / el, 2))
    cat(dn, sz$N, round(el, 1), "s\n")
  }
}
res <- do.call(rbind, rows)
# per-parameter results for validation/vignette_data.R
utils::write.csv(res, sprintf("validation/bench/models_%s.csv", format(Sys.Date())),
                 row.names = FALSE)
out <- sprintf("validation/bench/models_%s.txt", format(Sys.Date()))
sink(out)
cat(sprintf("panelspatialprobit %s, R %s; one chain on one core, times include warm-up\n\n",
            packageVersion("panelspatialprobit"), getRversion()))
for (key in unique(paste(res$design, res$N))) {
  r <- res[paste(res$design, res$N) == key, ]
  cat(sprintf("[%s] N = %d, T = %d (n = %d): %d iterations in %.1f s (%.2f ms/iter); ESS/s %.2f-%.2f\n",
              r$design[1], r$N[1], r$T[1], r$n[1], r$iterations[1], r$seconds[1],
              r$ms_per_iter[1], min(r$ess_per_sec), max(r$ess_per_sec)))
  print(r[, c("variable", "ess_bulk", "ess_per_sec")], row.names = FALSE)
  cat("\n")
}
sink()
cat(readLines(out), sep = "\n")
