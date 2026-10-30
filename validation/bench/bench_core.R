# Phase 1 benchmark: SAR random-effects binary, N = 1000 (25 x 40 rook grid), T = 10.
# Usage: Rscript validation/bench/bench_core.R [lib]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args)) args[1])
W <- psp_weights(rook_weights(25, 40), style = "W")
d <- sim_psprobit(W, T = 10, beta = c(0, 1, -1, 0.5, 0.5), rho = 0.5, sigma_alpha = 0.7, seed = 1)
t0 <- proc.time()
fit <- suppressWarnings(psprobit(y ~ x1 + x2 + x3 + x4, d, c("id", "time"), W, impacts = FALSE,
                                 control = psp_control(iter = 5000, warmup = 5000, chains = 1, seed = 1)))
el <- (proc.time() - t0)[["elapsed"]]
t1 <- proc.time(); imp <- psp_impacts(fit, draws = 200); el_imp <- (proc.time() - t1)[["elapsed"]]
s <- fit$diagnostics
cat(sprintf("N = 1000, T = 10, k = 5: 10000 iterations in %.1f s (%.2f ms/iter); impacts (200 draws) %.1f s\n",
            el, 1000 * el / 10000, el_imp))
print(data.frame(variable = s$variable, ess_bulk = round(s$ess_bulk), ess_per_sec = round(s$ess_bulk / el, 1)))
