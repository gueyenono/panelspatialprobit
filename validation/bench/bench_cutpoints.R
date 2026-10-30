# Roadmap 3.4: ESS per second of the two cut-point samplers ("transform" =
# Albert-Chib step + transformation move, "albert_chib" = Gibbs step only) on
# ordered SAR random-effects models, N = 400 (20 x 20 rook grid), T = 5.
# Writes validation/bench/cutpoints_<date>.txt (report) and .csv (per parameter).
# Usage: Rscript validation/bench/bench_cutpoints.R [lib]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args)) args[1])
W <- psp_weights(rook_weights(20, 20), style = "W")
designs <- list(J3 = 0.8, J5 = c(0.5, 1.0, 1.6))
rows <- list()
for (dn in names(designs)) {
  d <- sim_psprobit(W, T = 5, beta = c(0.3, 1, -0.5), rho = 0.5, sigma_alpha = 0.7,
                    cutpoints = designs[[dn]], seed = 10)
  for (smp in c("transform", "albert_chib")) {
    t0 <- proc.time()
    fit <- suppressWarnings(psoprobit(
      y ~ x1 + x2, d, c("id", "time"), W, impacts = FALSE,
      control = psp_control(iter = 3000, warmup = 1000, chains = 4, seed = 1,
                            cutpoint_sampler = smp, progress = FALSE)))
    el <- (proc.time() - t0)[["elapsed"]]
    s <- fit$diagnostics
    rows[[length(rows) + 1]] <- data.frame(design = dn, sampler = smp, variable = s$variable,
                                           ess_bulk = round(s$ess_bulk),
                                           ess_per_sec = round(s$ess_bulk / el, 1),
                                           seconds = round(el, 1))
  }
}
res <- do.call(rbind, rows)
# per-parameter results for validation/vignette_data.R
utils::write.csv(res, sprintf("validation/bench/cutpoints_%s.csv", format(Sys.Date())),
                 row.names = FALSE)
out <- sprintf("validation/bench/cutpoints_%s.txt", format(Sys.Date()))
sink(out)
cat("Cut-point samplers: N = 400, T = 5, 4 chains x 3000 draws (1000 warm-up)\n\n")
for (dn in names(designs)) {
  a <- res[res$design == dn & res$sampler == "transform", ]
  b <- res[res$design == dn & res$sampler == "albert_chib", ]
  cat(sprintf("[%s] seconds: transform %.1f, albert_chib %.1f\n", dn, a$seconds[1], b$seconds[1]))
  print(data.frame(variable = a$variable, ess_transform = a$ess_bulk,
                   ess_ac = b$ess_bulk, ess_s_transform = a$ess_per_sec,
                   ess_s_ac = b$ess_per_sec, ratio = round(a$ess_per_sec / b$ess_per_sec, 1)),
        row.names = FALSE)
  cat("\n")
}
sink()
cat(readLines(out), sep = "\n")
