# ESS per second on the mixing cases of plan/02 §2.6 ("Mixing experiments"),
# for the gate of a sampler change (CLAUDE.md "Mixing changes"; plan/06 R6):
#   balanced:     static SAR with random effects, N = 100 (10 x 10 rook), T = 5;
#   rare:         the same with about 5% ones;
#   rare_ordered: the same, ordered (J = 3), about 5% in the top category;
#   latent_ar:    latent AR(1), 4 x 4 grid, T = 4, gamma = 1.5, near-separated.
# 4,000 draws (after 1,000 warm-up) x 2 chains on one core, `seeds` data sets
# per case. Run it once per build and compare the two .csv files.
# Writes validation/bench/mixing_<label>_<date>.csv and prints a summary.
# Usage: Rscript validation/bench/bench_mixing.R <lib> <label> [seeds]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = args[1])
label <- args[2]
seeds <- if (length(args) >= 3) as.integer(args[3]) else 5L
W10 <- psp_weights(rook_weights(10, 10), style = "W")
W4 <- psp_weights(rook_weights(4, 4), style = "W")
ctl <- function(s) psp_control(iter = 4000, warmup = 1000, chains = 2, cores = 1, seed = s,
                               progress = FALSE)
cases <- list(
  balanced = function(s) list(fun = psprobit, W = W10, args = list(), d = sim_psprobit(
    W10, T = 5, beta = c(0, 1), rho = 0.5, sigma_alpha = 0.7, seed = s)),
  rare = function(s) list(fun = psprobit, W = W10, args = list(), d = sim_psprobit(
    W10, T = 5, beta = c(-1.45, 1), rho = 0.5, sigma_alpha = 0.7, seed = s)),
  rare_ordered = function(s) list(fun = psoprobit, W = W10, args = list(), d = sim_psprobit(
    W10, T = 5, beta = c(-0.3, 1), rho = 0.5, sigma_alpha = 0.7, cutpoints = 2.2, seed = s)),
  latent_ar = function(s) list(fun = psprobit, W = W4,
                               args = list(dynamic = psp_dynamic(latent = TRUE)),
                               d = sim_psprobit(W4, T = 4, beta = c(0, 2), rho = 0.3,
                                                sigma_alpha = 0.5, gamma_z = 1.5, seed = s)))
rows <- list()
for (cn in names(cases)) {
  for (s in seq_len(seeds)) {
    cs <- cases[[cn]](100 + s)
    t0 <- proc.time()
    fit <- suppressWarnings(do.call(cs$fun, c(list(y ~ x1, cs$d, c("id", "time"), cs$W,
                                                   impacts = FALSE, control = ctl(s)), cs$args)))
    el <- (proc.time() - t0)[["elapsed"]]
    dg <- fit$diagnostics
    rows[[length(rows) + 1]] <- data.frame(case = cn, seed = s, share_top = mean(cs$d$y == max(cs$d$y)),
                                           seconds = el, variable = dg$variable,
                                           ess_bulk = dg$ess_bulk, ess_per_sec = dg$ess_bulk / el)
  }
  cat(cn, "done\n")
}
res <- do.call(rbind, rows)
out <- sprintf("validation/bench/mixing_%s_%s.csv", label, format(Sys.Date()))
utils::write.csv(res, out, row.names = FALSE)
agg <- do.call(rbind, lapply(split(res, res$case), function(x) {
  per_seed <- split(x, x$seed)
  data.frame(case = x$case[1], seconds = round(mean(sapply(per_seed, function(p) p$seconds[1])), 2),
             share_top = round(mean(sapply(per_seed, function(p) p$share_top[1])), 3),
             min_ess = round(mean(sapply(per_seed, function(p) min(p$ess_bulk)))),
             min_ess_per_sec = round(mean(sapply(per_seed, function(p) min(p$ess_per_sec))), 1))
}))
cat(sprintf("panelspatialprobit %s (%s), means over %d data sets\n",
            packageVersion("panelspatialprobit"), label, seeds))
print(agg, row.names = FALSE)
