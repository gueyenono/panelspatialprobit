# Roadmap 3.6: cross-check the ordered sampler against MCMCpack::MCMCoprobit on
# a non-spatial pooled ordered probit (same parameterisation: intercept and
# first cut-point 0). Near-flat priors on both sides.
# Usage: Rscript validation/crosscheck/mcmcoprobit.R [lib]
args <- commandArgs(trailingOnly = TRUE)
library(panelspatialprobit, lib.loc = if (length(args)) args[1])
suppressPackageStartupMessages(library(MCMCpack))
W <- psp_weights(rook_weights(15, 20), style = "W")
d <- sim_psprobit(W, T = 4, beta = c(0.4, 1, -0.7), rho = 0, sigma_alpha = 0,
                  cutpoints = c(0.9, 1.7), seed = 3)
fit <- suppressWarnings(psoprobit(
  y ~ x1 + x2, d, c("id", "time"), W, model = "none", effects = "pooled", impacts = FALSE,
  priors = psp_priors(beta = prior_normal(0, 100), cutpoints = prior_flat()),
  control = psp_control(iter = 5000, warmup = 1000, chains = 4, seed = 1, progress = FALSE)))
mc <- MCMCoprobit(y ~ x1 + x2, data = d, burnin = 2000, mcmc = 40000, thin = 2, seed = 1,
                  tune = 0.3)
ours <- posterior::summarise_draws(fit$draws, "mean", "sd")
theirs <- data.frame(mean = colMeans(mc), sd = apply(mc, 2, sd))
cmp <- data.frame(parameter = ours$variable, psp_mean = ours$mean, psp_sd = ours$sd,
                  mcmcpack_mean = theirs$mean, mcmcpack_sd = theirs$sd)
cmp$diff_in_sd <- (cmp$psp_mean - cmp$mcmcpack_mean) / cmp$mcmcpack_sd
out <- sprintf("validation/crosscheck/results/mcmcoprobit_%s.txt", format(Sys.Date()))
sink(out)
cat("psoprobit(model = 'none', effects = 'pooled') vs MCMCpack::MCMCoprobit\n")
cat("N = 300, T = 4 (1200 obs), J = 4; MCMCpack columns: (Intercept), x1, x2, gamma2, gamma3\n\n")
print(cmp, digits = 3, row.names = FALSE)
sink()
cat(readLines(out), sep = "\n")
