# Cross-check of the ordered SAR model against spatialprobit::sarorderedprobit
# (T = 1, pooled SAR ordered probit, J = 4, near-diffuse priors on both
# sides). Same parameterisation: intercept, first cut-point fixed at 0.
# See plan/04 s4.4. spatialprobit's sampler has the issues listed in plan/02
# s2.10 (rho kernel, block order, restarted latent draws), and its cut-points
# use the Albert-Chib step, which mixes slowly, so it gets long chains.
# Usage: Rscript validation/crosscheck/spatialprobit_ordered.R [lib] [cores]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else NULL
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
library(panelspatialprobit, lib.loc = lib)
suppressPackageStartupMessages(library(spatialprobit))

W <- psp_weights(rook_weights(20, 20), style = "W")
truth <- c("(Intercept)" = 0.3, x1 = 1, x2 = -0.7, rho = 0.5, "cut[2|3]" = 0.9,
           "cut[3|4]" = 1.7)
d <- sim_psprobit(W, T = 1, beta = truth[1:3], rho = truth[["rho"]], sigma_alpha = 0,
                  cutpoints = truth[5:6], seed = 3)
d$yi <- as.integer(d$y)

# spatialprobit: 4 chains from the non-spatial ordered probit estimates
# (converted to the intercept / first cut-point 0 parameterisation)
pl <- MASS::polr(factor(yi) ~ x1 + x2, data = d, method = "probit")
start <- list(rho = 0, beta = c(-pl$zeta[[1]], coef(pl)),
              phi = c(-Inf, pl$zeta - pl$zeta[[1]], Inf))
nd <- 50000
t0 <- Sys.time()
sp <- parallel::mclapply(1:4, function(ch) {
  set.seed(ch)
  suppressWarnings(sarorderedprobit(yi ~ x1 + x2, W = as(W, "CsparseMatrix"), data = d,
                                    ndraw = nd, burn.in = 5000, start = start,
                                    computeMarginalEffects = FALSE))$B
}, mc.cores = cores)
t_sp <- Sys.time() - t0
stopifnot(!anyNA(unlist(sp)))
# columns beta_1..3, rho, y>=2 (fixed at 0), y>=3, y>=4; array [draw, chain, variable]
sp_arr <- posterior::as_draws_array(aperm(array(
  unlist(lapply(sp, function(B) B[, c(1:4, 6:7)])), c(nd, 6, 4)), c(1, 3, 2)))
posterior::variables(sp_arr) <- names(truth)

t0 <- Sys.time()
ps <- psoprobit(yi ~ x1 + x2, data = d, index = c("id", "time"), W = W, effects = "pooled",
                impacts = FALSE,
                priors = psp_priors(beta = prior_normal(0, 1e4), rho = prior_beta(1, 1),
                                    cutpoints = prior_flat()),
                control = psp_control(iter = 10000, warmup = 2000, chains = 4, cores = cores,
                                      seed = 1, progress = FALSE))
t_ps <- Sys.time() - t0

a <- posterior::summarise_draws(ps$draws, "mean", "sd", "rhat", "ess_bulk")
b <- posterior::summarise_draws(sp_arr, "mean", "sd", "rhat", "ess_bulk")
tab <- data.frame(parameter = a$variable, truth = unname(truth[a$variable]),
                  psp_mean = a$mean, psp_sd = a$sd, psp_rhat = a$rhat, psp_ess = a$ess_bulk,
                  sp_mean = b$mean[match(a$variable, b$variable)],
                  sp_sd = b$sd[match(a$variable, b$variable)],
                  sp_rhat = b$rhat[match(a$variable, b$variable)],
                  sp_ess = b$ess_bulk[match(a$variable, b$variable)])
tab$diff_in_sd <- (tab$psp_mean - tab$sp_mean) / tab$psp_sd
tab$diff_in_mcse <- (tab$psp_mean - tab$sp_mean) /
  sqrt(tab$psp_sd^2 / tab$psp_ess + tab$sp_sd^2 / tab$sp_ess)

out <- sprintf("validation/crosscheck/results/spatialprobit_ordered_%s.txt", format(Sys.Date()))
sink(out)
cat("psoprobit(model = 'sar', effects = 'pooled') vs spatialprobit::sarorderedprobit",
    sprintf("(spatialprobit %s)\n", packageVersion("spatialprobit")))
cat(sprintf("20 x 20 rook grid, T = 1, J = 4; category shares %s\n",
            paste(round(prop.table(table(d$yi)), 2), collapse = " / ")))
cat("Priors: beta N(0, 1e4^2) here, N(0, 1e12) there; rho uniform on (-1, 1); flat cut-points.\n")
cat(sprintf("Draws: here 4 x 8,000 after 2,000 warm-up; there 4 x %s after 5,000 burn-in.\n",
            format(nd, big.mark = ",")))
cat("diff_in_mcse: difference of the posterior means over its Monte Carlo standard error.\n\n")
print(tab, digits = 3, row.names = FALSE)
cat(sprintf("\nTime: spatialprobit %s (4 chains in parallel), panelspatialprobit %s\n",
            format(t_sp, digits = 3), format(t_ps, digits = 3)))
sink()
cat(readLines(out), sep = "\n")
