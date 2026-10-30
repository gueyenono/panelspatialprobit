# SBC experiment on spatialprobit::sar_probit_mcmc (plan/04 s4.2, roadmap 1.9),
# with panelspatialprobit on the identical design as a control.
# Cross-section (T = 1), 8 x 8 rook grid, pooled SAR probit,
# priors beta ~ N(0, I), rho ~ Beta(2, 2) on (-1, 1).
# Usage: Rscript validation/sbc/sbc_spatialprobit.R [reps] [cores] [lib]
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores <- if (length(args) >= 2) as.integer(args[2]) else 4L
lib <- if (length(args) >= 3) args[3] else NULL
library(panelspatialprobit, lib.loc = lib)
suppressPackageStartupMessages(library(spatialprobit))

W <- psp_weights(rook_weights(8, 8), style = "W")
L <- 400

one_rep <- function(rep) {
  set.seed(7e5 + rep)
  repeat {
    beta <- rnorm(3)
    rho <- 2 * rbeta(1, 2, 2) - 1
    d <- sim_psprobit(W, T = 1, beta = beta, rho = rho, sigma_alpha = 0)
    if (length(unique(d$y)) == 2) break
  }
  truth <- c(beta, rho)
  X <- cbind(1, d$x1, d$x2)
  sp <- tryCatch(suppressWarnings(sar_probit_mcmc(
    d$y, X, W, ndraw = L, burn.in = 1000, thinning = 10,
    prior = list(a1 = 2, a2 = 2, c = rep(0, 3), T = diag(3), lflag = 0),
    start = list(rho = 0, beta = rep(0, 3)), computeMarginalEffects = FALSE)),
    error = function(e) NULL)
  r_sp <- if (is.null(sp) || anyNA(sp$B)) rep(NA, 4) else colSums(sweep(sp$B, 2, truth, "<"))
  fit <- suppressWarnings(psprobit(y ~ x1 + x2, d, c("id", "time"), W, effects = "pooled",
                                   impacts = FALSE,
                                   priors = psp_priors(beta = prior_normal(0, 1),
                                                       rho = prior_beta(2, 2)),
                                   control = psp_control(iter = 4000, warmup = 1000, thin = 10,
                                                         chains = 1, seed = rep)))
  dm <- unclass(posterior::as_draws_matrix(fit$draws))
  r_psp <- colSums(sweep(dm, 2, truth, "<"))
  c(sp = r_sp, psp = r_psp)
}

out <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = cores)
res <- do.call(rbind, out[!vapply(out, inherits, logical(1), "try-error")])
dir.create("validation/sbc/results", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date())
saveRDS(res, sprintf("validation/sbc/results/sbc_spatialprobit_%s.rds", stamp))
pars <- c("(Intercept)", "x1", "x2", "rho")
report <- sprintf("SBC spatialprobit vs panelspatialprobit %s: %d reps, %d draws, N = 64, T = 1",
                  stamp, nrow(res), L)
report <- c(report, sprintf("spatialprobit failed (NaN/error) in %d reps", sum(is.na(res[, 1]))))
for (pkg in c("sp", "psp")) {
  for (j in seq_along(pars)) {
    x <- res[, j + 4 * (pkg == "psp")]
    x <- x[!is.na(x)]
    bins <- cut(x, breaks = seq(-0.5, L + 0.5, length.out = 11), labels = FALSE)
    pv <- chisq.test(tabulate(bins, 10))$p.value
    report <- c(report, sprintf("  %-4s %-12s chi2 p = %.4f  mean rank/L = %.3f  bins: %s", pkg,
                                pars[j], pv, mean(x) / L, paste(tabulate(bins, 10), collapse = " ")))
  }
}
writeLines(report, sprintf("validation/sbc/results/sbc_spatialprobit_%s.txt", stamp))
cat(report, sep = "\n")
