# Cross-check against spatialprobit::sarprobit on the Katrina data (T = 1,
# pooled SAR probit, near-diffuse priors). See plan/04 s4.4.
# Usage: Rscript validation/crosscheck/spatialprobit_katrina.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else NULL
library(panelspatialprobit, lib.loc = lib)
suppressPackageStartupMessages({
  library(spatialprobit)
  library(spdep)
})
data(Katrina, package = "spatialprobit")
nb <- suppressWarnings(knn2nb(knearneigh(cbind(Katrina$lat, Katrina$long), k = 11)))
W <- as(spatialreg::as_dgRMatrix_listw(nb2listw(nb, style = "W")), "CsparseMatrix")
Katrina$id <- seq_len(nrow(Katrina))
Katrina$time <- 1L
W <- psp_weights(W, ids = Katrina$id)
f <- y1 ~ flood_depth + log_medinc + small_size + large_size + low_status_customers +
  high_status_customers + owntype_sole_proprietor + owntype_national_chain

# spatialprobit's default start values (beta = 0, rho = 0.75) give NaN latent
# draws on these data (spatialprobit 1.0.4, tmvtnorm 1.7, Matrix 1.7); start it
# from the non-spatial probit estimates instead.
b0 <- coef(glm(f, family = binomial("probit"), data = Katrina))
set.seed(1)
t0 <- Sys.time()
sp <- suppressWarnings(sarprobit(f, W = W, data = Katrina, ndraw = 6000, burn.in = 1000,
                                 thinning = 1, computeMarginalEffects = FALSE,
                                 start = list(rho = 0.3, beta = b0)))
stopifnot(!anyNA(sp$B))
t_sp <- Sys.time() - t0

t0 <- Sys.time()
ps <- suppressWarnings(psprobit(f, data = Katrina, index = c("id", "time"), W = W,
                                effects = "pooled", impacts = FALSE,
                                priors = psp_priors(beta = prior_normal(0, 1e4),
                                                    rho = prior_beta(1, 1)),
                                control = psp_control(iter = 3000, warmup = 1000, chains = 2,
                                                      seed = 1)))
t_ps <- Sys.time() - t0

s <- posterior::summarise_draws(ps$draws, "mean", "sd", "ess_bulk")
tab <- data.frame(variable = s$variable, panelspatialprobit = s$mean, sd_psp = s$sd,
                  ess_psp = s$ess_bulk, spatialprobit = colMeans(sp$B),
                  sd_sp = apply(sp$B, 2, sd))
tab$diff_in_sd <- (tab$panelspatialprobit - tab$spatialprobit) / tab$sd_psp
print(tab, digits = 3, row.names = FALSE)
cat(sprintf("\nTime: spatialprobit %s, panelspatialprobit %s\n", format(t_sp, digits = 3),
            format(t_ps, digits = 3)))
