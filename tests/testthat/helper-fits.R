small_data <- function(seed = 1, cutpoints = NULL, T = 3, rho = 0.4, sigma_alpha = 0.5,
                       theta = NULL) {
  W <- psp_weights(rook_weights(4, 4), style = "W")
  d <- sim_psprobit(W, T = T, beta = c(0, 1, -0.5), rho = rho, sigma_alpha = sigma_alpha,
                    theta = theta, cutpoints = cutpoints, seed = seed)
  list(d = d, W = W)
}

fast_control <- function(...) psp_control(iter = 200, warmup = 100, chains = 2, seed = 3, ...)

# Time-varying weights: rook contiguity on a 4 x 4 grid with a different
# random tenth of the links removed in each period, row-standardised.
tv_weights <- function(T, seed = 1) {
  B <- rook_weights(4, 4)
  set.seed(seed)
  lapply(seq_len(T), function(t) {
    Bt <- as.matrix(B)
    links <- which(upper.tri(Bt) & Bt > 0)
    drop <- sample(links, round(length(links) / 10))
    Bt[drop] <- 0
    Bt[lower.tri(Bt)] <- t(Bt)[lower.tri(Bt)]
    psp_weights(Bt, style = "W")
  })
}

tv_data <- function(seed = 1, T = 3, rho = 0.4, sigma_alpha = 0.5, lambda = 0,
                    cutpoints = NULL) {
  W <- tv_weights(T, seed)
  d <- sim_psprobit(W, T = T, beta = c(0, 1, -0.5), rho = rho, sigma_alpha = sigma_alpha,
                    lambda = lambda, cutpoints = cutpoints, seed = seed)
  list(d = d, W = W)
}
