# Compact summaries of the validation results for vignette("validation"),
# which cannot read validation/ (excluded from the build).
# Usage: Rscript validation/vignette_data.R
# Writes vignettes/validation-sbc.csv, validation-mc.csv, validation-wk.csv,
# validation-bench.csv and validation-cutpoints.csv.

# SBC: every run, labelled by its file, as 10-bin rank histograms (ranks are in
# 0..400) with the chi-square p-value used in the text reports.
L <- 400
files <- list.files("validation/sbc/results", "^sbc_2.*\\.rds$", full.names = TRUE)
sbc <- do.call(rbind, lapply(files, function(f) {
  res <- readRDS(f)
  offset <- sub(".*_offset([0-9]+)\\.rds$", "\\1", f)
  offset <- if (offset == f) 0L else as.integer(offset)
  do.call(rbind, lapply(names(res), function(vn) {
    r <- res[[vn]]
    do.call(rbind, lapply(setdiff(colnames(r), "max_rhat"), function(p) {
      bins <- tabulate(cut(r[, p], breaks = seq(-0.5, L + 0.5, length.out = 11),
                           labels = FALSE), 10)
      data.frame(variant = vn, seed_offset = offset, parameter = p, reps = nrow(r),
                 median_rhat = round(stats::median(r[, "max_rhat"]), 4),
                 p = signif(suppressWarnings(stats::chisq.test(bins)$p.value), 3),
                 t(stats::setNames(bins, paste0("b", 1:10))))
    }))
  }))
}))
# a variant run twice with the same seeds is the same result: keep one copy
sbc <- sbc[!duplicated(sbc[, c("variant", "seed_offset", "parameter")]), ]
utils::write.csv(sbc, "vignettes/validation-sbc.csv", row.names = FALSE)

# Monte Carlo: bias, RMSE, 90% coverage and width per design and parameter.
mc_file <- tail(sort(list.files("validation/montecarlo/results", "^mc_.*\\.rds$",
                                full.names = TRUE)), 1)
res <- readRDS(mc_file)
mc <- do.call(rbind, lapply(names(res), function(dn) {
  ok <- Filter(function(r) is.list(r) && !inherits(r, "try-error"), res[[dn]])
  do.call(rbind, lapply(names(ok[[1]]), function(p) {
    m <- do.call(rbind, lapply(ok, `[[`, p))
    err <- m[, "mean"] - m[, "true"]
    data.frame(design = dn, parameter = p, reps = length(ok),
               true = round(mean(m[, "true"]), 3), bias = round(mean(err), 3),
               rmse = round(sqrt(mean(err^2)), 3),
               cover90 = round(mean(m[, "lo"] <= m[, "true"] & m[, "true"] <= m[, "hi"]), 3),
               width = round(mean(m[, "hi"] - m[, "lo"]), 3))
  }))
}))
utils::write.csv(mc, "vignettes/validation-mc.csv", row.names = FALSE)
# Wang & Kockelman (2009) replication: per cell and parameter, the average
# posterior mean, RMSE and 90% coverage, next to the paper's Table 1 averages.
wk_file <- tail(sort(list.files("validation/wang_kockelman/results", "^wk_.*\\.rds$",
                                full.names = TRUE)), 1)
wk_res <- readRDS(wk_file)
wk_res$cover <- wk_res$lo90 <= wk_res$true & wk_res$true <= wk_res$hi90
wk <- do.call(rbind, lapply(split(wk_res, list(wk_res$parameter, wk_res$rho_true,
                                               wk_res$lambda_true), drop = TRUE),
                            function(x) {
  data.frame(lambda = x$lambda_true[1], rho = x$rho_true[1], parameter = x$parameter[1],
             true = x$true[1], mean = round(mean(x$mean), 3),
             rmse = round(sqrt(mean((x$mean - x$true)^2)), 3), cover90 = round(mean(x$cover), 2),
             rhat_share_1.05 = round(mean(x$max_rhat >= 1.05), 2))
}))
# the paper's averages, from its Table 1
paper <- utils::read.csv("validation/wang_kockelman/wk_table1.csv", comment.char = "#")
paper <- reshape(paper, direction = "long", varying = names(paper)[-(1:2)],
                 v.names = "paper", timevar = "parameter", times = names(paper)[-(1:2)],
                 idvar = c("cell_lambda", "cell_rho"))
names(paper)[1:2] <- c("lambda", "rho")
wk <- merge(wk, paper, by = c("lambda", "rho", "parameter"))
wk <- wk[order(wk$lambda, wk$rho, match(wk$parameter, unique(paper$parameter))), ]
utils::write.csv(wk, "vignettes/validation-wk.csv", row.names = FALSE)

# Speed and mixing (validation/bench/bench_models.R): per model and size, the
# time per iteration and the range of ESS per second over the parameters.
bench_file <- tail(sort(list.files("validation/bench", "^models_.*\\.csv$", full.names = TRUE)), 1)
bres <- utils::read.csv(bench_file, check.names = FALSE)
bench <- do.call(rbind, lapply(split(bres, list(bres$design, bres$N), drop = TRUE), function(x) {
  data.frame(design = x$design[1], N = x$N[1], T = x$T[1], n = x$n[1],
             iterations = x$iterations[1], ms_per_iter = x$ms_per_iter[1],
             ess_per_sec_min = min(x$ess_per_sec), ess_per_sec_max = max(x$ess_per_sec),
             slowest = x$variable[which.min(x$ess_per_sec)])
}))
bench <- bench[order(bench$N, match(bench$design, c("sar_binary", "sar_ordered", "sem_binary"))), ]
utils::write.csv(bench, "vignettes/validation-bench.csv", row.names = FALSE)

# Cut-point samplers (validation/bench/bench_cutpoints.R): ESS per second of
# the transformation move and of the Albert-Chib step alone, per parameter.
cut_file <- tail(sort(list.files("validation/bench", "^cutpoints_.*\\.csv$", full.names = TRUE)), 1)
cres <- utils::read.csv(cut_file, check.names = FALSE)
tr <- cres[cres$sampler == "transform", ]
ac <- cres[cres$sampler == "albert_chib", ]
ac <- ac[match(paste(tr$design, tr$variable), paste(ac$design, ac$variable)), ]
cuts <- data.frame(design = tr$design, parameter = tr$variable,
                   ess_transform = tr$ess_bulk, ess_albert_chib = ac$ess_bulk,
                   ess_s_transform = tr$ess_per_sec, ess_s_albert_chib = ac$ess_per_sec,
                   ratio = round(tr$ess_per_sec / ac$ess_per_sec, 1))
utils::write.csv(cuts, "vignettes/validation-cutpoints.csv", row.names = FALSE)
cat(sprintf("%d SBC rows, %d Monte Carlo rows, %d Wang-Kockelman rows, %d benchmark rows, %d cut-point rows\n",
            nrow(sbc), nrow(mc), nrow(wk), nrow(bench), nrow(cuts)))
