# 02e_se_ratio_robust.R --------------------------------------------------------
# Robust SE ratio for the large-replicate correct-specification scenarios
# P2-01, P1-01 and P2-02 (Parts 1 and 2) and the Part 3 scenario without a
# treatment -> covariate effect (P3-01). The ratio of the mean bootstrap SE
# to the MC SD is below 1 in several of them although coverage is near 0.95;
# a few outlying replicates inflate the MC SD. This script reports the MC
# SD, a robust SD (IQR / 1.349), the mean and median bootstrap SE, both
# ratios, coverage, the kurtosis of the estimates and the largest
# standardised deviation.
#
# Inputs:  output/sim/sim_results_<run>.csv of the runs that supply the ids
#          (analysis/sim_paper_results.R)
# Outputs: output/sim/sim_se_ratio_robust.csv
#
# Run from the project root:  Rscript analysis/02e_se_ratio_robust.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))

cols <- c("id", "run", "scenario", "theta_hat", "se_boot", "covered")
res <- read_paper_results(c("P2-01", "P1-01", "P2-02", "P3-01"))[, cols]
res$key <- paste(res$run, res$scenario)

out <- do.call(rbind, lapply(unique(res$key), function(k) {
  x <- res[res$key == k, ]
  e <- x$theta_hat
  mc_sd <- stats::sd(e)
  rob_sd <- stats::IQR(e) / 1.349
  z <- (e - mean(e)) / mc_sd
  data.frame(
    id = x$id[1], run = x$run[1], scenario = x$scenario[1], n_rep = nrow(x),
    mc_sd = mc_sd, robust_sd = rob_sd,
    mean_se_boot = mean(x$se_boot), median_se_boot = stats::median(x$se_boot),
    se_ratio = mean(x$se_boot) / mc_sd,
    se_ratio_robust = mean(x$se_boot) / rob_sd,
    coverage = mean(x$covered),
    kurtosis = mean(z^4), max_abs_z = max(abs(z))
  )
}))

utils::write.csv(out, "output/sim/sim_se_ratio_robust.csv", row.names = FALSE)
print(out, digits = 3)
