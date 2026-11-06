# 02e_se_ratio_robust.R --------------------------------------------------------
# Robust SE ratio for the large-replicate correct-specification scenarios
# P2-01, P1-01 and P2-02 (Parts 1 and 2) and the Part 3 scenario without a
# treatment -> covariate effect (P3-01), for each of the two contrasts
# g2g1 (theta = risk(g2) - risk(g1)) and g3g2 (theta = risk(g3) - risk(g2)).
# The ratio of the mean bootstrap SE to the MC SD can be below 1 although
# coverage is near 0.95 when a few outlying replicates inflate the MC SD.
# This script reports the MC SD, a robust SD (IQR / 1.349), the mean and
# median bootstrap SE, both ratios, coverage, the kurtosis of the estimates
# and the largest standardised deviation.
#
# Inputs:  <sim_dir>/sim_results_<run>.csv of the runs that supply the ids
#          (analysis/sim_paper_results.R)
# Outputs: <sim_dir>/sim_se_ratio_robust.csv  one row per id x contrast
#          (column `contrast`: "g2g1" or "g3g2")
# CCW_SIM_DIR (default "output/sim") sets sim_dir.
#
# Run from the project root:  Rscript analysis/02e_se_ratio_robust.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))

sim_dir <- Sys.getenv("CCW_SIM_DIR", "")
if (!nzchar(sim_dir)) sim_dir <- file.path("output", "sim")
contrasts <- c("g2g1", "g3g2")

res <- read_paper_results(c("P2-01", "P1-01", "P2-02", "P3-01"),
                          sim_dir = sim_dir)
res$key <- paste(res$run, res$scenario)

out <- do.call(rbind, lapply(unique(res$key), function(k) {
  x <- res[res$key == k, ]
  do.call(rbind, lapply(contrasts, function(cc) {
    e <- x[[paste0("theta_hat_", cc)]]
    se <- x[[paste0("se_boot_", cc)]]
    cov <- x[[paste0("covered_", cc)]]
    mc_sd <- if (sum(!is.na(e)) >= 2) stats::sd(e, na.rm = TRUE) else
      NA_real_
    rob_sd <- stats::IQR(e, na.rm = TRUE) / 1.349
    z <- (e - mean(e, na.rm = TRUE)) / mc_sd
    data.frame(
      id = x$id[1], run = x$run[1], scenario = x$scenario[1],
      contrast = cc, n_rep = nrow(x),
      mc_sd = mc_sd, robust_sd = rob_sd,
      mean_se_boot = mean(se, na.rm = TRUE),
      median_se_boot = stats::median(se, na.rm = TRUE),
      se_ratio = mean(se, na.rm = TRUE) / mc_sd,
      se_ratio_robust = mean(se, na.rm = TRUE) / rob_sd,
      coverage = mean(cov, na.rm = TRUE),
      kurtosis = mean(z^4, na.rm = TRUE),
      max_abs_z = suppressWarnings(max(abs(z), na.rm = TRUE))
    )
  }))
}))

dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, file.path(sim_dir, "sim_se_ratio_robust.csv"),
                 row.names = FALSE)
print(out, digits = 3)
