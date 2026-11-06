# 02e_se_ratio_robust.R --------------------------------------------------------
# Robust SE ratio (MC SD vs IQR / 1.349 vs bootstrap SE), coverage and
# kurtosis for P2-01, P1-01, P2-02 and P3-01, per contrast g2g1 and g3g2.
# Inputs:  output/sim/sim_results_<run>.csv (analysis/sim_paper_results.R)
# Outputs: output/sim/sim_se_ratio_robust.csv
# Run:     Rscript analysis/02e_se_ratio_robust.R (CCW_SIM_DIR)

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
