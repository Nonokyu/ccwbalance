# 02c_simulation_tables_manycov.R ---------------------------------------------
# Calibration figure of the run manycov (see
# ccw_scenarios("manycov")): exceedance rate of the conditional parametric
# null reference (rrb_null()) of interval RRB with many correlated baseline
# covariates, small n, rare artificial censoring and extreme weights, and
# detection of an omitted covariate X2 correlated (rho = 0.5) or not
# (rho = 0) with X. All six scenarios of the run are shown; the manuscript
# reports only P2-05 (many_cov_n300) and P2-06 (extreme_weights), whose
# summaries are in output/sim/sim_paper_scenarios.csv
# (analysis/02b_simulation_tables_paper.R).
#
# Inputs:  output/sim/sim_results_manycov.csv, analysis/scenarios_manycov.csv
# Outputs: output/sim/sim_manycov_calibration.csv  scenario x history
#                                                  function exceedance rate,
#                                                  exact CI (id NA for the
#                                                  unreported scenarios)
#          output/figures/sim_part2_calibration_manycov.png  exceedance-rate
#                                                  plot
# CCW_OUT_SUFFIX (default "manycov") selects sim_results_<suffix>.csv; for a
# suffix other than "manycov" the outputs carry "_<suffix>" (e.g. smoke
# tests).
#
# Run from the project root:  Rscript analysis/02c_simulation_tables_manycov.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
library(ggplot2)

SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", "manycov")
out_dir <- file.path("output", "sim")
fig_dir <- file.path("output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tag <- if (SUFFIX == "manycov") "" else paste0("_", SUFFIX)
out_csv <- file.path(out_dir, sprintf("sim_manycov_calibration%s.csv", tag))
out_png <- file.path(fig_dir,
                     sprintf("sim_part2_calibration_manycov%s.png", tag))

res <- read.csv(file.path(out_dir, sprintf("sim_results_%s.csv", SUFFIX)))
scen <- read.csv(file.path("analysis", "scenarios_manycov.csv"))
scen <- scen[scen$scenario %in% unique(res$scenario), ]

feats <- sub("^null_exceed_", "",
             grep("^null_exceed_", names(res), value = TRUE))

# ---------------------------------------------------------------------------
# 1. calibration / detection: exceedance rate of the null-reference q95,
#    overall and per history function, with exact binomial 95% CI
cal <- do.call(rbind, lapply(scen$scenario, function(s) {
  x <- res[res$scenario == s, ]
  cols <- c(all = "null_exceed",
            stats::setNames(paste0("null_exceed_", feats), feats))
  do.call(rbind, lapply(names(cols), function(f) {
    v <- x[[cols[[f]]]]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NULL)
    ci <- stats::binom.test(sum(v), length(v))$conf.int
    data.frame(id = scen$id[scen$scenario == s], scenario = s, feature = f,
               n_rep = length(v), exceed_rate = mean(v), ci_lower = ci[1],
               ci_upper = ci[2])
  }))
}))
write.csv(cal, out_csv, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. figure
e_lab <- c(many_cov = "Many correlated covariates, n = 2000",
           many_cov_n300 = "Many correlated covariates, n = 300",
           rare_censoring = "Rare artificial censoring",
           extreme_weights = "Extreme weights",
           many_cov_omit_X2 = "Omitted X2 (rho = 0.5)",
           many_cov_omit_X2_rho0 = "Omitted X2 (rho = 0)")
e_names <- names(e_lab)[names(e_lab) %in% cal$scenario]
calp <- cal[cal$scenario %in% e_names, ]
calp$feature <- factor(ifelse(calp$feature == "all", "All history functions",
                              calp$feature),
                       levels = c("All history functions", feats))
# first scenario in grid order at the top
calp$scenario <- factor(e_lab[calp$scenario], levels = rev(e_lab[e_names]))
p_cal <- ggplot(calp, aes(exceed_rate, scenario)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.3) +
  facet_wrap(~feature, ncol = 4) +
  labs(x = paste("Rate of interval-RRB MaxIB (g1, g3) above the",
                 "null-reference 95th percentile (95% exact CI)"),
       y = NULL, caption = "Dashed line: nominal 0.05.") +
  theme_bw()
ggsave(out_png, p_cal, width = 12, height = 9, dpi = 150)

cat("Written:", out_csv, out_png, sep = "\n  ")
