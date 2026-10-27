# 02d_size_adjusted_detection.R ------------------------------------------------
# Size-adjusted comparison of the detection rates of interval RRB and
# censoring balance (CB) in the misspecification factorial (Part 1,
# P1-02 to P1-13 of paper_scenarios()).
# The two null references have different empirical false-positive rates under
# the correctly specified model (RRB about 0.023, CB about 0.050 at n = 2000,
# P1-01), so their nominal-5% detection rates are not
# comparable. Here the per-replicate threshold c * q95 (q95 = the replicate's
# own null 95th percentile) is calibrated so that the empirical false-positive
# rate under correct specification (P1-01, n = 2000, 1,000
# replicates) equals 0.05 for each statistic, and the same multiplier is
# applied to the Part 1 scenarios (n = 2000, 200 replicates).
#
# Inputs:  output/sim/sim_results_<run>.csv of the runs that supply P1-01 to
#          P1-13 (analysis/sim_paper_results.R)
# Outputs: output/sim/sim_part1_detection_size_adjusted.csv
#          output/sim/sim_part1_detection_size_multipliers.csv
# Run from the project root:  Rscript analysis/02d_size_adjusted_detection.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))
out_dir <- file.path("output", "sim")
scen <- paper_scenarios()
p1_ids <- scen$id[scen$part == 1]
d <- read_paper_results(p1_ids)

cal <- d[d$id == "P1-01", ]
# multiplier c such that mean(obs > c * q95) = 0.05: the 95th percentile of
# the ratio obs / q95 under the correct model
ratio_rrb <- cal$null_maxib_obs / cal$null_q95
ratio_cb <- cal$null_cb_maxib_obs / cal$null_cb_q95
c_rrb <- stats::quantile(ratio_rrb, 0.95, names = FALSE)
c_cb <- stats::quantile(ratio_cb, 0.95, names = FALSE)

rate <- function(x, c, obs, q95) {
  v <- as.numeric(x[[obs]] > c * x[[q95]])
  c(rate = mean(v), se = sqrt(mean(v) * (1 - mean(v)) / length(v)), v = list(v))
}
rows <- lapply(p1_ids, function(s) {
  x <- d[d$id == s, ]
  r <- rate(x, c_rrb, "null_maxib_obs", "null_q95")
  b <- rate(x, c_cb, "null_cb_maxib_obs", "null_cb_q95")
  dv <- r$v - b$v
  info <- scen[scen$id == s, ]
  data.frame(id = s, label = info$label, run = info$run,
             scenario = info$scenario,
             misspec = info$misspec, cens_strength = info$cens_strength,
             outcome_strength = info$outcome_strength, n_rep = nrow(x),
             nominal_rrb = mean(x$null_exceed),
             nominal_cb = mean(x$null_cb_exceed),
             adj_rrb = r$rate, adj_rrb_se = r$se, adj_cb = b$rate,
             adj_cb_se = b$se, diff_rrb_minus_cb = mean(dv),
             diff_se = stats::sd(dv) / sqrt(length(dv)))
})
res <- do.call(rbind, rows)
attr_row <- data.frame(c_rrb = c_rrb, c_cb = c_cb)
write.csv(res, file.path(out_dir, "sim_part1_detection_size_adjusted.csv"),
          row.names = FALSE)
write.csv(attr_row,
          file.path(out_dir, "sim_part1_detection_size_multipliers.csv"),
          row.names = FALSE)
cat(sprintf(paste("Multipliers c (obs > c * q95 gives FPR 0.05 under P1-01):",
                  "RRB %.3f, CB %.3f\n"), c_rrb, c_cb))
print(res[, setdiff(names(res), c("label", "run", "scenario"))],
      digits = 3, row.names = FALSE)
