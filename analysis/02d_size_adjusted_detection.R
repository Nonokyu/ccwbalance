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
# The comparison is made for the pooled statistics (strategy "all": maximum
# over all strategies, columns null_maxib_obs, null_q95, ...) and separately
# for each strategy g in g1, g2, g3 (columns null_maxib_obs_<g>,
# null_q95_<g>, ...). Each strategy has its own multipliers c, calibrated
# from that strategy's P1-01 ratios obs_<g> / q95_<g>.
# The same calibration is made for the studentised maximum of interval RRB
# (null_z_max[_<g>] against c * null_z_q95[_<g>]) and for the studentised
# maximum of CB restricted to the censorable rows (null_cb_cens_z_max[_<g>]
# against c * null_cb_cens_z_q95[_<g>]), the exceedance standard of the
# manuscript; these columns (prefixes stud_ and cb_cens_stud_) are appended
# after the raw MaxIB columns and are NA when the results lack them.
#
# Inputs:  <sim_dir>/sim_results_<run>.csv of the runs that supply P1-01 to
#          P1-13 (analysis/sim_paper_results.R)
# Outputs: <sim_dir>/sim_part1_detection_size_adjusted.csv
#            one row per strategy ("all", g1, g2, g3) x scenario id;
#            column strategy follows id; rows ordered by strategy, then id;
#            raw MaxIB: nominal_rrb, nominal_cb, adj_rrb(_se), adj_cb(_se),
#            diff_rrb_minus_cb, diff_se; then studentised: nominal_stud_rrb,
#            nominal_cb_cens_stud, stud_adj_rrb(_se), cb_cens_stud_adj
#            (_se), stud_diff_rrb_minus_cb_cens, stud_diff_se
#          <sim_dir>/sim_part1_detection_size_multipliers.csv
#            one row: c_rrb, c_cb (pooled), then c_rrb_<g>, c_cb_<g>; then
#            c_stud_rrb, c_cb_cens_stud (pooled) and c_stud_rrb_<g>,
#            c_cb_cens_stud_<g>
# CCW_SIM_DIR (default "output/sim") sets sim_dir.
# Run from the project root:  Rscript analysis/02d_size_adjusted_detection.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))
default_sim_dir <- file.path("output", "sim")
sim_dir <- Sys.getenv("CCW_SIM_DIR", "")
if (!nzchar(sim_dir)) sim_dir <- default_sim_dir
out_dir <- sim_dir
scen <- paper_scenarios()
p1_ids <- scen$id[scen$part == 1]
d <- read_paper_results(p1_ids, sim_dir = sim_dir)
strategies <- c("all", "g1", "g2", "g3")

# column names of the null-reference statistics for one strategy ("all" =
# the pooled columns without suffix)
null_cols <- function(g) {
  sfx <- if (g == "all") "" else paste0("_", g)
  stats::setNames(
    paste0(c("null_maxib_obs", "null_q95", "null_exceed", "null_cb_maxib_obs",
             "null_cb_q95", "null_cb_exceed"), sfx),
    c("obs", "q95", "exceed", "cb_obs", "cb_q95", "cb_exceed")
  )
}
# the same for the studentised maxima (RRB; CB on the censorable rows)
stud_cols <- function(g) {
  sfx <- if (g == "all") "" else paste0("_", g)
  stats::setNames(
    paste0(c("null_z_max", "null_z_q95", "null_z_exceed",
             "null_cb_cens_z_max", "null_cb_cens_z_q95",
             "null_cb_cens_z_exceed"), sfx),
    c("obs", "q95", "exceed", "cb_obs", "cb_q95", "cb_exceed")
  )
}
# a column, or NA when it is absent (results before 8a2bd6a)
col_get <- function(x, nm) if (nm %in% names(x)) x[[nm]] else
  rep(NA_real_, nrow(x))
q95_na <- function(v) {
  if (all(is.na(v))) NA_real_ else
    stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE)
}
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

cal <- d[d$id == "P1-01", ]
# multiplier c such that mean(obs > c * q95) = 0.05: the 95th percentile of
# the ratio obs / q95 under the correct model
multipliers <- lapply(strategies, function(g) {
  cl <- null_cols(g)
  ratio_rrb <- cal[[cl[["obs"]]]] / cal[[cl[["q95"]]]]
  ratio_cb <- cal[[cl[["cb_obs"]]]] / cal[[cl[["cb_q95"]]]]
  c(rrb = stats::quantile(ratio_rrb, 0.95, names = FALSE),
    cb = stats::quantile(ratio_cb, 0.95, names = FALSE))
})
names(multipliers) <- strategies
multipliers_stud <- lapply(strategies, function(g) {
  cl <- stud_cols(g)
  c(rrb = q95_na(col_get(cal, cl[["obs"]]) / col_get(cal, cl[["q95"]])),
    cb = q95_na(col_get(cal, cl[["cb_obs"]]) / col_get(cal, cl[["cb_q95"]])))
})
names(multipliers_stud) <- strategies

# NA-tolerant version for the studentised columns: replicates with NA are
# dropped (separately per statistic; the paired difference uses the
# replicates where both are available)
rate_na <- function(x, c, obs, q95) {
  v <- as.numeric(col_get(x, obs) > c * col_get(x, q95))
  n <- sum(!is.na(v))
  m <- mean_na(v)
  c(rate = m, se = if (n > 0) sqrt(m * (1 - m) / n) else NA_real_,
    v = list(v))
}
rate <- function(x, c, obs, q95) {
  v <- as.numeric(x[[obs]] > c * x[[q95]])
  c(rate = mean(v), se = sqrt(mean(v) * (1 - mean(v)) / length(v)), v = list(v))
}
detection_rows <- function(g) {
  cl <- null_cols(g)
  c_g <- multipliers[[g]]
  rows <- lapply(p1_ids, function(s) {
    x <- d[d$id == s, ]
    r <- rate(x, c_g[["rrb"]], cl[["obs"]], cl[["q95"]])
    b <- rate(x, c_g[["cb"]], cl[["cb_obs"]], cl[["cb_q95"]])
    dv <- r$v - b$v
    info <- scen[scen$id == s, ]
    data.frame(id = s, strategy = g, label = info$label, run = info$run,
               scenario = info$scenario,
               misspec = info$misspec, cens_strength = info$cens_strength,
               outcome_strength = info$outcome_strength, n_rep = nrow(x),
               nominal_rrb = mean(x[[cl[["exceed"]]]]),
               nominal_cb = mean(x[[cl[["cb_exceed"]]]]),
               adj_rrb = r$rate, adj_rrb_se = r$se, adj_cb = b$rate,
               adj_cb_se = b$se, diff_rrb_minus_cb = mean(dv),
               diff_se = stats::sd(dv) / sqrt(length(dv)),
               stud_row(x, g))
  })
  do.call(rbind, rows)
}
# studentised columns of one strategy x scenario (one-row data frame)
stud_row <- function(x, g) {
  cl <- stud_cols(g)
  c_g <- multipliers_stud[[g]]
  r <- rate_na(x, c_g[["rrb"]], cl[["obs"]], cl[["q95"]])
  b <- rate_na(x, c_g[["cb"]], cl[["cb_obs"]], cl[["cb_q95"]])
  dv <- r$v - b$v
  dv <- dv[!is.na(dv)]
  data.frame(nominal_stud_rrb = mean_na(col_get(x, cl[["exceed"]])),
             nominal_cb_cens_stud = mean_na(col_get(x, cl[["cb_exceed"]])),
             stud_adj_rrb = r$rate, stud_adj_rrb_se = r$se,
             cb_cens_stud_adj = b$rate, cb_cens_stud_adj_se = b$se,
             stud_diff_rrb_minus_cb_cens = mean_na(dv),
             stud_diff_se = if (length(dv) >= 2) {
               stats::sd(dv) / sqrt(length(dv))
             } else NA_real_)
}
res <- do.call(rbind, lapply(strategies, detection_rows))
rownames(res) <- NULL

attr_row <- data.frame(c_rrb = multipliers$all[["rrb"]],
                       c_cb = multipliers$all[["cb"]])
for (g in setdiff(strategies, "all")) {
  attr_row[[paste0("c_rrb_", g)]] <- multipliers[[g]][["rrb"]]
  attr_row[[paste0("c_cb_", g)]] <- multipliers[[g]][["cb"]]
}
attr_row$c_stud_rrb <- multipliers_stud$all[["rrb"]]
attr_row$c_cb_cens_stud <- multipliers_stud$all[["cb"]]
for (g in setdiff(strategies, "all")) {
  attr_row[[paste0("c_stud_rrb_", g)]] <- multipliers_stud[[g]][["rrb"]]
  attr_row[[paste0("c_cb_cens_stud_", g)]] <- multipliers_stud[[g]][["cb"]]
}
write.csv(res, file.path(out_dir, "sim_part1_detection_size_adjusted.csv"),
          row.names = FALSE)
write.csv(attr_row,
          file.path(out_dir, "sim_part1_detection_size_multipliers.csv"),
          row.names = FALSE)

cat("Multipliers c (obs > c * q95 gives FPR 0.05 under P1-01):\n")
for (g in strategies) {
  cat(sprintf("  %-3s  RRB %.3f, CB %.3f\n", g, multipliers[[g]][["rrb"]],
              multipliers[[g]][["cb"]]))
}
cat("Multipliers c, studentised (RRB; CB on the censorable rows):\n")
for (g in strategies) {
  cat(sprintf("  %-3s  RRB %.3f, CB_cens %.3f\n", g,
              multipliers_stud[[g]][["rrb"]],
              multipliers_stud[[g]][["cb"]]))
}
show_cols <- c("id", "misspec", "cens_strength", "outcome_strength", "n_rep",
               "nominal_rrb", "nominal_cb", "adj_rrb", "adj_cb",
               "diff_rrb_minus_cb", "diff_se", "nominal_stud_rrb",
               "nominal_cb_cens_stud", "stud_adj_rrb", "cb_cens_stud_adj",
               "stud_diff_rrb_minus_cb_cens", "stud_diff_se")
for (g in strategies) {
  cat(sprintf("\nStrategy %s%s:\n", g,
              if (g == "all") " (pooled over strategies)" else ""))
  print(res[res$strategy == g, show_cols], digits = 3, row.names = FALSE)
}
