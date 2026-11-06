# 02b_simulation_tables_paper.R -----------------------------------------------
# Tables and figures of the simulation study of the manuscript (Parts 1-3
# and the Supporting Information). One row per scenario id of paper_scenarios()
# (analysis/sim_scenarios.R; P1-01 to P3-04 in the row order of Table
# tab:grid), each taken from the run that supplies it (main, manycov,
# largerep, variants; analysis/sim_paper_results.R): cell means with Monte
# Carlo (MC) SEs, detection / false-positive rates against the conditional
# parametric null reference (rrb_null()) of interval RRB and of censoring
# balance (CB), and within-scenario Spearman correlations (secondary). No
# pooled Spearman / AUROC.
#
# Strategies: g1 initiates in intervals 0-2, g2 in intervals 0-8 (nested
# longer window), g3 does not initiate through interval 8. Two contrasts
# (suffix <cc>): g2g1, theta = risk(g2) - risk(g1); g3g2, theta = risk(g3) -
# risk(g2). The pooled diagnostics are aggregated over g1, g2 and g3;
# end-of-grace SMDs are per contrast (g2g1 at interval 3, g3g2 at interval
# 9). The censoring model is fitted separately per strategy (g1 and g2: one
# censoring interval each; g3: pooled logistic model over intervals 0-8), so
# the per-strategy summaries (suffix _<g>, g = g1, g2, g3) judge each
# strategy against the null reference of its own model.
#
# Inputs:  <sim_dir>/sim_results_<run>.csv of the runs main, manycov,
#          largerep and variants
# Outputs: analysis/scenarios_paper.csv        paper_scenarios(): id, label,
#                                              run and spec (default
#                                              sim_dir only)
#          <sim_dir>/sim_paper_scenarios.csv   one row per id: performance
#                                              per contrast (theta_true_<cc>,
#                                              bias_<cc>, mc_se_bias_<cc>,
#                                              mc_sd_<cc>,
#                                              abs_bias_over_mc_sd_<cc>,
#                                              rmse_<cc>, coverage_<cc>,
#                                              mcse_coverage_<cc>,
#                                              se_ratio_<cc>), diagnostic
#                                              means and MC SEs,
#                                              null-reference exceedance
#                                              (RRB overall and per history
#                                              function, CB), bootstrap-
#                                              interval exclusion rates of
#                                              cumulative RRB of X (g1, g2,
#                                              g3 x eog1, eog2, h);
#                                              per strategy g:
#                                              rrb_all_maxib_<g>_mean/_mcse,
#                                              exceed_rate_<g>(_se),
#                                              null_maxib_obs_<g>_mean
#                                              (model weights),
#                                              null_q95_<g>_mean, mc_q95_<g>
#                                              (95th pct over replicates of
#                                              null_maxib_obs_<g>),
#                                              jackson_cb_max_<g>_mean/_mcse,
#                                              null_cb_exceed_rate_<g>,
#                                              exceed_rate_<f>_<g>;
#                                              appended after these, pooled
#                                              (no suffix) and per strategy
#                                              (suffix _<g>): studentised
#                                              maximum (the exceedance
#                                              standard) stud_exceed_rate(_se),
#                                              stud_exceed_rate_<f>,
#                                              null_z_max_<g>_mean/_mcse,
#                                              null_z_q95_mean, mc_z_q95 (95th
#                                              pct over replicates of
#                                              null_z_max),
#                                              null_z_n_excluded_mean;
#                                              CB on the censorable rows:
#                                              jackson_cb_cens_max_mean/_mcse,
#                                              null_cb_cens_exceed_rate(_se),
#                                              null_cb_cens_maxib_obs_mean,
#                                              null_cb_cens_q95_mean,
#                                              null_cb_cens_stud_exceed_rate
#                                              (_se), null_cb_cens_z_max_mean,
#                                              null_cb_cens_z_q95_mean,
#                                              null_cb_cens_mc_z_q95,
#                                              null_cb_cens_exceed_rate_<f>
#                                              (raw only; no studentised CB
#                                              per feature)
#          <sim_dir>/sim_paper_calibration.csv  id x strategy ("all" =
#                                              pooled over g1-g3, or g1, g2,
#                                              g3) x diagnostic (RRB, CB) x
#                                              history function: exceedance
#                                              rate of the null q95, exact
#                                              95% CI (filter strategy ==
#                                              "all" for the pooled rows);
#                                              column statistic: "raw" (MaxIB,
#                                              rows of RRB and CB first),
#                                              "studentised" (RRB overall and
#                                              per function; CB_cens overall),
#                                              then "raw" rows of CB_cens (CB
#                                              on the censorable rows)
#          <sim_dir>/sim_paper_within_spearman.csv  diagnostic x contrast:
#                                              diagnostic vs |bias_<cc>|
#                                              within scenario, Fisher-z SE
#          <sim_dir>/sim_part1_factorial.csv   Part 1 2 x 2 cell means and
#                                              interaction contrasts
#                                              (quantities bias_g2g1,
#                                              bias_g3g2, the primary
#                                              diagnostics and
#                                              rrb_all_maxib_<g>,
#                                              jackson_cb_max_<g>,
#                                              jackson_cb_cens_max_<g>,
#                                              null_z_max_<g>)
#          <sim_dir>/sim_part1_detection_rrb_vs_cb.csv  Part 1: exceedance
#                                              of the null q95 by interval-RRB
#                                              MaxIB vs CB MaxIB, paired
#                                              difference and SEs
#          <fig_dir>/sim_part1_factorial.png   Part 1 factorial plot
#                                              (omitted term x {bias of
#                                              g3 - g2, interval-RRB MaxIB
#                                              of g1, g2, g3})
#          <fig_dir>/sim_part2_calibration.png  Part 2 null-reference
#                                              calibration (base mechanism)
#                                              per strategy g1, g2, g3,
#                                              overall and per history
#                                              function: exceedance of the
#                                              studentised maximum
#          <fig_dir>/sim_part2_calibration_raw.png  the same for the raw
#                                              interval-RRB MaxIB (Web
#                                              Appendix)
#
# Environment: CCW_SIM_DIR (default "output/sim") is the directory of the
# results (sim_dir) and of the CSV outputs. Figures go to output/figures
# (fig_dir) for the default and to <CCW_SIM_DIR>/figures otherwise;
# analysis/scenarios_paper.csv is written only for the default.
#
# Run from the project root:  Rscript analysis/02b_simulation_tables_paper.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))
library(ggplot2)

default_sim_dir <- file.path("output", "sim")
sim_dir <- Sys.getenv("CCW_SIM_DIR", "")
if (!nzchar(sim_dir)) sim_dir <- default_sim_dir
is_default_dir <- identical(normalizePath(sim_dir, mustWork = FALSE),
                            normalizePath(default_sim_dir, mustWork = FALSE))
out_dir <- sim_dir
fig_dir <- if (is_default_dir) file.path("output", "figures") else
  file.path(sim_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x) file.path(out_dir, sprintf("sim_%s.csv", x))
out_png <- function(x) file.path(fig_dir, sprintf("sim_%s.png", x))

scen <- paper_scenarios()
scen_csv <- file.path("analysis", "scenarios_paper.csv")
if (is_default_dir) write.csv(scen, scen_csv, row.names = FALSE)
res <- read_paper_results(sim_dir = sim_dir)

# contrasts: g2g1 = risk(g2) - risk(g1), g3g2 = risk(g3) - risk(g2)
contrasts <- c("g2g1", "g3g2")
primary <- c("rrb_all_maxib", "ks_max", "mhb_max", "jackson_cb_max",
             "cumrrb_base_maxib", "eog_smd_base_max_g2g1",
             "eog_smd_base_max_g3g2")
strategies <- c("g1", "g2", "g3")
# history functions: null_exceed_<f> without the per-strategy suffix _<g>
feats <- grep("^null_exceed_", names(res), value = TRUE)
feats <- sub("^null_exceed_", "",
             grep("_g[1-3]$", feats, value = TRUE, invert = TRUE))
labels_math <- c(
  bias_g2g1 = "'Bias ('*g[2]*' vs. '*g[1]*')'",
  bias_g3g2 = "'Bias ('*g[3]*' vs. '*g[2]*')'",
  rrb_all_maxib = "'Interval RRB (all history functions)'",
  rrb_all_maxib_g1 = "'MaxIB, '*g[1]",
  rrb_all_maxib_g2 = "'MaxIB, '*g[2]",
  rrb_all_maxib_g3 = "'MaxIB, '*g[3]",
  ks_max = "KS^RR", mhb_max = "M^RR",
  jackson_cb_max = "'Censoring balance'",
  cumrrb_base_maxib = "'Cumulative RRB (baseline)'",
  eog_smd_base_max_g2g1 = "'End-of-grace SMD of X ('*g[2]*' vs. '*g[1]*')'",
  eog_smd_base_max_g3g2 = "'End-of-grace SMD of X ('*g[3]*' vs. '*g[2]*')'"
)
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
mcse <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2) NA_real_ else stats::sd(v) / sqrt(length(v))
}
# mean of a column; NA when the column is absent or all NA
col_mean <- function(x, nm) if (nm %in% names(x)) mean_na(x[[nm]]) else
  NA_real_
# a column, or NA when it is absent
col_get <- function(x, nm) if (nm %in% names(x)) x[[nm]] else
  rep(NA_real_, nrow(x))
q95_na <- function(v) {
  if (all(is.na(v))) NA_real_ else
    stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE)
}
sd_na <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2) NA_real_ else stats::sd(v)
}
# binomial rate of a 0/1 column and its SE; NA when absent or all NA
rate_se <- function(x, nm) {
  v <- col_get(x, nm)
  v <- v[!is.na(v)]
  if (length(v) == 0) return(c(NA_real_, NA_real_))
  c(mean(v), sqrt(mean(v) * (1 - mean(v)) / length(v)))
}
cum_cells <- as.vector(outer(c("g1", "g2", "g3"), c("eog1", "eog2", "h"),
                             paste, sep = "_"))
key_cols <- c("id", "part", "label", "run", "scenario")

# ---------------------------------------------------------------------------
# 1. per-scenario summary
summ <- do.call(rbind, lapply(scen$id, function(s) {
  x <- res[res$id == s, ]
  n_rep <- nrow(x)
  exc <- mean_na(x$null_exceed)
  n_exc <- sum(!is.na(x$null_exceed))
  out <- data.frame(id = s, n_rep = n_rep)
  for (cc in contrasts) {
    th <- x[[paste0("theta_hat_", cc)]]
    b <- x[[paste0("bias_", cc)]]
    mc_sd <- sd_na(th)
    n_th <- sum(!is.na(th))
    cov <- mean_na(x[[paste0("covered_", cc)]])
    n_cov <- sum(!is.na(x[[paste0("covered_", cc)]]))
    out[[paste0("theta_true_", cc)]] <- x[[paste0("theta_true_", cc)]][1]
    out[[paste0("bias_", cc)]] <- mean_na(b)
    out[[paste0("mc_se_bias_", cc)]] <- mc_sd / sqrt(n_th)
    out[[paste0("mc_sd_", cc)]] <- mc_sd
    out[[paste0("abs_bias_over_mc_sd_", cc)]] <- abs(mean_na(b)) / mc_sd
    out[[paste0("rmse_", cc)]] <- sqrt(mean_na(b^2))
    out[[paste0("coverage_", cc)]] <- cov
    out[[paste0("mcse_coverage_", cc)]] <- sqrt(cov * (1 - cov) / n_cov)
    out[[paste0("se_ratio_", cc)]] <-
      mean_na(x[[paste0("se_boot_", cc)]]) / mc_sd
  }
  for (dg in primary) {
    out[[paste0(dg, "_mean")]] <- mean_na(x[[dg]])
    out[[paste0(dg, "_mcse")]] <- mcse(x[[dg]])
  }
  out$jackson_cb_est_max_mean <- col_mean(x, "jackson_cb_est_max")
  out$ess_ratio_min_mean <- col_mean(x, "ess_ratio_min")
  out$w99_max_mean <- col_mean(x, "w99_max")
  out$pi_min_mean <- col_mean(x, "pi_min")
  for (cc in contrasts) {
    out[[sprintf("strat_smd_L_max_%s_mean", cc)]] <-
      col_mean(x, paste0("strat_smd_L_max_", cc))
    out[[sprintf("strat_smd_X_max_%s_mean", cc)]] <-
      col_mean(x, paste0("strat_smd_X_max_", cc))
  }
  out$exceed_rate <- exc
  out$exceed_rate_se <- if (n_exc > 0) sqrt(exc * (1 - exc) / n_exc) else
    NA_real_
  out$null_q95_mean <- col_mean(x, "null_q95")
  out$mc_q95 <- q95_na(x$null_maxib_obs)
  out$null_maxib_obs_mean <- col_mean(x, "null_maxib_obs")
  out$null_z_max_mean <- col_mean(x, "null_z_max")
  out$null_frac_cells_p05_mean <- col_mean(x, "null_frac_cells_p05")
  out$null_n_cells_mean <- col_mean(x, "null_n_cells")
  for (f in feats) {
    out[[paste0("exceed_rate_", f)]] <- mean_na(x[[paste0("null_exceed_", f)]])
  }
  out$ltfu_cb_max_mean <- col_mean(x, "ltfu_cb_max")
  out$ltfu_rrb_max_mean <- col_mean(x, "ltfu_rrb_max")
  # null reference of censoring balance (same draws)
  e_cb <- if ("null_cb_exceed" %in% names(x)) x$null_cb_exceed else NA_real_
  n_cb <- sum(!is.na(e_cb))
  out$null_cb_exceed_rate <- mean_na(e_cb)
  out$null_cb_exceed_rate_se <- if (n_cb > 0) {
    sqrt(out$null_cb_exceed_rate * (1 - out$null_cb_exceed_rate) / n_cb)
  } else NA_real_
  out$null_cb_q95_mean <- col_mean(x, "null_cb_q95")
  out$null_cb_mc_q95 <- if ("null_cb_maxib_obs" %in% names(x)) {
    q95_na(x$null_cb_maxib_obs)
  } else NA_real_
  out$null_cb_maxib_obs_mean <- col_mean(x, "null_cb_maxib_obs")
  # bootstrap interval of the cumulative RRB of X
  for (cc in cum_cells) {
    out[[paste0("cumrrb_X_", cc, "_obs_mean")]] <-
      col_mean(x, sprintf("cumrrb_X_%s_obs", cc))
    out[[paste0("cumrrb_X_", cc, "_boot_se_mean")]] <-
      col_mean(x, sprintf("cumrrb_X_%s_boot_se", cc))
    out[[paste0("cumrrb_boot_excl0_rate_", cc)]] <-
      col_mean(x, sprintf("cumrrb_X_%s_excl0", cc))
  }
  out$cumrrb_boot_excl0_rate_any <- col_mean(x, "cumrrb_boot_excl0_any")
  out$cumrrb_base_eog1_mean <- col_mean(x, "cumrrb_base_eog1")
  out$cumrrb_base_eog2_mean <- col_mean(x, "cumrrb_base_eog2")
  for (cc in contrasts) {
    out[[sprintf("cumrrb_base_eog_diff_%s_mean", cc)]] <-
      col_mean(x, paste0("cumrrb_base_eog_diff_", cc))
  }
  # per strategy, each against the null reference of its own model
  for (g in strategies) {
    nm <- function(stem, end = "") paste0(stem, "_", g, end)
    rrb <- col_get(x, nm("rrb_all_maxib"))
    cbm <- col_get(x, nm("jackson_cb_max"))
    e_g <- col_get(x, nm("null_exceed"))
    r_g <- mean_na(e_g)
    n_g <- sum(!is.na(e_g))
    out[[nm("rrb_all_maxib", "_mean")]] <- mean_na(rrb)
    out[[nm("rrb_all_maxib", "_mcse")]] <- mcse(rrb)
    out[[nm("exceed_rate")]] <- r_g
    out[[nm("exceed_rate", "_se")]] <- if (n_g > 0) {
      sqrt(r_g * (1 - r_g) / n_g)
    } else NA_real_
    out[[nm("null_maxib_obs", "_mean")]] <- col_mean(x, nm("null_maxib_obs"))
    out[[nm("null_q95", "_mean")]] <- col_mean(x, nm("null_q95"))
    out[[nm("mc_q95")]] <- q95_na(col_get(x, nm("null_maxib_obs")))
    out[[nm("jackson_cb_max", "_mean")]] <- mean_na(cbm)
    out[[nm("jackson_cb_max", "_mcse")]] <- mcse(cbm)
    out[[nm("null_cb_exceed_rate")]] <- col_mean(x, nm("null_cb_exceed"))
    for (f in feats) {
      out[[sprintf("exceed_rate_%s_%s", f, g)]] <-
        col_mean(x, sprintf("null_exceed_%s_%s", f, g))
    }
  }
  # studentised maximum (the exceedance standard) and censoring balance
  # restricted to the censorable rows; appended after the columns above,
  # pooled (no suffix) and per strategy (suffix _<g>)
  for (g in c("all", strategies)) {
    sfx <- if (g == "all") "" else paste0("_", g)
    nm <- function(stem, end = "") paste0(stem, sfx, end)
    zr <- rate_se(x, nm("null_z_exceed"))
    out[[nm("stud_exceed_rate")]] <- zr[1]
    out[[nm("stud_exceed_rate", "_se")]] <- zr[2]
    for (f in feats) {
      out[[paste0("stud_exceed_rate_", f, sfx)]] <-
        col_mean(x, paste0("null_z_exceed_", f, sfx))
    }
    zmax <- col_get(x, nm("null_z_max"))
    if (g != "all") {
      # null_z_max_mean (pooled) is already written above
      out[[nm("null_z_max", "_mean")]] <- mean_na(zmax)
      out[[nm("null_z_max", "_mcse")]] <- mcse(zmax)
    }
    out[[nm("null_z_q95", "_mean")]] <- col_mean(x, nm("null_z_q95"))
    out[[nm("mc_z_q95")]] <- q95_na(zmax)
    out[[nm("null_z_n_excluded", "_mean")]] <-
      col_mean(x, nm("null_z_n_excluded"))
    cbc <- col_get(x, nm("jackson_cb_cens_max"))
    out[[nm("jackson_cb_cens_max", "_mean")]] <- mean_na(cbc)
    out[[nm("jackson_cb_cens_max", "_mcse")]] <- mcse(cbc)
    cr <- rate_se(x, nm("null_cb_cens_exceed"))
    out[[nm("null_cb_cens_exceed_rate")]] <- cr[1]
    out[[nm("null_cb_cens_exceed_rate", "_se")]] <- cr[2]
    out[[nm("null_cb_cens_maxib_obs", "_mean")]] <-
      col_mean(x, nm("null_cb_cens_maxib_obs"))
    out[[nm("null_cb_cens_q95", "_mean")]] <-
      col_mean(x, nm("null_cb_cens_q95"))
    czr <- rate_se(x, nm("null_cb_cens_z_exceed"))
    out[[nm("null_cb_cens_stud_exceed_rate")]] <- czr[1]
    out[[nm("null_cb_cens_stud_exceed_rate", "_se")]] <- czr[2]
    czmax <- col_get(x, nm("null_cb_cens_z_max"))
    out[[nm("null_cb_cens_z_max", "_mean")]] <- mean_na(czmax)
    out[[nm("null_cb_cens_z_q95", "_mean")]] <-
      col_mean(x, nm("null_cb_cens_z_q95"))
    out[[nm("null_cb_cens_mc_z_q95")]] <- q95_na(czmax)
    # per feature only the raw censorable-CB exceedance exists
    for (f in feats) {
      out[[paste0("null_cb_cens_exceed_rate_", f, sfx)]] <-
        col_mean(x, paste0("null_cb_cens_exceed_", f, sfx))
    }
  }
  out
}))
summ <- merge(scen[, c(key_cols, "misspec", "cens_strength",
                       "outcome_strength", "n", "n_boot", "n_extra",
                       "truncate")],
              summ, by = "id", sort = FALSE)
summ <- summ[match(scen$id, summ$id), ]
rownames(summ) <- NULL
# LTFU summaries are only meaningful where there is loss to follow-up
no_ltfu <- summ$id != "P3-03"
summ$ltfu_cb_max_mean[no_ltfu] <- NA
summ$ltfu_rrb_max_mean[no_ltfu] <- NA
write.csv(summ, out_csv("paper_scenarios"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. calibration / detection: exceedance rate of the null-reference q95 for
#    interval-RRB MaxIB and CB MaxIB, overall and per history function, with
#    exact binomial 95% CI; pooled over g1-g3 (strategy "all") and per
#    strategy (columns with suffix _<g>)
#    Column statistic: "raw" = MaxIB on the SMD scale (the rows of the
#    earlier versions, first), "studentised" = studentised maximum (the
#    exceedance standard; RRB overall and per history function, diagnostic
#    CB_cens = CB on the censorable rows, overall only). The raw
#    censorable-CB rows (diagnostic CB_cens, statistic raw) follow the
#    studentised rows.
cal_cols_raw <- do.call(rbind, lapply(c("all", strategies), function(g) {
  sfx <- if (g == "all") "" else paste0("_", g)
  rbind(
    data.frame(strategy = g, diagnostic = "RRB", feature = c("all", feats),
               col = paste0(c("null_exceed", paste0("null_exceed_", feats)),
                            sfx)),
    data.frame(strategy = g, diagnostic = "CB", feature = c("all", feats),
               col = paste0(c("null_cb_exceed",
                              paste0("null_cb_exceed_", feats)), sfx))
  )
}))
cal_cols_raw$statistic <- "raw"
cal_cols_stud <- do.call(rbind, lapply(c("all", strategies), function(g) {
  sfx <- if (g == "all") "" else paste0("_", g)
  rbind(
    data.frame(strategy = g, diagnostic = "RRB", feature = c("all", feats),
               col = paste0(c("null_z_exceed",
                              paste0("null_z_exceed_", feats)), sfx)),
    data.frame(strategy = g, diagnostic = "CB_cens", feature = "all",
               col = paste0("null_cb_cens_z_exceed", sfx))
  )
}))
cal_cols_stud$statistic <- "studentised"
cal_cols_cbc <- do.call(rbind, lapply(c("all", strategies), function(g) {
  sfx <- if (g == "all") "" else paste0("_", g)
  data.frame(strategy = g, diagnostic = "CB_cens",
             feature = c("all", feats),
             col = paste0(c("null_cb_cens_exceed",
                            paste0("null_cb_cens_exceed_", feats)), sfx))
}))
cal_cols_cbc$statistic <- "raw"
cal_table <- function(cols) {
  out <- do.call(rbind, lapply(scen$id, function(s) {
    x <- res[res$id == s, ]
    do.call(rbind, lapply(seq_len(nrow(cols)), function(k) {
      if (!cols$col[k] %in% names(x)) return(NULL)
      v <- x[[cols$col[k]]]
      v <- v[!is.na(v)]
      if (length(v) == 0) return(NULL)
      ci <- stats::binom.test(sum(v), length(v))$conf.int
      data.frame(id = s, strategy = cols$strategy[k],
                 diagnostic = cols$diagnostic[k],
                 feature = cols$feature[k], n_rep = length(v),
                 exceed_rate = mean(v), ci_lower = ci[1], ci_upper = ci[2],
                 statistic = cols$statistic[k])
    }))
  }))
  if (is.null(out)) return(NULL)
  cbind(scen[match(out$id, scen$id), key_cols], out[, -1],
        row.names = NULL)
}
cal <- rbind(cal_table(cal_cols_raw), cal_table(cal_cols_stud),
             cal_table(cal_cols_cbc))
write.csv(cal, out_csv("paper_calibration"), row.names = FALSE)

# 2b. Part 1 detection: interval-RRB MaxIB vs censoring-balance MaxIB, each
#     against its own conditional null q95 (same null draws). Paired by
#     replicate: SE of the mean of d_i = exceed_RRB_i - exceed_CB_i.
bin_se <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) NA_real_ else sqrt(mean(v) * (1 - mean(v)) / length(v))
}
p1_ids <- scen$id[scen$part == 1]
det <- do.call(rbind, lapply(p1_ids, function(s) {
  x <- res[res$id == s, ]
  e_r <- x$null_exceed
  e_c <- if ("null_cb_exceed" %in% names(x)) x$null_cb_exceed else
    rep(NA_real_, nrow(x))
  ok <- !is.na(e_r) & !is.na(e_c)
  dd <- e_r[ok] - e_c[ok]
  data.frame(
    id = s, n_rep = nrow(x), n_paired = sum(ok),
    exceed_rate = mean_na(e_r), exceed_rate_se = bin_se(e_r),
    null_cb_exceed_rate = mean_na(e_c), null_cb_exceed_rate_se = bin_se(e_c),
    diff_rrb_minus_cb = if (any(ok)) mean(dd) else NA_real_,
    diff_se = if (sum(ok) >= 2) stats::sd(dd) / sqrt(sum(ok)) else NA_real_
  )
}))
det <- merge(scen[, c(key_cols, "misspec", "cens_strength",
                      "outcome_strength")], det, by = "id", sort = FALSE)
det <- det[match(p1_ids, det$id), ]
write.csv(det, out_csv("part1_detection_rrb_vs_cb"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. within-scenario Spearman (secondary): diagnostic x contrast, against
#    |bias_<cc>|
spear <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
wsp <- do.call(rbind, lapply(scen$id, function(s) {
  x <- res[res$id == s, ]
  grid <- expand.grid(diagnostic = primary, contrast = contrasts,
                      stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
    dg <- grid$diagnostic[k]
    cc <- grid$contrast[k]
    ab <- x[[paste0("abs_bias_", cc)]]
    ok <- !is.na(x[[dg]]) & !is.na(ab)
    n_ok <- sum(ok)
    rho <- spear(x[[dg]], ab)
    se_z <- if (n_ok > 3) 1 / sqrt(n_ok - 3) else NA_real_
    data.frame(
      id = s, diagnostic = dg, contrast = cc, n_rep = n_ok, spearman = rho,
      se_fisher_z = se_z, se = (1 - rho^2) * se_z,
      ci_lower = tanh(atanh(rho) - 1.96 * se_z),
      ci_upper = tanh(atanh(rho) + 1.96 * se_z)
    )
  }))
}))
wsp <- cbind(scen[match(wsp$id, scen$id), key_cols], wsp[, -1],
             row.names = NULL)
write.csv(wsp, out_csv("paper_within_spearman"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Part 1 factorial: 2 x 2 cell means and interaction contrasts
#    interaction = (strong cens, strong outcome) - (strong, none)
#                  - (weak, strong) + (weak, none)
s1 <- summ[summ$part == 1 & summ$misspec != "none", ]
bias_q <- paste0("bias_", contrasts)
quant <- c(bias_q, primary, paste0("rrb_all_maxib_", strategies),
           paste0("jackson_cb_max_", strategies),
           paste0("jackson_cb_cens_max_", strategies),
           paste0("null_z_max_", strategies))
# cell mean and its MC SE; NA when the cell or the column is absent
cell_mean <- function(q, row) {
  get1 <- function(nm) {
    if (nrow(row) == 0 || !nm %in% names(row)) NA_real_ else
      as.numeric(row[[nm]][1])
  }
  if (q %in% bias_q) {
    cc <- sub("^bias_", "", q)
    c(get1(q), get1(paste0("mc_se_bias_", cc)))
  } else {
    c(get1(paste0(q, "_mean")), get1(paste0(q, "_mcse")))
  }
}
fact <- do.call(rbind, lapply(unique(s1$misspec), function(ms) {
  sub <- s1[s1$misspec == ms, ]
  get_row <- function(cs, os) {
    sub[sub$cens_strength == cs & sub$outcome_strength == os, ]
  }
  do.call(rbind, lapply(quant, function(q) {
    cells <- list(wn = cell_mean(q, get_row("weak", "none")),
                  ws = cell_mean(q, get_row("weak", "strong")),
                  sn = cell_mean(q, get_row("strong", "none")),
                  ss = cell_mean(q, get_row("strong", "strong")))
    m <- vapply(cells, `[`, 0, 1)
    se <- vapply(cells, `[`, 0, 2)
    data.frame(
      misspec = ms,
      ids = paste(sub$id[order(sub$id)], collapse = ","),
      quantity = q,
      mean_cw_on = m[["wn"]], mean_cw_os = m[["ws"]],
      mean_cs_on = m[["sn"]], mean_cs_os = m[["ss"]],
      mcse_cw_on = se[["wn"]], mcse_cw_os = se[["ws"]],
      mcse_cs_on = se[["sn"]], mcse_cs_os = se[["ss"]],
      main_cens = (m[["sn"]] + m[["ss"]] - m[["wn"]] - m[["ws"]]) / 2,
      main_outcome = (m[["ws"]] + m[["ss"]] - m[["wn"]] - m[["sn"]]) / 2,
      interaction = m[["ss"]] - m[["sn"]] - m[["ws"]] + m[["wn"]],
      interaction_mcse = sqrt(sum(se^2))
    )
  }))
}))
write.csv(fact, out_csv("part1_factorial"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. figures
misspec_lab <- c(omit_L = "'Omitted '*L[k]",
                 omit_Lsq = "'Omitted '*L[k]^2",
                 omit_lag = "'Omitted '*L[k-1]")
# the figure shows the bias of g3 - g2 (where the bias concentrates) and
# interval RRB per strategy (each strategy has its own censoring model); the
# bias of g2 - g1, pooled RRB, censoring balance, KS^RR and M^RR are kept in
# the CSV outputs (Table 2 / web appendix) but not plotted
fq <- c("bias_g3g2", paste0("rrb_all_maxib_", strategies))
fl <- do.call(rbind, lapply(fq, function(q) {
  cm <- vapply(seq_len(nrow(s1)), function(k) cell_mean(q, s1[k, ]),
               numeric(2))
  data.frame(misspec = s1$misspec, cens = s1$cens_strength,
             outcome = s1$outcome_strength, quantity = q,
             mean = cm[1, ], mcse = cm[2, ])
}))
fl$quantity <- factor(labels_math[fl$quantity], levels = labels_math[fq])
fl$misspec <- factor(misspec_lab[fl$misspec], levels = misspec_lab)
fl$cens <- factor(fl$cens, levels = c("weak", "strong"))
fl$outcome <- factor(fl$outcome, levels = c("none", "strong"))
ref_correct <- summ[summ$id == "P1-01", ]
ref <- data.frame(
  quantity = factor(labels_math[fq], levels = labels_math[fq]),
  value = vapply(fq, function(q) cell_mean(q, ref_correct)[1], 0)
)
# outcome effect encoded by symbol, line type and grey level (readable in
# monochrome and for colour-vision deficiency); the reference is dotted
outcome_lab <- "Outcome strength"
p_fact <- ggplot(fl, aes(cens, mean, colour = outcome, shape = outcome,
                         linetype = outcome, group = outcome)) +
  geom_hline(data = ref, aes(yintercept = value), linetype = "dotted",
             colour = "grey40") +
  geom_line(position = position_dodge(width = 0.2)) +
  geom_pointrange(aes(ymin = mean - 1.96 * mcse, ymax = mean + 1.96 * mcse),
                  position = position_dodge(width = 0.2), size = 0.45,
                  linetype = "solid", fill = "white") +
  # free scales: the bias column is an order of magnitude below RRB, and
  # the interaction (bias only with both factors strong) must stay visible
  facet_grid(misspec ~ quantity, scales = "free",
             labeller = label_parsed) +
  scale_colour_manual(values = c(none = "grey45", strong = "black"),
                      name = outcome_lab) +
  scale_shape_manual(values = c(none = 21, strong = 17),
                     name = outcome_lab) +
  scale_linetype_manual(values = c(none = "dashed", strong = "solid"),
                        name = outcome_lab) +
  labs(x = "Censoring strength",
       y = "Mean over replicates (95% MC interval)",
       caption = "Dotted line: correctly specified scenario.") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(out_png("part1_factorial"), p_fact, width = 9.5, height = 7,
       dpi = 150)

# Part 2 calibration (base mechanism): exceedance rate of the null-reference
# q95 of interval RRB per strategy (each against the null reference of its
# own censoring model), overall and per history function (P1-01, the
# n = 2000 correct scenario, is the n = 2000 member of Part 2); the pooled
# rows are in the CSV only
s2_lab <- c(`P2-01` = "Correct specification, n = 300",
            `P1-01` = "Correct specification, n = 2,000",
            `P2-02` = "Correct specification, n = 8,000",
            `P2-03` = "Near-violation of positivity",
            `P2-04` = "Weight truncation, 95th percentile")
s2_ids <- names(s2_lab)
# strategies by shape and Okabe-Ito colour; the legends of colour, fill
# and shape merge (same name and labels)
strat_col <- c(g1 = "#0072B2", g2 = "#D55E00", g3 = "#009E73")
strat_lab <- c(expression(g[1]), expression(g[2]), expression(g[3]))
dodge_y <- position_dodge(width = 0.6, orientation = "y", reverse = TRUE)
# one calibration plot for statistic "studentised" (Fig 3, the standard)
# or "raw" (interval-RRB MaxIB on the SMD scale; Web Appendix)
plot_calibration <- function(stat, xlab) {
  calp <- cal[cal$diagnostic == "RRB" & cal$statistic == stat &
                cal$id %in% s2_ids & cal$strategy %in% strategies, ]
  if (nrow(calp) == 0) {
    stop("No ", stat, " per-strategy exceedance columns in the results ",
         "(null_z_exceed_<g> needs the 02_simulation_study.R of 8a2bd6a).")
  }
  calp <- calp[order(match(calp$id, s2_ids)), ]
  calp$strategy <- factor(calp$strategy, levels = strategies)
  calp$feature <- factor(ifelse(calp$feature == "all",
                                "All history functions", calp$feature),
                         levels = c("All history functions", feats))
  calp$scenario <- factor(s2_lab[calp$id], levels = s2_lab[s2_ids])
  ggplot(calp, aes(exceed_rate, scenario, colour = strategy,
                   fill = strategy, shape = strategy)) +
    geom_vline(xintercept = 0.05, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.3,
                    position = dodge_y) +
    facet_wrap(~feature, ncol = 3) +
    scale_colour_manual(values = strat_col, labels = strat_lab,
                        name = "Strategy") +
    scale_fill_manual(values = strat_col, labels = strat_lab,
                      name = "Strategy") +
    scale_shape_manual(values = c(g1 = 21, g2 = 24, g3 = 22),
                       labels = strat_lab, name = "Strategy") +
    labs(x = xlab, y = NULL, caption = "Dashed line: nominal 0.05.") +
    theme_bw() + theme(legend.position = "bottom")
}
p_cal <- plot_calibration(
  "studentised",
  paste("Exceedance rate of the studentised maximum of interval RRB",
        "over its null 95th percentile (95% exact CI)")
)
ggsave(out_png("part2_calibration"), p_cal, width = 11, height = 7.5,
       dpi = 150)
p_cal_raw <- plot_calibration(
  "raw",
  paste("Rate of interval-RRB MaxIB above the null-reference",
        "95th percentile of the strategy's model (95% exact CI)")
)
ggsave(out_png("part2_calibration_raw"), p_cal_raw, width = 11,
       height = 7.5, dpi = 150)

perf_cols <- as.vector(t(outer(
  c("bias", "mc_sd", "abs_bias_over_mc_sd", "coverage", "se_ratio"),
  contrasts, paste, sep = "_"
)))
print(summ[, c("id", "label", "n_rep", perf_cols)], digits = 3,
      row.names = FALSE)
print(summ[, c("id", "rrb_all_maxib_mean", "exceed_rate", "null_q95_mean",
               "mc_q95")], digits = 3, row.names = FALSE)
cat("Per-strategy null-reference exceedance (interval RRB, CB):\n")
print(summ[, c("id", paste0("exceed_rate_", strategies),
               paste0("null_cb_exceed_rate_", strategies))], digits = 3,
      row.names = FALSE)
cat("Per-strategy studentised exceedance (interval RRB, censorable CB):\n")
print(summ[, c("id", paste0("stud_exceed_rate_", strategies),
               paste0("null_cb_cens_stud_exceed_rate_", strategies))],
      digits = 3, row.names = FALSE)
print(det[, c("id", "label", "exceed_rate", "null_cb_exceed_rate",
              "diff_rrb_minus_cb", "diff_se")], digits = 3,
      row.names = FALSE)
cat("Written:", if (is_default_dir) scen_csv,
    out_csv("paper_scenarios"), out_csv("paper_calibration"),
    out_csv("part1_detection_rrb_vs_cb"), out_csv("paper_within_spearman"),
    out_csv("part1_factorial"), out_png("part1_factorial"),
    out_png("part2_calibration"), out_png("part2_calibration_raw"),
    sep = "\n  ")
