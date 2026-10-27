# 02b_simulation_tables_paper.R -----------------------------------------------
# Tables and figures of the simulation study of the manuscript (Parts 1-3,
# Sections 5-6). One row per scenario id of paper_scenarios()
# (analysis/sim_scenarios.R; P1-01 to P3-04 in the row order of Table
# 2), each taken from the run that supplies it (main, manycov,
# largerep, variants; analysis/sim_paper_results.R): cell means with Monte
# Carlo (MC) SEs, detection / false-positive rates against the conditional
# parametric null reference (rrb_null()) of interval RRB and of censoring
# balance (CB), and within-scenario Spearman correlations (secondary). No
# pooled Spearman / AUROC.
#
# Inputs:  output/sim/sim_results_<run>.csv of the runs main, manycov, largerep
#          and variants
# Outputs: analysis/scenarios_paper.csv        paper_scenarios(): id, label,
#                                              run and spec
#          output/sim/sim_paper_scenarios.csv  one row per id: performance,
#                                              diagnostic means and MC SEs,
#                                              null-reference exceedance
#                                              (RRB overall and per history
#                                              function, CB), bootstrap-
#                                              interval exclusion rates of
#                                              cumulative RRB of X
#          output/sim/sim_paper_calibration.csv  id x diagnostic (RRB, CB) x
#                                              history function: exceedance
#                                              rate of the null q95, exact
#                                              95% CI
#          output/sim/sim_paper_within_spearman.csv  diagnostic vs |bias|
#                                              within scenario, Fisher-z SE
#          output/sim/sim_part1_factorial.csv  Part 1 2 x 2 cell means and
#                                              interaction contrasts
#          output/sim/sim_part1_detection_rrb_vs_cb.csv  Part 1: exceedance
#                                              of the null q95 by interval-RRB
#                                              MaxIB vs CB MaxIB, paired
#                                              difference and SEs
#          output/figures/sim_part1_factorial.png  Part 1 factorial plot
#                                              (bias, interval RRB, CB)
#          output/figures/sim_part2_calibration.png  Part 2 null-reference
#                                              calibration (base mechanism)
#
# Run from the project root:  Rscript analysis/02b_simulation_tables_paper.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))
source(file.path("analysis", "sim_paper_results.R"))
library(ggplot2)

out_dir <- file.path("output", "sim")
fig_dir <- file.path("output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x) file.path(out_dir, sprintf("sim_%s.csv", x))
out_png <- function(x) file.path(fig_dir, sprintf("sim_%s.png", x))

scen <- paper_scenarios()
write.csv(scen, file.path("analysis", "scenarios_paper.csv"),
          row.names = FALSE)
res <- read_paper_results()

primary <- c("rrb_all_maxib", "ks_max", "mhb_max", "jackson_cb_max",
             "cumrrb_base_maxib", "eog_smd_base_max")
feats <- sub("^null_exceed_", "",
             grep("^null_exceed_", names(res), value = TRUE))
labels_math <- c(
  bias = "'Bias'",
  rrb_all_maxib = "'Interval RRB (all history functions)'",
  ks_max = "KS^RR", mhb_max = "M^RR",
  jackson_cb_max = "'Censoring balance'",
  cumrrb_base_maxib = "'Cumulative RRB (baseline)'",
  eog_smd_base_max = "'End-of-grace SMD of X'"
)
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
mcse <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2) NA_real_ else stats::sd(v) / sqrt(length(v))
}
# mean of a column; NA when the column is absent or all NA
col_mean <- function(x, nm) if (nm %in% names(x)) mean_na(x[[nm]]) else
  NA_real_
q95_na <- function(v) {
  if (all(is.na(v))) NA_real_ else
    stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE)
}
cum_cells <- c("g1_eog", "g3_eog", "g1_h", "g3_h")
key_cols <- c("id", "part", "label", "run", "scenario")

# ---------------------------------------------------------------------------
# 1. per-scenario summary
summ <- do.call(rbind, lapply(scen$id, function(s) {
  x <- res[res$id == s, ]
  n_rep <- nrow(x)
  mc_sd <- stats::sd(x$theta_hat)
  cov <- mean(x$covered)
  exc <- mean(x$null_exceed, na.rm = TRUE)
  n_exc <- sum(!is.na(x$null_exceed))
  out <- data.frame(
    id = s, n_rep = n_rep, theta_true = x$theta_true[1],
    bias = mean(x$bias), mc_se_bias = mc_sd / sqrt(n_rep), mc_sd = mc_sd,
    abs_bias_over_mc_sd = abs(mean(x$bias)) / mc_sd,
    rmse = sqrt(mean(x$bias^2)),
    coverage = cov, mcse_coverage = sqrt(cov * (1 - cov) / n_rep),
    se_ratio = mean(x$se_boot) / mc_sd
  )
  for (dg in primary) {
    out[[paste0(dg, "_mean")]] <- mean_na(x[[dg]])
    out[[paste0(dg, "_mcse")]] <- mcse(x[[dg]])
  }
  out$jackson_cb_est_max_mean <- col_mean(x, "jackson_cb_est_max")
  out$ess_ratio_min_mean <- col_mean(x, "ess_ratio_min")
  out$w99_max_mean <- col_mean(x, "w99_max")
  out$pi_min_mean <- col_mean(x, "pi_min")
  out$strat_smd_L_max_mean <- col_mean(x, "strat_smd_L_max")
  out$strat_smd_X_max_mean <- col_mean(x, "strat_smd_X_max")
  out$exceed_rate <- exc
  out$exceed_rate_se <- sqrt(exc * (1 - exc) / n_exc)
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
  out$cumrrb_base_eog_mean <- col_mean(x, "cumrrb_base_eog")
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
#    exact binomial 95% CI
cal <- do.call(rbind, lapply(scen$id, function(s) {
  x <- res[res$id == s, ]
  cols <- rbind(
    data.frame(diagnostic = "RRB", feature = c("all", feats),
               col = c("null_exceed", paste0("null_exceed_", feats))),
    data.frame(diagnostic = "CB", feature = c("all", feats),
               col = c("null_cb_exceed", paste0("null_cb_exceed_", feats)))
  )
  do.call(rbind, lapply(seq_len(nrow(cols)), function(k) {
    if (!cols$col[k] %in% names(x)) return(NULL)
    v <- x[[cols$col[k]]]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NULL)
    ci <- stats::binom.test(sum(v), length(v))$conf.int
    data.frame(id = s, diagnostic = cols$diagnostic[k],
               feature = cols$feature[k], n_rep = length(v),
               exceed_rate = mean(v), ci_lower = ci[1], ci_upper = ci[2])
  }))
}))
cal <- cbind(scen[match(cal$id, scen$id), key_cols], cal[, -1],
             row.names = NULL)
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
# 3. within-scenario Spearman (secondary)
spear <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman",
                                                   use = "complete.obs"))
wsp <- do.call(rbind, lapply(scen$id, function(s) {
  x <- res[res$id == s, ]
  do.call(rbind, lapply(primary, function(dg) {
    ok <- !is.na(x[[dg]])
    n_ok <- sum(ok)
    rho <- spear(x[[dg]], x$abs_bias)
    se_z <- if (n_ok > 3) 1 / sqrt(n_ok - 3) else NA_real_
    data.frame(
      id = s, diagnostic = dg, n_rep = n_ok, spearman = rho,
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
quant <- c("bias", primary)
cell_mean <- function(q, row) {
  if (q == "bias") c(row$bias, row$mc_se_bias) else
    c(row[[paste0(q, "_mean")]], row[[paste0(q, "_mcse")]])
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
misspec_lab <- c(omit_L = "'Omitted '*L[t]",
                 omit_Lsq = "'Omitted '*L[t]^2",
                 omit_lag = "'Omitted '*L[t-1]")
# KS^RR and M^RR are kept in the CSV outputs but not plotted
fig_diags <- c("rrb_all_maxib", "jackson_cb_max")
fq <- c("bias", fig_diags)
fl <- do.call(rbind, lapply(fq, function(q) {
  data.frame(misspec = s1$misspec, cens = s1$cens_strength,
             outcome = s1$outcome_strength, quantity = q,
             mean = if (q == "bias") s1$bias else s1[[paste0(q, "_mean")]],
             mcse = if (q == "bias") s1$mc_se_bias else
               s1[[paste0(q, "_mcse")]])
}))
fl$quantity <- factor(labels_math[fl$quantity], levels = labels_math[fq])
fl$misspec <- factor(misspec_lab[fl$misspec], levels = misspec_lab)
fl$cens <- factor(fl$cens, levels = c("weak", "strong"))
fl$outcome <- factor(fl$outcome, levels = c("none", "strong"))
ref_correct <- summ[summ$id == "P1-01", ]
ref <- data.frame(
  quantity = factor(labels_math[fq], levels = labels_math[fq]),
  value = c(ref_correct$bias,
            vapply(fig_diags, function(q)
              ref_correct[[paste0(q, "_mean")]], 0))
)
p_fact <- ggplot(fl, aes(cens, mean, colour = outcome, group = outcome)) +
  geom_hline(data = ref, aes(yintercept = value), linetype = 3,
             colour = "grey40") +
  geom_line(position = position_dodge(width = 0.15)) +
  geom_pointrange(aes(ymin = mean - 1.96 * mcse, ymax = mean + 1.96 * mcse),
                  position = position_dodge(width = 0.15), size = 0.3) +
  facet_grid(misspec ~ quantity, scales = "free_y",
             labeller = label_parsed) +
  scale_colour_manual(values = c(none = "#1b9e77", strong = "#d95f02")) +
  labs(x = "Strength of the omitted term in the censoring model",
       y = "Mean over replicates (95% MC interval)",
       colour = "Outcome effect of the omitted term",
       caption = "Dotted line: correctly specified scenario.") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(out_png("part1_factorial"), p_fact, width = 9, height = 7.5,
       dpi = 150)

# Part 2 calibration (base mechanism): exceedance rate of the null-reference
# q95 of interval RRB, overall and per history function (P1-01, the
# n = 2000 correct scenario, is the n = 2000 member of Part 2)
s2_lab <- c(`P2-01` = "Correct, n = 300", `P1-01` = "Correct, n = 2000",
            `P2-02` = "Correct, n = 8000",
            `P2-03` = "Near positivity violation",
            `P2-04` = "Truncation at 95th pct")
s2_ids <- names(s2_lab)
calp <- cal[cal$diagnostic == "RRB" & cal$id %in% s2_ids, ]
calp <- calp[order(match(calp$id, s2_ids)), ]
calp$feature <- factor(ifelse(calp$feature == "all", "All history functions",
                              calp$feature),
                       levels = c("All history functions", feats))
calp$scenario <- factor(s2_lab[calp$id], levels = s2_lab[s2_ids])
p_cal <- ggplot(calp, aes(exceed_rate, scenario)) +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper), size = 0.3) +
  facet_wrap(~feature, ncol = 3) +
  labs(x = paste("Rate of interval-RRB MaxIB (g1, g3) above the",
                 "null-reference 95th percentile (95% exact CI)"),
       y = NULL, caption = "Dashed line: nominal 0.05.") +
  theme_bw()
ggsave(out_png("part2_calibration"), p_cal, width = 11, height = 7.5,
       dpi = 150)

print(summ[, c("id", "label", "n_rep", "bias", "mc_sd",
               "abs_bias_over_mc_sd", "coverage", "se_ratio",
               "rrb_all_maxib_mean", "exceed_rate", "null_q95_mean",
               "mc_q95")], digits = 3, row.names = FALSE)
print(det[, c("id", "label", "exceed_rate", "null_cb_exceed_rate",
              "diff_rrb_minus_cb", "diff_se")], digits = 3,
      row.names = FALSE)
cat("Written:", file.path("analysis", "scenarios_paper.csv"),
    out_csv("paper_scenarios"), out_csv("paper_calibration"),
    out_csv("part1_detection_rrb_vs_cb"), out_csv("paper_within_spearman"),
    out_csv("part1_factorial"), out_png("part1_factorial"),
    out_png("part2_calibration"), sep = "\n  ")
