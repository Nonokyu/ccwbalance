# 02_simulation_study.R -------------------------------------------------------
# Simulation study (report: 検証計画・再現可能資源・推奨プロトコル).
#
# Strategies of the DGM (simulate_ccw(), window = c(w1, w2), default 3, 9):
#   g1  initiate in intervals 0, ..., w1 - 1 (forced at w1 - 1); artificial
#       censoring only at t = w1 - 1 (phase "grace_end");
#   g2  the nested longer window: initiate in 0, ..., w2 - 1 (forced at
#       w2 - 1); artificial censoring only at t = w2 - 1 ("grace_end");
#   g3  no initiation through w2 - 1; censored on initiation in
#       0, ..., w2 - 1 ("no_init").
# NOTE: g2 was redefined (formerly it had a "no_init" phase before w1). The
# saved outputs in output/sim/ and analysis/scenarios_<run>.csv predate this
# redesign (and have the old one-contrast columns) until the runs are redone.
#
# For every scenario in ccw_scenarios(<run>) (analysis/sim_scenarios.R) and
# replicate:
#   * simulate CCW data, fit the scenario's artificial-censoring model,
#   * estimate strategy risks at the horizon and the two contrasts (CONTRASTS)
#       g2g1  theta = risk(g2) - risk(g1)  (longer vs shorter window),
#       g3g2  theta = risk(g3) - risk(g2)  (no initiation vs initiation
#             within the window);
#     |theta_hat - theta_true| with the truth from true_risk_ccw(),
#   * bootstrap SE (B resamples of individuals, censoring model refitted;
#     both contrasts from the same resamples) and a normal-approximation
#     95% CI,
#   * record every diagnostic (proposed + benchmark), aggregated over the
#     three strategies (DIAG_STRATEGIES = g1, g2, g3).
# End-of-grace intervals: eog1 = w1 (first interval after g1's grace
# period) and eog2 = w2 (first interval after g2's grace period and after
# the end of g3's censoring). Contrast cc is evaluated at t_eog(g2g1) = w1
# and t_eog(g3g2) = w2.
#
# Columns of sim_results_<run>.csv:
#   scenario, rep, n, risk_g1, risk_g2, risk_g3;
#   per contrast cc in g2g1, g3g2: theta_hat_<cc>, theta_true_<cc>,
#     se_boot_<cc>, abs_bias_<cc>, bias_<cc> (theta_hat - theta_true) and
#     covered_<cc> (|bias| <= 1.96 se_boot);
#   over g1, g2, g3 (cells with n_censored > 0): rrb_{raw,tr,all}_{iai,maxib},
#     ks_iai, ks_max, mhb_iai, mhb_max, cumrrb_base_iai, cumrrb_base_maxib,
#     jackson_cb_max, jackson_cb_mean, ess_ratio_min, w99_max, pi_min,
#     pr_pi_small_max;
#   end-of-grace strategy SMD per contrast (strategy_balance() with the
#     contrast's first strategy as reference, rows of the second, at
#     t_eog(cc)): eog_smd_base_max_<cc> (baseline X only; the practice),
#     eog_smd_post_max_<cc> (post-baseline L only), eog_smd_max_<cc> (max of
#     both); descriptive max |SMD| over follow-up: strat_smd_L_max_<cc>,
#     strat_smd_X_max_<cc>;
#   deterministic additions (computed after the RNG-using steps):
#     jackson_cb_est_max  Jackson-type censoring balance (max |SMD|) with
#                         the estimator's weights (censoring_balance(...,
#                         weights = "estimator"); responds to truncation),
#     cumrrb_base_eog1    max |cumulative RRB| of baseline X over g1, g2, g3
#                         at t = w1,
#     cumrrb_base_eog2    the same at t = w2,
#     cumrrb_base_eog_diff_g2g1  |cumRRB_X(g2, w1) - cumRRB_X(g1, w1)|,
#     cumrrb_base_eog_diff_g3g2  |cumRRB_X(g3, w2) - cumRRB_X(g2, w2)|,
#     ltfu_cb_max         max |SMD| of unweighted LTFU censoring balance
#                         (censored = "ltfu_next"), NA without LTFU,
#     ltfu_rrb_max        max |unweighted interval RRB| for LTFU over cells
#                         with n_censored > 0, NA without LTFU;
#   individual-bootstrap interval of the cumulative RRB of X (CUM_CELLS),
#     computed in the same bootstrap resamples as se_boot_<cc> (no
#     additional random numbers), <g> in g1, g2, g3 and <tp> in eog1
#     (t = w1), eog2 (t = w2), h (last interval): cumrrb_X_<g>_<tp>_obs,
#     cumrrb_X_<g>_<tp>_boot_se, cumrrb_X_<g>_<tp>_excl0 (1 if the 2.5-97.5%
#     percentile interval excludes 0) and cumrrb_boot_excl0_any (any of the
#     nine cells); NA when B = 0;
#   conditional parametric null reference rrb_null() of interval RRB
#     (CCW_N_NULL_DRAWS draws), over g1, g2, g3: MaxIB over all history
#     functions (null_maxib_obs, null_q95, null_p, null_exceed), the same
#     per history function (suffix _<feature>), null_z_max (max |z| over the
#     cells), null_frac_cells_p05 (fraction of cells with p_two < 0.05),
#     null_n_cells; the null reference of the Jackson-type censoring balance
#     (CB) from the same rrb_null() object: MaxIB of |CB| over cells with a
#     defined observed CB and all history functions (null_cb_maxib_obs,
#     null_cb_q95, null_cb_p, null_cb_exceed, null_cb_n_cells) and the
#     per-feature exceedance null_cb_exceed_<feature>.
#   per strategy <g> in g1, g2, g3 (each strategy has its own censoring
#     model, so each is judged against the null reference of its own
#     cells; appended after all columns above, including bias_<cc> and
#     covered_<cc>; NA when the strategy has no cell): rrb_all_maxib_<g>,
#     jackson_cb_max_<g>, null_maxib_obs_<g>, null_q95_<g>, null_p_<g>,
#     null_exceed_<g>, per history function null_{maxib_obs,q95,p,exceed}_
#     <feature>_<g>, null_z_max_<g>, null_n_cells_<g>, null_cb_maxib_obs_<g>,
#     null_cb_q95_<g>, null_cb_p_<g>, null_cb_exceed_<g>, null_cb_n_cells_<g>
#     and null_cb_exceed_<feature>_<g>. The pooled columns are the maxima
#     over <g> (null_n_cells the sum); their exceedance uses the pooled q95.
#   studentised maximum and censorable-row censoring balance (appended
#     after all columns above; from the same rrb_null() object, so no
#     additional random numbers), pooled over g1, g2, g3 (no suffix) and per
#     strategy (suffix _<g>), [_<g>] below:
#     studentised RRB (per cell |RRB| / null SD, maximum over the cells of
#       null_z_max[_<g>], cells with a non-positive null SD excluded):
#       null_z_q95[_<g>], null_z_p[_<g>], null_z_exceed[_<g>] (null_z_max >
#       null_z_q95), null_z_n_excluded[_<g>], per history function
#       null_z_max_<feature>[_<g>], null_z_q95_<feature>[_<g>],
#       null_z_exceed_<feature>[_<g>];
#     censoring balance over the censorable rows (censor_possible == TRUE;
#       rrb_null()'s cb_cens_* elements; cells with a defined observed
#       value): raw MaxIB null_cb_cens_maxib_obs[_<g>], null_cb_cens_q95
#       [_<g>], null_cb_cens_p[_<g>], null_cb_cens_exceed[_<g>],
#       null_cb_cens_n_cells[_<g>], null_cb_cens_exceed_<feature>[_<g>];
#       studentised null_cb_cens_z_max[_<g>], null_cb_cens_z_q95[_<g>],
#       null_cb_cens_z_exceed[_<g>];
#     jackson_cb_cens_max[_<g>]  max |SMD| of censoring_balance(weighted =
#       TRUE, rows = "censorable") with the model weights.
# Then: estimator performance per scenario and contrast (theta_true_<cc>,
# mean_theta_hat_<cc>, bias_<cc>, rmse_<cc>, mc_sd_<cc>, mean_se_boot_<cc>,
# coverage_<cc>, mean_abs_bias_<cc>), mean diagnostics and the
# null-reference exceedance rate, plus the per-strategy means
# rrb_all_maxib_<g> and null_exceed_rate_<g>, then the exceedance rates
# null_z_exceed_rate[_<g>], null_cb_cens_exceed_rate[_<g>] and
# null_cb_cens_z_exceed_rate[_<g>] appended as the last columns
# (sim_performance_<run>.csv); the tables and
# figures of the manuscript are made by
# analysis/02b_simulation_tables_paper.R.
#
# Runs (CCW_RUN; the four runs that supply the manuscript, see
# analysis/sim_scenarios.R):
#   main  (default) Parts 1-3 ccw_scenarios("main") (Part 1 factorial,
#      Part 2 null-reference calibration, Part 3); n and B per scenario
#      (scenario `n`, `n_boot`), 200 replicates (100 for n = 8000).
#   manycov  ccw_scenarios("manycov"): extended calibration of the null
#      reference (many correlated baseline covariates X2-X9, n = 300 with
#      many covariates, rare artificial censoring, extreme weights) and two
#      detection scenarios with an omitted X2 correlated (rho = 0.5) or not
#      (rho = 0) with X. Run like main (200 replicates, scenario n / n_boot);
#      the scenario's extra covariates (`extra_feats`) are added to the
#      monitored history functions, so the per-feature null columns
#      _X2 ... _X9 exist only for those scenarios (NA elsewhere).
#   largerep  ccw_scenarios("largerep"): large-replicate calibration of the
#      null reference for RRB and censoring balance under the correct model
#      (n = 300 / 2000 / 8000), early events in g1's grace period (end-of-grace
#      SMD of X non-zero under correct weights through outcome depletion;
#      cumulative RRB of X stays ~0), and the omitted-L detection scenario
#      with truncation at the 95th percentile or near-positivity violation.
#      The number of replicates is the scenario's `n_rep` (1000 / 500 / 200).
#   variants  ccw_scenarios("variants"): omitted L_t alone (beta2 = beta5
#      = 0), the L_t^2 cells without outcome effect at gamma0 = -4.5, and
#      early_events with 500 bootstrap resamples (scenario n_rep = 200).
# Outputs are suffixed by the run: <CCW_SIM_DIR>/sim_results_<run>.csv,
# sim_performance_<run>.csv, sim_settings_<run>.csv and (only with the
# default CCW_SIM_DIR) analysis/scenarios_<run>.csv.
#
# Configuration (environment variables):
#   CCW_RUN     main | manycov | largerep | variants              [main]
#   CCW_N_REPS  replicates per scenario    [main, manycov: 200, n = 8000:
#                                           100; largerep, variants: n_rep]
#   CCW_N_SUBJ  individuals per replicate (unless the scenario sets n)
#                                                                 [2000]
#   CCW_N_BOOT  bootstrap resamples per replicate for the CI
#                                                       [scenario n_boot]
#   CCW_N_NULL_DRAWS  null-reference draws per replicate          [500]
#   CCW_N_MC    Monte Carlo size for the true risks             [5e5]
#   CCW_OUT_SUFFIX  output suffix                               [run]
#   CCW_CORES   parallel workers (parallel::mclapply)       [detectCores()-1]
#   CCW_SIM_DIR output directory of the sim_*.csv files     [output/sim]
#               (any other value, e.g. a scratch directory for a smoke
#               test, also leaves analysis/scenarios_<run>.csv untouched)
# Random numbers: one L'Ecuyer-CMRG stream per task, assigned in task order
# (heavy scenarios first, by n x (B + 10)), so changing a scenario's n or B
# (or the number of replicates) reorders the tasks and changes the draws of
# every replicate. Data generation, bootstrap and null draws share the
# replicate's stream, except that rrb_null() is reseeded with a seed drawn
# from that stream after the bootstrap. The seed of a run is 20260923 plus
# its fixed offset (SEED_OFFSET; the offsets of the former script version,
# kept so that the draws do not change). The g2 redesign changes the
# simulated data (and hence the results) but not this stream design.
#
# Run from the project root:  Rscript analysis/02_simulation_study.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) as.numeric(Sys.getenv(x, as.character(default)))
RUN <- Sys.getenv("CCW_RUN", "main")
# seed offsets of the runs (fixed; do not change: they set the random
# numbers of the saved results)
SEED_OFFSET <- c(main = 3, manycov = 4, largerep = 5, variants = 7)
stopifnot(RUN %in% names(SEED_OFFSET))
SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", RUN)
# NA = per-scenario default (see reps_of(), nboot_of())
N_REPS <- env_num("CCW_N_REPS", NA)
N_SUBJ <- env_num("CCW_N_SUBJ", 2000)
N_BOOT <- env_num("CCW_N_BOOT", NA)
N_NULL_DRAWS <- env_num("CCW_N_NULL_DRAWS", 500)
N_MC <- env_num("CCW_N_MC", 5e5)
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
# contrasts theta_<cc> = risk(cc[1]) - risk(cc[2])
CONTRASTS <- list(g2g1 = c("g2", "g1"), g3g2 = c("g3", "g2"))
CC <- names(CONTRASTS)
DIAG_STRATEGIES <- c("g1", "g2", "g3")  # diagnostics aggregated over these
SEED <- 20260923

DEFAULT_SIM_DIR <- file.path("output", "sim")
out_dir <- Sys.getenv("CCW_SIM_DIR", DEFAULT_SIM_DIR)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- ccw_scenarios(RUN)
if (identical(out_dir, DEFAULT_SIM_DIR)) {
  write.csv(scenario_table(scenarios),
            file.path("analysis", sprintf("scenarios_%s.csv", RUN)),
            row.names = FALSE)
}

raw_feats <- c("X", "L")
tr_feats <- c("L_sq", "L_lagprod", "L_diff", "L_cummean", "L_gt0", "L_x_X")
all_feats <- c(raw_feats, tr_feats)
base_feats <- c("X", "L_base")
# history functions of a scenario: all_feats plus its extra baseline
# covariates (run manycov; identical to all_feats otherwise)
feats_of <- function(sc) c(all_feats, as.character(sc$extra_feats))

# per-scenario replicates and bootstrap size (run main defaults; the
# scenarios of runs largerep and variants set `n_rep`)
n_of <- function(sc) if (is.null(sc$params[["n"]])) N_SUBJ else sc$params[["n"]]
reps_of <- function(sc) {
  if (!is.na(N_REPS)) return(N_REPS)
  if (!is.null(sc$n_rep)) return(sc$n_rep)
  if (n_of(sc) >= 8000) 100 else 200
}
nboot_of <- function(sc) {
  if (!is.na(N_BOOT)) return(N_BOOT)
  if (!is.null(sc$n_boot)) sc$n_boot else 100
}

# ---------------------------------------------------------------------------
# helpers

# named vector of the contrasts at the horizon
contrast_of <- function(est) {
  est <- est[est$time == max(est$time), ]
  r <- stats::setNames(est$risk, as.character(est$strategy))
  vapply(CONTRASTS, function(cc) unname(r[cc[1]] - r[cc[2]]), numeric(1))
}

# end-of-grace times: eog1 = w1 (after g1's grace), eog2 = w2 (after g2's
# grace and g3's censoring); t_eog of contrast g2g1 is w1, of g3g2 w2
eog_times <- function(pars) c(eog1 = pars$window[1], eog2 = pars$window[2])
T_EOG_CC <- c(g2g1 = "eog1", g3g2 = "eog2")

max_or_na <- function(v) {
  v <- abs(v[!is.na(v)])
  if (length(v) == 0) NA_real_ else max(v)
}

# bootstrap SE of the contrasts; in the same resamples (no additional random
# numbers) the cumulative RRB of X for g1, g2, g3 at eog1, eog2 and the
# horizon. `cum_t` is the named vector c(eog1, eog2, h) of times. Returns
# list(se, cum): `se` named over the contrasts, `cum` a data frame
# (strategy, tp, boot_se, lcl, ucl) or NULL when B <= 0.
CUM_CELLS <- expand.grid(strategy = DIAG_STRATEGIES,
                         tp = c("eog1", "eog2", "h"),
                         stringsAsFactors = FALSE)
boot_se <- function(d, sc, B, cum_t) {
  if (B <= 0) {
    return(list(se = stats::setNames(rep(NA_real_, length(CC)), CC),
                cum = NULL))
  }
  idx <- split(seq_len(nrow(d)), d$id)
  ids <- names(idx)
  base_cols <- c("id", "strategy", "time", "censored_next", "censor_possible",
                 "phase", "event_next", "ltfu_next", "X", "L", "L_lag",
                 grep("^X[0-9]+$", names(d), value = TRUE))  # run manycov
  d <- d[, base_cols]
  cell_t <- unname(cum_t[CUM_CELLS$tp])
  n_cc <- length(CC)
  th <- vapply(seq_len(B), function(b) {
    s <- sample(ids, length(ids), replace = TRUE)
    rows <- unlist(idx[s], use.names = FALSE)
    db <- d[rows, ]
    db$id <- rep(seq_along(s), lengths(idx[s]))
    wb <- fit_censoring_model(db, sc$formula, truncate = sc$truncate)
    # deterministic given the resample: no random numbers drawn below
    rcb <- rrb(wb, "X", scope = "cumulative")
    k <- match(paste(CUM_CELLS$strategy, cell_t),
               paste(as.character(rcb$strategy), rcb$time))
    c(contrast_of(estimate_risk(wb)), rcb$rrb[k])
  }, numeric(n_cc + nrow(CUM_CELLS)))
  th <- matrix(th, nrow = n_cc + nrow(CUM_CELLS))
  q <- function(v, p) {
    if (all(is.na(v))) NA_real_ else
      stats::quantile(v, p, na.rm = TRUE, names = FALSE)
  }
  thc <- th[-seq_len(n_cc), , drop = FALSE]
  cum <- CUM_CELLS
  cum$boot_se <- apply(thc, 1, function(v) {
    if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE)
  })
  cum$lcl <- apply(thc, 1, q, 0.025)
  cum$ucl <- apply(thc, 1, q, 0.975)
  se <- stats::setNames(apply(th[seq_len(n_cc), , drop = FALSE], 1,
                              stats::sd), CC)
  list(se = se, cum = cum)
}

agg_max <- function(x, value, feats) {
  x <- x[as.character(x$strategy) %in% DIAG_STRATEGIES &
           as.character(x$feature) %in% feats & x$n_censored > 0, ]
  if (nrow(x) == 0) return(c(iai = NA, maxib = NA))
  a <- aggregate_balance(x, value = value, by = c("strategy", "feature"))
  c(iai = max(a$iai), maxib = max(a$maxib))
}

# MaxIB per diagnostic strategy (named over DIAG_STRATEGIES; the same cells
# as agg_max() restricted to one strategy, NA when it has none)
agg_max_by_strategy <- function(x, value, feats) {
  vapply(DIAG_STRATEGIES, function(g) {
    a <- agg_max(x[as.character(x$strategy) == g, ], value, feats)
    as.numeric(a[["maxib"]])
  }, numeric(1))
}

one_replicate <- function(sc, rep, truth) {
  pars <- do.call(ccw_sim_params, sc$params)
  if (is.null(sc$params[["n"]])) pars$n <- N_SUBJ
  d <- simulate_ccw(pars)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, sc$formula, truncate = sc$truncate)

  est <- estimate_risk(w)
  est_h <- est[est$time == max(est$time), ]
  risk <- stats::setNames(est_h$risk, as.character(est_h$strategy))
  theta_hat <- contrast_of(est)
  t_eog <- eog_times(pars)
  cum_t <- c(t_eog, h = max(d$time))
  bs <- boot_se(d, sc, nboot_of(sc), cum_t)
  se <- bs$se
  # seed of the null reference, drawn from the replicate's stream
  null_seed <- sample.int(.Machine$integer.max, 1)
  feats_sc <- feats_of(sc)

  r <- rrb(w, feats_sc)
  k <- ks_rr(w, feats_sc)
  m <- mahalanobis_rr(w, feats_sc, ridge = 0.01)
  m <- m[as.character(m$strategy) %in% DIAG_STRATEGIES & m$n_censored > 0, ]
  rc <- rrb(w, base_feats, scope = "cumulative")
  ws <- weight_support(w)
  ws <- ws[as.character(ws$strategy) %in% DIAG_STRATEGIES, ]
  cb <- censoring_balance(w, feats_sc, weighted = TRUE)
  cb <- cb[as.character(cb$strategy) %in% DIAG_STRATEGIES, ]

  a_raw <- agg_max(r, "rrb", raw_feats)
  a_tr <- agg_max(r, "rrb", tr_feats)
  a_all <- agg_max(r, "rrb", feats_sc)
  a_ks <- agg_max(k, "ks", feats_sc)
  a_cum <- agg_max(rc, "rrb", base_feats)
  mh_iai <- if (nrow(m) > 0) {
    max(tapply(m$mahalanobis, as.character(m$strategy), mean, na.rm = TRUE))
  } else NA

  # estimator: risks and, per contrast, estimate, truth, bootstrap SE
  out <- list(scenario = sc$name, rep = rep, n = pars$n,
              risk_g1 = risk[["g1"]], risk_g2 = risk[["g2"]],
              risk_g3 = risk[["g3"]])
  for (cc in CC) {
    out[[paste0("theta_hat_", cc)]] <- theta_hat[[cc]]
    out[[paste0("theta_true_", cc)]] <- truth[[cc]]
    out[[paste0("se_boot_", cc)]] <- se[[cc]]
    out[[paste0("abs_bias_", cc)]] <- abs(theta_hat[[cc]] - truth[[cc]])
  }
  out <- c(out, list(
    rrb_raw_iai = a_raw[["iai"]], rrb_raw_maxib = a_raw[["maxib"]],
    rrb_tr_iai = a_tr[["iai"]], rrb_tr_maxib = a_tr[["maxib"]],
    rrb_all_iai = a_all[["iai"]], rrb_all_maxib = a_all[["maxib"]],
    ks_iai = a_ks[["iai"]], ks_max = a_ks[["maxib"]],
    mhb_iai = mh_iai, mhb_max = max(m$mahalanobis, na.rm = TRUE),
    cumrrb_base_iai = a_cum[["iai"]], cumrrb_base_maxib = a_cum[["maxib"]],
    jackson_cb_max = max(abs(cb$smd), na.rm = TRUE),
    jackson_cb_mean = mean(abs(cb$smd), na.rm = TRUE),
    ess_ratio_min = min(ws$ess_ratio, na.rm = TRUE),
    w99_max = max(ws$w99, na.rm = TRUE),
    pi_min = min(ws$pi_min, na.rm = TRUE),
    pr_pi_small_max = max(ws$pr_pi_small, na.rm = TRUE)
  ))
  # end-of-grace strategy SMD per contrast (reference cc[1], rows cc[2]) at
  # t_eog(cc): baseline covariate only (the practice), post-baseline L only
  # (extension), and the maximum over both; descriptive strategy balance of
  # L and X over follow-up
  for (cc in CC) {
    sb <- strategy_balance(w, raw_feats, reference = CONTRASTS[[cc]][1])
    sb <- sb[as.character(sb$strategy) == CONTRASTS[[cc]][2], ]
    eog <- sb[sb$time == t_eog[[T_EOG_CC[[cc]]]], ]
    out[[paste0("eog_smd_max_", cc)]] <- max_or_na(eog$smd)
    out[[paste0("eog_smd_base_max_", cc)]] <-
      max_or_na(eog$smd[eog$feature == "X"])
    out[[paste0("eog_smd_post_max_", cc)]] <-
      max_or_na(eog$smd[eog$feature == "L"])
    out[[paste0("strat_smd_L_max_", cc)]] <-
      max_or_na(sb$smd[sb$feature == "L"])
    out[[paste0("strat_smd_X_max_", cc)]] <-
      max_or_na(sb$smd[sb$feature == "X"])
  }

  # additional deterministic diagnostics (no random numbers are drawn here)
  # Jackson-type censoring balance with the estimator's (truncated) weights
  cbe <- censoring_balance(w, feats_sc, weighted = TRUE,
                           weights = "estimator")
  cbe <- cbe[as.character(cbe$strategy) %in% DIAG_STRATEGIES, ]
  # cumulative RRB of baseline X at the end-of-grace intervals
  rc_x <- rc[rc$feature == "X" &
               as.character(rc$strategy) %in% DIAG_STRATEGIES, ]
  cum_x <- function(g, tt) {
    v <- rc_x$rrb[as.character(rc_x$strategy) == g & rc_x$time == tt]
    if (length(v) == 1) v else NA_real_
  }
  # loss to follow-up: censoring balance and RRB (unweighted; NA if no LTFU)
  cbl <- censoring_balance(w, feats_sc, weighted = FALSE,
                           censored = "ltfu_next")
  cbl <- cbl[as.character(cbl$strategy) %in% DIAG_STRATEGIES, ]
  rl <- rrb(w, feats_sc, method = "unweighted", censored = "ltfu_next")
  rl <- rl[as.character(rl$strategy) %in% DIAG_STRATEGIES &
             rl$n_censored > 0, ]
  out <- c(out, list(
    jackson_cb_est_max = max(abs(cbe$smd), na.rm = TRUE),
    cumrrb_base_eog1 = max_or_na(rc_x$rrb[rc_x$time == t_eog[["eog1"]]]),
    cumrrb_base_eog2 = max_or_na(rc_x$rrb[rc_x$time == t_eog[["eog2"]]])
  ))
  for (cc in CC) {
    tt <- t_eog[[T_EOG_CC[[cc]]]]
    out[[paste0("cumrrb_base_eog_diff_", cc)]] <-
      abs(cum_x(CONTRASTS[[cc]][1], tt) - cum_x(CONTRASTS[[cc]][2], tt))
  }
  out$ltfu_cb_max <- max_or_na(cbl$smd)
  out$ltfu_rrb_max <- max_or_na(rl$rrb)

  # cumulative RRB of X (g1, g2, g3; eog1, eog2, horizon): observed value
  # and individual-bootstrap SE / percentile interval (from boot_se())
  excl <- logical(0)
  for (k in seq_len(nrow(CUM_CELLS))) {
    g <- CUM_CELLS$strategy[k]
    tp <- CUM_CELLS$tp[k]
    nm <- sprintf("cumrrb_X_%s_%s", g, tp)
    out[[paste0(nm, "_obs")]] <- cum_x(g, cum_t[[tp]])
    bk <- if (is.null(bs$cum)) NULL else bs$cum[k, ]
    out[[paste0(nm, "_boot_se")]] <- if (is.null(bk)) NA_real_ else
      bk$boot_se
    ex <- if (is.null(bk) || is.na(bk$lcl) || is.na(bk$ucl)) NA_real_ else
      as.numeric(bk$lcl > 0 | bk$ucl < 0)
    out[[paste0(nm, "_excl0")]] <- ex
    excl <- c(excl, ex)
  }
  out$cumrrb_boot_excl0_any <- if (all(is.na(excl))) NA_real_ else
    as.numeric(any(excl == 1, na.rm = TRUE))
  # per strategy: MaxIB of interval RRB and max |CB| (moved to the end of
  # the results with the per-strategy null columns)
  a_all_g <- agg_max_by_strategy(r, "rrb", feats_sc)
  for (g in DIAG_STRATEGIES) {
    out[[paste0("rrb_all_maxib_", g)]] <- a_all_g[[g]]
    out[[paste0("jackson_cb_max_", g)]] <-
      max_or_na(cb$smd[as.character(cb$strategy) == g])
  }
  # censorable-row censoring balance with the model weights (deterministic)
  cbc <- censoring_balance(w, feats_sc, weighted = TRUE, rows = "censorable")
  cbc <- cbc[as.character(cbc$strategy) %in% DIAG_STRATEGIES, ]
  # null reference (last step: rrb_null() reseeds the RNG)
  nr <- rrb_null(w, feats_sc, n_draws = N_NULL_DRAWS, seed = null_seed)
  ext <- null_summary_ext(nr, feats_sc, cbc)
  res <- cbind(as.data.frame(out), null_summary(nr, feats_sc),
               as.data.frame(ext))
  # columns added after the four runs were first defined; moved behind all
  # other columns when the results are combined
  attr(res, "ext_cols") <- names(ext)
  res
}

# conditional parametric null reference of interval RRB, MaxIB over the
# diagnostic strategies (g1, g2, g3) and all / each history function; the
# same for the censoring balance (CB). `nr` is the rrb_null() object of the
# replicate
null_summary <- function(nr, feats = all_feats) {
  cells <- nr$cells
  sims <- matrix(nr$draws$rrb, nrow = nrow(cells))  # cells x draws
  in_s <- as.character(cells$strategy) %in% DIAG_STRATEGIES
  maxib_ref <- function(sel) {
    if (!any(sel)) return(c(obs = NA, q95 = NA, p = NA, exceed = NA))
    obs <- max(abs(cells$rrb_obs[sel]))
    m <- apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
    q95 <- stats::quantile(m, 0.95, na.rm = TRUE, names = FALSE)
    c(obs = obs, q95 = q95,
      p = (1 + sum(m >= obs, na.rm = TRUE)) / (1 + sum(!is.na(m))),
      exceed = as.numeric(obs > q95))
  }
  all_ref <- maxib_ref(in_s)
  res <- list(null_maxib_obs = all_ref[["obs"]], null_q95 = all_ref[["q95"]],
              null_p = all_ref[["p"]], null_exceed = all_ref[["exceed"]])
  for (f in feats) {
    fr <- maxib_ref(in_s & as.character(cells$feature) == f)
    res[[paste0("null_maxib_obs_", f)]] <- fr[["obs"]]
    res[[paste0("null_q95_", f)]] <- fr[["q95"]]
    res[[paste0("null_p_", f)]] <- fr[["p"]]
    res[[paste0("null_exceed_", f)]] <- fr[["exceed"]]
  }
  z <- abs(cells$z[in_s])
  res$null_z_max <- if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  res$null_frac_cells_p05 <- mean(cells$p_two[in_s] < 0.05)
  res$null_n_cells <- sum(in_s)
  # censoring balance (CB) from the same draws: MaxIB of |CB| over the
  # g1/g2/g3 cells with a defined observed CB (as in rrb_null()'s maxib_cb)
  sims_cb <- matrix(nr$draws$cb, nrow = nrow(cells))
  cb_ok <- in_s & !is.na(cells$cb_obs)
  cb_ref <- function(sel) {
    if (!any(sel)) return(c(obs = NA, q95 = NA, p = NA, exceed = NA))
    obs <- max(abs(cells$cb_obs[sel]))
    m <- apply(abs(sims_cb[sel, , drop = FALSE]), 2, function(v) {
      if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    })
    if (all(is.na(m))) return(c(obs = obs, q95 = NA, p = NA, exceed = NA))
    q95 <- stats::quantile(m, 0.95, na.rm = TRUE, names = FALSE)
    c(obs = obs, q95 = q95,
      p = (1 + sum(m >= obs, na.rm = TRUE)) / (1 + sum(!is.na(m))),
      exceed = as.numeric(obs > q95))
  }
  cb_all <- cb_ref(cb_ok)
  res$null_cb_maxib_obs <- cb_all[["obs"]]
  res$null_cb_q95 <- cb_all[["q95"]]
  res$null_cb_p <- cb_all[["p"]]
  res$null_cb_exceed <- cb_all[["exceed"]]
  res$null_cb_n_cells <- sum(cb_ok)
  for (f in feats) {
    res[[paste0("null_cb_exceed_", f)]] <-
      cb_ref(cb_ok & as.character(cells$feature) == f)[["exceed"]]
  }
  # per strategy (each has its own censoring model): the same summaries
  # over the cells of strategy g only
  is_f <- function(f) as.character(cells$feature) == f
  for (g in DIAG_STRATEGIES) {
    in_g <- in_s & as.character(cells$strategy) == g
    gr <- maxib_ref(in_g)
    res[[paste0("null_maxib_obs_", g)]] <- gr[["obs"]]
    res[[paste0("null_q95_", g)]] <- gr[["q95"]]
    res[[paste0("null_p_", g)]] <- gr[["p"]]
    res[[paste0("null_exceed_", g)]] <- gr[["exceed"]]
    for (f in feats) {
      fr <- maxib_ref(in_g & is_f(f))
      res[[sprintf("null_maxib_obs_%s_%s", f, g)]] <- fr[["obs"]]
      res[[sprintf("null_q95_%s_%s", f, g)]] <- fr[["q95"]]
      res[[sprintf("null_p_%s_%s", f, g)]] <- fr[["p"]]
      res[[sprintf("null_exceed_%s_%s", f, g)]] <- fr[["exceed"]]
    }
    z <- abs(cells$z[in_g])
    res[[paste0("null_z_max_", g)]] <- if (all(is.na(z))) NA_real_ else
      max(z, na.rm = TRUE)
    res[[paste0("null_n_cells_", g)]] <- sum(in_g)
    cb_g <- cb_ok & as.character(cells$strategy) == g
    cr <- cb_ref(cb_g)
    res[[paste0("null_cb_maxib_obs_", g)]] <- cr[["obs"]]
    res[[paste0("null_cb_q95_", g)]] <- cr[["q95"]]
    res[[paste0("null_cb_p_", g)]] <- cr[["p"]]
    res[[paste0("null_cb_exceed_", g)]] <- cr[["exceed"]]
    res[[paste0("null_cb_n_cells_", g)]] <- sum(cb_g)
    for (f in feats) {
      res[[sprintf("null_cb_exceed_%s_%s", f, g)]] <-
        cb_ref(cb_g & is_f(f))[["exceed"]]
    }
  }
  as.data.frame(res)
}

# observed vs null maximum of |value| over the cells `sel` (`obs_v` one
# value per cell, `sims` cells x draws). With `scale` (one value per cell,
# the null SD) the studentised maximum: values and draws divided by
# `scale`, cells with NA or non-positive `scale` excluded (`n_excluded`),
# as in rrb_null()'s maxib_stud
max_ref <- function(obs_v, sims, sel, scale = NULL) {
  n_exc <- 0
  if (!is.null(scale)) {
    ok_s <- !is.na(scale) & scale > 0
    n_exc <- sum(sel & !ok_s)
    sel <- sel & ok_s
  }
  na <- c(obs = NA, q95 = NA, p = NA, exceed = NA, n_excluded = n_exc)
  if (!any(sel)) return(na)
  ov <- obs_v[sel]
  sv <- sims[sel, , drop = FALSE]
  if (!is.null(scale)) {
    ov <- ov / scale[sel]
    sv <- sv / scale[sel]
  }
  obs <- max(abs(ov))
  m <- apply(abs(sv), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  if (all(is.na(m))) {
    na[["obs"]] <- obs
    return(na)
  }
  q95 <- stats::quantile(m, 0.95, na.rm = TRUE, names = FALSE)
  c(obs = obs, q95 = q95,
    p = (1 + sum(m >= obs, na.rm = TRUE)) / (1 + sum(!is.na(m))),
    exceed = as.numeric(obs > q95), n_excluded = n_exc)
}

# columns added for the studentised maximum and the censorable-row censoring
# balance (appended after all other columns), pooled over g1, g2, g3 (no
# suffix) and per strategy (suffix _<g>), from the same rrb_null() object as
# null_summary() (no additional random numbers). `cbc` is
# censoring_balance(rows = "censorable") with the model weights, restricted
# to DIAG_STRATEGIES
null_summary_ext <- function(nr, feats, cbc) {
  cells <- nr$cells
  sims <- matrix(nr$draws$rrb, nrow = nrow(cells))
  sims_cc <- matrix(nr$draws$cb_cens, nrow = nrow(cells))
  in_s <- as.character(cells$strategy) %in% DIAG_STRATEGIES
  is_f <- function(f) as.character(cells$feature) == f
  cc_ok <- !is.na(cells$cb_cens_obs)
  block <- function(sel, sfx, cb_smd) {
    nm <- function(x) paste0(x, sfx)
    r <- list()
    zr <- max_ref(cells$rrb_obs, sims, sel, cells$null_sd)
    r[[nm("null_z_q95")]] <- zr[["q95"]]
    r[[nm("null_z_p")]] <- zr[["p"]]
    r[[nm("null_z_exceed")]] <- zr[["exceed"]]
    r[[nm("null_z_n_excluded")]] <- zr[["n_excluded"]]
    for (f in feats) {
      fr <- max_ref(cells$rrb_obs, sims, sel & is_f(f), cells$null_sd)
      r[[nm(paste0("null_z_max_", f))]] <- fr[["obs"]]
      r[[nm(paste0("null_z_q95_", f))]] <- fr[["q95"]]
      r[[nm(paste0("null_z_exceed_", f))]] <- fr[["exceed"]]
    }
    sel_c <- sel & cc_ok
    cr <- max_ref(cells$cb_cens_obs, sims_cc, sel_c)
    r[[nm("null_cb_cens_maxib_obs")]] <- cr[["obs"]]
    r[[nm("null_cb_cens_q95")]] <- cr[["q95"]]
    r[[nm("null_cb_cens_p")]] <- cr[["p"]]
    r[[nm("null_cb_cens_exceed")]] <- cr[["exceed"]]
    r[[nm("null_cb_cens_n_cells")]] <- sum(sel_c)
    for (f in feats) {
      r[[nm(paste0("null_cb_cens_exceed_", f))]] <-
        max_ref(cells$cb_cens_obs, sims_cc, sel_c & is_f(f))[["exceed"]]
    }
    cz <- max_ref(cells$cb_cens_obs, sims_cc, sel_c, cells$cb_cens_null_sd)
    r[[nm("null_cb_cens_z_max")]] <- cz[["obs"]]
    r[[nm("null_cb_cens_z_q95")]] <- cz[["q95"]]
    r[[nm("null_cb_cens_z_exceed")]] <- cz[["exceed"]]
    r[[nm("jackson_cb_cens_max")]] <- max_or_na(cb_smd)
    r
  }
  res <- block(in_s, "", cbc$smd)
  for (g in DIAG_STRATEGIES) {
    res <- c(res, block(in_s & as.character(cells$strategy) == g,
                        paste0("_", g),
                        cbc$smd[as.character(cbc$strategy) == g]))
  }
  res
}

# ---------------------------------------------------------------------------
# truth and task list (one independent L'Ecuyer-CMRG stream per task)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + SEED_OFFSET[[RUN]])
t_start <- Sys.time()
# one true_risk_ccw() call per scenario; named vector over the contrasts
truths <- lapply(scenarios, function(sc) {
  tr <- true_risk_ccw(do.call(ccw_sim_params, sc$params), n_mc = N_MC)
  r <- stats::setNames(tr$risk, as.character(tr$strategy))
  vapply(CONTRASTS, function(cc) unname(r[cc[1]] - r[cc[2]]), numeric(1))
})

# heavy scenarios (large n / B) first for better load balancing
tasks <- do.call(rbind, lapply(seq_along(scenarios), function(s) {
  data.frame(s = s, rep = seq_len(reps_of(scenarios[[s]])))
}))
cost <- vapply(scenarios, function(sc) n_of(sc) * (nboot_of(sc) + 10),
               numeric(1))
tasks <- tasks[order(-cost[tasks$s], tasks$rep), ]
stream <- .Random.seed
seeds <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  stream <- parallel::nextRNGStream(stream)
  seeds[[i]] <- stream
}

message(sprintf(paste("Running %d tasks (%d scenarios), per-scenario n / B /",
                      "reps, %d null draws, %d cores"),
                nrow(tasks), length(scenarios), N_NULL_DRAWS, N_CORES))
res_list <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
  assign(".Random.seed", seeds[[i]], envir = globalenv())
  sc <- scenarios[[tasks$s[i]]]
  tryCatch(one_replicate(sc, tasks$rep[i], truths[[sc$name]]),
           error = function(e) {
             data.frame(scenario = sc$name, rep = tasks$rep[i],
                        error = conditionMessage(e))
           })
}, mc.cores = N_CORES, mc.preschedule = FALSE)
runtime <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

failed <- Filter(function(x) "error" %in% names(x), res_list)
if (length(failed) > 0) {
  warning(length(failed), " replicate(s) failed; first: ", failed[[1]]$error)
}
ok_list <- Filter(function(x) !("error" %in% names(x)), res_list)
ext_cols <- unique(unlist(lapply(ok_list, attr, "ext_cols")))
# runs manycov, largerep, variants: per-feature null columns may differ
# between scenarios -> bind_rows (fills NA); run main keeps rbind
res <- if (RUN %in% c("manycov", "largerep", "variants")) {
  as.data.frame(dplyr::bind_rows(ok_list))
} else {
  do.call(rbind, ok_list)
}
for (cc in CC) {
  b <- res[[paste0("theta_hat_", cc)]] - res[[paste0("theta_true_", cc)]]
  res[[paste0("bias_", cc)]] <- b
  res[[paste0("covered_", cc)]] <- abs(b) <= 1.96 * res[[paste0("se_boot_",
                                                                cc)]]
}
# per-strategy columns (suffix _g1, _g2, _g3; not risk_<g>) last, so that
# the pre-existing columns keep their positions
# pre-existing columns keep their positions; the studentised and
# censorable-row columns (ext_cols, see null_summary_ext()) come last
old_cols <- setdiff(names(res), ext_cols)
new_cols <- setdiff(
  grep(sprintf("_(%s)$", paste(DIAG_STRATEGIES, collapse = "|")),
       old_cols, value = TRUE),
  paste0("risk_", DIAG_STRATEGIES))
res <- res[, c(setdiff(old_cols, new_cols), new_cols,
               intersect(ext_cols, names(res)))]
write.csv(res, file.path(out_dir, sprintf("sim_results_%s.csv", SUFFIX)),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# performance by scenario (per contrast) and mean diagnostics
perf <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  out <- data.frame(scenario = x$scenario[1], n_rep = nrow(x))
  for (cc in CC) {
    v <- function(nm) x[[paste0(nm, "_", cc)]]
    out[[paste0("theta_true_", cc)]] <- v("theta_true")[1]
    out[[paste0("mean_theta_hat_", cc)]] <- mean(v("theta_hat"))
    out[[paste0("bias_", cc)]] <- mean(v("bias"))
    out[[paste0("rmse_", cc)]] <- sqrt(mean(v("bias")^2))
    out[[paste0("mc_sd_", cc)]] <- stats::sd(v("theta_hat"))
    out[[paste0("mean_se_boot_", cc)]] <- mean(v("se_boot"))
    out[[paste0("coverage_", cc)]] <- mean(v("covered"))
    out[[paste0("mean_abs_bias_", cc)]] <- mean(v("abs_bias"))
  }
  out$rrb_all_iai <- mean(x$rrb_all_iai)
  out$rrb_all_maxib <- mean(x$rrb_all_maxib)
  out$mhb_max <- mean(x$mhb_max)
  out$cumrrb_base_maxib <- mean(x$cumrrb_base_maxib)
  for (cc in CC) {
    nm <- paste0("strat_smd_L_max_", cc)
    out[[nm]] <- mean(x[[nm]], na.rm = TRUE)
  }
  out$ess_ratio_min <- mean(x$ess_ratio_min)
  out$w99_max <- mean(x$w99_max)
  out$jackson_cb_est_max <- mean(x$jackson_cb_est_max)
  out$cumrrb_base_eog1 <- mean(x$cumrrb_base_eog1)
  out$ltfu_cb_max <- mean(x$ltfu_cb_max, na.rm = TRUE)
  out$ltfu_rrb_max <- mean(x$ltfu_rrb_max, na.rm = TRUE)
  out
}))
# null-reference summaries; tables and figures of the manuscript in 02b
perf_d <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  data.frame(scenario = x$scenario[1],
             null_exceed_rate = mean(x$null_exceed, na.rm = TRUE),
             null_q95 = mean(x$null_q95, na.rm = TRUE),
             null_frac_cells_p05 = mean(x$null_frac_cells_p05,
                                        na.rm = TRUE))
}))
perf <- merge(perf, perf_d, by = "scenario", sort = FALSE)
# per-strategy means, appended after the pooled columns
perf_g <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  out <- data.frame(scenario = x$scenario[1])
  for (g in DIAG_STRATEGIES) {
    out[[paste0("rrb_all_maxib_", g)]] <-
      mean(x[[paste0("rrb_all_maxib_", g)]], na.rm = TRUE)
  }
  for (g in DIAG_STRATEGIES) {
    out[[paste0("null_exceed_rate_", g)]] <-
      mean(x[[paste0("null_exceed_", g)]], na.rm = TRUE)
  }
  out
}))
perf <- merge(perf, perf_g, by = "scenario", sort = FALSE)
# exceedance rates of the studentised maximum and the censorable-row
# censoring balance (raw and studentised), pooled and per strategy
perf_ext <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  out <- data.frame(scenario = x$scenario[1])
  for (sfx in c("", paste0("_", DIAG_STRATEGIES))) {
    for (st in c("null_z", "null_cb_cens", "null_cb_cens_z")) {
      out[[paste0(st, "_exceed_rate", sfx)]] <-
        mean(x[[paste0(st, "_exceed", sfx)]], na.rm = TRUE)
    }
  }
  out
}))
perf <- merge(perf, perf_ext, by = "scenario", sort = FALSE)
perf <- perf[match(names(scenarios), perf$scenario), ]
perf <- perf[!is.na(perf$scenario), ]
write.csv(perf, file.path(out_dir, sprintf("sim_performance_%s.csv", SUFFIX)),
          row.names = FALSE)
reps_d <- vapply(scenarios, reps_of, numeric(1))
boot_d <- vapply(scenarios, nboot_of, numeric(1))
settings <- data.frame(
  RUN,
  N_REPS = paste(sprintf("%s=%d", names(reps_d), reps_d), collapse = ";"),
  N_BOOT = paste(sprintf("%s=%d", names(boot_d), boot_d), collapse = ";"),
  N_NULL_DRAWS, N_MC, N_CORES, n_tasks = nrow(tasks),
  n_failed = length(failed), runtime_min = round(runtime, 2))
write.csv(settings, file.path(out_dir, sprintf("sim_settings_%s.csv", SUFFIX)),
          row.names = FALSE)
print(settings[, c("RUN", "N_NULL_DRAWS", "N_MC", "n_tasks", "n_failed",
                   "runtime_min")])
print(perf[, c("scenario",
               paste0(c("theta_true_", "bias_", "coverage_"), "g2g1"),
               paste0(c("theta_true_", "bias_", "coverage_"), "g3g2"),
               "rrb_all_maxib", "null_exceed_rate")],
      digits = 3, row.names = FALSE)
