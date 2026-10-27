# 02_simulation_study.R -------------------------------------------------------
# Simulation study of the manuscript (Sections 5-6).
#
# For every scenario in ccw_scenarios(<run>) (analysis/sim_scenarios.R) and
# replicate:
#   * simulate CCW data, fit the scenario's artificial-censoring model,
#   * estimate strategy risks at the horizon and the contrast
#     theta = risk(g3) - risk(g1); |theta_hat - theta_true| with the truth from
#     true_risk_ccw(),
#   * bootstrap SE (B resamples of individuals, censoring model refitted) and
#     a normal-approximation 95% CI,
#   * record every diagnostic (proposed + benchmark), aggregated over the two
#     contrast strategies (g1, g3).
#     End-of-grace strategy SMD (g1 vs g3 at the first interval after g1's
#     grace): eog_smd_base_max (baseline X only; the practice),
#     eog_smd_post_max (post-baseline L only) and eog_smd_max (max of both).
#     Added columns (deterministic, computed after the RNG-using steps):
#       jackson_cb_est_max  Jackson-type censoring balance (max |SMD|) with
#                           the estimator's weights (censoring_balance(...,
#                           weights = "estimator"); responds to truncation),
#       cumrrb_base_eog     max |cumulative RRB| of baseline X over g1, g3 at
#                           the end-of-grace interval (time = window[1]),
#       cumrrb_base_eog_diff  |cumRRB_g3 - cumRRB_g1| of X at that interval,
#       ltfu_cb_max         max |SMD| of unweighted LTFU censoring balance
#                           (censored = "ltfu_next"), NA without LTFU,
#       ltfu_rrb_max        max |unweighted interval RRB| for LTFU over cells
#                           with n_censored > 0, NA without LTFU.
# Then: estimator performance per scenario (bias, RMSE, MC SD, coverage)
# and the null-reference exceedance rate; the tables and figures of the
# manuscript are made by analysis/02b_simulation_tables_paper.R and
# 02c_simulation_tables_manycov.R.
#
# Runs (CCW_RUN; the four runs that supply the manuscript, see
# analysis/sim_scenarios.R):
#   main  (default) Parts 1-3 ccw_scenarios("main") (Part 1
#      factorial, Part 2 null-reference calibration, Part 3); n and B per
#      scenario (scenario `n`,
#      `n_boot`), 200 replicates (100 for n = 8000). Additionally records the
#      conditional parametric null reference rrb_null() of interval RRB
#      (CCW_N_NULL_DRAWS draws): MaxIB over g1/g3 and all history functions
#      (null_maxib_obs, null_q95, null_p, null_exceed), the same per history
#      function (suffix _<feature>), null_z_max (max |z| over g1/g3 cells) and
#      null_frac_cells_p05 (fraction of g1/g3 cells with p_two < 0.05).
#      See analysis/02b_simulation_tables_paper.R.
#   manycov  ccw_scenarios("manycov"): extended
#      calibration of the null
#      reference (many correlated baseline covariates X2-X9, n = 300 with
#      many covariates, rare artificial censoring, extreme weights) and two
#      detection scenarios with an omitted X2 correlated (rho = 0.5) or not
#      (rho = 0) with X. Run like main (200 replicates, scenario n / n_boot,
#      null reference); the scenario's extra covariates (`extra_feats`) are
#      added to the monitored history functions, so the per-feature null
#      columns _X2 ... _X9 exist only for those scenarios (NA elsewhere). See
#      analysis/02b_simulation_tables_paper.R and
#      02c_simulation_tables_manycov.R.
#   largerep  ccw_scenarios("largerep"): large-replicate
#      calibration of the null
#      reference for RRB and censoring balance under the correct model
#      (n = 300 / 2000 / 8000), early events in g1's grace period (end-of-grace
#      SMD of X non-zero under correct weights through outcome depletion;
#      cumulative RRB of X stays ~0), and the omitted-L detection scenario
#      with truncation at the 95th percentile or near-positivity violation.
#      Run like manycov (bind_rows, null reference); the number of replicates
#      is the scenario's `n_rep` (1000 / 500 / 200). See
#      analysis/02b_simulation_tables_paper.R.
#   variants  ccw_scenarios("variants"): omitted L_t
#      alone (beta2 = beta5
#      = 0), the L_t^2 cells without outcome effect at gamma0 = -4.5, and
#      early_events with 500 bootstrap resamples. Run like largerep (scenario
#      n_rep = 200, null reference).
#   Added columns (all runs; no additional random numbers):
#      * the null reference of the Jackson-type censoring balance (CB)
#        from the same rrb_null() object: MaxIB of |CB| over g1/g3 cells
#        with a defined observed CB and all history functions
#        (null_cb_maxib_obs, null_cb_q95, null_cb_p, null_cb_exceed,
#        null_cb_n_cells) and the per-feature exceedance
#        null_cb_exceed_<feature>;
#      * an individual-bootstrap interval of the cumulative RRB
#        of X for strategies g1, g3 at the end-of-grace interval
#        (eog = window[1]) and at the horizon (h = last interval), computed
#        in the same bootstrap resamples as se_boot (no additional random
#        numbers):
#        cumrrb_X_<g>_<t>_obs, cumrrb_X_<g>_<t>_boot_se and
#        cumrrb_X_<g>_<t>_excl0 (1 if the 2.5-97.5% percentile interval
#        excludes 0), <g> in g1, g3 and <t> in eog, h, and
#        cumrrb_boot_excl0_any (any of the four); NA when B = 0.
# Outputs are suffixed by the run: output/sim/sim_results_<run>.csv,
# sim_performance_<run>.csv, sim_settings_<run>.csv and
# analysis/scenarios_<run>.csv.
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
# Random numbers: one L'Ecuyer-CMRG stream per task, assigned in task order
# (heavy scenarios first, by n x (B + 10)), so changing a scenario's n or B
# (or the number of replicates) reorders the tasks and changes the draws of
# every replicate. Data generation, bootstrap and null
# draws share the replicate's stream, except that rrb_null() is reseeded
# with a seed drawn from that stream after the bootstrap. The seed of a run
# is 20260923 plus its fixed offset (SEED_OFFSET; the offsets of the former
# script version, kept so that the draws do not change).
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
CONTRAST <- c("g3", "g1")            # theta = risk(g3) - risk(g1)
DIAG_STRATEGIES <- c("g1", "g3")     # diagnostics aggregated over these
SEED <- 20260923

out_dir <- file.path("output", "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- ccw_scenarios(RUN)
write.csv(scenario_table(scenarios),
          file.path("analysis", sprintf("scenarios_%s.csv", RUN)),
          row.names = FALSE)

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
contrast_of <- function(est) {
  est <- est[est$time == max(est$time), ]
  r <- stats::setNames(est$risk, as.character(est$strategy))
  unname(r[CONTRAST[1]] - r[CONTRAST[2]])
}

# bootstrap SE of theta_hat; in the same resamples (no additional random
# numbers) the cumulative RRB of X for g1, g3 at the end-of-grace interval
# `t_eog` and at the horizon `t_h`. Returns list(se, cum) with `cum` a data
# frame (strategy, tp in eog/h, boot_se, lcl, ucl) or NULL when B <= 0.
CUM_CELLS <- expand.grid(strategy = DIAG_STRATEGIES, tp = c("eog", "h"),
                         stringsAsFactors = FALSE)
boot_se <- function(d, sc, B, t_eog, t_h) {
  if (B <= 0) return(list(se = NA_real_, cum = NULL))
  idx <- split(seq_len(nrow(d)), d$id)
  ids <- names(idx)
  base_cols <- c("id", "strategy", "time", "censored_next", "censor_possible",
                 "phase", "event_next", "ltfu_next", "X", "L", "L_lag",
                 grep("^X[0-9]+$", names(d), value = TRUE))  # run manycov
  d <- d[, base_cols]
  cum_t <- ifelse(CUM_CELLS$tp == "eog", t_eog, t_h)
  th <- vapply(seq_len(B), function(b) {
    s <- sample(ids, length(ids), replace = TRUE)
    rows <- unlist(idx[s], use.names = FALSE)
    db <- d[rows, ]
    db$id <- rep(seq_along(s), lengths(idx[s]))
    wb <- fit_censoring_model(db, sc$formula, truncate = sc$truncate)
    # deterministic given the resample: no random numbers drawn below
    rcb <- rrb(wb, "X", scope = "cumulative")
    k <- match(paste(CUM_CELLS$strategy, cum_t),
               paste(as.character(rcb$strategy), rcb$time))
    c(contrast_of(estimate_risk(wb)), rcb$rrb[k])
  }, numeric(1 + nrow(CUM_CELLS)))
  th <- matrix(th, nrow = 1 + nrow(CUM_CELLS))
  q <- function(v, p) {
    if (all(is.na(v))) NA_real_ else
      stats::quantile(v, p, na.rm = TRUE, names = FALSE)
  }
  cum <- CUM_CELLS
  cum$boot_se <- apply(th[-1, , drop = FALSE], 1, function(v) {
    if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE)
  })
  cum$lcl <- apply(th[-1, , drop = FALSE], 1, q, 0.025)
  cum$ucl <- apply(th[-1, , drop = FALSE], 1, q, 0.975)
  list(se = stats::sd(th[1, ]), cum = cum)
}

agg_max <- function(x, value, feats) {
  x <- x[as.character(x$strategy) %in% DIAG_STRATEGIES &
           as.character(x$feature) %in% feats & x$n_censored > 0, ]
  if (nrow(x) == 0) return(c(iai = NA, maxib = NA))
  a <- aggregate_balance(x, value = value, by = c("strategy", "feature"))
  c(iai = max(a$iai), maxib = max(a$maxib))
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
  t_eog <- pars$window[1]
  t_h <- max(d$time)
  bs <- boot_se(d, sc, nboot_of(sc), t_eog, t_h)
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
  sb <- strategy_balance(w, raw_feats, reference = CONTRAST[1])
  sb <- sb[as.character(sb$strategy) == CONTRAST[2], ]
  eog <- sb[sb$time == pars$window[1], ]  # first interval after g1's grace
  # end-of-grace SMD: baseline covariate only (the practice), post-baseline
  # L only (extension), and the maximum over both (earlier definition)
  eog_base <- eog$smd[eog$feature == "X"]
  eog_post <- eog$smd[eog$feature == "L"]

  a_raw <- agg_max(r, "rrb", raw_feats)
  a_tr <- agg_max(r, "rrb", tr_feats)
  a_all <- agg_max(r, "rrb", feats_sc)
  a_ks <- agg_max(k, "ks", feats_sc)
  a_cum <- agg_max(rc, "rrb", base_feats)
  mh_iai <- if (nrow(m) > 0) {
    max(tapply(m$mahalanobis, as.character(m$strategy), mean, na.rm = TRUE))
  } else NA

  # additional deterministic diagnostics (no random numbers are drawn here)
  # Jackson-type censoring balance with the estimator's (truncated) weights
  cbe <- censoring_balance(w, feats_sc, weighted = TRUE,
                           weights = "estimator")
  cbe <- cbe[as.character(cbe$strategy) %in% DIAG_STRATEGIES, ]
  # cumulative RRB of baseline X at the end-of-grace interval
  rc_eog <- rc[rc$feature == "X" & rc$time == pars$window[1] &
                 as.character(rc$strategy) %in% DIAG_STRATEGIES, ]
  rc_g <- stats::setNames(rc_eog$rrb, as.character(rc_eog$strategy))
  rc_diff <- if (all(CONTRAST %in% names(rc_g))) {
    abs(rc_g[[CONTRAST[1]]] - rc_g[[CONTRAST[2]]])
  } else NA_real_
  # loss to follow-up: censoring balance and RRB (unweighted; NA if no LTFU)
  max_or_na <- function(v) {
    v <- abs(v[!is.na(v)])
    if (length(v) == 0) NA_real_ else max(v)
  }
  cbl <- censoring_balance(w, feats_sc, weighted = FALSE,
                           censored = "ltfu_next")
  cbl <- cbl[as.character(cbl$strategy) %in% DIAG_STRATEGIES, ]
  rl <- rrb(w, feats_sc, method = "unweighted", censored = "ltfu_next")
  rl <- rl[as.character(rl$strategy) %in% DIAG_STRATEGIES &
             rl$n_censored > 0, ]

  out <- data.frame(
    scenario = sc$name, rep = rep, n = pars$n,
    risk_g1 = risk[["g1"]], risk_g2 = risk[["g2"]], risk_g3 = risk[["g3"]],
    theta_hat = theta_hat, theta_true = truth, se_boot = se,
    abs_bias = abs(theta_hat - truth),
    rrb_raw_iai = a_raw[["iai"]], rrb_raw_maxib = a_raw[["maxib"]],
    rrb_tr_iai = a_tr[["iai"]], rrb_tr_maxib = a_tr[["maxib"]],
    rrb_all_iai = a_all[["iai"]], rrb_all_maxib = a_all[["maxib"]],
    ks_iai = a_ks[["iai"]], ks_max = a_ks[["maxib"]],
    mhb_iai = mh_iai, mhb_max = max(m$mahalanobis, na.rm = TRUE),
    cumrrb_base_iai = a_cum[["iai"]], cumrrb_base_maxib = a_cum[["maxib"]],
    eog_smd_max = max(abs(eog$smd), na.rm = TRUE),
    eog_smd_base_max = max(abs(eog_base), na.rm = TRUE),
    eog_smd_post_max = max(abs(eog_post), na.rm = TRUE),
    jackson_cb_max = max(abs(cb$smd), na.rm = TRUE),
    jackson_cb_mean = mean(abs(cb$smd), na.rm = TRUE),
    ess_ratio_min = min(ws$ess_ratio, na.rm = TRUE),
    w99_max = max(ws$w99, na.rm = TRUE),
    pi_min = min(ws$pi_min, na.rm = TRUE),
    pr_pi_small_max = max(ws$pr_pi_small, na.rm = TRUE),
    # descriptive: strategy-balance of L (g3 vs g1) over follow-up
    strat_smd_L_max = max(abs(sb$smd[sb$feature == "L"]), na.rm = TRUE),
    strat_smd_X_max = max(abs(sb$smd[sb$feature == "X"]), na.rm = TRUE),
    jackson_cb_est_max = max(abs(cbe$smd), na.rm = TRUE),
    cumrrb_base_eog = max_or_na(rc_eog$rrb),
    cumrrb_base_eog_diff = rc_diff,
    ltfu_cb_max = max_or_na(cbl$smd),
    ltfu_rrb_max = max_or_na(rl$rrb)
  )
  # cumulative RRB of X (g1, g3; end of grace and horizon): observed value
  # and individual-bootstrap SE / percentile interval (from boot_se())
  cum_out <- list()
  excl <- logical(0)
  for (k in seq_len(nrow(CUM_CELLS))) {
    g <- CUM_CELLS$strategy[k]
    tp <- CUM_CELLS$tp[k]
    nm <- sprintf("cumrrb_X_%s_%s", g, tp)
    tt <- if (tp == "eog") t_eog else t_h
    v <- rc$rrb[rc$feature == "X" & as.character(rc$strategy) == g &
                  rc$time == tt]
    cum_out[[paste0(nm, "_obs")]] <- if (length(v) == 1) v else NA_real_
    bk <- if (is.null(bs$cum)) NULL else bs$cum[k, ]
    cum_out[[paste0(nm, "_boot_se")]] <- if (is.null(bk)) NA_real_ else
      bk$boot_se
    ex <- if (is.null(bk) || is.na(bk$lcl) || is.na(bk$ucl)) NA_real_ else
      as.numeric(bk$lcl > 0 | bk$ucl < 0)
    cum_out[[paste0(nm, "_excl0")]] <- ex
    excl <- c(excl, ex)
  }
  cum_out$cumrrb_boot_excl0_any <- if (all(is.na(excl))) NA_real_ else
    as.numeric(any(excl == 1, na.rm = TRUE))
  out <- cbind(out, as.data.frame(cum_out))
  cbind(out, null_summary(w, null_seed, feats_sc))
}

# conditional parametric null reference of interval RRB, MaxIB over the
# contrast strategies (g1, g3) and all / each history function; the same for
# the censoring balance (CB)
null_summary <- function(w, seed, feats = all_feats) {
  nr <- rrb_null(w, feats, n_draws = N_NULL_DRAWS, seed = seed)
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
  # censoring balance (CB) from the same draws: MaxIB of |CB| over the g1/g3
  # cells with a defined observed CB (as in rrb_null()'s maxib_cb)
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
  as.data.frame(res)
}

# ---------------------------------------------------------------------------
# truth and task list (one independent L'Ecuyer-CMRG stream per task)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + SEED_OFFSET[[RUN]])
t_start <- Sys.time()
truths <- vapply(scenarios, function(sc) {
  tr <- true_risk_ccw(do.call(ccw_sim_params, sc$params), n_mc = N_MC)
  tr$risk[tr$strategy == CONTRAST[1]] - tr$risk[tr$strategy == CONTRAST[2]]
}, numeric(1))

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
# runs manycov, largerep, variants: per-feature null columns may differ
# between scenarios -> bind_rows (fills NA); run main keeps rbind
res <- if (RUN %in% c("manycov", "largerep", "variants")) {
  as.data.frame(dplyr::bind_rows(ok_list))
} else {
  do.call(rbind, ok_list)
}
res$bias <- res$theta_hat - res$theta_true
res$covered <- abs(res$bias) <= 1.96 * res$se_boot
write.csv(res, file.path(out_dir, sprintf("sim_results_%s.csv", SUFFIX)), row.names = FALSE)

# ---------------------------------------------------------------------------
# performance by scenario
perf <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1], n_rep = nrow(x),
    theta_true = x$theta_true[1], mean_theta_hat = mean(x$theta_hat),
    bias = mean(x$bias), rmse = sqrt(mean(x$bias^2)),
    mc_sd = stats::sd(x$theta_hat), mean_se_boot = mean(x$se_boot),
    coverage = mean(x$covered),
    mean_abs_bias = mean(x$abs_bias),
    rrb_all_iai = mean(x$rrb_all_iai), rrb_all_maxib = mean(x$rrb_all_maxib),
    mhb_max = mean(x$mhb_max), cumrrb_base_maxib = mean(x$cumrrb_base_maxib),
    strat_smd_L_max = mean(x$strat_smd_L_max),
    ess_ratio_min = mean(x$ess_ratio_min), w99_max = mean(x$w99_max),
    jackson_cb_est_max = mean(x$jackson_cb_est_max),
    cumrrb_base_eog = mean(x$cumrrb_base_eog),
    ltfu_cb_max = mean(x$ltfu_cb_max, na.rm = TRUE),
    ltfu_rrb_max = mean(x$ltfu_rrb_max, na.rm = TRUE)
  )
}))
# null-reference summaries; tables and figures of the manuscript in 02b
# and 02c
perf_d <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  data.frame(scenario = x$scenario[1],
             null_exceed_rate = mean(x$null_exceed, na.rm = TRUE),
             null_q95 = mean(x$null_q95, na.rm = TRUE),
             null_frac_cells_p05 = mean(x$null_frac_cells_p05,
                                        na.rm = TRUE))
}))
perf <- merge(perf, perf_d, by = "scenario", sort = FALSE)
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
print(perf[, c("scenario", "theta_true", "bias", "mc_sd", "mean_se_boot",
               "coverage", "rrb_all_maxib", "null_exceed_rate")],
      digits = 3, row.names = FALSE)
