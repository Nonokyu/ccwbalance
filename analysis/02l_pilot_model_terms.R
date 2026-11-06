# 02l_pilot_model_terms.R -----------------------------------------------------
# Pilot: does leaving the censoring-model terms out of the studentised
# maximum in "single-cell fit groups" restore its calibration?
#
# Background (an earlier pilot, 02k_pilot_studentised.R, kept in the
# development repository only; its results are not in the manuscript;
# correct model): for the
# initiate-by-grace-end strategies g1 and g2 the studentised maximum
# max |RRB / null SD| is conservative and the per-function exceedance is
# lowest for X, L, L_sq and L_x_X, the terms of the censoring model
# ~ L + I(L^2) + X + L:X. fit_censoring_model() fits one logistic model per
# strategy x phase over the censor_possible rows; for g1 and g2 the
# grace_end group lies in a single time cell, so the score equations of
# the fit (in the data and in every refitted null draw) act on exactly one
# diagnostic cell. In g3 the no_init group pools several cells.
#
# The simulation, the censoring-model fit and the rrb_null() call are those
# of 02k (same helpers, same random-number use), so raw, stud, cb_cens and
# cb_cens_stud and the per-function raw/stud exceedances reproduce 02k for
# the same task seeds (see Random numbers).
#
# Definitions (per replicate and strategy g; entries are the strategy-g
# cell x history-function pairs of rrb_null()'s cells with a non-NA
# observed value; studentised statistics also drop entries whose null SD is
# NA or <= 0; for every statistic obs = max over the used entries of the
# |value| (divided by the null SD if studentised), the null maxima are the
# same maxima per draw, q95 = their 95th percentile, p = (1 + #{null >=
# obs}) / (1 + #{non-NA null}), exceed = obs > q95, entries = number of
# entries used):
#   raw                interval RRB, all entries (= 02k maxib_obs)
#   stud               RRB / null_sd, all entries (= 02k stud_obs)
#   stud_rule          as stud, but in every "single-fit-group cell" (below)
#                      the model-term functions are dropped
#   stud_nomodel       as stud, model-term functions dropped in all cells
#   raw_rule           as raw, with the exclusion of stud_rule
#   cb_cens            censorable-row censoring balance, raw (= 02k cb_cens)
#   cb_cens_stud       cb_cens / cb_cens_null_sd (= 02k cb_cens_stud)
#   cb_cens_stud_rule  as cb_cens_stud with the exclusion of stud_rule
#   raw_ldiff01, raw_ldiff10, stud_ldiff01, stud_ldiff10
#                      raw / stud recomputed with the L_diff entries of the
#                      observed RRB and of the null draws multiplied by 0.1
#                      / 10 (post hoc; RRB and its draws are linear in the
#                      feature, the null SD is recomputed from the scaled
#                      draws). The studentised maximum must be invariant:
#                      the script stops if any stud_ldiff* differs from stud
#                      (relative 1e-8; exceed must be equal).
# Single-fit-group cell: the censoring-model fit groups are the levels of
# the `by` columns of fit_censoring_model() (strategy x phase) over the
# censor_possible rows. A diagnostic cell (strategy x time) is flagged when
# every fit group with censorable rows in the cell has all of its
# censorable rows in that one cell (determined from the data of each
# replicate, not by strategy). Model-term functions of a cell: the terms of
# the censoring formula of its fit group(s) (same lookup as
# fit_censoring_model() if the formula is a list), mapped by
#   X -> X, L -> L, I(L^2) -> L_sq, L:X (or X:L) -> L_x_X,
#   any other term that is itself a monitored function name -> itself
#   (e.g. X2, ..., X9 of the manycov run);
# any other term is an error, as is a mapped function that is not
# monitored.
#
# Outputs (<CCW_SIM_DIR>, suffix "_<CCW_OUT_SUFFIX>" if set):
#   pilot_model_terms_reps.csv     one row per scenario x rep x strategy:
#     null_seed, n, <stat>_{obs,q95,p,exceed,entries} for the statistics
#     above, raw_exceed_<f> / stud_exceed_<f> per history function f (02k's
#     per-feature sets), n_flagged (flagged diagnostic cells),
#     flagged_times, model_feats, check_maxdiff (max abs difference of
#     raw / stud / cb_cens / cb_cens_stud from rrb_null()'s own tables; the
#     script stops if > 1e-10), scale_maxdiff, runtime_sec
#   pilot_model_terms.csv          scenario x strategy x statistic x feature
#     ("(all)" for the statistics above; per-function rows for raw/stud):
#     n_rep (successful replicates), n_failed (failed tasks of the
#     scenario), n_eval (non-NA exceed), rate, se (binomial), mean_obs,
#     mean_q95, mean_entries
#   pilot_model_terms_mcnemar.csv  scenario x strategy x comparison (a vs
#     b: stud_rule vs raw, stud_rule vs stud, cb_cens_stud_rule vs
#     cb_cens_stud): n_pairs, n_both, n_a_only, n_b_only, rate_a, rate_b,
#     p_mcnemar (exact two-sided binomial test of the discordant pairs; 1
#     if there are none)
#
# Configuration (environment variables):
#   CCW_N_REPS_CORRECT  replicates of scenario "correct"             [1000]
#   CCW_N_REPS_OMIT     replicates of every other scenario            [300]
#   CCW_N_NULL_DRAWS    null-reference draws per replicate            [500]
#   CCW_SCENARIOS       comma-separated scenario names of run main
#                         [correct,omit_lag_cw_on,omit_lag_cs_on]
#   CCW_N_SUBJ          individuals per replicate (unless the scenario
#                       sets n)                                      [2000]
#   CCW_CORES           parallel workers (parallel::mclapply)
#                                                      [detectCores() - 1]
#   CCW_SIM_DIR         output directory                      [output/sim]
#   CCW_OUT_SUFFIX      output suffix                                   [""]
#   CCW_SMOKE           1 = 3 replicates per scenario, 20 draws, suffix
#                       "smoke"                                         [0]
#   CCW_REPRO           1 = reproduce 02k: the three default scenarios,
#                       CCW_REPRO_BLOCK replicates each, suffix "repro";
#                       with 100 replicates, 500 draws and n = 2000 the
#                       summary is compared with CCW_REPRO_REF            [0]
#   CCW_REPRO_BLOCK     replicates per scenario of 02k's seed layout (02k's
#                       CCW_N_REPS)                                    [100]
#   CCW_REPRO_REF       02k summary (not distributed here; the
#                       comparison is skipped when the file is missing)
#                         [output/sim/pilot_studentised_summary.csv]
#
# Random numbers (one L'Ecuyer-CMRG stream per scenario x replicate task):
# RNGkind("L'Ecuyer-CMRG"); set.seed(20261102). Scenarios are indexed in the
# canonical order k = 1..K: correct, omit_lag_cw_on, omit_lag_cs_on, then
# the other scenarios of ccw_scenarios("main") in their order (K = 20).
# Block 1: seeds1 <- sample.int(.Machine$integer.max, B * K) with B =
# CCW_REPRO_BLOCK; replicate r <= B of scenario k gets seeds1[(k - 1) * B +
# r]. This is 02k's layout (tasks expand.grid(rep = 1:B, s), seeds
# sample.int(.Machine$integer.max, nrow(tasks))): sample.int() draws its
# values sequentially, so the first 3B values equal 02k's draw for its
# three default scenarios and replicate r of those scenarios uses 02k's
# seed for replicate r of a CCW_N_REPS = B run (for any subset or order of
# CCW_SCENARIOS). Block 2 (replicates r > B), drawn next from the same
# stream: seeds2 <- sample.int(.Machine$integer.max, (R - B) * K) with R the
# largest replicate count; replicate r of scenario k gets seeds2[(r - B - 1)
# * K + k] (replicate-major, so a seed does not depend on R). Inside a task
# set.seed(seed) is followed by exactly 02k's calls (simulate_ccw(),
# make_history_features(), fit_censoring_model(), the rrb_null() seed drawn
# by sample.int(), rrb_null()); nothing below draws random numbers.
#
# Run from the project root:
#   CCW_CORES=11 Rscript analysis/02l_pilot_model_terms.R
#   CCW_SMOKE=1 CCW_CORES=2 Rscript analysis/02l_pilot_model_terms.R
#   CCW_REPRO=1 CCW_CORES=11 Rscript analysis/02l_pilot_model_terms.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) {
  as.numeric(Sys.getenv(x, as.character(default)))
}
SMOKE <- Sys.getenv("CCW_SMOKE", "0") == "1"
REPRO <- Sys.getenv("CCW_REPRO", "0") == "1"
if (SMOKE && REPRO) stop("CCW_SMOKE and CCW_REPRO cannot both be set.")
DEFAULT_SC <- c("correct", "omit_lag_cw_on", "omit_lag_cs_on")
BLOCK <- env_num("CCW_REPRO_BLOCK", 100)
N_REPS_CORRECT <- env_num("CCW_N_REPS_CORRECT", 1000)
N_REPS_OMIT <- env_num("CCW_N_REPS_OMIT", 300)
N_NULL_DRAWS <- env_num("CCW_N_NULL_DRAWS", 500)
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
N_SUBJ <- env_num("CCW_N_SUBJ", 2000)
SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", "")
sc_names <- trimws(strsplit(
  Sys.getenv("CCW_SCENARIOS", paste(DEFAULT_SC, collapse = ",")), ",")[[1]])
if (SMOKE) {
  N_REPS_CORRECT <- 3
  N_REPS_OMIT <- 3
  N_NULL_DRAWS <- 20
  if (SUFFIX == "") SUFFIX <- "smoke"
}
if (REPRO) {
  if (!setequal(sc_names, DEFAULT_SC)) {
    stop("CCW_REPRO=1 uses the scenarios ", toString(DEFAULT_SC), ".")
  }
  N_REPS_CORRECT <- BLOCK
  N_REPS_OMIT <- BLOCK
  if (SUFFIX == "") SUFFIX <- "repro"
}
REPRO_REF <- Sys.getenv("CCW_REPRO_REF", file.path(
  "output", "sim", "pilot_studentised_summary.csv"))
out_dir <- Sys.getenv("CCW_SIM_DIR", file.path("output", "sim"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(stem) {
  file.path(out_dir, paste0(stem, if (SUFFIX != "") paste0("_", SUFFIX),
                            ".csv"))
}
STRATEGIES <- c("g1", "g2", "g3")

all_scen <- ccw_scenarios("main")
names(all_scen) <- vapply(all_scen, function(sc) sc$name, character(1))
miss <- setdiff(sc_names, names(all_scen))
if (length(miss) > 0) stop("Unknown scenario(s): ", toString(miss))
scenarios <- all_scen[sc_names]
canon <- c(DEFAULT_SC, setdiff(names(all_scen), DEFAULT_SC))
n_reps_of <- function(nm) if (nm == "correct") N_REPS_CORRECT else N_REPS_OMIT

# monitored history functions (as in 02_simulation_study.R and 02k)
all_feats <- c("X", "L", "L_sq", "L_lagprod", "L_diff", "L_cummean",
               "L_gt0", "L_x_X")
feats_of <- function(sc) c(all_feats, as.character(sc$extra_feats))

# model terms -> history functions ---------------------------------------------
TERM_MAP <- c(X = "X", L = "L", "I(L^2)" = "L_sq", "L:X" = "L_x_X",
              "X:L" = "L_x_X")
model_feats <- function(formula, feats) {
  labs <- labels(stats::terms(formula))
  out <- vapply(labs, function(lab) {
    f <- if (lab %in% names(TERM_MAP)) {
      TERM_MAP[[lab]]
    } else if (lab %in% feats) {
      lab
    } else {
      stop("Censoring-model term '", lab, "' cannot be mapped to a ",
           "history function.", call. = FALSE)
    }
    if (!f %in% feats) {
      stop("Censoring-model term '", lab, "' maps to '", f, "', which is ",
           "not monitored.", call. = FALSE)
    }
    f
  }, character(1))
  unique(unname(out))
}

# formula of a fit group (lookup of fit_censoring_model())
pick_formula <- function(formula, label, last) {
  if (!is.list(formula)) return(formula)
  for (nm in c(label, last, "default")) {
    if (!is.null(formula[[nm]])) return(formula[[nm]])
  }
  stop("No formula for censoring-model group '", label, "'.", call. = FALSE)
}

# per strategy x time cell with censorable rows: fit groups present, flag
# "every fit group of the cell lies entirely in this cell", model functions
rule_cells <- function(w, feats) {
  wt <- attr(w, "weighting")
  by <- wt$by
  poss <- as.logical(w$censor_possible) %in% TRUE
  cp <- as.data.frame(w[poss, c("strategy", "time", by), drop = FALSE])
  grp <- if (length(by) > 0) {
    as.character(interaction(cp[by], drop = TRUE))
  } else {
    rep("all", nrow(cp))
  }
  key <- paste(as.character(cp$strategy), cp$time, sep = "\r")
  g_ncell <- tapply(key, grp, function(k) length(unique(k)))
  g_feats <- lapply(split(seq_len(nrow(cp)), grp), function(ix) {
    last <- if (length(by) > 0) as.character(cp[[by[length(by)]]][ix[1]])
    model_feats(pick_formula(wt$formula, grp[ix[1]], last), feats)
  })
  ukey <- unique(key)
  do.call(rbind, lapply(ukey, function(k) {
    ix <- which(key == k)
    gs <- sort(unique(grp[ix]))
    data.frame(strategy = as.character(cp$strategy[ix[1]]),
               time = cp$time[ix[1]], groups = paste(gs, collapse = "+"),
               n_censorable = length(ix),
               flagged = all(g_ncell[gs] == 1),
               model_feats = paste(unique(unlist(g_feats[gs])),
                                   collapse = ","))
  }))
}

# maximum statistic over the selected entries (logic of .null_maxib())
max_stat <- function(obs, sims, sel, scale = NULL) {
  sel <- sel & !is.na(obs)
  if (!is.null(scale)) {
    ok <- !is.na(scale) & scale > 0
    sel <- sel & ok
    sc <- ifelse(ok, scale, 1)
    obs <- obs / sc
    sims <- sims / sc
  }
  if (!any(sel)) {
    return(c(obs = NA_real_, q95 = NA_real_, p = NA_real_,
             exceed = NA_real_, entries = 0))
  }
  m_obs <- max(abs(obs[sel]))
  m_null <- apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  q95 <- stats::quantile(m_null, 0.95, na.rm = TRUE, names = FALSE)
  c(obs = m_obs, q95 = q95,
    p = (1 + sum(m_null >= m_obs, na.rm = TRUE)) / (1 + sum(!is.na(m_null))),
    exceed = as.numeric(m_obs > q95), entries = sum(sel))
}

# row of a maxib table (as 02k's pick())
pick <- function(tab, g, f, obs_col) {
  r <- tab[tab$strategy == g & tab$feature == f, , drop = FALSE]
  if (nrow(r) == 0 || r$n_cells == 0) {
    return(c(obs = NA_real_, q95 = NA_real_, exceed = NA_real_))
  }
  c(obs = r[[obs_col]], q95 = r$null_q95,
    exceed = as.numeric(r[[obs_col]] > r$null_q95))
}

diff_max <- function(a, b) {
  d <- abs(a - b)
  d[is.na(a) & is.na(b)] <- 0
  d[xor(is.na(a), is.na(b))] <- Inf
  max(d)
}

one_replicate <- function(sc, rep) {
  t0 <- Sys.time()
  # --- exactly 02k's path ---
  pars <- do.call(ccw_sim_params, sc$params)
  if (is.null(sc$params[["n"]])) pars$n <- N_SUBJ
  d <- simulate_ccw(pars)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, sc$formula, truncate = sc$truncate)
  null_seed <- sample.int(.Machine$integer.max, 1)
  feats <- feats_of(sc)
  nr <- rrb_null(w, feats, n_draws = N_NULL_DRAWS, seed = null_seed)
  # --- post hoc statistics ---
  cells <- nr$cells
  nc <- nrow(cells)
  sims <- matrix(nr$draws$rrb, nc, N_NULL_DRAWS)
  sims_cbc <- matrix(nr$draws$cb_cens, nc, N_NULL_DRAWS)
  stopifnot(identical(as.character(nr$draws$feature[seq_len(nc)]),
                      as.character(cells$feature)),
            identical(nr$draws$time[seq_len(nc)], cells$time))
  rc <- rule_cells(w, feats)
  ckey <- paste(as.character(cells$strategy), cells$time, sep = "\r")
  m <- match(ckey, paste(rc$strategy, rc$time, sep = "\r"))
  if (anyNA(m)) stop("Diagnostic cell without censorable rows.")
  flag <- rc$flagged[m]
  is_model <- mapply(function(f, mf) f %in% strsplit(mf, ",")[[1]],
                     as.character(cells$feature), rc$model_feats[m])
  drop_rule <- flag & is_model
  feat <- as.character(cells$feature)
  ld <- feat == "L_diff"
  rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rows <- lapply(STRATEGIES, function(g) {
    sg <- as.character(cells$strategy) == g
    st <- list(
      raw = max_stat(cells$rrb_obs, sims, sg),
      stud = max_stat(cells$rrb_obs, sims, sg, cells$null_sd),
      stud_rule = max_stat(cells$rrb_obs, sims, sg & !drop_rule,
                           cells$null_sd),
      stud_nomodel = max_stat(cells$rrb_obs, sims, sg & !is_model,
                              cells$null_sd),
      raw_rule = max_stat(cells$rrb_obs, sims, sg & !drop_rule),
      cb_cens = max_stat(cells$cb_cens_obs, sims_cbc, sg),
      cb_cens_stud = max_stat(cells$cb_cens_obs, sims_cbc, sg,
                              cells$cb_cens_null_sd),
      cb_cens_stud_rule = max_stat(cells$cb_cens_obs, sims_cbc,
                                   sg & !drop_rule, cells$cb_cens_null_sd)
    )
    scales <- c(ldiff01 = 0.1, ldiff10 = 10)
    for (nm in names(scales)) {
      fac <- ifelse(ld, scales[[nm]], 1)
      o_k <- cells$rrb_obs * fac
      s_k <- sims * fac
      sd_k <- apply(s_k, 1, stats::sd, na.rm = TRUE)
      st[[paste0("raw_", nm)]] <- max_stat(o_k, s_k, sg)
      st[[paste0("stud_", nm)]] <- max_stat(o_k, s_k, sg, sd_k)
    }
    # checks: against rrb_null()'s own tables and scale invariance
    chk <- c(
      st$raw[c("obs", "q95", "exceed")] -
        pick(nr$maxib, g, "(all)", "maxib_obs"),
      st$stud[c("obs", "q95", "exceed")] -
        pick(nr$maxib_stud, g, "(all)", "stud_obs"),
      st$cb_cens[c("obs", "q95", "exceed")] -
        pick(nr$maxib_cb_cens, g, "(all)", "maxib_obs"),
      st$cb_cens_stud[c("obs", "q95", "exceed")] -
        pick(nr$maxib_cb_cens_stud, g, "(all)", "stud_obs")
    )
    chk[is.na(chk)] <- 0
    rel <- function(a, b) {
      v <- c(abs(a[c("obs", "q95")] - b[c("obs", "q95")]) /
               pmax(abs(b[c("obs", "q95")]), 1e-300),
             abs(a[c("exceed", "entries")] - b[c("exceed", "entries")]) *
               1e10)
      v[is.na(v)] <- 0
      max(v)
    }
    out <- data.frame(scenario = sc$name, rep = rep, strategy = g,
                      n = pars$n, null_seed = null_seed)
    for (s in names(st)) {
      for (v in names(st[[s]])) out[[paste0(s, "_", v)]] <- st[[s]][[v]]
    }
    for (f in feats) {
      out[[paste0("raw_exceed_", f)]] <-
        pick(nr$maxib, g, f, "maxib_obs")[["exceed"]]
      out[[paste0("stud_exceed_", f)]] <-
        pick(nr$maxib_stud, g, f, "stud_obs")[["exceed"]]
    }
    rg <- rc[rc$strategy == g & rc$flagged, , drop = FALSE]
    rg <- rg[paste(rg$strategy, rg$time, sep = "\r") %in% ckey, ,
             drop = FALSE]
    out$n_flagged <- nrow(rg)
    out$flagged_times <- paste(rg$time, collapse = ",")
    out$model_feats <- paste(unique(rc$model_feats[rc$strategy == g]),
                             collapse = "|")
    out$check_maxdiff <- max(abs(chk))
    out$scale_maxdiff <- max(rel(st$stud_ldiff01, st$stud),
                             rel(st$stud_ldiff10, st$stud))
    out$runtime_sec <- rt
    out
  })
  do.call(rbind, rows)
}

# tasks and seeds (see Random numbers) ---------------------------------------
tasks <- do.call(rbind, lapply(sc_names, function(nm) {
  data.frame(scenario = nm, k = match(nm, canon),
             rep = seq_len(n_reps_of(nm)))
}))
K <- length(canon)
R_MAX <- max(tasks$rep)
RNGkind("L'Ecuyer-CMRG")
set.seed(20261102)
seeds1 <- sample.int(.Machine$integer.max, BLOCK * K)
seeds2 <- if (R_MAX > BLOCK) {
  sample.int(.Machine$integer.max, (R_MAX - BLOCK) * K)
} else {
  integer(0)
}
tasks$seed <- ifelse(tasks$rep <= BLOCK,
                     seeds1[(tasks$k - 1) * BLOCK + pmin(tasks$rep, BLOCK)],
                     seeds2[pmax((tasks$rep - BLOCK - 1) * K + tasks$k, 1)])

message(sprintf(paste("Pilot (model terms): %d tasks (%s), %d null draws,",
                      "n = %d, %d cores, block %d%s"),
                nrow(tasks),
                paste(sprintf("%s x %d", sc_names,
                              vapply(sc_names, n_reps_of, numeric(1))),
                      collapse = ", "),
                N_NULL_DRAWS, N_SUBJ, N_CORES, BLOCK,
                if (REPRO) ", REPRO" else ""))
t_start <- Sys.time()
res_list <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
  set.seed(tasks$seed[i])
  sc <- scenarios[[tasks$scenario[i]]]
  tryCatch(one_replicate(sc, tasks$rep[i]),
           error = function(e) {
             data.frame(scenario = sc$name, rep = tasks$rep[i],
                        error = conditionMessage(e))
           })
}, mc.cores = N_CORES, mc.set.seed = TRUE, mc.preschedule = FALSE)
wall <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

is_fail <- vapply(res_list, function(x) "error" %in% names(x), logical(1))
fail_df <- if (any(is_fail)) do.call(rbind, res_list[is_fail]) else
  data.frame(scenario = character(0), rep = integer(0), error = character(0))
if (nrow(fail_df) > 0) {
  print(fail_df)
  warning(nrow(fail_df), " task(s) failed; first: ", fail_df$error[1])
}
res <- do.call(rbind, res_list[!is_fail])
if (is.null(res)) stop("All tasks failed.")
if (max(res$check_maxdiff) > 1e-10) {
  stop("Statistics differ from rrb_null()'s tables (max diff ",
       format(max(res$check_maxdiff)), ").")
}
if (max(res$scale_maxdiff) > 1e-8) {
  stop("Studentised maximum not invariant to the L_diff scale (max rel. ",
       "diff ", format(max(res$scale_maxdiff)), ").")
}
message(sprintf(paste("Checks passed: max diff vs rrb_null() tables %.2g;",
                      "max rel. diff of stud under L_diff scaling %.2g."),
                max(res$check_maxdiff), max(res$scale_maxdiff)))
write.csv(res, out_file("pilot_model_terms_reps"), row.names = FALSE)

# summary ----------------------------------------------------------------------
STATS <- c("raw", "stud", "stud_rule", "stud_nomodel", "raw_rule",
           "cb_cens", "cb_cens_stud", "cb_cens_stud_rule", "raw_ldiff01",
           "raw_ldiff10", "stud_ldiff01", "stud_ldiff10")
grp <- unique(res[, c("scenario", "strategy")])
summ <- do.call(rbind, lapply(seq_len(nrow(grp)), function(j) {
  x <- res[res$scenario == grp$scenario[j] &
             res$strategy == grp$strategy[j], , drop = FALSE]
  n_failed <- sum(fail_df$scenario == grp$scenario[j])
  one <- function(stat, feature, ex, obs = NULL, q95 = NULL, ent = NULL) {
    v <- ex[!is.na(ex)]
    p <- if (length(v) > 0) mean(v) else NA_real_
    mn <- function(z) if (is.null(z)) NA_real_ else mean(z, na.rm = TRUE)
    data.frame(scenario = grp$scenario[j], strategy = grp$strategy[j],
               statistic = stat, feature = feature, n_rep = nrow(x),
               n_failed = n_failed, n_eval = length(v), rate = p,
               se = if (length(v) > 0) sqrt(p * (1 - p) / length(v)) else NA,
               mean_obs = mn(obs), mean_q95 = mn(q95), mean_entries = mn(ent))
  }
  a <- lapply(STATS, function(s) {
    one(s, "(all)", x[[paste0(s, "_exceed")]], x[[paste0(s, "_obs")]],
        x[[paste0(s, "_q95")]], x[[paste0(s, "_entries")]])
  })
  fcols <- grep("^(raw|stud)_exceed_", names(x), value = TRUE)
  b <- lapply(fcols, function(cl) {
    one(sub("_exceed_.*$", "", cl), sub("^(raw|stud)_exceed_", "", cl),
        x[[cl]])
  })
  do.call(rbind, c(a, b))
}))
summ <- summ[order(match(summ$scenario, sc_names),
                   match(summ$strategy, STRATEGIES),
                   summ$feature != "(all)",
                   match(summ$statistic, STATS)), ]
write.csv(summ, out_file("pilot_model_terms"), row.names = FALSE)

# paired comparisons (exact McNemar) ------------------------------------------
COMP <- list(c("stud_rule", "raw"), c("stud_rule", "stud"),
             c("cb_cens_stud_rule", "cb_cens_stud"))
mcn <- do.call(rbind, lapply(seq_len(nrow(grp)), function(j) {
  x <- res[res$scenario == grp$scenario[j] &
             res$strategy == grp$strategy[j], , drop = FALSE]
  do.call(rbind, lapply(COMP, function(ab) {
    ea <- x[[paste0(ab[1], "_exceed")]]
    eb <- x[[paste0(ab[2], "_exceed")]]
    ok <- !is.na(ea) & !is.na(eb)
    ea <- ea[ok]
    eb <- eb[ok]
    n10 <- sum(ea == 1 & eb == 0)
    n01 <- sum(ea == 0 & eb == 1)
    data.frame(scenario = grp$scenario[j], strategy = grp$strategy[j],
               a = ab[1], b = ab[2], n_pairs = sum(ok),
               n_both = sum(ea == 1 & eb == 1), n_a_only = n10,
               n_b_only = n01, rate_a = mean(ea), rate_b = mean(eb),
               p_mcnemar = if (n10 + n01 > 0) {
                 stats::binom.test(n10, n10 + n01, 0.5)$p.value
               } else {
                 1
               })
  }))
}))
write.csv(mcn, out_file("pilot_model_terms_mcnemar"), row.names = FALSE)

# printout ---------------------------------------------------------------------
cat("\nFlagged cells and model-term functions (distinct values over",
    "replicates):\n")
fl <- unique(res[, c("scenario", "strategy", "n_flagged", "flagged_times",
                     "model_feats")])
print(fl, row.names = FALSE)
cat("\nExceedance rates (obs > null q95), all functions:\n")
wide <- reshape(summ[summ$feature == "(all)",
                     c("scenario", "strategy", "statistic", "rate")],
                idvar = c("scenario", "strategy"), timevar = "statistic",
                direction = "wide")
names(wide) <- sub("^rate\\.", "", names(wide))
print(wide, digits = 3, row.names = FALSE)
cat("\nExact McNemar comparisons:\n")
print(mcn, digits = 3, row.names = FALSE)

# reproduction of 02k ---------------------------------------------------------
if (REPRO) {
  if (BLOCK == 100 && N_NULL_DRAWS == 500 && N_SUBJ == 2000 &&
        nrow(fail_df) == 0 && file.exists(REPRO_REF)) {
    ref <- read.csv(REPRO_REF)
    get <- function(s, f, col) {
      r <- summ[summ$statistic == s & summ$feature == f, ]
      r[match(paste(ref$scenario, ref$strategy),
              paste(r$scenario, r$strategy)), col]
    }
    map <- list(
      exceed = get("raw", "(all)", "rate"),
      stud_exceed = get("stud", "(all)", "rate"),
      cb_cens_exceed = get("cb_cens", "(all)", "rate"),
      cb_cens_stud_exceed = get("cb_cens_stud", "(all)", "rate"),
      exceed_L = get("raw", "L", "rate"),
      exceed_L_diff = get("raw", "L_diff", "rate"),
      stud_exceed_L = get("stud", "L", "rate"),
      stud_exceed_L_diff = get("stud", "L_diff", "rate"),
      maxib_obs = get("raw", "(all)", "mean_obs"),
      null_q95 = get("raw", "(all)", "mean_q95"),
      stud_obs = get("stud", "(all)", "mean_obs"),
      stud_q95 = get("stud", "(all)", "mean_q95"),
      cb_cens_obs = get("cb_cens", "(all)", "mean_obs"),
      cb_cens_q95 = get("cb_cens", "(all)", "mean_q95"),
      cb_cens_stud_obs = get("cb_cens_stud", "(all)", "mean_obs"),
      cb_cens_stud_q95 = get("cb_cens_stud", "(all)", "mean_q95")
    )
    dif <- vapply(names(map), function(cl) diff_max(map[[cl]], ref[[cl]]),
                  numeric(1))
    cat("\nReproduction of", REPRO_REF, "(max abs difference per column):\n")
    print(signif(dif, 3))
    if (max(dif) > 1e-8) stop("02k summary NOT reproduced.")
    message("02k summary reproduced (max abs difference ",
            format(max(dif), digits = 3), ").")
  } else {
    message("REPRO: comparison with ", REPRO_REF, " skipped (needs 100 ",
            "replicates, 500 draws, n = 2000, no failures).")
  }
}

message(sprintf(paste("Done: %d of %d tasks, wall time %.2f min, mean",
                      "per-replicate runtime %.1f s."),
                nrow(res) / length(STRATEGIES), nrow(tasks), wall,
                mean(res$runtime_sec[!duplicated(res[, c("scenario",
                                                         "rep")])])))
message("Written: ", out_file("pilot_model_terms_reps"), ", ",
        out_file("pilot_model_terms"), ", ",
        out_file("pilot_model_terms_mcnemar"))
