# 02l_pilot_model_terms.R -----------------------------------------------------
# Per-strategy null exceedance of the raw and studentised maxima, with and
# without the censoring-model terms in cells that form a single fit group
# (Appendix S2). Writes output/sim/pilot_model_terms{,_reps,_mcnemar}.csv.
# Env vars: CCW_N_REPS_CORRECT (1000), CCW_N_REPS_OMIT (300),
#   CCW_N_NULL_DRAWS (500), CCW_SCENARIOS, CCW_N_SUBJ (2000), CCW_CORES,
#   CCW_SIM_DIR, CCW_OUT_SUFFIX, CCW_SMOKE (1 = quick test run).
# Run: Rscript analysis/02l_pilot_model_terms.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) {
  as.numeric(Sys.getenv(x, as.character(default)))
}
SMOKE <- Sys.getenv("CCW_SMOKE", "0") == "1"
DEFAULT_SC <- c("correct", "omit_lag_cw_on", "omit_lag_cs_on")
BLOCK <- 100 # replicates per scenario in the first seed block
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

# monitored history functions (as in 02_simulation_study.R)
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

# per strategy x time cell: flag cells whose fit groups lie entirely in the
# cell, and the model-term functions of the cell
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

# row of a maxib table
pick <- function(tab, g, f, obs_col) {
  r <- tab[tab$strategy == g & tab$feature == f, , drop = FALSE]
  if (nrow(r) == 0 || r$n_cells == 0) {
    return(c(obs = NA_real_, q95 = NA_real_, exceed = NA_real_))
  }
  c(obs = r[[obs_col]], q95 = r$null_q95,
    exceed = as.numeric(r[[obs_col]] > r$null_q95))
}

one_replicate <- function(sc, rep) {
  t0 <- Sys.time()
  # --- statistics returned by rrb_null() ---
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

# tasks and seeds: block 1 one seed per replicate up to BLOCK,
# block 2 (replicate-major) for the remaining replicates ---------------------
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
                N_NULL_DRAWS, N_SUBJ, N_CORES, BLOCK, ""))
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

message(sprintf(paste("Done: %d of %d tasks, wall time %.2f min, mean",
                      "per-replicate runtime %.1f s."),
                nrow(res) / length(STRATEGIES), nrow(tasks), wall,
                mean(res$runtime_sec[!duplicated(res[, c("scenario",
                                                         "rep")])])))
message("Written: ", out_file("pilot_model_terms_reps"), ", ",
        out_file("pilot_model_terms"), ", ",
        out_file("pilot_model_terms_mcnemar"))
