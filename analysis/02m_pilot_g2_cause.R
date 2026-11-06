# 02m_pilot_g2_cause.R ---------------------------------------------------------
# Per-strategy decomposition of the null reference of the interval RRB
# (Supporting Table S1): exceedance of the raw and studentised maxima under
# fitted vs true censoring probabilities, with and without refitting, on the
# whole risk set ("all") or the censorable rows ("cens"); test pairs S1-S5.
# Writes output/sim/pilot_g2_cause{,_summary,_mcnemar}.csv.
# Env vars: CCW_N_REPS_{CORRECT,N300,POWER,SCALE}, CCW_N_NULL_DRAWS (500),
#   CCW_SCENARIOS, CCW_CORES, CCW_OUT_SUFFIX, CCW_CHECK_PKG, CCW_SMOKE.
# Run: Rscript analysis/02m_pilot_g2_cause.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) {
  as.numeric(Sys.getenv(x, as.character(default)))
}
SMOKE <- Sys.getenv("CCW_SMOKE", "0") == "1"
SCALE <- env_num("CCW_N_REPS_SCALE", 1)
N_REPS <- c(
  correct = env_num("CCW_N_REPS_CORRECT", 1000),
  omit_lag_cs_on = env_num("CCW_N_REPS_POWER", 300),
  omit_L_cs_on = env_num("CCW_N_REPS_POWER", 300),
  correct_n300 = env_num("CCW_N_REPS_N300", 500)
)
N_REPS <- round(N_REPS * SCALE)
N_DRAWS <- env_num("CCW_N_NULL_DRAWS", 500)
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", if (SMOKE) "smoke" else "")
CHECK_PKG <- Sys.getenv("CCW_CHECK_PKG", if (SMOKE) "1" else "0") == "1"
if (SMOKE) {
  N_REPS[] <- 3
  N_DRAWS <- 20
}
SEED <- 20261103
PI_FLOOR <- 1e-6
STRATEGIES <- c("g1", "g2", "g3")
VARIANTS <- c("all", "cens")
SCHEMES <- c("fit_refit", "fit_norefit", "true_refit", "true_norefit")
PAIRS <- data.frame(
  pair = c("S1", "S2", "S3", "S4", "S5"),
  scheme = c(SCHEMES, "true_norefit"),
  observed = c("fit", "fit", "fit", "fit", "true"),
  stringsAsFactors = FALSE
)
out_dir <- file.path("output", "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x) {
  file.path(out_dir, paste0(x, if (SUFFIX != "") paste0("_", SUFFIX), ".csv"))
}

all_feats <- c("X", "L", "L_sq", "L_lagprod", "L_diff", "L_cummean", "L_gt0",
               "L_x_X")
main_sc <- ccw_scenarios("main")
names(main_sc) <- vapply(main_sc, function(sc) sc$name, character(1))
scenarios <- list(  # n = 2000 first (task order)
  correct = list(scenario = paper_scenario("P1-01"), calib = TRUE),
  omit_lag_cs_on = list(scenario = main_sc[["omit_lag_cs_on"]],
                        calib = FALSE),
  omit_L_cs_on = list(scenario = main_sc[["omit_L_cs_on"]], calib = FALSE),
  correct_n300 = list(scenario = paper_scenario("P2-01"), calib = TRUE)
)
sc_names <- trimws(strsplit(
  Sys.getenv("CCW_SCENARIOS", paste(names(scenarios), collapse = ",")),
  ",")[[1]])
miss <- setdiff(sc_names, names(scenarios))
if (length(miss) > 0) stop("Unknown scenario(s): ", toString(miss))
sc_names <- intersect(names(scenarios), sc_names)

# true P(remain uncensored | H_t) from the DGM (simulate_ccw())
true_pi <- function(d, pars) {
  b <- pars$beta
  lp <- b[["b0"]] + b[["b1"]] * d$L + b[["b2"]] * d$L^2 + b[["b3"]] * d$X +
    b[["b5"]] * d$L * d$X + b[["b6"]] * d$L_lag
  if (pars$n_extra > 0) {
    xe <- as.matrix(d[, paste0("X", seq_len(pars$n_extra) + 1), drop = FALSE])
    lp <- lp + as.vector(xe %*% pars$beta_extra)
  }
  p_init <- stats::plogis(lp)
  ph <- as.character(d$phase)
  ifelse(ph == "no_init", 1 - p_init, ifelse(ph == "grace_end", p_init, 1))
}

# maximum statistic over the selected entries (logic of .null_maxib()); also
# returns the number of excluded entries
max_stat <- function(obs, sims, sel, scale = NULL) {
  sel <- sel & !is.na(obs)
  n_exc <- 0
  if (!is.null(scale)) {
    ok <- !is.na(scale) & scale > 0
    n_exc <- sum(sel & !ok)
    sel <- sel & ok
    sc <- ifelse(ok, scale, 1)
    obs <- obs / sc
    sims <- sims / sc
  }
  if (!any(sel)) {
    return(c(obs = NA_real_, q95 = NA_real_, p = NA_real_,
             exceed = NA_real_, entries = 0, excluded = n_exc))
  }
  m_obs <- max(abs(obs[sel]))
  m_null <- apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  q95 <- stats::quantile(m_null, 0.95, na.rm = TRUE, names = FALSE)
  c(obs = m_obs, q95 = q95,
    p = (1 + sum(m_null >= m_obs, na.rm = TRUE)) / (1 + sum(!is.na(m_null))),
    exceed = as.numeric(m_obs > q95), entries = sum(sel), excluded = n_exc)
}

diff_max <- function(a, b) {
  d <- abs(a - b)
  d[is.na(a) & is.na(b)] <- 0
  d[xor(is.na(a), is.na(b))] <- Inf
  if (length(d) == 0) 0 else max(d)
}

one_rep <- function(nm, rep) {
  t0 <- Sys.time()
  sc <- scenarios[[nm]]$scenario
  calib <- scenarios[[nm]]$calib
  pars <- do.call(ccw_sim_params, sc$params)
  d <- simulate_ccw(pars)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, sc$formula, truncate = sc$truncate)
  null_seed <- sample.int(.Machine$integer.max, 1)
  feats <- c(all_feats, as.character(sc$extra_feats))
  data <- as.data.frame(w)
  n <- nrow(data)
  poss <- as.logical(data$censor_possible) %in% TRUE
  dd <- data$censored_next
  pt <- if (calib) true_pi(data, pars) else NULL
  wt <- attr(w, "weighting")
  by <- if (!is.null(wt$by)) wt$by else c("strategy", "phase")

  # diagnostic entries and package observed RRB (as rrb_null())
  data$.w_post_null <- data$w_ref / data$pi_hat
  ob <- rrb(data, feats, method = "weighted", scope = "interval",
            w_ref = "w_ref", w_post = ".w_post_null")
  ob <- ob[ob$n_censored >= 1 & !is.na(ob$rrb) &
             as.character(ob$strategy) %in% STRATEGIES, , drop = FALSE]
  strat <- as.character(ob$strategy)
  feat <- as.character(ob$feature)
  ckey <- paste(strat, ob$time, sep = "\r")
  preps <- list(
    all = .null_prep(data, feats, "strategy", "time", "w_ref"),
    cens = .null_prep(data[poss, , drop = FALSE], feats, "strategy", "time",
                      "w_ref")
  )
  rows_of <- list(all = seq_len(n), cens = which(poss))
  pos <- lapply(preps, function(p) cbind(match(ckey, p$cells),
                                         match(feat, feats)))
  if (anyNA(pos$all) || anyNA(pos$cens)) stop("Unmatched diagnostic cell.")
  rrb_v <- function(v, dvec, pvec) {
    ix <- rows_of[[v]]
    .null_rrb(preps[[v]], dvec[ix], pvec[ix])[pos[[v]]]
  }
  obs <- list(fit = lapply(stats::setNames(VARIANTS, VARIANTS), rrb_v,
                           dvec = dd, pvec = data$pi_hat))
  if (calib) {
    obs$true <- lapply(stats::setNames(VARIANTS, VARIANTS), rrb_v,
                       dvec = dd, pvec = pt)
  }
  check_pkg <- diff_max(obs$fit$all, ob$rrb)
  if (check_pkg > 1e-10) {
    stop("Observed RRB differs from rrb() (max diff ", format(check_pkg), ").")
  }

  # cell structure (rows of the diagnostic cells)
  rkey <- paste(as.character(data$strategy), data$time, sep = "\r")
  cell_rows <- lapply(stats::setNames(STRATEGIES, STRATEGIES), function(g) {
    which(rkey %in% unique(ckey[strat == g]))
  })
  n_risk_c <- table(factor(rkey, levels = unique(ckey)))
  n_poss_c <- table(factor(rkey[poss], levels = unique(ckey)))
  full_cell <- (n_risk_c == n_poss_c)[ckey]  # no pi = 1 rows

  # null draws (common random numbers across the schemes)
  groups <- .null_groups(data, sc$formula, by, "censor_possible")
  key <- paste(data$id[poss], data$time[poss], sep = "\r")
  grp <- match(key, unique(key))
  n_grp <- max(c(0L, grp))
  on0 <- .null_orient(data, "strategy", "censored_next", "censor_possible",
                      "A")[poss]
  schemes <- if (calib) SCHEMES else "fit_refit"
  p_src <- list(fit = data$pi_hat, true = pt)
  sims <- lapply(stats::setNames(VARIANTS, VARIANTS), function(v) {
    lapply(stats::setNames(schemes, schemes), function(s) {
      matrix(NA_real_, nrow(ob), N_DRAWS)
    })
  })
  set.seed(null_seed)
  for (b in seq_len(N_DRAWS)) {
    u <- stats::runif(n_grp)  # shared by the schemes
    for (s in schemes) {
      p <- p_src[[sub("_.*$", "", s)]]
      d_star <- integer(n)
      d_star[poss] <- .null_draw(p[poss], on0, grp, u)
      pi_star <- if (grepl("_refit$", s)) {
        .null_refit(groups, 1 - d_star, n)
      } else {
        p
      }
      pi_star <- pmax(pi_star, PI_FLOOR)
      for (v in VARIANTS) sims[[v]][[s]][, b] <- rrb_v(v, d_star, pi_star)
    }
  }

  # all = cens in cells without pi = 1 rows (observed and draws)
  check_cens <- vapply(STRATEGIES, function(g) {
    sel <- strat == g & full_cell
    if (!any(sel)) return(NA_real_)
    max(diff_max(obs$fit$all[sel], obs$fit$cens[sel]),
        diff_max(sims$all$fit_refit[sel, ], sims$cens$fit_refit[sel, ]))
  }, numeric(1))
  if (any(check_cens > 1e-10, na.rm = TRUE)) {
    stop("cens and all differ in cells without pi = 1 rows.")
  }

  # optional: S1 all against rrb_null() with the same null seed
  check_null <- NA_real_
  t_check <- 0
  if (CHECK_PKG) {
    t1 <- Sys.time()
    nr <- rrb_null(w, feats, n_draws = N_DRAWS, seed = null_seed,
                   pi_floor = PI_FLOOR)
    s1 <- sims$all$fit_refit
    sd1 <- apply(s1, 1, stats::sd, na.rm = TRUE)
    check_null <- max(vapply(STRATEGIES, function(g) {
      sel <- strat == g
      a <- max_stat(obs$fit$all, s1, sel)
      z <- max_stat(obs$fit$all, s1, sel, sd1)
      r1 <- nr$maxib[nr$maxib$strategy == g & nr$maxib$feature == "(all)", ]
      r2 <- nr$maxib_stud[nr$maxib_stud$strategy == g &
                            nr$maxib_stud$feature == "(all)", ]
      if (nrow(r1) == 0) return(if (any(sel)) Inf else 0)
      diff_max(c(a[c("obs", "q95", "p")], z[c("obs", "q95", "p")],
                 z[["excluded"]]),
               c(r1$maxib_obs, r1$null_q95, r1$p, r2$stud_obs, r2$null_q95,
                 r2$p, r2$n_excluded))
    }, numeric(1)))
    if (check_null > 1e-10) {
      stop("S1 (all) differs from rrb_null() (max diff ",
           format(check_null), ").")
    }
    t_check <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  }

  # statistics
  pr <- PAIRS[PAIRS$scheme %in% schemes & PAIRS$observed %in% names(obs), ]
  out <- list()
  for (v in VARIANTS) {
    for (k in seq_len(nrow(pr))) {
      o <- obs[[pr$observed[k]]][[v]]
      sm <- sims[[v]][[pr$scheme[k]]]
      sdv <- apply(sm, 1, stats::sd, na.rm = TRUE)
      for (g in STRATEGIES) {
        sel <- strat == g
        raw <- max_stat(o, sm, sel)
        stud <- max_stat(o, sm, sel, sdv)
        r <- data.frame(
          scenario = nm, rep = rep, n = pars$n, null_seed = null_seed,
          ref = v, pair = pr$pair[k], scheme = pr$scheme[k],
          observed = pr$observed[k], strategy = g,
          n_cells = length(unique(ckey[sel])),
          entries_raw = raw[["entries"]], maxib_obs = raw[["obs"]],
          maxib_q95 = raw[["q95"]], maxib_p = raw[["p"]],
          exceed_raw = raw[["exceed"]], entries_stud = stud[["entries"]],
          n_excluded = stud[["excluded"]], stud_obs = stud[["obs"]],
          stud_q95 = stud[["q95"]], stud_p = stud[["p"]],
          exceed_stud = stud[["exceed"]]
        )
        for (f in feats) {
          r[[paste0("exceed_stud_", f)]] <-
            max_stat(o, sm, sel & feat == f, sdv)[["exceed"]]
        }
        ix <- cell_rows[[g]]
        ret <- ix[poss[ix] & dd[ix] == 0]
        r$n_risk <- length(ix)
        r$n_censorable <- sum(poss[ix])
        r$n_censored <- sum(dd[ix] == 1, na.rm = TRUE)
        r$n_retained_cens <- length(ret)
        r$mean_inv_pihat_ret <- if (length(ret) > 0) {
          mean(1 / data$pi_hat[ret])
        } else {
          NA_real_
        }
        r$check_pkg <- check_pkg
        r$check_cens_eq <- check_cens[[g]]
        r$check_null <- check_null
        out[[length(out) + 1]] <- r
      }
    }
  }
  out <- do.call(rbind, out)
  # runtime without the optional rrb_null() check
  out$runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs")) -
    t_check
  out
}

# tasks and streams ------------------------------------------------------------
tasks <- do.call(rbind, lapply(names(scenarios), function(nm) {
  data.frame(scenario = nm, rep = seq_len(N_REPS[[nm]]),
             stringsAsFactors = FALSE)
}))
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
stream <- .Random.seed
seeds <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  stream <- parallel::nextRNGStream(stream)
  seeds[[i]] <- stream
}
run_idx <- which(tasks$scenario %in% sc_names)  # n = 2000 scenarios first
message(sprintf(
  "Pilot (g2 cause): %d tasks (%s), %d null draws, %d cores%s",
  length(run_idx),
  paste(sprintf("%s x %d", sc_names, N_REPS[sc_names]), collapse = ", "),
  N_DRAWS, N_CORES, if (CHECK_PKG) ", rrb_null() check" else ""
))
t_start <- Sys.time()
res_list <- parallel::mclapply(run_idx, function(i) {
  assign(".Random.seed", seeds[[i]], envir = globalenv())
  tryCatch(one_rep(tasks$scenario[i], tasks$rep[i]), error = function(e) {
    data.frame(scenario = tasks$scenario[i], rep = tasks$rep[i],
               error = conditionMessage(e))
  })
}, mc.cores = N_CORES, mc.preschedule = FALSE)
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
res <- res[order(match(res$scenario, names(scenarios)), res$rep,
                 match(res$ref, VARIANTS), res$pair,
                 match(res$strategy, STRATEGIES)), ]
rownames(res) <- NULL
write.csv(res, out_csv("pilot_g2_cause"), row.names = FALSE)

# summary ----------------------------------------------------------------------
rate_se <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(c(NA_real_, NA_real_))
  c(mean(v), sqrt(mean(v) * (1 - mean(v)) / length(v)))
}
q95 <- function(v) stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE)
fcols <- grep("^exceed_stud_", names(res), value = TRUE)
scols <- c("n_risk", "n_censorable", "n_censored", "n_retained_cens",
           "mean_inv_pihat_ret")
summ <- do.call(rbind, lapply(
  split(res, list(res$scenario, res$ref, res$pair, res$strategy), drop = TRUE),
  function(x) {
    er <- rate_se(x$exceed_raw)
    es <- rate_se(x$exceed_stud)
    out <- data.frame(
      scenario = x$scenario[1], ref = x$ref[1], pair = x$pair[1],
      scheme = x$scheme[1], observed = x$observed[1],
      strategy = x$strategy[1], n_rep = nrow(x),
      n_failed = sum(fail_df$scenario == x$scenario[1]),
      exceed_raw = er[1], exceed_raw_se = er[2],
      exceed_stud = es[1], exceed_stud_se = es[2],
      maxib_obs_mean = mean(x$maxib_obs, na.rm = TRUE),
      mc_q95_maxib_obs = q95(x$maxib_obs),
      maxib_q95_mean = mean(x$maxib_q95, na.rm = TRUE),
      stud_obs_mean = mean(x$stud_obs, na.rm = TRUE),
      mc_q95_stud_obs = q95(x$stud_obs),
      stud_q95_mean = mean(x$stud_q95, na.rm = TRUE),
      n_excluded_mean = mean(x$n_excluded)
    )
    for (cl in fcols) out[[cl]] <- mean(x[[cl]], na.rm = TRUE)
    for (cl in scols) out[[paste0(cl, "_mean")]] <- mean(x[[cl]], na.rm = TRUE)
    out
  }))
summ <- summ[order(match(summ$scenario, names(scenarios)),
                   match(summ$ref, VARIANTS), summ$pair,
                   match(summ$strategy, STRATEGIES)), ]
rownames(summ) <- NULL
write.csv(summ, out_csv("pilot_g2_cause_summary"), row.names = FALSE)

# paired comparisons (exact McNemar) ------------------------------------------
mcnemar <- function(xa, xb, lab) {
  k <- paste(xa$scenario, xa$rep, xa$strategy)
  m <- match(k, paste(xb$scenario, xb$rep, xb$strategy))
  ea <- xa$exceed_stud
  eb <- xb$exceed_stud[m]
  ok <- !is.na(ea) & !is.na(eb)
  ea <- ea[ok]
  eb <- eb[ok]
  n10 <- sum(ea == 1 & eb == 0)
  n01 <- sum(ea == 0 & eb == 1)
  data.frame(lab, n_pairs = sum(ok), n_both = sum(ea == 1 & eb == 1),
             n_a_only = n10, n_b_only = n01, rate_a = mean(ea),
             rate_b = mean(eb),
             p_mcnemar = if (n10 + n01 > 0) {
               stats::binom.test(n10, n10 + n01, 0.5)$p.value
             } else {
               1
             })
}
grp <- unique(res[, c("scenario", "strategy")])
mcn <- do.call(rbind, lapply(seq_len(nrow(grp)), function(j) {
  x <- res[res$scenario == grp$scenario[j] &
             res$strategy == grp$strategy[j], , drop = FALSE]
  sub_x <- function(v, p) x[x$ref == v & x$pair == p, , drop = FALSE]
  rows <- list()
  for (v in VARIANTS) {
    for (p in setdiff(unique(x$pair), "S1")) {
      rows[[length(rows) + 1]] <- mcnemar(
        sub_x(v, "S1"), sub_x(v, p),
        data.frame(scenario = grp$scenario[j], strategy = grp$strategy[j],
                   ref_a = v, pair_a = "S1", ref_b = v, pair_b = p))
    }
  }
  rows[[length(rows) + 1]] <- mcnemar(
    sub_x("cens", "S1"), sub_x("all", "S1"),
    data.frame(scenario = grp$scenario[j], strategy = grp$strategy[j],
               ref_a = "cens", pair_a = "S1", ref_b = "all", pair_b = "S1"))
  do.call(rbind, rows)
}))
mcn <- mcn[order(match(mcn$scenario, names(scenarios)),
                 match(mcn$strategy, STRATEGIES)), ]
rownames(mcn) <- NULL
write.csv(mcn, out_csv("pilot_g2_cause_mcnemar"), row.names = FALSE)

# printout ---------------------------------------------------------------------
cat("\nCalibration scenarios (exceedance of the null q95):\n")
cal <- summ[summ$scenario %in% names(scenarios)[vapply(
  scenarios, function(s) s$calib, logical(1))], ]
print(cal[, c("scenario", "ref", "pair", "strategy", "n_rep", "exceed_raw",
              "exceed_stud", "exceed_stud_se", "stud_obs_mean",
              "mc_q95_stud_obs", "stud_q95_mean", "n_excluded_mean")],
      digits = 3, row.names = FALSE)
cat("\nPower scenarios (S1):\n")
print(summ[!summ$scenario %in% cal$scenario,
           c("scenario", "ref", "pair", "strategy", "n_rep", "exceed_raw",
             "exceed_stud")], digits = 3, row.names = FALSE)
cat("\ng2 cell structure (means):\n")
print(unique(summ[summ$strategy == "g2" & summ$pair == "S1" &
                    summ$ref == "all",
                  c("scenario", paste0(scols, "_mean"))]),
      digits = 4, row.names = FALSE)
cat(sprintf(paste0("\nChecks: max |obs_fit - rrb()| %.2g; max |all - cens| ",
                   "in cells without pi = 1 rows %.2g; rrb_null() check %s\n"),
            max(res$check_pkg), max(c(0, res$check_cens_eq), na.rm = TRUE),
            if (CHECK_PKG) sprintf("%.2g", max(res$check_null)) else "off"))
first <- !duplicated(res[, c("scenario", "rep")])
cat(sprintf(paste0("Runtime: %.2f min wall; mean per-task %.1f s; ",
                   "failed: %d of %d\n"),
            wall, mean(res$runtime_sec[first]), nrow(fail_df),
            length(run_idx)))
print(tapply(res$runtime_sec[first], res$scenario[first], mean))
cat("Written:", out_csv("pilot_g2_cause"), out_csv("pilot_g2_cause_summary"),
    out_csv("pilot_g2_cause_mcnemar"), sep = "\n  ")
