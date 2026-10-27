# 02g_null_coupling_check.R ---------------------------------------------------
# Check of the within-individual coupling of the conditional parametric null
# reference (rrb_null()). The default redraw gives the censorable clones of
# one individual at one time (id x time, a "tied group") one shared uniform
# variate U, and clone g is "treated in t" iff U < p_g, where
# p_g = 1 - pi_hat (clone censored on A_t = 1) or pi_hat (censored on
# A_t = 0) is the initiation probability implied by its strategy-specific
# censoring model; it is censored iff treated (censored on A_t = 1) or iff
# untreated (censored on A_t = 0). The redrawn states of two clones are mutually
# inconsistent (one treated, the other untreated) with probability
# |p_g - p_g'|, i.e. exactly when the strategy-specific models disagree.
#
# For each replicate (or the cohort) the script reports
#   * tied groups: number, fraction of censorable rows in tied groups, and
#     the within-group range of p (max p - min p): mean, median, 95th pct,
#     separately for groups whose clones share the orientation (g2 / g3 at
#     t < w1) and groups with opposite orientations (g1 grace_end vs g2 / g3
#     no_init at t = w1 - 1; g2 grace_end vs g3 at t = w2 - 1);
#   * expected inconsistency rate (mean range over tied groups) and the
#     empirical rate under the shared-U and the independent redraw;
#   * the MaxIB null reference (95th percentile, exceedance of the observed
#     MaxIB) under three redraw schemes:
#       independent  one Bernoulli(1 - pi_hat) per censorable row
#                    (rrb_null(coupling = "independent")),
#       shared       one U per id x time (rrb_null() default, .null_draw()),
#       common       one A* ~ Bernoulli(p_common) per id x time from a single
#                    initiation model P(A_t = 1 | H_t) fitted to the
#                    untreated person-time; D* = A* (censored on A_t = 1) or
#                    1 - A* (censored on A_t = 0);
#     each followed by the package's refit (.null_refit()) and interval RRB
#     (.null_rrb()); same cell selection as rrb_null();
#   * the null SD of the interval RRB of X in the g1 grace-end cell and the
#     g3 cell at t = w1 - 1 under each scheme.
#
# Part 1: simulation, correctly specified censoring model, n = 2000 and 300,
#         MaxIB over g1, g3 and eight history functions (as in the paper).
# Part 2: Webster-Clark SynPUF cohort (as analysis/03b, section 1), MaxIB
#         over all strategies and per strategy, five covariates. Skipped
#         with a message if the data are absent.
#
# Outputs: output/sim/null_coupling_check.csv     one row per replicate
#          output/sim/null_coupling_summary.csv   one row per n
#          output/empirical/wc_null_coupling.csv  scheme x strategy
#
# Configuration: CCW_CORES (default detectCores() - 1), N_REPS (200),
#                N_DRAWS (300; the Webster-Clark part always uses 500),
#                COUPLING_SIM (default 1; 0 skips part 1 and leaves the
#                output/sim files unchanged).
# Run from the project root:  Rscript analysis/02g_null_coupling_check.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))

env_num <- function(x, default) as.numeric(Sys.getenv(x, as.character(default)))
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
N_REPS <- env_num("N_REPS", 200)
N_DRAWS <- env_num("N_DRAWS", 300)
N_DRAWS_WC <- 500
RUN_SIM <- Sys.getenv("COUPLING_SIM", "1") == "1"
SCHEMES <- c("independent", "shared", "common")
FEATS <- c("X", "L", "L_sq", "L_lagprod", "L_diff", "L_cummean", "L_gt0",
           "L_x_X")
SIM_FORMULA <- ~ L + I(L^2) + X + L:X
COMMON_SIM <- A ~ L + I(L^2) + X + L:X + factor(time)
sim_dir <- file.path("output", "sim")
emp_dir <- file.path("output", "empirical")
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

# ---------------------------------------------------------------------------
# helpers

# Fixed pieces of the null reference, exactly as in rrb_null()
null_setup <- function(w, features) {
  wt <- attr(w, "weighting")
  by <- if (!is.null(wt$by)) wt$by else c("strategy", "phase")
  data <- as.data.frame(w)
  data$.w_post_null <- data$w_ref / data$pi_hat
  obs <- rrb(data, features, method = "weighted", scope = "interval",
             w_ref = "w_ref", w_post = ".w_post_null")
  obs <- obs[obs$n_censored >= 1 & !is.na(obs$rrb), , drop = FALSE]
  prep <- .null_prep(data, features, "strategy", "time", "w_ref")
  groups <- .null_groups(data, wt$formula, by, "censor_possible")
  poss <- as.logical(data$censor_possible)
  key <- paste(data$id[poss], data$time[poss], sep = "\r")
  grp <- match(key, unique(key))
  on0 <- .null_orient(data, "strategy", "censored_next", "censor_possible",
                      "A")[poss]
  pi <- data$pi_hat[poss]
  pos <- cbind(match(paste(as.character(obs$strategy), obs$time, sep = "\r"),
                     prep$cells),
               match(as.character(obs$feature), features))
  list(data = data, obs = obs, prep = prep, groups = groups, poss = poss,
       grp = grp, n_grp = max(c(0L, grp)), on0 = on0, pi = pi,
       p_init = ifelse(on0, pi, 1 - pi), pos = pos)
}

# Tied-group statistics of the implied initiation probabilities
tie_stats <- function(s) {
  sz <- tabulate(s$grp, s$n_grp)
  rng <- as.vector(tapply(s$p_init, s$grp, max) -
                     tapply(s$p_init, s$grp, min))
  n0 <- as.vector(rowsum(as.numeric(s$on0), s$grp))
  tied <- sz >= 2
  opp <- tied & n0 > 0 & n0 < sz
  same <- tied & !opp
  q <- function(v, pr) {
    if (length(v) == 0) NA_real_ else
      stats::quantile(v, pr, names = FALSE)
  }
  list(stats = data.frame(
    n_censorable = length(s$grp), n_groups = s$n_grp,
    n_tied = sum(tied), n_tied_same = sum(same), n_tied_opp = sum(opp),
    frac_rows_tied = sum(sz[tied]) / length(s$grp),
    range_mean = mean(rng[tied]), range_median = q(rng[tied], 0.5),
    range_q95 = q(rng[tied], 0.95),
    range_mean_same = if (any(same)) mean(rng[same]) else NA_real_,
    range_mean_opp = if (any(opp)) mean(rng[opp]) else NA_real_,
    incons_expected = mean(rng[tied])
  ), tied = tied, sz = sz)
}

# Common initiation model: probability of initiation per id x time group
common_p <- function(s, fit_rows, formula) {
  fit <- stats::glm(formula, family = stats::binomial(), data = fit_rows)
  first <- which(s$poss)[!duplicated(s$grp)]
  nd <- s$data[first, , drop = FALSE]
  # groups are numbered in order of first appearance, so `first` is sorted
  # by group
  unname(stats::predict(fit, newdata = nd, type = "response"))
}

# One redraw of D* on the censorable rows under a scheme
draw_d <- function(scheme, s, p_common) {
  switch(scheme,
    independent = stats::rbinom(length(s$pi), 1, 1 - s$pi),
    shared = .null_draw(s$pi, s$on0, s$grp, stats::runif(s$n_grp)),
    common = {
      a <- stats::rbinom(s$n_grp, 1, p_common)[s$grp]
      as.integer(ifelse(s$on0, 1 - a, a))
    })
}

# Null draws of the interval RRB (cells x draws) under one scheme, plus the
# empirical inconsistency rate over tied groups
run_scheme <- function(scheme, s, ts, n_draws, seed, p_common = NULL,
                       pi_floor = 1e-6) {
  set.seed(seed)
  n <- nrow(s$data)
  sims <- matrix(NA_real_, nrow(s$obs), n_draws)
  incons <- numeric(n_draws)
  for (b in seq_len(n_draws)) {
    d_poss <- draw_d(scheme, s, p_common)
    treated <- ifelse(s$on0, 1 - d_poss, d_poss)
    k <- as.vector(rowsum(treated, s$grp))
    incons[b] <- mean((k > 0 & k < ts$sz)[ts$tied])
    d_star <- integer(n)
    d_star[s$poss] <- d_poss
    pi_star <- pmax(.null_refit(s$groups, 1 - d_star, n), pi_floor)
    sims[, b] <- .null_rrb(s$prep, d_star, pi_star)[s$pos]
  }
  list(sims = sims, incons = mean(incons))
}

# MaxIB observed vs null over the selected cells
maxib_ref <- function(s, sims, sel) {
  obs <- max(abs(s$obs$rrb[sel]))
  m <- apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
  q95 <- stats::quantile(m, 0.95, na.rm = TRUE, names = FALSE)
  data.frame(n_cells = sum(sel), maxib_obs = obs,
             null_mean = mean(m, na.rm = TRUE), null_q95 = q95,
             p = (1 + sum(m >= obs, na.rm = TRUE)) / (1 + sum(!is.na(m))),
             exceed = as.numeric(obs > q95))
}

cell_sd <- function(s, sims, g, tt, f = "X") {
  i <- which(as.character(s$obs$strategy) == g & s$obs$time == tt &
               as.character(s$obs$feature) == f)
  if (length(i) != 1) NA_real_ else stats::sd(sims[i, ], na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# 1. simulation
one_rep <- function(n, rep) {
  pars <- ccw_sim_params(n = n)
  w1 <- pars$window[1]
  d <- simulate_ccw(pars, seed = 7000 + rep)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, SIM_FORMULA)
  s <- null_setup(w, FEATS)
  ts <- tie_stats(s)
  fit_rows <- s$data[s$data$strategy == "g3" & s$data$phase == "no_init", ]
  pc <- common_p(s, fit_rows, COMMON_SIM)
  sel <- as.character(s$obs$strategy) %in% c("g1", "g3")
  out <- cbind(data.frame(n = n, rep = rep), ts$stats)
  for (sc in SCHEMES) {
    r <- run_scheme(sc, s, ts, N_DRAWS, 100 + rep, pc)
    m <- maxib_ref(s, r$sims, sel)
    out$maxib_obs <- m$maxib_obs
    out$n_cells <- m$n_cells
    out[[paste0("q95_", sc)]] <- m$null_q95
    out[[paste0("p_", sc)]] <- m$p
    out[[paste0("exceed_", sc)]] <- m$exceed
    out[[paste0("incons_emp_", sc)]] <- r$incons
    out[[paste0("sd_x_g1_", sc)]] <- cell_sd(s, r$sims, "g1", w1 - 1)
    out[[paste0("sd_x_g3_", sc)]] <- cell_sd(s, r$sims, "g3", w1 - 1)
  }
  out
}

if (RUN_SIM) {
  tasks <- expand.grid(rep = seq_len(N_REPS), n = c(2000, 300))
  tasks <- tasks[order(-tasks$n, tasks$rep), ]
  message(sprintf("Simulation: %d tasks, %d draws x %d schemes, %d cores",
                  nrow(tasks), N_DRAWS, length(SCHEMES), N_CORES))
  res_list <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    tryCatch(one_rep(tasks$n[i], tasks$rep[i]), error = function(e) {
      data.frame(n = tasks$n[i], rep = tasks$rep[i],
                 error = conditionMessage(e))
    })
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  errs <- vapply(res_list, function(x) "error" %in% names(x), logical(1))
  if (any(errs)) {
    print(do.call(rbind, res_list[errs]))
    stop(sum(errs), " replicate(s) failed.")
  }
  res <- do.call(rbind, res_list)
  write.csv(res, file.path(sim_dir, "null_coupling_check.csv"),
            row.names = FALSE)

  summ <- do.call(rbind, lapply(c(2000, 300), function(nn) {
    x <- res[res$n == nn, ]
    num <- setdiff(names(x), c("n", "rep"))
    out <- data.frame(n = nn, n_rep = nrow(x))
    for (v in num) out[[paste0(v, "_mean")]] <- mean(x[[v]], na.rm = TRUE)
    for (sc in SCHEMES) {
      e <- x[[paste0("exceed_", sc)]]
      e <- e[!is.na(e)]
      out[[paste0("exceed_rate_", sc)]] <- mean(e)
      out[[paste0("exceed_rate_se_", sc)]] <- sqrt(mean(e) * (1 - mean(e)) /
                                                     length(e))
    }
    q <- x[, paste0("q95_", SCHEMES)]
    out$q95_max_abs_diff_mean <- mean(apply(q, 1, function(v) diff(range(v))))
    out$q95_scheme_mean_max_diff <- diff(range(colMeans(q)))
    out$mc_q95_maxib_obs <- stats::quantile(x$maxib_obs, 0.95, names = FALSE)
    out
  }))
  write.csv(summ, file.path(sim_dir, "null_coupling_summary.csv"),
            row.names = FALSE)
  cat("\n=== Simulation: coupling check (means over replicates) ===\n")
  print(t(summ), digits = 4)
} else {
  message("COUPLING_SIM = 0: simulation part skipped.")
}

# ---------------------------------------------------------------------------
# 2. Webster-Clark cohort (data preparation as analysis/03b, section 1)
data_file <- file.path("data-raw", "webster_clark", "Data",
                       "baseccwcohort.sas7bdat")
if (!file.exists(data_file)) {
  message("Webster-Clark data not found (", data_file, "); part 2 skipped.")
} else {
  dir.create(emp_dir, recursive = TRUE, showWarnings = FALSE)
  DELTA <- 10
  HORIZON <- 180
  K <- HORIZON / DELTA
  raw <- haven::read_sas(data_file)
  cohort <- data.frame(
    id = seq_len(nrow(raw)), ID = raw$ID,
    age = raw$baseline_age,
    female = as.numeric(raw$sex == "2"),
    renal = as.numeric(raw$renal == "Y"),
    race_black = as.numeric(raw$race == "2"),
    race_other = as.numeric(raw$race %in% c("3", "5")),
    rx = raw$RX_start, fu = raw$followup, outcome = raw$outcome
  )
  day_to_int <- function(d) floor(d / DELTA)
  cohort$init_time <- ifelse(is.na(cohort$rx) | cohort$rx >= HORIZON, Inf,
                             day_to_int(cohort$rx))
  cohort$event_time <- ifelse(cohort$outcome == 1 & cohort$fu <= HORIZON,
                              pmin(day_to_int(cohort$fu), K - 1), Inf)
  cohort$ltfu_time <- ifelse(cohort$outcome == 0 & cohort$fu < HORIZON,
                             day_to_int(cohort$fu), Inf)

  covs <- c("age", "female", "renal", "race_black", "race_other")
  w30 <- day_to_int(30)
  w90 <- day_to_int(90)
  strategies <- list("0_30" = c(0, w30 - 1), "0_90" = c(0, w90 - 1),
                     "30_90" = c(w30, w90 - 1), "no_0_90" = c(NA, w90 - 1))
  long <- ccw_clone(cohort, strategies, n_time = K, covariates = covs)
  formulas <- list(
    no_init = ~ age + female + renal + time + I(time^2),
    grace_end = ~ age + female + renal
  )
  w <- fit_censoring_model(long, formulas, stabilized = TRUE)

  s <- null_setup(w, covs)
  ts <- tie_stats(s)
  fit_rows <- s$data[s$data$strategy == "no_0_90" &
                       s$data$phase == "no_init", ]
  pc <- common_p(s, fit_rows, A ~ age + female + renal + factor(time))
  strat <- as.character(s$obs$strategy)
  levs <- c("(all)", intersect(levels(s$obs$strategy), strat))
  wc <- do.call(rbind, lapply(SCHEMES, function(sc) {
    r <- run_scheme(sc, s, ts, N_DRAWS_WC, 20260925, pc)
    do.call(rbind, lapply(levs, function(g) {
      sel <- if (g == "(all)") rep(TRUE, length(strat)) else strat == g
      cbind(data.frame(scheme = sc, strategy = g), maxib_ref(s, r$sims, sel),
            incons_emp = r$incons,
            sd_age_0_30_t2 = cell_sd(s, r$sims, "0_30", w30 - 1, "age"),
            sd_age_no_0_90_t2 = cell_sd(s, r$sims, "no_0_90", w30 - 1, "age"))
    }))
  }))
  wc <- cbind(wc, ts$stats[rep(1, nrow(wc)), ], row.names = NULL)
  write.csv(wc, file.path(emp_dir, "wc_null_coupling.csv"), row.names = FALSE)
  cat("\n=== Webster-Clark: coupling check ===\n")
  print(ts$stats, digits = 4, row.names = FALSE)
  print(wc[, c("scheme", "strategy", "n_cells", "maxib_obs", "null_q95", "p",
               "exceed", "incons_emp", "sd_age_0_30_t2", "sd_age_no_0_90_t2")],
        digits = 4, row.names = FALSE)
}

cat(sprintf("\nTotal runtime: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
