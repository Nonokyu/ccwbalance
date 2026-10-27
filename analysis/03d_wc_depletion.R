# 03d_wc_depletion.R --------------------------------------------------------
# Depletion term of the end-of-grace SMD in the Webster-Clark SynPUF example.
#
# PURPOSE: the end-of-grace between-strategy SMD of a baseline covariate
#   compares the weighted individuals retained (and event-free) under two
#   strategies; after outcome events it contains a depletion term
#   (depletion equation of Section 3.2 of the manuscript), whereas
#   cumulative RRB carries individuals with earlier events forward (weights
#   rescaled by the stabilising numerators) and targets the baseline
#   population. For the baseline covariates at day 30 and day 90 (interval
#   starting on that day, `w30`/`w90` of 03_empirical_webster_clark.R) this
#   script reports
#     * the weighted and unweighted between-strategy SMD (strategy_balance()),
#     * the cumulative RRB of each strategy (rrb(scope = "cumulative")),
#     * the within-strategy depletion component
#         dep_g = (mu_B,g - mu_cum,g) / s0
#       (survivors' mean minus the cumulative-reconstruction mean, on the
#       baseline-SD scale of cumulative RRB), so that mu_B,g - mu0 =
#       s0 * (cum_rrb_g + dep_g), and the exact decomposition of the SMD
#         SMD = s0 / s_ref * {(dep_g - dep_g') + (cum_rrb_g - cum_rrb_g')},
#     * the cohort-level depletion (unweighted mean among individuals
#       observed event-free at that day minus the baseline mean, / s0),
#       with the crude incidence F and the mean among events mu_E.
#   Data construction, censoring model and weights are identical to
#   analysis/03_empirical_webster_clark.R (DELTA = 10, age + sex + renal).
#
# INPUTS:  data-raw/webster_clark/Data/baseccwcohort.sas7bdat (git-ignored;
#          git clone https://github.com/mawcpharmdphd/CCW_tutorial
#          data-raw/webster_clark)
# OUTPUTS: output/empirical/wc_depletion.csv (tidy: time, time_day,
#          covariate, strategy, pair, weighted, quantity, value)
#
# Run from the project root:  Rscript analysis/03d_wc_depletion.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
data_file <- file.path("data-raw", "webster_clark", "Data",
                       "baseccwcohort.sas7bdat")
if (!file.exists(data_file)) {
  stop("Data not found: ", data_file, "\nRun:\n  git clone ",
       "https://github.com/mawcpharmdphd/CCW_tutorial data-raw/webster_clark")
}
t_start <- Sys.time()
suppressMessages(pkgload::load_all(".", quiet = TRUE))

DELTA <- 10
HORIZON <- 180
K <- HORIZON / DELTA
out_dir <- file.path("output", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Cohort, long data and weights (as in 03_empirical_webster_clark.R) -----
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

# cross-check against the saved outputs of 03_empirical_webster_clark.R
est <- estimate_risk(w)
saved_est <- utils::read.csv(file.path(out_dir, "wc_estimates.csv"))
r_cum <- rrb(w, covs, method = c("unweighted", "weighted"),
             scope = "cumulative")
saved_rrb <- utils::read.csv(file.path(out_dir, "wc_rrb.csv"))
saved_cum <- saved_rrb[saved_rrb$scope == "cumulative", ]
chk <- merge(r_cum, saved_cum,
             by = c("strategy", "time", "feature", "method", "scope"))
cat(sprintf("Cross-check: max |risk - saved| = %.2e; ",
            max(abs(est$risk[est$time == K - 1] -
                      saved_est$estimate[1:4]))),
    sprintf("max |cumulative RRB - saved| = %.2e (%d rows)\n",
            max(abs(chk$rrb.x - chk$rrb.y), na.rm = TRUE), nrow(chk)))

# 2. Quantities at day 30 and day 90 ----------------------------------------
times <- c(w30, w90)
strat <- names(strategies)
pairs <- list(c("0_90", "0_30"), c("0_30", "no_0_90"),
              c("0_90", "no_0_90"), c("30_90", "0_30"),
              c("30_90", "no_0_90"), c("0_90", "30_90"))
row_out <- function(tt, cv, g, pr, wt, q, v) {
  data.frame(time = tt, time_day = tt * DELTA, covariate = cv, strategy = g,
             pair = pr, weighted = wt, quantity = q, value = v)
}
ret <- w[w$censored_next == 0, ]
res <- list()
for (wt in c(TRUE, FALSE)) {
  mth <- if (wt) "weighted" else "unweighted"
  rc <- r_cum[r_cum$method == mth, ]
  # package strategy balance (unweighted: all at-risk rows), all references
  sb <- do.call(rbind, lapply(strat, function(r) {
    strategy_balance(w, covs, reference = r, weighted = wt)
  }))
  for (tt in times) {
    for (cv in covs) {
      # within-strategy survivor means (retained rows) and components
      comp <- lapply(strat, function(g) {
        x <- ret[ret$strategy == g & ret$time == tt, ]
        wi <- if (wt) x$w_post else rep(1, nrow(x))
        cr <- rc[rc$strategy == g & rc$time == tt & rc$feature == cv, ]
        mu_b <- .wmean(x[[cv]], wi)
        dep <- (mu_b - cr$mean_post) / cr$sd_pre
        # weighted share of earlier events in the cumulative set (F_g)
        ev <- w[w$strategy == g & w$censored_next == 0 & w$time < tt &
                  w$event_next == 1, ]
        list(g = g, mu_b = mu_b, var_b = .wvar(x[[cv]], wi), dep = dep,
             cum = cr$rrb, mu0 = cr$mean_pre, s0 = cr$sd_pre,
             mu_cum = cr$mean_post, n_ret = nrow(x), n_ev = nrow(ev),
             ess = sum(wi)^2 / sum(wi^2))
      })
      names(comp) <- strat
      for (g in strat) {
        cg <- comp[[g]]
        res[[length(res) + 1]] <- rbind(
          row_out(tt, cv, g, NA, wt, "mean_survivors", cg$mu_b),
          row_out(tt, cv, g, NA, wt, "mean_cumulative", cg$mu_cum),
          row_out(tt, cv, g, NA, wt, "mean_baseline", cg$mu0),
          row_out(tt, cv, g, NA, wt, "sd_baseline", cg$s0),
          row_out(tt, cv, g, NA, wt, "cum_rrb", cg$cum),
          row_out(tt, cv, g, NA, wt, "depletion", cg$dep),
          row_out(tt, cv, g, NA, wt, "survivor_vs_baseline",
                  cg$cum + cg$dep),
          row_out(tt, cv, g, NA, wt, "n_retained", cg$n_ret),
          row_out(tt, cv, g, NA, wt, "n_earlier_events", cg$n_ev),
          row_out(tt, cv, g, NA, wt, "ess_retained", cg$ess)
        )
      }
      for (p in pairs) {
        a <- comp[[p[1]]]
        b <- comp[[p[2]]]
        s_ref <- sqrt((a$var_b + b$var_b) / 2)
        k <- a$s0 / s_ref
        pr <- paste(p, collapse = " vs ")
        sbp <- sb[sb$strategy == p[1] & sb$reference == p[2] &
                    sb$time == tt & sb$feature == cv, "smd"]
        res[[length(res) + 1]] <- rbind(
          row_out(tt, cv, NA, pr, wt, "smd_package", sbp),
          row_out(tt, cv, NA, pr, wt, "smd_retained", (a$mu_b - b$mu_b) /
                    s_ref),
          row_out(tt, cv, NA, pr, wt, "smd_depletion_part",
                  k * (a$dep - b$dep)),
          row_out(tt, cv, NA, pr, wt, "smd_cumrrb_part",
                  k * (a$cum - b$cum))
        )
      }
    }
  }
}

# 3. Cohort-level (strategy-free, unweighted) depletion ---------------------
for (tt in times) {
  free <- cohort$event_time >= tt & cohort$ltfu_time >= tt
  evd <- cohort$event_time < tt
  for (cv in covs) {
    x <- cohort[[cv]]
    mu0 <- mean(x)
    s0 <- sqrt(mean((x - mu0)^2))
    res[[length(res) + 1]] <- rbind(
      row_out(tt, cv, "cohort", NA, FALSE, "cohort_F_crude", mean(evd)),
      row_out(tt, cv, "cohort", NA, FALSE, "cohort_mean_events",
              mean(x[evd])),
      row_out(tt, cv, "cohort", NA, FALSE, "cohort_n_eventfree", sum(free)),
      row_out(tt, cv, "cohort", NA, FALSE, "cohort_depletion",
              (mean(x[free]) - mu0) / s0),
      row_out(tt, cv, "cohort", NA, FALSE, "cohort_events_vs_baseline",
              (mean(x[evd]) - mu0) / s0)
    )
  }
}
out <- do.call(rbind, res)
rownames(out) <- NULL
write.csv(out, file.path(out_dir, "wc_depletion.csv"), row.names = FALSE)

# 4. Console summary ---------------------------------------------------------
show <- function(q, wt, by = "pair") {
  x <- out[out$quantity == q & out$weighted == wt, ]
  x$key <- paste0(x$time_day, "d ", x[[by]])
  m <- stats::xtabs(value ~ key + covariate, data = x)
  print(round(unclass(m)[, covs, drop = FALSE], 3))
}
for (wt in c(TRUE, FALSE)) {
  cat(sprintf("\n==== %s ====\n", if (wt) "WEIGHTED" else "UNWEIGHTED"))
  for (q in c("smd_package", "smd_retained", "smd_depletion_part",
              "smd_cumrrb_part")) {
    cat("\n--", q, "\n")
    show(q, wt)
  }
  for (q in c("cum_rrb", "depletion", "survivor_vs_baseline",
              "ess_retained")) {
    cat("\n--", q, "\n")
    show(q, wt, "strategy")
  }
}
cat("\n==== COHORT (unweighted) ====\n")
for (q in c("cohort_F_crude", "cohort_depletion",
            "cohort_events_vs_baseline")) {
  cat("\n--", q, "\n")
  show(q, FALSE, "strategy")
}
cat(sprintf("\nRun time: %.1f s\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
