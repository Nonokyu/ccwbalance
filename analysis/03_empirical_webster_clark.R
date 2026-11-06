# 03_empirical_webster_clark.R -----------------------------------------------
# Empirical example: Webster-Clark et al. (2025, Pharmacoepidemiol Drug Saf,
# doi:10.1002/pds.70223) CCW tutorial on a cohort built from the CMS
# Synthetic Medicare Public Use Files (SynPUF): treatment initiation after
# myocardial infarction (MI) discharge, outcome within 180 days.
#
# DATA (public, no registration; ~200 KB):
#   git clone https://github.com/mawcpharmdphd/CCW_tutorial data-raw/webster_clark
#   -> data-raw/webster_clark/Data/baseccwcohort.sas7bdat   (git-ignored)
# The repository also contains the authors' SAS and R code (R_code/), which
# this script re-implements in the ccwbalance long-data contract.
# Time-varying covariates come from the cohort's CMS DE-SynPUF claims:
#   bash data-raw/get_synpuf_claims.sh -> data-raw/synpuf_claims/ (git-ignored)
# Data preparation is shared with 03b, 03c, 03d and 02g (analysis/wc_data.R):
#   visit (carrier claim), outpt (hospital outpatient claim) and fill (Part D
#   dispensing other than the index antiplatelet) in the 10 days before each
#   interval. Three history-depth functions of the same dispensings
#   (fill_lag: fill in the previous window; fill_30: number in the 30 days
#   before the interval; fill_cum: number since discharge) are monitored by
#   the diagnostics but enter no censoring model.
#
# Variables used (one row per patient, n = 1590):
#   ID, Discharge_date, birthdate, baseline_age, sex (SynPUF: 1 = male,
#   2 = female), renal ("Y" = end-stage renal disease), race (1 white,
#   2 black, 3 other, 5 Hispanic), RX_start (days from discharge to first
#   treatment; NA = never), followup (days), outcome (1 = event at `followup`).
#
# Design (as in the tutorial, strategies over a 180-day horizon):
#   "0_30"  initiate within days 0-30      "0_90"  initiate within days 0-90
#   "no_0_90" no initiation during days 0-89 (intervals 0-8); initiation
#             from day 90 onwards is allowed and not censored. (Earlier
#             versions used "none" = never initiate over the 180 days.)
#
# DEVIATIONS from the tutorial (discrete-time re-implementation):
#   * Time is discretised into DELTA-day intervals [k*DELTA, (k+1)*DELTA).
#     Hence "0_30" = initiation on days 0-29 and "0_90" ends on day 89;
#     patients starting exactly on day 30 or day 90 are handled differently
#     from the tutorial (counts printed below).
#   * The tutorial fits the grace-end censoring model among "recent starters"
#     (last 7 days) and non-starters; here the model at the grace-end interval
#     is fitted among all clones still uninitiated at the start of that
#     interval (starters within the last DELTA days vs non-starters).
#   * An outcome (or loss to follow-up) in the grace-end interval precedes the
#     "not started by day 30/90" deviation (ccw_clone() default), as in the
#     tutorial; the early-initiation deviation of "no_0_90" is applied at the
#     start of the interval of initiation.
#   * Loss to follow-up (followup < 180 without outcome) is treated as
#     non-informative (as in the tutorial) and is not weighted.
#   * Risks: weighted discrete-time Kaplan-Meier (estimate_risk()); the
#     tutorial uses the baseline survival of a weighted null Cox model. A
#     direct port of the tutorial's continuous-time "0_30"/"0_90" analyses is
#     included below as a cross-check.
#
# Diagnostics: interval RRB, M^RR (interval), strategy and censoring balance
# over the eleven history functions (five baseline covariates, three
# time-varying claims indicators, three depth functions); cumulative RRB
# over the baseline covariates. Censoring model as in the tutorial (age, sex,
# renal); 03c adds race and/or the time-varying covariates.
#
# Outputs (output/empirical/, git-ignored) include wc_rrb.csv,
# wc_timevarying.csv (id x time covariates and depth functions) and
# wc_timevarying_assoc.csv (column `model`: "indicators" = initiation and
# outcome models with visit, outpt, fill; "indicators + depth" = initiation
# model with fill_lag, fill_30, fill_cum added).
#
# Run from the project root:  Rscript analysis/03_empirical_webster_clark.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "wc_data.R"))
library(ggplot2)

DELTA <- as.numeric(Sys.getenv("WC_DELTA", "10"))
HORIZON <- 180
N_BOOT <- as.numeric(Sys.getenv("WC_N_BOOT", "200"))
fig_dir <- file.path("output", "figures")
out_dir <- file.path("output", "empirical")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_fig <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height,
         dpi = 150)
}

# 1. Cohort, long data and time-varying covariates (analysis/wc_data.R) -------
wc <- wc_prepare(delta = DELTA, horizon = HORIZON)
raw <- wc$raw
cohort <- wc$cohort
long <- wc$long
covs <- wc$covs
tvf <- wc$tvf
depth <- wc$depth
feats <- wc$feats
strategies <- wc$strategies
K <- wc$K
w30 <- wc$w30
w90 <- wc$w90
cat(sprintf("n = %d; initiators within 180 d: %d; outcomes within 180 d: %d;",
            nrow(cohort), sum(is.finite(cohort$init_time)),
            sum(is.finite(cohort$event_time))),
    sprintf("LTFU before 180 d: %d\n", sum(is.finite(cohort$ltfu_time))))
cat(sprintf("Starters exactly on day 30: %d, on day 90: %d, on day 0: %d\n",
            sum(cohort$rx == 30, na.rm = TRUE),
            sum(cohort$rx == 90, na.rm = TRUE),
            sum(cohort$rx == 0, na.rm = TRUE)))
write.csv(wc$tv, file.path(out_dir, "wc_timevarying.csv"), row.names = FALSE)
cat("\nMeans of the time-varying covariates and depth functions by",
    "interval:\n")
print(round(t(sapply(c(tvf, depth),
                     function(v) tapply(wc$tv[[v]], wc$tv$time, mean))), 2))
cat("Overall means:\n")
print(round(colMeans(wc$tv[, c(tvf, depth)]), 3))

# initiation and outcome associations of the time-varying covariates
# (no_0_90 clones in the no-initiation phase, intervals 0-8)
pp <- long[long$strategy == "no_0_90" & long$phase == "no_init", ]
# model "indicators + depth": initiation model with the depth functions
# added, to check whether they add to the prediction of initiation
assoc_fit <- function(y, terms, label) {
  f <- stats::as.formula(paste(y, "~ age + female + renal + time +",
                               paste(terms, collapse = " + ")))
  m <- stats::glm(f, family = stats::binomial(), data = pp)
  cf <- summary(m)$coef
  data.frame(model = label,
             response = ifelse(y == "censored_next", "initiation",
                               "outcome"),
             term = terms, log_or = cf[terms, 1], se = cf[terms, 2],
             or = exp(cf[terms, 1]), p = cf[terms, 4],
             n_rows = stats::nobs(m), n_y = sum(m$y))
}
assoc <- rbind(
  assoc_fit("censored_next", tvf, "indicators"),
  assoc_fit("event_next", tvf, "indicators"),
  assoc_fit("censored_next", c(tvf, depth), "indicators + depth")
)
rownames(assoc) <- NULL
cat("\nInitiation (censored_next) and outcome (event_next) associations of",
    "the time-varying covariates\n(logistic models on the no_0_90",
    "no-initiation rows, adjusted for age, sex, renal, time; model",
    "'indicators + depth' adds fill_lag, fill_30, fill_cum):\n")
print(assoc, digits = 3, row.names = FALSE)
write.csv(assoc, file.path(out_dir, "wc_timevarying_assoc.csv"),
          row.names = FALSE)

# 2. Censoring model and IPCW --------------------------------------------------
# as in the tutorial: age, sex, renal; time terms for the no-initiation phase
formulas <- wc$formulas$base
w <- fit_censoring_model(long, formulas, stabilized = TRUE)
est <- estimate_risk(w)
risk_h <- est[est$time == K - 1, c("strategy", "risk")]
rd <- function(e) {
  r <- stats::setNames(e$risk[e$time == max(e$time)],
                       as.character(e$strategy[e$time == max(e$time)]))
  c(rd_0_90_vs_0_30 = r[["0_90"]] - r[["0_30"]])
}
rd_hat <- rd(est)

# bootstrap (individuals resampled, censoring model refitted)
set.seed(20260923)
boot <- t(vapply(seq_len(N_BOOT), function(b) {
  s <- sample(cohort$id, replace = TRUE)
  cb <- cohort[s, ]
  cb$id <- seq_len(nrow(cb))
  lb <- ccw_clone(cb, strategies, n_time = K, covariates = covs)
  wb <- suppressWarnings(fit_censoring_model(lb, formulas))
  eb <- estimate_risk(wb)
  c(stats::setNames(eb$risk[eb$time == K - 1],
                    as.character(eb$strategy[eb$time == K - 1])), rd(eb))
}, numeric(4)))
res_tab <- data.frame(
  quantity = c(paste0("risk_", names(strategies)), names(rd_hat)),
  estimate = c(risk_h$risk, rd_hat),
  boot_se = apply(boot, 2, stats::sd)
)
res_tab$lcl <- res_tab$estimate - 1.96 * res_tab$boot_se
res_tab$ucl <- res_tab$estimate + 1.96 * res_tab$boot_se
cat(sprintf("\n180-day risks and risk differences (bootstrap B = %d):\n",
            N_BOOT))
print(res_tab, digits = 3, row.names = FALSE)
write.csv(res_tab, file.path(out_dir, "wc_estimates.csv"), row.names = FALSE)

# unweighted (naive CCW) risks for comparison
est0 <- estimate_risk(w, weight = NULL)
cat("\nUnweighted CCW risks at 180 d:\n")
print(est0[est0$time == K - 1, ], digits = 3, row.names = FALSE)

# 3. Cross-check: direct port of the tutorial's 0_30 and 0_90 analysis --------
tutorial_risk <- function(window, recent_from) {
  r <- raw
  nonstart <- is.na(r$RX_start) | r$RX_start > window
  short <- r$followup < window + 1
  cf <- ifelse(short, r$followup, ifelse(nonstart, window, r$followup))
  co <- ifelse(short, r$outcome, ifelse(nonstart, 0, r$outcome))
  cn <- ifelse(short, 0, ifelse(nonstart, 1, 0))
  age_at <- function(start) {
    as.integer(as.numeric(difftime(r$Discharge_date + start, r$birthdate,
                                   units = "days")) / 365.25)
  }
  row0 <- data.frame(id = seq_len(nrow(r)), start = 0,
                     end = pmin(cf, window),
                     y = ifelse(cf <= window, co, 0), w = 1)
  keep2 <- !(cf <= window & cn == 0)
  r2 <- data.frame(
    id = seq_len(nrow(r)), start = window,
    end = ifelse(cn == 1, window, pmin(cf, HORIZON)),
    y = ifelse(cn == 1, 0, ifelse(cf > HORIZON, 0, co)),
    cens = cn, age = age_at(window), female = r$sex == "1",
    renal = r$renal,
    recent = ifelse(cn == 1, 1, as.numeric(r$RX_start >= recent_from))
  )[keep2, ]
  m <- stats::glm(cens ~ age + female + renal, family = stats::binomial(),
                  data = r2[r2$recent == 1, ])
  p_unc <- 1 - stats::predict(m, newdata = r2, type = "response")
  r2$w <- ifelse(r2$cens == 1, 0, ifelse(r2$recent == 1, 1 / p_unc, 1))
  dd <- rbind(row0, r2[, names(row0)])
  dd <- dd[dd$w > 0 & dd$end > dd$start, ]
  fit <- survival::coxph(survival::Surv(start, end, y) ~ 1, data = dd,
                         weights = w, ties = "efron")
  sf <- survival::survfit(fit)
  1 - utils::tail(sf$surv[sf$time <= HORIZON], 1)
}
xcheck <- data.frame(
  strategy = c("0_30", "0_90"),
  tutorial_port = c(tutorial_risk(30, 23), tutorial_risk(90, 83)),
  this_script = risk_h$risk[match(c("0_30", "0_90"), risk_h$strategy)]
)
cat("\nCross-check vs port of the tutorial code (180-day risk):\n")
print(xcheck, digits = 3, row.names = FALSE)
write.csv(xcheck, file.path(out_dir, "wc_crosscheck.csv"), row.names = FALSE)

# risk curves
est$method <- "IPCW"
est0$method <- "unweighted"
p_risk <- ggplot(rbind(est, est0), aes((time + 1) * DELTA, risk,
                                       colour = strategy,
                                       linetype = method)) +
  geom_step() +
  labs(x = "Days since discharge", y = "Cumulative risk",
       title = "Webster-Clark SynPUF cohort: CCW risk curves") +
  theme_bw()
save_fig(p_risk, "wc_risk_curves.png", 8, 5)

# 4. Conventional diagnostics -------------------------------------------------
# (1) weight distribution
ret <- w[w$censored_next == 0, ]
wq <- do.call(rbind, lapply(split(ret, ret$strategy), function(x) {
  q <- stats::quantile(x$w_post, c(0, 0.5, 0.95, 0.99, 1))
  data.frame(strategy = x$strategy[1], mean = mean(x$w_post), min = q[1],
             median = q[2], p95 = q[3], p99 = q[4], max = q[5])
}))
cat("\nWeight distribution (retained clone-intervals, stabilised):\n")
print(wq, digits = 3, row.names = FALSE)
write.csv(wq, file.path(out_dir, "wc_weight_distribution.csv"),
          row.names = FALSE)
p_w <- ggplot(ret, aes(factor(time), w_post)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~strategy, labeller = label_both) + scale_y_log10() +
  labs(x = "Interval", y = "Cumulative stabilised weight (log scale)",
       title = "Weight distribution over follow-up") +
  theme_bw()
save_fig(p_w, "wc_weight_distribution.png", 11, 6)

# (2) end-of-grace SMD across strategies (weighted retained clones at the
#     first interval after the 0-90 grace period), reference 0_30
sb <- strategy_balance(w, feats, reference = "0_30")
eog <- sb[sb$time == w90, ]
cat("\nEnd-of-grace strategy SMD (interval", w90, ", reference 0_30):\n")
print(eog[, c("strategy", "feature", "smd")], digits = 3, row.names = FALSE)
p_sb <- ggplot(sb, aes(time, smd, colour = feature)) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = 2, colour = "grey50") +
  geom_line() + geom_point(size = 1) +
  facet_wrap(~strategy, labeller = label_both) +
  labs(x = "Interval", y = "Weighted SMD vs 0_30",
       title = "Strategy-balance (conventional; descriptive)") +
  theme_bw()
save_fig(p_sb, "wc_strategy_balance.png", 10, 4)

# (3) Jackson-type censoring balance
cb <- rbind(censoring_balance(w, feats, weighted = FALSE),
            censoring_balance(w, feats, weighted = TRUE))
p_cb <- ggplot(cb, aes(time, abs(smd), colour = feature)) +
  geom_line() + geom_point(size = 1) +
  facet_grid(weighted ~ strategy, labeller = label_both) +
  labs(x = "Interval", y = "|SMD| censored vs retained",
       title = "Jackson-type censoring balance") +
  theme_bw()
save_fig(p_cb, "wc_censoring_balance.png", 11, 5)

# 5. Proposed diagnostics -----------------------------------------------------
# interval RRB over the eleven history functions (baseline covariates, the
# three time-varying claims indicators and the three depth functions);
# cumulative RRB over the baseline covariates only
r_int <- rrb(w, feats, method = c("unweighted", "weighted"))
r_cum <- rrb(w, covs, method = c("unweighted", "weighted"),
             scope = "cumulative")
write.csv(rbind(r_int, r_cum), file.path(out_dir, "wc_rrb.csv"),
          row.names = FALSE)
save_fig(plot_love_dynamic(r_int), "wc_love_dynamic.png", 12, 6)
save_fig(plot_love_dynamic(r_cum), "wc_love_cumulative.png", 12, 6)
save_fig(plot_balance_heatmap(r_int[r_int$method == "weighted", ]),
         "wc_heatmap_rrb.png", 10, 9)
# interval scope over all eleven history functions; cumulative scope over the
# baseline covariates only (0-targeted only for baseline-fixed features)
mh <- rbind(
  mahalanobis_rr(w, feats, method = c("unweighted", "weighted"),
                 scope = "interval", ridge = 0.01),
  mahalanobis_rr(w, covs, method = c("unweighted", "weighted"),
                 scope = "cumulative", ridge = 0.01)
)
save_fig(plot_mahalanobis(mh), "wc_mahalanobis.png", 10, 4)
ws <- weight_support(w)
write.csv(ws, file.path(out_dir, "wc_weight_support.csv"), row.names = FALSE)
save_fig(plot_weight_support(ws), "wc_weight_support.png", 12, 3.5)
db <- plot_dashboard(w, feats, ecdf_feature = "age", ecdf_strategy = "0_30",
                     title = "Webster-Clark SynPUF cohort: CCW diagnostics",
                     ridge = 0.01)
save_fig(db, "wc_dashboard.png", 16, 17)

agg <- aggregate_balance(r_int[r_int$n_censored > 0, ], "rrb")
cat("\nInterval RRB summaries (IAI / MaxIB), weighted vs unweighted:\n")
print(agg[order(agg$method, agg$strategy, -agg$maxib), ], digits = 3,
      row.names = FALSE)
aggc <- aggregate_balance(r_cum[r_cum$n_censored > 0, ], "rrb")
cat("\nCumulative RRB summaries (baseline covariates):\n")
print(aggc[order(aggc$method, aggc$strategy, -aggc$maxib), ], digits = 3,
      row.names = FALSE)
cat("\nWeight support at grace-end intervals:\n")
print(ws[ws$time %in% c(w30 - 1, w90 - 1),
         c("strategy", "time", "n_risk", "n_retained", "ess_ratio",
           "entropy", "pi_min", "pr_pi_small", "w99")],
      digits = 3, row.names = FALSE)
