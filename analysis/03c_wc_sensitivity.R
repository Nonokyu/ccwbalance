# 03c_wc_sensitivity.R --------------------------------------------------------
# Webster-Clark SynPUF example (analysis/03): uncertainty of the RRB profile
# and a sensitivity analysis that adds race to the censoring model.
#
# Data preparation and the base censoring model are identical to
# analysis/03_empirical_webster_clark.R and 03b_wc_null_reference.R.
#
# 1. Base model (age, sex, renal; time terms in the no-initiation phase):
#    individual bootstrap of the weighted RRB (rrb_boot(), B = WC_RRB_BOOT,
#    censoring model refitted per resample), cumulative and interval scope.
# 2. Sensitivity model: race_black + race_other added to both phase formulas.
#    Interval and cumulative RRB (unweighted, weighted), weight support, the
#    conditional parametric null reference (rrb_null(), NULL_DRAWS draws) and
#    the 180-day risks / risk differences vs 0_30 with the individual
#    bootstrap SE of analysis/03 (same code, same seed 20260923).
# 3. Comparison table base vs race model (printed and saved).
#
# Outputs (git-ignored, output/empirical/):
#   wc_rrb_boot_cumulative.csv, wc_rrb_boot_interval.csv   (base model)
#   wc_race_rrb.csv, wc_race_weight_support.csv,
#   wc_race_null_maxib.csv, wc_race_null_cells.csv, wc_race_estimates.csv,
#   wc_race_comparison.csv
#
# Run from the project root:  Rscript analysis/03c_wc_sensitivity.R
#   WC_RRB_BOOT (default 200), WC_N_BOOT (default 200), NULL_DRAWS (500)

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
data_file <- file.path("data-raw", "webster_clark", "Data",
                       "baseccwcohort.sas7bdat")
if (!file.exists(data_file)) {
  stop("Data not found: ", data_file, "\nRun:\n  git clone ",
       "https://github.com/mawcpharmdphd/CCW_tutorial data-raw/webster_clark")
}
suppressMessages(pkgload::load_all(".", quiet = TRUE))

RRB_BOOT <- as.numeric(Sys.getenv("WC_RRB_BOOT", "200"))
N_BOOT <- as.numeric(Sys.getenv("WC_N_BOOT", "200"))
N_DRAWS <- as.numeric(Sys.getenv("NULL_DRAWS", "500"))
DELTA <- 10
HORIZON <- 180
K <- HORIZON / DELTA
out_dir <- file.path("output", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}
t_start <- Sys.time()

# 1. Cohort and base censoring model (as in analysis/03) ----------------------
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

# base model: individual bootstrap of the weighted RRB
for (sc in c("cumulative", "interval")) {
  t0 <- Sys.time()
  bt <- rrb_boot(w, covs, B = RRB_BOOT, scope = sc, seed = 20261018)
  save_csv(bt, sprintf("wc_rrb_boot_%s.csv", sc))
  cat(sprintf("rrb_boot (%s, B = %d): %.1f s\n", sc, RRB_BOOT,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# 2. Sensitivity model with race ----------------------------------------------
formulas_race <- list(
  no_init = ~ age + female + renal + race_black + race_other + time +
    I(time^2),
  grace_end = ~ age + female + renal + race_black + race_other
)
w_race <- fit_censoring_model(long, formulas_race, stabilized = TRUE)
r_race <- rbind(
  rrb(w_race, covs, method = c("unweighted", "weighted")),
  rrb(w_race, covs, method = c("unweighted", "weighted"),
      scope = "cumulative")
)
save_csv(r_race, "wc_race_rrb.csv")
ws_race <- weight_support(w_race)
save_csv(ws_race, "wc_race_weight_support.csv")
t0 <- Sys.time()
nr_race <- rrb_null(w_race, covs, n_draws = N_DRAWS, seed = 20260925)
cat(sprintf("rrb_null (race model, %d draws): %.1f s\n", N_DRAWS,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
save_csv(nr_race$maxib, "wc_race_null_maxib.csv")
save_csv(nr_race$cells, "wc_race_null_cells.csv")

# 180-day risks and risk differences with the bootstrap of analysis/03
rd <- function(e) {
  r <- stats::setNames(e$risk[e$time == max(e$time)],
                       as.character(e$strategy[e$time == max(e$time)]))
  c(rd_0_90_vs_0_30 = r[["0_90"]] - r[["0_30"]],
    rd_30_90_vs_0_30 = r[["30_90"]] - r[["0_30"]])
}
estimates <- function(wf, fml) {
  est <- estimate_risk(wf)
  risk_h <- est[est$time == K - 1, c("strategy", "risk")]
  rd_hat <- rd(est)
  set.seed(20260923)
  boot <- t(vapply(seq_len(N_BOOT), function(b) {
    s <- sample(cohort$id, replace = TRUE)
    cb <- cohort[s, ]
    cb$id <- seq_len(nrow(cb))
    lb <- ccw_clone(cb, strategies, n_time = K, covariates = covs)
    wb <- suppressWarnings(fit_censoring_model(lb, fml))
    eb <- estimate_risk(wb)
    c(stats::setNames(eb$risk[eb$time == K - 1],
                      as.character(eb$strategy[eb$time == K - 1])), rd(eb))
  }, numeric(6)))
  tab <- data.frame(
    quantity = c(paste0("risk_", names(strategies)), names(rd_hat)),
    estimate = c(risk_h$risk, rd_hat),
    boot_se = apply(boot, 2, stats::sd)
  )
  tab$lcl <- tab$estimate - 1.96 * tab$boot_se
  tab$ucl <- tab$estimate + 1.96 * tab$boot_se
  tab
}
est_race <- estimates(w_race, formulas_race)
save_csv(est_race, "wc_race_estimates.csv")
cat(sprintf("\nRace model: 180-day risks and risk differences (B = %d):\n",
            N_BOOT))
print(est_race, digits = 3, row.names = FALSE)
# base model, identical to analysis/03's wc_estimates.csv (same seed)
est_base <- estimates(w, formulas)

# 3. Comparison base vs race model --------------------------------------------
r_base <- rbind(
  rrb(w, covs, method = "weighted"),
  rrb(w, covs, method = "weighted", scope = "cumulative")
)
pick_rrb <- function(r, sc, g, tt, f = "race_other") {
  v <- r$rrb[r$method == "weighted" & r$scope == sc &
               as.character(r$strategy) == g & r$time == tt & r$feature == f]
  if (length(v) == 1) v else NA_real_
}
ess_min <- function(ws, tt) min(ws$ess_ratio[ws$time == tt], na.rm = TRUE)
models <- list(base = list(r = r_base, ws = weight_support(w), est = est_base),
               race = list(r = r_race, ws = ws_race, est = est_race))
rows <- list()
for (g in c("0_90", "30_90")) {
  rows[[length(rows) + 1]] <- data.frame(
    quantity = sprintf("interval RRB race_other, %s, time 8 (day 90)", g),
    base = pick_rrb(models$base$r, "interval", g, w90 - 1),
    race = pick_rrb(models$race$r, "interval", g, w90 - 1))
}
for (g in names(strategies)) {
  rows[[length(rows) + 1]] <- data.frame(
    quantity = sprintf("cumulative RRB race_other, %s, time %d (horizon)", g,
                       K - 1),
    base = pick_rrb(models$base$r, "cumulative", g, K - 1),
    race = pick_rrb(models$race$r, "cumulative", g, K - 1))
}
rows[[length(rows) + 1]] <- data.frame(
  quantity = "smallest ESS/N at time 8 (over strategies)",
  base = ess_min(models$base$ws, w90 - 1),
  race = ess_min(models$race$ws, w90 - 1))
for (q in c("rd_0_90_vs_0_30", "rd_30_90_vs_0_30")) {
  for (v in c("estimate", "boot_se")) {
    rows[[length(rows) + 1]] <- data.frame(
      quantity = paste(q, v),
      base = models$base$est[[v]][models$base$est$quantity == q],
      race = models$race$est[[v]][models$race$est$quantity == q])
  }
}
cmp <- do.call(rbind, rows)
save_csv(cmp, "wc_race_comparison.csv")
cat("\nBase model (age, sex, renal) vs race model",
    "(+ race_black, race_other):\n")
print(cmp, digits = 3, row.names = FALSE)
cat("\nRace model, null reference MaxIB (observed vs conditional null):\n")
print(nr_race$maxib[nr_race$maxib$feature %in% c("race_other", "(all)"), ],
      digits = 3, row.names = FALSE)
cat(sprintf("\nTotal runtime: %.1f s\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
