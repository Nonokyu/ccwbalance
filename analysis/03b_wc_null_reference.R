# 03b_wc_null_reference.R ----------------------------------------------------
# Conditional parametric null reference for the interval RRB / MaxIB
# (rrb_null()) applied to the Webster-Clark SynPUF example (analysis/03).
#
# Data preparation and censoring models are identical to
# analysis/03_empirical_webster_clark.R (sections 1-2; no bootstrap).
#
# Outputs (git-ignored):
#   output/empirical/wc_rrb_null_cells.csv, wc_rrb_null_maxib.csv
#
# Run from the project root:  Rscript analysis/03b_wc_null_reference.R
#   NULL_DRAWS (default 500)

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
data_file <- file.path("data-raw", "webster_clark", "Data",
                       "baseccwcohort.sas7bdat")
if (!file.exists(data_file)) {
  stop("Data not found: ", data_file, "\nRun:\n  git clone ",
       "https://github.com/mawcpharmdphd/CCW_tutorial data-raw/webster_clark")
}
suppressMessages(pkgload::load_all(".", quiet = TRUE))

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

# 1. Webster-Clark cohort and censoring model (as in analysis/03) -------------
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

# 2. Null reference ------------------------------------------------------------
t0 <- Sys.time()
r0 <- rrb_null(w, covs, n_draws = N_DRAWS, seed = 20260925)
t_wc <- difftime(Sys.time(), t0, units = "secs")
save_csv(r0$cells, "wc_rrb_null_cells.csv")
save_csv(r0$maxib, "wc_rrb_null_maxib.csv")

show_cells <- function(r, g, tt) {
  x <- r$cells[as.character(r$cells$strategy) == g & r$cells$time == tt,
               c("strategy", "time", "feature", "n_risk", "n_censored",
                 "rrb_obs", "null_sd", "z", "p_two")]
  print(x, digits = 3, row.names = FALSE)
}
frac_z <- function(r) {
  z <- r$cells$z[!is.na(r$cells$z)]
  sprintf("%d / %d = %.3f", sum(abs(z) > 1.96), length(z),
          mean(abs(z) > 1.96))
}

cat(sprintf("\n=== Webster-Clark: rrb_null(), %d draws (%.1f s) ===\n",
            N_DRAWS, as.numeric(t_wc)))
cat("\nInterval ending day 90 (time 8), strategy 0_90:\n")
show_cells(r0, "0_90", 8)
cat("\nInterval ending day 90 (time 8), strategy 30_90:\n")
show_cells(r0, "30_90", 8)
cat("\nInterval ending day 30 (time 2), strategy 0_30:\n")
show_cells(r0, "0_30", 2)
cat("\nMaxIB (observed vs conditional null):\n")
print(r0$maxib, digits = 3, row.names = FALSE)
cat("\nCells (all strategies/features, with censoring) with |z| > 1.96:",
    frac_z(r0), "\n")
cat("\nAll cells with |z| > 2:\n")
big <- r0$cells[!is.na(r0$cells$z) & abs(r0$cells$z) > 2, ]
print(big[, c("strategy", "time", "feature", "n_risk", "n_censored",
              "rrb_obs", "null_sd", "z", "p_two")],
      digits = 3, row.names = FALSE)

cat(sprintf("\nTotal runtime: %.1f s\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
