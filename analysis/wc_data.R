# wc_data.R --------------------------------------------------------------------
# Shared data preparation for the Webster-Clark SynPUF case study
# (analysis/03_empirical_webster_clark.R, 03b, 03c, 03d and part 2 of
# 02g_null_coupling_check.R). Source after pkgload::load_all(".").
#
# INPUTS (git-ignored):
#   data-raw/webster_clark/Data/baseccwcohort.sas7bdat
#     git clone https://github.com/mawcpharmdphd/CCW_tutorial
#       data-raw/webster_clark
#   data-raw/synpuf_claims/{carrier,outpatient,pde}.csv
#     bash data-raw/get_synpuf_claims.sh   (CMS DE-SynPUF claims of the
#     1,590 cohort members; first column DESYNPUF_ID = cohort ID)
#
# wc_prepare(delta = 10, horizon = 180) returns a list with
#   raw         the tutorial cohort (haven::read_sas)
#   cohort      one row per individual: id, ID, age, female, renal,
#               race_black, race_other, rx, fu, outcome, init_time,
#               event_time, ltfu_time (intervals of `delta` days)
#   long        ccw_clone() long data with the time-varying covariates and
#               the depth functions merged by (id, time); row order
#               strategy, id, time
#   covs        the five baseline covariates
#   tvf         the three time-varying covariates c("visit", "outpt",
#               "fill"); for interval k the window is the `delta` days
#               before the start of the interval (days delta*k - delta to
#               delta*k - 1 after discharge; for k = 0 the last days of the
#               index admission):
#                 visit  >= 1 carrier (physician/supplier) claim
#                        (CLM_FROM_DT) in the window
#                 outpt  >= 1 hospital outpatient claim (CLM_FROM_DT)
#                 fill   >= 1 Part D dispensing (SRVC_DT), excluding the
#                        NDCs the patient filled on the initiation date
#                        RX_date (the index antiplatelet)
#   depth       three history-depth functions of the same non-index Part D
#               dispensings, monitored by the diagnostics but never entered
#               into a censoring model, c("fill_lag", "fill_30",
#               "fill_cum"):
#                 fill_lag  `fill` of the window of interval k - 1 (0 for
#                           k = 0)
#                 fill_30   number of dispensings in the 30 days before the
#                           interval (days delta*k - 30 to delta*k - 1; for
#                           small k partly in the index admission)
#                 fill_cum  number of dispensings from discharge (day 0) to
#                           day delta*k - 1 (0 for k = 0)
#   feats       c(covs, tvf, depth), the eleven monitored history functions
#   strategies  0_30, 0_90, no_0_90
#   K, w30, w90, delta, horizon
#   tv          id x time table of the time-varying covariates and the
#               depth functions (with ID)
#   formulas    censoring-model formula lists base, race, use, full (terms
#               from covs and tvf only; the depth functions enter none)
#
# wc_boot_long(wc, s) builds the long data of a bootstrap resample.

wc_prepare <- function(delta = 10, horizon = 180) {
  data_file <- file.path("data-raw", "webster_clark", "Data",
                         "baseccwcohort.sas7bdat")
  if (!file.exists(data_file)) {
    stop("Data not found: ", data_file, "\nRun:\n  git clone ",
         "https://github.com/mawcpharmdphd/CCW_tutorial ",
         "data-raw/webster_clark")
  }
  cl_dir <- file.path("data-raw", "synpuf_claims")
  cl_files <- file.path(cl_dir, c("carrier.csv", "outpatient.csv",
                                  "pde.csv"))
  if (!all(file.exists(cl_files))) {
    stop("SynPUF claims not found: ",
         paste(cl_files[!file.exists(cl_files)], collapse = ", "),
         "\nRun:\n  bash data-raw/get_synpuf_claims.sh")
  }
  k_int <- horizon / delta

  # 1. Cohort -----------------------------------------------------------------
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
  day_to_int <- function(d) floor(d / delta)
  cohort$init_time <- ifelse(is.na(cohort$rx) | cohort$rx >= horizon, Inf,
                             day_to_int(cohort$rx))
  cohort$event_time <- ifelse(cohort$outcome == 1 & cohort$fu <= horizon,
                              pmin(day_to_int(cohort$fu), k_int - 1), Inf)
  cohort$ltfu_time <- ifelse(cohort$outcome == 0 & cohort$fu < horizon,
                             day_to_int(cohort$fu), Inf)

  covs <- c("age", "female", "renal", "race_black", "race_other")
  tvf <- c("visit", "outpt", "fill")
  depth <- c("fill_lag", "fill_30", "fill_cum")
  w30 <- day_to_int(30)
  w90 <- day_to_int(90)
  strategies <- list("0_30" = c(0, w30 - 1), "0_90" = c(0, w90 - 1),
                     "no_0_90" = c(NA, w90 - 1))
  long <- ccw_clone(cohort, strategies, n_time = k_int, covariates = covs)

  # 2. Time-varying covariates from the SynPUF claims -------------------------
  read_claims <- function(f, cols) {
    x <- utils::read.csv(file.path(cl_dir, f), colClasses = "character",
                         check.names = FALSE)
    x[, cols, drop = FALSE]
  }
  ymd <- function(x) as.Date(x, "%Y%m%d")
  # days since discharge; a claim on day d falls in the window of interval
  # k = floor(d / delta) + 1 (days delta*k - delta .. delta*k - 1)
  window_of <- function(id, dt) {
    d <- as.numeric(dt - raw$Discharge_date[match(id, raw$ID)])
    floor(d / delta) + 1
  }
  car <- read_claims("carrier.csv", c("DESYNPUF_ID", "CLM_FROM_DT"))
  opd <- read_claims("outpatient.csv", c("DESYNPUF_ID", "CLM_FROM_DT"))
  pde <- read_claims("pde.csv", c("DESYNPUF_ID", "SRVC_DT", "PROD_SRVC_ID"))
  pde_dt <- ymd(pde$SRVC_DT)
  # index antiplatelet: NDCs dispensed to the patient on RX_date
  rx_date <- raw$RX_date[match(pde$DESYNPUF_ID, raw$ID)]
  key <- paste(pde$DESYNPUF_ID, pde$PROD_SRVC_ID)
  idx_key <- unique(key[!is.na(rx_date) & pde_dt == rx_date])
  keep <- !(key %in% idx_key)
  any_in <- function(id, k) {
    ok <- !is.na(k) & k >= 0 & k <= k_int - 1
    hit <- unique(paste(id[ok], k[ok]))
    as.numeric(paste(tv$ID, tv$time) %in% hit)
  }
  tv <- data.frame(
    id = rep(cohort$id, k_int), ID = rep(cohort$ID, k_int),
    time = rep(seq_len(k_int) - 1, each = nrow(cohort))
  )
  tv$visit <- any_in(car$DESYNPUF_ID,
                     window_of(car$DESYNPUF_ID, ymd(car$CLM_FROM_DT)))
  tv$outpt <- any_in(opd$DESYNPUF_ID,
                     window_of(opd$DESYNPUF_ID, ymd(opd$CLM_FROM_DT)))
  tv$fill <- any_in(pde$DESYNPUF_ID[keep],
                    window_of(pde$DESYNPUF_ID[keep], pde_dt[keep]))
  tv <- tv[order(tv$id, tv$time), ]
  rownames(tv) <- NULL

  # history depth of the non-index dispensings (days since discharge)
  fill_day <- as.numeric(pde_dt[keep] -
                           raw$Discharge_date[match(pde$DESYNPUF_ID[keep],
                                                    raw$ID)])
  fill_id <- pde$DESYNPUF_ID[keep]
  # number of dispensings on days lo(k) .. delta*k - 1, per tv row
  n_fills <- function(lo) {
    out <- numeric(nrow(tv))
    for (k in seq_len(k_int) - 1) {
      rows <- which(tv$time == k)
      sel <- fill_day >= lo(k) & fill_day <= delta * k - 1
      out[rows] <- as.numeric(table(factor(fill_id[sel],
                                           levels = tv$ID[rows])))
    }
    out
  }
  n10 <- n_fills(function(k) delta * k - delta)
  if (!all(tv$fill == as.numeric(n10 > 0))) {
    stop("fill differs from the indicator of the window count.")
  }
  prev <- match(paste(tv$id, tv$time - 1), paste(tv$id, tv$time))
  tv$fill_lag <- ifelse(is.na(prev), 0, tv$fill[prev])
  tv$fill_30 <- n_fills(function(k) delta * k - 30)
  tv$fill_cum <- n_fills(function(k) 0)

  # merge onto the long data, keeping the ccw_clone() row order
  pos <- match(paste(long$id, long$time), paste(tv$id, tv$time))
  if (anyNA(pos)) stop("Unmatched id x time rows in the long data.")
  for (v in c(tvf, depth)) long[[v]] <- tv[[v]][pos]

  # 3. Censoring models -------------------------------------------------------
  base <- "age + female + renal"
  race <- "race_black + race_other"
  use <- "visit + outpt + fill"
  fml <- function(terms) {
    rhs <- paste(terms, collapse = " + ")
    list(no_init = stats::as.formula(paste("~", rhs, "+ time + I(time^2)")),
         grace_end = stats::as.formula(paste("~", rhs)))
  }
  formulas <- list(
    base = list(no_init = ~ age + female + renal + time + I(time^2),
                grace_end = ~ age + female + renal),
    race = fml(c(base, race)),
    use = fml(c(base, use)),
    full = fml(c(base, race, use))
  )

  list(raw = raw, cohort = cohort, long = long, covs = covs, tvf = tvf,
       depth = depth, feats = c(covs, tvf, depth), strategies = strategies,
       K = k_int, w30 = w30, w90 = w90, delta = delta, horizon = horizon,
       tv = tv, formulas = formulas)
}

# wc_boot_long(wc, s): long data of an individual bootstrap resample `s`
# (original ids, with replacement), built as in analysis/03 (cohort rows
# resampled, new ids 1..n, ccw_clone()) with the time-varying covariates of
# the original individuals (and the depth functions) attached by (original
# id, time).
wc_boot_long <- function(wc, s) {
  cb <- wc$cohort[s, ]
  cb$orig_id <- cb$id
  cb$id <- seq_len(nrow(cb))
  lb <- ccw_clone(cb, wc$strategies, n_time = wc$K,
                  covariates = c(wc$covs, "orig_id"))
  pos <- match(paste(lb$orig_id, lb$time), paste(wc$tv$id, wc$tv$time))
  for (v in c(wc$tvf, wc$depth)) lb[[v]] <- wc$tv[[v]][pos]
  lb$orig_id <- NULL
  lb
}
