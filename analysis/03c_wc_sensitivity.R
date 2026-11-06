# 03c_wc_sensitivity.R --------------------------------------------------------
# Case study: individual bootstrap of the RRB under the tutorial (base) model,
# and RRB, weight support, null reference and risk difference under the
# censoring models race, use and full; comparison table (Tables 5 and S7).
# Inputs: the case-study data (analysis/wc_data.R; see analysis/README.md).
# Outputs: output/empirical/wc_rrb_boot_*.csv, wc_<model>_*.csv,
#   wc_model_comparison.csv
# Run: Rscript analysis/03c_wc_sensitivity.R
#   (WC_RRB_BOOT, WC_N_BOOT, NULL_DRAWS)

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "wc_data.R"))

RRB_BOOT <- as.numeric(Sys.getenv("WC_RRB_BOOT", "200"))
N_BOOT <- as.numeric(Sys.getenv("WC_N_BOOT", "200"))
N_DRAWS <- as.numeric(Sys.getenv("NULL_DRAWS", "500"))
out_dir <- file.path("output", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}
secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))
t_start <- Sys.time()

# 1. Data and base censoring model --------------------------------------------
wc <- wc_prepare(delta = 10, horizon = 180)
long <- wc$long
covs <- wc$covs
feats <- wc$feats
strategies <- wc$strategies
K <- wc$K
w30 <- wc$w30
w90 <- wc$w90
w <- fit_censoring_model(long, wc$formulas$base, stabilized = TRUE)

# base model: individual bootstrap of the weighted RRB
for (sc in c("cumulative", "interval")) {
  t0 <- Sys.time()
  fs <- if (sc == "interval") feats else covs
  bt <- rrb_boot(w, fs, B = RRB_BOOT, scope = sc, seed = 20261018)
  save_csv(bt, sprintf("wc_rrb_boot_%s.csv", sc))
  cat(sprintf("rrb_boot (%s, B = %d): %.1f s\n", sc, RRB_BOOT, secs(t0)))
}

# 2. Diagnostics and estimates per censoring model ----------------------------
rd <- function(e) {
  r <- stats::setNames(e$risk[e$time == max(e$time)],
                       as.character(e$strategy[e$time == max(e$time)]))
  c(rd_0_90_vs_0_30 = r[["0_90"]] - r[["0_30"]])
}
# 180-day risks and risk difference; the seed of analysis/03 gives the same
# resamples for every model
estimates <- function(wf, fml) {
  est <- estimate_risk(wf)
  risk_h <- est[est$time == K - 1, c("strategy", "risk")]
  rd_hat <- rd(est)
  set.seed(20260923)
  boot <- t(vapply(seq_len(N_BOOT), function(b) {
    s <- sample(wc$cohort$id, replace = TRUE)
    lb <- wc_boot_long(wc, s)
    wb <- suppressWarnings(fit_censoring_model(lb, fml))
    eb <- estimate_risk(wb)
    c(stats::setNames(eb$risk[eb$time == K - 1],
                      as.character(eb$strategy[eb$time == K - 1])), rd(eb))
  }, numeric(4)))
  tab <- data.frame(
    quantity = c(paste0("risk_", names(strategies)), names(rd_hat)),
    estimate = c(risk_h$risk, rd_hat),
    boot_se = apply(boot, 2, stats::sd)
  )
  tab$lcl <- tab$estimate - 1.96 * tab$boot_se
  tab$ucl <- tab$estimate + 1.96 * tab$boot_se
  tab
}
run_model <- function(m, save = TRUE) {
  t0 <- Sys.time()
  fml <- wc$formulas[[m]]
  wm <- if (m == "base") w else fit_censoring_model(long, fml,
                                                     stabilized = TRUE)
  r <- rbind(
    rrb(wm, feats, method = c("unweighted", "weighted")),
    rrb(wm, covs, method = c("unweighted", "weighted"), scope = "cumulative")
  )
  ws <- weight_support(wm)
  nr <- rrb_null(wm, feats, n_draws = N_DRAWS, seed = 20260925)
  est <- estimates(wm, fml)
  if (save) {
    save_csv(r, sprintf("wc_%s_rrb.csv", m))
    save_csv(ws, sprintf("wc_%s_weight_support.csv", m))
    save_csv(nr$maxib, sprintf("wc_%s_null_maxib.csv", m))
    save_csv(nr$maxib_stud, sprintf("wc_%s_null_maxib_stud.csv", m))
    save_csv(nr$cells, sprintf("wc_%s_null_cells.csv", m))
    save_csv(est, sprintf("wc_%s_estimates.csv", m))
  }
  cat(sprintf("model %s (null %d draws, bootstrap B = %d): %.1f s\n", m,
              N_DRAWS, N_BOOT, secs(t0)))
  list(r = r, ws = ws, nr = nr, est = est)
}
model_names <- c("base", "race", "use", "full")
# base model: outputs already written by analysis/03 and 03b
res <- lapply(stats::setNames(model_names, model_names), function(m) {
  run_model(m, save = m != "base")
})

# 3. Comparison table ---------------------------------------------------------
last_t <- c("0_30" = w30 - 1, "0_90" = w90 - 1, "no_0_90" = w90 - 1)
pick_rrb <- function(r, sc, g, tt, f) {
  v <- r$rrb[r$method == "weighted" & r$scope == sc &
               as.character(r$strategy) == g & r$time == tt & r$feature == f]
  if (length(v) == 1) v else NA_real_
}
pick_maxib <- function(nr, f, col) {
  nr$maxib[[col]][nr$maxib$strategy == "(all)" & nr$maxib$feature == f]
}
pick_stud <- function(nr, f, col) {
  nr$maxib_stud[[col]][nr$maxib_stud$strategy == "(all)" &
                         nr$maxib_stud$feature == f]
}
pick_est <- function(est, q, col) est[[col]][est$quantity == q]
row <- function(label, fun) {
  cbind(data.frame(quantity = label),
        as.data.frame(lapply(res, fun)))
}
rows <- list()
for (f in feats) {
  rows[[length(rows) + 1]] <- row(
    sprintf("interval RRB %s, 0_30, time %d (day 30)", f, w30 - 1),
    function(x) pick_rrb(x$r, "interval", "0_30", w30 - 1, f))
  rows[[length(rows) + 1]] <- row(
    sprintf("interval RRB %s, 0_90, time %d (day 90)", f, w90 - 1),
    function(x) pick_rrb(x$r, "interval", "0_90", w90 - 1, f))
}
for (g in names(strategies)) {
  rows[[length(rows) + 1]] <- row(
    sprintf("cumulative RRB race_other, %s, time %d (grace end)", g,
            last_t[[g]]),
    function(x) pick_rrb(x$r, "cumulative", g, last_t[[g]], "race_other"))
}
for (f in c(feats, "(all)")) {
  lab <- if (f != "(all)") sprintf("MaxIB %s (all strategies)", f) else
    "MaxIB all features and strategies"
  for (col in c("maxib_obs", "null_q95", "p")) {
    rows[[length(rows) + 1]] <- row(
      paste(lab, sub("maxib_obs", "observed", col)),
      function(x) pick_maxib(x$nr, f, col))
  }
}
rows[[length(rows) + 1]] <- row(
  "smallest ESS/|R| at time 8 (over strategies)",
  function(x) min(x$ws$ess_ratio[x$ws$time == w90 - 1], na.rm = TRUE))
rows[[length(rows) + 1]] <- row(
  "ESS/|R| at time 8, 0_90",
  function(x) x$ws$ess_ratio[x$ws$time == w90 - 1 &
                               x$ws$strategy == "0_90"])
rows[[length(rows) + 1]] <- row(
  "w99 at time 8, 0_90",
  function(x) x$ws$w99[x$ws$time == w90 - 1 & x$ws$strategy == "0_90"])
for (g in names(strategies)) {
  q <- paste0("risk_", g)
  rows[[length(rows) + 1]] <- row(paste(q, "(180 d)"),
                                  function(x) pick_est(x$est, q, "estimate"))
}
q <- "rd_0_90_vs_0_30"
for (col in c("estimate", "boot_se", "lcl", "ucl")) {
  rows[[length(rows) + 1]] <- row(paste(q, col),
                                  function(x) pick_est(x$est, q, col))
}
# studentised maximum (appended after the existing rows)
for (f in c(feats, "(all)")) {
  lab <- if (f != "(all)") {
    sprintf("Studentised max %s (all strategies)", f)
  } else {
    "Studentised max all features and strategies"
  }
  for (col in c("stud_obs", "null_q95", "p")) {
    rows[[length(rows) + 1]] <- row(
      paste(lab, sub("stud_obs", "observed", col)),
      function(x) pick_stud(x$nr, f, col))
  }
}
cmp <- do.call(rbind, rows)
rownames(cmp) <- NULL
save_csv(cmp, "wc_model_comparison.csv")
cat("\nCensoring models: base (age, sex, renal), race (+ race_black,",
    "race_other),\nuse (+ visit, outpt, fill), full (+ race and use):\n")
print(cmp, digits = 3, row.names = FALSE)

cat("\nWeighted interval RRB at the grace ends (0_30 time 2, 0_90 time 8):\n")
for (m in model_names) {
  r <- res[[m]]$r
  ge <- r[r$scope == "interval" & r$method == "weighted" &
            ((r$strategy == "0_30" & r$time == w30 - 1) |
               (r$strategy == "0_90" & r$time == w90 - 1)), ]
  tab <- stats::xtabs(rrb ~ feature + strategy, data = ge,
                      drop.unused.levels = TRUE)
  cat("--", m, "\n")
  print(round(unclass(tab)[feats, c("0_30", "0_90")], 3))
}
stud_feats <- c("fill", wc$depth, "race_other", "outpt", "(all)")
for (m in model_names) {
  cat(sprintf("\nModel %s: null reference MaxIB (all strategies):\n", m))
  nr <- res[[m]]$nr
  print(nr$maxib[nr$maxib$strategy == "(all)", ], digits = 3,
        row.names = FALSE)
  cat(sprintf("Model %s: studentised maximum max|z| (selected functions):\n",
              m))
  print(nr$maxib_stud[nr$maxib_stud$feature %in% stud_feats, ], digits = 3,
        row.names = FALSE)
  big <- nr$cells[!is.na(nr$cells$z) & abs(nr$cells$z) > 2, ]
  cat("cells with |z| > 2:\n")
  print(big[, c("strategy", "time", "feature", "n_censored", "rrb_obs",
                "null_sd", "z", "p_two")], digits = 3, row.names = FALSE)
}
cat(sprintf("\nTotal runtime: %.1f s\n", secs(t_start)))
