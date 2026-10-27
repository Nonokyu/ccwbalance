# 02h_null_oracle.R -----------------------------------------------------------
# Pure-oracle check of the conditional null reference of the interval-RRB
# MaxIB: no estimation anywhere. With the TRUE data-generating censoring
# probabilities pi = P(D_{t+1} = 0 | H_t) (DGM of simulate_ccw(), correct
# scenarios P2-01, P1-01, P2-02 of paper_scenarios(), i.e. correct_n300_1000,
# correct_1000, correct_n8000_500 of the run largerep),
#   * the reference weights are W_{t-1} = prod_{s < t} 1 / pi_s (variant
#     "true_w") or 1 (variant "unit_w");
#   * the observed interval RRB uses the post weight w_ref / pi;
#   * the null redraws D* from pi with the coupled (within-individual) draw
#     of rrb_null() and keeps pi as the denominator (no refit).
# Given the risk sets and weights, the per-cell null is then exact, so the
# exceedance rate of the observed MaxIB over its null 95th percentile shows
# the conservatism that remains when the plug-in error of pi_hat is absent
# (the coupling of the null scale to the realised risk sets). The
# studentised maximum max |RRB| / null SD is reported alongside the raw
# MaxIB. Cells: g1 and g3 with at least one censored row, 8 history
# functions (X, L, L_sq, L_lagprod, L_diff, L_cummean, L_gt0, L_x_X).
#
# Outputs:
#   output/sim/null_oracle_results[_suffix].csv  one row per n x replicate x
#     weight variant: maxib_obs, null_q95, exceed_raw, stud_obs, stud_q95,
#     exceed_stud, cell_rej (share of cells with per-cell two-sided
#     p = (1 + #{|RRB*| >= |RRB|}) / (1 + B) < 0.05), cell_rej_g1,
#     cell_rej_g3late (g3 cells with time >= 4), n_cells
#   output/sim/null_oracle_summary[_suffix].csv  per n x variant: n_rep,
#     exceedance rates (raw, studentised) and per-cell rejection rates with
#     SEs (binomial for the exceedances; across-replicate SD / sqrt(n_rep)
#     for the cell shares), MC q95 of maxib_obs vs mean null q95, and the
#     correlation over replicates of null q95 and observed MaxIB
#
# Configuration (environment variables):
#   CCW_CORES       parallel workers (parallel::mclapply)   [detectCores()-1]
#   CCW_OUT_SUFFIX  output suffix ("" -> files above)                    [""]
#   CCW_SMOKE       1 = smoke test: 3 replicates per n, 20 draws, suffix
#                   "smoke"                                              [0]
# Replicates: 1000 (n = 300), 1000 (n = 2000), 500 (n = 8000); 500 null
# draws. Random numbers: one L'Ecuyer-CMRG stream per task (n x replicate);
# within a task both weight variants use the same null draws.
#
# Run from the project root:  Rscript analysis/02h_null_oracle.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) as.numeric(Sys.getenv(x, as.character(default)))
SMOKE <- Sys.getenv("CCW_SMOKE", "0") == "1"
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", if (SMOKE) "smoke" else "")
N_DRAWS <- if (SMOKE) 20 else 500
SEED <- 20261024
DIAG_STRATEGIES <- c("g1", "g3")
sfx <- if (SUFFIX == "") "" else paste0("_", SUFFIX)
out_dir <- file.path("output", "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x) file.path(out_dir, paste0(x, sfx, ".csv"))

all_feats <- c("X", "L", "L_sq", "L_lagprod", "L_diff", "L_cummean", "L_gt0",
               "L_x_X")
mechanisms <- list(
  oracle_n300 = list(scenario = paper_scenario("P2-01"), n_rep = 1000),
  oracle_n2000 = list(scenario = paper_scenario("P1-01"), n_rep = 1000),
  oracle_n8000 = list(scenario = paper_scenario("P2-02"), n_rep = 500)
)
if (SMOKE) for (m in names(mechanisms)) mechanisms[[m]]$n_rep <- 3

# true P(remain uncensored | H_t) from the DGM (as in 02i)
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

colmax <- function(m) {
  apply(m, 2, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
}

one_rep <- function(mech, rep) {
  sc <- mechanisms[[mech]]$scenario
  pars <- do.call(ccw_sim_params, sc$params)
  d <- simulate_ccw(pars)
  d <- as.data.frame(make_history_features(d, interactions = list(c("L",
                                                                    "X"))))
  d <- d[order(d$strategy, d$id, d$time), ]
  pt <- true_pi(d, pars)
  n <- nrow(d)
  clone <- paste(d$strategy, d$id)
  d$w_true <- stats::ave(1 / pt, clone, FUN = cumprod) * pt  # W_{t-1}
  d$w_unit <- 1
  poss <- as.logical(d$censor_possible)
  dd <- d$censored_next
  key <- paste(d$id[poss], d$time[poss], sep = "\r")
  grp <- match(key, unique(key))
  n_grp <- max(grp)
  on0 <- .null_orient(d, "strategy", "censored_next", "censor_possible",
                      "A")[poss]
  null_seed <- sample.int(.Machine$integer.max, 1)
  out <- list()
  for (wv in c("true_w", "unit_w")) {
    prep <- .null_prep(d, all_feats, "strategy", "time",
                       if (wv == "true_w") "w_true" else "w_unit")
    ro <- .null_rrb(prep, dd, pt)
    ncens <- rowsum(dd, prep$cell, reorder = FALSE)[, 1]
    st <- sub("\r.*", "", prep$cells)
    tm <- as.numeric(sub(".*\r", "", prep$cells))
    cells <- expand.grid(cell = seq_along(prep$cells),
                         f = seq_along(all_feats))
    keep <- ncens[cells$cell] >= 1 & st[cells$cell] %in% DIAG_STRATEGIES &
      !is.na(ro[cbind(cells$cell, cells$f)])
    cells <- cells[keep, ]
    pos <- cbind(cells$cell, cells$f)
    sims <- matrix(NA_real_, nrow(cells), N_DRAWS)
    set.seed(null_seed)  # common random numbers across the variants
    for (b in seq_len(N_DRAWS)) {
      ds <- integer(n)
      ds[poss] <- .null_draw(pt[poss], on0, grp, stats::runif(n_grp))
      sims[, b] <- .null_rrb(prep, ds, pt)[pos]
    }
    obs <- ro[pos]
    sd_null <- apply(sims, 1, stats::sd, na.rm = TRUE)
    sd_null[!is.finite(sd_null) | sd_null <= 0] <- NA_real_
    m_null <- colmax(abs(sims))
    z_null <- colmax(abs(sims) / sd_null)
    q95 <- stats::quantile(m_null, 0.95, na.rm = TRUE, names = FALSE)
    z_obs <- max(abs(obs) / sd_null, na.rm = TRUE)
    z_q95 <- stats::quantile(z_null, 0.95, na.rm = TRUE, names = FALSE)
    p_cell <- (1 + rowSums(abs(sims) >= abs(obs), na.rm = TRUE)) /
      (1 + rowSums(!is.na(sims)))
    s_c <- st[cells$cell]
    late <- s_c == "g3" & tm[cells$cell] >= 4
    out[[wv]] <- data.frame(
      mechanism = mech, rep = rep, n = pars$n, weights = wv,
      n_cells = nrow(cells), maxib_obs = max(abs(obs)), null_q95 = q95,
      exceed_raw = as.numeric(max(abs(obs)) > q95),
      stud_obs = z_obs, stud_q95 = z_q95,
      exceed_stud = as.numeric(z_obs > z_q95),
      cell_rej = mean(p_cell < 0.05),
      cell_rej_g1 = mean(p_cell[s_c == "g1"] < 0.05),
      cell_rej_g3late = if (any(late)) mean(p_cell[late] < 0.05) else NA_real_
    )
  }
  do.call(rbind, out)
}

tasks <- do.call(rbind, lapply(names(mechanisms), function(m) {
  data.frame(mech = m, rep = seq_len(mechanisms[[m]]$n_rep))
}))
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
stream <- .Random.seed
seeds <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  stream <- parallel::nextRNGStream(stream)
  seeds[[i]] <- stream
}
message(sprintf("Running %d tasks, %d draws, %d cores", nrow(tasks), N_DRAWS,
                N_CORES))
t_start <- Sys.time()
ord <- order(-match(tasks$mech, names(mechanisms)), tasks$rep)  # n = 8000 first
res_list <- parallel::mclapply(ord, function(i) {
  assign(".Random.seed", seeds[[i]], envir = globalenv())
  tryCatch(one_rep(tasks$mech[i], tasks$rep[i]), error = function(e) {
    data.frame(mechanism = tasks$mech[i], rep = tasks$rep[i],
               error = conditionMessage(e))
  })
}, mc.cores = N_CORES, mc.preschedule = FALSE)
failed <- Filter(function(x) "error" %in% names(x), res_list)
if (length(failed) > 0) {
  warning(length(failed), " replicate(s) failed; first: ", failed[[1]]$error)
}
res <- do.call(rbind, Filter(function(x) !("error" %in% names(x)), res_list))
res <- res[order(match(res$mechanism, names(mechanisms)), res$rep,
                 res$weights), ]
utils::write.csv(res, out_csv("null_oracle_results"), row.names = FALSE)

rate_se <- function(v) c(mean(v), sqrt(mean(v) * (1 - mean(v)) / length(v)))
mean_se <- function(v) {
  v <- v[!is.na(v)]
  c(mean(v), stats::sd(v) / sqrt(length(v)))
}
summ <- do.call(rbind, lapply(split(res, list(res$mechanism, res$weights),
                                    drop = TRUE), function(x) {
  er <- rate_se(x$exceed_raw)
  es <- rate_se(x$exceed_stud)
  c0 <- mean_se(x$cell_rej)
  c1 <- mean_se(x$cell_rej_g1)
  c3 <- mean_se(x$cell_rej_g3late)
  data.frame(
    mechanism = x$mechanism[1], n = x$n[1], weights = x$weights[1],
    n_rep = nrow(x), exceed_raw = er[1], exceed_raw_se = er[2],
    exceed_stud = es[1], exceed_stud_se = es[2],
    cell_rej = c0[1], cell_rej_se = c0[2],
    cell_rej_g1 = c1[1], cell_rej_g1_se = c1[2],
    cell_rej_g3late = c3[1], cell_rej_g3late_se = c3[2],
    mc_q95_obs = stats::quantile(x$maxib_obs, 0.95, names = FALSE),
    null_q95_mean = mean(x$null_q95),
    cor_q95_obs = stats::cor(x$null_q95, x$maxib_obs),
    mc_q95_stud = stats::quantile(x$stud_obs, 0.95, names = FALSE),
    stud_q95_mean = mean(x$stud_q95)
  )
}))
summ <- summ[order(summ$n, summ$weights), ]
rownames(summ) <- NULL
utils::write.csv(summ, out_csv("null_oracle_summary"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
cat(sprintf("Runtime: %.1f min; failed: %d\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            length(failed)))
cat("Written:", out_csv("null_oracle_results"), out_csv("null_oracle_summary"),
    sep = "\n  ")
