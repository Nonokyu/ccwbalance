# 02i_null_decomposition.R ----------------------------------------------------
# Why is the conditional parametric null reference (rrb_null()) conservative
# at n = 300? Decomposes the null 95th percentile (q95) of the interval-RRB
# MaxIB (g1, g3; all history functions) into the contribution of
#   * the censoring probabilities the null indicators are drawn from:
#     the fitted pi_hat (the package null) or the TRUE data-generating
#     pi = P(D_{t+1} = 0 | H_t) (oracle), and
#   * refitting the censoring model on each draw (refit) or keeping the
#     drawing probabilities as the weights' denominators (no refit),
# i.e. four schemes {fit, true} x {refit, norefit}; "fit_refit" is the
# package null (rrb_null(), default within-individual coupling).
#
# Mechanisms (specs of the manuscript scenarios, paper_scenario() in
# analysis/sim_scenarios.R; own replicates):
#   many_cov_n300  P2-05 (8 extra baseline covariates X2-X9 in the model;
#                  16 history functions; n = 300)
#   correct_n300   P2-01 (base DGM, correct model; 8 history functions;
#                  n = 300)
#   correct, correct_n8000  P1-01, P2-02: the same mechanism at n = 2000 and
#                  8000
# (The mechanisms were taken from the run main's correct_n300, correct and
# correct_n8000 before the scenario ids existed; these have the same
# parameters, censoring model, truncation and monitored functions as P2-01,
# P1-01 and P2-02, so the random numbers and results are unchanged.)
# For every replicate: simulate, fit the scenario's censoring model, record
# the observed MaxIB (with pi_hat, as in the simulation study) and the null
# q95 under each scheme. All four schemes use the same null seed (common
# random numbers) and the same coupled redraw as rrb_null(). The true pi is
# computed from the DGM equations of simulate_ccw() (ccw_sim_params()): with
# p_init = plogis(beta0 + beta1 L_t + beta2 L_t^2 + beta3 X + beta5 L_t X +
# beta6 L_{t-1} + beta_extra' X_extra), pi = 1 - p_init in "no_init" rows
# (censored on initiation) and pi = p_init in "grace_end" rows (simulate_ccw()
# does not return the true probability).
#
# Outputs:
#   output/sim/null_decomposition.csv          one row per mechanism x
#                                              replicate x scheme
#   output/sim/null_decomposition_summary.csv  per mechanism x scheme: mean
#     null q95, exceedance rate of the observed MaxIB against it (with
#     binomial SE), MC q95 of the observed MaxIB, and the mean over
#     replicates of the SD of logit(pi_hat) and of logit(pi_true) over the
#     censorable rows
#
# Configuration (environment variables):
#   CCW_N_REPS        replicates per mechanism                 [200]
#   CCW_N_NULL_DRAWS  null draws per replicate and scheme      [300]
#   CCW_CORES         parallel workers (parallel::mclapply)    [detectCores()-1]
#   CCW_OUT_SUFFIX    output suffix ("" -> files above)        [""]
# Random numbers: one L'Ecuyer-CMRG stream per task (mechanism x replicate).
#
# Run from the project root:  Rscript analysis/02i_null_decomposition.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "sim_scenarios.R"))

env_num <- function(x, default) as.numeric(Sys.getenv(x, as.character(default)))
N_REPS <- env_num("CCW_N_REPS", 200)
N_DRAWS <- env_num("CCW_N_NULL_DRAWS", 300)
N_CORES <- env_num("CCW_CORES", max(1, parallel::detectCores() - 1))
SUFFIX <- Sys.getenv("CCW_OUT_SUFFIX", "")
SEED <- 20261018
DIAG_STRATEGIES <- c("g1", "g3")
out_dir <- file.path("output", "sim")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- function(x) {
  file.path(out_dir, sprintf("%s%s.csv", x,
                             if (SUFFIX == "") "" else paste0("_", SUFFIX)))
}

all_feats <- c("X", "L", "L_sq", "L_lagprod", "L_diff", "L_cummean", "L_gt0",
               "L_x_X")
mechanisms <- list(
  many_cov_n300 = paper_scenario("P2-05"),
  correct_n300 = paper_scenario("P2-01"),
  # appended 2026-10-19 (larger n; the task order and hence the RNG streams
  # of the two mechanisms above are unchanged)
  correct = paper_scenario("P1-01"),
  correct_n8000 = paper_scenario("P2-02")
)
feats_of <- function(sc) c(all_feats, as.character(sc$extra_feats))

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

# null q95 of MaxIB (g1/g3 cells, all features) for the four schemes, with
# the internals of rrb_null() (coupling = "individual", pi_floor = 1e-6)
null_schemes <- function(w, feats, formula, pi_true, n_draws, seed) {
  data <- as.data.frame(w)
  data$.w_post_null <- data$w_ref / data$pi_hat
  obs <- rrb(data, feats, method = "weighted", scope = "interval",
             w_ref = "w_ref", w_post = ".w_post_null")
  obs <- obs[obs$n_censored >= 1 & !is.na(obs$rrb), , drop = FALSE]
  sel <- as.character(obs$strategy) %in% DIAG_STRATEGIES
  prep <- .null_prep(data, feats, "strategy", "time", "w_ref")
  groups <- .null_groups(data, formula, c("strategy", "phase"),
                         "censor_possible")
  poss <- as.logical(data$censor_possible)
  n <- nrow(data)
  key <- paste(data$id[poss], data$time[poss], sep = "\r")
  grp_poss <- match(key, unique(key))
  n_grp <- max(c(0L, grp_poss))
  on0_poss <- .null_orient(data, "strategy", "censored_next",
                           "censor_possible", "A")[poss]
  pos <- cbind(match(paste(as.character(obs$strategy), obs$time, sep = "\r"),
                     prep$cells),
               match(as.character(obs$feature), feats))[sel, , drop = FALSE]
  p_src <- list(fit = data$pi_hat, true = pi_true)
  schemes <- expand.grid(refit = c("refit", "norefit"), src = c("fit", "true"),
                         stringsAsFactors = FALSE)
  m <- matrix(NA_real_, nrow(schemes), n_draws)
  set.seed(seed)
  for (b in seq_len(n_draws)) {
    u <- stats::runif(n_grp)  # shared by the four schemes
    for (k in seq_len(nrow(schemes))) {
      p <- p_src[[schemes$src[k]]]
      d_star <- integer(n)
      d_star[poss] <- .null_draw(p[poss], on0_poss, grp_poss, u)
      pi_star <- if (schemes$refit[k] == "refit") {
        .null_refit(groups, 1 - d_star, n)
      } else p
      r <- .null_rrb(prep, d_star, pmax(pi_star, 1e-6))[pos]
      m[k, b] <- if (all(is.na(r))) NA_real_ else max(abs(r), na.rm = TRUE)
    }
  }
  data.frame(scheme = paste(schemes$src, schemes$refit, sep = "_"),
             obs = max(abs(obs$rrb[sel])),
             q95 = apply(m, 1, stats::quantile, 0.95, na.rm = TRUE,
                         names = FALSE),
             n_cells = sum(sel))
}

one_rep <- function(mech, rep) {
  sc <- mechanisms[[mech]]
  pars <- do.call(ccw_sim_params, sc$params)
  d <- simulate_ccw(pars)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, sc$formula, truncate = sc$truncate)
  pt <- true_pi(w, pars)
  seed <- sample.int(.Machine$integer.max, 1)
  out <- null_schemes(w, feats_of(sc), sc$formula, pt, N_DRAWS, seed)
  poss <- as.logical(w$censor_possible)
  lg <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  cbind(data.frame(mechanism = mech, rep = rep, n = pars$n,
                   n_feats = length(feats_of(sc))),
        out,
        sd_logit_pihat = stats::sd(lg(w$pi_hat[poss])),
        sd_logit_pitrue = stats::sd(lg(pt[poss])),
        min_pihat = min(w$pi_hat[poss]), min_pitrue = min(pt[poss]))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
tasks <- expand.grid(rep = seq_len(N_REPS), mech = names(mechanisms),
                     stringsAsFactors = FALSE)
stream <- .Random.seed
seeds <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  stream <- parallel::nextRNGStream(stream)
  seeds[[i]] <- stream
}
message(sprintf(
  "Running %d tasks (%d mechanisms x %d reps), %d draws, %d cores",
  nrow(tasks), length(mechanisms), N_REPS, N_DRAWS, N_CORES
))
t_start <- Sys.time()
res_list <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
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
res$exceed <- as.numeric(res$obs > res$q95)
write.csv(res, out_csv("null_decomposition"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, list(res$mechanism, res$scheme),
                                    drop = TRUE), function(x) {
  e <- mean(x$exceed)
  data.frame(
    mechanism = x$mechanism[1], scheme = x$scheme[1], n_rep = nrow(x),
    n_feats = x$n_feats[1], null_q95_mean = mean(x$q95),
    exceed_rate = e, exceed_rate_se = sqrt(e * (1 - e) / nrow(x)),
    mc_q95_obs = stats::quantile(x$obs, 0.95, names = FALSE),
    obs_mean = mean(x$obs),
    sd_logit_pihat_mean = mean(x$sd_logit_pihat),
    sd_logit_pitrue_mean = mean(x$sd_logit_pitrue),
    min_pihat_mean = mean(x$min_pihat), min_pitrue_mean = mean(x$min_pitrue)
  )
}))
summ$scheme <- factor(summ$scheme, levels = c("fit_refit", "fit_norefit",
                                              "true_refit", "true_norefit"))
summ <- summ[order(summ$mechanism, summ$scheme), ]
write.csv(summ, out_csv("null_decomposition_summary"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
cat(sprintf("Runtime: %.1f min; failed: %d\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            length(failed)))
cat("Written:", out_csv("null_decomposition"),
    out_csv("null_decomposition_summary"), sep = "\n  ")
