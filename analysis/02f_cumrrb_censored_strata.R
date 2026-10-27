# 02f_cumrrb_censored_strata.R ------------------------------------------------
# Cumulative RRB of the baseline functions (X, L_0) summarised over the strata
# with artificial censoring (strategy x interval with at least one
# artificially censored individual; g1, g3), the same strata as interval RRB.
# 02_simulation_study.R took the cumulative-RRB maximum over every interval
# from the first artificial censoring to the end of follow-up. Without loss to
# follow-up the two maxima coincide, because the weights do not change after
# the last artificial censoring; with unweighted loss to follow-up they differ.
#
# The script reproduces the replicate data of the run main exactly (same
# L'Ecuyer-CMRG streams as 02_simulation_study.R; data simulation is the
# first use of each stream), refits the censoring model
# and recomputes both maxima for extra_ltfu (P3-03) and, as a check, the
# first CCW_CHECK_REPS replicates of truncation_95 (P2-04) and of the run's
# own correct scenario `correct` (the manuscript reports P1-01 from the run
# largerep). It checks that the old maximum equals cumrrb_base_maxib in
# output/sim/sim_results_main.csv. The column `scenario` of the outputs holds
# the run-main scenario names.
#
# Inputs:  output/sim/sim_results_main.csv
# Outputs: output/sim/sim_paper_cumrrb_censored_strata.csv (per replicate),
#          output/sim/sim_paper_cumrrb_censored_strata_summary.csv
#
# Run from the project root:  Rscript analysis/02f_cumrrb_censored_strata.R
# (CCW_CORES sets the number of cores; CCW_CHECK_REPS the replicates of the
# scenarios without loss to follow-up that are recomputed as a check)

pkgload::load_all(".", quiet = TRUE)
source(file.path("analysis", "sim_scenarios.R"))

# seed offset of the run main (equals SEED_OFFSET[["main"]] in
# 02_simulation_study.R)
SEED_OFFSET_MAIN <- 3
SEED <- 20260923
N_MC <- 5e5
N_SUBJ <- 2000
CONTRAST <- c("g3", "g1")
DIAG_STRATEGIES <- c("g1", "g3")
base_feats <- c("X", "L_base")
n_cores <- as.numeric(Sys.getenv("CCW_CORES",
                                 max(1, parallel::detectCores() - 1)))
check_reps <- as.numeric(Sys.getenv("CCW_CHECK_REPS", 20))

scenarios <- ccw_scenarios("main")
n_of <- function(sc) if (is.null(sc$params[["n"]])) N_SUBJ else sc$params[["n"]]
reps_of <- function(sc) {
  if (!is.null(sc$n_rep)) return(sc$n_rep)
  if (n_of(sc) >= 8000) 100 else 200
}
nboot_of <- function(sc) if (!is.null(sc$n_boot)) sc$n_boot else 100

# streams as in 02_simulation_study.R (the truths consume the initial stream)
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED + SEED_OFFSET_MAIN)
invisible(lapply(scenarios, function(sc) {
  true_risk_ccw(do.call(ccw_sim_params, sc$params), n_mc = N_MC)
}))
tk <- do.call(rbind, lapply(seq_along(scenarios), function(s) {
  data.frame(s = s, rep = seq_len(reps_of(scenarios[[s]])))
}))
cost <- vapply(scenarios, function(sc) n_of(sc) * (nboot_of(sc) + 10),
               numeric(1))
tasks <- tk[order(-cost[tk$s], tk$rep), ]
stream <- .Random.seed
seeds <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  stream <- parallel::nextRNGStream(stream)
  seeds[[i]] <- stream
}

scen_names <- vapply(scenarios, function(sc) sc$name, character(1))
tasks$scenario <- scen_names[tasks$s]
run <- tasks$scenario == "extra_ltfu" |
  (tasks$scenario %in% c("correct", "truncation_95") &
     tasks$rep <= check_reps)

max_abs <- function(x) {
  x <- x[as.character(x$strategy) %in% DIAG_STRATEGIES &
           as.character(x$feature) %in% base_feats & !is.na(x$rrb), ]
  if (nrow(x) == 0) NA_real_ else max(abs(x$rrb))
}

res <- parallel::mclapply(which(run), function(i) {
  assign(".Random.seed", seeds[[i]], envir = globalenv())
  sc <- scenarios[[tasks$s[i]]]
  pars <- do.call(ccw_sim_params, sc$params)
  if (is.null(sc$params[["n"]])) pars$n <- N_SUBJ
  d <- simulate_ccw(pars)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  w <- fit_censoring_model(d, sc$formula, truncate = sc$truncate)
  rc <- rrb(w, base_feats, scope = "cumulative")
  ri <- rrb(w, "X")
  cens <- paste(ri$strategy, ri$time)[ri$n_censored > 0]
  rc_cens <- rc[paste(rc$strategy, rc$time) %in% cens, ]
  data.frame(scenario = sc$name, rep = tasks$rep[i],
             cumrrb_base_maxib_followup = max_abs(rc[rc$n_censored > 0, ]),
             cumrrb_base_maxib = max_abs(rc_cens))
}, mc.cores = n_cores, mc.preschedule = FALSE)
res <- do.call(rbind, res)

old <- utils::read.csv("output/sim/sim_results_main.csv")
k <- match(paste(res$scenario, res$rep), paste(old$scenario, old$rep))
res$stored <- old$cumrrb_base_maxib[k]
dev <- max(abs(res$cumrrb_base_maxib_followup - res$stored))
message(sprintf("max |recomputed - stored| (old definition): %.2e", dev))
stopifnot(dev < 1e-10)
res$stored <- NULL

summ <- do.call(rbind, lapply(split(res, res$scenario), function(x) {
  data.frame(
    scenario = x$scenario[1], n_rep = nrow(x),
    cumrrb_base_maxib_mean = mean(x$cumrrb_base_maxib),
    cumrrb_base_maxib_mcse = stats::sd(x$cumrrb_base_maxib) / sqrt(nrow(x)),
    cumrrb_base_maxib_followup_mean = mean(x$cumrrb_base_maxib_followup),
    max_abs_diff = max(abs(x$cumrrb_base_maxib -
                             x$cumrrb_base_maxib_followup))
  )
}))
print(summ, digits = 4)
utils::write.csv(res, "output/sim/sim_paper_cumrrb_censored_strata.csv",
                 row.names = FALSE)
utils::write.csv(summ,
                 "output/sim/sim_paper_cumrrb_censored_strata_summary.csv",
                 row.names = FALSE)
