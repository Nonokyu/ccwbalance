# sim_paper_results.R ----------------------------------------------------------
# Replicate-level results of the scenarios of the manuscript. Each id of
# paper_scenarios() (analysis/sim_scenarios.R) points to the run scenario
# whose results the manuscript reports, so that the later runs replace the
# superseded rows of earlier runs by construction (see the column
# `supersedes` of paper_scenarios()).
# Sourced by the paper table scripts 02b, 02d and 02e after
#   pkgload::load_all("."); source("analysis/sim_scenarios.R")

# Read output/sim/sim_results_<run>.csv of the runs that supply the
# requested ids (default: all) and return the rows of their run scenarios,
# in id order and in file order within an id, with the columns id, part,
# label and run in front (dplyr::bind_rows, so that columns absent
# from a run, e.g. the per-feature columns of X2-X9, are NA). sim_dir
# defaults to the environment variable CCW_SIM_DIR, else output/sim.
read_paper_results <- function(ids = NULL, sim_dir = NULL) {
  if (is.null(sim_dir)) {
    sim_dir <- Sys.getenv("CCW_SIM_DIR", "")
    if (!nzchar(sim_dir)) sim_dir <- file.path("output", "sim")
  }
  ps <- paper_scenarios()
  if (!is.null(ids)) {
    stopifnot(all(ids %in% ps$id))
    ps <- ps[match(ids, ps$id), ]
  }
  runs <- unique(ps$run)
  raw <- stats::setNames(lapply(runs, function(r) {
    utils::read.csv(file.path(sim_dir, sprintf("sim_results_%s.csv", r)))
  }), runs)
  out <- lapply(seq_len(nrow(ps)), function(k) {
    x <- raw[[ps$run[k]]]
    x <- x[x$scenario == ps$scenario[k], , drop = FALSE]
    if (nrow(x) == 0) {
      stop("No results for ", ps$id[k], " (", ps$scenario[k],
           ") in sim_results_", ps$run[k], ".csv")
    }
    key <- data.frame(id = ps$id[k], part = ps$part[k], label = ps$label[k],
                      run = ps$run[k])
    cbind(key[rep(1, nrow(x)), ], x, row.names = NULL)
  })
  as.data.frame(dplyr::bind_rows(out))
}
