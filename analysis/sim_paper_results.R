# sim_paper_results.R ----------------------------------------------------------
# Replicate-level results of the manuscript scenarios, read from the runs
# that supply them (see paper_scenarios() in analysis/sim_scenarios.R).
# Sourced by 02b, 02d and 02e after
#   pkgload::load_all("."); source("analysis/sim_scenarios.R")

# Rows of the requested ids (default: all) from sim_results_<run>.csv in
# `sim_dir` (default: CCW_SIM_DIR, else output/sim), with id, part, label
# and run in front; columns absent from a run are NA.
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
