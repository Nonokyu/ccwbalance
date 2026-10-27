#' Build history features f(H_k)
#'
#' Adds transformed history features `f(H_k)` for the diagnostics: squares
#' `L_k^2`, lag products `L_k L_{k-1}`, differences `L_k - L_{k-1}`,
#' cumulative means, threshold indicators `I(L_k > c)`, baseline values
#' `L_0` and pairwise interactions. Baseline values are baseline-fixed, so
#' their cumulative RRB is zero-targeted (see [rrb()]).
#'
#' Lags are taken within individual x strategy after ordering by `time`. An
#' existing column `<v>_lag` (as returned by [simulate_ccw()]) is used as
#' `L_{k-1}`; otherwise the lag is computed and set to the current value in
#' the first row, so that `L_0 - L_{-1} = 0`.
#'
#' @param data Long data (see [ccwbalance-package]).
#' @param vars Time-varying covariates to transform.
#' @param thresholds Numeric vector of cut-offs `c` for `I(L_k > c)`.
#' @param interactions List of character pairs, e.g. `list(c("L", "X"))`.
#' @param baseline Add `<v>_base`, the value of `<v>` in the first row of the
#'   individual under the strategy?
#' @param id,strategy,time Column names.
#' @return `data`, sorted by strategy, individual and interval, with added
#'   columns `<v>_sq`, `<v>_lagprod`, `<v>_diff`, `<v>_cummean`, `<v>_gt<c>`,
#'   `<v>_base` and `<a>_x_<b>`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 100), seed = 1)
#' head(make_history_features(d, interactions = list(c("L", "X"))))
make_history_features <- function(data, vars = "L", thresholds = 0,
                                   interactions = list(), baseline = TRUE,
                                   id = "id",
                                   strategy = "strategy", time = "time") {
  .check_cols(data, c(id, strategy, time, vars))
  ord <- order(data[[strategy]], data[[id]], data[[time]])
  data <- data[ord, , drop = FALSE]
  grp <- interaction(data[[strategy]], data[[id]], drop = TRUE)
  for (v in vars) {
    x <- data[[v]]
    lag_col <- paste0(v, "_lag")
    if (lag_col %in% names(data)) {
      xl <- data[[lag_col]]
    } else {
      xl <- stats::ave(x, grp, FUN = function(z) c(z[1], z[-length(z)]))
    }
    data[[paste0(v, "_sq")]] <- x^2
    data[[paste0(v, "_lagprod")]] <- x * xl
    data[[paste0(v, "_diff")]] <- x - xl
    data[[paste0(v, "_cummean")]] <- stats::ave(
      x, grp, FUN = function(z) cumsum(z) / seq_along(z)
    )
    if (baseline) {
      data[[paste0(v, "_base")]] <- stats::ave(x, grp,
                                               FUN = function(z) rep(z[1],
                                                                     length(z)))
    }
    for (cc in thresholds) {
      data[[paste0(v, "_gt", format(cc))]] <- as.numeric(x > cc)
    }
  }
  for (pr in interactions) {
    .check_cols(data, pr)
    data[[paste0(pr[1], "_x_", pr[2])]] <- data[[pr[1]]] * data[[pr[2]]]
  }
  rownames(data) <- NULL
  data
}
