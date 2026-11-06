#' Build history features f(H_t) for balance diagnostics
#'
#' Adds transformed history features used as `f(H_t)` in the reconstruction
#' diagnostics: squares `L_t^2`, lag products `L_t L_{t-1}`, differences
#' `Delta L_t = L_t - L_{t-1}`, cumulative means, threshold indicators
#' `I(L_t > c)`, baseline values `L_0` (for cumulative RRB, which is
#' 0-targeted only for baseline-fixed features) and pairwise interactions.
#'
#' Lags are computed within clone (`id` x `strategy`) after ordering by time.
#' If a `<var>_lag` column already exists (as produced by [simulate_ccw()]) it
#' is used as `L_{t-1}`; otherwise the lag is computed and, at the first row,
#' set equal to the current value (so `Delta L_0 = 0`).
#'
#' @param data Long clone data.
#' @param vars Time-varying covariates to transform.
#' @param thresholds Numeric vector of cut-offs `c` for `I(L_t > c)`.
#' @param interactions List of character pairs, e.g. `list(c("L", "X"))`.
#' @param baseline Add `<v>_base`, the value of `<v>` at the clone's first
#'   time (baseline-fixed feature)?
#' @param id,strategy,time Column names.
#' @return `data` with added columns `<v>_sq`, `<v>_lagprod`, `<v>_diff`,
#'   `<v>_cummean`, `<v>_gt<c>`, `<v>_base` and `<a>_x_<b>`.
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
