#' Individual-level bootstrap of RRB
#'
#' Nonparametric bootstrap of the RRB trajectory (cumulative scope by
#' default). Individuals are resampled with replacement, with all their rows
#' under every strategy; each resampled individual gets a new identifier,
#' so an individual drawn twice contributes two independent histories. The
#' censoring model is refitted on each resample with the settings stored by
#' [fit_censoring_model()] (formula, numerator, model groups, stabilisation,
#' truncation) and [rrb()] is recomputed.
#'
#' This is an uncertainty summary of the observed RRB, mainly for the
#' cumulative scope, where the RRB of time-varying features is descriptive.
#' For a zero-targeted reference of the interval scope see [rrb_null()].
#'
#' @param data Long data returned by [fit_censoring_model()] (its
#'   `"weighting"` attribute supplies the model settings).
#' @param features Character vector of feature columns `f(H_k)`.
#' @param B Number of bootstrap resamples.
#' @param scope `"cumulative"` (default) or `"interval"` (see [rrb()]).
#' @param method One RRB weighting method (see [rrb()]).
#' @param seed Optional random seed.
#' @param id,strategy,time Column names.
#' @param ... Further arguments passed to [rrb()] (e.g. `censored`,
#'   `event`, `carry_events`).
#' @return A data frame with one row per observed strategy x interval x
#'   feature cell: `strategy, time, feature, n_risk, rrb_obs, boot_mean,
#'   boot_se, lcl, ucl` (percentile 2.5% and 97.5% limits), `n_boot_ok`
#'   (resamples in which the cell's RRB is defined), `B` and `scope`. A cell
#'   absent from a resample is ignored in that resample.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' bt <- rrb_boot(w, c("X", "L"), B = 5, seed = 1)
#' head(bt)
rrb_boot <- function(data, features, B = 200,
                     scope = c("cumulative", "interval"),
                     method = "weighted", seed = NULL, id = "id",
                     strategy = "strategy", time = "time", ...) {
  scope <- match.arg(scope)
  stopifnot(length(method) == 1, B >= 1)
  wt <- attr(data, "weighting")
  if (is.null(wt) || is.null(wt$formula)) {
    stop("`data` has no \"weighting\" attribute: use data returned by ",
         "fit_censoring_model().", call. = FALSE)
  }
  .check_cols(data, c(id, strategy, time, features))
  data <- as.data.frame(data)

  obs <- rrb(data, features, method = method, scope = scope,
             strategy = strategy, time = time, ...)
  key_obs <- paste(obs$strategy, obs$time, obs$feature, sep = "\r")

  base <- data[, setdiff(names(data), .boot_weight_cols()), drop = FALSE]
  idx <- split(seq_len(nrow(base)), base[[id]])
  n_id <- length(idx)

  if (!is.null(seed)) set.seed(seed)
  sims <- matrix(NA_real_, nrow(obs), B)
  for (b in seq_len(B)) {
    s <- sample.int(n_id, n_id, replace = TRUE)
    rb <- .boot_rrb(base, idx[s], wt, features, method, scope, id, strategy,
                    time, ...)
    k <- match(key_obs, paste(rb$strategy, rb$time, rb$feature, sep = "\r"))
    sims[, b] <- rb$rrb[k]
  }

  n_ok <- rowSums(!is.na(sims))
  q <- function(p) {
    apply(sims, 1, function(v) {
      if (all(is.na(v))) NA_real_ else
        stats::quantile(v, p, na.rm = TRUE, names = FALSE)
    })
  }
  data.frame(
    strategy = obs$strategy, time = obs$time, feature = obs$feature,
    n_risk = obs$n_risk, rrb_obs = obs$rrb,
    boot_mean = ifelse(n_ok > 0, rowMeans(sims, na.rm = TRUE), NA_real_),
    boot_se = apply(sims, 1, function(v) {
      if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE)
    }),
    lcl = q(0.025), ucl = q(0.975), n_boot_ok = n_ok, B = B, scope = scope
  )
}

# Columns added by fit_censoring_model(), dropped before refitting
.boot_weight_cols <- function() {
  c("pi_hat", "q_hat", "weight_unstab", "weight_stab", "w_prev_unstab",
    "w_prev_stab", "weight", "w_prev", "w_post", "w_ref",
    "weights_stabilized", ".uncens")
}

# RRB of one resample; `rows` holds one row-index vector per resampled
# individual, which gets its position as the new identifier
.boot_rrb <- function(base, rows, wt, features, method, scope, id, strategy,
                      time, ...) {
  db <- base[unlist(rows, use.names = FALSE), , drop = FALSE]
  db[[id]] <- rep(seq_along(rows), lengths(rows))
  wb <- fit_censoring_model(db, wt$formula, numerator = wt$numerator,
                            by = wt$by, stabilized = wt$stabilized,
                            truncate = wt$truncate, id = id,
                            strategy = strategy, time = time)
  rrb(wb, features, method = method, scope = scope, strategy = strategy,
      time = time, ...)
}
