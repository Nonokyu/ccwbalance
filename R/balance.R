.check_methods <- function(method) {
  bad <- setdiff(method, .rr_methods)
  if (length(bad) > 0) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Use one of: ", paste(.rr_methods, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Risk-set Reconstruction Balance (RRB)
#'
#' For each strategy `g`, feature `j` and time `t`, compares a reference
#' ("pre") distribution with the IPCW-weighted distribution of the clones
#' retained after artificial censoring ("post"), standardised by the
#' (reference-weighted) reference SD:
#' \deqn{RRB_{j,g}(t) = \frac{\mu^{post,w}_{j,g}(t) - \mu^{pre}_{j,g}(t)}{s^{pre}_{j,g}(t)}.}
#'
#' ## Scope
#' * `scope = "interval"` (default; local reconstruction): the reference is
#'   the risk set at `t` just before artificial censoring in `[t, t+1)`, the
#'   post set the clones retained after it. Under a correctly specified
#'   censoring model and positivity the IPW identity
#'   `E[I(D = 0) f(H_t) / pi(H_t)] = E[f(H_t)]` gives a target of 0 for
#'   **every** history feature `f(H_t)`, including post-baseline `L_t`.
#' * `scope = "cumulative"` (reconstruction of the baseline population): the
#'   reference is the cloned population at baseline (first time of the
#'   strategy, unweighted), the post set the clones retained at `t` weighted
#'   by the full cumulative weight `W_t`, plus clones that had the outcome
#'   while uncensored before `t` (carried forward with their weight and
#'   feature values at the event; `carry_events = TRUE`). The target is 0
#'   **only for baseline-fixed features** (`X`, `L_0` via `<v>_base` from
#'   [make_history_features()]) and only if there is no loss to follow-up that
#'   is left unweighted. For post-baseline features such as `L_t`,
#'   cumulative RRB is **not** 0-targeted: `L_t` evolves over time and is
#'   affected by treatment (`A_{t-1} -> L_t`), so a non-zero value is expected
#'   even under correct weighting; it is descriptive only.
#'
#' ## Methods (weights on the reference / post side)
#' * `"unweighted"`: 1 / 1 (the raw selection induced by artificial
#'   censoring);
#' * `"unstabilized"`: interval scope `W_{t-1}` / `W_{t-1} / pi_t`; cumulative
#'   scope 1 / `W_t` (unstabilised);
#' * `"stabilized"`: interval scope `W_{t-1} q_t` / `W_{t-1} q_t / pi_t`;
#'   cumulative scope 1 / `W_t^s` (carried event clones
#'   are extended with the cell numerators `q_k`, which requires the default
#'   numerator of [fit_censoring_model()]; a warning is given otherwise);
#' * `"weighted"`: columns `w_ref` / `w_post` as produced by
#'   [fit_censoring_model()] (the chosen scheme, including truncation;
#'   cumulative scope uses `w_post` only).
#'
#' @param data Long clone data with weights (see [fit_censoring_model()]).
#' @param features Character vector of feature columns `f(H_t)`.
#' @param method One or more of `"weighted"`, `"unweighted"`,
#'   `"unstabilized"`, `"stabilized"`.
#' @param scope One or more of `"interval"`, `"cumulative"`.
#' @param strategy,time,censored Column names.
#' @param w_ref,w_post,pi_hat,q_hat,w_prev_unstab,w_prev_stab Weight column
#'   names.
#' @param event Outcome column used to carry event clones forward in the
#'   cumulative scope.
#' @param carry_events Carry forward clones with the outcome (cumulative
#'   scope)?
#' @return A data frame with `strategy, time, feature, method, scope, n_risk,
#'   n_censored, mean_pre, mean_post, sd_pre, rrb`. For the cumulative scope
#'   `n_risk` is the baseline number of clones and `n_censored` the
#'   cumulative number artificially censored up to and including `t`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' head(rrb(w, c("X", "L"), method = c("unweighted", "weighted")))
#' tail(rrb(w, "X", scope = "cumulative"))
rrb <- function(data, features, method = "weighted", scope = "interval",
                strategy = "strategy", time = "time",
                censored = "censored_next", w_ref = "w_ref",
                w_post = "w_post", pi_hat = "pi_hat", q_hat = "q_hat",
                w_prev_unstab = "w_prev_unstab", w_prev_stab = "w_prev_stab",
                event = "event_next", carry_events = TRUE) {
  .check_methods(method)
  .check_scope(scope)
  .check_cols(data, c(strategy, time, censored, features))
  grid <- expand.grid(method = method, scope = scope,
                      stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(grid)), function(k) {
    m <- grid$method[k]
    sc <- grid$scope[k]
    fr <- .rr_frame(data, features, m, sc, strategy, time, censored, w_ref,
                    w_post, pi_hat, q_hat, w_prev_unstab, w_prev_stab, event,
                    carry_events)
    out <- .to_long(fr$frame, features) |>
      dplyr::group_by(.data$strategy, .data$time, .data$feature) |>
      dplyr::summarise(
        mean_pre = .wmean(.data$.value, .data$.ref_w),
        mean_post = .wmean(.data$.value, .data$.post_w),
        sd_pre = sqrt(.wvar(.data$.value, .data$.ref_w)),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        method = m, scope = sc,
        rrb = ifelse(!is.na(.data$sd_pre) & .data$sd_pre > 1e-10,
                     (.data$mean_post - .data$mean_pre) / .data$sd_pre,
                     NA_real_)
      )
    .join_counts(out, fr$counts)
  })
  out <- dplyr::bind_rows(res)
  out <- out[, c("strategy", "time", "feature", "method", "scope", "n_risk",
                 "n_censored", "mean_pre", "mean_post", "sd_pre", "rrb")]
  as.data.frame(out)
}

.check_scope <- function(scope) {
  bad <- setdiff(scope, c("interval", "cumulative"))
  if (length(bad) > 0) {
    stop("Unknown scope(s): ", paste(bad, collapse = ", "),
         ". Use \"interval\" and/or \"cumulative\".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Weighted-ECDF Kolmogorov-Smirnov reconstruction distance (KS^RR)
#'
#' \deqn{KS^{RR}_{j,g}(t) = \sup_x |F^{post,w}_{j,g,t}(x) - F^{pre}_{j,g,t}(x)|}
#' using the same reference / post sets and weights as [rrb()] (see there for
#' the `"interval"` and `"cumulative"` scopes and which features are
#' 0-targeted in each).
#'
#' @inheritParams rrb
#' @return A data frame with `strategy, time, feature, method, scope, n_risk,
#'   n_censored, ks`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' head(ks_rr(w, "L"))
ks_rr <- function(data, features, method = "weighted", scope = "interval",
                  strategy = "strategy", time = "time",
                  censored = "censored_next", w_ref = "w_ref",
                  w_post = "w_post", pi_hat = "pi_hat", q_hat = "q_hat",
                  w_prev_unstab = "w_prev_unstab",
                  w_prev_stab = "w_prev_stab", event = "event_next",
                  carry_events = TRUE) {
  .check_methods(method)
  .check_scope(scope)
  .check_cols(data, c(strategy, time, censored, features))
  grid <- expand.grid(method = method, scope = scope,
                      stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(grid)), function(k) {
    m <- grid$method[k]
    sc <- grid$scope[k]
    fr <- .rr_frame(data, features, m, sc, strategy, time, censored, w_ref,
                    w_post, pi_hat, q_hat, w_prev_unstab, w_prev_stab, event,
                    carry_events)
    out <- .to_long(fr$frame, features) |>
      dplyr::group_by(.data$strategy, .data$time, .data$feature) |>
      dplyr::summarise(
        ks = {
          kp <- .data$.post_w > 0
          kr <- .data$.ref_w > 0
          .wks(.data$.value[kp], .data$.post_w[kp],
               .data$.value[kr], .data$.ref_w[kr])
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(method = m, scope = sc)
    .join_counts(out, fr$counts)
  })
  out <- dplyr::bind_rows(res)
  as.data.frame(out[, c("strategy", "time", "feature", "method", "scope",
                        "n_risk", "n_censored", "ks")])
}

#' Mahalanobis global reconstruction balance (M^RR)
#'
#' \deqn{M^{RR}_g(t) = \Delta_g(t)^\top \Sigma^{-1}_{g,t,ref} \Delta_g(t),
#'   \quad \Delta_g(t) = \mu^{post,w}_g(t) - \mu^{pre}_g(t),}
#' where `Sigma` is the reference-weighted covariance of the reference set
#' (see [rrb()] for the `"interval"` / `"cumulative"` scopes). Features with
#' (near) zero reference variance at a given time are dropped. Optional ridge
#' regularisation uses `Sigma + ridge * mean(diag(Sigma)) * I`. In the
#' cumulative scope the metric is 0-targeted only when all features are
#' baseline-fixed.
#'
#' @inheritParams rrb
#' @param ridge Non-negative ridge constant (0 = none).
#' @return A data frame with `strategy, time, method, scope, n_risk,
#'   n_censored, n_features, mahalanobis`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' head(mahalanobis_rr(w, c("X", "L", "L_lag")))
mahalanobis_rr <- function(data, features, method = "weighted", ridge = 0,
                           scope = "interval", strategy = "strategy",
                           time = "time", censored = "censored_next",
                           w_ref = "w_ref", w_post = "w_post",
                           pi_hat = "pi_hat", q_hat = "q_hat",
                           w_prev_unstab = "w_prev_unstab",
                           w_prev_stab = "w_prev_stab", event = "event_next",
                           carry_events = TRUE) {
  .check_methods(method)
  .check_scope(scope)
  .check_cols(data, c(strategy, time, censored, features))
  stopifnot(ridge >= 0)
  grid <- expand.grid(method = method, scope = scope,
                      stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(grid)), function(k) {
    m <- grid$method[k]
    sc <- grid$scope[k]
    fr <- .rr_frame(data, features, m, sc, strategy, time, censored, w_ref,
                    w_post, pi_hat, q_hat, w_prev_unstab, w_prev_stab, event,
                    carry_events)
    d <- fr$frame
    d <- d[stats::complete.cases(d[, features, drop = FALSE]), , drop = FALSE]
    keys <- unique(d[, c("strategy", "time")])
    keys <- keys[order(keys$strategy, keys$time), , drop = FALSE]
    rows <- lapply(seq_len(nrow(keys)), function(i) {
      sel <- d$strategy == keys$strategy[i] & d$time == keys$time[i]
      sub <- d[sel, , drop = FALSE]
      Fm <- as.matrix(sub[, features, drop = FALSE])
      wr <- sub$.ref_w
      wp <- sub$.post_w
      val <- NA_real_
      nf <- 0L
      if (sum(wr) > 0 && sum(wp) > 0) {
        mu_ref <- colSums(Fm * wr) / sum(wr)
        mu_post <- colSums(Fm * wp) / sum(wp)
        cen <- sweep(Fm, 2, mu_ref)
        S <- crossprod(cen * sqrt(wr)) / sum(wr)
        ok <- diag(S) > 1e-10
        nf <- sum(ok)
        if (nf > 0) {
          S <- S[ok, ok, drop = FALSE]
          if (ridge > 0) S <- S + ridge * mean(diag(S)) * diag(nf)
          delta <- (mu_post - mu_ref)[ok]
          val <- tryCatch(drop(t(delta) %*% solve(S, delta)),
                          error = function(e) NA_real_)
        }
      }
      data.frame(strategy = keys$strategy[i], time = keys$time[i],
                 method = m, scope = sc, n_features = nf, mahalanobis = val)
    })
    .join_counts(do.call(rbind, rows), fr$counts)
  })
  out <- do.call(rbind, res)
  out <- out[, c("strategy", "time", "method", "scope", "n_risk",
                 "n_censored", "n_features", "mahalanobis")]
  rownames(out) <- NULL
  out
}

#' Strategy-balance (descriptive dynamic strategy-risk-set contrast)
#'
#' \deqn{B^S_j(t) = \frac{\mu^w_{j,g}(t) - \mu^w_{j,g'}(t)}{s_{j,ref}}}
#' comparing each strategy with a reference strategy at each time, where
#' `s_ref` is the square root of the average of the two (weighted) variances.
#' With `weighted = TRUE` the means are over clones retained after artificial
#' censoring at `t`, weighted by `w_post`; with `weighted = FALSE` they are
#' unweighted means over the at-risk rows. For post-baseline covariates
#' affected by earlier treatment this is **descriptive only** and is not a
#' zero-balance target.
#'
#' @param data Long clone data.
#' @param features Feature columns.
#' @param reference Reference strategy (default: first level).
#' @param weighted Logical.
#' @param strategy,time,censored,w_post Column names.
#' @return A data frame with `strategy, reference, time, feature, weighted,
#'   mean, mean_ref, smd`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' head(strategy_balance(d, c("X", "L"), weighted = FALSE))
strategy_balance <- function(data, features, reference = NULL,
                             weighted = TRUE, strategy = "strategy",
                             time = "time", censored = "censored_next",
                             w_post = "w_post") {
  .check_cols(data, c(strategy, time, features))
  if (weighted) {
    .check_cols(data, c(censored, w_post))
    d <- data[data[[censored]] == 0, , drop = FALSE]
    d$.w <- d[[w_post]]
  } else {
    d <- data
    d$.w <- 1
  }
  lv <- if (is.factor(data[[strategy]])) levels(data[[strategy]]) else
    sort(unique(data[[strategy]]))
  if (is.null(reference)) reference <- lv[1]
  long <- data.frame(strategy = as.character(d[[strategy]]),
                     time = d[[time]], .w = d$.w,
                     d[, features, drop = FALSE], check.names = FALSE) |>
    tidyr::pivot_longer(dplyr::all_of(features), names_to = "feature",
                        values_to = ".value")
  long <- long[!is.na(long$.value), , drop = FALSE]
  summ <- long |>
    dplyr::group_by(.data$strategy, .data$time, .data$feature) |>
    dplyr::summarise(mean = .wmean(.data$.value, .data$.w),
                     var = .wvar(.data$.value, .data$.w), .groups = "drop")
  ref <- summ[summ$strategy == reference, c("time", "feature", "mean", "var")]
  names(ref) <- c("time", "feature", "mean_ref", "var_ref")
  out <- dplyr::inner_join(summ[summ$strategy != reference, ], ref,
                           by = c("time", "feature")) |>
    dplyr::mutate(
      reference = reference,
      weighted = weighted,
      smd = (.data$mean - .data$mean_ref) /
        sqrt((.data$var + .data$var_ref) / 2),
      smd = ifelse(is.finite(.data$smd), .data$smd, NA_real_),
      feature = factor(.data$feature, levels = features),
      strategy = factor(.data$strategy, levels = lv)
    ) |>
    dplyr::arrange(.data$strategy, .data$feature, .data$time)
  as.data.frame(out[, c("strategy", "reference", "time", "feature",
                        "weighted", "mean", "mean_ref", "smd")])
}

#' Censoring balance (after Jackson)
#'
#' Within each strategy x time, compares clones censored in the next interval
#' (`D_{t+1} = 1`) with those retained (`D_{t+1} = 0`):
#' \deqn{B^C_{j,g}(t) = \frac{\mu_{j,g,t}(D_{t+1}=1) - \mu_{j,g,t}(D_{t+1}=0)}{s_{j,g,t}}.}
#' With `weighted = TRUE` each clone is weighted by the inverse probability of
#' its observed censoring status, so that balance is expected under a correct
#' model; `s` is the reference-weighted SD of the at-risk set. Two weight
#' versions are available (`weights`):
#' * `"model"` (default): censored rows get `W_{t-1} / (1 - pi_t)` and
#'   retained rows `W_{t-1} / pi_t`, with the untruncated unstabilised
#'   `W_{t-1}` (`w_prev`); `s` is `W_{t-1}`-weighted. This checks the fitted
#'   censoring model only and cannot respond to weight truncation.
#' * `"estimator"`: retained rows get `w_post` (the possibly truncated
#'   cumulative weight `W_t` actually used by the estimator, i.e. the same
#'   post weights as `rrb(method = "weighted")`), censored rows
#'   `min(w_ref / (1 - pi_t), cap)`, where `w_ref` is the truncated
#'   `W_{t-1}` (times `q_t` when stabilised) and `cap` the truncation cap of
#'   [fit_censoring_model()]; `s` is `w_ref`-weighted (the reference weights
#'   of [rrb()]).
#'
#' Without truncation the two versions coincide: with unstabilised weights
#' the weights are identical, and with stabilised weights (default numerator,
#' common to all clones of a strategy x time cell) they differ only by the
#' cell-constant factor `prod_{k <= t} q_k`, which cancels in the weighted
#' means and SD.
#'
#' With `rows = "censorable"` only rows at which artificial censoring is
#' possible (`censor_possible == TRUE`) enter the censored and the retained
#' group; rows that cannot be censored in the next interval (e.g. clones of
#' individuals who already initiated, with `pi_t = 1`) are dropped from the
#' retained group. The standardising SD `s` is still computed over the whole
#' at-risk set (with the weights described above), so that `"all"` and
#' `"censorable"` share the same scale.
#'
#' Only strategy x time cells with at least one censored and one retained
#' clone (among the rows selected by `rows`) are returned. The same function
#' can be applied to loss to follow-up via `censored = "ltfu_next"` (unweighted).
#'
#' @param data Long clone data.
#' @param features Feature columns.
#' @param weighted Logical.
#' @param weights `"model"` or `"estimator"` (used when `weighted = TRUE`).
#' @param strategy,time,censored,pi_hat,w_prev,w_ref,w_post Column names
#'   (`w_prev` is used by `weights = "model"`, `w_ref` and `w_post` by
#'   `weights = "estimator"`).
#' @param rows `"all"` (default; every row of the risk set) or
#'   `"censorable"` (only rows with `censor_possible == TRUE` enter the two
#'   groups; see Details).
#' @param censor_possible Column flagging rows at which artificial
#'   censoring is possible (used by `rows = "censorable"`).
#' @param cap Truncation cap for the censored-clone weights of
#'   `weights = "estimator"`; `NULL` (default) takes
#'   `attr(data, "weighting")$cap` (`Inf` if absent).
#' @return A data frame with `strategy, time, feature, weighted, weights,
#'   n_censored, n_retained, mean_censored, mean_retained, sd, smd`
#'   (`weights` is `"model"` or `"estimator"`, `NA` when unweighted).
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' head(censoring_balance(d, c("X", "L")))
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X, truncate = 0.9)
#' head(censoring_balance(w, c("X", "L"), weighted = TRUE,
#'                        weights = "estimator"))
censoring_balance <- function(data, features, weighted = FALSE,
                              weights = c("model", "estimator"),
                              strategy = "strategy", time = "time",
                              censored = "censored_next", pi_hat = "pi_hat",
                              w_prev = "w_prev_unstab", w_ref = "w_ref",
                              w_post = "w_post", cap = NULL,
                              rows = c("all", "censorable"),
                              censor_possible = "censor_possible") {
  weights <- match.arg(weights)
  rows <- match.arg(rows)
  .check_cols(data, c(strategy, time, censored, features,
                      if (rows == "censorable") censor_possible))
  cen <- data[[censored]] == 1
  poss <- if (rows == "censorable") {
    as.logical(data[[censor_possible]]) %in% TRUE
  } else {
    rep(TRUE, nrow(data))
  }
  if (weighted && weights == "model") {
    .check_cols(data, c(pi_hat, w_prev))
    p <- data[[pi_hat]]
    wprev <- data[[w_prev]]
    data$.w <- ifelse(cen, wprev / pmax(1 - p, 1e-6), wprev / p)
    data$.w0 <- wprev
  } else if (weighted) {
    .check_cols(data, c(pi_hat, w_ref, w_post))
    if (is.null(cap)) {
      cap <- attr(data, "weighting")$cap
      if (is.null(cap)) cap <- Inf
    }
    p <- data[[pi_hat]]
    wr <- data[[w_ref]]
    data$.w <- ifelse(cen, pmin(wr / pmax(1 - p, 1e-6), cap), data[[w_post]])
    data$.w0 <- wr
  } else {
    data$.w <- 1
    data$.w0 <- 1
  }
  wlab <- if (weighted) weights else NA_character_
  long <- data.frame(strategy = data[[strategy]], time = data[[time]],
                     .cen = cen, .poss = poss, .w = data$.w,
                     .w0 = data$.w0, data[, features, drop = FALSE],
                     check.names = FALSE) |>
    tidyr::pivot_longer(dplyr::all_of(features), names_to = "feature",
                        values_to = ".value")
  long <- long[!is.na(long$.value), , drop = FALSE]
  out <- long |>
    dplyr::group_by(.data$strategy, .data$time, .data$feature) |>
    dplyr::filter(any(.data$.cen & .data$.poss),
                  any(!.data$.cen & .data$.poss)) |>
    dplyr::summarise(
      n_censored = sum(.data$.cen & .data$.poss),
      n_retained = sum(!.data$.cen & .data$.poss),
      mean_censored = .wmean(.data$.value[.data$.cen & .data$.poss],
                             .data$.w[.data$.cen & .data$.poss]),
      mean_retained = .wmean(.data$.value[!.data$.cen & .data$.poss],
                             .data$.w[!.data$.cen & .data$.poss]),
      sd = sqrt(.wvar(.data$.value, .data$.w0)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      weighted = weighted,
      weights = wlab,
      smd = ifelse(.data$sd > 1e-10,
                   (.data$mean_censored - .data$mean_retained) / .data$sd,
                   NA_real_),
      feature = factor(.data$feature, levels = features)
    ) |>
    dplyr::arrange(.data$strategy, .data$feature, .data$time)
  as.data.frame(out[, c("strategy", "time", "feature", "weighted", "weights",
                        "n_censored", "n_retained", "mean_censored",
                        "mean_retained", "sd", "smd")])
}

#' Time aggregation of balance trajectories: IAI and MaxIB
#'
#' \deqn{IAI = \frac{\sum_t q(t) |B(t)|}{\sum_t q(t)}, \qquad
#'   MaxIB = \max_t |B(t)|}
#' with equal time weights `q(t) = 1` by default (recommended primary
#' analysis); a column of time weights (e.g. risk-set size or ESS) can be
#' supplied for secondary analyses.
#'
#' @param x Balance trajectory table (e.g. output of [rrb()]).
#' @param value Column holding the balance value (e.g. `"rrb"`, `"ks"`,
#'   `"smd"`, `"mahalanobis"`).
#' @param by Grouping columns (those absent from `x` are ignored).
#' @param time Time column.
#' @param time_weight `NULL` (equal weights) or a column name.
#' @param times Optional subset of times to aggregate over.
#' @return A data frame with the grouping columns and `n_times, iai, maxib,
#'   time_at_max`.
#' @export
#' @examples
#' x <- data.frame(feature = "L", time = 0:3, rrb = c(0, 0.1, -0.3, 0.2))
#' aggregate_balance(x)
aggregate_balance <- function(x, value = "rrb",
                              by = c("strategy", "feature", "method",
                                     "scope"),
                              time = "time", time_weight = NULL,
                              times = NULL) {
  .check_cols(x, c(value, time))
  by <- intersect(by, names(x))
  if (!is.null(times)) x <- x[x[[time]] %in% times, , drop = FALSE]
  x$.v <- abs(x[[value]])
  x$.q <- if (is.null(time_weight)) 1 else x[[time_weight]]
  x$.t <- x[[time]]
  x <- x[!is.na(x$.v), , drop = FALSE]
  out <- x |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_times = dplyr::n(),
      iai = sum(.data$.q * .data$.v) / sum(.data$.q),
      maxib = max(.data$.v),
      time_at_max = .data$.t[which.max(.data$.v)],
      .groups = "drop"
    )
  as.data.frame(out)
}

#' Weight-support and positivity diagnostics by strategy x time
#'
#' For each strategy and time: risk-set size, number retained, effective
#' sample size `ESS = (sum w)^2 / sum w^2` of the retained clones, ESS ratio
#' `ESS / N_risk`, normalised weight entropy
#' `-sum p log p / log N` (`p = w / sum w`), `pi_hat` summaries over the
#' at-risk rows (min, 1st and 5th percentiles, `Pr(pi_hat < 0.05)`), and
#' weight quantiles (mean, 95th, 99th, max). These are weight-concentration /
#' support metrics, not balance metrics; they can signal but not prove
#' positivity failure.
#'
#' @param data Weighted long clone data.
#' @param weight Weight column of retained clones.
#' @param strategy,time,censored,pi_hat Column names.
#' @param pi_threshold Threshold for `Pr(pi_hat < threshold)`.
#' @return A data frame of diagnostics.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' head(weight_support(w))
weight_support <- function(data, weight = "w_post", strategy = "strategy",
                           time = "time", censored = "censored_next",
                           pi_hat = "pi_hat", pi_threshold = 0.05) {
  .check_cols(data, c(strategy, time, censored, weight, pi_hat))
  d <- data.frame(strategy = data[[strategy]], time = data[[time]],
                  .cen = data[[censored]] == 1, .w = data[[weight]],
                  .pi = data[[pi_hat]])
  out <- d |>
    dplyr::group_by(.data$strategy, .data$time) |>
    dplyr::summarise(
      n_risk = dplyr::n(),
      n_retained = sum(!.data$.cen),
      ess = {
        w <- .data$.w[!.data$.cen]
        if (length(w) == 0) NA_real_ else sum(w)^2 / sum(w^2)
      },
      entropy = {
        w <- .data$.w[!.data$.cen]
        w <- w[w > 0]
        if (length(w) < 2) NA_real_ else {
          p <- w / sum(w)
          -sum(p * log(p)) / log(length(w))
        }
      },
      pi_min = min(.data$.pi),
      pi_p01 = stats::quantile(.data$.pi, 0.01, names = FALSE),
      pi_p05 = stats::quantile(.data$.pi, 0.05, names = FALSE),
      pr_pi_small = mean(.data$.pi < pi_threshold),
      w_mean = .safe(.data$.w[!.data$.cen], mean),
      w95 = .safe(.data$.w[!.data$.cen], stats::quantile, 0.95,
                  names = FALSE),
      w99 = .safe(.data$.w[!.data$.cen], stats::quantile, 0.99,
                  names = FALSE),
      w_max = .safe(.data$.w[!.data$.cen], max),
      .groups = "drop"
    ) |>
    dplyr::mutate(ess_ratio = .data$ess / .data$n_risk, .after = "ess")
  as.data.frame(out)
}
