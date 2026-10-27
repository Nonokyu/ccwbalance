.check_methods <- function(method) {
  bad <- setdiff(method, .rr_methods)
  if (length(bad) > 0) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Use one of: ", paste(.rr_methods, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Risk-set reconstruction balance (RRB)
#'
#' For each strategy `g`, history feature `f_j(H_k)` and interval `k`,
#' compares the weighted mean of the post set with the mean of the reference
#' set, standardised by the reference-weighted standard deviation:
#' \deqn{RRB_{j,g}(k) = (\mu^{post}_{j,g}(k) - \mu^{ref}_{j,g}(k)) /
#'   s^{ref}_{j,g}(k).}
#'
#' ## Scope
#' * `scope = "interval"` (default): the reference set is the risk set of
#'   the stratum just before artificial censoring at the end of interval
#'   `k`, and the post set is the rows retained (`D_k = 0`). Under a
#'   correctly specified censoring model and positivity the target is 0 for
#'   every history function `f(H_k)`, including post-baseline covariates.
#' * `scope = "cumulative"`: the reference set is the unweighted risk set at
#'   the first interval of the strategy (the baseline population). The post
#'   set is the rows retained at `k`, weighted by the cumulative weight
#'   `W_k`, together with the individuals who had the outcome while retained
#'   in an earlier interval `tau < k` (`carry_events = TRUE`). These enter
#'   with their feature values at `tau` and the weight `W_tau Q_k / Q_tau`,
#'   that is, the weight at the event rescaled by the stratum numerators
#'   `q_l`, `tau < l <= k` (with unstabilised weights the factor is 1). The
#'   target is 0 only for baseline-fixed features (`X`, `<v>_base` from
#'   [make_history_features()]) and only when loss to follow-up is absent or
#'   weighted. For post-baseline features such as `L_k` the cumulative RRB
#'   is descriptive.
#'
#' ## Methods (reference / post weights)
#' * `"unweighted"`: 1 / 1;
#' * `"unstabilized"`: interval `W^u_{k-1}` / `W^u_{k-1} / pi_k`;
#'   cumulative 1 / `W^u_k`;
#' * `"stabilized"`: interval `W_{k-1} q_k` / `W_{k-1} q_k / pi_k`;
#'   cumulative 1 / `W_k`. The rescaling of carried rows needs `q_hat` to be
#'   constant within strata (the default numerator of
#'   [fit_censoring_model()]); a warning is given otherwise;
#' * `"weighted"`: the columns `w_ref` / `w_post` of [fit_censoring_model()]
#'   (cumulative scope: 1 / `w_post`). With truncated weights this gives the
#'   RRB of the weights used by the estimator. The interval RRB of the
#'   untruncated model weight (post weight `w_ref / pi_hat`) is returned by
#'   [rrb_null()] as `cells$rrb_obs`.
#'
#' @param data Long data with weights (see [fit_censoring_model()]).
#' @param features Character vector of feature columns `f(H_k)`.
#' @param method One or more of `"weighted"`, `"unweighted"`,
#'   `"unstabilized"`, `"stabilized"`.
#' @param scope One or more of `"interval"`, `"cumulative"`.
#' @param strategy,time,censored Column names.
#' @param w_ref,w_post,pi_hat,q_hat,w_prev_unstab,w_prev_stab Weight column
#'   names (see [fit_censoring_model()]).
#' @param event Outcome column used to carry individuals with the outcome
#'   forward (cumulative scope).
#' @param carry_events Carry individuals with the outcome forward
#'   (cumulative scope)?
#' @return A data frame with `strategy, time, feature, method, scope, n_risk,
#'   n_censored, mean_pre, mean_post, sd_pre, rrb` (`mean_pre` and `sd_pre`
#'   refer to the reference set). In the interval scope `n_risk` and
#'   `n_censored` are the stratum counts; in the cumulative scope `n_risk` is
#'   the baseline size and `n_censored` the number artificially censored up
#'   to and including interval `k`. `rrb` is `NA` when the reference
#'   standard deviation is 0.
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

#' Weighted Kolmogorov-Smirnov reconstruction distance (KS^RR)
#'
#' \deqn{KS^{RR}_{j,g}(k) = \sup_x |F^{post}_{j,g,k}(x) - F^{ref}_{j,g,k}(x)|,}
#' the largest distance between the weighted empirical distribution
#' functions of the post and reference sets of [rrb()] (same scopes, methods
#' and zero targets).
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

#' Mahalanobis reconstruction balance (M^RR)
#'
#' \deqn{M^{RR}_g(k) = \Delta_g(k)^\top \Sigma_{g,k}^{-1} \Delta_g(k),}
#' where `Delta_g(k)` is the vector of differences between the post and
#' reference means of the features and `Sigma_{g,k}` the reference-weighted
#' covariance matrix of the reference set (sets and methods as in [rrb()]).
#' Features with reference variance below `1e-10` in a stratum are dropped
#' there. With `ridge > 0`, `Sigma + ridge * mean(diag(Sigma)) * I` is
#' inverted. In the cumulative scope the target is 0 only when all features
#' are baseline-fixed.
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

#' Between-strategy balance (standardised mean differences)
#'
#' \deqn{SMD_{j}(g, k) = (\mu_{j,g}(k) - \mu_{j,g_0}(k)) /
#'   \sqrt{(s^2_{j,g}(k) + s^2_{j,g_0}(k)) / 2}}
#' compares each strategy `g` with a reference strategy `g_0` in every
#' interval. With `weighted = TRUE` the means and variances are over the
#' rows retained in interval `k` (`D_k = 0`), weighted by `w_post`; with
#' `weighted = FALSE` they are unweighted over the at-risk rows. The
#' end-of-grace SMD is the weighted value in the first interval after the
#' grace period. For post-baseline covariates affected by earlier treatment
#' this contrast is descriptive and has no zero target.
#'
#' @param data Long data (with weights if `weighted = TRUE`).
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
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' sb <- strategy_balance(w, "X", reference = "g1")
#' sb[sb$time == 3, ]
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
#' Within each stratum, compares the rows censored at the end of interval
#' `k` (`D_k = 1`) with the retained rows (`D_k = 0`):
#' \deqn{CB_{j,g}(k) = (\mu_{j,g,k}(D_k = 1) - \mu_{j,g,k}(D_k = 0)) /
#'   s_{j,g,k}.}
#' With `weighted = FALSE` the means and the standard deviation `s` of the
#' stratum are unweighted. With `weighted = TRUE` each row is weighted by
#' the inverse probability of its observed censoring status:
#' * `weights = "model"` (default): censored rows get
#'   `W^u_{k-1} / (1 - pi_k)` and retained rows `W^u_{k-1} / pi_k`, where
#'   `W^u_{k-1}` is the untruncated, unstabilised earlier weight
#'   (`w_prev_unstab`); `s` is `W^u_{k-1}`-weighted. This version checks the
#'   fitted censoring model and is not affected by truncation.
#' * `weights = "estimator"`: retained rows get `w_post` (the possibly
#'   truncated `W_k` of the estimator) and censored rows
#'   `min(w_ref / (1 - pi_k), cap)`; `s` is `w_ref`-weighted.
#'
#' Without truncation both versions give the same values. `1 - pi_k` is
#' floored at `1e-6`. Only strata with at least one censored and one
#' retained row are returned. Loss to follow-up can be examined with
#' `censored = "ltfu_next"` (unweighted).
#'
#' @param data Long data (with weights if `weighted = TRUE`).
#' @param features Feature columns.
#' @param weighted Logical.
#' @param weights `"model"` or `"estimator"` (used when `weighted = TRUE`).
#' @param strategy,time,censored,pi_hat,w_prev,w_ref,w_post Column names
#'   (`w_prev` is used by `weights = "model"`, `w_ref` and `w_post` by
#'   `weights = "estimator"`).
#' @param cap Truncation cap for the censored-row weights of
#'   `weights = "estimator"`; `NULL` (default) takes
#'   `attr(data, "weighting")$cap` (`Inf` if absent).
#' @return A data frame with `strategy, time, feature, weighted, weights,
#'   n_censored, n_retained, mean_censored, mean_retained, sd, smd`
#'   (`weights` is `NA` when unweighted).
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
                              w_post = "w_post", cap = NULL) {
  weights <- match.arg(weights)
  .check_cols(data, c(strategy, time, censored, features))
  cen <- data[[censored]] == 1
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
                     .cen = cen, .w = data$.w, .w0 = data$.w0,
                     data[, features, drop = FALSE], check.names = FALSE) |>
    tidyr::pivot_longer(dplyr::all_of(features), names_to = "feature",
                        values_to = ".value")
  long <- long[!is.na(long$.value), , drop = FALSE]
  out <- long |>
    dplyr::group_by(.data$strategy, .data$time, .data$feature) |>
    dplyr::filter(any(.data$.cen), any(!.data$.cen)) |>
    dplyr::summarise(
      n_censored = sum(.data$.cen),
      n_retained = sum(!.data$.cen),
      mean_censored = .wmean(.data$.value[.data$.cen], .data$.w[.data$.cen]),
      mean_retained = .wmean(.data$.value[!.data$.cen],
                             .data$.w[!.data$.cen]),
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

#' Aggregate balance over intervals: IAI and MaxIB
#'
#' \deqn{IAI = \sum_k v_k |B(k)| / \sum_k v_k, \qquad MaxIB = \max_k |B(k)|,}
#' with equal interval weights `v_k = 1` by default or a column of weights
#' (e.g. risk-set size). The summaries are taken over all rows of `x` with a
#' non-missing value. To obtain MaxIB over the strata with artificial
#' censoring, filter `x` first, e.g. `x[x$n_censored > 0, ]`.
#'
#' @param x Balance table, e.g. the output of [rrb()].
#' @param value Column holding the balance value (e.g. `"rrb"`, `"ks"`,
#'   `"smd"`, `"mahalanobis"`).
#' @param by Grouping columns (those absent from `x` are ignored).
#' @param time Interval column.
#' @param time_weight `NULL` (equal weights) or a column name.
#' @param times Optional subset of intervals.
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

#' Weight support and positivity summaries
#'
#' For each stratum: risk-set size, number retained, effective sample size
#' `ESS = (sum w)^2 / sum w^2` of the retained rows, `ESS / n_risk`,
#' normalised entropy `-sum p log p / log N` of the positive retained
#' weights (`p = w / sum w`), and the mean, 95th and 99th percentiles and
#' maximum of the retained weights. The `pi_hat` summaries (minimum, 1st
#' and 5th percentiles and the proportion below `pi_threshold`) are over all
#' at-risk rows, including rows where censoring is not possible
#' (`pi_hat = 1`). These summaries can signal, but not establish, a lack of
#' positivity.
#'
#' @param data Long data with weights (see [fit_censoring_model()]).
#' @param weight Weight column of the retained rows.
#' @param strategy,time,censored,pi_hat Column names.
#' @param pi_threshold Threshold for the proportion of small `pi_hat`.
#' @return A data frame with `strategy, time, n_risk, n_retained, ess,
#'   ess_ratio, entropy, pi_min, pi_p01, pi_p05, pr_pi_small, w_mean, w95,
#'   w99, w_max`.
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
