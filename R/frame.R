# Reference and post sets of the reconstruction diagnostics.
#
# Returns a list with `frame` (strategy, time, .ref_w, .post_w and the
# features; a row may belong to the reference set, the post set or both)
# and `counts` (strategy, time, n_risk, n_censored).
#
# Interval scope: reference = risk set of the stratum, post = rows retained
# after artificial censoring (see .rr_weights()).
# Cumulative scope: reference = unweighted risk set at the first interval of
# the strategy, post = rows retained at k weighted by W_k, plus rows with
# the outcome while retained in an earlier interval tau (feature values of
# that row; weight W_tau, times Q_k / Q_tau for stabilised weights).

.rr_frame <- function(data, features, method, scope, strategy, time,
                      censored, w_ref, w_post, pi_hat, q_hat, w_prev_unstab,
                      w_prev_stab, event = "event_next",
                      carry_events = TRUE) {
  if (scope == "interval") {
    d <- .rr_weights(data, method, censored, w_ref, w_post, pi_hat, q_hat,
                     w_prev_unstab, w_prev_stab)
    frame <- data.frame(strategy = d[[strategy]], time = d[[time]],
                        .ref_w = d$.ref_w, .post_w = d$.post_w,
                        d[, features, drop = FALSE], check.names = FALSE)
    counts <- stats::aggregate(
      list(n_risk = rep(1L, nrow(d)), n_censored = d[[censored]] == 1),
      by = list(strategy = d[[strategy]], time = d[[time]]), FUN = sum)
    return(list(frame = frame, counts = counts))
  }
  if (scope != "cumulative") {
    stop("`scope` must be \"interval\" or \"cumulative\".", call. = FALSE)
  }
  .cum_frame(data, features, method, strategy, time, censored, w_post,
             pi_hat, q_hat, w_prev_unstab, w_prev_stab, event, carry_events)
}

.cum_frame <- function(data, features, method, strategy, time, censored,
                       w_post, pi_hat, q_hat, w_prev_unstab, w_prev_stab,
                       event, carry_events) {
  .check_cols(data, c(strategy, time, censored))
  n <- nrow(data)
  ret <- data[[censored]] == 0
  # cumulative weight W_k of each row and whether it is stabilised
  stab <- FALSE
  if (method == "unweighted") {
    W <- rep(1, n)
  } else if (method == "unstabilized") {
    .check_cols(data, c(pi_hat, w_prev_unstab))
    W <- data[[w_prev_unstab]] / data[[pi_hat]]
  } else if (method == "stabilized") {
    .check_cols(data, c(pi_hat, q_hat, w_prev_stab))
    W <- data[[w_prev_stab]] * data[[q_hat]] / data[[pi_hat]]
    stab <- TRUE
  } else if (method == "weighted") {
    .check_cols(data, w_post)
    W <- data[[w_post]]
    stab <- "weights_stabilized" %in% names(data) &&
      isTRUE(any(as.logical(data$weights_stabilized)))
  } else {
    stop("Unknown method: ", method, call. = FALSE)
  }
  sv <- data[[strategy]]
  tv <- data[[time]]
  # stratum numerators q_k, used to rescale carried stabilised weights
  if (stab) {
    .check_cols(data, q_hat)
    qrange <- stats::ave(data[[q_hat]], sv, tv,
                         FUN = function(z) {
                           if (length(z) > 0) diff(range(z)) else 0
                         })
    if (any(qrange > 1e-8)) {
      warning("q_hat varies within strategy x time cells: the stabilised ",
              "pseudo-population is then not the baseline population and ",
              "cumulative RRB is not 0-targeted. Use the default numerator ",
              "or method = \"unstabilized\".", call. = FALSE)
    }
  }
  has_event <- carry_events && event %in% names(data)
  ev <- if (has_event) {
    ret & !is.na(data[[event]]) & data[[event]] == 1
  } else {
    rep(FALSE, n)
  }
  strategies <- if (is.factor(sv)) levels(droplevels(sv)) else
    sort(unique(sv))
  frames <- list()
  counts <- list()
  for (g in strategies) {
    ig <- which(as.character(sv) == as.character(g))
    times <- sort(unique(tv[ig]))
    t0 <- times[1]
    base <- ig[tv[ig] == t0]
    qbar <- NULL
    if (stab) {
      qbar <- tapply(data[[q_hat]][ig], tv[ig], mean)
    }
    cens_by_t <- tapply(data[[censored]][ig] == 1, tv[ig], sum)
    for (tt in times) {
      post_rows <- ig[tv[ig] == tt & ret[ig]]
      post_w <- W[post_rows]
      carry <- ig[ev[ig] & tv[ig] < tt]
      carry_w <- W[carry]
      if (stab && length(carry) > 0) {
        # W_tau Q_k / Q_tau
        cq <- cumprod(qbar)
        tn <- names(qbar)
        idx_t <- match(as.character(tt), tn)
        idx_s <- match(as.character(tv[carry]), tn)
        carry_w <- carry_w * cq[idx_t] / cq[idx_s]
      }
      rows <- c(base, post_rows, carry)
      frames[[length(frames) + 1]] <- data.frame(
        strategy = sv[rows],
        time = tt,
        .ref_w = c(rep(1, length(base)), rep(0, length(post_rows) +
                                               length(carry))),
        .post_w = c(rep(0, length(base)), post_w, carry_w),
        data[rows, features, drop = FALSE],
        check.names = FALSE
      )
      counts[[length(counts) + 1]] <- data.frame(
        strategy = sv[base[1]], time = tt, n_risk = length(base),
        n_censored = sum(cens_by_t[as.character(times[times <= tt])])
      )
    }
  }
  list(frame = do.call(rbind, frames), counts = do.call(rbind, counts))
}

.to_long <- function(frame, features) {
  long <- tidyr::pivot_longer(frame, cols = dplyr::all_of(features),
                              names_to = "feature", values_to = ".value")
  long <- long[!is.na(long$.value), , drop = FALSE]
  long$feature <- factor(long$feature, levels = features)
  long
}

.join_counts <- function(x, counts) {
  counts$strategy <- as.character(counts$strategy)
  lv <- if (is.factor(x$strategy)) levels(x$strategy) else NULL
  x$strategy <- as.character(x$strategy)
  out <- dplyr::left_join(x, counts, by = c("strategy", "time"))
  if (!is.null(lv)) out$strategy <- factor(out$strategy, levels = lv)
  out
}
