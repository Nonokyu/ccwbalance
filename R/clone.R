#' Build long data under each strategy with artificial censoring
#'
#' Evaluates the observed treatment history of every individual against
#' every strategy and returns long data in the package format (see
#' [ccwbalance-package]), with artificial censoring in the interval of the
#' first deviation. Two input modes are available:
#'
#' * Initiation-window strategies (default): one row per individual with
#'   the intervals of initiation, outcome and loss to follow-up.
#' * Sustained strategies: person-period data (supply `time` and
#'   `treatment`) and strategies created by [ccw_sustained()], such as
#'   "treat in every interval" or "never treat".
#'
#' Window rules, with `tau` the initiation interval (`Inf` if none):
#' * window `c(a, b)` ("initiate in intervals `a, ..., b`"): initiation
#'   before `a` is a deviation in interval `tau` (phase `"no_init"` for
#'   `k < a`); no initiation by the end of `b` is a deviation in interval
#'   `b` (phase `"grace_end"` at `k = b` among those not yet initiated);
#' * window `c(NA, b)` ("do not initiate through interval `b`"):
#'   initiation at `tau <= b` is a deviation in interval `tau` (phase
#'   `"no_init"` for `k <= b`).
#'
#' Within an interval the order is artificial censoring, loss to follow-up,
#' outcome. The exception is controlled by `outcome_first_at_grace_end`:
#' non-initiation by the end of interval `b` is a deviation only at the end
#' of that interval, so by default an individual whose outcome or loss to
#' follow-up falls in interval `b` is not artificially censored there; such
#' rows get phase `"none"` and do not enter the grace-end censoring model.
#' Rows stop at the first of artificial censoring, outcome, loss to
#' follow-up or interval `n_time - 1`.
#'
#' @param data One row per individual (window mode) or one row per
#'   individual x interval (sustained mode).
#' @param strategies Named list of windows, e.g.
#'   `list(g1 = c(0, 2), g2 = c(3, 8), g3 = c(NA, 8))`, or (sustained mode)
#'   of [ccw_sustained()] objects.
#' @param n_time Number of intervals of follow-up (window mode; optional in
#'   sustained mode, where rows with `time >= n_time` are dropped).
#' @param id,init_time,event_time,ltfu_time Column names; times are interval
#'   indices (`0, 1, ...`), `Inf` or `NA` if the event does not occur.
#' @param covariates Columns copied to every row.
#' @param outcome_first_at_grace_end See Details. With `FALSE` artificial
#'   censoring at the end of the window comes first, as in [simulate_ccw()].
#' @param time,treatment Sustained mode: interval index and binary treatment
#'   columns of the person-period data. Supplying `time` selects sustained
#'   mode.
#' @param event,ltfu Sustained mode: outcome indicator of the interval (`NA`
#'   allowed when unobserved) and optional loss-to-follow-up indicator
#'   (lost before the outcome of the interval). Individuals whose records
#'   end without an outcome leave the risk set.
#' @return Long data with `id, strategy, time, at_risk, A, A_lag,
#'   censor_possible, phase, censored_next, ltfu_next, event_next` and the
#'   `covariates` (sustained mode adds `run_off`, see [ccw_sustained()]).
#' @export
#' @examples
#' obs <- data.frame(id = 1:3, tau = c(1, 5, Inf), ev = c(Inf, 4, 2),
#'                   lt = Inf, age = c(70, 80, 65))
#' ccw_clone(obs, list(early = c(0, 2), never = c(NA, 6)), n_time = 7,
#'           init_time = "tau", event_time = "ev", ltfu_time = "lt",
#'           covariates = "age")
#' pp <- data.frame(id = rep(1:2, each = 3), t = rep(0:2, 2),
#'                  a = c(1, 1, 0, 0, 0, 0), y = c(0, 0, 0, 0, 1, NA),
#'                  L = c(0.2, -0.1, 0.4, 1.1, 0.3, -0.5))
#' ccw_clone(pp, list(always = ccw_sustained(1), never = ccw_sustained(0)),
#'           time = "t", treatment = "a", event = "y", covariates = "L")
ccw_clone <- function(data, strategies, n_time = NULL, id = "id",
                      init_time = "init_time", event_time = "event_time",
                      ltfu_time = "ltfu_time", covariates = character(),
                      outcome_first_at_grace_end = TRUE, time = NULL,
                      treatment = NULL, event = NULL, ltfu = NULL) {
  if (is.null(names(strategies)) || any(names(strategies) == "")) {
    stop("`strategies` must be a named list.", call. = FALSE)
  }
  if (!is.null(time)) {
    return(.ccw_clone_sustained(data, strategies, n_time, id, time,
                                treatment, event, ltfu, covariates))
  }
  if (is.null(n_time)) stop("`n_time` is required.", call. = FALSE)
  if (any(vapply(strategies, inherits, logical(1), "ccw_sustained"))) {
    stop("Sustained strategies need person-period data: supply `time` and ",
         "`treatment`.", call. = FALSE)
  }
  .check_cols(data, c(id, init_time, event_time, ltfu_time, covariates))
  fix <- function(x) ifelse(is.na(x), Inf, x)
  tau <- fix(data[[init_time]])
  ev <- fix(data[[event_time]])
  lt <- fix(data[[ltfu_time]])
  n <- nrow(data)
  K <- n_time
  base_t <- rep(seq_len(K) - 1, times = n)
  base_i <- rep(seq_len(n), each = K)
  out <- lapply(names(strategies), function(g) {
    w <- strategies[[g]]
    a <- w[1]
    b <- w[2]
    if (is.na(b) || b > K - 1 || (!is.na(a) && a > b)) {
      stop("Invalid window for strategy ", g, call. = FALSE)
    }
    cens <- if (is.na(a)) {
      ifelse(tau <= b, tau, Inf)
    } else {
      grace <- ifelse(tau > b, b, Inf)
      if (outcome_first_at_grace_end) {
        grace <- ifelse(grace == b & (ev == b | lt == b), Inf, grace)
      }
      ifelse(tau < a, tau, grace)
    }
    last <- pmin(cens, ev, lt, K - 1)
    keep <- base_t <= last[base_i]
    i <- base_i[keep]
    t <- base_t[keep]
    not_init <- tau[i] >= t
    phase <- if (is.na(a)) {
      ifelse(t <= b & not_init, "no_init", "none")
    } else {
      ge <- t == b & not_init
      if (outcome_first_at_grace_end) {
        ge <- ge & !(ev[i] == b | lt[i] == b)
      }
      ifelse(t < a & not_init, "no_init", ifelse(ge, "grace_end", "none"))
    }
    cn <- as.integer(t == cens[i])
    lf <- ifelse(cn == 1, 0L, as.integer(t == lt[i]))
    d <- data.frame(
      id = data[[id]][i], strategy = g, time = t, at_risk = 1L,
      A = as.integer(tau[i] <= t), A_lag = as.integer(tau[i] < t),
      censor_possible = phase != "none", phase = phase,
      censored_next = cn, ltfu_next = lf,
      event_next = ifelse(cn == 1 | lf == 1, NA_integer_,
                          as.integer(t == ev[i]))
    )
    cbind(d, data[i, covariates, drop = FALSE])
  })
  out <- do.call(rbind, out)
  out$strategy <- factor(out$strategy, levels = names(strategies))
  out <- out[order(out$strategy, out$id, out$time), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Sustained treatment strategy for ccw_clone()
#'
#' Defines the strategy "receive treatment level `value` in every interval",
#' optionally allowing up to `grace` consecutive non-adherent intervals. An
#' individual deviates (`censored_next = 1`) in the first interval `k` in
#' which its run of consecutive non-adherent intervals reaches `grace + 1`.
#' Artificial censoring is therefore possible (phase `"sustain"`) in rows
#' whose run of non-adherence up to interval `k - 1` (column `run_off`)
#' equals `grace`; rows with a shorter run get phase `"none"` and
#' `censor_possible = FALSE`. With `grace = 0` every non-adherent interval
#' is a deviation. The censoring model of phase `"sustain"` is a model for
#' `P(A_k = value | H_k)`.
#'
#' Within an interval the order is `L_k`, `A_k`, artificial censoring, loss
#' to follow-up, outcome; a censored row has `event_next = NA`.
#'
#' @param value Treatment level to sustain (1 = always treat, 0 = never
#'   treat).
#' @param grace Number of consecutive non-adherent intervals allowed.
#' @return An object of class `"ccw_sustained"`.
#' @export
#' @examples
#' ccw_sustained(1)
#' ccw_sustained(1, grace = 1)
ccw_sustained <- function(value = 1, grace = 0) {
  stopifnot(length(value) == 1, value %in% c(0, 1), length(grace) == 1,
            grace >= 0, grace == round(grace))
  structure(list(value = value, grace = grace), class = "ccw_sustained")
}

#' @export
print.ccw_sustained <- function(x, ...) {
  cat(sprintf("<ccw_sustained: A = %d in every interval, grace = %d>\n",
              x$value, x$grace))
  invisible(x)
}

.ccw_clone_sustained <- function(data, strategies, n_time, id, time,
                                 treatment, event, ltfu, covariates) {
  if (!all(vapply(strategies, inherits, logical(1), "ccw_sustained"))) {
    stop("In person-period mode all strategies must be ccw_sustained().",
         call. = FALSE)
  }
  if (is.null(treatment) || is.null(event)) {
    stop("Sustained mode needs `treatment` and `event` columns.",
         call. = FALSE)
  }
  .check_cols(data, c(id, time, treatment, event, ltfu, covariates))
  data <- as.data.frame(data)
  if (!is.null(n_time)) data <- data[data[[time]] < n_time, , drop = FALSE]
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  if (anyDuplicated(data[, c(id, time)])) {
    stop("Duplicate id x time rows.", call. = FALSE)
  }
  idv <- data[[id]]
  tv <- data[[time]]
  a <- data[[treatment]]
  y <- data[[event]]
  lf <- if (is.null(ltfu)) rep(0L, nrow(data)) else data[[ltfu]]
  lf[is.na(lf)] <- 0L
  first <- !duplicated(idv)
  a_lag <- c(0, a[-length(a)])
  a_lag[first] <- 0
  # first outcome / loss-to-follow-up interval (later rows are dropped)
  first_time <- function(flag) {
    ft <- tapply(ifelse(!is.na(flag) & flag == 1, tv, Inf), idv, min)
    unname(ft[as.character(idv)])
  }
  ev_t <- first_time(y)
  lt_t <- first_time(lf)
  out <- lapply(names(strategies), function(g) {
    st <- strategies[[g]]
    nonadh <- as.integer(a != st$value)
    # run of consecutive non-adherent intervals up to and including k
    run <- stats::ave(nonadh, idv, FUN = function(z) {
      r <- integer(length(z))
      acc <- 0L
      for (k in seq_along(z)) {
        acc <- if (z[k] == 1L) acc + 1L else 0L
        r[k] <- acc
      }
      r
    })
    run_prev <- c(0L, run[-length(run)])
    run_prev[first] <- 0L
    dev <- run > st$grace
    cens_t <- tapply(ifelse(dev, tv, Inf), idv, min)
    cens_t <- unname(cens_t[as.character(idv)])
    last <- pmin(cens_t, ev_t, lt_t)
    keep <- tv <= last
    cn <- as.integer(tv == cens_t)
    lfn <- ifelse(cn == 1, 0L, as.integer(lf == 1))
    d <- data.frame(
      id = idv, strategy = g, time = tv, at_risk = 1L,
      A = a, A_lag = a_lag, run_off = run_prev,
      censor_possible = run_prev == st$grace,
      phase = ifelse(run_prev == st$grace, "sustain", "none"),
      censored_next = cn, ltfu_next = lfn,
      event_next = ifelse(cn == 1 | lfn == 1, NA_integer_,
                          as.integer(y))
    )
    cbind(d, data[, covariates, drop = FALSE])[keep, , drop = FALSE]
  })
  out <- do.call(rbind, out)
  out$strategy <- factor(out$strategy, levels = names(strategies))
  out <- out[order(out$strategy, out$id, out$time), , drop = FALSE]
  rownames(out) <- NULL
  out
}
