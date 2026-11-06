#' Fit an artificial-censoring model and compute IPCW weights
#'
#' Estimates `pi_hat = P(D_{t+1} = 0 | R_t = 1, H_t, G = g)` by logistic
#' regression among rows where artificial censoring is possible
#' (`censor_possible`), separately within the groups given by `by`
#' (default: strategy x censoring phase). Rows where censoring cannot occur get
#' `pi_hat = q_hat = 1`. Misspecification scenarios are obtained simply by
#' supplying different right-hand-side formulas.
#'
#' The stabilising numerator `q_t` must be common to all clones at risk in a
#' strategy x time cell (otherwise the target population changes, e.g. clones
#' that already initiated would be up-weighted relative to those facing
#' censoring). It is therefore estimated on **all** at-risk rows of each
#' strategy: with `numerator = NULL` (default) `q_t` is the observed proportion
#' of clones remaining uncensored in the strategy x time cell; with a formula
#' it is a logistic model fitted per strategy on all at-risk rows (it should
#' contain only baseline covariates and time terms).
#'
#' Cumulative weights are `W_t = prod_{k <= t} 1 / pi_k` (unstabilised) and
#' `W_t = prod_{k <= t} q_k / pi_k` (stabilised), computed within clone.
#' Optional truncation caps the chosen cumulative weight at the given quantile
#' of its distribution among retained rows (pooled over strategies and times).
#'
#' @param data Long clone data (see [simulate_ccw()]).
#' @param formula One-sided formula for the censoring (denominator) model, e.g.
#'   `~ L + I(L^2) + X + L:X`, or a named list of formulas chosen per model
#'   group: an element named by the full group label (e.g.
#'   `"g2.grace_end"`) is used first, then one named by the last `by`
#'   component (e.g. the phase `"no_init"`), then an element named
#'   `"default"`.
#' @param numerator `NULL` (default; nonparametric strategy x time proportion)
#'   or a one-sided formula for the stabilising numerator `q`, e.g.
#'   `~ factor(time) + X`.
#' @param by Grouping columns within which separate models are fitted.
#' @param stabilized Use stabilised weights for `w_ref`/`w_post`?
#' @param truncate `NULL` or a quantile in (0, 1), e.g. `0.99`.
#' @param id,strategy,time,censored,censor_possible Column names.
#' @param pi_floor Lower bound applied to `pi_hat` for numerical stability.
#' @return `data` with added columns `pi_hat, q_hat, weight_unstab,
#'   weight_stab, w_prev_unstab, w_prev_stab, weight, w_prev, w_ref, w_post,
#'   weights_stabilized`.
#'   `w_post` is the (possibly truncated) cumulative weight `W_t` of the
#'   chosen scheme; `w_ref` is `W_{t-1}` (times `q_t` if stabilised), i.e. the
#'   weight defining the pre-censoring reference risk set. The fitted models
#'   are attached as attribute `"models"`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' summary(w$w_post)
fit_censoring_model <- function(data, formula, numerator = NULL,
                                by = c("strategy", "phase"),
                                stabilized = TRUE, truncate = NULL,
                                id = "id", strategy = "strategy",
                                time = "time", censored = "censored_next",
                                censor_possible = "censor_possible",
                                pi_floor = 1e-6) {
  .check_cols(data, c(id, strategy, time, censored, censor_possible, by))
  data <- data[order(data[[strategy]], data[[id]], data[[time]]), ,
               drop = FALSE]
  rownames(data) <- NULL
  data$.uncens <- 1 - data[[censored]]
  pick_formula <- function(label, last) {
    if (!is.list(formula)) return(formula)
    for (nm in c(label, last, "default")) {
      if (!is.null(formula[[nm]])) return(formula[[nm]])
    }
    stop("No formula for censoring-model group '", label,
         "' (supply a \"default\" element).", call. = FALSE)
  }

  data$pi_hat <- 1
  data$q_hat <- 1
  poss <- as.logical(data[[censor_possible]])
  grp <- if (length(by) > 0) {
    interaction(data[by], drop = TRUE)
  } else {
    factor(rep("all", nrow(data)))
  }
  models <- list()
  for (g in levels(droplevels(grp[poss]))) {
    idx <- which(poss & grp == g)
    sub <- data[idx, , drop = FALSE]
    last <- if (length(by) > 0) as.character(sub[[by[length(by)]]][1]) else g
    f_den <- stats::update(pick_formula(g, last), .uncens ~ .)
    den <- .fit_logit(f_den, sub)
    data$pi_hat[idx] <- den$p
    models[[g]] <- list(denominator = den$fit)
  }
  data$pi_hat <- pmax(data$pi_hat, pi_floor)

  # stabilising numerator: common to all at-risk clones of a strategy x time
  if (is.null(numerator)) {
    data$q_hat <- stats::ave(data$.uncens, data[[strategy]], data[[time]],
                             FUN = mean)
  } else {
    f_num <- stats::update(numerator, .uncens ~ .)
    for (g in unique(as.character(data[[strategy]]))) {
      idx <- which(as.character(data[[strategy]]) == g)
      num <- .fit_logit(f_num, data[idx, , drop = FALSE])
      data$q_hat[idx] <- num$p
      models[[paste0(g, ".numerator")]] <- list(numerator = num$fit)
    }
  }

  clone <- interaction(data[[strategy]], data[[id]], drop = TRUE)
  cp <- function(x) stats::ave(x, clone, FUN = cumprod)
  lag1 <- function(x) stats::ave(x, clone, FUN = function(z) c(1, z[-length(z)]))
  data$weight_unstab <- cp(1 / data$pi_hat)
  data$weight_stab <- cp(data$q_hat / data$pi_hat)
  data$w_prev_unstab <- lag1(data$weight_unstab)
  data$w_prev_stab <- lag1(data$weight_stab)

  w <- if (stabilized) data$weight_stab else data$weight_unstab
  wp <- if (stabilized) data$w_prev_stab else data$w_prev_unstab
  cap <- Inf
  if (!is.null(truncate)) {
    stopifnot(truncate > 0, truncate < 1)
    cap <- stats::quantile(w[data[[censored]] == 0], truncate, names = FALSE)
  }
  data$weight <- pmin(w, cap)
  data$w_prev <- pmin(wp, cap)
  data$w_post <- data$weight
  data$w_ref <- data$w_prev * (if (stabilized) data$q_hat else 1)
  data$weights_stabilized <- stabilized
  data$.uncens <- NULL
  attr(data, "models") <- models
  attr(data, "weighting") <- list(formula = formula, numerator = numerator,
                                  stabilized = stabilized, truncate = truncate,
                                  cap = cap, by = by)
  data
}

# Internal: logistic fit returning fitted probabilities (constant if the
# outcome does not vary)
.fit_logit <- function(f, sub) {
  if (length(unique(sub$.uncens)) < 2) {
    return(list(fit = NULL, p = rep(mean(sub$.uncens), nrow(sub))))
  }
  fit <- suppressWarnings(stats::glm(f, family = stats::binomial(), data = sub))
  list(fit = fit, p = unname(stats::predict(fit, type = "response")))
}

#' Weighted discrete-time risk estimate by strategy
#'
#' Weighted Kaplan-Meier-type estimator: the interval hazard is the weighted
#' proportion of events among clones retained (not artificially censored and
#' not lost to follow-up) in `[t, t+1)`; risk is `1 - prod(1 - h_t)`.
#'
#' @param data Weighted long clone data.
#' @param weight Weight column (`NULL` for unweighted).
#' @param strategy,time,censored,event,ltfu Column names (`ltfu` may be
#'   `NULL`).
#' @return Data frame with `strategy, time, hazard, risk`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' tail(estimate_risk(w))
estimate_risk <- function(data, weight = "w_post", strategy = "strategy",
                          time = "time", censored = "censored_next",
                          event = "event_next", ltfu = "ltfu_next") {
  .check_cols(data, c(strategy, time, censored, event))
  keep <- data[[censored]] == 0
  if (!is.null(ltfu) && ltfu %in% names(data)) keep <- keep & data[[ltfu]] == 0
  d <- data[keep & !is.na(data[[event]]), , drop = FALSE]
  d$.w <- if (is.null(weight)) 1 else d[[weight]]
  d$.y <- d[[event]]
  out <- d |>
    dplyr::group_by(strategy = .data[[strategy]], time = .data[[time]]) |>
    dplyr::summarise(hazard = sum(.data$.w * .data$.y) / sum(.data$.w),
                     .groups = "drop") |>
    dplyr::arrange(.data$strategy, .data$time) |>
    dplyr::mutate(risk = 1 - cumprod(1 - .data$hazard), .by = "strategy")
  as.data.frame(out)
}
