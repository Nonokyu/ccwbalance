#' Fit the artificial-censoring model and compute the weights
#'
#' Estimates `pi_k = P(D_k = 0 | R_k = 1, H_k)` by logistic regression on the
#' rows where artificial censoring is possible (`censor_possible`),
#' separately within each group defined by `by` (default: strategy x
#' phase). Rows where censoring cannot occur get `pi_hat = 1`.
#'
#' The stabilising numerator `q_k` must be common to all individuals at risk
#' in a stratum (strategy x interval). It is estimated on all at-risk rows of
#' each strategy: with `numerator = NULL` (default) `q_k` is the observed
#' proportion remaining uncensored in the stratum; with a formula it is a
#' logistic model fitted per strategy on all at-risk rows, which should then
#' contain only time terms (a covariate-dependent numerator is not common
#' to the stratum).
#'
#' Cumulative weights are formed within individual x strategy:
#' `W^u_k = prod_{l <= k} 1 / pi_l` (unstabilised) and
#' `W_k = prod_{l <= k} q_l / pi_l` (stabilised). With `truncate`, the
#' chosen cumulative weight is capped at that quantile of its distribution
#' over the retained rows (`D_k = 0`), pooled over strategies and intervals;
#' the cap is applied to both `w_post` and `w_prev`.
#'
#' @param data Long data (see [ccwbalance-package]).
#' @param formula One-sided formula for the censoring model, e.g.
#'   `~ L + I(L^2) + X + L:X`, or a named list of formulas. For each model
#'   group the list element named by the full group label (e.g.
#'   `"g2.grace_end"`) is used, else the element named by the value of the
#'   last `by` column (e.g. `"no_init"`), else the element `"default"`.
#' @param numerator `NULL` (default; observed stratum proportion) or a
#'   one-sided formula for `q_k`, e.g. `~ factor(time)`.
#' @param by Columns defining the groups in which separate censoring models
#'   are fitted.
#' @param stabilized Use stabilised weights for `w_ref` and `w_post`?
#' @param truncate `NULL` (no truncation) or a quantile in (0, 1), e.g.
#'   `0.99`.
#' @param id,strategy,time,censored,censor_possible Column names.
#' @param pi_floor Lower bound applied to `pi_hat`.
#' @return `data`, sorted by strategy, individual and interval, with added
#'   columns
#'   * `pi_hat` (`pi_k`) and `q_hat` (`q_k`);
#'   * `weight_unstab` (`W^u_k`), `weight_stab` (`W_k`) and their values in
#'     the previous interval, `w_prev_unstab` (`W^u_{k-1}`) and
#'     `w_prev_stab` (`W_{k-1}`), all untruncated;
#'   * `weight` and `w_prev`: the chosen scheme's `W_k` and `W_{k-1}`,
#'     truncated if requested;
#'   * `w_post` (= `weight`): the estimator weight of a retained row;
#'   * `w_ref`: the reference weight `w_prev * q_hat` (stabilised) or
#'     `w_prev` (unstabilised);
#'   * `weights_stabilized`: the value of `stabilized`.
#'
#'   The fitted models are attached as attribute `"models"` and the settings
#'   (formula, numerator, stabilisation, truncation, cap, groups) as
#'   attribute `"weighting"`.
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

  # stabilising numerator: common to all at-risk rows of a stratum
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
  lag1 <- function(x) stats::ave(x, clone,
                                 FUN = function(z) c(1, z[-length(z)]))
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

# Logistic fit returning fitted probabilities (constant if the outcome does
# not vary)
.fit_logit <- function(f, sub) {
  if (length(unique(sub$.uncens)) < 2) {
    return(list(fit = NULL, p = rep(mean(sub$.uncens), nrow(sub))))
  }
  fit <- suppressWarnings(stats::glm(f, family = stats::binomial(), data = sub))
  list(fit = fit, p = unname(stats::predict(fit, type = "response")))
}

#' Weighted discrete-time risk by strategy
#'
#' Weighted Kaplan-Meier-type estimator. The hazard of interval `k` is the
#' weighted proportion of outcomes among the rows retained in interval `k`
#' (`D_k = 0` and, if `ltfu` is present, not lost to follow-up), and the
#' risk by the end of interval `k` is `1 - prod_{l <= k} (1 - h_l)`.
#'
#' @param data Long data with weights (see [fit_censoring_model()]).
#' @param weight Weight column; `NULL` for unweighted.
#' @param strategy,time,censored,event,ltfu Column names (`ltfu` is ignored
#'   if `NULL` or absent).
#' @return A data frame with `strategy, time, hazard, risk`.
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
