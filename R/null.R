#' Conditional parametric null reference for interval RRB and censoring
#' balance
#'
#' Even under a correctly specified censoring model the interval RRB of a
#' finite sample is not exactly 0: its signed value is centred at 0, but its
#' absolute value and its maximum over strata (MaxIB) are positive, of order
#' `ESS^{-1/2}`. `rrb_null()` gives the distribution of the interval RRB,
#' the censoring balance (CB) and their maxima that would arise if the
#' fitted censoring model were true, so that observed values can be judged
#' against calibrated thresholds.
#'
#' ## Construction
#' For each draw `b = 1, ..., n_draws`:
#' 1. on the rows where artificial censoring is possible, draw
#'    `D*_k ~ Bernoulli(1 - pi_hat)` (other rows keep `D*_k = 0`; see
#'    Coupling);
#' 2. refit the censoring model to `D*_k` with the same formula(s) and model
#'    groups as [fit_censoring_model()], giving `pi*_k`;
#' 3. keep the observed reference weights `w_ref` (`W_{k-1} q_k`) and give
#'    the rows retained under `D*_k` the weight `w_ref / pi*_k`;
#' 4. compute the interval RRB and CB of every strategy x interval x
#'    feature cell.
#'
#' The observed values are computed in the same way from `D_k` and
#' `pi_hat`: the observed RRB has post weight `w_ref / pi_hat`, the model
#' weight of the current interval applied to the reference weights. With
#' truncated weights it therefore differs from `rrb(method = "weighted")`,
#' which gives the estimator-weight RRB; without truncation the two agree.
#' The CB of a cell is `(mu_C - mu_U) / s_R`, with censored rows weighted
#' by `w_ref / (1 - pi)` (`1 - pi` floored at `1e-6`), retained rows by
#' `w_ref / pi` and `s_R` the `w_ref`-weighted standard deviation of the
#' stratum; without truncation the observed CB equals
#' `censoring_balance(weighted = TRUE, weights = "model")`.
#'
#' Only strata with at least one artificially censored row in the observed
#' data are used, so the maxima are taken over the same cells in the data
#' and in the draws. A cell without censored rows in a draw contributes
#' `RRB* = 0` and `CB* = NA`; cells with an undefined observed value are
#' left out of the maxima.
#'
#' The reference is conditional on the observed risk sets and reference
#' weights: the redrawn censoring is not propagated to later risk sets or to
#' `W_{k-1}`. It is therefore defined for the interval scope only. It
#' reflects the variability of the censoring indicators and of the
#' censoring-model fit, but not that of the covariates or of the earlier
#' weights.
#'
#' ## Coupling
#' The strategies of one individual are censored by the same treatment
#' event: at the end of a grace period, for example, the individual is
#' censored under "initiate within the window" if it does not initiate and
#' under "do not initiate" if it does. With `coupling = "individual"`
#' (default) the censorable rows of one individual in one interval (`id` x
#' `time`) share a uniform variate `U`, and `U < P(A_k = 1 | H_k)`
#' represents treatment in interval `k`. A row of a group censored on
#' `A_k = 1` (where `P(A_k = 1) = 1 - pi_hat`) is censored if
#' `U < 1 - pi_hat`, and a row of a group censored on `A_k = 0` if
#' `U >= pi_hat`. Each row keeps its marginal probability `1 - pi_hat`;
#' when the models of two strategies agree (`pi_1 = 1 - pi_2`) the
#' indicators are exactly complementary. A row that is the only censorable
#' row of its `id` x `time` group is an independent draw. With
#' `coupling = "independent"` every censorable row gets its own Bernoulli
#' draw.
#'
#' The orientation of each strategy x phase group (censored on `A_k = 1` or
#' on `A_k = 0`) is the value of the `treatment` column among its censored
#' rows. A group without censored rows, and every group when the `treatment`
#' column is absent, falls back to the phase rule: `"grace_end"` is censored
#' on `A_k = 0`, any other phase on `A_k = 1`. This rule fits the
#' initiation-window strategies of [ccw_clone()] and [simulate_ccw()] but
#' not a sustained "always treat" strategy, which needs the `treatment`
#' column. If the censored rows of a group have both treatment values,
#' `rrb_null()` stops; use `coupling = "independent"`.
#'
#' @param data Long data returned by [fit_censoring_model()]; its
#'   `"weighting"` attribute supplies `formula` and `by` unless given.
#' @param features Character vector of feature columns `f(H_k)`.
#' @param n_draws Number of null draws.
#' @param seed Optional random seed.
#' @param strategy,time,censored,censor_possible,pi_hat,w_ref Column names.
#' @param by Model-group columns (default: those of [fit_censoring_model()],
#'   else `c(strategy, "phase")`).
#' @param formula Censoring-model formula or named list of formulas
#'   (default: the one of [fit_censoring_model()]).
#' @param pi_floor Lower bound applied to the refitted `pi*_k`.
#' @param coupling `"individual"` (default) or `"independent"`; see
#'   Coupling.
#' @param id Individual identifier column (used by
#'   `coupling = "individual"`).
#' @param treatment Binary treatment column `A_k` used to orient the
#'   coupling; if absent the phase rule is used.
#' @return An object of class `"rrb_null"`, a list with
#'   * `cells`: one row per cell with `strategy, time, feature, n_risk,
#'     n_censored, rrb_obs, null_mean, null_sd, null_q025, null_q975, z,
#'     p_two`, where `z = rrb_obs / null_sd` and
#'     `p_two = (1 + #{|RRB*| >= |RRB_obs|}) / (1 + n_draws)`, followed by
#'     `cb_obs, cb_null_mean, cb_null_sd, cb_z, cb_p_two`, the same
#'     summaries for CB (`NA` where the observed CB is undefined);
#'   * `maxib`: MaxIB of interval RRB per strategy x feature, per strategy,
#'     per feature and overall (`"(all)"` labels), with `n_cells, maxib_obs,
#'     null_mean, null_q95, p`, where
#'     `p = (1 + #{MaxIB* >= MaxIB_obs}) / (1 + n_draws)`;
#'   * `maxib_cb`: as `maxib`, for CB;
#'   * `draws`: long data frame `draw, strategy, time, feature, rrb, cb`;
#'   * `n_draws`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' nr <- rrb_null(w, c("X", "L"), n_draws = 20, seed = 1)
#' nr
#' head(nr$cells)
rrb_null <- function(data, features, n_draws = 200, seed = NULL,
                     strategy = "strategy", time = "time",
                     censored = "censored_next",
                     censor_possible = "censor_possible", pi_hat = "pi_hat",
                     w_ref = "w_ref", by = NULL, formula = NULL,
                     pi_floor = 1e-6,
                     coupling = c("individual", "independent"),
                     id = "id", treatment = "A") {
  coupling <- match.arg(coupling)
  wt <- attr(data, "weighting")
  if (is.null(formula)) formula <- wt$formula
  if (is.null(formula)) {
    stop("No censoring-model formula: supply `formula` or use data returned ",
         "by fit_censoring_model().", call. = FALSE)
  }
  if (is.null(by)) by <- if (!is.null(wt$by)) wt$by else c(strategy, "phase")
  .check_cols(data, c(strategy, time, censored, censor_possible, pi_hat,
                      w_ref, features, by,
                      if (coupling == "individual") id))
  stopifnot(n_draws >= 1)
  data <- as.data.frame(data)

  # observed interval RRB (post weight w_ref / pi_hat)
  data$.w_post_null <- data[[w_ref]] / data[[pi_hat]]
  obs <- rrb(data, features, method = "weighted", scope = "interval",
             strategy = strategy, time = time, censored = censored,
             w_ref = w_ref, w_post = ".w_post_null")
  obs <- obs[obs$n_censored >= 1 & !is.na(obs$rrb), , drop = FALSE]
  if (nrow(obs) == 0) stop("No cells with artificial censoring.", call. = FALSE)

  prep <- .null_prep(data, features, strategy, time, w_ref)
  groups <- .null_groups(data, formula, by, censor_possible)
  poss <- as.logical(data[[censor_possible]])
  p_obs <- data[[pi_hat]]
  n <- nrow(data)
  if (coupling == "individual") {
    key <- paste(data[[id]][poss], data[[time]][poss], sep = "\r")
    grp_poss <- match(key, unique(key))
    n_grp <- max(c(0L, grp_poss))
    on0_poss <- .null_orient(data, strategy, censored, censor_possible,
                             treatment)[poss]
  }

  cell_obs <- match(paste(as.character(obs$strategy), obs$time, sep = "\r"),
                    prep$cells)
  feat_obs <- match(as.character(obs$feature), features)
  pos <- cbind(cell_obs, feat_obs)

  if (!is.null(seed)) set.seed(seed)
  sims <- matrix(NA_real_, nrow(obs), n_draws)
  sims_cb <- sims
  # one null draw: redraw censoring from `p`, refit, return RRB* and CB*
  world <- function(p) {
    d_star <- integer(n)
    d_star[poss] <- if (coupling == "individual") {
      .null_draw(p[poss], on0_poss, grp_poss, stats::runif(n_grp))
    } else {
      stats::rbinom(sum(poss), 1, 1 - p[poss])
    }
    pi_star <- pmax(.null_refit(groups, 1 - d_star, n), pi_floor)
    list(r = .null_rrb(prep, d_star, pi_star)[pos],
         cb = .null_cb(prep, d_star, pi_star)[pos])
  }
  for (b in seq_len(n_draws)) {
    wb <- world(p_obs)
    sims[, b] <- wb$r
    sims_cb[, b] <- wb$cb
  }
  cb_obs <- .null_cb(prep, data[[censored]], pmax(p_obs, pi_floor))[pos]

  cells <- data.frame(
    strategy = obs$strategy, time = obs$time, feature = obs$feature,
    n_risk = obs$n_risk, n_censored = obs$n_censored, rrb_obs = obs$rrb,
    null_mean = rowMeans(sims, na.rm = TRUE),
    null_sd = apply(sims, 1, stats::sd, na.rm = TRUE),
    null_q025 = apply(sims, 1, stats::quantile, 0.025, na.rm = TRUE,
                      names = FALSE),
    null_q975 = apply(sims, 1, stats::quantile, 0.975, na.rm = TRUE,
                      names = FALSE)
  )
  cells$z <- ifelse(cells$null_sd > 0, cells$rrb_obs / cells$null_sd,
                    NA_real_)
  cells$p_two <- (1 + rowSums(abs(sims) >= abs(cells$rrb_obs),
                              na.rm = TRUE)) /
    (1 + rowSums(!is.na(sims)))

  cells$cb_obs <- cb_obs
  n_cb <- rowSums(!is.na(sims_cb))
  cells$cb_null_mean <- ifelse(n_cb > 0, rowMeans(sims_cb, na.rm = TRUE),
                               NA_real_)
  cells$cb_null_sd <- apply(sims_cb, 1, function(v) {
    if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE)
  })
  cells$cb_z <- ifelse(!is.na(cells$cb_null_sd) & cells$cb_null_sd > 0,
                       cb_obs / cells$cb_null_sd, NA_real_)
  cells$cb_p_two <- ifelse(
    is.na(cb_obs), NA_real_,
    (1 + rowSums(abs(sims_cb) >= abs(cb_obs), na.rm = TRUE)) / (1 + n_cb)
  )

  maxib <- .null_maxib(cells, sims)
  maxib_cb <- .null_maxib(cells, sims_cb, value = "cb_obs")
  draws <- data.frame(
    draw = rep(seq_len(n_draws), each = nrow(obs)),
    strategy = rep(obs$strategy, n_draws), time = rep(obs$time, n_draws),
    feature = rep(obs$feature, n_draws), rrb = as.vector(sims),
    cb = as.vector(sims_cb)
  )
  structure(list(cells = cells, maxib = maxib, maxib_cb = maxib_cb,
                 draws = draws, n_draws = n_draws),
            class = "rrb_null")
}

#' @export
print.rrb_null <- function(x, ...) {
  cat(sprintf(paste0("Conditional parametric null for interval RRB and ",
                     "censoring balance (%d draws)\n"), x$n_draws))
  cat("\nMaxIB, Interval RRB (observed vs null):\n")
  print(x$maxib, digits = 3, row.names = FALSE)
  if (!is.null(x$maxib_cb)) {
    cat("\nMaxIB, Censoring balance (observed vs null):\n")
    print(x$maxib_cb, digits = 3, row.names = FALSE)
  }
  big <- x$cells[!is.na(x$cells$z) & abs(x$cells$z) > 2, , drop = FALSE]
  cat(sprintf("\nCells with |z| > 2: %d of %d\n", nrow(big), nrow(x$cells)))
  if (nrow(big) > 0) {
    print(big[, c("strategy", "time", "feature", "n_risk", "n_censored",
                  "rrb_obs", "null_sd", "z", "p_two")],
          digits = 3, row.names = FALSE)
  }
  invisible(x)
}

# Coupled redraw of censoring indicators. `cens_on0` flags rows censored on
# A_k = 0, `grp` indexes the id x time group of each row and `u` holds one
# uniform variate per group. Treatment is U < P(A_k = 1), with
# P(A_k = 1) = pi if censored on A_k = 0 and 1 - pi otherwise, so that
# P(D* = 1) = 1 - pi for every row.
.null_draw <- function(pi, cens_on0, grp, u) {
  uu <- u[grp]
  as.integer(ifelse(cens_on0, uu >= pi, uu < 1 - pi))
}

# Orientation of the coupling: TRUE where artificial censoring occurs on
# A_k = 0, from the treatment value of the censored rows of each
# strategy x phase group (fallback: grace_end -> A_k = 0, else A_k = 1).
.null_orient <- function(data, strategy, censored, censor_possible,
                         treatment) {
  has_phase <- "phase" %in% names(data)
  phase <- if (has_phase) as.character(data$phase) else rep("", nrow(data))
  out <- phase == "grace_end"
  if (is.null(treatment) || !treatment %in% names(data)) return(out)
  poss <- as.logical(data[[censor_possible]])
  grp <- interaction(data[c(strategy, if (has_phase) "phase")], drop = TRUE,
                     sep = ".")
  a <- data[[treatment]]
  cens <- data[[censored]] %in% 1
  for (g in levels(droplevels(grp[poss]))) {
    idx <- which(poss & grp == g)
    v <- unique(a[idx][cens[idx]])
    v <- v[!is.na(v)]
    if (length(v) == 0) next
    if (length(v) > 1 || !v %in% c(0, 1)) {
      stop(sprintf(paste0("Artificial censoring in group '%s' is not a ",
                          "function of `%s`; use coupling = 'independent'."),
                   g, treatment), call. = FALSE)
    }
    out[idx] <- v == 0
  }
  out
}

# Draw-invariant pieces of the interval RRB
.null_prep <- function(data, features, strategy, time, w_ref) {
  key <- paste(as.character(data[[strategy]]), data[[time]], sep = "\r")
  cells <- unique(key)
  cell <- match(key, cells)
  fm <- as.matrix(data[, features, drop = FALSE])
  storage.mode(fm) <- "double"
  ok <- !is.na(fm)
  f0 <- fm
  f0[!ok] <- 0
  wr <- data[[w_ref]]
  wr_ok <- wr * ok
  sw <- rowsum(wr_ok, cell, reorder = FALSE)
  mu <- rowsum(wr_ok * f0, cell, reorder = FALSE) / sw
  dev2 <- (f0 - mu[cell, , drop = FALSE])^2
  sd_pre <- sqrt(rowsum(wr_ok * dev2, cell, reorder = FALSE) / sw)
  sd_pre[!is.finite(sd_pre) | sd_pre <= 1e-10] <- NA_real_
  list(cells = cells, cell = cell, f0 = f0, ok = ok, wr = wr, mu = mu,
       sd_pre = sd_pre)
}

# Interval RRB matrix (cells x features) for censoring D and pi
.null_rrb <- function(prep, d, pi) {
  wp <- (1 - d) * prep$wr / pi * prep$ok
  mu_post <- rowsum(wp * prep$f0, prep$cell, reorder = FALSE) /
    rowsum(wp, prep$cell, reorder = FALSE)
  r <- (mu_post - prep$mu) / prep$sd_pre
  r[!is.finite(r)] <- NA_real_
  no_cens <- rowsum(d, prep$cell, reorder = FALSE)[, 1] == 0
  r[no_cens, ] <- ifelse(is.na(prep$sd_pre[no_cens, ]), NA_real_, 0)
  r
}

# Censoring-balance matrix (cells x features): (mu_C - mu_U) / s_R; NA where
# a cell has no censored or no retained row with an observed value
.null_cb <- function(prep, d, pi) {
  wk <- (1 - d) * prep$wr / pi * prep$ok
  wc <- d * prep$wr / pmax(1 - pi, 1e-6) * prep$ok
  mu_k <- rowsum(wk * prep$f0, prep$cell, reorder = FALSE) /
    rowsum(wk, prep$cell, reorder = FALSE)
  mu_c <- rowsum(wc * prep$f0, prep$cell, reorder = FALSE) /
    rowsum(wc, prep$cell, reorder = FALSE)
  r <- (mu_c - mu_k) / prep$sd_pre
  r[!is.finite(r)] <- NA_real_
  r
}

# Design matrices of the censoring-model groups of fit_censoring_model()
.null_groups <- function(data, formula, by, censor_possible) {
  poss <- as.logical(data[[censor_possible]])
  grp <- if (length(by) > 0) {
    interaction(data[by], drop = TRUE)
  } else {
    factor(rep("all", nrow(data)))
  }
  pick <- function(label, last) {
    if (!is.list(formula)) return(formula)
    for (nm in c(label, last, "default")) {
      if (!is.null(formula[[nm]])) return(formula[[nm]])
    }
    stop("No formula for censoring-model group '", label, "'.", call. = FALSE)
  }
  lapply(levels(droplevels(grp[poss])), function(g) {
    idx <- which(poss & grp == g)
    sub <- data[idx, , drop = FALSE]
    last <- if (length(by) > 0) as.character(sub[[by[length(by)]]][1]) else g
    tt <- stats::delete.response(stats::terms(pick(g, last)))
    x <- stats::model.matrix(tt, stats::model.frame(tt, sub,
                                                    na.action = stats::na.pass))
    if (anyNA(x)) {
      stop("Missing covariate values in censoring-model group '", g, "'.",
           call. = FALSE)
    }
    list(idx = idx, x = x)
  })
}

# Refit the censoring models to a new uncensored indicator
.null_refit <- function(groups, uncens, n) {
  p <- rep(1, n)
  for (g in groups) {
    y <- uncens[g$idx]
    p[g$idx] <- if (length(unique(y)) < 2) {
      rep(mean(y), length(y))
    } else {
      suppressWarnings(stats::glm.fit(g$x, y,
                                      family = stats::binomial()))$fitted.values
    }
  }
  p
}

# Cell selections of the MaxIB pooling levels (strategy x feature,
# strategy, feature, overall); `sel` excludes cells with an NA observed value
.null_maxib_sets <- function(cells, value = "rrb_obs") {
  obs_v <- cells[[value]]
  strat <- as.character(cells$strategy)
  feat <- as.character(cells$feature)
  lv_s <- if (is.factor(cells$strategy)) {
    intersect(levels(cells$strategy), strat)
  } else {
    sort(unique(strat))
  }
  lv_f <- if (is.factor(cells$feature)) {
    intersect(levels(cells$feature), feat)
  } else {
    unique(feat)
  }
  grid <- rbind(
    expand.grid(strategy = lv_s, feature = lv_f, stringsAsFactors = FALSE),
    data.frame(strategy = lv_s, feature = "(all)"),
    data.frame(strategy = "(all)", feature = lv_f),
    data.frame(strategy = "(all)", feature = "(all)")
  )
  sets <- lapply(seq_len(nrow(grid)), function(k) {
    sel <- (grid$strategy[k] == "(all)" | strat == grid$strategy[k]) &
      (grid$feature[k] == "(all)" | feat == grid$feature[k])
    if (!any(sel)) return(NULL)
    list(strategy = grid$strategy[k], feature = grid$feature[k],
         sel = sel & !is.na(obs_v))
  })
  sets[!vapply(sets, is.null, logical(1))]
}

# Column-wise maximum of |sims| over the selected cells
.null_colmax <- function(sims, sel) {
  apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
}

# Observed vs null MaxIB at each pooling level; `value` names the observed
# column of `cells`
.null_maxib <- function(cells, sims, value = "rrb_obs") {
  obs_v <- cells[[value]]
  rows <- lapply(.null_maxib_sets(cells, value), function(s) {
    sel <- s$sel
    if (!any(sel)) {
      return(data.frame(
        strategy = s$strategy, feature = s$feature,
        n_cells = 0L, maxib_obs = NA_real_, null_mean = NA_real_,
        null_q95 = NA_real_, p = NA_real_
      ))
    }
    m_obs <- max(abs(obs_v[sel]))
    m_null <- .null_colmax(sims, sel)
    data.frame(
      strategy = s$strategy, feature = s$feature,
      n_cells = sum(sel), maxib_obs = m_obs,
      null_mean = mean(m_null, na.rm = TRUE),
      null_q95 = stats::quantile(m_null, 0.95, na.rm = TRUE, names = FALSE),
      p = (1 + sum(m_null >= m_obs, na.rm = TRUE)) / (1 + sum(!is.na(m_null)))
    )
  })
  do.call(rbind, rows)
}
