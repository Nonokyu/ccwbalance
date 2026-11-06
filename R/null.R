#' Conditional parametric null reference for interval RRB, censoring balance
#' and their maxima
#'
#' Even under a correctly specified censoring model the interval RRB of a
#' finite sample is not exactly 0: its signed value is centred at 0, but its
#' absolute value and the maximum over cells (MaxIB) are positive, of order
#' `ESS^{-1/2}`. `rrb_null()` builds a null reference that says how large the
#' interval RRB and MaxIB would be if the fitted censoring model were the true
#' one, so that observed values can be judged against it (calibrated MaxIB
#' thresholds and z-type statistics). The same draws give a null reference
#' for the Jackson-type censoring balance (CB) of each cell, so that the
#' proposed diagnostic and censoring balance can be judged on the same
#' calibrated scale.
#'
#' ## Construction (conditional parametric null)
#' For `b = 1, ..., n_draws`:
#' 1. on rows where artificial censoring is possible, redraw the censoring
#'    indicator with `P(D*_{t+1} = 1) = 1 - pi_hat` from the fitted model
#'    (rows where censoring is impossible keep `D* = 0`). With the default
#'    `coupling = "individual"` the censorable rows of one individual at one
#'    time (`id` x `time`, i.e. the clones of that individual) share one
#'    uniform variate `U ~ Uniform(0, 1)`, and `U < P(A_t = 1 | H_t)` plays
#'    the role of the shared treatment event "treated in interval `t`". A row
#'    of a strategy x phase group that is censored on `A_t = 1` (where
#'    `P(A_t = 1) = 1 - pi_hat`) is censored iff `U < 1 - pi_hat`; a row of a
#'    group censored on `A_t = 0` (`P(A_t = 1) = pi_hat`) iff `U >= pi_hat`
#'    (see Within-individual dependence for how the orientation is found).
#'    With `coupling = "independent"` every censorable row gets its own
#'    Bernoulli draw;
#' 2. refit the censoring model to `D*` with the same formula(s) and model
#'    groups (`by`) as [fit_censoring_model()], giving `pi*`;
#' 3. keep the reference weights `w_ref` (`W_{t-1} q_t`) of the observed data
#'    fixed and give the clones retained under `D*` the weight
#'    `w_ref / pi*`;
#' 4. compute the interval RRB `RRB*` and the censoring balance `CB*` of
#'    every strategy x time x feature cell.
#'
#' ## Censoring balance
#' The Jackson-type censoring balance of a cell is
#' `CB = (mu_C - mu_K) / s_R`, where `mu_K` is the mean of `f(H_t)` among the
#' retained rows (`D = 0`) weighted by `w_ref / pi`, `mu_C` the mean among the
#' censored rows (`D = 1`) weighted by `w_ref / (1 - pi)` (with `1 - pi`
#' floored at `1e-6`), and `s_R` the `w_ref`-weighted SD of `f(H_t)` over the
#' whole cell (the RRB denominator). The observed CB uses `D` and `pi_hat`
#' and equals `censoring_balance(weighted = TRUE, weights = "model")` (the
#' common factor `q_t` of `w_ref` cancels within a cell, untruncated
#' weights); `CB*` uses `D*` and `pi*` of the same draw (no additional
#' random numbers are drawn, so the RRB output does not depend on whether CB
#' is computed). CB is undefined in a cell without censored rows, so such a
#' cell contributes `CB* = NA` in that draw (ignored by the maxima and the
#' null summaries).
#'
#' The censoring balance is returned for two row sets. The `cb_*` elements
#' use every row of the risk set. The `cb_cens_*` elements use only the
#' rows at which artificial censoring is possible (`censor_possible ==
#' TRUE`), so that rows that cannot be artificially censored (for example,
#' individuals who already initiated) do not enter the retained group; their
#' observed value equals `censoring_balance(weighted = TRUE, rows =
#' "censorable")`. Both versions share the standardising SD `s_R` of the
#' whole risk set and are computed from the same draws.
#'
#' The observed RRB is computed in the same way from the observed `D` and
#' `pi_hat` (post weight `w_ref / pi_hat`; this equals `rrb(method =
#' "weighted")` when the weights are not truncated). Only cells with at least
#' one artificially censored clone in the observed data are used, so MaxIB is
#' taken over the same cells in the data and in the draws; a cell without any
#' censored clone in a draw contributes `RRB* = 0` (it would not be a
#' diagnostic cell), and cells whose observed RRB is `NA` (zero reference
#' variance) are dropped.
#'
#' The reference is **conditional** on the observed risk sets and reference
#' weights: redrawn censoring is not propagated to later risk sets or to
#' `W_{t-1}`. It is therefore defined for the **interval scope only** (the
#' 0-targeted local reconstruction; see [rrb()]). It captures the sampling
#' variability of the censoring indicators and of the censoring-model fit, but
#' not that of covariates or of earlier-interval weights.
#'
#' ## Within-individual dependence
#' The clones of one individual are artificially censored by functions of
#' the same treatment history: at the end of a grace period, for example,
#' the clone of an "initiate within the window" strategy is censored if the
#' individual does not initiate and the clone of a "do not initiate"
#' strategy is censored if the individual does. Redrawing every censorable
#' row independently would ignore this dependence, so the null distribution
#' of statistics that combine strategies (MaxIB over strategies) would not
#' reflect the joint behaviour of the clones. The default
#' `coupling = "individual"` preserves each row's marginal censoring
#' probability `1 - pi_hat` exactly (cell-level null marginals are
#' unchanged) and ties the clones of an individual together through a
#' shared uniform variate: when the strategy-specific models agree
#' (`pi_grace_end = 1 - pi_no_init`) the indicators are exactly
#' complementary, and otherwise the coupling is the maximally dependent one
#' compatible with both marginals. A row that is the only censorable row of
#' its `id` x `time` group gets its own variate, which is the independent
#' draw. `coupling = "independent"` reproduces the earlier behaviour
#' (independent Bernoulli draws per row).
#'
#' The orientation of each strategy x phase group (censored on `A_t = 1` or
#' on `A_t = 0`) is read from the data: it is the value of the `treatment`
#' column among the artificially censored rows of the group. In window mode
#' ([ccw_clone()], [simulate_ccw()]) `"no_init"` rows are censored on
#' `A_t = 1` and `"grace_end"` rows on `A_t = 0`; in sustained mode the
#' "always treat" strategy is censored on `A_t = 0` and the "never treat"
#' strategy on `A_t = 1`, so their clones are coupled as complementary
#' events. A group without censored rows falls back to the phase rule
#' (`"grace_end"`: `A_t = 0`; any other phase: `A_t = 1`), and so does every
#' group when the `treatment` column is absent (which is correct for window
#' mode, but not for sustained "always treat" strategies). If the censored
#' rows of a group have both treatment values, censoring is not a function
#' of the treatment and `rrb_null()` stops; use `coupling = "independent"`.
#'
#' @param data Long clone data returned by [fit_censoring_model()] (the
#'   `"weighting"` attribute supplies `formula` and `by` unless given).
#' @param features Character vector of feature columns `f(H_t)`.
#' @param n_draws Number of null draws.
#' @param seed Optional random seed.
#' @param strategy,time,censored,censor_possible,pi_hat,w_ref Column names.
#' @param by Model-group columns (default: those used by
#'   [fit_censoring_model()], else `c("strategy", "phase")`).
#' @param formula Censoring-model formula or named list of formulas (default:
#'   the one used by [fit_censoring_model()]).
#' @param pi_floor Lower bound applied to refitted `pi*`.
#' @param coupling How censoring indicators are redrawn: `"individual"`
#'   (default; clones of one individual at one time share a uniform variate,
#'   see Within-individual dependence) or `"independent"` (one independent
#'   Bernoulli draw per censorable row).
#' @param id Individual identifier column (used by `coupling =
#'   "individual"`).
#' @param treatment Binary treatment column `A_t` used to orient the coupling
#'   (`coupling = "individual"`); if absent the phase rule is used.
#' @return An object of class `"rrb_null"`: a list with
#'   * `cells`: `strategy, time, feature, n_risk, n_censored, rrb_obs,
#'     null_mean, null_sd, null_q025, null_q975, z, p_two`, where
#'     `z = rrb_obs / null_sd` and `p_two` is the two-sided exceedance
#'     probability `(1 + #{|RRB*| >= |RRB_obs|}) / (1 + n_draws)`;
#'   * `maxib`: MaxIB per strategy x feature, per strategy (all features),
#'     per feature (all strategies) and overall (`"(all)"` labels), with
#'     `n_cells, maxib_obs, null_mean, null_q95, p` where
#'     `p = (1 + #{MaxIB* >= MaxIB_obs}) / (1 + n_draws)`;
#'   * `cells` also holds `cb_obs, cb_null_mean, cb_null_sd, cb_z,
#'     cb_p_two`: the same summaries for the censoring balance (CB, see
#'     Censoring balance; `NA` where the observed CB is undefined);
#'   * `maxib_cb`: as `maxib`, for the censoring balance (`n_cells` counts
#'     the cells with a defined observed CB);
#'   * `maxib_stud`, `maxib_cb_stud`: the studentised maximum for the same
#'     pooling levels as `maxib` / `maxib_cb`. Every cell's observed value
#'     and draws are divided by the cell's null SD (`null_sd` /
#'     `cb_null_sd`, estimated from the same draws); cells whose null SD is
#'     `NA` or not positive are left out (`n_excluded`). Columns
#'     `n_cells, n_excluded, stud_obs, null_mean, null_q95, p`, with
#'     `stud_obs` the maximum of `|z|` (`|cb_z|`) over the included cells
#'     and `p` defined as for `maxib`;
#'   * the censorable-row censoring balance (see Censoring balance), from
#'     the same draws: `cells` columns `cb_cens_obs, cb_cens_null_mean,
#'     cb_cens_null_sd, cb_cens_z, cb_cens_p_two` and the maxima
#'     `maxib_cb_cens` (as `maxib_cb`) and `maxib_cb_cens_stud` (as
#'     `maxib_cb_stud`);
#'   * `draws`: long data frame `draw, strategy, time, feature, rrb, cb,
#'     cb_cens` of the null draws;
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
  cens_mask <- poss %in% TRUE
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
  sims_cb_cens <- sims
  # one parametric-bootstrap world: censoring drawn from `p` (all rows),
  # refitted; returns RRB* and CB* (all rows, censorable rows) of the
  # diagnostic cells
  world <- function(p) {
    d_star <- integer(n)
    d_star[poss] <- if (coupling == "individual") {
      .null_draw(p[poss], on0_poss, grp_poss, stats::runif(n_grp))
    } else {
      stats::rbinom(sum(poss), 1, 1 - p[poss])
    }
    pi_star <- pmax(.null_refit(groups, 1 - d_star, n), pi_floor)
    list(r = .null_rrb(prep, d_star, pi_star)[pos],
         cb = .null_cb(prep, d_star, pi_star)[pos],
         cb_cens = .null_cb(prep, d_star, pi_star, cens_mask)[pos])
  }
  for (b in seq_len(n_draws)) {
    wb <- world(p_obs)
    sims[, b] <- wb$r
    sims_cb[, b] <- wb$cb
    sims_cb_cens[, b] <- wb$cb_cens
  }
  p_obs_fl <- pmax(p_obs, pi_floor)
  cb_obs <- .null_cb(prep, data[[censored]], p_obs_fl)[pos]
  cb_cens_obs <- .null_cb(prep, data[[censored]], p_obs_fl, cens_mask)[pos]

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

  cb_sum <- .null_cb_summary(cb_obs, sims_cb)
  cells[names(cb_sum)] <- cb_sum
  cb_cens_sum <- .null_cb_summary(cb_cens_obs, sims_cb_cens)
  names(cb_cens_sum) <- sub("^cb_", "cb_cens_", names(cb_cens_sum))
  cells[names(cb_cens_sum)] <- cb_cens_sum

  maxib <- .null_maxib(cells, sims)
  maxib_cb <- .null_maxib(cells, sims_cb, value = "cb_obs")
  maxib_stud <- .null_maxib(cells, sims, scale = cells$null_sd)
  maxib_cb_stud <- .null_maxib(cells, sims_cb, value = "cb_obs",
                               scale = cells$cb_null_sd)
  maxib_cb_cens <- .null_maxib(cells, sims_cb_cens, value = "cb_cens_obs")
  maxib_cb_cens_stud <- .null_maxib(cells, sims_cb_cens,
                                    value = "cb_cens_obs",
                                    scale = cells$cb_cens_null_sd)
  draws <- data.frame(
    draw = rep(seq_len(n_draws), each = nrow(obs)),
    strategy = rep(obs$strategy, n_draws), time = rep(obs$time, n_draws),
    feature = rep(obs$feature, n_draws), rrb = as.vector(sims),
    cb = as.vector(sims_cb), cb_cens = as.vector(sims_cb_cens)
  )
  structure(list(cells = cells, maxib = maxib, maxib_cb = maxib_cb,
                 maxib_stud = maxib_stud, maxib_cb_stud = maxib_cb_stud,
                 maxib_cb_cens = maxib_cb_cens,
                 maxib_cb_cens_stud = maxib_cb_cens_stud,
                 draws = draws, n_draws = n_draws),
            class = "rrb_null")
}

# Internal: per-cell null summaries of a censoring balance (observed values
# `obs`, draws `sims`, cells x draws): columns cb_obs, cb_null_mean,
# cb_null_sd, cb_z, cb_p_two
.null_cb_summary <- function(obs, sims) {
  n_cb <- rowSums(!is.na(sims))
  null_sd <- apply(sims, 1, function(v) {
    if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE)
  })
  data.frame(
    cb_obs = obs,
    cb_null_mean = ifelse(n_cb > 0, rowMeans(sims, na.rm = TRUE), NA_real_),
    cb_null_sd = null_sd,
    cb_z = ifelse(!is.na(null_sd) & null_sd > 0, obs / null_sd, NA_real_),
    cb_p_two = ifelse(
      is.na(obs), NA_real_,
      (1 + rowSums(abs(sims) >= abs(obs), na.rm = TRUE)) / (1 + n_cb)
    )
  )
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
  if (!is.null(x$maxib_stud)) {
    cat("\nStudentised maximum, Interval RRB (max |z|, observed vs null):\n")
    print(x$maxib_stud, digits = 3, row.names = FALSE)
  }
  if (!is.null(x$maxib_cb_stud)) {
    cat(paste0("\nStudentised maximum, Censoring balance ",
               "(max |cb_z|, observed vs null):\n"))
    print(x$maxib_cb_stud, digits = 3, row.names = FALSE)
  }
  if (!is.null(x$maxib_cb_cens)) {
    cat(paste0("\nMaxIB, Censoring balance over censorable rows ",
               "(observed vs null):\n"))
    print(x$maxib_cb_cens, digits = 3, row.names = FALSE)
  }
  if (!is.null(x$maxib_cb_cens_stud)) {
    cat(paste0("\nStudentised maximum, Censoring balance over censorable ",
               "rows (max |cb_cens_z|, observed vs null):\n"))
    print(x$maxib_cb_cens_stud, digits = 3, row.names = FALSE)
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

# Internal: coupled redraw of censoring indicators. `pi` is the probability
# of remaining uncensored, `cens_on0` flags rows censored on A_t = 0 (else
# censored on A_t = 1), `grp` indexes the id x time group of each row and `u`
# holds one uniform variate per group. The shared event is
# "treated in t" <=> U < P(A_t = 1), with P(A_t = 1) = pi if censored on
# A_t = 0 and 1 - pi if censored on A_t = 1, so each row keeps
# P(D* = 1) = 1 - pi exactly.
.null_draw <- function(pi, cens_on0, grp, u) {
  uu <- u[grp]
  as.integer(ifelse(cens_on0, uu >= pi, uu < 1 - pi))
}

# Internal: orientation of the coupling. Returns a logical vector (all rows)
# that is TRUE where artificial censoring occurs on A_t = 0. Determined per
# strategy x phase group from the treatment value of the censored rows; the
# phase rule (grace_end -> A_t = 0, otherwise A_t = 1) is the fallback.
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

# Internal: fixed (draw-invariant) pieces of the interval RRB --------------
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

# Internal: interval RRB matrix (cells x features) for censoring D and pi
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

# Internal: censoring-balance matrix (cells x features) for censoring D and
# pi: (mu_C - mu_K) / s_R with retained rows weighted by w_ref / pi and
# censored rows by w_ref / (1 - pi); NA where a cell has no censored (or no
# retained) row with an observed feature value. `mask` (logical, one per
# row; NULL = all rows) restricts both groups to the masked rows; s_R is
# always the w_ref-weighted SD of the whole risk set (prep$sd_pre)
.null_cb <- function(prep, d, pi, mask = NULL) {
  ok <- if (is.null(mask)) prep$ok else prep$ok * mask
  wk <- (1 - d) * prep$wr / pi * ok
  wc <- d * prep$wr / pmax(1 - pi, 1e-6) * ok
  mu_k <- rowsum(wk * prep$f0, prep$cell, reorder = FALSE) /
    rowsum(wk, prep$cell, reorder = FALSE)
  mu_c <- rowsum(wc * prep$f0, prep$cell, reorder = FALSE) /
    rowsum(wc, prep$cell, reorder = FALSE)
  r <- (mu_c - mu_k) / prep$sd_pre
  r[!is.finite(r)] <- NA_real_
  r
}

# Internal: design matrices of the censoring-model groups (as in
# fit_censoring_model(); the response is supplied per draw)
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

# Internal: refit the censoring models to a new uncensored indicator
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

# Internal: cell selections of the MaxIB pooling levels (per strategy x
# feature, per strategy, per feature, overall), in the row order of
# .null_maxib(); levels without any cell are omitted and `sel` excludes
# cells whose observed value (`value` column of `cells`) is NA
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

# Internal: column-wise maximum of |sims| over the selected cells (NA for a
# column without any non-NA selected value)
.null_colmax <- function(sims, sel) {
  apply(abs(sims[sel, , drop = FALSE]), 2, function(v) {
    if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  })
}

# Internal: observed vs null MaxIB at several pooling levels; `value` names
# the observed-value column of `cells` (cells with NA observed values are
# left out of both the observed and the null maxima). With `scale` (one
# value per cell, e.g. the null SD) the studentised maximum is returned:
# observed values and draws are divided by `scale`, cells with NA or
# non-positive `scale` are excluded (counted in `n_excluded`) and the
# statistic column is `stud_obs` instead of `maxib_obs`
.null_maxib <- function(cells, sims, value = "rrb_obs", scale = NULL) {
  obs_v <- cells[[value]]
  stud <- !is.null(scale)
  if (stud) {
    ok_s <- !is.na(scale) & scale > 0
    sc <- ifelse(ok_s, scale, 1)
    obs_v <- obs_v / sc
    sims <- sims / sc
  }
  rows <- lapply(.null_maxib_sets(cells, value), function(s) {
    sel <- s$sel
    n_exc <- 0L
    if (stud) {
      n_exc <- sum(sel & !ok_s)
      sel <- sel & ok_s
    }
    if (!any(sel)) {
      out <- data.frame(
        strategy = s$strategy, feature = s$feature,
        n_cells = 0L, maxib_obs = NA_real_, null_mean = NA_real_,
        null_q95 = NA_real_, p = NA_real_
      )
    } else {
      m_obs <- max(abs(obs_v[sel]))
      m_null <- .null_colmax(sims, sel)
      out <- data.frame(
        strategy = s$strategy, feature = s$feature,
        n_cells = sum(sel), maxib_obs = m_obs,
        null_mean = mean(m_null, na.rm = TRUE),
        null_q95 = stats::quantile(m_null, 0.95, na.rm = TRUE,
                                   names = FALSE),
        p = (1 + sum(m_null >= m_obs, na.rm = TRUE)) /
          (1 + sum(!is.na(m_null)))
      )
    }
    if (stud) {
      names(out)[names(out) == "maxib_obs"] <- "stud_obs"
      out <- data.frame(out[1:3], n_excluded = n_exc, out[-(1:3)])
    }
    out
  })
  do.call(rbind, rows)
}
