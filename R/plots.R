# Rows of groups with at least two points (for geom_line)
.multi_point <- function(x, by) {
  k <- interaction(x[by], drop = TRUE)
  x[stats::ave(seq_along(k), k, FUN = length) > 1, , drop = FALSE]
}

#' Dynamic Love plot of RRB
#'
#' Plots `|RRB_{j,g}(k)|` against the interval, one line per feature, with
#' one panel column per strategy and one panel row per method (method x
#' scope when `x` has both scopes).
#'
#' @param x Output of [rrb()].
#' @param threshold Height of a reference line (default 0.1); `NULL` for
#'   none.
#' @param censored_times_only Keep only strata with artificial censoring
#'   (`n_censored > 0`)?
#' @return A ggplot object.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' plot_love_dynamic(rrb(w, c("X", "L"), method = c("unweighted", "weighted")))
plot_love_dynamic <- function(x, threshold = 0.1, censored_times_only = TRUE) {
  .check_cols(x, c("strategy", "time", "feature", "method", "rrb"), "x")
  if (censored_times_only && "n_censored" %in% names(x)) {
    x <- x[x$n_censored > 0, , drop = FALSE]
  }
  x <- x[!is.na(x$rrb), , drop = FALSE]
  x$abs_rrb <- abs(x$rrb)
  x$method <- factor(x$method, levels = intersect(
    c("unweighted", "unstabilized", "stabilized", "weighted"),
    unique(x$method)))
  scopes <- if ("scope" %in% names(x)) unique(x$scope) else "interval"
  if (length(scopes) > 1) {
    x$panel <- interaction(x$method, x$scope, sep = " / ", lex.order = TRUE)
    facet <- ggplot2::facet_grid(panel ~ strategy,
                                 labeller = ggplot2::labeller(
                                   strategy = ggplot2::label_both))
    grp <- c("strategy", "method", "scope", "feature")
  } else {
    facet <- ggplot2::facet_grid(method ~ strategy,
                                 labeller = ggplot2::label_both)
    grp <- c("strategy", "method", "feature")
  }
  ttl <- if (identical(scopes, "cumulative")) {
    "Dynamic Love plot: cumulative RRB"
  } else {
    "Dynamic Love plot: risk-set reconstruction balance"
  }
  p <- ggplot2::ggplot(x, ggplot2::aes(.data$time, .data$abs_rrb,
                                       colour = .data$feature,
                                       group = .data$feature)) +
    ggplot2::geom_line(data = .multi_point(x, grp), alpha = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    facet +
    ggplot2::labs(x = "Interval k", y = "|RRB(k)|",
                  colour = "Feature f(H_k)", title = ttl) +
    ggplot2::theme_bw()
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = 2,
                                 colour = "grey40")
  }
  p
}

#' Balance heatmap with IAI and MaxIB margins
#'
#' Shows the absolute balance value by feature and interval, one panel per
#' strategy, with two margin columns holding IAI and MaxIB computed by
#' [aggregate_balance()] over the displayed intervals.
#'
#' @param x Balance table with a single method and a single scope, e.g. the
#'   output of [rrb()] or [ks_rr()].
#' @param value Column to display in absolute value (`"rrb"` or `"ks"`).
#' @param censored_times_only See [plot_love_dynamic()].
#' @return A ggplot object.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' plot_balance_heatmap(rrb(w, c("X", "L", "L_lag")))
plot_balance_heatmap <- function(x, value = "rrb", censored_times_only = TRUE) {
  .check_cols(x, c("strategy", "time", "feature", value), "x")
  if (censored_times_only && "n_censored" %in% names(x)) {
    x <- x[x$n_censored > 0, , drop = FALSE]
  }
  if ("method" %in% names(x) && length(unique(x$method)) > 1) {
    stop("Supply a single method to plot_balance_heatmap().", call. = FALSE)
  }
  if ("scope" %in% names(x) && length(unique(x$scope)) > 1) {
    stop("Supply a single scope to plot_balance_heatmap().", call. = FALSE)
  }
  ttl <- if ("scope" %in% names(x) && identical(unique(x$scope),
                                                "cumulative")) {
    "Cumulative reconstruction balance heatmap"
  } else {
    "Reconstruction balance heatmap"
  }
  agg <- aggregate_balance(x, value = value, by = c("strategy", "feature"))
  times <- sort(unique(x$time))
  lv <- c(as.character(times), "IAI", "MaxIB")
  cells <- data.frame(strategy = x$strategy, feature = x$feature,
                      col = as.character(x$time), v = abs(x[[value]]),
                      margin = FALSE)
  marg <- rbind(
    data.frame(strategy = agg$strategy, feature = agg$feature, col = "IAI",
               v = agg$iai, margin = TRUE),
    data.frame(strategy = agg$strategy, feature = agg$feature, col = "MaxIB",
               v = agg$maxib, margin = TRUE)
  )
  df <- rbind(cells, marg)
  df$col <- factor(df$col, levels = lv)
  df$feature <- factor(df$feature, levels = rev(levels(factor(x$feature))))
  ggplot2::ggplot(df, ggplot2::aes(.data$col, .data$feature, fill = .data$v)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(data = df[df$margin, , drop = FALSE],
                       ggplot2::aes(label = sprintf("%.2f", .data$v)),
                       size = 2.4) +
    ggplot2::geom_vline(xintercept = length(times) + 0.5, linewidth = 0.8) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#cb181d",
                                 name = paste0("|", value, "|"),
                                 na.value = "grey85") +
    ggplot2::facet_wrap(~strategy, ncol = 1,
                        labeller = ggplot2::label_both) +
    ggplot2::labs(x = "Interval k  |  summaries",
                  y = "Feature f(H_k)",
                  title = ttl) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Weight-support plot
#'
#' Plots `ESS / n_risk`, the 95th and 99th weight percentiles and
#' `min(pi_hat)` against the interval, by strategy.
#'
#' @param x Output of [weight_support()].
#' @return A ggplot object.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' plot_weight_support(weight_support(w))
plot_weight_support <- function(x) {
  .check_cols(x, c("strategy", "time", "ess_ratio", "w95", "w99", "pi_min"),
              "x")
  long <- tidyr::pivot_longer(
    x[, c("strategy", "time", "ess_ratio", "w95", "w99", "pi_min")],
    cols = c("ess_ratio", "w95", "w99", "pi_min"),
    names_to = "metric", values_to = "value"
  )
  long$metric <- factor(long$metric,
                        levels = c("ess_ratio", "w95", "w99", "pi_min"),
                        labels = c("ESS / N", "w95", "w99", "min pi-hat"))
  ggplot2::ggplot(long, ggplot2::aes(.data$time, .data$value,
                                     colour = .data$strategy)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
    ggplot2::labs(x = "Interval k", y = NULL, colour = "Strategy",
                  title = "Weight support / positivity") +
    ggplot2::theme_bw()
}

#' ECDF before and after artificial censoring
#'
#' For one feature, strategy and selected intervals, shows the weighted ECDF
#' of the risk set before censoring (reference weights), the unweighted ECDF
#' of the retained rows and the weighted ECDF of the retained rows (post
#' weights), as in the interval scope of [rrb()].
#'
#' @param data Long data with weights (see [fit_censoring_model()]).
#' @param feature Feature column.
#' @param strategy_level Strategy to display.
#' @param times Intervals to display (default: the first four intervals
#'   with artificial censoring).
#' @param method Weighting method (see [rrb()]).
#' @param strategy,time,censored Column names.
#' @param ... Weight column names (`w_ref`, `w_post`, `pi_hat`, `q_hat`,
#'   `w_prev_unstab`, `w_prev_stab`), if different from the defaults.
#' @return A ggplot object.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' plot_ecdf_pre_post(w, "L", "g3", times = c(0, 4))
plot_ecdf_pre_post <- function(data, feature, strategy_level, times = NULL,
                               method = "weighted", strategy = "strategy",
                               time = "time", censored = "censored_next",
                               ...) {
  .check_cols(data, c(strategy, time, censored, feature))
  .check_methods(method)
  cols <- utils::modifyList(
    list(w_ref = "w_ref", w_post = "w_post", pi_hat = "pi_hat",
         q_hat = "q_hat", w_prev_unstab = "w_prev_unstab",
         w_prev_stab = "w_prev_stab"), list(...))
  d <- data[data[[strategy]] == strategy_level, , drop = FALSE]
  if (is.null(times)) {
    ct <- sort(unique(d[[time]][d[[censored]] == 1]))
    times <- utils::head(ct, 4)
  }
  d <- d[d[[time]] %in% times & !is.na(d[[feature]]), , drop = FALSE]
  d <- .rr_weights(d, method, censored, cols$w_ref, cols$w_post, cols$pi_hat,
                   cols$q_hat, cols$w_prev_unstab, cols$w_prev_stab)
  curves <- lapply(sort(unique(d[[time]])), function(tt) {
    s <- d[d[[time]] == tt, , drop = FALSE]
    x <- s[[feature]]
    ret <- s[[censored]] == 0
    at <- sort(unique(x))
    rbind(
      data.frame(time = tt, x = at, F = .wecdf(x, s$.ref_w, at),
                 set = "pre-censoring risk set"),
      data.frame(time = tt, x = at, F = .wecdf(x[ret], rep(1, sum(ret)), at),
                 set = "retained (unweighted)"),
      data.frame(time = tt, x = at, F = .wecdf(x[ret], s$.post_w[ret], at),
                 set = "retained (IPCW)")
    )
  })
  cdf <- do.call(rbind, curves)
  cdf$set <- factor(cdf$set, levels = c("pre-censoring risk set",
                                        "retained (unweighted)",
                                        "retained (IPCW)"))
  ggplot2::ggplot(cdf, ggplot2::aes(.data$x, .data$F, colour = .data$set,
                                    linetype = .data$set)) +
    ggplot2::geom_step() +
    ggplot2::facet_wrap(~time, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = c("black", "#E69F00", "#0072B2")) +
    ggplot2::scale_linetype_manual(values = c(1, 2, 1)) +
    ggplot2::labs(x = feature, y = "ECDF", colour = NULL, linetype = NULL,
                  title = sprintf("ECDF of %s, strategy %s: pre vs post",
                                  feature, strategy_level)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Mahalanobis reconstruction balance plot
#'
#' Plots `M^RR_g(k)` against the interval, by strategy and method (one
#' panel per scope when `x` has both scopes).
#'
#' @param x Output of [mahalanobis_rr()].
#' @param censored_times_only See [plot_love_dynamic()].
#' @return A ggplot object.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' plot_mahalanobis(mahalanobis_rr(w, c("X", "L"), c("unweighted", "weighted")))
plot_mahalanobis <- function(x, censored_times_only = TRUE) {
  .check_cols(x, c("strategy", "time", "method", "mahalanobis"), "x")
  if (censored_times_only && "n_censored" %in% names(x)) {
    x <- x[x$n_censored > 0, , drop = FALSE]
  }
  x <- x[!is.na(x$mahalanobis), , drop = FALSE]
  multi_scope <- "scope" %in% names(x) && length(unique(x$scope)) > 1
  grp <- c("strategy", "method", if ("scope" %in% names(x)) "scope")
  p <- ggplot2::ggplot(x, ggplot2::aes(.data$time, .data$mahalanobis,
                                       colour = .data$strategy,
                                       linetype = .data$method)) +
    ggplot2::geom_line(data = .multi_point(x, grp)) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(x = "Interval k", y = expression(M^RR * (k)),
                  colour = "Strategy", linetype = "Method",
                  title = "Mahalanobis reconstruction balance") +
    ggplot2::theme_bw()
  if (multi_scope) {
    p <- p + ggplot2::facet_wrap(~scope, scales = "free_y",
                                 labeller = ggplot2::label_both)
  }
  p
}

#' Diagnostic dashboard
#'
#' Combines the dynamic Love plot (unweighted and weighted interval RRB),
#' the heatmap of weighted interval RRB, the Mahalanobis plot, the
#' weight-support plot and an ECDF plot. The panels are arranged with
#' patchwork if installed, else with gridExtra; without either, the list of
#' plots is returned.
#'
#' @param data Long data with weights (see [fit_censoring_model()]).
#' @param features Feature columns for RRB and the Mahalanobis metric.
#' @param ecdf_feature Feature of the ECDF plot.
#' @param ecdf_strategy Strategy of the ECDF plot (default: last level).
#' @param ecdf_times Intervals of the ECDF plot.
#' @param title Optional title.
#' @param ridge Ridge constant passed to [mahalanobis_rr()].
#' @return A patchwork object, a gtable (gridExtra) or a named list of
#'   ggplot objects; the plots are also attached as attribute `"plots"`.
#' @export
#' @examples
#' \donttest{
#' d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
#' d <- make_history_features(d)
#' w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
#' db <- plot_dashboard(w, c("X", "L", "L_sq"))
#' }
plot_dashboard <- function(data, features, ecdf_feature = features[1],
                           ecdf_strategy = NULL, ecdf_times = NULL,
                           title = NULL, ridge = 0) {
  if (is.null(ecdf_strategy)) {
    s <- data$strategy
    ecdf_strategy <- if (is.factor(s)) utils::tail(levels(s), 1) else
      utils::tail(sort(unique(s)), 1)
  }
  r <- rrb(data, features, method = c("unweighted", "weighted"))
  plots <- list(
    love = plot_love_dynamic(r),
    heatmap = plot_balance_heatmap(r[r$method == "weighted", ]),
    mahalanobis = plot_mahalanobis(
      mahalanobis_rr(data, features, method = c("unweighted", "weighted"),
                     ridge = ridge)),
    support = plot_weight_support(weight_support(data)),
    ecdf = plot_ecdf_pre_post(data, ecdf_feature, ecdf_strategy,
                              times = ecdf_times)
  )
  if (requireNamespace("patchwork", quietly = TRUE)) {
    out <- (plots$love | plots$heatmap) / (plots$mahalanobis | plots$ecdf) /
      plots$support
    out <- out + patchwork::plot_annotation(title = title)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    out <- gridExtra::arrangeGrob(
      plots$love, plots$heatmap, plots$mahalanobis, plots$ecdf,
      plots$support,
      layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)), top = title
    )
  } else {
    out <- plots
  }
  attr(out, "plots") <- plots
  out
}
