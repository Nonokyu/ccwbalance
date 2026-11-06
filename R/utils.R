# Internal helpers ----------------------------------------------------------

.check_cols <- function(data, cols, arg = "data") {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf("`%s` is missing required column(s): %s", arg,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.wmean <- function(x, w) {
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(NA_real_)
  sum(w * x) / sw
}

.wvar <- function(x, w) {
  m <- .wmean(x, w)
  if (is.na(m)) return(NA_real_)
  sum(w * (x - m)^2) / sum(w)
}

# Weighted ECDF evaluated at `at`
.wecdf <- function(x, w, at) {
  o <- order(x)
  xs <- x[o]
  cw <- cumsum(w[o]) / sum(w)
  idx <- findInterval(at, xs)
  c(0, cw)[idx + 1]
}

# Weighted two-sample KS distance
.wks <- function(x1, w1, x2, w2) {
  if (length(x1) == 0 || length(x2) == 0 || sum(w1) <= 0 || sum(w2) <= 0) {
    return(NA_real_)
  }
  at <- sort(unique(c(x1, x2)))
  max(abs(.wecdf(x1, w1, at) - .wecdf(x2, w2, at)))
}

# Build reference-side and post-censoring weights for a reconstruction
# comparison. Returns data with `.ref_w` (weight of every at-risk row in the
# pre-censoring reference) and `.post_w` (weight of retained rows; 0 for rows
# censored in [t, t+1)).
.rr_weights <- function(data, method, censored, w_ref, w_post, pi_hat, q_hat,
                        w_prev_unstab, w_prev_stab) {
  retained <- as.numeric(data[[censored]] == 0)
  if (method == "unweighted") {
    .check_cols(data, censored)
    data$.ref_w <- 1
    data$.post_w <- retained
  } else if (method == "weighted") {
    .check_cols(data, c(censored, w_ref, w_post))
    data$.ref_w <- data[[w_ref]]
    data$.post_w <- retained * data[[w_post]]
  } else if (method == "unstabilized") {
    .check_cols(data, c(censored, pi_hat, w_prev_unstab))
    data$.ref_w <- data[[w_prev_unstab]]
    data$.post_w <- retained * data[[w_prev_unstab]] / data[[pi_hat]]
  } else if (method == "stabilized") {
    .check_cols(data, c(censored, pi_hat, q_hat, w_prev_stab))
    data$.ref_w <- data[[w_prev_stab]] * data[[q_hat]]
    data$.post_w <- retained * data[[w_prev_stab]] * data[[q_hat]] /
      data[[pi_hat]]
  } else {
    stop("Unknown method: ", method, call. = FALSE)
  }
  data
}

.rr_methods <- c("weighted", "unweighted", "unstabilized", "stabilized")

# Apply a summary function, returning NA for an empty vector
.safe <- function(x, f, ...) if (length(x) == 0) NA_real_ else f(x, ...)
