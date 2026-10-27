#' Parameters of the simulation data-generating mechanism
#'
#' Returns the parameter list of [simulate_ccw()] and [true_risk_ccw()].
#' Any element can be overridden through `...`; the coefficient vectors
#' `alpha`, `beta`, `gamma` and `delta` can be overridden partially by
#' name, e.g. `ccw_sim_params(beta = c(b2 = 0.4))`.
#'
#' In intervals `k = 0, ..., n_time - 1`:
#' \deqn{L_k = \alpha_0 + \alpha_1 L_{k-1} + \alpha_2 A_{k-1} + \alpha_3 X +
#'   \epsilon_k,}
#' \deqn{logit P(A_k = 1) = \beta_0 + \beta_1 L_k + \beta_2 L_k^2 +
#'   \beta_3 X + \beta_4 A_{k-1} + \beta_5 L_k X + \beta_6 L_{k-1} +
#'   \beta_{extra}' X_{extra},}
#' \deqn{logit P(C_k = 1) = \delta_0 + \delta_1 L_k,}
#' \deqn{logit P(Y_k = 1) = \gamma_0 + \gamma_1 L_k + \gamma_2 A_k +
#'   \gamma_3 X + \gamma_4 L_k^2 + \gamma_5 L_k X + \gamma_6 L_{k-1} +
#'   \gamma_{extra}' X_{extra},}
#' with `X ~ N(0, 1)`, `epsilon_k ~ N(0, sigma_L^2)`, `L_{-1} = 0`,
#' `A_{-1} = 0` and `C_k` loss to follow-up (`delta0 = -Inf`, the default,
#' gives none). With `absorbing = TRUE` (default) treatment stays on once
#' initiated and `beta4` is unused. The terms `gamma4` to `gamma6` (default
#' 0) let terms omitted from a censoring model also affect the outcome.
#'
#' Extra baseline covariates (`n_extra > 0`, default 0):
#' `X_extra = rho_extra X + sqrt(1 - rho_extra^2) Z`, `Z ~ N(0, I)`, so that
#' each is standard normal with correlation `rho_extra` with `X`. They enter
#' the treatment model with coefficients `beta_extra` (default
#' `rep(c(0.2, -0.2), length.out = n_extra)`) and the outcome model with
#' `gamma_extra` (default `rep(0.2, n_extra)`). With `n_extra = 0` no
#' additional random numbers are drawn.
#'
#' The strategies are defined by `window = c(w1, w2)`:
#' * `g1`: initiate in intervals `0, ..., w1 - 1`;
#' * `g2`: initiate in intervals `w1, ..., w2 - 1`;
#' * `g3`: do not initiate in intervals `0, ..., w2 - 1`.
#'
#' @param ... Named overrides.
#' @return A named list of parameters.
#' @export
#' @examples
#' str(ccw_sim_params(n = 500))
ccw_sim_params <- function(...) {
  defaults <- list(
    n = 2000,
    n_time = 12,
    window = c(3, 9),
    alpha = c(a0 = 0, a1 = 0.5, a2 = -0.5, a3 = 0.3),
    sigma_L = 0.6,
    beta = c(b0 = -2.2, b1 = 0.4, b2 = 0.1, b3 = 0.3, b4 = 0, b5 = 0.15,
             b6 = 0),
    gamma = c(g0 = -4, g1 = 0.5, g2 = -0.5, g3 = 0.3, g4 = 0, g5 = 0, g6 = 0),
    delta = c(d0 = -Inf, d1 = 0),
    absorbing = TRUE,
    n_extra = 0,
    rho_extra = 0.5,
    beta_extra = NULL,
    gamma_extra = NULL
  )
  args <- list(...)
  vec_pars <- c("alpha", "beta", "gamma", "delta")
  for (nm in intersect(names(args), vec_pars)) {
    v <- args[[nm]]
    if (is.null(names(v)) || !all(names(v) %in% names(defaults[[nm]]))) {
      stop(sprintf("`%s` must be a named vector with names in: %s", nm,
                   paste(names(defaults[[nm]]), collapse = ", ")),
           call. = FALSE)
    }
    defaults[[nm]][names(v)] <- v
    args[[nm]] <- NULL
  }
  p <- utils::modifyList(defaults, args)
  .resolve_extra(p)
}

# Resolve and validate the extra-covariate parameters
.resolve_extra <- function(p) {
  k <- p$n_extra
  if (length(k) != 1 || is.na(k) || k < 0 || k != round(k)) {
    stop("`n_extra` must be a single non-negative integer.", call. = FALSE)
  }
  p$n_extra <- as.integer(k)
  if (length(p$rho_extra) != 1 || abs(p$rho_extra) > 1) {
    stop("`rho_extra` must be a single number in [-1, 1].", call. = FALSE)
  }
  if (is.null(p$beta_extra) && k > 0) {
    p$beta_extra <- rep(c(0.2, -0.2), length.out = k)
  }
  if (is.null(p$gamma_extra) && k > 0) {
    p$gamma_extra <- rep(0.2, k)
  }
  for (nm in c("beta_extra", "gamma_extra")) {
    if (length(p[[nm]]) != k) {
      stop(sprintf("`%s` must have length n_extra = %d (got %d).", nm, k,
                   length(p[[nm]])), call. = FALSE)
    }
  }
  p
}

# Person-period matrices. `regime = NULL` gives the natural course; "g1",
# "g2" or "g3" sets treatment by the strategy (no loss to follow-up).
.sim_paths <- function(p, regime = NULL) {
  n <- p$n
  K <- p$n_time
  a <- p$alpha
  b <- p$beta
  gm <- p$gamma
  d <- p$delta
  w1 <- p$window[1]
  w2 <- p$window[2]

  X <- stats::rnorm(n)
  k_extra <- if (is.null(p$n_extra)) 0L else p$n_extra
  lp_extra <- y_extra <- 0
  Xe <- NULL
  if (k_extra > 0) {
    rho <- p$rho_extra
    Z <- matrix(stats::rnorm(n * k_extra), n, k_extra)
    Xe <- rho * X + sqrt(1 - rho^2) * Z
    lp_extra <- as.vector(Xe %*% p$beta_extra)
    y_extra <- as.vector(Xe %*% p$gamma_extra)
  }
  Lm <- Am <- Pm <- Ym <- Cm <- matrix(NA_real_, n, K)
  L_prev <- rep(0, n)
  A_prev <- rep(0, n)
  for (t in seq_len(K) - 1) {
    eps <- stats::rnorm(n, 0, p$sigma_L)
    L <- a[["a0"]] + a[["a1"]] * L_prev + a[["a2"]] * A_prev + a[["a3"]] * X +
      eps
    lp <- b[["b0"]] + b[["b1"]] * L + b[["b2"]] * L^2 + b[["b3"]] * X +
      b[["b5"]] * L * X + b[["b6"]] * L_prev + lp_extra
    if (p$absorbing) {
      p_init <- ifelse(A_prev == 1, 1, stats::plogis(lp))
    } else {
      p_init <- stats::plogis(lp + b[["b4"]] * A_prev)
    }
    A <- stats::rbinom(n, 1, p_init)
    if (!is.null(regime)) {
      ever <- A_prev == 1
      A <- switch(regime,
        g1 = if (t < w1 - 1) A else if (t == w1 - 1) rep(1, n) else
          if (p$absorbing) rep(1, n) else A,
        g2 = if (t < w1) rep(0, n) else if (t < w2 - 1) A else
          if (t == w2 - 1) rep(1, n) else if (p$absorbing) rep(1, n) else A,
        g3 = if (t <= w2 - 1) rep(0, n) else A,
        stop("Unknown regime: ", regime, call. = FALSE)
      )
      if (regime %in% c("g1", "g2") && p$absorbing) A[ever] <- 1
    }
    ltfu <- if (is.null(regime)) {
      stats::rbinom(n, 1, stats::plogis(d[["d0"]] + d[["d1"]] * L))
    } else {
      rep(0, n)
    }
    ly <- gm[["g0"]] + gm[["g1"]] * L + gm[["g2"]] * A + gm[["g3"]] * X +
      gm[["g4"]] * L^2 + gm[["g5"]] * L * X + gm[["g6"]] * L_prev +
      y_extra
    y <- stats::rbinom(n, 1, stats::plogis(ly))
    Lm[, t + 1] <- L
    Am[, t + 1] <- A
    Pm[, t + 1] <- L_prev
    Ym[, t + 1] <- y
    Cm[, t + 1] <- ltfu
    L_prev <- L
    A_prev <- A
  }
  list(X = X, L = Lm, A = Am, L_lag = Pm, Y = Ym, C = Cm, X_extra = Xe)
}

.first_one <- function(m) {
  apply(m, 1, function(r) {
    w <- which(r == 1)
    if (length(w) == 0) Inf else w[1] - 1
  })
}

#' Simulate long data from the data-generating mechanism
#'
#' Simulates observational data from the mechanism of [ccw_sim_params()],
#' evaluates every individual under the three strategies and applies
#' artificial censoring in the interval of the first deviation. Returns long
#' data in the package format (see [ccwbalance-package]).
#'
#' Within interval `k` the order is `L_k`, `A_k`, artificial censoring
#' (`D_k`), loss to follow-up, outcome. A row with `censored_next = 1` has
#' `event_next = NA`.
#'
#' @param params Parameter list from [ccw_sim_params()].
#' @param seed Optional random seed.
#' @return A data frame with columns `id, strategy, time, at_risk, X, L,
#'   L_lag, A, A_lag, censor_possible, phase, censored_next, ltfu_next,
#'   event_next` (with `n_extra > 0`, the columns `X2, ..., X{n_extra + 1}`
#'   follow `X`).
#'   `phase` is `"no_init"` (censored if initiating), `"grace_end"`
#'   (censored if not initiating) or `"none"`. The parameters are attached
#'   as attribute `"params"`.
#' @export
#' @examples
#' d <- simulate_ccw(ccw_sim_params(n = 200), seed = 1)
#' head(d)
simulate_ccw <- function(params = ccw_sim_params(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  K <- p$n_time
  w1 <- p$window[1]
  w2 <- p$window[2]
  if (!(w1 >= 1 && w2 > w1 && w2 <= K)) {
    stop("`window` must satisfy 1 <= w1 < w2 <= n_time.", call. = FALSE)
  }
  s <- .sim_paths(p)
  n <- p$n

  tau <- .first_one(s$A) # initiation interval
  ev <- .first_one(s$Y) # event interval
  lt <- .first_one(s$C) # LTFU interval

  cens <- list(
    g1 = ifelse(tau > w1 - 1, w1 - 1, Inf),
    g2 = ifelse(tau < w1, tau, ifelse(tau > w2 - 1, w2 - 1, Inf)),
    g3 = ifelse(tau <= w2 - 1, tau, Inf)
  )

  times <- seq_len(K) - 1
  base <- data.frame(
    id = rep(seq_len(n), each = K),
    time = rep(times, times = n),
    X = rep(s$X, each = K),
    L = as.vector(t(s$L)),
    L_lag = as.vector(t(s$L_lag)),
    A = as.vector(t(s$A)),
    Y = as.vector(t(s$Y)),
    C = as.vector(t(s$C))
  )
  extra_names <- character(0)
  if (!is.null(s$X_extra)) {
    extra_names <- paste0("X", seq_len(ncol(s$X_extra)) + 1)
    for (j in seq_along(extra_names)) {
      base[[extra_names[j]]] <- rep(s$X_extra[, j], each = K)
    }
  }
  base$A_lag <- stats::ave(base$A, base$id,
                           FUN = function(a) c(0, a[-length(a)]))

  out <- lapply(names(cens), function(g) {
    ct <- cens[[g]][base$id]
    last <- pmin(ct, ev[base$id], lt[base$id], K - 1)
    d <- base[base$time <= last, , drop = FALSE]
    ctd <- ct[base$time <= last]
    d$strategy <- g
    d$censored_next <- as.integer(d$time == ctd)
    not_init <- d$A_lag == 0
    d$phase <- switch(g,
      g1 = ifelse(d$time == w1 - 1 & not_init, "grace_end", "none"),
      g2 = ifelse(d$time < w1 & not_init, "no_init",
                  ifelse(d$time == w2 - 1 & not_init, "grace_end", "none")),
      g3 = ifelse(d$time <= w2 - 1 & not_init, "no_init", "none")
    )
    d$censor_possible <- d$phase != "none"
    d$ltfu_next <- ifelse(d$censored_next == 1, 0L, as.integer(d$C))
    d$event_next <- ifelse(d$censored_next == 1 | d$ltfu_next == 1,
                           NA_integer_, as.integer(d$Y))
    d
  })
  out <- do.call(rbind, out)
  out$at_risk <- 1L
  out$strategy <- factor(out$strategy, levels = c("g1", "g2", "g3"))
  out <- out[order(out$strategy, out$id, out$time),
             c("id", "strategy", "time", "at_risk", "X", extra_names, "L",
               "L_lag", "A", "A_lag", "censor_possible", "phase",
               "censored_next", "ltfu_next", "event_next")]
  rownames(out) <- NULL
  attr(out, "params") <- p
  out
}

#' True risk under each simulated strategy
#'
#' Monte Carlo evaluation of the risk of the outcome by the end of interval
#' `n_time - 1` under the regimes identified by the weighted analysis of
#' [simulate_ccw()] data, without loss to follow-up: natural initiation
#' within the window with initiation forced in its last interval (`g1`;
#' `g2`, which also prevents initiation before `w1`), or no initiation up
#' to interval `w2 - 1` (`g3`).
#'
#' @param params Parameter list from [ccw_sim_params()]; `n` is replaced by
#'   `n_mc`.
#' @param n_mc Monte Carlo sample size.
#' @param seed Optional random seed.
#' @return A data frame with `strategy` and `risk`.
#' @export
#' @examples
#' true_risk_ccw(ccw_sim_params(), n_mc = 2000, seed = 1)
true_risk_ccw <- function(params = ccw_sim_params(), n_mc = 2e5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  p$n <- n_mc
  risk <- vapply(c("g1", "g2", "g3"), function(g) {
    s <- .sim_paths(p, regime = g)
    mean(rowSums(s$Y) > 0)
  }, numeric(1))
  data.frame(strategy = factor(names(risk), levels = c("g1", "g2", "g3")),
             risk = unname(risk))
}
