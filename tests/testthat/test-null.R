null_data <- function(formula = ~ L + I(L^2) + X + L:X) {
  d <- simulate_ccw(ccw_sim_params(n = 1000), seed = 3)
  fit_censoring_model(d, formula)
}
null_features <- c("X", "L", "L_lag")

test_that("rrb_null returns the documented structure", {
  w <- small_weighted(n = 300)
  nr <- rrb_null(w, c("X", "L"), n_draws = 20, seed = 1)
  expect_s3_class(nr, "rrb_null")
  expect_equal(nr$n_draws, 20)
  expect_true(all(c("strategy", "time", "feature", "n_risk", "n_censored",
                    "rrb_obs", "null_sd", "null_q025", "null_q975", "z",
                    "p_two") %in% names(nr$cells)))
  expect_true(all(nr$cells$n_censored >= 1))
  expect_true(all(nr$cells$p_two > 0 & nr$cells$p_two <= 1))
  expect_equal(nrow(nr$draws), 20 * nrow(nr$cells))
  expect_true(all(c("strategy", "feature", "n_cells", "maxib_obs",
                    "null_mean", "null_q95", "p") %in% names(nr$maxib)))
  # per strategy x feature, per strategy, per feature and overall
  n_s <- length(unique(nr$cells$strategy))
  expect_equal(nrow(nr$maxib), n_s * 2 + n_s + 2 + 1)
  all_row <- nr$maxib$strategy == "(all)" & nr$maxib$feature == "(all)"
  expect_equal(nr$maxib$maxib_obs[all_row], max(abs(nr$cells$rrb_obs)))
  # observed values agree with rrb() (untruncated weights)
  ref <- rrb(w, c("X", "L"))
  ref <- ref[ref$n_censored >= 1 & !is.na(ref$rrb), ]
  expect_equal(nr$cells$rrb_obs, ref$rrb)
  expect_output(print(nr), "MaxIB")
})

test_that("internal refit and RRB engine reproduce the observed fit", {
  w <- small_weighted(n = 300)
  wt <- attr(w, "weighting")
  expect_equal(wt$by, c("strategy", "phase"))
  gr <- .null_groups(w, wt$formula, wt$by, "censor_possible")
  p <- pmax(.null_refit(gr, 1 - w$censored_next, nrow(w)), 1e-6)
  expect_equal(p, w$pi_hat)
  prep <- .null_prep(w, "L", "strategy", "time", "w_ref")
  r <- .null_rrb(prep, w$censored_next, w$pi_hat)
  ref <- rrb(w, "L")
  k <- match(paste(ref$strategy, ref$time, sep = "\r"), prep$cells)
  ok <- ref$n_censored >= 1
  expect_equal(unname(r[k[ok], 1]), ref$rrb[ok])
})

test_that("observed RRB is not extreme under a correct censoring model", {
  nr <- rrb_null(null_data(), null_features, n_draws = 50, seed = 2)
  per_feature <- nr$maxib[nr$maxib$strategy == "(all)" &
                            nr$maxib$feature != "(all)", ]
  expect_gte(sum(per_feature$p >= 0.05), 2)
  expect_lt(abs(mean(nr$cells$z, na.rm = TRUE)), 0.5)
})

test_that("omitting L from the censoring model is flagged for L", {
  nr <- rrb_null(null_data(~ X), null_features, n_draws = 50, seed = 2)
  m <- nr$maxib
  # L-driven censoring accumulates over the no-initiation phase of g3
  expect_lt(m$p[m$strategy == "g3" & m$feature == "L"], 0.05)
})

test_that("rrb_null is reproducible with a seed", {
  w <- small_weighted(n = 300)
  a <- rrb_null(w, "X", n_draws = 10, seed = 7)
  b <- rrb_null(w, "X", n_draws = 10, seed = 7)
  expect_identical(a$draws, b$draws)
  expect_identical(a$maxib, b$maxib)
})

test_that("coupled redraw keeps each row's marginal censoring probability", {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
  w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
  poss <- as.logical(w$censor_possible)
  key <- paste(w$id[poss], w$time[poss])
  grp <- match(key, unique(key))
  expect_lt(max(grp), sum(poss)) # some id x time groups have several rows
  ge <- w$phase[poss] == "grace_end"
  pi <- w$pi_hat[poss]
  set.seed(11)
  dm <- vapply(seq_len(400), function(b) {
    .null_draw(pi, ge, grp, stats::runif(max(grp)))
  }, integer(sum(poss)))
  m <- rowMeans(dm)
  expect_gt(stats::cor(m, 1 - pi), 0.95)
  expect_lt(mean(abs(m - (1 - pi))), 0.08)
})

test_that("coupled redraw ties grace-end and no-initiation clones", {
  u <- seq(0.005, 0.995, by = 0.01)
  n_u <- length(u)
  grp <- rep(seq_len(n_u), each = 2)
  ge <- rep(c(TRUE, FALSE), n_u)
  # mirror-image models: pi_ge = 1 - pi_ni -> exactly complementary
  d <- .null_draw(rep(c(0.3, 0.7), n_u), ge, grp, u)
  expect_true(all(d[ge] + d[!ge] == 1))
  expect_equal(mean(d[ge]), 0.7, tolerance = 0.02)
  # pi_ge >= 1 - pi_ni: D*_ni = 1 implies D*_ge = 0 (never both censored)
  d <- .null_draw(rep(c(0.6, 0.8), n_u), ge, grp, u)
  expect_true(all(d[ge][d[!ge] == 1] == 0))
  expect_true(any(d[ge] + d[!ge] == 0))
  # pi_ge < 1 - pi_ni: D*_ge = 0 implies D*_ni = 1 (at least one censored)
  d <- .null_draw(rep(c(0.2, 0.5), n_u), ge, grp, u)
  expect_true(all(d[!ge][d[ge] == 0] == 1))
  expect_true(any(d[ge] + d[!ge] == 2))
  # two no-initiation clones are censored by the same event
  d <- .null_draw(rep(0.6, 2 * n_u), rep(FALSE, 2 * n_u), grp, u)
  expect_identical(d[c(TRUE, FALSE)], d[c(FALSE, TRUE)])
  # a singleton group is a plain Bernoulli(1 - pi) draw of its own U
  expect_identical(.null_draw(0.4, FALSE, 1L, 0.5), 1L)
  expect_identical(.null_draw(0.4, TRUE, 1L, 0.5), 1L)
  expect_identical(.null_draw(0.4, TRUE, 1L, 0.3), 0L)
})

test_that("coupling = 'independent' reproduces the earlier draws", {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
  w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
  nr <- rrb_null(w, c("X", "L"), n_draws = 20, seed = 1,
                 coupling = "independent")
  expect_equal(nr$maxib$null_q95,
               c(0.216874, 0.269155, 0.178523, 0.088963, 0.131392, 0.177494,
                 0.216874, 0.269155, 0.191723, 0.269155, 0.190505, 0.269155),
               tolerance = 1e-5)
  expect_equal(head(nr$cells$null_sd, 5),
               c(0.0898459, 0.0640627, 0.1060405, 0.0622723, 0.0194842),
               tolerance = 1e-5)
})

test_that("rrb_null default coupling runs and checks the id column", {
  w <- small_weighted(n = 300)
  nr <- rrb_null(w, c("X", "L"), n_draws = 10, seed = 1)
  ni <- rrb_null(w, c("X", "L"), n_draws = 10, seed = 1,
                 coupling = "independent")
  expect_equal(nr$cells$rrb_obs, ni$cells$rrb_obs)
  expect_false(identical(nr$draws$rrb, ni$draws$rrb))
  expect_output(print(nr), "MaxIB")
  w2 <- w
  names(w2)[names(w2) == "id"] <- "pid"
  expect_error(rrb_null(w2, "X", n_draws = 2), "missing required column")
  expect_s3_class(rrb_null(w2, "X", n_draws = 2, id = "pid"), "rrb_null")
  expect_s3_class(rrb_null(w2, "X", n_draws = 2, coupling = "independent"),
                  "rrb_null")
})

test_that("window-mode orientation from A equals the phase rule", {
  w <- small_weighted(n = 300)
  on0 <- .null_orient(w, "strategy", "censored_next", "censor_possible", "A")
  poss <- w$censor_possible
  expect_identical(on0[poss], w$phase[poss] == "grace_end")
  # without a treatment column the phase rule is used
  w2 <- w[, setdiff(names(w), "A")]
  expect_identical(
    .null_orient(w2, "strategy", "censored_next", "censor_possible", "A"),
    w$phase == "grace_end"
  )
})

test_that("sustained always / never clones are coupled as complements", {
  set.seed(5)
  n <- 300
  pp <- data.frame(id = rep(seq_len(n), each = 4), t = rep(0:3, n),
                   L = stats::rnorm(4 * n))
  pp$a <- stats::rbinom(nrow(pp), 1, stats::plogis(0.8 * pp$L))
  pp$y <- stats::rbinom(nrow(pp), 1, 0.03)
  d <- ccw_clone(pp, list(always = ccw_sustained(1),
                          never = ccw_sustained(0)),
                 time = "t", treatment = "a", event = "y", covariates = "L")
  w <- fit_censoring_model(d, ~ L)
  on0 <- .null_orient(w, "strategy", "censored_next", "censor_possible", "A")
  poss <- w$censor_possible
  expect_true(all(on0[poss & w$strategy == "always"]))
  expect_false(any(on0[poss & w$strategy == "never"]))
  # pairs of always / never clones censorable at the same id x time
  key <- paste(w$id, w$time)
  ia <- which(poss & w$strategy == "always")
  inv <- which(poss & w$strategy == "never")
  inv <- inv[match(key[ia], key[inv])]
  ok <- !is.na(inv)
  ia <- ia[ok]
  inv <- inv[ok]
  expect_gt(length(ia), 100)
  u <- seq(0.005, 0.995, by = 0.01)
  for (k in seq_len(5)) {
    pa <- w$pi_hat[ia[k]]
    # mirror probabilities: exactly complementary
    d1 <- .null_draw(rep(c(pa, 1 - pa), length(u)),
                     rep(c(on0[ia[k]], on0[inv[k]]), length(u)),
                     rep(seq_along(u), each = 2), u)
    expect_true(all(d1[c(TRUE, FALSE)] + d1[c(FALSE, TRUE)] == 1))
    # fitted probabilities: counter-monotone, never both of the ordered
    # outcomes that would indicate comonotone coupling
    pn <- w$pi_hat[inv[k]]
    d2 <- .null_draw(rep(c(pa, pn), length(u)),
                     rep(c(on0[ia[k]], on0[inv[k]]), length(u)),
                     rep(seq_along(u), each = 2), u)
    da <- d2[c(TRUE, FALSE)]
    dn <- d2[c(FALSE, TRUE)]
    if (pa >= 1 - pn) expect_true(all(da + dn <= 1))
    if (pa <= 1 - pn) expect_true(all(da + dn >= 1))
  }
  nr <- rrb_null(w, "L", n_draws = 5, seed = 1)
  expect_s3_class(nr, "rrb_null")
})

test_that("non-orientable censoring groups are rejected", {
  toy <- data.frame(strategy = "g", phase = "no_init",
                    censor_possible = TRUE, censored_next = c(1, 1, 0),
                    A = c(1, 0, 0))
  expect_error(.null_orient(toy, "strategy", "censored_next",
                            "censor_possible", "A"),
               "not a function of `A`")
  w <- small_weighted(n = 300)
  k <- which(w$censored_next == 1 & w$phase == "no_init" &
               w$strategy == "g3")[1]
  w$A[k] <- 1 - w$A[k]
  expect_error(rrb_null(w, "X", n_draws = 2), "coupling = 'independent'")
  expect_s3_class(rrb_null(w, "X", n_draws = 2, coupling = "independent"),
                  "rrb_null")
})

# Censoring balance in the null reference ------------------------------------
cb_fixture <- function() {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
}

test_that("adding censoring balance leaves the RRB output unchanged", {
  # stored rrb_null() output computed without the censoring-balance columns
  ref <- readRDS(test_path("ref-rrb-null.rds"))
  w <- cb_fixture()
  f <- c("X", "L", "L_lag")
  for (cp in c("individual", "independent")) {
    nr <- rrb_null(w, f, n_draws = 30, seed = 4, coupling = cp)
    expect_identical(nr$cells[, c("rrb_obs", "null_mean", "null_sd", "z")],
                     ref[[cp]]$cells)
    expect_identical(nr$maxib, ref[[cp]]$maxib)
    expect_identical(nr$draws$rrb, ref[[cp]]$rrb)
  }
})

test_that("observed CB equals censoring_balance(weights = 'model')", {
  w <- cb_fixture()
  f <- c("X", "L", "L_lag")
  nr <- rrb_null(w, f, n_draws = 5, seed = 1)
  expect_true(all(c("cb_obs", "cb_null_mean", "cb_null_sd", "cb_z",
                    "cb_p_two") %in% names(nr$cells)))
  cb <- censoring_balance(w, f, weighted = TRUE, weights = "model")
  k <- match(paste(nr$cells$strategy, nr$cells$time, nr$cells$feature),
             paste(cb$strategy, cb$time, cb$feature))
  ok <- !is.na(k)
  expect_gt(sum(ok), 10)
  expect_equal(nr$cells$cb_obs[ok], cb$smd[k[ok]], tolerance = 1e-8)
  # cells without a matching censoring_balance row have undefined CB
  expect_true(all(is.na(nr$cells$cb_obs[!ok])))
  expect_true("cb" %in% names(nr$draws))
  expect_output(print(nr), "Censoring balance")
  expect_output(print(nr), "Interval RRB")
})

test_that("CB null runs under both couplings and has a positive spread", {
  w <- cb_fixture()
  for (cp in c("individual", "independent")) {
    nr <- rrb_null(w, c("X", "L"), n_draws = 20, seed = 2, coupling = cp)
    expect_identical(nr$maxib_cb[, c("strategy", "feature")],
                     nr$maxib[, c("strategy", "feature")])
    expect_true(all(c("n_cells", "maxib_obs", "null_mean", "null_q95",
                      "p") %in% names(nr$maxib_cb)))
    all_row <- nr$maxib_cb$strategy == "(all)" &
      nr$maxib_cb$feature == "(all)"
    expect_equal(nr$maxib_cb$maxib_obs[all_row],
                 max(abs(nr$cells$cb_obs), na.rm = TRUE))
    has_cb <- !is.na(nr$cells$cb_obs)
    expect_true(all(nr$cells$cb_null_sd[has_cb] > 0, na.rm = TRUE))
    expect_gt(mean(!is.na(nr$cells$cb_null_sd[has_cb])), 0.8)
    p <- nr$cells$cb_p_two[has_cb]
    expect_true(all(p > 0 & p <= 1))
  }
})

test_that("internal CB engine is NA without censored rows", {
  w <- cb_fixture()
  prep <- .null_prep(w, "L", "strategy", "time", "w_ref")
  r <- .null_cb(prep, integer(nrow(w)), w$pi_hat)
  expect_true(all(is.na(r)))
})
