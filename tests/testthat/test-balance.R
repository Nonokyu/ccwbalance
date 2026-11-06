test_that("IPW identity holds exactly on a toy example", {
  toy <- toy_ipw()
  # Horvitz-Thompson form
  ht <- sum((toy$censored_next == 0) / toy$pi_hat * toy$f)
  expect_equal(ht, sum(toy$f))
  r <- rrb(toy, "f", method = c("unstabilized", "stabilized"))
  expect_equal(r$mean_pre, rep(mean(toy$f), 2))
  expect_equal(r$mean_post, rep(mean(toy$f), 2))
  expect_equal(r$rrb, c(0, 0))
  # unweighted selection is not balanced
  ru <- rrb(toy, "f", method = "unweighted")
  expect_gt(abs(ru$rrb), 0.1)
  # KS and Mahalanobis are 0 as well
  expect_equal(ks_rr(toy, "f", method = "unstabilized")$ks, 0)
  expect_equal(mahalanobis_rr(toy, "f", method = "unstabilized")$mahalanobis,
               0)
})

test_that("RRB is near 0 under a correctly specified censoring model", {
  set.seed(2024)
  d <- simulate_ccw(ccw_sim_params(n = 5000), seed = 2024)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  feats <- c("X", "L", "L_sq", "L_x_X")
  w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
  r <- rrb(w, feats, method = c("unweighted", "weighted"))
  r <- r[r$n_censored > 0, ]
  max_w <- max(abs(r$rrb[r$method == "weighted"]), na.rm = TRUE)
  max_u <- max(abs(r$rrb[r$method == "unweighted"]), na.rm = TRUE)
  expect_lt(max_w, 0.15)
  expect_gt(max_u, 0.2)
  # omitting L from the censoring model leaves more residual imbalance in L
  # (compared on the time-integrated IAI, which is less noisy than MaxIB)
  wm <- fit_censoring_model(d, ~ X)
  rm <- rrb(wm, feats)
  rm <- rm[rm$n_censored > 0, ]
  iai_ok <- aggregate_balance(r[r$method == "weighted", ], by = "feature")
  iai_bad <- aggregate_balance(rm, by = "feature")
  expect_gt(iai_bad$iai[iai_bad$feature == "L"],
            1.5 * iai_ok$iai[iai_ok$feature == "L"])
})

test_that("rrb validates inputs", {
  toy <- toy_ipw()
  expect_error(rrb(toy, "f", method = "nope"), "Unknown method")
  expect_error(rrb(toy, "missing_col", method = "unweighted"), "missing")
})

test_that("ks_rr matches a hand-computed weighted KS", {
  toy <- data.frame(strategy = "g", time = 0, f = c(1, 2, 3, 4),
                    censored_next = c(0, 0, 1, 1))
  k <- ks_rr(toy, "f", method = "unweighted")
  # pre ECDF at 2 = 0.5, post ECDF at 2 = 1
  expect_equal(k$ks, 0.5)
})

test_that("mahalanobis_rr equals squared RRB in one dimension", {
  w <- small_weighted()
  r <- rrb(w, "L")
  m <- mahalanobis_rr(w, "L")
  mm <- merge(r, m, by = c("strategy", "time", "method"))
  mm <- mm[mm$n_censored.x > 0, ]
  expect_equal(mm$mahalanobis, mm$rrb^2, tolerance = 1e-8)
  # ridge shrinks the metric
  m2 <- mahalanobis_rr(w, c("X", "L", "L_sq"), ridge = 1)
  m0 <- mahalanobis_rr(w, c("X", "L", "L_sq"), ridge = 0)
  ok <- !is.na(m0$mahalanobis) & m0$mahalanobis > 0
  expect_true(all(m2$mahalanobis[ok] <= m0$mahalanobis[ok] + 1e-12))
})

test_that("strategy_balance is exactly 0 at baseline for unweighted clones", {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 3)
  sb <- strategy_balance(d, c("X", "L"), weighted = FALSE)
  expect_true(all(abs(sb$smd[sb$time == 0]) < 1e-12))
  expect_equal(sort(unique(as.character(sb$strategy))), c("g2", "g3"))
  w <- small_weighted()
  sbw <- strategy_balance(w, "X", reference = "g3")
  expect_true(all(sbw$reference == "g3"))
})

test_that("censoring_balance compares censored vs retained", {
  toy <- data.frame(strategy = "g", time = 0, f = c(0, 0, 2, 2),
                    censored_next = c(0, 0, 1, 1))
  cb <- censoring_balance(toy, "f")
  expect_equal(cb$mean_censored, 2)
  expect_equal(cb$mean_retained, 0)
  expect_equal(cb$smd, 2)
  w <- small_weighted()
  cbw <- censoring_balance(w, c("X", "L"), weighted = TRUE)
  cbu <- censoring_balance(w, c("X", "L"), weighted = FALSE)
  expect_lt(mean(abs(cbw$smd)), mean(abs(cbu$smd)))
})

test_that("aggregate_balance computes IAI and MaxIB", {
  x <- data.frame(feature = "L", time = 0:3, rrb = c(0, 0.1, -0.3, 0.2))
  a <- aggregate_balance(x)
  expect_equal(a$iai, 0.15)
  expect_equal(a$maxib, 0.3)
  expect_equal(a$time_at_max, 2)
  x$q <- c(0, 0, 1, 1)
  expect_equal(aggregate_balance(x, time_weight = "q")$iai, 0.25)
  expect_equal(aggregate_balance(x, times = 0:1)$maxib, 0.1)
})

test_that("weight_support gives ESS = N and entropy 1 for equal weights", {
  d <- data.frame(strategy = "g", time = 0, censored_next = c(0, 0, 0, 1),
                  w_post = 2, pi_hat = c(1, 0.5, 0.04, 0.2))
  ws <- weight_support(d)
  expect_equal(ws$ess, 3)
  expect_equal(ws$ess_ratio, 0.75)
  expect_equal(ws$entropy, 1)
  expect_equal(ws$pi_min, 0.04)
  expect_equal(ws$pr_pi_small, 0.25)
  w <- small_weighted()
  ws2 <- weight_support(w)
  expect_true(all(ws2$ess_ratio > 0 & ws2$ess_ratio <= 1 + 1e-12))
})

test_that("cumulative IPW identity holds exactly on a toy example", {
  toy <- toy_cum()
  r <- rrb(toy, "f", method = c("unstabilized", "stabilized"),
           scope = "cumulative")
  r1 <- r[r$time == 1, ]
  expect_equal(r1$mean_pre, c(2, 2))
  expect_equal(r1$mean_post, c(2, 2))
  expect_equal(r1$rrb, c(0, 0))
  expect_equal(r1$n_censored, c(3, 3))
  expect_equal(r1$n_risk, c(10, 10))
  # without carrying the event clone forward the identity fails
  rn <- rrb(toy, "f", method = "unstabilized", scope = "cumulative",
            carry_events = FALSE)
  expect_gt(abs(rn$rrb[rn$time == 1]), 0.1)
  # unweighted cumulative selection is not balanced
  ru <- rrb(toy, "f", method = "unweighted", scope = "cumulative")
  expect_gt(abs(ru$rrb[ru$time == 1]), 0.1)
  expect_equal(ks_rr(toy, "f", "unstabilized", scope = "cumulative")$ks,
               c(0, 0))
  m <- mahalanobis_rr(toy, "f", "stabilized", scope = "cumulative")
  expect_equal(m$mahalanobis, c(0, 0))
})

test_that("cumulative RRB of baseline X is near 0 under correct spec", {
  d <- simulate_ccw(ccw_sim_params(n = 5000), seed = 2025)
  d <- make_history_features(d)
  feats <- c("X", "L_base")
  w_s <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
  w_u <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X, stabilized = FALSE)
  rs <- rrb(w_s, feats, method = c("unweighted", "weighted", "stabilized"),
            scope = "cumulative")
  ru <- rrb(w_u, feats, method = "weighted", scope = "cumulative")
  rs <- rs[rs$n_censored > 0, ]
  ru <- ru[ru$n_censored > 0, ]
  expect_lt(max(abs(rs$rrb[rs$method != "unweighted"])), 0.1)
  expect_lt(max(abs(ru$rrb)), 0.1)
  expect_gt(max(abs(rs$rrb[rs$method == "unweighted"])), 0.2)
  # stabilized weights via rrb(method = "stabilized") and via w_post agree
  a <- rs[rs$method == "stabilized", "rrb"]
  b <- rs[rs$method == "weighted", "rrb"]
  expect_equal(a, b)
  # post-baseline L_t is NOT 0-targeted cumulatively (treatment affects L)
  rl <- rrb(w_s, "L", scope = "cumulative")
  expect_gt(max(abs(rl$rrb), na.rm = TRUE), 0.1)
})

test_that("cumulative scope warns for covariate-dependent numerators", {
  w <- small_weighted(numerator = ~ factor(time) + X)
  expect_warning(rrb(w, "X", method = "stabilized", scope = "cumulative"),
                 "varies within")
})

test_that("scope is validated and plots accept cumulative output", {
  w <- small_weighted()
  expect_error(rrb(w, "X", scope = "nope"), "Unknown scope")
  r <- rrb(w, c("X", "L_base"), method = c("unweighted", "weighted"),
           scope = c("interval", "cumulative"))
  expect_setequal(unique(r$scope), c("interval", "cumulative"))
  expect_s3_class(plot_love_dynamic(r), "ggplot")
  expect_s3_class(plot_love_dynamic(r[r$scope == "cumulative", ]), "ggplot")
  expect_s3_class(plot_balance_heatmap(
    r[r$scope == "cumulative" & r$method == "weighted", ]), "ggplot")
  expect_error(plot_balance_heatmap(r[r$method == "weighted", ]),
               "single scope")
  m <- mahalanobis_rr(w, c("X", "L_base"), scope = c("interval",
                                                     "cumulative"))
  expect_s3_class(plot_mahalanobis(m), "ggplot")
  a <- aggregate_balance(r)
  expect_true("scope" %in% names(a))
})

test_that("weight_support returns NA weight summaries when nobody is retained", {
  d <- data.frame(strategy = "g", time = 0, censored_next = c(1, 1),
                  w_post = 2, pi_hat = 0.3)
  ws <- expect_silent(weight_support(d))
  expect_true(is.na(ws$w_max) && is.na(ws$w99) && is.na(ws$ess))
  expect_equal(ws$n_retained, 0)
})

test_that("censoring_balance estimator weights match model without truncation", {
  d <- simulate_ccw(ccw_sim_params(n = 400), seed = 3)
  f <- ~ L + I(L^2) + X + L:X
  for (st in c(TRUE, FALSE)) {
    w <- fit_censoring_model(d, f, stabilized = st)
    cm <- censoring_balance(w, c("X", "L"), weighted = TRUE)
    ce <- censoring_balance(w, c("X", "L"), weighted = TRUE,
                            weights = "estimator")
    expect_true(all(cm$weights == "model"))
    expect_true(all(ce$weights == "estimator"))
    expect_equal(ce$smd, cm$smd, tolerance = 1e-10)
  }
  expect_true(all(is.na(censoring_balance(w, "X")$weights)))
  wt <- fit_censoring_model(d, f, truncate = 0.9)
  cm <- censoring_balance(wt, c("X", "L"), weighted = TRUE)
  ce <- censoring_balance(wt, c("X", "L"), weighted = TRUE,
                          weights = "estimator")
  expect_gt(max(abs(ce$smd - cm$smd), na.rm = TRUE), 1e-6)
  ci <- censoring_balance(wt, c("X", "L"), weighted = TRUE,
                          weights = "estimator", cap = Inf)
  expect_false(isTRUE(all.equal(ci$smd, ce$smd)))
})
