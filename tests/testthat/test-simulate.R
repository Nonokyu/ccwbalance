test_that("ccw_sim_params merges partial coefficient overrides", {
  p <- ccw_sim_params(beta = c(b2 = 0.9), n = 10)
  expect_equal(p$beta[["b2"]], 0.9)
  expect_equal(p$beta[["b1"]], ccw_sim_params()$beta[["b1"]])
  expect_equal(p$n, 10)
  expect_error(ccw_sim_params(beta = c(zz = 1)))
})

test_that("simulate_ccw returns the long-data contract", {
  d <- simulate_ccw(ccw_sim_params(n = 200), seed = 1)
  expect_true(all(c("id", "strategy", "time", "at_risk", "censored_next",
                    "censor_possible", "phase", "X", "L", "L_lag",
                    "event_next") %in% names(d)))
  expect_equal(levels(d$strategy), c("g1", "g2", "g3"))
  # baseline clones are identical across strategies
  b <- d[d$time == 0, ]
  expect_equal(nrow(b), 3 * 200)
  xs <- tapply(b$X, b$strategy, mean)
  expect_equal(unname(xs[1]), unname(xs[2]))
  expect_equal(unname(xs[1]), unname(xs[3]))
  # artificial censoring only where possible; each clone censored at most once
  expect_true(all(d$censored_next[!d$censor_possible] == 0))
  expect_true(all(tapply(d$censored_next, interaction(d$strategy, d$id),
                         sum) <= 1))
  # no rows after censoring
  last <- tapply(d$time, interaction(d$strategy, d$id), max)
  cen_t <- d$time[d$censored_next == 1]
  cen_k <- interaction(d$strategy, d$id)[d$censored_next == 1]
  expect_equal(as.vector(last[as.character(cen_k)]), cen_t)
  # outcome unobserved when censored
  expect_true(all(is.na(d$event_next[d$censored_next == 1])))
})

test_that("true_risk_ccw returns one risk per strategy", {
  tr <- true_risk_ccw(ccw_sim_params(), n_mc = 2000, seed = 2)
  expect_equal(nrow(tr), 3)
  expect_true(all(tr$risk > 0 & tr$risk < 1))
})

test_that("default DGM regression values (g2 = nested window)", {
  # fixed-seed values of the default DGM; any change flags a DGM change
  d <- simulate_ccw(ccw_sim_params(n = 200), seed = 1)
  expect_equal(nrow(d), 4404)
  expect_equal(sum(d$X * d$L * d$time), 7334.85215549933, tolerance = 1e-8)
  expect_false(any(grepl("^X[0-9]+$", names(d))))
  tr <- true_risk_ccw(ccw_sim_params(), n_mc = 5000, seed = 7)
  expect_equal(tr$risk, c(0.1158, 0.1400, 0.2248), tolerance = 1e-12)
})

test_that("g2 is the nested longer window containing g1", {
  p <- ccw_sim_params(n = 300)
  w1 <- p$window[1]
  w2 <- p$window[2]
  d <- simulate_ccw(p, seed = 1)
  g1 <- d[d$strategy == "g1", ]
  g2 <- d[d$strategy == "g2", ]
  # (a) g2 has no "no_init" phase; grace end only at w2 - 1
  expect_true(all(g2$phase %in% c("grace_end", "none")))
  expect_true(all(g2$time[g2$phase == "grace_end"] == w2 - 1))
  expect_true(all(g2$censored_next[g2$time < w2 - 1] == 0))
  # (b) every g1 row up to w1 - 1 also exists in g2
  k1 <- paste(g1$id, g1$time)[g1$time <= w1 - 1]
  expect_true(all(k1 %in% paste(g2$id, g2$time)))
  # (c) clones censored from g1 for non-initiation stay at risk in g2 at w1
  #     unless they had the outcome or LTFU by then
  cen1 <- g1$id[g1$censored_next == 1 & g1$phase == "grace_end"]
  expect_gt(length(cen1), 0)
  gone <- g2$id[g2$time <= w1 - 1 &
                  (g2$ltfu_next %in% 1 | g2$event_next %in% 1)]
  still <- setdiff(cen1, gone)
  expect_gt(length(still), 0)
  expect_true(all(still %in% g2$id[g2$time == w1]))
  # (d) treatment is protective (gamma2 < 0): longer window, higher risk
  tr <- true_risk_ccw(ccw_sim_params(), n_mc = 5000, seed = 7)
  expect_lt(tr$risk[1], tr$risk[2])
  expect_lt(tr$risk[2], tr$risk[3])
})

test_that("extra correlated baseline covariates are added", {
  p <- ccw_sim_params(n = 5000, n_extra = 3, rho_extra = 0.5)
  expect_equal(p$beta_extra, c(0.2, -0.2, 0.2))
  expect_equal(p$gamma_extra, rep(0.2, 3))
  d <- simulate_ccw(p, seed = 11)
  nm <- names(d)
  expect_equal(nm[match("X", nm) + 1:3], c("X2", "X3", "X4"))
  for (v in c("X2", "X3", "X4")) {
    expect_true(all(tapply(d[[v]], d$id, function(z) length(unique(z))) == 1))
  }
  b <- d[d$time == 0 & d$strategy == "g1", ]
  expect_lt(abs(stats::cor(b$X, b$X2) - 0.5), 0.1)
  expect_error(ccw_sim_params(n_extra = 3, beta_extra = c(1, 2)),
               "beta_extra")
  expect_error(ccw_sim_params(n_extra = 2, gamma_extra = 1), "gamma_extra")
  h <- make_history_features(d, interactions = list(c("L", "X")))
  expect_true(all(c("X2", "L_sq", "L_x_X") %in% names(h)))
  tr <- true_risk_ccw(p, n_mc = 2000, seed = 1)
  expect_true(all(tr$risk > 0 & tr$risk < 1))
})

test_that("outcome terms gamma4-6 default to 0 and enter the DGM", {
  p <- ccw_sim_params()
  expect_equal(unname(p$gamma[c("g4", "g5", "g6")]), c(0, 0, 0))
  r0 <- true_risk_ccw(ccw_sim_params(), n_mc = 5000, seed = 3)
  r1 <- true_risk_ccw(ccw_sim_params(gamma = c(g4 = 1)), n_mc = 5000,
                      seed = 3)
  expect_true(all(r1$risk > r0$risk))
})
