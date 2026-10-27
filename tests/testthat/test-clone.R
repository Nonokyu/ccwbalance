test_that("ccw_clone applies the censoring rules", {
  obs <- data.frame(id = 1:4, tau = c(1, 5, Inf, 0), ev = c(Inf, 4, 2, Inf),
                    lt = c(Inf, Inf, Inf, 3), age = c(70, 80, 65, 60))
  d <- ccw_clone(obs, list(early = c(0, 2), late = c(3, 5), never = c(NA, 6)),
                 n_time = 7, init_time = "tau", event_time = "ev",
                 ltfu_time = "lt", covariates = "age")
  get <- function(g, i) d[d$strategy == g & d$id == i, ]
  # early: id 1 initiates at 1 -> retained to the end
  expect_equal(max(get("early", 1)$time), 6)
  expect_equal(sum(get("early", 1)$censored_next), 0)
  # early: id 2 initiates at 5 -> censored at grace end (t = 2)
  e2 <- get("early", 2)
  expect_equal(max(e2$time), 2)
  expect_equal(e2$censored_next[e2$time == 2], 1)
  expect_equal(e2$phase[e2$time == 2], "grace_end")
  # early: id 3 has the event in the grace-end interval 2: by default the
  # outcome comes first (not censored, not in the grace-end model) ...
  e3 <- get("early", 3)
  expect_equal(e3$censored_next[e3$time == 2], 0)
  expect_equal(e3$event_next[e3$time == 2], 1)
  expect_equal(e3$phase[e3$time == 2], "none")
  # ... unless censoring is ordered first
  d2 <- ccw_clone(obs, list(early = c(0, 2)), n_time = 7, init_time = "tau",
                  event_time = "ev", ltfu_time = "lt",
                  outcome_first_at_grace_end = FALSE)
  e3b <- d2[d2$id == 3, ]
  expect_equal(e3b$censored_next[e3b$time == 2], 1)
  expect_true(is.na(e3b$event_next[e3b$time == 2]))
  # late: id 1 initiates at 1 < 3 -> censored at 1 (no_init phase)
  l1 <- get("late", 1)
  expect_equal(max(l1$time), 1)
  expect_equal(l1$phase, c("no_init", "no_init"))
  # late: id 2 initiates at 5 (in window), event at 4
  l2 <- get("late", 2)
  expect_equal(l2$event_next[l2$time == 4], 1)
  expect_equal(sum(l2$censored_next), 0)
  # never: id 4 initiates at 0 -> censored at 0
  n4 <- get("never", 4)
  expect_equal(nrow(n4), 1)
  expect_equal(n4$censored_next, 1)
  # early: id 4 LTFU at 3
  e4 <- get("early", 4)
  expect_equal(e4$ltfu_next[e4$time == 3], 1)
  expect_true(is.na(e4$event_next[e4$time == 3]))
  expect_equal(unique(get("never", 1)$age), 70)
  expect_error(ccw_clone(obs, list(c(0, 2)), 7, init_time = "tau",
                         event_time = "ev", ltfu_time = "lt"), "named")
})

test_that("ccw_sustained validates and prints", {
  expect_s3_class(ccw_sustained(1), "ccw_sustained")
  expect_error(ccw_sustained(2))
  expect_error(ccw_sustained(1, grace = -1))
  expect_output(print(ccw_sustained(0, 1)), "grace = 1")
})

test_that("ccw_clone sustained mode censors at every deviation", {
  pp <- data.frame(
    id = rep(1:3, each = 4), t = rep(0:3, 3),
    a = c(1, 1, 0, 1,   0, 0, 0, 0,   1, 0, 0, 1),
    y = c(0, 0, 0, 0,   0, 1, NA, NA, 0, 0, 0, 1),
    L = 1:12
  )
  pp <- pp[!(pp$id == 2 & pp$t > 1), ]   # rows end after the outcome
  d <- ccw_clone(pp, list(always = ccw_sustained(1),
                          never = ccw_sustained(0),
                          always_g1 = ccw_sustained(1, grace = 1)),
                 time = "t", treatment = "a", event = "y", covariates = "L")
  get <- function(g, i) d[d$strategy == g & d$id == i, ]
  # always: id 1 deviates at t = 2
  a1 <- get("always", 1)
  expect_equal(a1$time, 0:2)
  expect_equal(a1$censored_next, c(0, 0, 1))
  expect_true(all(a1$censor_possible))
  expect_true(is.na(a1$event_next[3]))
  # never: id 2 adherent, outcome at t = 1, rows stop
  n2 <- get("never", 2)
  expect_equal(n2$event_next, c(0, 1))
  expect_equal(sum(n2$censored_next), 0)
  # never: id 1 deviates at t = 0
  expect_equal(get("never", 1)$censored_next, 1)
  # always with grace 1: id 1 (one off interval) is not censored
  g1 <- get("always_g1", 1)
  expect_equal(nrow(g1), 4)
  expect_equal(sum(g1$censored_next), 0)
  # deviation only possible after one off interval (run_off == grace)
  expect_equal(g1$run_off, c(0, 0, 0, 1))
  expect_equal(g1$phase, c("none", "none", "none", "sustain"))
  expect_equal(g1$censor_possible, c(FALSE, FALSE, FALSE, TRUE))
  # id 3: two consecutive off intervals -> deviation at t = 2
  g3 <- get("always_g1", 3)
  expect_equal(max(g3$time), 2)
  expect_equal(g3$censored_next[g3$time == 2], 1)
  expect_equal(g3$L, 9:11)
  # sustained strategies need person-period data
  expect_error(ccw_clone(pp, list(a = ccw_sustained(1)), n_time = 4),
               "person-period")
  expect_error(ccw_clone(pp, list(a = c(0, 1)), time = "t", treatment = "a",
                         event = "y"), "ccw_sustained")
})

test_that("sustained clones feed fit_censoring_model and rrb", {
  set.seed(4)
  n <- 3000
  pp <- data.frame(id = rep(seq_len(n), each = 4), t = rep(0:3, n))
  pp$L <- stats::rnorm(nrow(pp))
  pp$a <- stats::rbinom(nrow(pp), 1, stats::plogis(1 + pp$L))
  pp$y <- stats::rbinom(nrow(pp), 1, 0.05)
  d <- ccw_clone(pp, list(always = ccw_sustained(1),
                          never = ccw_sustained(0)),
                 time = "t", treatment = "a", event = "y", covariates = "L")
  w <- fit_censoring_model(d, ~ L + time, by = "strategy")
  expect_true(all(w$pi_hat > 0 & w$pi_hat <= 1))
  r <- rrb(w, "L", method = c("unweighted", "weighted"))
  expect_true(all(r$n_censored > 0))
  r <- r[r$n_risk >= 200, ]
  expect_lt(max(abs(r$rrb[r$method == "weighted"]), na.rm = TRUE),
            max(abs(r$rrb[r$method == "unweighted"]), na.rm = TRUE))
})
