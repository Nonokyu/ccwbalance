test_that("make_history_features computes the transforms", {
  d <- data.frame(id = c(1, 1, 1), strategy = "g", time = 0:2,
                  L = c(1, 3, 2), X = 2)
  out <- make_history_features(d, thresholds = c(0, 2.5),
                               interactions = list(c("L", "X")))
  expect_equal(out$L_sq, c(1, 9, 4))
  expect_equal(out$L_diff, c(0, 2, -1))
  expect_equal(out$L_lagprod, c(1, 3, 6))
  expect_equal(out$L_cummean, c(1, 2, 2))
  expect_equal(out$L_gt0, c(1, 1, 1))
  expect_equal(out$L_gt2.5, c(0, 1, 0))
  expect_equal(out$L_x_X, c(2, 6, 4))
  expect_equal(out$L_base, c(1, 1, 1))
})

test_that("make_history_features uses an existing lag column", {
  d <- data.frame(id = 1, strategy = "g", time = 0:1, L = c(1, 3),
                  L_lag = c(0.5, 1))
  out <- make_history_features(d)
  expect_equal(out$L_diff, c(0.5, 2))
})

test_that("make_history_features handles several variables", {
  d <- data.frame(id = rep(1:2, each = 2), strategy = "g",
                  time = c(0, 1, 0, 1),
                  L1 = c(1, 2, 0, 4), L2 = c(3, 1, 2, 2))
  out <- make_history_features(d, vars = c("L1", "L2"))
  expect_equal(out$L1_sq, c(1, 4, 0, 16))
  expect_equal(out$L2_diff, c(0, -2, 0, 0))
  expect_equal(out$L1_cummean, c(1, 1.5, 0, 2))
  expect_equal(out$L2_base, c(3, 3, 2, 2))
})
