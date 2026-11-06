test_that("fit_censoring_model adds valid probabilities and weights", {
  w <- small_weighted()
  expect_true(all(w$pi_hat > 0 & w$pi_hat <= 1))
  expect_true(all(w$pi_hat[!w$censor_possible] == 1))
  # cumulative weight = previous weight * q / pi
  expect_equal(w$weight_stab, w$w_prev_stab * w$q_hat / w$pi_hat)
  expect_equal(w$weight_unstab, w$w_prev_unstab / w$pi_hat)
  expect_equal(w$w_ref, w$w_prev * w$q_hat)
  expect_true(all(w$w_prev[w$time == 0] == 1))
  expect_type(attr(w, "models"), "list")
  # default numerator is common within strategy x time
  qv <- tapply(w$q_hat, interaction(w$strategy, w$time, drop = TRUE),
               function(z) diff(range(z)))
  expect_true(all(qv < 1e-12))
  # formula numerator
  w2 <- small_weighted(numerator = ~ factor(time) + X)
  expect_true(all(w2$q_hat > 0 & w2$q_hat <= 1))
})

test_that("truncation caps weights", {
  w <- small_weighted(truncate = 0.9)
  cap <- attr(w, "weighting")$cap
  expect_true(is.finite(cap))
  expect_lte(max(w$w_post), cap + 1e-12)
})

test_that("estimate_risk returns monotone risks", {
  w <- small_weighted()
  r <- estimate_risk(w)
  expect_true(all(c("strategy", "time", "hazard", "risk") %in% names(r)))
  for (g in unique(r$strategy)) {
    expect_true(all(diff(r$risk[r$strategy == g]) >= -1e-12))
  }
})

test_that("fit_censoring_model accepts a named list of formulas", {
  d <- simulate_ccw(ccw_sim_params(n = 400), seed = 5)
  fl <- list(no_init = ~ L + X, grace_end = ~ X, "g2.grace_end" = ~ 1)
  w <- fit_censoring_model(d, fl)
  m <- attr(w, "models")
  expect_equal(names(coef(m[["g2.grace_end"]]$denominator)), "(Intercept)")
  expect_equal(names(coef(m[["g3.no_init"]]$denominator)),
               c("(Intercept)", "L", "X"))
  expect_equal(names(coef(m[["g1.grace_end"]]$denominator)),
               c("(Intercept)", "X"))
  expect_error(fit_censoring_model(d, list(no_init = ~ L)), "No formula")
})
