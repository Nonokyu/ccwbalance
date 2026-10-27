boot_fixture <- function() {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
  fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
}

test_that("rrb_boot returns one row per observed cell for both scopes", {
  w <- boot_fixture()
  for (sc in c("cumulative", "interval")) {
    bt <- rrb_boot(w, c("X", "L"), B = 5, scope = sc, seed = 1)
    obs <- rrb(w, c("X", "L"), scope = sc)
    expect_equal(nrow(bt), nrow(obs))
    expect_equal(bt$rrb_obs, obs$rrb)
    expect_true(all(c("strategy", "time", "feature", "n_risk", "rrb_obs",
                      "boot_mean", "boot_se", "lcl", "ucl", "n_boot_ok",
                      "B", "scope") %in% names(bt)))
    expect_true(all(bt$B == 5))
    expect_true(all(bt$scope == sc))
    expect_true(all(bt$n_boot_ok <= 5))
    # cells with a structurally zero RRB (baseline time in the cumulative
    # scope, no censoring in the interval scope) have zero spread
    nz <- abs(bt$rrb_obs) > 1e-12
    se <- bt$boot_se[nz & !is.na(bt$boot_se)]
    expect_gt(length(se), 0.8 * sum(nz))
    expect_true(all(se > 0))
    ok <- !is.na(bt$lcl)
    expect_true(all(bt$lcl[ok] <= bt$ucl[ok]))
  }
})

test_that("rrb_boot is reproducible with a seed", {
  w <- boot_fixture()
  a <- rrb_boot(w, "L", B = 3, seed = 9)
  b <- rrb_boot(w, "L", B = 3, seed = 9)
  expect_identical(a, b)
})

test_that("the identity resample reproduces the observed RRB", {
  w <- boot_fixture()
  base <- as.data.frame(w)[, setdiff(names(w), .boot_weight_cols())]
  idx <- split(seq_len(nrow(base)), base$id)
  rb <- .boot_rrb(base, idx, attr(w, "weighting"), c("X", "L"), "weighted",
                  "cumulative", "id", "strategy", "time")
  obs <- rrb(w, c("X", "L"), scope = "cumulative")
  expect_equal(rb$rrb, obs$rrb)
  expect_equal(rb$n_risk, obs$n_risk)
})

test_that("rrb_boot requires fit_censoring_model() output", {
  d <- simulate_ccw(ccw_sim_params(n = 100), seed = 1)
  expect_error(rrb_boot(d, "L", B = 2), "weighting")
})
