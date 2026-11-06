# Studentised maximum and censorable-row censoring balance -------------------
stud_fixture <- function() {
  d <- simulate_ccw(ccw_sim_params(n = 300), seed = 1)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  fit_censoring_model(d, ~ L + I(L^2) + X + L:X)
}

test_that("studentised maxima cover every MaxIB set", {
  w <- stud_fixture()
  nr <- rrb_null(w, c("X", "L", "L_diff"), n_draws = 30, seed = 1)
  for (nm in c("maxib_stud", "maxib_cb_stud")) {
    ms <- nr[[nm]]
    ref <- if (nm == "maxib_stud") nr$maxib else nr$maxib_cb
    expect_identical(ms[, c("strategy", "feature")],
                     ref[, c("strategy", "feature")])
    expect_true(all(c("n_cells", "n_excluded", "stud_obs", "null_mean",
                      "null_q95", "p") %in% names(ms)))
    expect_true(all(ms$n_excluded >= 0))
    has <- ms$n_cells > 0
    expect_true(all(is.finite(ms$null_q95[has])))
    expect_true(all(ms$p[has] > 0 & ms$p[has] <= 1))
    expect_equal(ms$n_cells + ms$n_excluded, ref$n_cells)
  }
  z <- nr$cells$z
  all_row <- nr$maxib_stud$strategy == "(all)" &
    nr$maxib_stud$feature == "(all)"
  expect_equal(nr$maxib_stud$stud_obs[all_row], max(abs(z), na.rm = TRUE))
  for (g in c("g1", "g2", "g3")) {
    r <- nr$maxib_stud$strategy == g & nr$maxib_stud$feature == "(all)"
    zg <- z[nr$cells$strategy == g]
    expect_equal(nr$maxib_stud$stud_obs[r], max(abs(zg), na.rm = TRUE))
  }
  cbz <- nr$cells$cb_z
  all_cb <- nr$maxib_cb_stud$strategy == "(all)" &
    nr$maxib_cb_stud$feature == "(all)"
  expect_equal(nr$maxib_cb_stud$stud_obs[all_cb], max(abs(cbz), na.rm = TRUE))
  expect_output(print(nr), "Studentised maximum")
})

test_that("studentised maximum excludes cells without a positive null SD", {
  cells <- data.frame(strategy = "g", time = 1:3, feature = "f",
                      rrb_obs = c(0.2, -0.6, 0.1))
  sims <- matrix(c(0.1, -0.1, 0.3, 0.2, -0.2, 0.1, 0, 0, 0), 3,
                 byrow = TRUE)
  sc <- c(0.5, NA, 0)
  m <- .null_maxib(cells, sims, scale = sc)
  r <- m[m$strategy == "g" & m$feature == "f", ]
  expect_equal(r$n_cells, 1L)
  expect_equal(r$n_excluded, 2L)
  expect_equal(r$stud_obs, 0.4)
  expect_equal(r$null_q95,
               stats::quantile(c(0.2, 0.2, 0.6), 0.95, names = FALSE))
  expect_equal(r$p, (1 + 1) / 4)
})

test_that("cb_cens_* elements give the censorable-row censoring balance", {
  w <- stud_fixture()
  f <- c("X", "L")
  nr <- rrb_null(w, f, n_draws = 20, seed = 3)
  cl <- nr$cells
  expect_true(all(c("cb_cens_obs", "cb_cens_null_mean", "cb_cens_null_sd",
                    "cb_cens_z", "cb_cens_p_two") %in% names(cl)))
  expect_identical(nr$maxib_cb_cens[, c("strategy", "feature")],
                   nr$maxib_cb[, c("strategy", "feature")])
  expect_true("stud_obs" %in% names(nr$maxib_cb_cens_stud))
  expect_true("cb_cens" %in% names(nr$draws))
  # cells whose risk set holds non-censorable rows
  nonposs <- as.data.frame(w)[!w$censor_possible, c("strategy", "time")]
  key_np <- unique(paste(nonposs$strategy, nonposs$time))
  key <- paste(cl$strategy, cl$time)
  pure <- !key %in% key_np
  expect_gt(sum(pure), 0)
  for (v in c("obs", "null_mean", "null_sd", "z", "p_two")) {
    expect_identical(cl[[paste0("cb_cens_", v)]][pure],
                     cl[[paste0("cb_", v)]][pure])
  }
  dk <- paste(nr$draws$strategy, nr$draws$time)
  expect_identical(nr$draws$cb_cens[!dk %in% key_np],
                   nr$draws$cb[!dk %in% key_np])
  mixed <- !pure & !is.na(cl$cb_obs)
  expect_gt(sum(mixed), 0)
  expect_true(all(abs(cl$cb_obs[mixed] - cl$cb_cens_obs[mixed]) > 1e-8))
  # observed censorable CB equals censoring_balance(rows = "censorable")
  cb <- censoring_balance(w, f, weighted = TRUE, rows = "censorable")
  k <- match(paste(cl$strategy, cl$time, cl$feature),
             paste(cb$strategy, cb$time, cb$feature))
  ok <- !is.na(k)
  expect_gt(sum(ok), 0)
  expect_equal(cl$cb_cens_obs[ok], cb$smd[k[ok]], tolerance = 1e-8)
  expect_true(all(is.na(cl$cb_cens_obs[!ok])))
  # the censorable maxima are those of the cb_cens columns
  all_row <- nr$maxib_cb_cens$strategy == "(all)" &
    nr$maxib_cb_cens$feature == "(all)"
  expect_equal(nr$maxib_cb_cens$maxib_obs[all_row],
               max(abs(cl$cb_cens_obs), na.rm = TRUE))
  all_st <- nr$maxib_cb_cens_stud$strategy == "(all)" &
    nr$maxib_cb_cens_stud$feature == "(all)"
  expect_equal(nr$maxib_cb_cens_stud$stud_obs[all_st],
               max(abs(cl$cb_cens_z), na.rm = TRUE))
  expect_output(print(nr), "censorable rows")
})

test_that("censoring_balance(rows = 'censorable') matches a hand computation", {
  toy <- data.frame(strategy = "g", time = 0,
                    f = c(0, 1, 2, 4, 5, 9),
                    censored_next = c(1, 0, 0, 1, 0, 0),
                    censor_possible = c(TRUE, TRUE, TRUE, TRUE, FALSE,
                                        FALSE))
  cb <- censoring_balance(toy, "f", rows = "censorable")
  expect_equal(cb$n_censored, 2)
  expect_equal(cb$n_retained, 2)
  expect_equal(cb$mean_censored, 2)
  expect_equal(cb$mean_retained, 1.5)
  # SD over the whole risk set (population form, as .wvar)
  s <- sqrt(mean((toy$f - mean(toy$f))^2))
  expect_equal(cb$sd, s)
  expect_equal(cb$smd, 0.5 / s)
  ca <- censoring_balance(toy, "f")
  expect_equal(ca$n_retained, 4)
  expect_equal(ca$mean_retained, 17 / 4)
  expect_equal(ca$sd, s)
  expect_error(censoring_balance(toy[, -5], "f", rows = "censorable"),
               "missing required column")
  # no censorable retained row: the cell is dropped
  toy2 <- toy
  toy2$censor_possible <- toy2$censored_next == 1
  expect_equal(nrow(censoring_balance(toy2, "f", rows = "censorable")), 0)
})
