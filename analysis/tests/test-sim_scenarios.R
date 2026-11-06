# Tests of the simulation scenario grids defined in analysis/sim_scenarios.R.
# Run from the project root:
#   Rscript -e 'pkgload::load_all("."); testthat::test_file(
#     "analysis/tests/test-sim_scenarios.R")'

local({
  root <- if (file.exists(file.path("analysis", "sim_scenarios.R"))) {
    "."
  } else {
    file.path("..", "..")
  }
  source(file.path(root, "analysis", "sim_scenarios.R"), local = globalenv())
})

test_that("ccw_scenarios() accepts only the four paper runs", {
  runs <- c("main", "manycov", "largerep", "variants")
  expect_identical(ccw_scenarios(), ccw_scenarios("main"))
  for (bad in c("rerun", "revision", "other")) {
    expect_error(ccw_scenarios(bad))
  }
  for (r in runs) {
    sc <- ccw_scenarios(r)
    tab <- scenario_table(sc)
    expect_equal(nrow(tab), length(sc))
    expect_true(all(tab$run == r))
    expect_false(anyNA(tab$misspec))
    for (s in sc) expect_type(do.call(ccw_sim_params, s$params), "list")
  }
  p <- do.call(ccw_sim_params, ccw_scenarios("main")$omit_lag_cs_os$params)
  d <- simulate_ccw(utils::modifyList(p, list(n = 100)), seed = 1)
  expect_s3_class(d, "data.frame")
})

test_that("run main has the three parts with the specified structure", {
  d <- ccw_scenarios("main")
  expect_length(d, 20)
  expect_false(anyDuplicated(names(d)) > 0)
  tab <- scenario_table(d)
  expect_equal(nrow(tab), 20)
  expect_true(all(c("run", "id", "part", "label", "misspec", "cens_strength",
                    "outcome_strength", "n_boot", "n") %in% names(tab)))
  expect_false("study" %in% names(tab))
  # ids of the scenarios the manuscript reports from this run
  expect_equal(as.vector(table(tab$part)[c("1", "2", "3")]), c(6, 2, 3))
  expect_equal(sum(is.na(tab$id)), 9)
  expect_true(all(tab$run == "main"))
  expect_setequal(unique(tab$misspec),
                  c("none", "omit_L", "omit_Lsq", "omit_lag"))
  expect_setequal(unique(tab$cens_strength), c("none", "weak", "strong"))
  expect_setequal(unique(tab$outcome_strength), c("none", "strong"))
  s1 <- tab[tab$misspec != "none", ]
  expect_equal(nrow(s1), 12)
  expect_true(all(table(s1$misspec, s1$cens_strength, s1$outcome_strength)
                  == 1))
  expect_true(all(s1$n_boot == 100))
  expect_equal(tab$n[tab$scenario == "correct_n300"], 300)
  expect_equal(tab$n[tab$scenario == "correct_n8000"], 8000)
  expect_equal(tab$truncate[tab$scenario == "truncation_95"], 0.95)
  # formulas
  fml <- function(nm) paste(deparse(d[[nm]]$formula), collapse = "")
  expect_equal(fml("omit_L_cs_os"), "~X")
  expect_equal(fml("omit_Lsq_cw_on"), "~L + X")
  expect_equal(fml("omit_lag_cs_os"), "~L + I(L^2) + X + L:X")
  expect_equal(fml("correct"), fml("extra_ltfu"))
  # coefficients
  pd <- function(nm) do.call(ccw_sim_params, d[[nm]]$params)
  expect_equal(pd("omit_Lsq_cw_on")$beta[["b5"]], 0)
  expect_equal(pd("omit_Lsq_cs_os")$beta[["b5"]], 0)
  expect_equal(pd("omit_L_cs_on")$beta[["b1"]], 1.0)
  expect_equal(pd("omit_L_cw_on")$gamma[["g1"]], 0)
  expect_equal(pd("omit_Lsq_cs_os")$gamma[["g4"]], 1.0)
  expect_equal(pd("omit_lag_cw_os")$beta[["b6"]], 0.6)
  expect_equal(pd("omit_L_cs_os")$gamma[["g1"]], 1.5)
  expect_equal(pd("omit_L_cw_os")$gamma[["g0"]], -4)
  expect_equal(pd("omit_Lsq_cw_os")$gamma[["g0"]], -4.5)
  expect_equal(pd("omit_lag_cs_os")$gamma[["g6"]], 1.2)
  expect_equal(pd("omit_lag_cw_os")$gamma[["g0"]], -4)
  expect_equal(pd("omit_lag_cw_on")$gamma[["g0"]], -4)
  expect_equal(pd("near_positivity")$beta[["b0"]], -2)
  expect_gt(pd("extra_ltfu")$delta[["d1"]], 0)
})

test_that("run manycov scenarios are well formed", {
  e <- ccw_scenarios("manycov")
  expect_length(e, 6)
  expect_equal(names(e), c("many_cov", "many_cov_n300", "rare_censoring",
                           "extreme_weights", "many_cov_omit_X2",
                           "many_cov_omit_X2_rho0"))
  expect_true(all(vapply(e, `[[`, "", "run") == "manycov"))
  expect_equal(e$many_cov_n300$id, "P2-05")
  expect_equal(e$extreme_weights$id, "P2-06")
  expect_null(e$many_cov$id)
  for (s in e) {
    expect_equal(length(s$formula), 2)  # one-sided
    p <- do.call(ccw_sim_params, s$params)
    expect_equal(s$extra_feats, if (p$n_extra > 0) {
      paste0("X", seq_len(p$n_extra) + 1)
    } else character(0))
    expect_true(all(all.vars(s$formula) %in% c("L", "X", s$extra_feats)))
  }
  expect_false("X2" %in% all.vars(e$many_cov_omit_X2$formula))
  expect_true("X9" %in% all.vars(e$many_cov_omit_X2$formula))
  expect_equal(do.call(ccw_sim_params, e$many_cov_omit_X2_rho0$params)$
                 rho_extra, 0)
  expect_equal(do.call(ccw_sim_params, e$many_cov_omit_X2$params)$
                 beta_extra[1], 1)
  tab <- scenario_table(e)
  expect_true(all(c("n_extra", "extra_feats") %in% names(tab)))
  expect_equal(tab$n_extra, c(8, 8, 0, 0, 8, 8))
  expect_equal(tab$n[tab$scenario == "many_cov_n300"], 300)
  expect_equal(tab$extra_feats[tab$scenario == "rare_censoring"], "")
  expect_equal(tab$extra_feats[1], paste(paste0("X", 2:9), collapse = ","))
  expect_true(all(scenario_table(ccw_scenarios("main"))$n_extra == 0))
})

test_that("run largerep scenarios are well formed", {
  f <- ccw_scenarios("largerep")
  expect_equal(names(f), c("correct_1000", "correct_n300_1000",
                           "correct_n8000_500", "early_events",
                           "omit_L_ss_trunc", "omit_L_ss_nearpos"))
  expect_true(all(vapply(f, `[[`, "", "run") == "largerep"))
  expect_equal(f$correct_1000$id, "P1-01")
  expect_null(f$early_events$id)
  tab <- scenario_table(f)
  expect_equal(nrow(tab), 6)
  expect_equal(tab$n, c(2000, 300, 8000, 2000, 2000, 2000))
  expect_equal(tab$n_boot, c(500, 500, 200, 100, 100, 100))
  expect_equal(tab$n_rep, c(1000, 1000, 500, 200, 200, 200))
  expect_equal(tab$truncate[tab$scenario == "omit_L_ss_trunc"], 0.95)
  pf <- function(nm) do.call(ccw_sim_params, f[[nm]]$params)
  expect_equal(pf("omit_L_ss_nearpos")$beta[["b1"]], 1)
  expect_equal(pf("omit_L_ss_nearpos")$beta[["b0"]], -2)
  expect_equal(pf("omit_L_ss_trunc")$gamma[["g1"]], 1.5)
  expect_equal(deparse(f$omit_L_ss_trunc$formula), "~X")
  expect_lt(pf("early_events")$gamma[["g2"]], -1)
  # runs main and manycov carry no n_rep (their replicate counts are
  # unchanged)
  expect_true(all(is.na(scenario_table(ccw_scenarios("main"))$n_rep)))
  expect_true(all(is.na(scenario_table(ccw_scenarios("manycov"))$n_rep)))
})

test_that("run variants scenarios are well formed", {
  h <- ccw_scenarios("variants")
  expect_equal(names(h), c("omit_L1_cw_on", "omit_L1_cw_os", "omit_L1_cs_on",
                           "omit_L1_cs_os", "omit_Lsq_cw_on_g45",
                           "omit_Lsq_cs_on_g45", "early_events_b500"))
  expect_true(all(vapply(h, `[[`, "", "run") == "variants"))
  expect_equal(unname(vapply(h, `[[`, "", "id")),
               c(sprintf("P1-%02d", c(2:6, 8)), "P3-04"))
  expect_true(all(vapply(h, `[[`, 0, "n_rep") == 200))
  ph <- function(nm) do.call(ccw_sim_params, h[[nm]]$params)
  expect_equal(ph("omit_L1_cs_os")$beta[c("b1", "b2", "b5")],
               c(b1 = 1, b2 = 0, b5 = 0))
  expect_equal(ph("omit_L1_cs_os")$gamma[["g1"]], 1.5)
  expect_equal(deparse(h$omit_L1_cw_on$formula), "~X")
  expect_equal(ph("omit_Lsq_cw_on_g45")$gamma[c("g0", "g4")],
               c(g0 = -4.5, g4 = 0))
  expect_equal(h$early_events_b500$n_boot, 500)
  expect_equal(nrow(scenario_table(h)), 7)
})

test_that("paper_scenarios() follows Table tab:grid of the manuscript", {
  ps <- paper_scenarios()
  ids <- c(sprintf("P1-%02d", 1:13), sprintf("P2-%02d", 1:6),
           sprintf("P3-%02d", 1:4))
  # 23 scenarios in the row order of Table tab:grid; no E ids
  expect_equal(ps$id, ids)
  expect_false(anyDuplicated(ps$id) > 0)
  expect_false(any(grepl("^E-", ps$id)))
  expect_equal(ps$part, rep(1:3, c(13, 6, 4)))
  fac <- c("weak / none", "weak / strong", "strong / none",
           "strong / strong")
  expect_equal(ps$label, c(
    "Correct specification, n = 2000", paste0("Omitted L_k, ", fac),
    paste0("Omitted L_k^2, ", fac), paste0("Omitted L_{k-1}, ", fac),
    "Correct specification, n = 300", "Correct specification, n = 8000",
    "Near-violation of positivity", "Weight truncation, 95th percentile",
    "Many correlated covariates", "Extreme weights",
    "No treatment->covariate effect", "Strong treatment->covariate effect",
    "Unweighted loss to follow-up", "Many early events"))
  expect_true(all(ps$run %in% c("main", "manycov", "largerep", "variants")))
  # every id resolves to an existing run scenario carrying the same id
  for (k in seq_len(nrow(ps))) {
    sc <- ccw_scenarios(ps$run[k])
    expect_true(ps$scenario[k] %in% names(sc))
    expect_equal(sc[[ps$scenario[k]]]$id, ps$id[k])
    expect_equal(sc[[ps$scenario[k]]]$part, ps$part[k])
    expect_equal(sc[[ps$scenario[k]]]$label, ps$label[k])
    expect_identical(paper_scenario(ps$id[k]), sc[[ps$scenario[k]]])
  }
  expect_error(paper_scenario("P9-99"))
  # no other run scenario carries an id
  n_id <- sum(vapply(c("main", "manycov", "largerep", "variants"), function(r) {
    sum(!is.na(scenario_table(ccw_scenarios(r))$id))
  }, 0))
  expect_equal(n_id, nrow(ps))
  # specs of Table tab:grid: n, replicates, bootstrap resamples
  expect_equal(ps$n, c(rep(2000, 13), 300, 8000, 2000, 2000, 300, 2000,
                       rep(2000, 4)))
  expect_equal(ps$n_rep, c(1000, rep(200, 12), 1000, 500, rep(200, 8)))
  expect_equal(ps$n_boot, c(500, rep(100, 12), 500, 200, 500, 500, 100, 100,
                            100, 100, 100, 500))
  expect_equal(ps$truncate[ps$id == "P2-04"], 0.95)
  expect_equal(ps$censoring_formula[2:5], rep("~X", 4))
  expect_equal(ps$censoring_formula[6:9], rep("~L + X", 4))
  pp <- function(id) do.call(ccw_sim_params, paper_scenario(id)$params)
  expect_equal(unname(pp("P1-05")$beta[c("b1", "b2", "b5")]), c(1, 0, 0))
  expect_equal(pp("P1-06")$gamma[["g0"]], -4.5)
  expect_equal(pp("P1-07")$gamma[["g0"]], -4.5)
  expect_equal(pp("P3-04")$gamma[["g2"]], -3.5)
  # Part 1 factor levels (omitted L_k cells say "omit_L", not "omit_L1")
  p1 <- ps[ps$part == 1 & ps$misspec != "none", ]
  expect_equal(unique(p1$misspec), c("omit_L", "omit_Lsq", "omit_lag"))
  expect_equal(p1$cens_strength, rep(rep(c("weak", "strong"), each = 2), 3))
  expect_equal(p1$outcome_strength, rep(c("none", "strong"), 6))
  # superseded rows of earlier runs
  expect_equal(ps$supersedes[ps$id %in% c("P1-01", "P2-01", "P2-02", "P3-04")],
               c("main:correct", "main:correct_n300", "main:correct_n8000",
                 "largerep:early_events"))
  # the committed table matches
  root <- if (file.exists(file.path("analysis", "scenarios_paper.csv"))) {
    "."
  } else {
    file.path("..", "..")
  }
  csv <- utils::read.csv(file.path(root, "analysis", "scenarios_paper.csv"),
                         na.strings = character(0))
  expect_equal(csv$id, ps$id)
  expect_equal(csv$run, ps$run)
  expect_equal(csv$scenario, ps$scenario)
  expect_equal(csv$supersedes, ps$supersedes)
  expect_equal(csv$param_overrides, ps$param_overrides)
})
