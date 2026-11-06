# Simulation scenario grids for the analysis scripts.
#
# Defines ccw_scenarios(), scenario_table(), paper_scenarios(),
# paper_scenario() and the .scenarios_* helpers (moved here from the
# package's former R/scenarios.R). These are plain functions, not part of
# ccwbalance. Load the package first, e.g.
#   pkgload::load_all(".")
#   source("analysis/sim_scenarios.R")
# The only package function called is ccw_sim_params() (exported, so it
# resolves with library(ccwbalance) as well as with load_all()).
#
# Runs and paper scenarios
# ------------------------
# A run is one execution of analysis/02_simulation_study.R
# (CCW_RUN=<run>) over the grid ccw_scenarios(<run>); its raw outputs are
# output/sim/sim_results_<run>.csv, sim_performance_<run>.csv,
# sim_settings_<run>.csv and analysis/scenarios_<run>.csv. The four runs
# supply the manuscript:
#   main      Parts 1-3 of the manuscript
#   manycov   many correlated covariates, extreme weights (and four
#             unreported scenarios)
#   largerep  large-replicate runs of the correct specification (1,000 /
#             1,000 / 500 replicates; and unreported scenarios)
#   variants  variants of main's scenarios (omitted L_t alone, L_t^2 cells
#             at gamma0 = -4.5, early events with 500 bootstrap resamples)
#
# Each scenario of the manuscript (Table "tab:grid" of the manuscript, in
# row order) has an id: P1-01 to P1-13
# (Part 1), P2-01 to P2-06 (Part 2), P3-01 to P3-04 (Part 3). Each id points
# to the run scenario whose results the manuscript reports; the column
# `supersedes` of paper_scenarios() records the rows of earlier runs that
# it replaces. Scenarios reported only in the
# Supporting Information or used only by the supplementary checks would
# carry ids E-01, E-02, ...; there are none: the supplementary checks
# (02d, 02e, 02f, 02g, 02j, 02l, 02m) and every table of the
# Supporting Information concern the mechanisms of P ids.
# paper_scenarios() returns the mapping table (written to
# analysis/scenarios_paper.csv) and paper_scenario(id) the run scenario.

# Simulation scenario grid of a run
#
# Returns the scenario grid of a run as a list of scenario definitions.
# Each scenario holds DGM overrides for [ccw_sim_params()], the censoring
# model formula used in the analysis, truncation, and the expected diagnostic
# pattern.
#
# * `run = "main"` (default): the scenarios of Parts 1-3
#   (20 scenarios; `f_correct = ~ L + I(L^2) + X + L:X`).
#   - Part 1 (n = 2000): the correct model plus, for each omitted term
#     (`omit_L`: analysis model `~ X`; `omit_Lsq`: DGM without the `L X`
#     term (`beta5 = 0`) and analysis model `~ L + X`, so that exactly
#     `L^2` is omitted; `omit_lag`: `f_correct` while `L_{t-1}` enters the
#     DGM), a 2 x 2
#     factorial of the omitted term's strength in the censoring (treatment)
#     model (`weak`/`strong`: `beta1` 0.4/1.0, `beta2` 0.3/0.6, `beta6`
#     0.6/1.2) and in the outcome model (`none`/`strong`: `gamma1` 0/1.5,
#     `gamma4` 0/1.0, `gamma6` 0/1.2). Names are
#     `<misspec>_c<w|s>_o<n|s>`. The strong outcome coefficients were chosen
#     (pilot: 100 replicates at n = 2000, no bootstrap) as the smallest step
#     giving |bias| >= 2 MC SD in the strong/strong cell (`omit_L`:
#     `gamma1` 1.0 -> 1.5; `omit_Lsq`: `gamma4 = 1`; `omit_lag`: 1.2 as
#     is). `gamma0` stays at -4
#     except where the true g3 horizon risk would leave 0.15-0.40
#     (`true_risk_ccw()`, `n_mc = 1e5`): `omit_Lsq` strong outcome uses
#     `gamma0 = -4.5` (g3 risk 0.375). g3 risks: `correct` 0.219,
#     `omit_L` 0.197 (none) / 0.33 (strong), `omit_Lsq` 0.22 / 0.375,
#     `omit_lag` 0.22-0.23 / 0.33. The manuscript reports `omit_Lsq_c?_os`
#     and `omit_lag_*` from this run (P1-07, P1-09, P1-10 to P1-13); the
#     run `variants` replaces `omit_L_*` and `omit_Lsq_c?_on`, and the run
#     `largerep` replaces `correct`.
#   - Part 2 (correct model; null-reference calibration): `correct_n300`,
#     `correct_n8000` (replaced by the run `largerep`), `near_positivity`
#     (P2-03), `truncation_95` (P2-04).
#   - Part 3: `no_tx_cov_effect` (P3-01), `strong_tx_cov_effect` (P3-02),
#     `extra_ltfu` (P3-03; unweighted informative loss to follow-up).
#
#   Scenarios of this run carry `misspec`, `cens_strength`,
#   `outcome_strength` and `n_boot` (bootstrap resamples).
#
# * `run = "manycov"`: an extended calibration check of
#   the conditional parametric null reference [rrb_null()] (default
#   within-individual coupling of the redrawn censoring indicators;
#   6 scenarios, n = 2000 unless stated, `n_boot = 100`). `many_cov` adds
#   8 extra baseline covariates `X2`-`X9` (`n_extra = 8`, correlation
#   `rho_extra = 0.5` with `X`, default `beta_extra` / `gamma_extra`) that
#   are all in the censoring model
#   `~ L + I(L^2) + X + L:X + X2 + ... + X9`; `many_cov_n300` (P2-05) is the
#   same at n = 300; `rare_censoring` (`beta0 = -3.5`) leaves few censored
#   clones in the no-initiation cells; `extreme_weights` (P2-06;
#   `beta0 = -2`, `beta1 = 2`, `beta2 = 0.5`, `beta3 = 1`) gives very small
#   ESS. Two detection
#   scenarios omit `X2` (`beta_extra[1] = 1`, strong in the censoring model)
#   from the many-covariate model, with `rho_extra = 0.5`
#   (`many_cov_omit_X2`: `X2` partly captured by the modelled `X`, `X3`-`X9`)
#   or `rho_extra = 0` (`many_cov_omit_X2_rho0`). `many_cov`,
#   `rare_censoring` and the two `X2` scenarios are not in the manuscript.
#   Scenarios of this run carry the fields of run main plus `extra_feats`,
#   the extra covariates monitored as history functions (`character(0)` when
#   none; `X2` is monitored even where it is omitted from the model).
#
# * `run = "largerep"` (6 scenarios, `f_correct`, every
#   scenario sets `n`, `n_boot` and the number of simulation replicates
#   `n_rep`). Large-replicate calibration of the null reference of interval
#   RRB and censoring balance under the correct model: `correct_1000`
#   (P1-01; n = 2000, 1000 replicates), `correct_n300_1000` (P2-01; n = 300,
#   1000 replicates), `correct_n8000_500` (P2-02; n = 8000, 500
#   replicates). `early_events` (correct
#   model; `gamma0 = -2`, `gamma2 = -3.5`, `gamma3 = 2.5`; replaced by the
#   run `variants`): many outcome
#   events during g1's grace period and a strongly prognostic `X`, so that
#   the end-of-grace strategy SMD of `X` (g2 vs g1) is clearly non-zero under
#   correct weights through outcome depletion, while the cumulative RRB of
#   the baseline features stays near 0 (pilot: 20 replicates per setting,
#   n = 2000: mean |end-of-grace SMD of X| about 0.25, mean cumulative-RRB
#   MaxIB of `X`, `L_base` about 0.05, true g3 horizon risk about 0.63).
#   `omit_L_ss_trunc`
#   (run main's `omit_L_cs_os`, `beta1 = 1`, `gamma1 = 1.5`, analysis model
#   `~ X`, weights truncated at the 95th percentile) and `omit_L_ss_nearpos`
#   (the same with the near-positivity coefficients `beta0 = -2`,
#   `beta2 = 0.4`, `beta3 = 0.6` of run main's `near_positivity`) are not in
#   the manuscript. Scenarios of this run carry the fields of run manycov
#   plus `n_rep`.
#
# * `run = "variants"`: variants of run main's scenarios (7 scenarios,
#   fields of run largerep, 200 replicates). `omit_L1_*` (P1-02 to P1-05)
#   repeat run main's `omit_L_*` factorial
#   (`beta1` 0.4 / 1, `gamma1` 0 / 1.5, analysis model `~ X`) with
#   `beta2 = beta5 = 0`, so that the model `~ X` omits `L_t` alone.
#   `omit_Lsq_cw_on_g45` (P1-06) and `omit_Lsq_cs_on_g45` (P1-08) repeat run
#   main's `omit_Lsq_*_on` with `gamma0 = -4.5`, the intercept of the strong
#   outcome side, so that all four cells of the `L_t^2` factorial share
#   one intercept. `early_events_b500` (P3-04) repeats run largerep's
#   `early_events` with `n_boot = 500` for the percentile interval of
#   cumulative RRB.
#
# Argument `run`: `"main"`, `"manycov"`, `"largerep"` or `"variants"`.
# Value: A named list; each element has `name`, `description`, `params`
#   (list of overrides), `formula`, `truncate`, `expected`, `run`,
#   `misspec`, `cens_strength`, `outcome_strength` and `n_boot` (runs
#   manycov, largerep, variants: also `extra_feats`; runs largerep,
#   variants: also `n_rep`). Scenarios of the manuscript also carry `id`,
#   `part` and `label` (see paper_scenarios()).
ccw_scenarios <- function(run = c("main", "manycov", "largerep", "variants")) {
  run <- match.arg(run)
  sc <- switch(run,
               main = .scenarios_main(), manycov = .scenarios_manycov(),
               largerep = .scenarios_largerep(),
               variants = .scenarios_variants())
  .annotate_paper(sc, run)
}

.scenarios_main <- function() {
  f_correct <- ~ L + I(L^2) + X + L:X
  sc_d <- function(name, description, params = list(), formula = f_correct,
                   truncate = NULL, expected, misspec = "none",
                   cens_strength = "none", outcome_strength = "none",
                   n_boot = 100) {
    list(name = name, description = description, params = params,
         formula = formula, truncate = truncate, expected = expected,
         run = "main", misspec = misspec,
         cens_strength = cens_strength, outcome_strength = outcome_strength,
         n_boot = n_boot)
  }
  # Part 1: omitted term x censoring strength x outcome strength
  s1_spec <- list(
    omit_L = list(
      formula = ~ X, beta = "b1", cens = c(weak = 0.4, strong = 1.0),
      gamma = "g1", outc = c(none = 0, strong = 1.5), g0 = NULL,
      label = "L_t omitted", expected = "RRB of L up; bias if gamma1 != 0"),
    omit_Lsq = list(
      formula = ~ L + X, beta = "b2", cens = c(weak = 0.3, strong = 0.6),
      gamma = "g4", outc = c(none = 0, strong = 1.0), g0 = -4.5,
      beta_fixed = c(b5 = 0), label = "L_t^2 omitted, beta5 = 0",
      expected = "L_sq / KS up; bias if gamma4 != 0"),
    omit_lag = list(
      formula = f_correct, beta = "b6", cens = c(weak = 0.6, strong = 1.2),
      gamma = "g6", outc = c(none = 0, strong = 1.2), g0 = NULL,
      label = "L_{t-1} omitted",
      expected = "lag / slope features up; bias if gamma6 != 0")
  )
  s1 <- list(sc_d("correct", "Correct specification (base DGM)",
                  expected = "RRB/KS/MHB ~ 0, small bias", n_boot = 500))
  for (ms in names(s1_spec)) {
    sp <- s1_spec[[ms]]
    for (cs in names(sp$cens)) {
      for (os in names(sp$outc)) {
        beta <- c(sp$beta_fixed, stats::setNames(sp$cens[[cs]], sp$beta))
        gamma <- stats::setNames(sp$outc[[os]], sp$gamma)
        desc <- sprintf("%s (%s = %s, %s = %s)", sp$label,
                        sub("b", "beta", sp$beta), sp$cens[[cs]],
                        sub("g", "gamma", sp$gamma), sp$outc[[os]])
        if (os == "strong" && !is.null(sp$g0)) {
          gamma <- c(g0 = sp$g0, gamma)
          desc <- sprintf("%s; gamma0 = %s keeps the true g3 risk <= 0.40",
                          desc, sp$g0)
        }
        s1[[length(s1) + 1]] <- sc_d(
          sprintf("%s_c%s_o%s", ms, substr(cs, 1, 1), substr(os, 1, 1)),
          desc, params = list(beta = beta, gamma = gamma),
          formula = sp$formula, expected = sp$expected, misspec = ms, cens_strength = cs, outcome_strength = os,
          n_boot = 100)
      }
    }
  }
  # Part 2
  s2 <- list(
    sc_d("correct_n300", "Correct specification, n = 300",
         params = list(n = 300), expected = "null reference calibrated",
         n_boot = 500),
    sc_d("correct_n8000", "Correct specification, n = 8000",
         params = list(n = 8000), expected = "null reference calibrated",
         n_boot = 200),
    sc_d("near_positivity", "Extreme adherence-model coefficients",
         params = list(beta = c(b0 = -2, b1 = 1.2, b2 = 0.4, b3 = 0.6)),
         expected = "ESS down, large weights; null reference calibrated",
         n_boot = 500),
    sc_d("truncation_95", "Truncate weights at 95th pct", truncate = 0.95,
         expected = "variance down; interval RRB (untruncated) calibrated",
         n_boot = 500)
  )
  # Part 3
  s3 <- list(
    sc_d("no_tx_cov_effect", "alpha2 = 0", params = list(alpha = c(a2 = 0)),
         expected = "strategy-balance also near 0"),
    sc_d("strong_tx_cov_effect", "alpha2 = -1.5",
         params = list(alpha = c(a2 = -1.5)),
         expected = "strategy SMD of L up but RRB ~ 0, small bias"),
    sc_d("extra_ltfu",
         "Informative LTFU (d0 = -2, d1 = 1.5, gamma1 = 1) not weighted",
         params = list(delta = c(d0 = -2, d1 = 1.5), gamma = c(g1 = 1.0)),
         expected = "RRB ~ 0 but LTFU censoring balance off; bias")
  )
  sc <- c(s1, s2, s3)
  stats::setNames(sc, vapply(sc, `[[`, "", "name"))
}

.scenarios_manycov <- function() {
  f_correct <- ~ L + I(L^2) + X + L:X
  f_many <- ~ L + I(L^2) + X + L:X + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9
  f_omit_x2 <- ~ L + I(L^2) + X + L:X + X3 + X4 + X5 + X6 + X7 + X8 + X9
  extras <- paste0("X", 2:9)
  sc_e <- function(name, description, params = list(), formula = f_correct,
                   expected, misspec = "none", cens_strength = "none",
                   outcome_strength = "none", extra_feats = character(0)) {
    list(name = name, description = description, params = params,
         formula = formula, truncate = NULL, expected = expected,
         run = "manycov", misspec = misspec,
         cens_strength = cens_strength, outcome_strength = outcome_strength,
         n_boot = 100, extra_feats = extra_feats)
  }
  many <- list(n_extra = 8, rho_extra = 0.5)
  beta_x2 <- c(1.0, rep(c(0.2, -0.2), length.out = 7))
  sc <- list(
    sc_e("many_cov",
         "8 extra baseline covariates (rho = 0.5 with X), all modelled",
         params = many, formula = f_many,
         expected = "null reference calibrated with many correlated covariates",
         extra_feats = extras),
    sc_e("many_cov_n300",
         "8 extra baseline covariates (rho = 0.5), all modelled, n = 300",
         params = c(many, list(n = 300)), formula = f_many,
         expected = "null reference calibrated with many correlated covariates",
         extra_feats = extras),
    sc_e("rare_censoring", "Rare artificial censoring (beta0 = -3.5)",
         params = list(beta = c(b0 = -3.5)),
         expected = paste("few censored clones per g3 cell, heavy grace-end",
                          "censoring in g1; null reference calibrated")),
    sc_e("extreme_weights",
         "Extreme weights (beta0 = -2, beta1 = 2, beta2 = 0.5, beta3 = 1)",
         params = list(beta = c(b0 = -2, b1 = 2.0, b2 = 0.5, b3 = 1.0)),
         expected = paste("very small ESS and pi-hat tails; null reference",
                          "calibrated")),
    sc_e("many_cov_omit_X2",
         "X2 (beta_extra1 = 1, rho = 0.5) omitted from the censoring model",
         params = c(many, list(beta_extra = beta_x2)), formula = f_omit_x2,
         expected = paste("omitted X2 partly captured by correlated X;",
                          "detection reduced"),
         misspec = "omit_X2", cens_strength = "strong",
         outcome_strength = "weak", extra_feats = extras),
    sc_e("many_cov_omit_X2_rho0",
         "X2 (beta_extra1 = 1, rho = 0) omitted from the censoring model",
         params = list(n_extra = 8, rho_extra = 0, beta_extra = beta_x2),
         formula = f_omit_x2,
         expected = paste("omitted X2 uncorrelated with modelled covariates;",
                          "detection"),
         misspec = "omit_X2", cens_strength = "strong",
         outcome_strength = "weak", extra_feats = extras)
  )
  stats::setNames(sc, vapply(sc, `[[`, "", "name"))
}

.scenarios_largerep <- function() {
  f_correct <- ~ L + I(L^2) + X + L:X
  sc_f <- function(name, description, params, formula = f_correct,
                   truncate = NULL, expected, misspec = "none",
                   cens_strength = "none", outcome_strength = "none",
                   n_boot, n_rep) {
    list(name = name, description = description, params = params,
         formula = formula, truncate = truncate, expected = expected,
         run = "largerep", misspec = misspec,
         cens_strength = cens_strength, outcome_strength = outcome_strength,
         n_boot = n_boot, extra_feats = character(0), n_rep = n_rep)
  }
  omit_l <- list(beta = c(b1 = 1.0), gamma = c(g1 = 1.5))
  sc <- list(
    sc_f("correct_1000", "Correct specification, n = 2000, 1000 replicates",
         params = list(n = 2000),
         expected = "RRB and CB null references calibrated (precise)",
         n_boot = 500, n_rep = 1000),
    sc_f("correct_n300_1000",
         "Correct specification, n = 300, 1000 replicates",
         params = list(n = 300),
         expected = "RRB and CB null references calibrated at small n",
         n_boot = 500, n_rep = 1000),
    sc_f("correct_n8000_500",
         "Correct specification, n = 8000, 500 replicates",
         params = list(n = 8000),
         expected = "RRB and CB null references calibrated at large n",
         n_boot = 200, n_rep = 500),
    sc_f("early_events",
         paste("Early events in g1's grace period (gamma0 = -2,",
               "gamma2 = -3.5, gamma3 = 2.5), correct model"),
         params = list(n = 2000, gamma = c(g0 = -2, g2 = -3.5, g3 = 2.5)),
         expected = paste("end-of-grace SMD of X non-zero (outcome",
                          "depletion) while cumulative RRB of X ~ 0 and",
                          "interval RRB calibrated; small bias"),
         outcome_strength = "strong", n_boot = 100, n_rep = 200),
    sc_f("omit_L_ss_trunc",
         paste("L_t omitted (beta1 = 1, gamma1 = 1.5), weights truncated",
               "at the 95th pct"),
         params = c(list(n = 2000), omit_l), formula = ~ X, truncate = 0.95,
         expected = "RRB of L up; bias; truncation does not mask RRB",
         misspec = "omit_L", cens_strength = "strong",
         outcome_strength = "strong", n_boot = 100, n_rep = 200),
    sc_f("omit_L_ss_nearpos",
         paste("L_t omitted (beta1 = 1, gamma1 = 1.5) with near-positivity",
               "coefficients (beta0 = -2, beta2 = 0.4, beta3 = 0.6)"),
         params = list(n = 2000,
                       beta = c(b0 = -2, b1 = 1.0, b2 = 0.4, b3 = 0.6),
                       gamma = c(g1 = 1.5)),
         formula = ~ X,
         expected = "RRB of L up with small ESS; bias",
         misspec = "omit_L", cens_strength = "strong",
         outcome_strength = "strong", n_boot = 100, n_rep = 200)
  )
  stats::setNames(sc, vapply(sc, `[[`, "", "name"))
}

.scenarios_variants <- function() {
  f_correct <- ~ L + I(L^2) + X + L:X
  sc_h <- function(name, description, params, formula = f_correct,
                   expected, misspec = "none", cens_strength = "none",
                   outcome_strength = "none", n_boot = 100) {
    list(name = name, description = description, params = params,
         formula = formula, truncate = NULL, expected = expected,
         run = "variants", misspec = misspec,
         cens_strength = cens_strength, outcome_strength = outcome_strength,
         n_boot = n_boot, extra_feats = character(0), n_rep = 200)
  }
  sc <- list()
  for (cs in c("weak", "strong")) {
    for (os in c("none", "strong")) {
      b1 <- c(weak = 0.4, strong = 1.0)[[cs]]
      g1 <- c(none = 0, strong = 1.5)[[os]]
      sc[[length(sc) + 1]] <- sc_h(
        sprintf("omit_L1_c%s_o%s", substr(cs, 1, 1), substr(os, 1, 1)),
        sprintf(paste("L_t alone omitted (beta1 = %s, gamma1 = %s;",
                      "beta2 = beta5 = 0)"), b1, g1),
        params = list(beta = c(b1 = b1, b2 = 0, b5 = 0), gamma = c(g1 = g1)),
        formula = ~ X, expected = "RRB of L up; bias if gamma1 != 0",
        misspec = "omit_L1", cens_strength = cs,
        outcome_strength = os)
    }
  }
  for (cs in c("weak", "strong")) {
    b2 <- c(weak = 0.3, strong = 0.6)[[cs]]
    sc[[length(sc) + 1]] <- sc_h(
      sprintf("omit_Lsq_c%s_on_g45", substr(cs, 1, 1)),
      sprintf(paste("L_t^2 omitted, beta5 = 0 (beta2 = %s, gamma4 = 0);",
                    "gamma0 = -4.5"), b2),
      params = list(beta = c(b5 = 0, b2 = b2), gamma = c(g0 = -4.5, g4 = 0)),
      formula = ~ L + X, expected = "L_sq / KS up; no bias",
      misspec = "omit_Lsq", cens_strength = cs,
      outcome_strength = "none")
  }
  sc[[length(sc) + 1]] <- sc_h(
    "early_events_b500",
    paste("Early events in g1's grace period (gamma0 = -2, gamma2 = -3.5,",
          "gamma3 = 2.5), correct model, 500 bootstrap resamples"),
    params = list(n = 2000, gamma = c(g0 = -2, g2 = -3.5, g3 = 2.5)),
    expected = "as largerep early_events, more precise percentile interval",
    outcome_strength = "strong", n_boot = 500)
  stats::setNames(sc, vapply(sc, `[[`, "", "name"))
}

# Flatten a scenario list to a data frame (e.g. for `scenarios_<run>.csv`)
#
# Argument `scenarios`: Output of [ccw_scenarios()].
# Value: A data frame with one row per scenario. `run` is the run;
#   `id`, `part` and `label` identify the scenarios of the manuscript (`NA`
#   for the others, see paper_scenarios()); `n` is the scenario's sample
#   size; `n_extra` is the number of extra baseline
#   covariates (0 when absent) and `extra_feats` the comma-separated extra
#   covariates monitored as history functions (`""` when none); `n_rep` is
#   the scenario's number of simulation replicates when the scenario sets
#   it (runs largerep and variants; `NA` otherwise).
scenario_table <- function(scenarios = ccw_scenarios()) {
  chr_or_na <- function(x) if (is.null(x)) NA_character_ else x
  flat <- function(p) {
    if (length(p) == 0) return("")
    paste(vapply(names(p), function(nm) {
      v <- p[[nm]]
      val <- if (!is.null(names(v))) {
        paste(paste0(names(v), "=", v), collapse = ";")
      } else {
        paste(v, collapse = ";")
      }
      paste0(nm, ":{", val, "}")
    }, ""), collapse = " ")
  }
  do.call(rbind, lapply(scenarios, function(s) {
    data.frame(
      scenario = s$name,
      description = s$description,
      param_overrides = flat(s$params),
      censoring_formula = paste(deparse(s$formula), collapse = ""),
      truncate = if (is.null(s$truncate)) NA_real_ else s$truncate,
      expected = s$expected,
      run = chr_or_na(s$run),
      id = chr_or_na(s$id),
      part = if (is.null(s$part)) NA_integer_ else s$part,
      label = chr_or_na(s$label),
      misspec = chr_or_na(s$misspec),
      cens_strength = chr_or_na(s$cens_strength),
      outcome_strength = chr_or_na(s$outcome_strength),
      n = if (is.null(s$params[["n"]])) ccw_sim_params()$n else
        s$params[["n"]],
      n_boot = if (is.null(s$n_boot)) NA_real_ else s$n_boot,
      n_extra = if (is.null(s$params$n_extra)) 0 else s$params$n_extra,
      extra_feats = paste(s$extra_feats, collapse = ","),
      n_rep = if (is.null(s$n_rep)) NA_real_ else s$n_rep
    )
  }))
}

# Scenarios of the manuscript, in the row order of Table tab:grid of the
# manuscript: id, part, label (manuscript wording, ASCII),
# the run and run scenario that supply the reported results, the factor
# levels of Part 1 (NA elsewhere: taken from the run scenario) and the
# "<run>:<scenario>" keys of the rows of earlier runs that the run scenario
# replaces.
.paper_map <- function() {
  fac <- c("weak / none", "weak / strong", "strong / none", "strong / strong")
  cells <- c("cw_on", "cw_os", "cs_on", "cs_os")
  p1 <- data.frame(
    id = sprintf("P1-%02d", 1:13), part = 1L,
    label = c("Correct specification, n = 2000",
              paste0("Omitted L_k, ", fac), paste0("Omitted L_k^2, ", fac),
              paste0("Omitted L_{k-1}, ", fac)),
    run = c("largerep", rep("variants", 4),
            "variants", "main", "variants", "main", rep("main", 4)),
    scenario = c("correct_1000", paste0("omit_L1_", cells),
                 "omit_Lsq_cw_on_g45", "omit_Lsq_cw_os",
                 "omit_Lsq_cs_on_g45", "omit_Lsq_cs_os",
                 paste0("omit_lag_", cells)),
    misspec = c("none", rep(c("omit_L", "omit_Lsq", "omit_lag"), each = 4)),
    cens_strength = c("none", rep(rep(c("weak", "strong"), each = 2), 3)),
    outcome_strength = c("none", rep(c("none", "strong"), 6)),
    supersedes = c("main:correct", paste0("main:omit_L_", cells),
                   "main:omit_Lsq_cw_on", "", "main:omit_Lsq_cs_on", "",
                   rep("", 4))
  )
  p23 <- data.frame(
    id = c(sprintf("P2-%02d", 1:6), sprintf("P3-%02d", 1:4)),
    part = rep(2:3, c(6, 4)),
    label = c("Correct specification, n = 300",
              "Correct specification, n = 8000",
              "Near-violation of positivity",
              "Weight truncation, 95th percentile",
              "Many correlated covariates", "Extreme weights",
              "No treatment->covariate effect",
              "Strong treatment->covariate effect",
              "Unweighted loss to follow-up", "Many early events"),
    run = c("largerep", "largerep", "main", "main", "manycov", "manycov",
            "main", "main", "main", "variants"),
    scenario = c("correct_n300_1000", "correct_n8000_500", "near_positivity",
                 "truncation_95", "many_cov_n300", "extreme_weights",
                 "no_tx_cov_effect", "strong_tx_cov_effect", "extra_ltfu",
                 "early_events_b500"),
    misspec = NA_character_, cens_strength = NA_character_,
    outcome_strength = NA_character_,
    supersedes = c("main:correct_n300", "main:correct_n8000", rep("", 7),
                   "largerep:early_events")
  )
  rbind(p1, p23)
}

# Add id / part / label to the scenarios of a run that the manuscript
# reports
.annotate_paper <- function(sc, run) {
  map <- .paper_map()
  map <- map[map$run == run, ]
  for (k in seq_len(nrow(map))) {
    nm <- map$scenario[k]
    stopifnot(nm %in% names(sc))
    sc[[nm]]$id <- map$id[k]
    sc[[nm]]$part <- map$part[k]
    sc[[nm]]$label <- map$label[k]
  }
  sc
}

# Scenarios of the manuscript and the runs that supply them
#
# Value: A data frame with one row per scenario of the manuscript, in the
#   row order of Table tab:grid (P1-01 to P1-13, P2-01 to P2-06, P3-01 to
#   P3-04; E ids would follow, see the header): `id`, `part`, `label`,
#   `run` and `scenario` (the run scenario whose results the manuscript
#   reports), `supersedes` ("<run>:<scenario>" keys of the rows of earlier
#   runs that the run scenario replaces, "" if none), `misspec`,
#   `cens_strength`, `outcome_strength` (Part 1 factor levels; for the
#   omitted-L_k cells "omit_L", while the run scenario says "omit_L1") and
#   the scenario spec needed to rerun it: `description`, `param_overrides`,
#   `censoring_formula`, `truncate`, `n`, `n_rep` (replicates; the run
#   default of 200, or 100 at n >= 8000, unless the scenario sets it),
#   `n_boot`, `n_extra`, `extra_feats` and `expected`.
paper_scenarios <- function() {
  map <- .paper_map()
  runs <- unique(map$run)
  tabs <- stats::setNames(lapply(runs, function(r) {
    scenario_table(ccw_scenarios(r))
  }), runs)
  spec <- do.call(rbind, lapply(seq_len(nrow(map)), function(k) {
    tb <- tabs[[map$run[k]]]
    tb[tb$scenario == map$scenario[k], ]
  }))
  fill <- function(a, b) ifelse(is.na(a), b, a)
  n_rep <- ifelse(is.na(spec$n_rep), ifelse(spec$n >= 8000, 100, 200),
                  spec$n_rep)
  out <- data.frame(
    id = map$id, part = map$part, label = map$label, run = map$run,
    scenario = map$scenario, supersedes = map$supersedes,
    misspec = fill(map$misspec, spec$misspec),
    cens_strength = fill(map$cens_strength, spec$cens_strength),
    outcome_strength = fill(map$outcome_strength, spec$outcome_strength),
    description = spec$description, param_overrides = spec$param_overrides,
    censoring_formula = spec$censoring_formula, truncate = spec$truncate,
    n = spec$n, n_rep = n_rep, n_boot = spec$n_boot, n_extra = spec$n_extra,
    extra_feats = spec$extra_feats, expected = spec$expected
  )
  rownames(out) <- NULL
  out
}

# Run scenario (element of ccw_scenarios(run)) of a manuscript scenario id
paper_scenario <- function(id) {
  map <- .paper_map()
  k <- match(id, map$id)
  if (is.na(k)) stop("Unknown scenario id: ", id)
  ccw_scenarios(map$run[k])[[map$scenario[k]]]
}
