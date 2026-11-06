# sim_scenarios.R --------------------------------------------------------------
# Scenario grids of the simulation study (plain functions, not part of the
# package). ccw_scenarios(run) returns the grid of one run of
# 02_simulation_study.R (main, manycov, largerep, variants);
# paper_scenarios() maps the manuscript ids P1-01 to P3-04 (Table 2, row
# order) to the run scenario that supplies each (saved by 02b as
# analysis/scenarios_paper.csv). Usage:
#   pkgload::load_all("."); source("analysis/sim_scenarios.R")

# Scenario grid of a run: a named list of scenarios, each with `name`,
# `description`, `params` (overrides of ccw_sim_params()), `formula`
# (censoring model), `truncate`, `expected` and run-specific fields
# (`n_boot`, `extra_feats`, `n_rep`; `id`, `part`, `label` for the
# scenarios of the manuscript).
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

# Flatten a scenario list to a data frame (one row per scenario)
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

# Manuscript scenarios in Table 2 row order: id, part, label, the run and
# run scenario that supply the results, Part 1 factor levels and
# `supersedes` ("<run>:<scenario>" keys of the rows it replaces)
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

# Add id / part / label to the manuscript scenarios of a run
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

# Manuscript scenarios and the runs that supply them: one row per id with
# the run scenario and the spec needed to rerun it
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
