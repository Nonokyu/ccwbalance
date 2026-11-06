# 02j_depletion_incidence.R ----------------------------------------------------
# Outcome incidence during the grace periods (cited in the text) for the base
# mechanism (P1-01) and early_events (P3-04), and for the Webster-Clark
# SynPUF cohort (crude and Kaplan-Meier).
# Inputs:  paper_scenario(); data-raw/webster_clark/Data/baseccwcohort.sas7bdat
# Outputs: output/sim/sim_depletion_incidence.csv,
#          output/empirical/wc_grace_incidence.csv
# Run:     Rscript analysis/02j_depletion_incidence.R

pkgload::load_all(".", quiet = TRUE)
source(file.path("analysis", "sim_scenarios.R"))

n_obs <- 2e5
n_regime <- 5e5
ids <- c(correct_1000 = "P1-01", early_events = "P3-04")

sim <- do.call(rbind, lapply(names(ids), function(nm) {
  pr <- paper_scenario(ids[[nm]])$params
  pr$n <- n_obs
  p <- do.call(ccw_sim_params, pr)
  set.seed(20261025)
  s <- .sim_paths(p)
  tau <- .first_one(s$A)
  ev <- .first_one(s$Y)
  r <- data.frame(
    scenario = nm,
    ci_int2 = mean(ev <= 2), ci_int8 = mean(ev <= 8),
    ci_int11 = mean(ev <= 11),
    ci_int2_noinit = mean(ev[tau > 2] <= 2),
    ci_int8_noinit = mean(ev[tau > 8] <= 8),
    ev_before_init_int01 = mean(ev <= 1 & tau > ev),
    init_by_int2 = mean(tau <= 2), init_by_int8 = mean(tau <= 8)
  )
  set.seed(20261026)
  pm <- p
  pm$n <- n_regime
  for (g in c("g1", "g2", "g3")) {
    eg <- .first_one(.sim_paths(pm, regime = g)$Y)
    r[[paste0(g, "_risk_int2")]] <- mean(eg <= 2)
    r[[paste0(g, "_risk_int8")]] <- mean(eg <= 8)
    r[[paste0(g, "_risk_end")]] <- mean(eg <= 11)
  }
  r
}))
utils::write.csv(sim, "output/sim/sim_depletion_incidence.csv",
                 row.names = FALSE)
print(t(sim), digits = 3)

raw <- haven::read_sas(file.path("data-raw", "webster_clark", "Data",
                                 "baseccwcohort.sas7bdat"))
fu <- raw$followup
y <- raw$outcome
rx <- raw$RX_start
km <- survival::survfit(survival::Surv(fu, y) ~ 1)
km_at <- function(d) 1 - summary(km, times = d, extend = TRUE)$surv
wc_row <- function(d) {
  ev <- y == 1 & fu < d            # outcome on days 0 to d - 1
  noinit <- is.na(rx) | rx >= d    # not initiated by day d - 1
  data.frame(day = d, n = nrow(raw), n_event = sum(ev), crude = mean(ev),
             km = km_at(d - 0.5), n_noinit = sum(noinit),
             crude_noinit = mean(ev[noinit]),
             prop_event_before_init = mean(ev & (is.na(rx) | rx > fu)))
}
wc <- rbind(wc_row(30), wc_row(90), wc_row(181))
wc$km[3] <- km_at(180)
dir.create("output/empirical", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(wc, "output/empirical/wc_grace_incidence.csv",
                 row.names = FALSE)
print(wc, digits = 3)
