# Shared small fixtures ------------------------------------------------------
small_weighted <- function(n = 400, seed = 11, formula = ~ L + I(L^2) + X + L:X,
                           ...) {
  d <- simulate_ccw(ccw_sim_params(n = n), seed = seed)
  d <- make_history_features(d, interactions = list(c("L", "X")))
  fit_censoring_model(d, formula, ...)
}

# Deterministic toy data for which the IPW identity holds exactly:
# stratum f = 1: 4 rows, 2 censored, pi = 0.5
# stratum f = 3: 5 rows, 1 censored, pi = 0.8
toy_ipw <- function() {
  data.frame(
    id = 1:9, strategy = "g", time = 0,
    f = c(rep(1, 4), rep(3, 5)),
    censored_next = c(1, 1, 0, 0, 1, 0, 0, 0, 0),
    pi_hat = c(rep(0.5, 4), rep(0.8, 5)),
    q_hat = 0.7,
    w_prev_unstab = 1, w_prev_stab = 1
  )
}

# Two-interval toy for the cumulative IPW identity.
# t = 0: 10 clones (5 with f = 1, 5 with f = 3), no censoring possible;
#        one f = 1 clone has the outcome at t = 0 (carried forward).
# t = 1: f = 1: 4 rows, 2 censored, pi = 0.5; f = 3: 5 rows, 1 censored,
#        pi = 0.8. Numerator q = 0.7 at t = 1 (common within the cell).
toy_cum <- function() {
  t0 <- data.frame(id = 1:10, strategy = "g", time = 0,
                   f = c(rep(1, 5), rep(3, 5)), censored_next = 0,
                   event_next = c(1, rep(0, 9)), pi_hat = 1, q_hat = 1,
                   w_prev_unstab = 1, w_prev_stab = 1)
  t1 <- data.frame(id = 2:10, strategy = "g", time = 1,
                   f = c(rep(1, 4), rep(3, 5)),
                   censored_next = c(1, 1, 0, 0, 1, 0, 0, 0, 0),
                   event_next = 0,
                   pi_hat = c(rep(0.5, 4), rep(0.8, 5)), q_hat = 0.7,
                   w_prev_unstab = 1, w_prev_stab = 1)
  t1$event_next[t1$censored_next == 1] <- NA
  rbind(t0, t1)
}
