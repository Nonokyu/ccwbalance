# ccwbalance

ccwbalance provides diagnostics for the artificial-censoring weights of
clone-censor-weight (CCW) analyses in target trial emulation. Its main
diagnostic is risk-set reconstruction balance (RRB). Interval RRB checks
whether the retained individuals, reweighted by the inverse probability of
remaining uncensored, reproduce the risk set just before artificial
censoring. Cumulative RRB checks whether they reproduce the baseline
population. The package also computes censoring balance, a conditional
parametric null reference that judges each strategy by its studentised
maximum imbalance, between-strategy standardised mean differences, weight
summaries and plots.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Nonokyu/ccwbalance")
```

## Data contract

All diagnostics take long data with one row per individual x strategy x
interval in which the individual is in the risk set of the strategy.
`ccw_clone()` builds such data from one row per individual (initiation-window
strategies) or from person-period data (sustained strategies), and
`simulate_ccw()` simulates it. Column names can be changed through function
arguments.

| Column | Content |
|---|---|
| `id` | Individual identifier |
| `strategy` | Strategy label |
| `time` | Interval index `k = 0, 1, ...` |
| `censored_next` | `D_k`: 1 if artificially censored at the end of interval `k` |
| `censor_possible` | Whether artificial censoring can occur in the row |
| `phase` | Kind of censoring possible (`"no_init"`, `"grace_end"`, `"sustain"`, `"none"`); defines the censoring-model groups with `strategy` |
| `event_next` | Outcome in interval `k` (`NA` when censored or lost to follow-up) |
| `ltfu_next` | Optional: loss to follow-up in interval `k` |
| `A` | Optional: treatment `A_k`, used by `rrb_null()` |
| `X`, `L`, ... | History features `f(H_k)` |

`fit_censoring_model()` adds the columns `pi_hat`, `q_hat`, `w_prev_unstab`,
`w_prev_stab`, `w_ref` and `w_post`.

## Example

```r
library(ccwbalance)

# simulated data with three strategies g1, g2, g3
d <- simulate_ccw(ccw_sim_params(n = 1000), seed = 1)
d <- make_history_features(d, interactions = list(c("L", "X")))
feats <- c("X", "L", "L_sq", "L_lag", "L_x_X")

# censoring model by strategy x phase; stabilised weights
w <- fit_censoring_model(d, ~ L + I(L^2) + X + L:X)

# interval RRB, unweighted and weighted
r <- rrb(w, feats, method = c("unweighted", "weighted"))
r_cens <- r[r$n_censored > 0, ]
aggregate_balance(r_cens, by = c("method", "feature"))

# null reference for interval RRB and censoring balance
nr <- rrb_null(w, feats, n_draws = 50, seed = 1)
# per stratum x history function: null SD and z = RRB / null SD
head(nr$cells[, c("strategy", "time", "feature", "rrb_obs", "null_sd", "z",
                  "cb_cens_z")])
# per strategy: studentised maximum (largest |z| over its strata and
# history functions) against its null 95th percentile
st <- nr$maxib_stud[nr$maxib_stud$feature == "(all)" &
                      nr$maxib_stud$strategy != "(all)", ]
st$exceeds <- st$stud_obs > st$null_q95
st
# the same for censoring balance over the censorable individuals
cb_st <- nr$maxib_cb_cens_stud
cb_st[cb_st$feature == "(all)" & cb_st$strategy != "(all)", ]
# raw maximum |RRB| (MaxIB), descriptive
nr$maxib[nr$maxib$feature == "(all)", ]

# cumulative RRB of baseline-fixed features, with bootstrap intervals
rc <- rrb(w, c("X", "L_base"), scope = "cumulative")
bt <- rrb_boot(w, c("X", "L_base"), B = 20, seed = 1)
head(bt)

# censoring balance, end-of-grace SMD, weight support
cb <- censoring_balance(w, feats, weighted = TRUE)
sb <- strategy_balance(w, "X", reference = "g1")
sb[sb$time == 3, ]
ws <- weight_support(w)
head(ws[, c("strategy", "time", "ess_ratio", "w99", "pi_min")])

# plots
plot_love_dynamic(r)
plot_weight_support(ws)
```

## Functions

| Purpose | Functions |
|---|---|
| Data construction | `ccw_clone()`, `ccw_sustained()`, `make_history_features()` |
| Weights and risk | `fit_censoring_model()`, `estimate_risk()` |
| Reconstruction diagnostics | `rrb()`, `ks_rr()`, `mahalanobis_rr()`, `aggregate_balance()` |
| Null reference and uncertainty | `rrb_null()`, `rrb_boot()` |
| Other balance and weight summaries | `censoring_balance()`, `strategy_balance()`, `weight_support()` |
| Plots | `plot_love_dynamic()`, `plot_balance_heatmap()`, `plot_mahalanobis()`, `plot_weight_support()`, `plot_ecdf_pre_post()`, `plot_dashboard()` |
| Simulation | `ccw_sim_params()`, `simulate_ccw()`, `true_risk_ccw()` |

## Notation

| Symbol | Meaning | Package |
|---|---|---|
| `k` | Interval index | column `time` |
| `D_k` | Artificial censoring at the end of interval `k` | column `censored_next` |
| `pi_k = P(D_k = 0 \| R_k = 1, H_k)` | Probability of remaining uncensored | column `pi_hat` |
| `q_k` | Stabilising numerator, common to the stratum | column `q_hat` |
| `W_{k-1} q_k` | Reference weight of interval RRB | column `w_ref` |
| `W_k`, `W^u_k` | Stabilised and unstabilised cumulative weights | columns `weight_stab`, `weight_unstab`; `w_post` holds the chosen, possibly truncated, weight |
| stratum | Strategy x interval | `strategy`, `time` |
| `f(H_k)` | History function | argument `features` |
| `RRB_{j,g}(k)` | Interval or cumulative RRB | `rrb(scope = )` |
| `CB_{j,g}(k)` | Censoring balance | `censoring_balance(weighted = TRUE)`, `rrb_null()$cells$cb_obs` |
| `z = RRB / SD_0` | RRB of a stratum in units of its null SD | `rrb_null()$cells$z` (`cb_cens_z` for censoring balance) |
| studentised maximum | Largest `\|z\|` of a strategy over its strata and history functions (primary) | `rrb_null()$maxib_stud`, `$maxib_cb_cens_stud` |
| MaxIB | Maximum absolute imbalance over strata (descriptive) | `aggregate_balance()`, `rrb_null()$maxib` |

## Zero targets

* Interval RRB has target 0 for every history function `f(H_k)`, including
  post-baseline covariates, when the censoring model is correct and
  positivity holds. `rrb_null()` gives its finite-sample null reference:
  for every stratum the null SD and `z = RRB / SD_0` of interval RRB and of
  censoring balance (over all rows and over the censorable individuals),
  and for every strategy the studentised maximum `max |z|` over its strata
  and history functions with its null 95th percentile and exceedance
  probability `p`. The studentised maximum is the primary summary; the raw
  maximum `max |RRB|` (MaxIB) is reported as a descriptive summary.
* Cumulative RRB has target 0 only for baseline-fixed features (`X`,
  `<v>_base`). For post-baseline covariates it is descriptive.
* Between-strategy balance of post-baseline covariates, including the
  end-of-grace SMD, is descriptive and has no zero target.

`aggregate_balance()` summarises all rows it is given; filter to the strata
with artificial censoring (`n_censored > 0`) to obtain MaxIB over the same
strata as `rrb_null()`.

## Reproducing the manuscript

The folder `analysis/` holds the simulation study, the supplementary checks
and the case study of the manuscript; `data-raw/get_synpuf_claims.sh`
extracts the claims of the case-study cohort; `output/` holds the saved
simulation outputs, so the tables and figures can be re-created without rerunning the
simulations. None of these folders is part of the installed package. See
[`analysis/README.md`](analysis/README.md) for the required packages, the
data, the run order and the outputs behind each table and figure.

## Citation

Nonomiya Y. Risk-set reconstruction balance: zero-targeted diagnostics for
artificial censoring in clone–censor–weight analyses. In preparation.

## License

MIT. See `LICENSE.md`.
