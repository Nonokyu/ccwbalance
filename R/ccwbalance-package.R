#' ccwbalance: artificial-censoring diagnostics for clone-censor-weight
#' analyses
#'
#' Diagnostics for the artificial-censoring weights of clone-censor-weight
#' (CCW) analyses: interval and cumulative risk-set reconstruction balance
#' (RRB) with distributional (KS) and multivariate (Mahalanobis) versions,
#' Jackson-type censoring balance, a conditional parametric null reference
#' that judges each strategy by the studentised maximum (the largest
#' |RRB| / null SD over its strata and history functions) as well as the
#' raw maximum absolute imbalance (MaxIB), an individual-level bootstrap, between-strategy (including end-of-grace) standardised mean
#' differences, weight-support summaries, plots, tools to build the input
#' data and a simulation data-generating mechanism.
#'
#' ## Notation
#' Follow-up is divided into intervals `I_k`, `k = 0, 1, ..., K - 1`
#' (column `time`). `D_k` indicates artificial censoring at the end of
#' interval `k` (column `censored_next`), `R_k` membership of the risk set
#' at the start of interval `k`, and `H_k` the history up to `L_k` and
#' `A_{k-1}`. `pi_k = P(D_k = 0 | R_k = 1, H_k)` is the probability of
#' remaining uncensored (column `pi_hat`) and `q_k` the stabilising
#' numerator, common to all individuals in a stratum (column `q_hat`). A
#' stratum is a strategy x interval combination. The cumulative weight is
#' `W_k = prod_{l <= k} q_l / pi_l = Q_k W^u_k`, with `Q_k = prod_{l <= k}
#' q_l` and the unstabilised weight `W^u_k = prod_{l <= k} 1 / pi_l`.
#'
#' ## Data contract
#' The diagnostic functions take long data with one row per individual x
#' strategy x interval in which the individual is in the risk set of the
#' strategy (`R_k = 1`). Column names can be changed through the function
#' arguments. Required columns:
#' * `id`: individual identifier;
#' * `strategy`: strategy label (factor or character);
#' * `time`: interval index `k` (`0, 1, ...`);
#' * `censored_next`: `D_k` (1 if artificially censored at the end of
#'   interval `k`, else 0);
#' * `censor_possible`: whether artificial censoring can occur in the row;
#' * `phase`: kind of artificial censoring possible in the row (e.g.
#'   `"no_init"`, `"grace_end"`, `"sustain"`, `"none"`); together with
#'   `strategy` it defines the censoring-model groups;
#' * `event_next`: outcome in interval `k` (`NA` when censored);
#' * history features `f(H_k)`, e.g. `X`, `L`, `L_lag`, and those added by
#'   [make_history_features()].
#'
#' Optional columns: `ltfu_next` (loss to follow-up in interval `k`) and
#' `A` (treatment `A_k`, used to orient the coupling of [rrb_null()]).
#' [fit_censoring_model()] adds the weight columns `pi_hat`, `q_hat`,
#' `w_prev_unstab`, `w_prev_stab`, `w_ref` and `w_post`. [ccw_clone()]
#' builds data in this format from one row per individual or from
#' person-period data, and [simulate_ccw()] simulates it.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data
NULL

utils::globalVariables(c(".ref_w", ".post_w", ".value", ".feature", ".grp"))
