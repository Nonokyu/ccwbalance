# 03b_wc_null_reference.R ----------------------------------------------------
# Conditional parametric null reference for the interval RRB / MaxIB
# (rrb_null()) applied to the Webster-Clark SynPUF example (analysis/03).
#
# Data preparation (analysis/wc_data.R) and the base censoring model (age,
# sex, renal) are identical to analysis/03_empirical_webster_clark.R; the
# null reference covers the eleven history functions (five baseline
# covariates, the three time-varying claims indicators and the three depth
# functions fill_lag, fill_30, fill_cum, which enter no censoring model).
# Two maxima per strategy x function: the raw MaxIB (max |RRB|, $maxib) and
# the studentised maximum (max |RRB| / null SD of the cell, $maxib_stud),
# each against the null 95th percentile of the same maximum.
#
# Outputs (git-ignored):
#   output/empirical/wc_rrb_null_cells.csv, wc_rrb_null_maxib.csv,
#   wc_rrb_null_maxib_stud.csv
#
# Run from the project root:  Rscript analysis/03b_wc_null_reference.R
#   NULL_DRAWS (default 500)

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("analysis", "wc_data.R"))

N_DRAWS <- as.numeric(Sys.getenv("NULL_DRAWS", "500"))
out_dir <- file.path("output", "empirical")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}
t_start <- Sys.time()

# 1. Webster-Clark data and base censoring model (as in analysis/03) ----------
wc <- wc_prepare(delta = 10, horizon = 180)
feats <- wc$feats
w <- fit_censoring_model(wc$long, wc$formulas$base, stabilized = TRUE)

# 2. Null reference ------------------------------------------------------------
t0 <- Sys.time()
r0 <- rrb_null(w, feats, n_draws = N_DRAWS, seed = 20260925)
t_wc <- difftime(Sys.time(), t0, units = "secs")
save_csv(r0$cells, "wc_rrb_null_cells.csv")
save_csv(r0$maxib, "wc_rrb_null_maxib.csv")
save_csv(r0$maxib_stud, "wc_rrb_null_maxib_stud.csv")

show_cells <- function(r, g, tt) {
  x <- r$cells[as.character(r$cells$strategy) == g & r$cells$time == tt,
               c("strategy", "time", "feature", "n_risk", "n_censored",
                 "rrb_obs", "null_sd", "z", "p_two")]
  print(x, digits = 3, row.names = FALSE)
}
frac_z <- function(r) {
  z <- r$cells$z[!is.na(r$cells$z)]
  sprintf("%d / %d = %.3f", sum(abs(z) > 1.96), length(z),
          mean(abs(z) > 1.96))
}

cat(sprintf("\n=== Webster-Clark: rrb_null(), %d draws (%.1f s) ===\n",
            N_DRAWS, as.numeric(t_wc)))
cat("\nMaxIB of fill, the depth functions, race_other and all features",
    "(observed vs null):\n")
print(r0$maxib[r0$maxib$feature %in% c("fill", wc$depth, "race_other",
                                       "(all)"), ],
      digits = 3, row.names = FALSE)
cat("\nStudentised maximum max|z| of the same functions and outpt",
    "(observed vs null):\n")
print(r0$maxib_stud[r0$maxib_stud$feature %in%
                      c("fill", wc$depth, "race_other", "outpt", "(all)"), ],
      digits = 3, row.names = FALSE)
cat("\nInterval ending day 90 (time 8), strategy 0_90:\n")
show_cells(r0, "0_90", 8)
cat("\nInterval ending day 30 (time 2), strategy 0_30:\n")
show_cells(r0, "0_30", 2)
cat("\nMaxIB (observed vs conditional null):\n")
print(r0$maxib, digits = 3, row.names = FALSE)
cat("\nCells (all strategies/features, with censoring) with |z| > 1.96:",
    frac_z(r0), "\n")
cat("\nAll cells with |z| > 2:\n")
big <- r0$cells[!is.na(r0$cells$z) & abs(r0$cells$z) > 2, ]
print(big[, c("strategy", "time", "feature", "n_risk", "n_censored",
              "rrb_obs", "null_sd", "z", "p_two")],
      digits = 3, row.names = FALSE)

cat(sprintf("\nTotal runtime: %.1f s\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
