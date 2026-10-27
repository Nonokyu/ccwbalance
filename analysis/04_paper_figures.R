# 04_paper_figures.R ----------------------------------------------------------
# Case-study figure of the manuscript (Figure 4, fig4_wc_rrb_profile.png),
# written to output/figures/. Reads only saved outputs. The simulation
# figures of the manuscript (Figure 2, sim_part1_factorial.png; Figure 3,
# sim_part2_calibration.png) are written to output/figures/ by
# analysis/02b_simulation_tables_paper.R.
#
# Webster-Clark SynPUF example (fig4_wc_rrb_profile.png): weighted interval
# and cumulative risk-set reconstruction balance (RRB) of the baseline
# covariates over follow-up (source: output/empirical/wc_rrb.csv, written by
# analysis/03_empirical_webster_clark.R). When the individual-bootstrap
# files output/empirical/wc_rrb_boot_interval.csv and
# wc_rrb_boot_cumulative.csv exist (analysis/03c_wc_sensitivity.R) and the
# environment variable CCW_FIG4_BANDS is "1", 95% percentile bootstrap bands
# are added for race_other (grey ribbons); by default (CCW_FIG4_BANDS unset)
# the figure has no bands, as in the manuscript.
#
# Run from the project root, after analysis/03_empirical_webster_clark.R:
#   Rscript analysis/04_paper_figures.R

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(ggplot2)

fig_dir <- file.path("output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Webster-Clark RRB profile ----------------------------------------------------
delta <- 10 # days per interval (analysis/03_empirical_webster_clark.R)
rrb_file <- file.path("output", "empirical", "wc_rrb.csv")
if (!file.exists(rrb_file)) {
  stop("Missing ", rrb_file, "; run analysis/03_empirical_webster_clark.R.")
}
wc <- read.csv(rrb_file)
wc <- wc[wc$method == "weighted", ]
# each strategy up to its last interval with artificial censoring (end of
# its grace period); the weights do not change afterwards
cens <- wc[wc$scope == "interval" & wc$n_censored > 0, ]
last_cens <- tapply(cens$time, cens$strategy, max)
wc <- wc[wc$time <= last_cens[as.character(wc$strategy)], ]

strategy_labels <- c(
  "0_30" = "Initiate days 0-29", "0_90" = "Initiate days 0-89",
  "30_90" = "Initiate days 30-89", "no_0_90" = "No initiation days 0-89"
)
wc$day <- (wc$time + 1) * delta # end of interval, where censoring occurs
wc$strategy_lab <- factor(strategy_labels[wc$strategy],
                          levels = strategy_labels)
wc$group <- ifelse(wc$feature == "race_other",
                   "Race: other (not in censoring model)",
                   "Other baseline covariates")
wc$group <- factor(wc$group, levels = c("Race: other (not in censoring model)",
                                        "Other baseline covariates"))
wc$scope_lab <- factor(
  ifelse(wc$scope == "interval",
         "A. Interval RRB (pre-censoring risk set)",
         "B. Cumulative RRB (baseline cloned population)"),
  levels = c("A. Interval RRB (pre-censoring risk set)",
             "B. Cumulative RRB (baseline cloned population)")
)
wc$line_id <- interaction(wc$strategy, wc$feature, wc$scope)

# 95% percentile bootstrap bands (race_other, strategies 0_90 and 30_90),
# matched to the weighted rows by strategy x time x feature x scope
boot_files <- file.path("output", "empirical",
                        sprintf("wc_rrb_boot_%s.csv",
                                c("interval", "cumulative")))
band <- NULL
show_bands <- identical(Sys.getenv("CCW_FIG4_BANDS", "0"), "1")
if (show_bands && all(file.exists(boot_files))) {
  bt <- do.call(rbind, lapply(boot_files, function(f) {
    read.csv(f)[, c("strategy", "time", "feature", "scope", "lcl", "ucl")]
  }))
  # strategies discussed in the text (the 0-29 band is too wide to read)
  band <- merge(wc[wc$feature == "race_other" &
                     as.character(wc$strategy) %in% c("0_90", "30_90"), ], bt,
                by = c("strategy", "time", "feature", "scope"))
  band <- band[order(band$scope, band$strategy, band$time), ]
} else if (show_bands) {
  message("Bootstrap files not found (run analysis/03c_wc_sensitivity.R); ",
          "figure without bands.")
}

p_wc <- ggplot(wc, aes(.data$day, .data$rrb, group = .data$line_id)) +
  geom_hline(yintercept = c(-0.1, 0.1), colour = "grey55", linewidth = 0.3,
             linetype = "longdash") +
  (if (!is.null(band)) {
    geom_ribbon(data = band, aes(x = .data$day, ymin = .data$lcl,
                                 ymax = .data$ucl, group = .data$line_id),
                inherit.aes = FALSE, fill = "grey", alpha = 0.2)
  }) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(data = wc[wc$group == "Other baseline covariates", ],
            aes(linetype = .data$strategy_lab), colour = "grey65",
            linewidth = 0.35) +
  geom_line(data = wc[wc$group != "Other baseline covariates", ],
            aes(linetype = .data$strategy_lab), colour = "black",
            linewidth = 0.7) +
  geom_point(data = wc[wc$group != "Other baseline covariates", ],
             aes(shape = .data$strategy_lab), colour = "black", size = 1.6) +
  facet_wrap(~scope_lab, nrow = 1) +
  # all four strategies (incl. "no_0_90") end their grace period by day 90
  scale_x_continuous(breaks = seq(0, 90, 10),
                     limits = c(0, max(90, max(wc$day)))) +
  scale_shape_manual(values = c(16, 17, 15, 1)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
  labs(x = "Days since discharge (end of 10-day interval)",
       y = "Weighted RRB (standardised mean difference)",
       linetype = "Strategy", shape = "Strategy") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0))

ggsave(file.path(fig_dir, "fig4_wc_rrb_profile.png"), p_wc, width = 8,
       height = 4.2, dpi = 300, bg = "white")

cat("Figure written to", file.path(fig_dir, "fig4_wc_rrb_profile.png"),
    "\n")
