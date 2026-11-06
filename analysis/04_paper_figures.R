# 04_paper_figures.R ----------------------------------------------------------
# Figure 4 of the manuscript (case study): A, weighted interval RRB under the
# tutorial model; B, studentised maximum per strategy and history function
# against its null 95th percentile, tutorial and revised (use) model;
# C, cumulative RRB of the baseline covariates under no initiation, with
# bootstrap intervals at day 90. Reads output/empirical/ (from 03, 03b, 03c);
# writes output/figures/fig4_wc_rrb_profile.png.
# Run: Rscript analysis/04_paper_figures.R (CCW_FIG4_BANDS=1 adds panel A bands)

if (!file.exists("DESCRIPTION")) stop("Run this script from the project root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(ggplot2)
library(patchwork)

fig_dir <- file.path("output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Webster-Clark RRB profile -------------------------------------------------
delta <- 10 # days per interval
emp_dir <- file.path("output", "empirical")
need <- function(f, hint) {
  path <- file.path(emp_dir, f)
  if (!file.exists(path)) stop("Missing ", path, "; run ", hint, ".")
  read.csv(path)
}
wc_all <- need("wc_rrb.csv", "analysis/03_empirical_webster_clark.R")
use_all <- need("wc_use_rrb.csv", "analysis/03c_wc_sensitivity.R")
wc <- wc_all[wc_all$method == "weighted", ]
# each strategy up to its last interval with artificial censoring
cens <- wc[wc$scope == "interval" & wc$n_censored > 0, ]
last_cens <- tapply(cens$time, cens$strategy, max)
wc <- wc[wc$time <= last_cens[as.character(wc$strategy)], ]

strategy_labels <- c(
  "0_30" = "Initiation,\ndays 0-29", "0_90" = "Initiation,\ndays 0-89",
  "no_0_90" = "No initiation\nbefore day 90"
)
wc$day <- (wc$time + 1) * delta # end of interval, where censoring occurs
wc$strategy_lab <- factor(strategy_labels[wc$strategy],
                          levels = strategy_labels)
# point groups of panel A: the dispensing indicator, its depth functions and
# all other functions (small grey)
depth <- c("fill_lag", "fill_30", "fill_cum")
grp_labels <- c(
  fill = "Dispensing indicator",
  fill_lag = "Preceding-window indicator",
  fill_30 = "Dispensings in 30 days",
  fill_cum = "Dispensings since discharge",
  other = "Other functions"
)
grp <- ifelse(wc$feature %in% c("fill", depth), wc$feature, "other")
wc$group <- factor(grp_labels[grp], levels = grp_labels)
hl <- wc$group != grp_labels[["other"]]
cat("\nPanels A-B (weighted RRB, highlighted functions):\n")
print(wc[hl, c("scope", "strategy", "time", "day", "feature", "rrb")],
      digits = 3, row.names = FALSE)
wc$line_id <- interaction(wc$strategy, wc$feature, wc$scope)

# 95% percentile bootstrap bands (race_other, strategy 0_90)
boot_files <- file.path(emp_dir, "wc_rrb_boot_interval.csv")
band <- NULL
show_bands <- identical(Sys.getenv("CCW_FIG4_BANDS", "0"), "1")
if (show_bands && all(file.exists(boot_files))) {
  bt <- do.call(rbind, lapply(boot_files, function(f) {
    read.csv(f)[, c("strategy", "time", "feature", "scope", "lcl", "ucl")]
  }))
  band <- merge(wc[wc$feature == "race_other" &
                     as.character(wc$strategy) == "0_90", ], bt,
                by = c("strategy", "time", "feature", "scope"))
  band <- band[order(band$scope, band$strategy, band$time), ]
} else if (show_bands) {
  message("Bootstrap files not found (run analysis/03c_wc_sensitivity.R); ",
          "figure without bands.")
}

wc_a <- wc[wc$scope == "interval" & wc$n_censored > 0, ]
# panel C: cumulative RRB under no initiation, tutorial censoring model
model_labels <- c(base = "Tutorial model", use = "Revised model")
base_feats <- c("age", "female", "renal", "race_black", "race_other")
base_labels <- c(age = "Age", female = "Female", renal = "Renal",
                 race_black = "Race Black",
                 race_other = "Race other")
wc_c <- wc_all[wc_all$method == "weighted" & wc_all$scope == "cumulative" &
                 wc_all$strategy == "no_0_90" &
                 wc_all$time <= last_cens[["no_0_90"]], ]
wc_c$day <- (wc_c$time + 1) * delta
wc_c$cov_lab <- factor(base_labels[wc_c$feature], levels = base_labels)
# bootstrap intervals at day 90, staggered to the right of day 90
boot_cum_file <- file.path(emp_dir, "wc_rrb_boot_cumulative.csv")
whisk <- NULL
if (file.exists(boot_cum_file)) {
  whisk <- read.csv(boot_cum_file)
  whisk <- whisk[whisk$strategy == "no_0_90" &
                   whisk$scope == "cumulative" &
                   whisk$time == last_cens[["no_0_90"]] &
                   whisk$feature %in% base_feats, ]
  whisk$cov_lab <- factor(base_labels[whisk$feature], levels = base_labels)
  # whiskers at x = 92, 94, ..., 100 in legend order
  whisk$x <- 92 + 2 * (as.integer(whisk$cov_lab) - 1)
  cat("\nPanel C: bootstrap intervals at day 90 (no_0_90, cumulative)\n")
  print(whisk[, c("feature", "rrb_obs", "boot_se", "lcl", "ucl", "B")],
        digits = 3, row.names = FALSE)
} else {
  message("Missing ", boot_cum_file, " (run analysis/03c_wc_sensitivity.R); ",
          "panel C without bootstrap intervals.")
}
hl_a <- wc_a$group != grp_labels[["other"]]
key_w <- grid::unit(1.4, "cm")
shape_fun <- c(16, 17, 15, 18, 16)
size_fun <- c(1.9, 1.9, 1.8, 2.4, 0.9)
# Okabe-Ito colours (colour-blind safe); shapes and line types are kept as
# a redundant encoding for greyscale printing
oi <- c(blue = "#0072B2", vermillion = "#D55E00", green = "#009E73",
        pink = "#CC79A7", orange = "#E69F00")
col_fun <- c(unname(oi[c("vermillion", "blue", "green", "pink")]), "grey60")
band_of <- function(sc) if (is.null(band)) NULL else band[band$scope == sc, ]
x_scale <- scale_x_continuous(breaks = seq(10, 90, 10), limits = c(5, 95),
                              expand = expansion(mult = 0.02))
# legends at the right so that patchwork aligns the panel bodies
panel_theme <- theme_minimal(base_size = 10) +
  theme(legend.position = "right", legend.title.position = "top",
        legend.margin = ggplot2::margin(0, 0, 0, 4),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 8.8, hjust = 0,
                                  margin = margin(b = 2)),
        strip.text = element_text(hjust = 0))

# initiation strategies censor in one interval each: discrete day axis and
# panel widths proportional to the number of cells
wc_a$day_f <- factor(wc_a$day, levels = sort(unique(wc_a$day)))
p_a <- ggplot(wc_a, aes(.data$day_f, .data$rrb)) +
  (if (!is.null(band_of("interval"))) {
    geom_linerange(data = transform(band_of("interval"),
                                    day_f = factor(day, levels(wc_a$day_f))),
                   aes(x = .data$day_f, ymin = .data$lcl, ymax = .data$ucl),
                   inherit.aes = FALSE, colour = "grey70", linewidth = 2)
  }) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  # shape, colour and size share a name so that they form one legend
  geom_point(data = wc_a[!hl_a, ],
             aes(shape = .data$group, colour = .data$group,
                 size = .data$group)) +
  geom_point(data = wc_a[hl_a, ],
             aes(shape = .data$group, colour = .data$group,
                 size = .data$group)) +
  facet_grid(~strategy_lab, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = expansion(add = 1.3)) +
  scale_shape_manual(values = shape_fun, limits = unname(grp_labels),
                     name = "History function") +
  scale_colour_manual(values = col_fun, limits = unname(grp_labels),
                      name = "History function") +
  scale_size_manual(values = size_fun, limits = unname(grp_labels),
                    name = "History function") +
  guides(shape = guide_legend(ncol = 1),
         colour = guide_legend(ncol = 1),
         size = guide_legend(ncol = 1)) +
  labs(x = "Day since discharge (end of the 10-day interval)",
       y = "Interval RRB") +
  panel_theme

stopifnot(nrow(wc_c) > 0, setequal(wc_c$feature, base_feats))
lt_cov <- c("solid", "dashed", "dotted", "dotdash", "longdash")
# open symbols so that the line type stays visible in the merged legend key
shape_cov <- c(21, 24, 22, 23, 25)
col_cov <- unname(oi[c("blue", "vermillion", "green", "pink", "orange")])
p_cum <- ggplot(wc_c, aes(.data$day, .data$rrb, group = .data$cov_lab)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(aes(linetype = .data$cov_lab, colour = .data$cov_lab),
            linewidth = 0.6) +
  geom_point(aes(shape = .data$cov_lab, colour = .data$cov_lab),
             fill = "white", size = 1.6) +
  (if (!is.null(whisk)) {
    list(
      geom_linerange(data = whisk,
                     aes(x = .data$x, ymin = .data$lcl, ymax = .data$ucl,
                         colour = .data$cov_lab),
                     inherit.aes = FALSE, linewidth = 0.5,
                     show.legend = FALSE),
      geom_point(data = whisk,
                 aes(x = .data$x, y = .data$rrb_obs, shape = .data$cov_lab,
                     colour = .data$cov_lab),
                 inherit.aes = FALSE, fill = "white", size = 2)
    )
  }) +
  scale_x_continuous(breaks = seq(10, 90, 10), limits = c(5, 101.5),
                     expand = expansion(mult = 0.02)) +
  scale_linetype_manual(values = lt_cov, limits = unname(base_labels),
                        name = "Baseline covariate") +
  scale_shape_manual(values = shape_cov, limits = unname(base_labels),
                     name = "Baseline covariate") +
  scale_colour_manual(values = col_cov, limits = unname(base_labels),
                      name = "Baseline covariate") +
  guides(linetype = guide_legend(keywidth = key_w, ncol = 1),
         shape = guide_legend(keywidth = key_w, ncol = 1),
         colour = guide_legend(keywidth = key_w, ncol = 1)) +
  labs(x = "Day since discharge (end of the 10-day interval)",
       y = "Cumulative RRB") +
  panel_theme

# Panel B: per-strategy studentised maximum of the weighted interval RRB --
feature_labels <- c(
  age = "Age", female = "Female", renal = "Renal", race_black = "Race Black",
  race_other = "Race other", visit = "Visit", outpt = "Outpatient",
  fill = "Dispensing indicator", fill_lag = "Preceding-window indicator",
  fill_30 = "Dispensings in 30 days",
  fill_cum = "Dispensings since discharge",
  "(all)" = "All eleven"
)
# censoring-model terms (formulas base and use of analysis/wc_data.R)
in_model <- list(base = c("age", "female", "renal"),
                 use = c("age", "female", "renal", "visit", "outpt", "fill"))
mx <- rbind(
  cbind(model = "base",
        need("wc_rrb_null_maxib_stud.csv",
             "analysis/03b_wc_null_reference.R")),
  cbind(model = "use",
        need("wc_use_null_maxib_stud.csv", "analysis/03c_wc_sensitivity.R"))
)
mx <- mx[mx$strategy != "(all)", ] # per-strategy rows only
mx$feature_lab <- factor(feature_labels[mx$feature],
                         levels = rev(feature_labels))
mx$model_lab <- factor(model_labels[mx$model],
                       levels = unname(model_labels))
mx$strategy_lab <- factor(strategy_labels[mx$strategy],
                          levels = strategy_labels)
mx$type <- factor(as.character(mx$model_lab),
                  levels = unname(model_labels))
cat("\nPanel B: observed studentised maximum and null q95 per strategy",
    "and model\n")
print(data.frame(model = mx$model, strategy = mx$strategy,
                 feature = mx$feature, n_cells = mx$n_cells,
                 stud_obs = round(mx$stud_obs, 3),
                 null_q95 = round(mx$null_q95, 3), p = round(mx$p, 3),
                 exceeds = mx$stud_obs > mx$null_q95),
      row.names = FALSE)
# one row per function: tutorial model in the upper half, revised model in
# the lower half; the half-row is shaded when the function is a model term
off <- c(base = 0.2, use = -0.2)
mx$y <- as.integer(mx$feature_lab) + off[mx$model]
term_rows <- do.call(rbind, lapply(names(in_model), function(m) {
  y <- match(feature_labels[in_model[[m]]], levels(mx$feature_lab))
  data.frame(ymin = if (m == "base") y else y - 0.5,
             ymax = if (m == "base") y + 0.5 else y)
}))
bar_col <- c(base = "#BDD7EC", use = "#F5C9B0")

p_c <- ggplot(mx, aes(y = .data$y)) +
  geom_rect(data = term_rows,
            aes(xmin = -Inf, xmax = Inf, ymin = .data$ymin,
                ymax = .data$ymax),
            inherit.aes = FALSE, fill = "grey92") +
  geom_linerange(aes(xmin = 0, xmax = .data$null_q95),
                 colour = bar_col[mx$model], linewidth = 1.6) +
  geom_point(aes(x = .data$stud_obs, shape = .data$type,
                 fill = .data$type, colour = .data$type),
             size = 1.6, stroke = 0.6) +
  facet_grid(~strategy_lab, scales = "free_x") +
  scale_shape_manual(values = c(21, 21), limits = levels(mx$type),
                     name = "Observed max|z|") +
  scale_colour_manual(values = c(oi[["blue"]], oi[["vermillion"]]),
                      limits = levels(mx$type), name = "Observed max|z|") +
  scale_fill_manual(values = c(oi[["blue"]], "white"),
                    limits = levels(mx$type), name = "Observed max|z|") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06)),
                     n.breaks = 4) +
  scale_y_continuous(breaks = seq_along(levels(mx$feature_lab)),
                     labels = levels(mx$feature_lab),
                     expand = expansion(add = c(0.6, 0.6))) +
  labs(x = paste("Studentised maximum max|z| (bars: null 95th percentile",
                 "of each model; shading: model term)"),
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing.x = grid::unit(1, "lines"),
        legend.position = "right", legend.title.position = "top",
        legend.margin = ggplot2::margin(0, 0, 0, 4),
        strip.text = element_text(hjust = 0))

# separate legends; y-axis titles freed from the alignment with panel B
fig4 <- free(p_a, type = "label", side = "l") /
  p_c /
  free(p_cum, type = "label", side = "l") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
fig4 <- fig4 +
  plot_layout(heights = c(1, 1.3, 0.9))
ggsave(file.path(fig_dir, "fig4_wc_rrb_profile.png"), fig4, width = 8,
       height = 9, dpi = 300, bg = "white")

cat("Figure written to", file.path(fig_dir, "fig4_wc_rrb_profile.png"),
    "\n")
