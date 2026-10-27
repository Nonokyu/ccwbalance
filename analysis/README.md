# Reproducing the manuscript

This folder holds the code that produces the tables and figures of

> Nonomiya Y. Risk-set reconstruction balance: zero-targeted diagnostics for
> artificial censoring in clone–censor–weight analyses. In preparation.

It covers the simulation study (Parts 1–3, Sections 5–6), the checks of the
null reference (Supporting Information, Appendix S2) and the case study on
the Webster-Clark Synthetic Medicare cohort (Section 7, Appendix S3). The
folders `analysis/`, `output/` and `data-raw/` are excluded from the package
build.

## Setup

Run every script from the repository root, which is the package root. Each
script loads the development version of the package with
`pkgload::load_all(".")`.

Packages needed in addition to the package's own imports:

| Package | Used by |
|---|---|
| `pkgload` | all scripts |
| `haven` | case study (`03*.R`), `02g`, `02j` (reads the SAS data file) |
| `survival` | `03_empirical_webster_clark.R`, `02j` |
| `patchwork` or `gridExtra` | dashboard plot of `03_empirical_webster_clark.R` |
| `testthat` | `tests/test-sim_scenarios.R` |

`parallel` (base R) runs the simulations with `parallel::mclapply()`; set the
number of workers with `CCW_CORES` (default: all cores but one). Results were
produced with R 4.4.3.

The case study needs the public cohort of the Webster-Clark CCW tutorial
(no registration, about 200 KB):

```sh
git clone https://github.com/mawcpharmdphd/CCW_tutorial data-raw/webster_clark
```

The scripts read `data-raw/webster_clark/Data/baseccwcohort.sas7bdat`. The
folder is git-ignored.

## Scenarios

`sim_scenarios.R` defines the scenario grids. A *run* is one execution of
`02_simulation_study.R` over the grid `ccw_scenarios(<run>)`. Four runs
(`main`, `manycov`, `largerep`, `variants`) supply the manuscript.
`paper_scenarios()` maps the scenario ids of Table 2 (P1-01 to P1-13,
P2-01 to P2-06, P3-01 to P3-04) to the run scenario whose results are
reported; the table is saved as `scenarios_paper.csv`, and
`paper_scenario(id)` returns the definition of one scenario.
`sim_paper_results.R` reads the replicate-level results of given ids. The
grids of the runs are saved as `scenarios_<run>.csv`.

## Saved outputs

`output/sim/` holds the replicate-level results of the four runs
(`sim_results_<run>.csv`, with `sim_performance_<run>.csv` and
`sim_settings_<run>.csv`) and every summary that the manuscript cites.
`output/figures/` holds the simulation figures. With these files, the table
and figure scripts run in seconds without rerunning any simulation.

## Run order

Fast (seconds; read only saved outputs):

```sh
Rscript analysis/02b_simulation_tables_paper.R   # tables and figures, Parts 1-3
Rscript analysis/02c_simulation_tables_manycov.R # calibration, run manycov
Rscript analysis/02d_size_adjusted_detection.R   # size-adjusted detection
Rscript analysis/02e_se_ratio_robust.R           # robust SE ratio
Rscript -e 'pkgload::load_all("."); testthat::test_file("analysis/tests/test-sim_scenarios.R")'
```

Case study (needs `data-raw/webster_clark/`; outputs go to the git-ignored
`output/empirical/`, and Figure 4 to `output/figures/fig4_wc_rrb_profile.png`,
which is therefore not included here):

```sh
Rscript analysis/03_empirical_webster_clark.R  # about 30 s
Rscript analysis/03b_wc_null_reference.R       # null reference
Rscript analysis/03c_wc_sensitivity.R          # bootstrap, race in model
Rscript analysis/03d_wc_depletion.R            # depletion term
Rscript analysis/04_paper_figures.R            # Figure 4 (after 03)
```

Simulations (long; runtimes with 40 cores in parentheses):

```sh
CCW_RUN=main     Rscript analysis/02_simulation_study.R  # (~2 h)
CCW_RUN=manycov  Rscript analysis/02_simulation_study.R  # (~25 min)
CCW_RUN=largerep Rscript analysis/02_simulation_study.R  # (~2.75 h)
CCW_RUN=variants Rscript analysis/02_simulation_study.R  # (~40 min)
Rscript analysis/02f_cumrrb_censored_strata.R  # needs sim_results_main.csv
Rscript analysis/02g_null_coupling_check.R   # coupling of the null reference
Rscript analysis/02h_null_oracle.R           # oracle check
Rscript analysis/02i_null_decomposition.R    # decomposition of the shortfall
Rscript analysis/02j_depletion_incidence.R   # grace-period incidence (+ data)
```

`02g` and `02j` also use the Webster-Clark data (`02g` skips its case-study
part when the data are absent). After a simulation rerun, rerun `02b`,
`02c`, `02d` and `02e`. The seeds are fixed in the scripts; the header of
each script lists its settings (environment variables) and outputs.

## Outputs behind the manuscript

| Manuscript | Output | Script |
|---|---|---|
| Table 2 (scenarios) | `analysis/scenarios_paper.csv` | `02b` |
| Table 3 (Part 1) | `output/sim/sim_paper_scenarios.csv`, `sim_part1_factorial.csv` | `02b` |
| Table 4 (Part 2) | `output/sim/sim_paper_scenarios.csv` | `02b` |
| Figure 2 (Part 1) | `output/figures/sim_part1_factorial.png` | `02b` |
| Figure 3 (Part 2) | `output/figures/sim_part2_calibration.png` | `02b` |
| Figure 4 (case study) | `output/figures/fig4_wc_rrb_profile.png` | `04` (from `03`) |
| Table S1 | `output/sim/null_coupling_summary.csv`, `output/empirical/wc_null_coupling.csv` | `02g` |
| Table S2 | `output/sim/null_oracle_summary.csv` | `02h` |
| Table S3 | `output/sim/null_decomposition_summary.csv` | `02i` |
| Tables S4, S6, S7 | `output/sim/sim_paper_scenarios.csv` | `02b` |
| Table S5 | `output/sim/sim_part1_detection_size_adjusted.csv`, `sim_part1_detection_size_multipliers.csv` | `02d` |
| Table S8 | `output/empirical/wc_race_comparison.csv`, `wc_*_null_maxib.csv`, `wc_*estimates.csv` | `03`, `03b`, `03c` |
| Table S9 | `output/empirical/wc_rrb.csv`, `wc_race_rrb.csv` | `03`, `03c` |

Numbers cited in the text also come from `sim_paper_calibration.csv`
(`02b`), `sim_se_ratio_robust.csv` (`02e`),
`sim_paper_cumrrb_censored_strata_summary.csv` (`02f`),
`sim_depletion_incidence.csv` and `output/empirical/wc_grace_incidence.csv`
(`02j`), and `output/empirical/wc_depletion.csv` (`03d`). Figure 1 is a
schematic and has no script. `02c` writes the calibration plot of all six
scenarios of the run `manycov`
(`output/figures/sim_part2_calibration_manycov.png`); the manuscript reports
two of them (P2-05, P2-06) through `sim_paper_scenarios.csv`.
