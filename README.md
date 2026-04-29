# Meta-Analysis of Hg and PCB Effects on Fish Reproductive Success

This repository contains the data and R scripts used in the meta-analysis
reported in **[Manuscript citation TBD]**.

The effect-size calculations themselves (Hedges' *d* and Fisher's *Z*,
per endpoint and per toxicant, under both *Conservative* and *Inclusive*
inclusion rules) were performed in `data/Fish_Toxi_Meta_Data.xlsx`.
The R scripts in this repo apply the Borenstein correction for shared
control groups, run publication-bias diagnostics, and produce the
manuscript figures.

## Repository structure

```
Hg_PCBs_Meta_Repo/
├── data/        # Master data sheet, analytic inputs, and intermediate CSVs
├── scripts/     # R scripts for pub-bias analysis and figure generation
├── figures/     # Output figures (created when scripts are run)
└── README.md
```

## Data

| File | Description |
|------|-------------|
| `data/Fish_Toxi_Meta_Data.xlsx` | **Master spreadsheet.** Contains the full study list (`Studies_Table`), search terms, results summary (`Results_all`), and per-endpoint effect-size calculations under both Conservative and Inclusive inclusion rules: `Con_Size`, `Inclusive Size`, `Con_Muscle`, `Inclusive_Muscle`, `Con_Hatching`, `Inclusive Hatching`, `Con_Mortality`, `Inclusive Mortality`, `Con_Growth`, `Inclusive_Growth`. |
| `data/ICC List for Meta.xlsx` | Flat per-study yi/vi table (sheet `Sheet1`) used as the analytic input to the R pub-bias pipeline. |
| `data/borenstein_corrected_hatching.csv` | Borenstein-aggregated hatching effect sizes (output of `borenstein_correction.R`). |
| `data/borenstein_corrected_mortality.csv` | Borenstein-aggregated mortality effect sizes. |
| `data/borenstein_corrected_growth.csv` | Borenstein-aggregated growth effect sizes. |
| `data/borenstein_corrected_muscle.csv` | Borenstein-aggregated muscle (Fisher's Z) effect sizes. |
| `data/borenstein_corrected_size.csv` | Borenstein-aggregated size (Fisher's Z) effect sizes. |
| `data/forest_plot_borenstein.csv` | Long-format table consumed by the forest-plot script. |
| `data/forest plot.xlsx` | Original forest-plot input (size and muscle endpoints). |

## Scripts

### Publication-bias pipeline

1. **`scripts/borenstein_correction.R`** — applies the Borenstein et al.
   (2009) correction for shared control groups to the per-study effect
   sizes in `ICC List for Meta.xlsx`, producing the
   `borenstein_corrected_*.csv` files used downstream.
2. **`scripts/publication_bias_analysis.R`** — runs Egger's regression test,
   funnel plots, and trim-and-fill correction on the Borenstein-corrected
   larval-success data and on the Fisher's-Z maternal data.

### Figure generation

3. **`scripts/forest_plots.R`** — generates forest plots for larval
   success and size/muscle endpoints from the corrected and original
   tables.
4. **`scripts/Publication_effect_sizes.R`** — produces the publication-
   quality summary plots from the final pooled effect sizes (values are
   hard-coded inside the script — these are the values reported in the
   manuscript and computed in `Fish_Toxi_Meta_Data.xlsx`).

## Dependencies

R (>= 4.1) with:

```r
install.packages(c("metafor", "ggplot2", "cowplot", "readxl",
                   "dplyr", "MASS", "gridExtra"))
```

## How to run

All paths in the scripts are relative to the repo root, so:

1. Set R's working directory to the repo root. In RStudio:
   *Session > Set Working Directory > To Project Directory* (or open the
   folder as a project). From the console:
   ```r
   setwd("path/to/Hg_PCBs_Meta_Repo")
   ```
2. Source the scripts in this order:
   ```r
   source("scripts/borenstein_correction.R")     # writes data/borenstein_corrected_*.csv
   source("scripts/publication_bias_analysis.R") # writes figures/funnel_*.png + publication_bias_summary.csv
   source("scripts/forest_plots.R")              # writes figures/*_forest_plot.png
   source("scripts/Publication_effect_sizes.R")  # writes figures/larval_success_plot.* and figures/maternal_characteristics_plot.*
   ```

Figures are written to `figures/`. Intermediate corrected-effect-size
tables are written back to `data/` so downstream scripts can pick them up.

## Citation

If you use these data or scripts, please cite:

> Sosa, B. et al. (YEAR). *[Manuscript title]*. *[Journal]*. doi: TBD

## Contact

Brandon Sosa — bsosa012@fiu.edu
