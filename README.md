# CRC scRNA-seq Atlas Explorer (R Shiny)

Interactive browser for the per-sample summary tables of the colorectal cancer
single-cell atlas: filter on any metadata field, then plot sample and cell
numbers per study and compare groups.

## Files

| file | purpose |
|---|---|
| `app.R` | the Shiny app (UI + server) |
| `samples.csv` | the input table — one row per sample, both cohorts, 102 columns |
| `crc_atlas_study_overview.csv` | static per-study summary (same content as the app's *Study overview* tab, unfiltered) |

## Run

```r
install.packages(c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", "scales", "DT"))
setwd("path/to/this/folder")
shiny::runApp()
```

The app reads `samples.csv` from its own folder. To point it elsewhere, use
either of:

```r
options(crc_atlas_file = "/path/to/samples.csv")   # before runApp()
Sys.setenv(CRC_ATLAS_FILE = "/path/to/samples.csv")
```

A table that still has `crc_atlas_samples.csv`'s name is picked up as a
fallback, so older copies keep working.

## Input handling on load

The app derives three columns and leaves the other 102 untouched:

- **`cohort`** — `Healthy` where `medical_condition == "healthy"`, otherwise
  `CRC`; `colorectal polyp` therefore counts as CRC. Filter on
  `medical_condition` to separate cancer (1,210 samples) from polyp (203) from
  healthy (186).
- **`n_cells`** — row sum of the 12 cell-type columns, recomputed live whenever
  the cell-type selection changes.
- **`age_numeric`** — `age` parsed to a number, range midpoint for entries like
  `20-25`.

It also repairs two artefacts of stacking per-cohort exports: rows that merely
repeat the column names are dropped (the supplied `samples.csv` has one, at
position 1413), and the numeric columns such a row turns into text are coerced
back to numbers. `True`/`False` text columns become logical.

## What the app does

**Sidebar — filters** (combine with AND across columns, OR within a column)

- *Quick filters*: `cohort`, `medical_condition`, `study_id`, `sample_type`,
  `sample_tissue`, `platform`.
- *More metadata filters*: pick any of the 76 categorical metadata columns
  (every categorical field with 2–120 levels, i.e. all but the free-text and
  per-sample identifier columns)
  (mutation status, MSI/MMR, TNM stage, CMS and TME subtype, histology,
  treatment, enrichment, country, processing lab, reference genome, …) and a
  value selector appears for it. Each column's label shows its number of levels
  and how many samples are unannotated.
- Columns with missing values offer the explicit choice `(not annotated)`, so
  unannotated samples can be kept or excluded deliberately rather than silently.
- *Numeric filters*: patient age (with a switch for samples whose age is
  unknown), minimum cells per sample, and a cell-type checkbox group that
  restricts **every** cell number in the app — including the minimum-cells
  threshold — to the selected cell types.
- Value boxes above the tabs always show what survived the current filters:
  samples, patients, studies, datasets, cells, median cells per sample.

**Overview tab** — one pair of selectors (`group by` / `colour by`) drives three
plots, so any categorical column can go on the x axis:

1. *Samples per group* — sample counts, stacked / dodged / normalised to 100%.
2. *Cells per group* — summed cell numbers, optional log10 axis.
3. *Cell-type composition* — stacked proportions, absolute counts, or a
   proportion heatmap of the 12 cell types.

Each plot has a PNG download, and the exact plotted numbers are shown underneath
in a table with CSV/Excel export.

**Compare groups tab** — choose a grouping variable (e.g. `cohort`,
`microsatellite_status`, `KRAS_status`, `treatment_status_before_resection`) and
optionally split panels by a second variable. Shows samples and cells per group,
the distribution of cells per sample (log10 boxplot + points), cell-type
composition per group, and a summary table of samples / patients / studies /
cells / median cells per sample.

**Study overview tab** — one row per study with datasets, samples, patients,
cells, median cells per sample, sample types, tissues, platforms, countries,
PMID and DOI; searchable per column and downloadable.

**Sample table tab** — the filtered samples with a user-chosen column set, plus
a sidebar download of the full filtered table.

## Data notes

- 1,599 samples, 42 studies, 69 datasets, 602 patients, 4,248,732 cells
  (1,413 CRC-cohort and 186 healthy-cohort samples).
- All 203 polyp samples come from a single study (`Chen_2021_Cell`), so any
  polyp-vs-tumour comparison is also a between-study comparison.
- `n_cells` is the row sum of the 12 cell-type count columns. `low_cnv_burden`
  and `high_cnv_burden` are kept as metadata and are **not** part of that total.
- `age` mixes exact ages and ranges (e.g. `20-25`); `age_numeric` uses the range
  midpoint, and 271 samples have no parsable age.
- Mutation-status, CMS and TME-subtype columns are annotated for a minority of
  samples, so filtering on them shrinks the cohort sharply — read the value
  boxes to see how much.
- `sample_type` mixes tumour, adjacent normal, healthy normal, polyp, lymph node,
  metastasis and blood; a filter on `medical_condition` alone does not separate
  tumour from adjacent-normal tissue.
