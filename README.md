# CRC scRNA-seq Atlas Explorer (R Shiny)

Interactive browser for the per-sample summary tables of the colorectal cancer
single-cell atlas: filter on any metadata field, then plot sample and cell
numbers per study and compare groups.

## Files

| file | purpose |
|---|---|
| `app.R` | the Shiny app (UI + server) |
| `samples.csv` | the input table — one row per sample, both cohorts, 102 columns |
| `.gitignore` | ignores R session files (`.Rhistory`, `.RData`, `.Rproj.user`, …) |
| `crc_atlas_study_overview.csv` | static per-study summary (same content as the app's *Study overview* tab, unfiltered) |
| `crc_atlas_column_overview.csv` | static column dictionary (same content as the app's *Metadata columns* table, unfiltered) |

The first three files are what a `runGitHub` repository needs; the two CSV
summaries are convenience exports and can be left out.

## Run straight from GitHub

Collaborators need no copy of the data or the code — only the packages:

```r
install.packages(c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", "scales", "DT"))
shiny::runGitHub("crc-atlas-explorer", "YOUR-GITHUB-USERNAME")
```

`runGitHub` downloads the repository to a temporary folder and launches the app
from it, so `samples.csv` must sit next to `app.R` in the repo root — it does.
Notes:

- **Branch.** shiny >= 1.5 defaults to `ref = "HEAD"` and resolves whatever the
  default branch is. On older shiny the default is `ref = "master"`, which fails
  on a `main`-branch repo — those collaborators pass `ref = "main"` explicitly.
- **Private repositories are not supported.** `runGitHub` fetches the archive
  anonymously, so a private repo returns a download error. Either make the repo
  public, or have collaborators clone it and run locally:
  `git clone <url> && R -e 'shiny::runApp("crc-atlas-explorer")'`.
- **Nothing is cached.** Each call re-downloads ~1.2 MB; fine interactively, but
  for repeated use cloning once is kinder.
- The app refuses to start with a named list of missing packages rather than a
  stack trace, so a collaborator who skipped the `install.packages` line gets a
  readable message.

## Run from a local folder

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
- *Treatment*: `treatment_status_before_resection`, `treatment_drug`,
  `treatment_response` and `RECIST`, each labelled with its level count and how
  many samples are unannotated (69 / 1,118 / 1,212 / 1,285 respectively).
- *More metadata filters*: pick any of the remaining 66 categorical metadata columns
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

**Metadata columns tab** — the in-app data dictionary. One row per column
(105 in total, including the three derived ones) giving its thematic group,
type, number of levels, how many samples are unannotated and a preview of the
most frequent values; searchable and sortable per field, and downloadable as a
CSV. Clicking a row inspects that column: categorical fields get a frequency
table (samples, % of samples, patients, studies and cells per level, with
`(not annotated)` included) and a bar chart of the 40 most frequent levels;
numeric fields get n / missing / min / Q1 / median / Q3 / max / mean / SD and a
histogram. A *Describe* switch chooses between the currently filtered selection
and the whole atlas, so you can see how a filter reshapes a column — e.g. with
the healthy cohort selected, `KRAS_status` reports 0 levels and 186/186
unannotated.

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
- Treatment annotation is internally inconsistent: of the 1,186 samples marked
  `naive` in `treatment_status_before_resection`, 139 carry a `treatment_drug`
  value — 51 of them the literal string `naive`, but 88 an actual regimen. Two
  `treated` samples have no drug recorded. `treatment_status_before_resection`
  is the more reliable of the two for splitting treated from untreated;
  `treatment_drug` also mixes regimen and cycle count in one free-text field
  (`4x FOLFOX`, `FOLFOX + bevacizumab x 6 cycles`), so its 34 levels are not a
  clean drug vocabulary.
- The 344 treated samples come from 10 studies. `RECIST` is annotated in 5
  studies (158 / 89 / 57 / 8 / 2 samples) and `treatment_response` in 7, both
  dominated by `Chen_2024_Cancer_Cell` — so response-stratified comparisons are
  driven by one or two cohorts and should be read per study (use *Split panels
  by* `study_id` on the Compare tab).
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
