# ─────────────────────────────────────────────────────────────────────────────
# CRC scRNA-seq Atlas Explorer
# Browse, filter and compare per-sample metadata and cell counts of a
# colorectal cancer single-cell atlas.
#
# Run with:  shiny::runApp()      (from this directory)
# Input:     samples.csv  next to this file -- one row per sample, both cohorts.
#            Override with  options(crc_atlas_file = "/path/to/samples.csv")  or
#            the CRC_ATLAS_FILE environment variable.
# ─────────────────────────────────────────────────────────────────────────────

.required <- c("shiny", "bslib", "dplyr", "tidyr", "ggplot2", "scales", "DT")
.missing <- .required[!vapply(.required, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("CRC Atlas Explorer needs these R packages: ", paste(.missing, collapse = ", "),
       "\nInstall them with:\n  install.packages(c(",
       paste0('"', .missing, '"', collapse = ", "), "))", call. = FALSE)

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(DT)

# ── Data ─────────────────────────────────────────────────────────────────────

CELLTYPES <- c("Myeloid cell", "Granulocyte", "B cell", "Plasma cell", "T cell",
               "NK cell", "ILC", "Cancer cell", "Epithelial cell", "Hepatocyte",
               "Stromal cell", "Schwann cell")

MISSING_LABEL <- "(not annotated)"

parse_age <- function(x) {
  s <- trimws(as.character(x))
  out <- rep(NA_real_, length(s))
  single <- grepl("^[0-9]+(\\.[0-9]+)?$", s)
  out[single] <- as.numeric(s[single])
  rng <- grepl("^[0-9]+\\s*-\\s*[0-9]+$", s)
  if (any(rng)) {
    parts <- strsplit(s[rng], "\\s*-\\s*")
    out[rng] <- vapply(parts, function(p) mean(as.numeric(p)), numeric(1))
  }
  out
}

DATA_CANDIDATES <- c("samples.csv", "crc_atlas_samples.csv")

find_data_file <- function() {
  f <- getOption("crc_atlas_file", Sys.getenv("CRC_ATLAS_FILE", ""))
  if (nzchar(f)) {
    if (!file.exists(f)) stop("Atlas file not found: ", f)
    return(f)
  }
  hit <- DATA_CANDIDATES[file.exists(DATA_CANDIDATES)]
  if (!length(hit))
    stop("No atlas table found. Put samples.csv next to app.R, or set ",
         "options(crc_atlas_file = \"/path/to/samples.csv\").")
  hit[1]
}

# Tables assembled by stacking per-cohort exports can carry the second file's
# header as a data row; drop any row that is mostly a copy of the column names.
drop_header_rows <- function(a) {
  is_hdr <- apply(a, 1, function(r) mean(as.character(r) == names(a), na.rm = TRUE) > 0.5)
  is_hdr[is.na(is_hdr)] <- FALSE
  if (any(is_hdr)) message("Dropped ", sum(is_hdr), " repeated header row(s) from the input table.")
  a[!is_hdr, , drop = FALSE]
}

load_atlas <- function() {
  path <- find_data_file()
  a <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("NA", "", "NaN", "nan"))
  a <- drop_header_rows(a)

  missing_ct <- setdiff(CELLTYPES, names(a))
  if (length(missing_ct)) stop("Missing cell-type columns: ", paste(missing_ct, collapse = ", "))

  # a stray header row forces whole columns to character on read -> re-type here
  num_cols <- intersect(c(CELLTYPES, "low_cnv_burden", "high_cnv_burden"), names(a))
  a[num_cols] <- lapply(a[num_cols], function(z) suppressWarnings(as.numeric(as.character(z))))
  a[CELLTYPES] <- lapply(a[CELLTYPES], function(z) { z[is.na(z)] <- 0; z })
  a[] <- lapply(a, function(z) {
    if (is.character(z) && all(na.omit(z) %in% c("True", "False"))) z == "True" else z
  })

  a$age_numeric <- parse_age(a$age)
  a$n_cells <- rowSums(a[CELLTYPES])

  # CRC vs healthy split, kept as its own column so it can be filtered and
  # compared directly ("colorectal polyp" belongs to the CRC cohort).
  if (!"cohort" %in% names(a))
    a$cohort <- ifelse(is.na(a$medical_condition), NA_character_,
                       ifelse(a$medical_condition == "healthy", "Healthy", "CRC"))

  attr(a, "source_file") <- path
  a
}

atlas <- load_atlas()
ATLAS_SOURCE <- attr(atlas, "source_file")

# Categorical columns usable as filters / grouping variables (drop ID-like and
# free-text fields with too many levels to be useful on an axis).
.cat_cols <- names(atlas)[vapply(atlas, function(z)
  is.character(z) || is.logical(z) || is.factor(z), logical(1))]
.nlev <- vapply(atlas[.cat_cols], function(z) length(unique(na.omit(as.character(z)))), integer(1))
.drop <- c("sample_id", "original_obs_names", "study_doi")
GROUP_COLS  <- setdiff(.cat_cols[.nlev >= 2 & .nlev <= 120], .drop)
FILTER_COLS <- GROUP_COLS

# Thematic grouping, used to navigate the metadata-column overview.
COL_GROUPS <- list(
  identity   = c("dataset", "study_id", "sample_id", "patient_id", "donor_id", "replicate",
                 "original_obs_names"),
  condition  = c("cohort", "medical_condition", "cancer_type", "sample_type", "tumor_source",
                 "sample_tissue", "anatomic_region", "anatomic_location", "disease",
                 "tissue", "tissue_type", "is_primary_data"),
  tumor      = c("tumor_stage", "tumor_stage_TNM", "tumor_stage_TNM_T", "tumor_stage_TNM_N",
                 "tumor_stage_TNM_M", "tumor_size", "tumor_dimensions", "tumor_grade",
                 "histological_type", "tumor_immune_type", "CMS_type_Guinney_2015_Nat_Med",
                 "TME_subtype_Bagaev_Cancer_Cell_2021"),
  molecular  = c("microsatellite_status", "mismatch_repair_deficiency_status",
                 "MLH1_promoter_methylation_status"),
  patient    = c("sex", "age", "age_numeric", "self_reported_ethnicity", "development_stage",
                 "organism"),
  treatment  = c("treatment_status_before_resection", "treatment_drug", "treatment_response",
                 "RECIST"),
  technical  = c("platform", "platform_fine", "assay", "cellranger_version", "reference_genome",
                 "matrix_type", "enrichment_cell_types", "tissue_cell_state", "suspension_type"),
  provenance = c("tissue_processing_lab", "hospital_location", "country",
                 "NCBI_BioProject_accession", "SRA_sample_accession", "GEO_sample_accession",
                 "ENA_sample_accession", "synapse_sample_accession", "study_doi", "study_pmid"),
  cell_counts = c(CELLTYPES, "low_cnv_burden", "high_cnv_burden", "n_cells")
)
# gene-level status columns, and the CELLxGENE ontology ids, are regular patterns
COL_GROUPS$molecular <- c(COL_GROUPS$molecular,
  setdiff(grep("_status$", names(atlas), value = TRUE),
          c("microsatellite_status", "mismatch_repair_deficiency_status",
            "MLH1_promoter_methylation_status", "treatment_status_before_resection")))
COL_GROUPS$ontology <- grep("_ontology_term_id$", names(atlas), value = TRUE)

col_group <- function(cols) {
  lookup <- unlist(lapply(names(COL_GROUPS), function(g)
    setNames(rep(g, length(COL_GROUPS[[g]])), COL_GROUPS[[g]])))
  out <- unname(lookup[cols])
  out[is.na(out)] <- "other"
  out
}

fmt_num <- function(x) {
  if (!is.finite(x)) return("-")
  if (abs(x - round(x)) < 1e-9) format(round(x), big.mark = ",", trim = TRUE)
  else format(round(x, 2), big.mark = ",", trim = TRUE)
}

# Columns offered as always-visible "quick" filters.
QUICK <- intersect(c("cohort", "medical_condition", "study_id", "sample_type",
                     "sample_tissue", "platform"), FILTER_COLS)
# Treatment fields get their own always-visible panel.
TREATMENT <- intersect(c("treatment_status_before_resection", "treatment_drug",
                         "treatment_response", "RECIST"), FILTER_COLS)
PINNED <- c(QUICK, TREATMENT)
EXTRA_DEFAULT <- character(0)

choices_for <- function(col) {
  v <- atlas[[col]]
  ch <- sort(unique(na.omit(as.character(v))))
  if (anyNA(v)) ch <- c(ch, MISSING_LABEL)
  ch
}

n_missing <- function(col) sum(is.na(atlas[[col]]))

# Label a column for the UI: name + number of levels + missingness.
col_label <- function(col) {
  sprintf("%s  (%d levels, %d NA)", col, length(setdiff(choices_for(col), MISSING_LABEL)),
          n_missing(col))
}

match_vals <- function(x, vals) {
  keep <- as.character(x) %in% setdiff(vals, MISSING_LABEL)
  if (MISSING_LABEL %in% vals) keep <- keep | is.na(x)
  keep
}

lab_na <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- MISSING_LABEL
  x
}

PALETTE <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860",
             "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD", "#1F4E79", "#A23B72",
             "#2E8B57", "#B8860B", "#6A5ACD", "#CD5C5C", "#20B2AA", "#708090")

fill_scale <- function(n, name) {
  if (n <= length(PALETTE)) {
    scale_fill_manual(values = setNames(PALETTE[seq_len(n)], NULL), name = name,
                      na.value = "grey80", drop = FALSE)
  } else {
    scale_fill_manual(values = hue_pal()(n), name = name, na.value = "grey80", drop = FALSE)
  }
}

base_theme <- theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title.position = "plot",
        legend.position = "right",
        strip.text = element_text(face = "bold"))

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# ── UI ───────────────────────────────────────────────────────────────────────

sidebar_ui <- sidebar(
  width = 340, open = TRUE,
  accordion(
    open = c("Quick filters", "Treatment", "Numeric filters"),
    accordion_panel(
      "Quick filters", icon = NULL,
      lapply(QUICK, function(col)
        selectizeInput(paste0("flt_", make.names(col)), col,
                       choices = choices_for(col), selected = NULL, multiple = TRUE,
                       options = list(placeholder = "all", plugins = list("remove_button")))
      )
    ),
    accordion_panel(
      "Treatment", icon = NULL,
      helpText("Only a subset of samples carries treatment annotation - the",
               "counts in each label show how many are unannotated."),
      lapply(TREATMENT, function(col)
        selectizeInput(paste0("flt_", make.names(col)), col_label(col),
                       choices = choices_for(col), selected = NULL, multiple = TRUE,
                       options = list(placeholder = "all", plugins = list("remove_button")))
      )
    ),
    accordion_panel(
      "More metadata filters",
      helpText("Pick any metadata column, then choose the values to keep."),
      selectizeInput("extra_cols", NULL, choices = setdiff(FILTER_COLS, PINNED),
                     selected = EXTRA_DEFAULT, multiple = TRUE,
                     options = list(placeholder = "add a column to filter on...",
                                    plugins = list("remove_button"))),
      uiOutput("extra_filters")
    ),
    accordion_panel(
      "Numeric filters",
      sliderInput("age", "Patient age (years)",
                  min = floor(min(atlas$age_numeric, na.rm = TRUE)),
                  max = ceiling(max(atlas$age_numeric, na.rm = TRUE)),
                  value = c(floor(min(atlas$age_numeric, na.rm = TRUE)),
                            ceiling(max(atlas$age_numeric, na.rm = TRUE))), step = 1),
      checkboxInput("age_na", "keep samples with unknown age", TRUE),
      sliderInput("mincells", "Minimum cells per sample",
                  min = 0, max = signif(max(atlas$n_cells), 2),
                  value = 0, step = 100),
      checkboxGroupInput("celltypes_keep", "Count only these cell types",
                         choices = CELLTYPES, selected = CELLTYPES, inline = FALSE)
    )
  ),
  hr(),
  actionButton("reset", "Reset all filters", class = "btn-outline-secondary btn-sm"),
  downloadButton("dl_filtered", "Download filtered table (CSV)", class = "btn-sm")
)

plot_tab <- function(id, title, extra = NULL) {
  nav_panel(
    title,
    card(
      full_screen = TRUE,
      card_header(
        div(class = "d-flex justify-content-between align-items-center",
            span(title),
            downloadButton(paste0("dl_", id, "_png"), "PNG", class = "btn-sm btn-light"))
      ),
      card_body(extra, uiOutput(paste0(id, "_ui")))
    ),
    card(card_header("Plotted values"),
         card_body(DTOutput(paste0(id, "_tbl"))))
  )
}

ui <- page_sidebar(
  title = "CRC scRNA-seq Atlas Explorer",
  theme = bs_theme(version = 5, preset = "flatly"),
  sidebar = sidebar_ui,
  layout_columns(
    fill = FALSE, col_widths = c(2, 2, 2, 2, 2, 2),
    value_box("Samples", textOutput("vb_samples"), theme = "primary"),
    value_box("Patients", textOutput("vb_patients"), theme = "secondary"),
    value_box("Studies", textOutput("vb_studies"), theme = "secondary"),
    value_box("Datasets", textOutput("vb_datasets"), theme = "secondary"),
    value_box("Cells", textOutput("vb_cells"), theme = "success"),
    value_box("Median cells/sample", textOutput("vb_median"), theme = "light")
  ),
  navset_card_tab(
    id = "tabs",
    nav_panel(
      "Overview",
      layout_columns(
        col_widths = c(12),
        card(card_header("Group by / colour by"),
             card_body(
               layout_columns(
                 col_widths = c(4, 4, 4),
                 selectizeInput("xvar", "Group samples by (x axis)", choices = GROUP_COLS,
                                selected = "study_id"),
                 selectizeInput("colvar", "Colour by", choices = c("(none)", GROUP_COLS),
                                selected = "sample_type"),
                 radioButtons("barmode", "Bars", inline = TRUE,
                              choices = c("stacked", "dodged", "proportion"),
                              selected = "stacked")
               ),
               layout_columns(
                 col_widths = c(4, 4, 4),
                 radioButtons("ordering", "Order x by", inline = TRUE,
                              choices = c("value", "name"), selected = "value"),
                 checkboxInput("logx", "log10 axis (cell counts)", FALSE),
                 checkboxInput("flip", "horizontal bars", TRUE)
               )
             ))
      ),
      navset_card_underline(
        plot_tab("samples", "Samples per group"),
        plot_tab("cells", "Cells per group"),
        plot_tab("comp", "Cell-type composition",
                 extra = radioButtons("compmode", NULL, inline = TRUE,
                                      choices = c("proportion", "absolute", "heatmap"),
                                      selected = "proportion"))
      )
    ),
    nav_panel(
      "Compare groups",
      card(card_body(layout_columns(
        col_widths = c(5, 5, 2),
        selectizeInput("cmpvar", "Compare by", choices = GROUP_COLS, selected = "cohort"),
        selectizeInput("cmpfacet", "Split panels by (optional)",
                       choices = c("(none)", GROUP_COLS), selected = "(none)"),
        checkboxInput("cmp_drop_na", "drop not-annotated", TRUE)
      ))),
      layout_columns(
        col_widths = c(6, 6),
        card(full_screen = TRUE, card_header("Samples and cells per group"),
             card_body(uiOutput("cmp_bar_ui"))),
        card(full_screen = TRUE, card_header("Cells per sample"),
             card_body(uiOutput("cmp_box_ui")))
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(full_screen = TRUE, card_header("Cell-type composition per group"),
             card_body(uiOutput("cmp_comp_ui"))),
        card(full_screen = TRUE, card_header("Group summary"),
             card_body(DTOutput("cmp_tbl")))
      )
    ),
    nav_panel(
      "Metadata columns",
      card(card_body(layout_columns(
        col_widths = c(5, 4, 3),
        selectizeInput("ov_col", "Inspect column", choices = NULL, selected = NULL),
        radioButtons("ov_scope", "Describe", inline = TRUE,
                     choices = c("current filters", "whole atlas"),
                     selected = "current filters"),
        div(class = "pt-4", downloadButton("dl_overview", "Column overview (CSV)",
                                           class = "btn-sm"))
      ))),
      layout_columns(
        col_widths = c(7, 5),
        card(full_screen = TRUE,
             card_header("All metadata columns - click a row to inspect it"),
             card_body(DTOutput("ov_tbl"))),
        card(full_screen = TRUE, card_header(textOutput("ov_detail_header")),
             card_body(uiOutput("ov_plot_ui"), DTOutput("ov_detail_tbl")))
      )
    ),
    nav_panel(
      "Study overview",
      card(full_screen = TRUE,
           card_header(div(class = "d-flex justify-content-between align-items-center",
                           span("One row per study"),
                           downloadButton("dl_studies", "CSV", class = "btn-sm btn-light"))),
           card_body(DTOutput("study_tbl")))
    ),
    nav_panel(
      "Sample table",
      card(card_body(selectizeInput(
        "tbl_cols", "Columns to show", choices = names(atlas), multiple = TRUE,
        selected = intersect(c("cohort", "study_id", "dataset", "sample_id", "patient_id",
                               "sample_type", "sample_tissue", "anatomic_region",
                               "tumor_stage", "microsatellite_status", "sex", "age",
                               "platform", "n_cells"), names(atlas)),
        options = list(plugins = list("remove_button"))))),
      card(full_screen = TRUE, card_header("Filtered samples"),
           card_body(DTOutput("sample_tbl")))
    ),
    nav_panel(
      "About",
      card(card_body(
        markdown(sprintf(
          "### CRC scRNA-seq Atlas Explorer

Loaded `%s`: **%d samples**, **%d studies**, **%d datasets**, **%d patients**,
**%s cells** across %d annotated cell types.

**`cohort`** is derived from `medical_condition`: `Healthy` for healthy samples,
`CRC` for both `colorectal cancer` and `colorectal polyp`. Use
`medical_condition` itself to separate cancer from polyp. Rows that merely
repeat the column names (an artefact of stacking per-cohort exports) are dropped
on load.

**Filtering.** Quick filters cover the most used fields; *More metadata filters*
lets you filter on any of the %d categorical metadata columns. Each filter lists
the value `%s` when the column has missing entries, so you can explicitly keep
or exclude unannotated samples. Filters combine with AND across columns and OR
within a column. *Count only these cell types* restricts every cell number in the
app (including the minimum-cells filter) to the selected cell types.

**Cell numbers** are per-sample counts of cells passing the atlas QC, summed from
the %d cell-type columns. `low_cnv_burden` / `high_cnv_burden` are kept in the
metadata but are not part of the cell-type totals.

**Caveats.** `age` mixes exact ages and ranges; the numeric age filter uses the
range midpoint and samples without a parsable age are controlled by the
*keep samples with unknown age* checkbox. Mutation-status and subtype columns are
annotated for a minority of samples, so filtering on them shrinks the cohort
sharply — the value boxes always show what survived.",
          basename(ATLAS_SOURCE),
          nrow(atlas), dplyr::n_distinct(atlas$study_id), dplyr::n_distinct(atlas$dataset),
          dplyr::n_distinct(atlas$patient_id), fmt_n(sum(atlas$n_cells)), length(CELLTYPES),
          length(FILTER_COLS), MISSING_LABEL, length(CELLTYPES)))
      ))
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # dynamic UI for user-chosen filter columns
  output$extra_filters <- renderUI({
    cols <- input$extra_cols
    if (is.null(cols) || !length(cols)) return(NULL)
    lapply(cols, function(col) {
      id <- paste0("flt_", make.names(col))
      selectizeInput(id, col_label(col), choices = choices_for(col),
                     selected = isolate(input[[id]]), multiple = TRUE,
                     options = list(placeholder = "all", plugins = list("remove_button")))
    })
  })

  observeEvent(input$reset, {
    for (col in union(PINNED, input$extra_cols))
      updateSelectizeInput(session, paste0("flt_", make.names(col)), selected = character(0))
    updateSelectizeInput(session, "extra_cols", selected = character(0))
    updateSliderInput(session, "age", value = c(floor(min(atlas$age_numeric, na.rm = TRUE)),
                                                ceiling(max(atlas$age_numeric, na.rm = TRUE))))
    updateCheckboxInput(session, "age_na", value = TRUE)
    updateSliderInput(session, "mincells", value = 0)
    updateCheckboxGroupInput(session, "celltypes_keep", selected = CELLTYPES)
  })

  keep_ct <- reactive({
    ct <- input$celltypes_keep
    if (is.null(ct) || !length(ct)) CELLTYPES else ct
  })

  # filtered data; n_cells recomputed over the selected cell types
  dat <- reactive({
    d <- atlas
    d$n_cells <- rowSums(d[keep_ct()], na.rm = TRUE)

    for (col in union(PINNED, input$extra_cols)) {
      v <- input[[paste0("flt_", make.names(col))]]
      if (!is.null(v) && length(v)) d <- d[match_vals(d[[col]], v), , drop = FALSE]
    }

    rng <- input$age
    if (!is.null(rng)) {
      keep <- !is.na(d$age_numeric) & d$age_numeric >= rng[1] & d$age_numeric <= rng[2]
      if (isTRUE(input$age_na)) keep <- keep | is.na(d$age_numeric)
      d <- d[keep, , drop = FALSE]
    }
    if (!is.null(input$mincells)) d <- d[d$n_cells >= input$mincells, , drop = FALSE]
    d
  })

  output$vb_samples  <- renderText(fmt_n(nrow(dat())))
  output$vb_patients <- renderText(fmt_n(n_distinct(dat()$patient_id)))
  output$vb_studies  <- renderText(fmt_n(n_distinct(dat()$study_id)))
  output$vb_datasets <- renderText(fmt_n(n_distinct(dat()$dataset)))
  output$vb_cells    <- renderText(fmt_n(sum(dat()$n_cells)))
  output$vb_median   <- renderText({
    d <- dat(); if (!nrow(d)) "-" else fmt_n(round(median(d$n_cells)))
  })

  # ---- aggregation shared by the two overview bar charts --------------------
  agg <- function(measure) {
    d <- dat()
    validate(need(nrow(d) > 0, "No samples match the current filters."))
    xv <- input$xvar; cv <- input$colvar
    d$.x <- lab_na(d[[xv]])
    d$.grp <- if (is.null(cv) || cv == "(none)") "all samples" else lab_na(d[[cv]])
    out <- d %>% group_by(.x, .grp) %>%
      summarise(value = if (measure == "samples") dplyr::n() else sum(n_cells),
                .groups = "drop")
    tot <- out %>% group_by(.x) %>% summarise(t = sum(value), .groups = "drop")
    lev <- if (identical(input$ordering, "value")) tot$.x[order(tot$t)] else sort(unique(out$.x), decreasing = TRUE)
    out$.x <- factor(out$.x, levels = lev)
    out$.grp <- factor(out$.grp, levels = sort(unique(out$.grp)))
    out
  }

  bar_plot <- function(measure) {
    out <- agg(measure)
    pos <- switch(input$barmode, stacked = "stack", dodged = position_dodge(preserve = "single"),
                  proportion = "fill")
    ylab <- if (measure == "samples") "Samples" else
      paste0("Cells (", length(keep_ct()), "/", length(CELLTYPES), " cell types)")
    if (identical(input$barmode, "proportion")) ylab <- paste0("Fraction of ", tolower(ylab))
    p <- ggplot(out, aes(x = .x, y = value, fill = .grp)) +
      geom_col(position = pos, width = 0.78) +
      fill_scale(nlevels(out$.grp), if (is.null(input$colvar) || input$colvar == "(none)") NULL else input$colvar) +
      labs(x = input$xvar, y = ylab) + base_theme
    if (identical(input$barmode, "proportion")) {
      p <- p + scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0, 0.02)))
    } else if (isTRUE(input$logx) && measure == "cells") {
      p <- p + scale_y_log10(labels = label_number(big.mark = ","), expand = expansion(c(0, 0.05)))
    } else {
      p <- p + scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(c(0, 0.05)))
    }
    if (isTRUE(input$flip)) p + coord_flip()
    else p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  plot_h <- function(n_lev) {
    if (isTRUE(input$flip)) max(320, min(1600, 60 + 20 * n_lev)) else 460
  }

  register_plot <- function(id, plot_fun, data_fun, nlev_fun) {
    output[[paste0(id, "_ui")]] <- renderUI({
      plotOutput(paste0(id, "_plot"), height = paste0(plot_h(nlev_fun()), "px"))
    })
    output[[paste0(id, "_plot")]] <- renderPlot(plot_fun(), res = 100)
    output[[paste0(id, "_tbl")]] <- renderDT(
      datatable(data_fun(), rownames = FALSE, extensions = "Buttons",
                options = list(pageLength = 10, dom = "Bfrtip", buttons = c("csv", "excel"),
                               scrollX = TRUE))
    )
    output[[paste0("dl_", id, "_png")]] <- downloadHandler(
      filename = function() paste0("crc_atlas_", id, "_", Sys.Date(), ".png"),
      content = function(file)
        ggsave(file, plot_fun(), width = 10,
               height = max(4, min(18, plot_h(nlev_fun()) / 80)), dpi = 300, bg = "white")
    )
  }

  pretty_agg <- function(measure) {
    out <- agg(measure)
    names(out) <- c(input$xvar, if (is.null(input$colvar) || input$colvar == "(none)") "group" else input$colvar,
                    if (measure == "samples") "n_samples" else "n_cells")
    out
  }

  register_plot("samples", function() bar_plot("samples"),
                function() pretty_agg("samples"),
                function() n_distinct(lab_na(dat()[[input$xvar]])))
  register_plot("cells", function() bar_plot("cells"),
                function() pretty_agg("cells"),
                function() n_distinct(lab_na(dat()[[input$xvar]])))

  # ---- cell-type composition ------------------------------------------------
  comp_data <- reactive({
    d <- dat()
    validate(need(nrow(d) > 0, "No samples match the current filters."))
    xv <- input$xvar
    d$.x <- lab_na(d[[xv]])
    out <- d %>% select(all_of(c(".x", keep_ct()))) %>%
      pivot_longer(all_of(keep_ct()), names_to = "cell_type", values_to = "n_cells") %>%
      group_by(.x, cell_type) %>% summarise(n_cells = sum(n_cells), .groups = "drop") %>%
      group_by(.x) %>%
      mutate(fraction = if (sum(n_cells) > 0) n_cells / sum(n_cells) else NA_real_) %>%
      ungroup()
    tot <- out %>% group_by(.x) %>% summarise(t = sum(n_cells), .groups = "drop")
    lev <- if (identical(input$ordering, "value")) tot$.x[order(tot$t)] else sort(unique(out$.x), decreasing = TRUE)
    out$.x <- factor(out$.x, levels = lev)
    out$cell_type <- factor(out$cell_type, levels = intersect(CELLTYPES, keep_ct()))
    out
  })

  comp_plot <- function() {
    out <- comp_data()
    if (identical(input$compmode, "heatmap")) {
      p <- ggplot(out, aes(x = cell_type, y = .x, fill = fraction)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        scale_fill_viridis_c(labels = percent_format(accuracy = 1), name = "fraction\nof cells",
                             na.value = "grey90", option = "mako", direction = -1) +
        labs(x = NULL, y = input$xvar) + base_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank())
      return(p)
    }
    yv <- if (identical(input$compmode, "proportion")) "fraction" else "n_cells"
    p <- ggplot(out, aes(x = .x, y = .data[[yv]], fill = cell_type)) +
      geom_col(width = 0.78) +
      fill_scale(nlevels(out$cell_type), "cell type") +
      labs(x = input$xvar,
           y = if (yv == "fraction") "Fraction of cells" else "Cells") +
      base_theme
    p <- if (yv == "fraction")
      p + scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0, 0.02)))
    else p + scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(c(0, 0.05)))
    if (isTRUE(input$flip)) p + coord_flip()
    else p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  output$comp_ui <- renderUI({
    n <- n_distinct(lab_na(dat()[[input$xvar]]))
    h <- if (identical(input$compmode, "heatmap")) max(320, min(1600, 120 + 18 * n)) else plot_h(n)
    plotOutput("comp_plot", height = paste0(h, "px"))
  })
  output$comp_plot <- renderPlot(comp_plot(), res = 100)
  output$comp_tbl <- renderDT({
    out <- comp_data()
    names(out)[1] <- input$xvar
    out$fraction <- round(out$fraction, 4)
    datatable(out, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 10, dom = "Bfrtip", buttons = c("csv", "excel"),
                             scrollX = TRUE))
  })
  output$dl_comp_png <- downloadHandler(
    filename = function() paste0("crc_atlas_composition_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, comp_plot(), width = 10, height = 8, dpi = 300, bg = "white")
  )

  # ---- compare groups -------------------------------------------------------
  cmp_data <- reactive({
    d <- dat()
    validate(need(nrow(d) > 0, "No samples match the current filters."))
    d$.g <- lab_na(d[[input$cmpvar]])
    d$.f <- if (identical(input$cmpfacet, "(none)")) "all" else lab_na(d[[input$cmpfacet]])
    if (isTRUE(input$cmp_drop_na)) d <- d[d$.g != MISSING_LABEL, , drop = FALSE]
    validate(need(nrow(d) > 0, "No annotated samples for this grouping variable."))
    d
  })

  cmp_summary <- reactive({
    cmp_data() %>% group_by(.g) %>%
      summarise(n_samples = dplyr::n(), n_patients = n_distinct(patient_id),
                n_studies = n_distinct(study_id),
                median_cells_per_sample = median(n_cells),
                n_cells = sum(n_cells), .groups = "drop") %>%
      relocate(n_cells, .before = median_cells_per_sample) %>%
      arrange(desc(n_samples))
  })

  cmp_bar_plot <- function() {
    s <- cmp_summary() %>%
      select(.g, n_samples, n_cells) %>%
      pivot_longer(c(n_samples, n_cells), names_to = "measure", values_to = "value") %>%
      mutate(measure = factor(measure, c("n_samples", "n_cells"),
                              c("samples", "cells")),
             .g = factor(.g, levels = rev(cmp_summary()$.g)))
    ggplot(s, aes(x = .g, y = value, fill = measure)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_text(aes(label = fmt_n(value)), hjust = -0.1, size = 3.2) +
      facet_wrap(~measure, scales = "free_x") +
      scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(c(0, 0.25))) +
      fill_scale(2, NULL) + coord_flip() +
      labs(x = input$cmpvar, y = NULL) + base_theme
  }

  cmp_box_plot <- function() {
    d <- cmp_data()
    d$.g <- factor(d$.g, levels = rev(cmp_summary()$.g))
    p <- ggplot(d, aes(x = .g, y = pmax(n_cells, 1), fill = .g)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, show.legend = FALSE) +
      geom_jitter(height = 0, width = 0.18, size = 0.6, alpha = 0.35, show.legend = FALSE) +
      guides(fill = "none") +
      scale_y_log10(labels = label_number(big.mark = ",")) +
      fill_scale(n_distinct(d$.g), NULL) +
      coord_flip() + labs(x = input$cmpvar, y = "Cells per sample (log10)") + base_theme
    if (!identical(input$cmpfacet, "(none)")) p <- p + facet_wrap(~.f)
    p
  }

  cmp_comp_plot <- function() {
    d <- cmp_data()
    out <- d %>% select(all_of(c(".g", ".f", keep_ct()))) %>%
      pivot_longer(all_of(keep_ct()), names_to = "cell_type", values_to = "n_cells") %>%
      group_by(.g, .f, cell_type) %>% summarise(n_cells = sum(n_cells), .groups = "drop") %>%
      group_by(.g, .f) %>%
      mutate(fraction = if (sum(n_cells) > 0) n_cells / sum(n_cells) else NA_real_) %>%
      ungroup() %>%
      mutate(cell_type = factor(cell_type, levels = intersect(CELLTYPES, keep_ct())),
             .g = factor(.g, levels = rev(cmp_summary()$.g)))
    p <- ggplot(out, aes(x = .g, y = fraction, fill = cell_type)) +
      geom_col(width = 0.75) +
      fill_scale(nlevels(out$cell_type), "cell type") +
      scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(c(0, 0.02))) +
      coord_flip() + labs(x = input$cmpvar, y = "Fraction of cells") + base_theme
    if (!identical(input$cmpfacet, "(none)")) p <- p + facet_wrap(~.f)
    p
  }

  cmp_h <- function(mult = 26, base = 140) {
    n <- n_distinct(cmp_data()$.g)
    f <- if (identical(input$cmpfacet, "(none)")) 1 else n_distinct(cmp_data()$.f)
    max(300, min(1600, base + mult * n * ceiling(f / 2)))
  }
  output$cmp_bar_ui  <- renderUI(plotOutput("cmp_bar",  height = paste0(cmp_h(), "px")))
  output$cmp_box_ui  <- renderUI(plotOutput("cmp_box",  height = paste0(cmp_h(), "px")))
  output$cmp_comp_ui <- renderUI(plotOutput("cmp_comp", height = paste0(cmp_h(), "px")))
  output$cmp_bar  <- renderPlot(cmp_bar_plot(),  res = 100)
  output$cmp_box  <- renderPlot(cmp_box_plot(),  res = 100)
  output$cmp_comp <- renderPlot(cmp_comp_plot(), res = 100)
  output$cmp_tbl <- renderDT({
    s <- cmp_summary(); names(s)[1] <- input$cmpvar
    datatable(s, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, dom = "Bfrtip", buttons = c("csv", "excel"),
                             scrollX = TRUE)) %>%
      formatCurrency(c("n_cells", "n_samples"), currency = "", digits = 0)
  })

  # ---- metadata column overview ---------------------------------------------
  ov_data <- reactive({
    if (identical(input$ov_scope, "whole atlas")) atlas else dat()
  })

  # one row per column: type, number of levels, missingness, and a value preview
  col_overview <- reactive({
    d <- ov_data()
    validate(need(nrow(d) > 0, "No samples match the current filters."))
    cols <- names(d)
    rows <- lapply(cols, function(cl) {
      v <- d[[cl]]
      n_miss <- sum(is.na(v))
      if (is.numeric(v)) {
        q <- suppressWarnings(quantile(v, c(0, 0.5, 1), na.rm = TRUE))
        summ <- if (all(is.na(v))) "all missing" else
          paste0("min ", fmt_num(q[[1]]), " | median ", fmt_num(q[[2]]),
                 " | max ", fmt_num(q[[3]]))
        n_lev <- length(unique(stats::na.omit(v)))
      } else {
        vc <- sort(table(as.character(v)), decreasing = TRUE)
        n_lev <- length(vc)
        top <- utils::head(vc, 3)
        summ <- if (!n_lev) "all missing" else
          paste0(paste0(names(top), " (", as.integer(top), ")", collapse = ", "),
                 if (n_lev > 3) paste0(", +", n_lev - 3, " more") else "")
      }
      data.frame(column = cl, group = col_group(cl),
                 type = if (is.numeric(v)) "numeric" else
                        if (is.logical(v)) "logical" else "categorical",
                 n_levels = n_lev, n_missing = n_miss,
                 pct_missing = round(100 * n_miss / nrow(d), 1),
                 values = summ, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })

  observeEvent(col_overview(), {
    cur <- isolate(input$ov_col)
    choices <- split(col_overview()$column, col_overview()$group)
    updateSelectizeInput(session, "ov_col", choices = choices,
                         selected = if (!is.null(cur) && nzchar(cur)) cur else "medical_condition",
                         server = FALSE)
  }, once = TRUE)

  observeEvent(input$ov_tbl_rows_selected, {
    updateSelectizeInput(session, "ov_col",
                         selected = col_overview()$column[input$ov_tbl_rows_selected])
  })

  ov_col <- reactive({
    cl <- input$ov_col
    if (is.null(cl) || !nzchar(cl) || !cl %in% names(atlas)) "medical_condition" else cl
  })

  output$ov_detail_header <- renderText({
    cl <- ov_col()
    row <- col_overview()[col_overview()$column == cl, ]
    paste0(cl, "  -  ", row$type, ", ", row$n_levels, " levels, ",
           row$n_missing, " of ", nrow(ov_data()), " samples not annotated")
  })

  ov_detail <- reactive({
    d <- ov_data(); cl <- ov_col(); v <- d[[cl]]
    if (is.numeric(v)) {
      s <- c(n = sum(!is.na(v)), missing = sum(is.na(v)),
             min = suppressWarnings(min(v, na.rm = TRUE)),
             q1 = unname(quantile(v, 0.25, na.rm = TRUE)),
             median = stats::median(v, na.rm = TRUE),
             q3 = unname(quantile(v, 0.75, na.rm = TRUE)),
             max = suppressWarnings(max(v, na.rm = TRUE)),
             mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE))
      data.frame(statistic = names(s), value = vapply(s, fmt_num, character(1)),
                 stringsAsFactors = FALSE, row.names = NULL)
    } else {
      out <- d %>% mutate(.v = lab_na(.data[[cl]])) %>% group_by(.v) %>%
        summarise(n_samples = dplyr::n(), n_patients = n_distinct(patient_id),
                  n_studies = n_distinct(study_id), n_cells = sum(n_cells), .groups = "drop") %>%
        mutate(pct_samples = round(100 * n_samples / sum(n_samples), 1)) %>%
        arrange(desc(n_samples)) %>%
        relocate(pct_samples, .after = n_samples)
      names(out)[1] <- cl
      as.data.frame(out)
    }
  })

  ov_plot <- function() {
    d <- ov_data(); cl <- ov_col(); v <- d[[cl]]
    if (is.numeric(v)) {
      validate(need(any(!is.na(v)), paste(cl, "has no values in this selection.")))
      return(ggplot(d[!is.na(d[[cl]]), , drop = FALSE], aes(x = .data[[cl]])) +
               geom_histogram(bins = 30, fill = PALETTE[1], colour = "white", linewidth = 0.2) +
               scale_y_continuous(expand = expansion(c(0, 0.05))) +
               labs(x = cl, y = "Samples") + base_theme)
    }
    tb <- ov_detail()
    names(tb)[1] <- ".v"
    capped <- utils::head(tb, 40)
    capped$.v <- factor(capped$.v, levels = rev(capped$.v))
    ggplot(capped, aes(x = .v, y = n_samples)) +
      geom_col(fill = PALETTE[1], width = 0.78) +
      geom_text(aes(label = n_samples), hjust = -0.15, size = 3) +
      scale_y_continuous(expand = expansion(c(0, 0.15))) +
      coord_flip() +
      labs(x = NULL, y = "Samples",
           subtitle = if (nrow(tb) > 40) paste("showing the 40 most frequent of",
                                               nrow(tb), "levels") else NULL) +
      base_theme
  }

  output$ov_plot_ui <- renderUI({
    v <- ov_data()[[ov_col()]]
    h <- if (is.numeric(v)) 300 else
      max(220, min(900, 70 + 20 * min(40, length(unique(lab_na(v))))))
    plotOutput("ov_plot", height = paste0(h, "px"))
  })
  output$ov_plot <- renderPlot(ov_plot(), res = 100)

  output$ov_tbl <- renderDT(
    datatable(col_overview(), rownames = FALSE, filter = "top",
              selection = list(mode = "single", selected = NULL),
              options = list(pageLength = 15, scrollX = TRUE,
                             columnDefs = list(list(width = "260px", targets = 6))))
  )
  output$ov_detail_tbl <- renderDT(
    datatable(ov_detail(), rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 10, dom = "Bfrtip", buttons = c("csv", "excel"),
                             scrollX = TRUE))
  )
  output$dl_overview <- downloadHandler(
    filename = function() paste0("crc_atlas_column_overview_", Sys.Date(), ".csv"),
    content = function(file) write.csv(col_overview(), file, row.names = FALSE)
  )

  # ---- study overview -------------------------------------------------------
  study_summary <- reactive({
    d <- dat()
    validate(need(nrow(d) > 0, "No samples match the current filters."))
    collapse_u <- function(x) paste(sort(unique(na.omit(as.character(x)))), collapse = ", ")
    d %>% group_by(study_id) %>%
      summarise(n_datasets = n_distinct(dataset), n_samples = dplyr::n(),
                n_patients = n_distinct(patient_id),
                median_cells_per_sample = median(n_cells),
                n_cells = sum(n_cells),
                sample_types = collapse_u(sample_type), tissues = collapse_u(sample_tissue),
                platforms = collapse_u(platform), countries = collapse_u(country),
                pmid = collapse_u(study_pmid), doi = collapse_u(study_doi),
                .groups = "drop") %>%
      relocate(n_cells, .before = median_cells_per_sample) %>%
      arrange(desc(n_cells))
  })
  output$study_tbl <- renderDT(
    datatable(study_summary(), rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  )
  output$dl_studies <- downloadHandler(
    filename = function() paste0("crc_atlas_study_overview_", Sys.Date(), ".csv"),
    content = function(file) write.csv(study_summary(), file, row.names = FALSE)
  )

  # ---- sample table ---------------------------------------------------------
  output$sample_tbl <- renderDT({
    d <- dat()
    cols <- intersect(input$tbl_cols, names(d))
    if (!length(cols)) cols <- c("study_id", "sample_id", "n_cells")
    datatable(d[, cols, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })
  output$dl_filtered <- downloadHandler(
    filename = function() paste0("crc_atlas_filtered_samples_", Sys.Date(), ".csv"),
    content = function(file) write.csv(dat(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
