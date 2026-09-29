Step 3: NSTI Quality Assessment
================

- [Introduction](#introduction)
  - [Purpose](#purpose)
  - [What Is NSTI?](#what-is-nsti)
  - [Prerequisites](#prerequisites)
  - [What This Notebook Does](#what-this-notebook-does)
  - [Statistics and Calculations](#statistics-and-calculations)
  - [Expected Input](#expected-input)
  - [Expected Output](#expected-output)
- [Environment Setup](#environment-setup)
  - [Load Required Packages](#load-packages)
- [Configuration](#configuration)
  - [Adjust Run Parameters](#adjust-run-parameters)
  - [Define Path Parameters](#define-paths)
- [Data Import](#data-import)
  - [Import Sample Metadata](#import-metadata)
  - [Import Per-ASV NSTI Values](#import-per-asv-nsti)
  - [Import Per-Sample Weighted NSTI Values](#import-weighted-nsti)
  - [Import Step 1 ASV Abundance Table](#import-asv-abundance)
  - [Validate Sample ID Consistency](#validate-sample-ids)
- [Cross-Reference Step 2’s NSTI Filtering Log](#cross-reference-log)
- [Compute QC Metrics](#compute-qc-metrics)
  - [Per-ASV NSTI Table](#compute-per-asv)
  - [Per-Sample NSTI Table](#compute-per-sample)
- [Visualizations](#visualizations)
  - [Boxplot: Weighted NSTI by Condition](#plot-boxplot)
  - [Histogram/Density: Per-ASV NSTI Distribution](#plot-histogram)
  - [Scatterplot: NSTI vs. Abundance](#plot-scatter)
- [Results Summary](#results-summary)
  - [Compute Filtering Log Summary Table](#summary-log-table)
  - [Display Per-Sample and Per-ASV Tables](#summary-tables)
  - [Export Summary to Excel](#export-summary)
  - [Document and Export Column Dictionary](#column-dictionary)
- [Output File Summary](#output-file-summary)
- [Output Interpretation](#output-interpretation)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [NSTI and PICRUSt2](#nsti-and-picrust2)
  - [R Packages](#r-packages)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [data/metadata.tsv Not Found](#datametadatatsv-not-found)
    - [Sample ID Mismatches](#sample-id-mismatches)
    - [combined_marker_predicted_and_nsti.tsv or `weighted_nsti.tsv` Not
      Found](#combined_marker_predicted_and_nstitsv-or-weighted_nstitsv-not-found)
    - [Log Cross-Check Count Mismatch](#log-cross-check-count-mismatch)
    - [Suspiciously High NSTI Values Across the Entire
      Dataset](#suspiciously-high-nsti-values-across-the-entire-dataset)
    - [`mean_relative_abundance` Shows `#NUM!` for (Nearly) Every ASV in
      Excel](#mean_relative_abundance-shows-num-for-nearly-every-asv-in-excel)
    - [BIOM Table Fails to Import as a Numeric
      Matrix](#biom-table-fails-to-import-as-a-numeric-matrix)

<!-- The hidden setup chunk below controls whether analysis code is executed. HTML reports run the workflow; GitHub Markdown remains a non-executing tutorial and code reference. -->
<!-- The hidden CSS chunk below affects only the rendered HTML report and is deliberately excluded from GitHub Markdown. -->

# Introduction

## Purpose

This notebook is **Step 3** of the
[PICRUSt2](https://github.com/picrust/picrust2/wiki) workflow. Before
drawing any biological conclusions from the functional predictions
produced in [Step 2](2_picrust2_pipeline.md), it is essential to check
*how confident* those predictions actually are. This notebook visualizes
and tabulates the **Nearest Sequenced Taxon Index (NSTI)** — PICRUSt2’s
own built-in measure of prediction reliability — for every ASV and every
sample, so that questionable data can be identified before it propagates
into downstream differential-abundance interpretation.

## What Is NSTI?

The **Nearest Sequenced Taxon Index (NSTI)** quantifies, for a single
input sequence (ASV), the phylogenetic distance between where that
sequence was placed in PICRUSt2’s reference tree and the nearest
reference genome with actually-sequenced gene content. Concretely, it is
the sum of branch lengths separating the ASV’s placement point from its
nearest sequenced relative in the tree.

- **Low NSTI** (close to 0): the ASV sits very close to a reference
  genome, so the inferred gene content is based on a close relative —
  high confidence.
- **High NSTI**: the nearest reference genome is phylogenetically
  distant, so the inferred gene content is an extrapolation over a
  larger evolutionary distance — lower confidence.

NSTI is computed automatically during [Step 2](2_picrust2_pipeline.md)’s
hidden-state prediction and is already used internally: Step 2 passes
its documented `max_nsti_threshold` explicitly to
`picrust2_pipeline.py`, excluding any ASV above that threshold from the
KO/EC metagenome predictions before this notebook ever sees the data.
The workflow default is 2.0 and can be changed consistently with the
`PICRUST2_MAX_NSTI` environment variable. This notebook does not
re-implement the PICRUSt2 filter — it visualizes and cross-checks what
the pipeline did using PICRUSt2’s own output files.

PICRUSt2 itself writes two related but distinct NSTI outputs, and this
notebook uses both rather than recomputing either from scratch:

<div class="nsti-source-table">

| File | Level | What it reflects |
|:---|:---|:---|
| [combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv) | Per-ASV | The NSTI value for **every** placed ASV, including any that [Step 2](2_picrust2_pipeline.md) subsequently excluded for being above the `max_nsti` cut-off. This is the complete, unfiltered picture. |
| [KO_metagenome_out/weighted_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/KO_metagenome_out/weighted_nsti.tsv) | Per-sample | A single abundance-weighted mean NSTI value per sample, computed by PICRUSt2 **after** already dropping ASVs above the `max_nsti` cut-off. |

</div>

<div class="alert alert-info">

**Why use the `combined_` table rather than a `bac_` or `arc_` table?**
PICRUSt2 places sequences against bacterial and archaeal reference trees
separately, writing
[bac_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/bac_marker_predicted_and_nsti.tsv)
and
[arc_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/arc_marker_predicted_and_nsti.tsv).
[combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv)
merges the selected placements into one non-redundant per-ASV table and
records the selected domain in `best_domain`. It is also the file read
by PICRUSt2’s metagenome-prediction stage, so it is the appropriate
source for assessing the ASVs that feed the community predictions.

</div>

<div class="alert alert-warning">

**This distinction matters for interpretation.** Because
`weighted_nsti.tsv` is computed only over the ASVs that *survived* [Step
2](2_picrust2_pipeline.md)’s filtering, a sample can show a reassuringly
low per-sample weighted NSTI even if a substantial fraction of its reads
belonged to ASVs that were excluded entirely. This notebook therefore
also reports, per sample, what fraction of reads were excluded by the
`max_nsti` filter — a number the official `weighted_nsti` file cannot
tell you on its own.

</div>

There is no universal, field-wide cut-off for what counts as a “bad”
*per-sample weighted* NSTI value (unlike the configured per-ASV
`max_nsti` filter) — acceptable values vary by environment. This
notebook therefore flags unusually high per-sample values using a
data-driven rule (values above the third quartile plus 1.5 times the
interquartile range by default), so the flag is relative to the study’s
own distribution rather than an invented universal threshold.

## Prerequisites

Before running this notebook, ensure that:

1.  Required R packages are installed — run
    [setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
    once per R environment if you have not already done so.
2.  **[Step 1](1_prepare_picrust2_inputs.md)** has been completed,
    producing
    [study_seqs.biom](../../results/1_prepare_picrust2_inputs/study_seqs.biom)
    in
    [results/1_prepare_picrust2_inputs/](../../results/1_prepare_picrust2_inputs/).
3.  **[Step 2](2_picrust2_pipeline.md)** has been completed, with
    `decompress_output_files <- TRUE` (the default), producing
    [combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv)
    and
    [KO_metagenome_out/weighted_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/KO_metagenome_out/weighted_nsti.tsv)
    in [its output
    directory](../../results/2_picrust2_pipeline/picrust2_out_pipeline/).
4.  **[data/metadata.tsv](../../data/metadata.tsv)** is present with the
    required `SampleID`, `Condition`, and `Reference` columns (see
    [data/README.md](../../data/README.md)). This file is new as of this
    step — it was not required for Steps
    [1](1_prepare_picrust2_inputs.md)-[2](2_picrust2_pipeline.md).

## What This Notebook Does

1.  **Imports sample metadata** and validates that every sample ID
    matches the PICRUSt2 output being analyzed.
2.  **Imports per-ASV NSTI values** directly from PICRUSt2’s own
    [combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv)
    output.
3.  **Imports per-sample weighted NSTI values** directly from PICRUSt2’s
    own `weighted_nsti.tsv` output.
4.  **Cross-references [Step 2](2_picrust2_pipeline.md)’s run log** to
    report how many ASVs the pipeline’s own `max_nsti` filter already
    excluded.
5.  **Computes a per-sample “reads excluded by NSTI filter” metric** —
    the complementary number the official `weighted_nsti` file cannot
    show on its own (see the warning above).
6.  **Produces three interactive (plotly) visualizations**: a per-sample
    weighted NSTI boxplot by condition, a per-ASV NSTI distribution
    histogram/density plot, and an NSTI-vs-abundance scatterplot.
7.  **Exports every table to a shared Excel workbook**, with each
    interactive plot additionally saved as a standalone, shareable HTML
    file.

## Statistics and Calculations

- **Per-ASV NSTI** is read directly from PICRUSt2; lower values indicate
  a closer sequenced reference. `exceeds_max_nsti_threshold` reproduces
  the configured per-ASV filtering rule as a transparent flag.
- **Per-sample weighted NSTI** is PICRUSt2’s abundance-weighted mean
  over ASVs retained by that filter; this notebook does not recompute or
  alter it.
- **Excluded-read percentage** is `100 ×` the within-sample relative
  abundance assigned to ASVs above the NSTI threshold. It complements
  weighted NSTI by quantifying reads that were omitted before the
  weighted mean was calculated.
- **High-value IQR-fence flag** uses Tukey’s upper fence,
  `Q3 + iqr_fence_multiplier × IQR` (default multiplier 1.5). It is a
  descriptive screening rule, not a hypothesis test, confidence
  interval, or p-value. No between-condition significance test is
  performed in this step.

## Expected Input

- [results/1_prepare_picrust2_inputs/study_seqs.biom](../../results/1_prepare_picrust2_inputs/study_seqs.biom)
  — ASV x sample count table from Step 1 (used to compute per-ASV mean
  relative abundance and per-sample excluded-read fractions).
- [results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv)
  — per-ASV NSTI from Step 2. This notebook expects it already
  decompressed (i.e. Step 2 was run with
  `decompress_output_files <- TRUE`, the default).
- [results/2_picrust2_pipeline/picrust2_out_pipeline/KO_metagenome_out/weighted_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/KO_metagenome_out/weighted_nsti.tsv)
  — per-sample weighted NSTI from Step 2 (also expected already
  decompressed).
- [data/metadata.tsv](../../data/metadata.tsv) — sample metadata with
  exactly three required columns: `SampleID`, `Condition`, and
  `Reference`.

## Expected Output

- [results/3_nsti_quality_assessment/nsti_quality_summary.xlsx](../../results/3_nsti_quality_assessment/nsti_quality_summary.xlsx)
  — workbook with `Per_Sample_NSTI`, `Per_ASV_NSTI`, and
  `NSTI_Filtering_Log_Summary` sheets, plus a trailing
  `Column_Dictionary` sheet documenting every column in the workbook.
- [results/3_nsti_quality_assessment/plots/](../../results/3_nsti_quality_assessment/plots/)
  — three standalone interactive HTML plots (weighted NSTI boxplot,
  per-ASV NSTI distribution, NSTI vs. abundance scatterplot).

------------------------------------------------------------------------

# Environment Setup

<div class="alert alert-info">

**Before running this notebook**: install the required R packages by
running
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
once per R environment. This only needs to be done once — not before
every run of this notebook.

</div>

## Load Required Packages

The chunk below loads every R package this notebook depends on, plus
three helper functions shared across the whole pipeline: Excel writing,
column-dictionary documentation, and clickable output links (each
sourced from [R/functions/](../functions/)).

``` r
# NOTE: This chunk assumes the packages below are already installed. If any
# library() call fails with "there is no package called ...", run
# setup/install_R_dependencies.R once to install every package this
# workflow needs, then re-run this notebook.

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# openxlsx: Excel file creation and manipulation
library(openxlsx)

# readr: Fast, consistent delimited-file import; transparently handles both
# plain-text and gzip-compressed (.gz) files without manual decompression
library(readr)

# dplyr: Data manipulation grammar (filter/mutate/summarise/join)
library(dplyr)

# tidyr: Reshaping data between wide and long formats
library(tidyr)

# biomformat: Bioconductor package for reading BIOM-format files
# Used here to read back the ASV x sample count table written in Step 1
library(biomformat)

# ggplot2: Grammar-of-graphics static plotting; the base layer converted to
# an interactive widget by plotly::ggplotly() below
library(ggplot2)

# plotly: Converts ggplot2 plots into interactive (zoomable, hoverable)
# HTML widgets
library(plotly)

# htmlwidgets: Saves plotly widgets as standalone, self-contained HTML files
library(htmlwidgets)

# DT: interactive, paginated/scrollable HTML tables — used below so
# per-sample/per-column tables are not dumped into the HTML report in full.
library(DT)

# Source the custom Excel utility function from the project's function library
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source the custom column-dictionary utility function from the project's
# function library. This function documents every column of a sheet being
# exported, so a trailing "Column_Dictionary" sheet can be added to the
# workbook (see the Document and Export Column Dictionary section below).
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source the custom output-links utility function from the project's function
# library. This function renders a chunk's output file path(s) as clickable
# Markdown links in the knitted report, for reproducibility across machines.
source(here("R", "functions", "render_output_links_function.R"))

# Source the custom render_output_tree utility function from the project's
# function library. This function prints this notebook's entire output
# folder as a clickable directory tree (used in "Output File Summary" below),
# scanned live from disk at knit time so it always matches what was actually
# produced on this run.
source(here("R", "functions", "render_output_tree_function.R"))
```

------------------------------------------------------------------------

# Configuration

## Adjust Run Parameters

Set the per-ASV NSTI flagging threshold and the per-sample outlier
sensitivity here before running the rest of the notebook – both default
to values already explained above and rarely need changing unless Step 2
was run with a non-default `--max_nsti`.

``` r
# ------------------------------------------------------------------
# Per-ASV NSTI cut-off used for flagging in this notebook's tables/plots.
# Uses the same environment-variable override as Step 2 and Step 5. The
# workflow default is 2.0; setting PICRUST2_MAX_NSTI keeps all three steps in
# sync without editing multiple notebooks.
# ------------------------------------------------------------------
max_nsti_text <- Sys.getenv("PICRUST2_MAX_NSTI", unset = "2.0")
max_nsti_threshold <- suppressWarnings(as.numeric(max_nsti_text))
if (length(max_nsti_threshold) != 1L || is.na(max_nsti_threshold) ||
    !is.finite(max_nsti_threshold) || max_nsti_threshold < 0) {
  stop("PICRUST2_MAX_NSTI must be one finite non-negative number; received: ",
       max_nsti_text)
}

# ------------------------------------------------------------------
# Multiplier used for the descriptive Tukey/IQR upper-fence rule applied
# to per-sample weighted NSTI values (see the Purpose section above for
# why no universal fixed threshold is used here). 1.5 is the conventional
# statistical default for a moderate outlier; increase it (e.g. to 3) to
# flag only more extreme outliers.
# ------------------------------------------------------------------
iqr_fence_multiplier <- 1.5  # <-- EDIT to make fence exceedance flagging stricter/looser
if (length(iqr_fence_multiplier) != 1L || is.na(iqr_fence_multiplier) ||
    !is.finite(iqr_fence_multiplier) || iqr_fence_multiplier < 0) {
  stop("iqr_fence_multiplier must be one finite non-negative number.")
}

cat(
  "NSTI QC Configuration:\n",
  "- Per-ASV max_nsti flagging threshold:", max_nsti_threshold, "\n",
  "- Per-sample Tukey/IQR fence multiplier:", iqr_fence_multiplier, "\n"
)
```

## Define Path Parameters

This chunk resolves every input and output path relative to the project
root via `here()`, creates this step’s own output folder(s) if they do
not already exist, and fails fast with an informative error if any
expected upstream file cannot be found – so a missing prerequisite is
caught immediately, before any computation begins.

``` r
# Base folders, mirroring the layout established in Steps 1-2. The example
# runner overrides both so its test data and outputs remain isolated.
data_folder    <- here(Sys.getenv("PICRUST2_DATA_DIR", unset = "data"))
results_folder <- here(Sys.getenv("PICRUST2_RESULTS_DIR", unset = "results"))

# Step 1's output (this notebook's ASV abundance input)
step1_output_folder <- here(results_folder, "1_prepare_picrust2_inputs")
input_biom_path      <- here(step1_output_folder, "study_seqs.biom")

# Step 2's single consolidated results folder (see 2_picrust2_pipeline.Rmd),
# with its own logs subfolder
step2_output_folder <- here(results_folder, "2_picrust2_pipeline")
run_output_dir       <- here(step2_output_folder, "picrust2_out_pipeline")
step2_log_folder     <- here(step2_output_folder, "logs")
run_log_path          <- here(step2_log_folder, "picrust2_run.log")

# This notebook's own output folder
output_folder <- here(results_folder, "3_nsti_quality_assessment")
plots_folder  <- here(output_folder, "plots")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_folder, recursive = TRUE, showWarnings = FALSE)

# Path to the metadata file (new as of this step; see data/README.md)
input_metadata_path <- here(data_folder, "metadata.tsv")

# Path to the Excel workbook summarizing this step's QC results
output_excel_path <- here(output_folder, "nsti_quality_summary.xlsx")

if (!dir.exists(run_output_dir)) {
  stop(
    "Step 2's run output directory was not found:\n  ", run_output_dir, "\n",
    "Run Step 2 (2_picrust2_pipeline.Rmd) first."
  )
}

# ------------------------------------------------------------------
# Helper: confirm an expected input file actually exists and return its
# path, or stop with a clear, actionable error if it does not. This
# notebook expects Step 2's output to already be decompressed (i.e.
# decompress_output_files <- TRUE, the default in Step 2's Configuration) —
# it does not look for .gz versions of these files.
# ------------------------------------------------------------------
resolve_existing_file <- function(expected_path) {
  if (!file.exists(expected_path)) {
    stop(
      "Could not find the expected file:\n  ", expected_path, "\n",
      "This notebook expects Step 2's output to already be decompressed ",
      "(decompress_output_files <- TRUE in Step 2's Configuration, the default). ",
      "If Step 2 was run with decompression turned off, either re-run it with ",
      "decompress_output_files <- TRUE, or decompress this specific file manually."
    )
  }
  expected_path
}

# PICRUSt2 places sequences against separate bacterial and archaeal
# reference trees, writing "bac_marker_predicted_and_nsti.tsv" and
# "arc_marker_predicted_and_nsti.tsv" respectively. "combined_..." merges
# these into one non-redundant per-ASV table. This is also the exact file
# PICRUSt2's own metagenome_pipeline.py reads before generating the KO/EC
# predictions, so using it here matches Step 2's actual ASV coverage.
input_per_asv_nsti_path <- resolve_existing_file(
  file.path(run_output_dir, "combined_marker_predicted_and_nsti.tsv")
)

input_weighted_nsti_path <- resolve_existing_file(
  file.path(run_output_dir, "KO_metagenome_out", "weighted_nsti.tsv")
)

if (!file.exists(input_metadata_path)) {
  stop(
    "Could not find the sample metadata file: ", input_metadata_path, "\n",
    "This file is required starting from Step 3 — see data/README.md for the ",
    "expected format (SampleID, Condition, Reference)."
  )
}

if (!file.exists(input_biom_path)) {
  stop(
    "Could not find the BIOM input file: ", input_biom_path, "\n",
    "Run Step 1 (1_prepare_picrust2_inputs.Rmd) first."
  )
}

cat(
  "Path Configuration:\n",
  "- Step 2 run directory used:", run_output_dir, "\n",
  "- Per-ASV NSTI file:", input_per_asv_nsti_path, "\n",
  "- Per-sample weighted NSTI file:", input_weighted_nsti_path, "\n",
  "- Metadata file:", input_metadata_path, "\n",
  "- Step 1 BIOM table:", input_biom_path, "\n"
)
```

------------------------------------------------------------------------

# Data Import

## Import Sample Metadata

This chunk reads [data/metadata.tsv](../../data/metadata.tsv) and
confirms that the simplified `SampleID`, `Condition`, and `Reference`
schema is valid before any results are joined to it.

``` r
# Read the metadata file (tab-separated; see data/README.md for the schema)
sample_metadata <- read_tsv(input_metadata_path, show_col_types = FALSE)

# Confirm the required columns are present before proceeding
required_metadata_columns <- c("SampleID", "Condition", "Reference")
missing_metadata_columns <- setdiff(required_metadata_columns, colnames(sample_metadata))
if (length(missing_metadata_columns) > 0) {
  stop(
    "data/metadata.tsv is missing required column(s): ",
    paste(missing_metadata_columns, collapse = ", "),
    ". See data/README.md for the expected format."
  )
}
unexpected_metadata_columns <- setdiff(
  colnames(sample_metadata), required_metadata_columns
)
if (length(unexpected_metadata_columns) > 0L) {
  stop(
    "data/metadata.tsv must use the simplified three-column schema exactly. ",
    "Unexpected column(s): ", paste(unexpected_metadata_columns, collapse = ", "),
    ". Expected: ", paste(required_metadata_columns, collapse = ", ")
  )
}

# Enforce one row per non-empty sample identifier before any joins. Duplicate
# IDs could multiply rows; missing conditions would make a group summary
# uninterpretable even if the join itself succeeded.
if (anyDuplicated(sample_metadata$SampleID)) {
  stop("data/metadata.tsv contains duplicate SampleID values; each sample must occur exactly once.")
}
if (any(is.na(sample_metadata$SampleID) | trimws(sample_metadata$SampleID) == "")) {
  stop("data/metadata.tsv contains an empty SampleID value.")
}
if (any(is.na(sample_metadata$Condition) | trimws(sample_metadata$Condition) == "")) {
  stop("data/metadata.tsv contains an empty Condition value.")
}

# Normalize common logical encodings and reject anything that cannot be
# interpreted unambiguously as TRUE or FALSE.
reference_flag <- if (is.logical(sample_metadata$Reference)) {
  sample_metadata$Reference
} else {
  as.logical(sample_metadata$Reference)
}
if (any(is.na(reference_flag))) {
  stop("data/metadata.tsv's Reference column must contain only TRUE or FALSE.")
}
sample_metadata$Reference <- reference_flag

# A reference designation belongs to an entire condition, not to selected
# samples within that condition. First verify within-condition consistency,
# then require exactly one complete reference condition.
reference_flags_by_condition <- sample_metadata %>%
  group_by(Condition) %>%
  summarise(n_flags = n_distinct(Reference), .groups = "drop")
if (any(reference_flags_by_condition$n_flags != 1)) {
  stop("Every sample in a Condition must have the same Reference value.")
}
reference_conditions <- unique(sample_metadata$Condition[sample_metadata$Reference])
if (length(reference_conditions) != 1) {
  stop("Exactly one Condition must be marked Reference = TRUE.")
}

cat(
  "Imported metadata for", nrow(sample_metadata), "samples.\n",
  "- Conditions found:", paste(unique(sample_metadata$Condition), collapse = ", "), "\n",
  "- Reference condition:", reference_conditions, "\n"
)
```

## Import Per-ASV NSTI Values

Reads PICRUSt2’s own
[combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv)
output directly — this is the complete, unfiltered per-ASV NSTI table
(i.e. it still includes ASVs that Step 2 subsequently excluded via the
`max_nsti` cut-off).

``` r
# read_tsv() reads the plain-text .tsv file directly.
per_asv_nsti_raw <- read_tsv(input_per_asv_nsti_path, show_col_types = FALSE)

# Confirm the columns this notebook depends on are present. PICRUSt2 writes
# the ASV/sequence identifier as "sequence" and the NSTI value as
# "metadata_NSTI" (confirmed against both PICRUSt2's own
# metagenome_pipeline.py source and a real combined_marker_predicted_and_nsti.tsv
# file, which also includes "16S_rRNA_Count", "best_domain", and
# "closest_reference_genome" columns this notebook does not need) — if
# either required column is missing, the file is not the expected format.
required_nsti_columns <- c("sequence", "metadata_NSTI")
missing_nsti_columns <- setdiff(required_nsti_columns, colnames(per_asv_nsti_raw))
if (length(missing_nsti_columns) > 0) {
  stop(
    "combined_marker_predicted_and_nsti.tsv is missing expected column(s): ",
    paste(missing_nsti_columns, collapse = ", "),
    ". Found columns: ", paste(colnames(per_asv_nsti_raw), collapse = ", ")
  )
}

# Keep only the columns this notebook needs, renamed for consistency with
# the ASV_ID naming convention used elsewhere in this repository
per_asv_nsti <- per_asv_nsti_raw %>%
  select(ASV_ID = sequence, NSTI = metadata_NSTI)

cat("Imported per-ASV NSTI values for", nrow(per_asv_nsti), "ASVs.\n")
```

## Import Per-Sample Weighted NSTI Values

Reads PICRUSt2’s own `weighted_nsti.tsv` output — the abundance-weighted
mean NSTI per sample, computed by PICRUSt2 *after* excluding ASVs above
the `max_nsti` cut-off.

``` r
weighted_nsti_raw <- read_tsv(input_weighted_nsti_path, show_col_types = FALSE)

# PICRUSt2 writes the sample identifier column as "sample" and the value
# column as "weighted_NSTI" (confirmed against calc_weighted_nsti() in
# PICRUSt2's own metagenome_pipeline.py source).
required_weighted_columns <- c("sample", "weighted_NSTI")
missing_weighted_columns <- setdiff(required_weighted_columns, colnames(weighted_nsti_raw))
if (length(missing_weighted_columns) > 0) {
  stop(
    "weighted_nsti file is missing expected column(s): ",
    paste(missing_weighted_columns, collapse = ", "),
    ". Found columns: ", paste(colnames(weighted_nsti_raw), collapse = ", ")
  )
}

weighted_nsti <- weighted_nsti_raw %>%
  rename(SampleID = sample, weighted_NSTI = weighted_NSTI)

cat("Imported per-sample weighted NSTI values for", nrow(weighted_nsti), "samples.\n")
```

## Import Step 1 ASV Abundance Table

Reads back the BIOM feature table written in [Step
1](1_prepare_picrust2_inputs.md), used below to compute each ASV’s mean
relative abundance (for the NSTI-vs-abundance scatterplot) and the
fraction of each sample’s reads excluded by the `max_nsti` filter.

``` r
# read_biom() returns a biom-class object. biom_data() is biomformat's own
# accessor for the underlying ASV (rows) x sample (columns) count matrix,
# matching the orientation written by make_biom() in Step 1 -- used here
# instead of as(biom_object, "matrix") because that S4 coercion is not
# reliably registered across biomformat versions and can silently return a
# non-numeric object (this is what produced the "'x' must be numeric or
# complex" error from colSums() below). biom_data() can itself return
# either a dense matrix or a sparse Matrix (dgTMatrix) depending on how the
# BIOM file was written, so as.matrix() is applied on top to guarantee a
# plain, dense, numeric matrix either way.
asv_count_matrix <- as.matrix(biom_data(read_biom(input_biom_path)))

# Defensive check: confirm the matrix actually came out numeric before it
# is used in any arithmetic below. A non-numeric result here would point
# to an unexpected BIOM file structure rather than a routine input problem,
# so this fails loudly with the actual storage mode reported.
if (!is.numeric(asv_count_matrix)) {
  stop(
    "The BIOM feature table did not import as a numeric matrix (storage mode: ",
    storage.mode(asv_count_matrix), "). Confirm study_seqs.biom from Step 1 ",
    "was written correctly -- re-running Step 1 may help."
  )
}

# ------------------------------------------------------------------
# Detect samples with zero total ASV counts (e.g. an empty negative/blank
# control that yielded no amplicon reads). Relative abundance is
# mathematically undefined (0/0) for such a sample -- if left unguarded,
# that single all-NaN sample column silently propagates into *every* ASV's
# mean_relative_abundance downstream (rowMeans() with the default
# na.rm = FALSE returns NaN for a row if *any* element in that row is NaN),
# which openxlsx then writes to the exported workbook as "#NUM!" for
# essentially the entire Per_ASV_NSTI column -- surfacing it here up front
# keeps the actual cause visible instead of an opaque spreadsheet error.
# ------------------------------------------------------------------
sample_total_counts <- colSums(asv_count_matrix)
zero_count_samples   <- names(sample_total_counts)[sample_total_counts == 0]

if (length(zero_count_samples) > 0) {
  warning(
    length(zero_count_samples), " sample(s) have zero total ASV counts in Step 1's ",
    "BIOM table, so relative abundance is undefined (0/0) for them: ",
    paste(zero_count_samples, collapse = ", "), ". ",
    "This is expected for an empty negative/blank control but would be unexpected ",
    "for a genuine biological sample -- confirm which case applies here. These ",
    "sample(s) are excluded from the mean_relative_abundance and ",
    "pct_reads_excluded_by_nsti_filter calculations below (rather than silently ",
    "corrupting every ASV's mean), but remain in Per_Sample_NSTI with PICRUSt2's ",
    "own reported weighted_NSTI for them."
  )
}

# Convert to relative abundance within each sample (columns sum to 1). The
# zero_count_samples column(s) identified above become entirely NaN here
# (0/0) -- this is expected and is guarded against explicitly downstream
# rather than fixed at the matrix level, so the raw NaN remains visible to
# anyone inspecting asv_relative_abundance_matrix directly.
asv_relative_abundance_matrix <- sweep(asv_count_matrix, MARGIN = 2, STATS = colSums(asv_count_matrix), FUN = "/")

cat(
  "Imported Step 1 ASV abundance table:\n",
  "- ASVs:", nrow(asv_count_matrix), "\n",
  "- Samples:", ncol(asv_count_matrix), "\n",
  "- Samples with zero total counts:", length(zero_count_samples),
  if (length(zero_count_samples) > 0) paste0(" (", paste(zero_count_samples, collapse = ", "), ")") else "", "\n"
)
```

## Validate Sample ID Consistency

Confirm that the metadata, the weighted-NSTI table, and the ASV
abundance table all describe the same set of samples before joining them
— silent mismatches here would otherwise produce misleading group
assignments downstream.

``` r
metadata_sample_ids   <- sample_metadata$SampleID
weighted_sample_ids    <- weighted_nsti$SampleID
abundance_sample_ids   <- colnames(asv_count_matrix)

# Check both directions: a PICRUSt2 sample missing from metadata is a hard
# error (every result needs a Condition assignment), while a metadata sample
# missing from the PICRUSt2 output is only a soft warning (that row is
# simply unused below, not a sign anything is broken).
samples_missing_from_metadata <- setdiff(weighted_sample_ids, metadata_sample_ids)
samples_missing_from_weighted  <- setdiff(metadata_sample_ids, weighted_sample_ids)
samples_missing_from_abundance <- setdiff(weighted_sample_ids, abundance_sample_ids)

if (length(samples_missing_from_metadata) > 0) {
  stop(
    "Found ", length(samples_missing_from_metadata), " sample(s) in the PICRUSt2 output ",
    "with no matching row in data/metadata.tsv:\n",
    paste("  -", samples_missing_from_metadata, collapse = "\n")
  )
}

if (length(samples_missing_from_abundance) > 0) {
  stop(
    "Found ", length(samples_missing_from_abundance), " sample(s) in the PICRUSt2 output ",
    "with no matching column in the Step 1 BIOM table:\n",
    paste("  -", samples_missing_from_abundance, collapse = "\n"),
    "\nThis should not normally happen if Steps 1-2 were run on the same input; ",
    "check that both were run from the same data/asv_count_table.csv."
  )
}

if (length(samples_missing_from_weighted) > 0) {
  warning(
    "Found ", length(samples_missing_from_weighted), " sample(s) in data/metadata.tsv ",
    "with no matching entry in the PICRUSt2 output (these will simply be excluded ",
    "from this notebook's plots/tables): ",
    paste(samples_missing_from_weighted, collapse = ", ")
  )
}

cat("Sample ID consistency check passed for", length(weighted_sample_ids), "samples.\n")
```

------------------------------------------------------------------------

# Cross-Reference Step 2’s NSTI Filtering Log

PICRUSt2 reports, in its own console output, exactly how many ASVs its
`max_nsti` filter removed before this notebook ever sees the data
(e.g. *“14 of 5178 ASVs were above the max NSTI cut-off of 2.0 and were
removed”*). Parsing that message from [Step 2](2_picrust2_pipeline.md)’s
saved log provides an independent cross-check against this notebook’s
own count of ASVs above `max_nsti_threshold` computed below.

``` r
# Regex matches PICRUSt2's own two possible log message formats:
#   "<N> of <M> ASVs were above the max NSTI cut-off of <X> and were removed..."
#   "All ASVs were below the max NSTI cut-off of <X> and so all were retained..."
nsti_filter_log_message <- NA_character_
nsti_filter_removed_count <- NA_integer_
nsti_filter_total_count <- NA_integer_

if (file.exists(run_log_path)) {
  # Read the saved combined stdout/stderr from Step 2. The log is supporting
  # provenance only; downstream calculations still use the actual NSTI table.
  log_lines <- readLines(run_log_path, warn = FALSE)

  removed_line <- grep("ASVs were above the max NSTI cut-off", log_lines, value = TRUE)
  retained_line <- grep("All ASVs were below the max NSTI cut-off", log_lines, value = TRUE)

  if (length(removed_line) > 0) {
    # Extract the two leading integers from PICRUSt2's removal message:
    # number removed and total ASVs evaluated.
    nsti_filter_log_message <- removed_line[[1]]
    parsed_numbers <- regmatches(
      nsti_filter_log_message,
      regexec("^(\\d+) of (\\d+) ASVs", nsti_filter_log_message)
    )[[1]]
    if (length(parsed_numbers) == 3) {
      nsti_filter_removed_count <- as.integer(parsed_numbers[[2]])
      nsti_filter_total_count   <- as.integer(parsed_numbers[[3]])
    }
  } else if (length(retained_line) > 0) {
    # The all-retained message states the removed count directly (zero) but
    # does not always expose the total in a stable, parseable form.
    nsti_filter_log_message <- retained_line[[1]]
    nsti_filter_removed_count <- 0L
  } else {
    message(
      "Could not find PICRUSt2's own NSTI filtering message in ", run_log_path, ". ",
      "This can happen with older/newer PICRUSt2 versions that word this message ",
      "differently — this is only a cross-check and is not required for the rest ",
      "of this notebook to run."
    )
  }
} else {
  message(
    "Step 2 log file not found at ", run_log_path, " — skipping the log cross-check. ",
    "This is only a supplementary validation step."
  )
}

cat(
  "Step 2 NSTI Filtering Log Cross-Check:\n",
  "- Log message found:", if (is.na(nsti_filter_log_message)) "(none)" else nsti_filter_log_message, "\n",
  "- ASVs removed (per Step 2's own log):", nsti_filter_removed_count, "\n"
)
```

------------------------------------------------------------------------

# Compute QC Metrics

## Per-ASV NSTI Table

Combines the per-ASV NSTI values with each ASV’s mean relative abundance
across samples, and flags ASVs above `max_nsti_threshold` (matching
[Step 2](2_picrust2_pipeline.md)’s own `max_nsti` filter).

``` r
# Mean relative abundance of each ASV across all samples.
# na.rm = TRUE: the only source of NaN in asv_relative_abundance_matrix is a
# zero-total-count sample column (identified and warned about in Data
# Import above), where every ASV's relative abundance is undefined (0/0).
# Excluding it here means each ASV's mean is taken over the samples where
# it was actually measurable, instead of every ASV's mean silently becoming
# NaN (and every Excel cell showing "#NUM!") because one sample column was
# entirely NaN.
asv_mean_relative_abundance <- data.frame(
  ASV_ID = rownames(asv_relative_abundance_matrix),
  mean_relative_abundance = rowMeans(asv_relative_abundance_matrix, na.rm = TRUE),
  row.names = NULL
)

# Join NSTI values with abundance; an inner join is used deliberately so
# that only ASVs present in both the NSTI table and the Step 1 abundance
# table are retained (mismatches would already have surfaced upstream if
# this notebook's earlier consistency checks matter for ASVs specifically —
# note those checks were at the sample level; a few ASV-level mismatches can
# still legitimately occur here, e.g. singleton ASVs PICRUSt2 could not place).
per_asv_nsti_table <- per_asv_nsti %>%
  inner_join(asv_mean_relative_abundance, by = "ASV_ID") %>%
  mutate(
    exceeds_max_nsti_threshold = NSTI > max_nsti_threshold
  ) %>%
  arrange(desc(NSTI))

n_asvs_unmatched <- nrow(per_asv_nsti) - nrow(per_asv_nsti_table)
if (n_asvs_unmatched > 0) {
  warning(
    n_asvs_unmatched, " ASV(s) with an NSTI value had no matching entry in the ",
    "Step 1 abundance table and were excluded from the per-ASV table/plots."
  )
}

# Cross-check this notebook's own count of flagged ASVs against Step 2's log
n_flagged_here <- sum(per_asv_nsti_table$exceeds_max_nsti_threshold)
cat(
  "Per-ASV NSTI Table:\n",
  "- ASVs summarized:", nrow(per_asv_nsti_table), "\n",
  "- ASVs above max_nsti_threshold (", max_nsti_threshold, "):", n_flagged_here, "\n"
)
if (!is.na(nsti_filter_removed_count) && n_flagged_here != nsti_filter_removed_count) {
  warning(
    "This notebook counted ", n_flagged_here, " ASV(s) above max_nsti_threshold, but ",
    "Step 2's own log reported ", nsti_filter_removed_count, " removed. A mismatch can ",
    "happen if max_nsti_threshold here differs from the --max_nsti value Step 2 was ",
    "actually run with — confirm both use the same cut-off if this matters to you."
  )
}
```

## Per-Sample NSTI Table

Combines the official per-sample weighted NSTI with metadata condition
assignment, the descriptive Tukey/IQR upper-fence flag, and the
complementary “fraction of reads excluded by the `max_nsti` filter”
metric described in [Purpose](#purpose).

``` r
# For each sample, compute what fraction of its total (relative) abundance
# belonged to ASVs above max_nsti_threshold — i.e. ASVs Step 2 excluded
# before ever computing the official weighted_NSTI value.
asvs_above_threshold <- per_asv_nsti_table$ASV_ID[per_asv_nsti_table$exceeds_max_nsti_threshold]

pct_reads_excluded_by_nsti_filter <- vapply(
  colnames(asv_relative_abundance_matrix),
  function(sample_id) {
    # A zero-total-count sample (flagged in Data Import above) has an
    # entirely NaN relative-abundance column. Report NA explicitly here
    # rather than letting it resolve to a misleading 0% (which only happens
    # to be harmless today because no ASV in this run actually exceeds
    # max_nsti_threshold) or propagate as NaN/"#NUM!" in a future run where
    # some ASVs do.
    if (sample_id %in% zero_count_samples) {
      return(NA_real_)
    }
    asvs_in_sample <- rownames(asv_relative_abundance_matrix)
    excluded_mask  <- asvs_in_sample %in% asvs_above_threshold
    100 * sum(asv_relative_abundance_matrix[excluded_mask, sample_id])
  },
  FUN.VALUE = numeric(1)
)

reads_excluded_df <- data.frame(
  SampleID = names(pct_reads_excluded_by_nsti_filter),
  pct_reads_excluded_by_nsti_filter = unname(pct_reads_excluded_by_nsti_filter),
  row.names = NULL
)

# Descriptive Tukey/IQR upper-fence rule for per-sample weighted NSTI —
# see the Purpose section above for why a fixed universal cut-off is not used.
weighted_nsti_quartiles <- quantile(weighted_nsti$weighted_NSTI, probs = c(0.25, 0.75), na.rm = TRUE)
weighted_nsti_iqr <- weighted_nsti_quartiles[[2]] - weighted_nsti_quartiles[[1]]
weighted_nsti_iqr_fence <- weighted_nsti_quartiles[[2]] + iqr_fence_multiplier * weighted_nsti_iqr

# Assemble the full per-sample table
metadata_columns_to_join <- c("SampleID", "Condition", "Reference")

per_sample_nsti_table <- weighted_nsti %>%
  inner_join(sample_metadata[, metadata_columns_to_join], by = "SampleID") %>%
  inner_join(reads_excluded_df, by = "SampleID") %>%
  mutate(
    exceeds_weighted_nsti_iqr_fence = weighted_NSTI > weighted_nsti_iqr_fence
  ) %>%
  arrange(desc(weighted_NSTI))

cat(
  "Per-Sample NSTI Table:\n",
  "- Samples summarized:", nrow(per_sample_nsti_table), "\n",
  "- Tukey/IQR upper fence (Q3 + ", iqr_fence_multiplier, " x IQR):",
  round(weighted_nsti_iqr_fence, 4), "\n",
  "- Samples exceeding the fence:", sum(per_sample_nsti_table$exceeds_weighted_nsti_iqr_fence), "\n"
)
```

------------------------------------------------------------------------

# Visualizations

## Boxplot: Weighted NSTI by Condition

Shows whether prediction quality (per-sample weighted NSTI) differs
visibly between metadata conditions. A condition with a higher
distribution may reflect a less well-characterized community in the
reference database or a technical/batch issue worth investigating. This
plot is descriptive; this notebook does not attach a between-condition
p-value to it.

``` r
# Build the static ggplot first: boxes summarize each condition, while jittered
# points preserve the individual samples that generated those summaries.
nsti_boxplot_static <- ggplot(
  per_sample_nsti_table,
  aes(x = Condition, y = weighted_NSTI, fill = Condition, text = SampleID)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  # The dashed line is the study-specific descriptive upper fence calculated
  # earlier; it is not a universal quality threshold or significance cutoff.
  geom_hline(
    yintercept = weighted_nsti_iqr_fence,
    linetype = "dashed", color = "#e74c3c"
  ) +
  labs(
    title = "Per-Sample Weighted NSTI by Condition",
    subtitle = paste0("Dashed line: descriptive Tukey/IQR upper fence (Q3 + ", iqr_fence_multiplier, " x IQR)"),
    x = "Condition",
    y = "Weighted NSTI (abundance-weighted mean, PICRUSt2 output)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Convert to Plotly only after the statistical layers and labels are fixed.
# Hover text exposes SampleID without crowding the static display.
nsti_boxplot_interactive <- ggplotly(nsti_boxplot_static, tooltip = c("x", "y", "text"))
nsti_boxplot_interactive
```

The next chunk saves the displayed widget as a self-contained HTML file
and adds a portable report link. Saving is separated from plot
construction so users running chunks interactively can inspect the plot
before writing it.

## Histogram/Density: Per-ASV NSTI Distribution

Shows the overall distribution of per-ASV NSTI values across the whole
dataset, with a reference line at `max_nsti_threshold` — the fraction of
the distribution to the right of that line is exactly what [Step
2](2_picrust2_pipeline.md)’s `max_nsti` filter already excluded.

``` r
# Overlay a normalized histogram and kernel-density estimate so both binned
# counts and the overall distribution shape are visible on the same scale.
nsti_histogram_static <- ggplot(per_asv_nsti_table, aes(x = NSTI)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#2c3e50", alpha = 0.75) +
  geom_density(color = "#f39c12", linewidth = 1) +
  # Mark the exact per-ASV threshold passed to PICRUSt2 in Step 2.
  geom_vline(
    xintercept = max_nsti_threshold,
    linetype = "dashed", color = "#e74c3c"
  ) +
  labs(
    title = "Distribution of Per-ASV NSTI Values",
    subtitle = paste0("Dashed line: max_nsti_threshold = ", max_nsti_threshold),
    x = "NSTI (nearest sequenced taxon index)",
    y = "Density"
  ) +
  theme_minimal(base_size = 13)

# The interactive form provides exact NSTI values on hover while retaining the
# same bins, density curve, and threshold as the static plot.
nsti_histogram_interactive <- ggplotly(nsti_histogram_static)
nsti_histogram_interactive
```

The following reporting chunk saves and links the interactive
distribution. The underlying per-ASV values remain available at full
precision in the Excel workbook.

## Scatterplot: NSTI vs. Abundance

Identifies ASVs that are simultaneously **abundant** and **poorly
supported** (high NSTI) — the ones of greatest concern, since they
represent a meaningful share of the community whose predicted function
rests on a distant reference genome. Points are colored by whether they
exceed `max_nsti_threshold`.

``` r
# Map each ASV's mean relative abundance to the x-axis and its placement
# distance to the y-axis. Color separates ASVs retained versus excluded by the
# configured per-ASV NSTI threshold.
nsti_scatter_static <- ggplot(
  per_asv_nsti_table,
  aes(
    x = mean_relative_abundance, y = NSTI,
    color = exceeds_max_nsti_threshold, text = ASV_ID
  )
) +
  geom_point(alpha = 0.7, size = 1.8) +
  # Reuse the same threshold line and colors as the distribution plot so the
  # visual quality language is consistent across the report.
  geom_hline(
    yintercept = max_nsti_threshold,
    linetype = "dashed", color = "#e74c3c"
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  scale_color_manual(
    values = c(`FALSE` = "#2c3e50", `TRUE` = "#e74c3c"),
    name = paste0("Exceeds max_nsti_threshold (", max_nsti_threshold, ")")
  ) +
  labs(
    title = "Per-ASV NSTI vs. Mean Relative Abundance",
    subtitle = "Top-right quadrant: abundant ASVs with unreliable functional predictions",
    x = "Mean relative abundance across samples",
    y = "NSTI (nearest sequenced taxon index)"
  ) +
  theme_minimal(base_size = 13)

# Retain ASV_ID in hover text so individual points can be traced back to the
# Per_ASV_NSTI workbook sheet without labeling every point on the canvas.
nsti_scatter_interactive <- ggplotly(nsti_scatter_static, tooltip = c("text", "x", "y"))
nsti_scatter_interactive
```

This final visualization chunk saves the scatterplot and adds its report
link, using the same self-contained format as the other two NSTI plots.

------------------------------------------------------------------------

# Results Summary

## Compute Filtering Log Summary Table

This chunk assembles a one-row summary table combining [Step
2](2_picrust2_pipeline.md)’s own log-reported filtering counts with this
notebook’s independently recomputed count, so the cross-check performed
above has a durable, exportable record.

``` r
nsti_filtering_log_summary <- data.frame(
  `Log File` = if (file.exists(run_log_path)) basename(run_log_path) else NA_character_,
  `Log Message` = nsti_filter_log_message,
  `ASVs Removed (per Step 2 log)` = nsti_filter_removed_count,
  `Total ASVs (per Step 2 log)` = nsti_filter_total_count,
  `ASVs Above Threshold (recomputed here)` = n_flagged_here,
  `max_nsti_threshold Used Here` = max_nsti_threshold,
  check.names = FALSE
)

# Display the summary using DT::datatable() (already loaded above), the
# same interactive table style used elsewhere in this notebook. The one-row
# summary is reshaped into a tidy two-column Field/Value table first, since
# DT::datatable() (unlike kable()) is not designed to transpose a
# single-row data frame directly.
nsti_filtering_log_summary_long <- data.frame(
  Field = names(nsti_filtering_log_summary),
  Value = unlist(nsti_filtering_log_summary, use.names = FALSE),
  check.names = FALSE
)
datatable(
  nsti_filtering_log_summary_long,
  options = list(dom = 't', pageLength = nrow(nsti_filtering_log_summary_long),
                 scrollY = "400px", scrollCollapse = TRUE, paging = FALSE,
                 ordering = FALSE, searching = FALSE),
  rownames = FALSE,
  caption = "NSTI filtering log summary"
)
```

## Display Per-Sample and Per-ASV Tables

This chunk presents the two principal quality-control tables as
horizontally scrollable views without sorting controls or search fields:
one row per sample for weighted NSTI and excluded-read burden, and one
row per ASV for placement distance and abundance context.

``` r
# per_sample_nsti_table has one row per sample (potentially hundreds) --
# display it as an interactive, paginated/scrollable table rather than a
# single static block. DT::formatRound() rounds only the on-screen display,
# leaving the underlying per_sample_nsti_table object (and its later Excel
# export) at full precision.
numeric_columns_per_sample <- names(per_sample_nsti_table)[sapply(per_sample_nsti_table, is.numeric)]
numeric_columns_per_sample <- setdiff(
  numeric_columns_per_sample,
  "pct_reads_excluded_by_nsti_filter"
)
datatable(
  per_sample_nsti_table,
  options = list(dom = 't', pageLength = 10, scrollX = TRUE,
                 scrollY = "400px", scrollCollapse = TRUE, paging = FALSE,
                 ordering = FALSE, searching = FALSE),
  rownames = FALSE,
  caption = "Per-Sample Weighted NSTI"
) %>%
  formatRound(columns = numeric_columns_per_sample, digits = 4) %>%
  formatRound(columns = "pct_reads_excluded_by_nsti_filter", digits = 1)

# per_asv_nsti_table is already bounded to the 20 highest values; displayed
# as an interactive DT::datatable() for consistency with the rest of this
# notebook, with formatRound() rounding only the on-screen numeric display.
per_asv_nsti_table_top20 <- head(per_asv_nsti_table, 20)
numeric_columns_per_asv <- which(sapply(per_asv_nsti_table_top20, is.numeric))
datatable(
  per_asv_nsti_table_top20,
  options = list(dom = 't', pageLength = 20, scrollX = TRUE,
                 scrollY = "400px", scrollCollapse = TRUE, paging = FALSE,
                 ordering = FALSE, searching = FALSE),
  rownames = FALSE,
  caption = "Per-ASV NSTI (20 highest values shown)"
) %>%
  formatRound(columns = numeric_columns_per_asv, digits = 4)
```

## Export Summary to Excel

This chunk writes the per-sample table, per-ASV table, and filtering log
summary computed above to the output workbook using the project’s shared
Excel-writing helper,
[add_sheet_to_excel_function.R](../functions/add_sheet_to_excel_function.R).

## Document and Export Column Dictionary

Every column written to this workbook is documented in a trailing
`Column_Dictionary` sheet (Sheet / Column / Explanation), so the file is
self-explanatory to anyone opening it outside of this notebook.

``` r
# Human-readable explanation for every column written to the Per_Sample_NSTI
# sheet, keyed by exact column name.
per_sample_nsti_column_descriptions <- c(
  `SampleID` =
    "Sample identifier, matching the column names in Step 1's ASV count table and the SampleID column in data/metadata.tsv.",
  `weighted_NSTI` =
    "Per-sample abundance-weighted mean Nearest Sequenced Taxon Index, taken directly from PICRUSt2's own weighted_nsti.tsv output. Computed by PICRUSt2 only over ASVs that survived the max_nsti filter (see the Purpose section above).",
  `Condition` =
    "Sample condition from data/metadata.tsv, used to compare NSTI quality descriptively across experimental conditions.",
  `Reference` =
    "TRUE when the sample belongs to the reference/baseline Condition; FALSE otherwise.",
  `pct_reads_excluded_by_nsti_filter` =
    "Percentage of the sample's total relative abundance that belonged to ASVs excluded by Step 2's max_nsti filter. NA for samples with zero total ASV counts, where relative abundance is undefined (0/0).",
  `exceeds_weighted_nsti_iqr_fence` =
    "TRUE if weighted_NSTI exceeds the descriptive Tukey/IQR upper fence (Q3 + iqr_fence_multiplier x IQR) computed across all samples in this run; this is a screening flag, not a hypothesis test or p-value."
)

# Human-readable explanation for every column written to the Per_ASV_NSTI
# sheet.
per_asv_nsti_column_descriptions <- c(
  `ASV_ID` =
    "ASV identifier (PICRUSt2's own \"sequence\" column in combined_marker_predicted_and_nsti.tsv), matching the IDs used throughout Steps 1-2.",
  `NSTI` =
    "Nearest Sequenced Taxon Index for this ASV, taken directly from PICRUSt2's combined_marker_predicted_and_nsti.tsv output (its metadata_NSTI column). Lower values indicate a closer reference genome and a higher-confidence functional prediction.",
  `mean_relative_abundance` =
    "Mean relative abundance of this ASV across all samples (0-1 scale), excluding any zero-total-count sample column (rowMeans(..., na.rm = TRUE)) so one empty sample cannot corrupt every ASV's mean.",
  `exceeds_max_nsti_threshold` =
    "TRUE if NSTI is greater than max_nsti_threshold (matches PICRUSt2's own --max_nsti default of 2.0) -- these are exactly the ASVs Step 2 already excluded from the KO/EC metagenome predictions."
)

# Human-readable explanation for every column written to the
# NSTI_Filtering_Log_Summary sheet.
nsti_filtering_log_summary_column_descriptions <- c(
  `Log File` =
    "File name of the Step 2 console log parsed for the NSTI filtering cross-check message.",
  `Log Message` =
    "The exact PICRUSt2 log line reporting how many ASVs were excluded (or that all were retained) by the max_nsti filter, as parsed from Step 2's log.",
  `ASVs Removed (per Step 2 log)` =
    "Number of ASVs PICRUSt2 itself reported removing for exceeding max_nsti, parsed directly from Step 2's log.",
  `Total ASVs (per Step 2 log)` =
    "Total ASV count PICRUSt2 reported in the same log line, when the \"N of M ASVs\" removal message format was found (NA if only the \"all ASVs retained\" message format was found).",
  `ASVs Above Threshold (recomputed here)` =
    "This notebook's own independent count of ASVs with NSTI > max_nsti_threshold, computed directly from combined_marker_predicted_and_nsti.tsv, for cross-checking against Step 2's log.",
  `max_nsti_threshold Used Here` =
    "The max_nsti_threshold value configured in this notebook's Configuration section, used for the independent recomputation above."
)

# Build the Sheet / Column / Explanation rows for each sheet written to this
# workbook and combine them into one table.
column_dictionary <- rbind(
  build_column_dictionary(
    sheet_name    = "Per_Sample_NSTI",
    data          = per_sample_nsti_table,
    descriptions  = per_sample_nsti_column_descriptions,
    workbook_path = output_excel_path,
    rownames      = FALSE
  ),
  build_column_dictionary(
    sheet_name    = "Per_ASV_NSTI",
    data          = per_asv_nsti_table,
    descriptions  = per_asv_nsti_column_descriptions,
    workbook_path = output_excel_path,
    rownames      = FALSE
  ),
  build_column_dictionary(
    sheet_name    = "NSTI_Filtering_Log_Summary",
    data          = nsti_filtering_log_summary,
    descriptions  = nsti_filtering_log_summary_column_descriptions,
    workbook_path = output_excel_path,
    rownames      = FALSE
  )
)

# The dictionary is intentionally not printed in the report. It is exported
# only as the trailing Column_Dictionary worksheet below.
```

This chunk writes the displayed dictionary to the workbook as its final
sheet and links the completed NSTI quality-assessment workbook.

------------------------------------------------------------------------

# Output File Summary

The tree below lists every file this notebook has written to its own
output folder,
[results/3_nsti_quality_assessment/](../../results/3_nsti_quality_assessment/),
as a clickable, portable link (relative to this notebook’s own location)
with a short description – built live from what is actually on disk at
knit time via the project’s shared
[render_output_tree_function.R](../functions/render_output_tree_function.R)
helper, so it always matches this run’s real output rather than a
hand-maintained list.

------------------------------------------------------------------------

# Output Interpretation

- **`Per_Sample_NSTI` sheet**: one row per sample, with the official
  PICRUSt2 `weighted_NSTI` value, its metadata `Condition` and
  `Reference` values, the `pct_reads_excluded_by_nsti_filter` metric,
  and an `exceeds_weighted_nsti_iqr_fence` screening flag. Samples with
  a high `pct_reads_excluded_by_nsti_filter` warrant a closer look even
  if their `weighted_NSTI` itself looks fine — remember `weighted_NSTI`
  is computed only over the ASVs that survived filtering.
- **`Per_ASV_NSTI` sheet**: one row per ASV, sorted by NSTI descending,
  with `mean_relative_abundance` and an `exceeds_max_nsti_threshold`
  flag. ASVs flagged `TRUE` are exactly the ones [Step
  2](2_picrust2_pipeline.md) already excluded from the KO/EC metagenome
  predictions. If any sample has zero total ASV counts (see the warning
  printed during [Data Import](#import-asv-abundance)), that sample is
  excluded from this mean (`na.rm = TRUE`) rather than being allowed to
  silently turn every ASV’s mean into `NaN`/`#NUM!` in Excel.
- **`NSTI_Filtering_Log_Summary` sheet**: a single-row cross-check
  comparing this notebook’s own count of flagged ASVs against the count
  PICRUSt2 itself reported in [Step 2](2_picrust2_pipeline.md)’s log.
  These should match; if they don’t, `max_nsti_threshold` in this
  notebook’s [Configuration](#adjust-run-parameters) likely differs from
  the `--max_nsti` value [Step 2](2_picrust2_pipeline.md) was actually
  run with.
- **`Column_Dictionary` sheet**: a trailing Sheet / Column / Explanation
  table documenting every column in the three sheets above, generated
  automatically from the actual exported data (see [Document and Export
  Column Dictionary](#column-dictionary)) so it cannot drift out of sync
  with what the workbook actually contains.
- **What counts as “good”** depends heavily on the environment being
  studied — well-characterized communities can show lower NSTI values
  overall than under-sequenced environments. Use the boxplot and
  scatterplot to look for *within-study* outliers and systematic
  differences between your own conditions, rather than comparing against
  an absolute external benchmark.
- If a substantial fraction of ASVs or reads are flagged, treat
  downstream differential-abundance results from [Step
  4](4_functional_differential_abundance_analysis.md) and taxon
  contribution results from [Step 5](5_taxon_contribution.md) with
  appropriate caution, and consider reporting these QC metrics alongside
  any biological conclusions drawn from this workflow.

------------------------------------------------------------------------

# Recommended Next Step

This notebook does not export a file that the later analysis notebooks
read as an input — both use [Step 2](2_picrust2_pipeline.md)’s
[picrust2_out_pipeline/](../../results/2_picrust2_pipeline/picrust2_out_pipeline/)
output directly. Its value is diagnostic: use the
`nsti_quality_summary.xlsx` workbook and the three plots above to check
for samples with a high `pct_reads_excluded_by_nsti_filter` or a
`weighted_NSTI` exceeding the descriptive Tukey/IQR upper fence, and for
abundant, high-NSTI ASVs in the scatterplot, and decide whether any of
those warrant exclusion or extra caution before interpreting results
downstream. Once you are satisfied with data quality, proceed to [Step 4
— Functional Differential Abundance
Analysis](4_functional_differential_abundance_analysis.md).

------------------------------------------------------------------------

# Session Information

Record the R environment for reproducibility.

------------------------------------------------------------------------

# References

## NSTI and PICRUSt2

- Langille MGI, et al. (2013). Predictive functional profiling of
  microbial communities using 16S rRNA marker gene sequences. *Nat
  Biotechnol* 31, 814-821. (original NSTI concept, PICRUSt1)
  <https://doi.org/10.1038/nbt.2676>
- Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for
  prediction of metagenome functions. *Nat Biotechnol* 38, 685-688.
  <https://doi.org/10.1038/s41587-020-0548-6>
- [PICRUSt2 Metagenome Prediction wiki
  page](https://github.com/picrust/picrust2/wiki/Metagenome-prediction)
  (describes `weighted_nsti.tsv.gz` and the `--max_nsti` filter)
- [PICRUSt2 Full Pipeline Script wiki
  page](https://github.com/picrust/picrust2/wiki/Full-pipeline-script)
  (`--max_nsti` default and other `picrust2_pipeline.py` options)
- [PICRUSt2 Key
  Limitations](https://github.com/picrust/picrust2/wiki/Key-Limitations)

## R Packages

- Sievert C (2020). *Interactive Web-Based Data Visualization with R,
  plotly, and shiny*. <https://plotly-r.com/> (`plotly`)
- Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
  Springer-Verlag New York. (`ggplot2`)
- McDonald D, Clemente JC, Kuczynski J, et al. (2012). The Biological
  Observation Matrix (BIOM) format or: how I learned to stop worrying
  and love the ome-ome. *GigaScience* 1, 7.
  <https://doi.org/10.1186/2047-217X-1-7> (`biomformat`, BIOM format)

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### [data/metadata.tsv](../../data/metadata.tsv) Not Found

**Error**:
`Could not find the sample metadata file: .../data/metadata.tsv`

**Cause**: This file is a new requirement starting at [Step
3](3_nsti_quality_assessment.md) — it was not needed for [Step
1](1_prepare_picrust2_inputs.md) or [Step 2](2_picrust2_pipeline.md).

**Actions**:

- Create [data/metadata.tsv](../../data/metadata.tsv) following the
  three-column schema documented in
  [data/README.md](../../data/README.md): `SampleID`, `Condition`, and
  `Reference`.
- Confirm the file is tab-separated (not comma-separated) and has a
  header row.

### Sample ID Mismatches

**Error**:
`Found N sample(s) in the PICRUSt2 output with no matching row in [data/metadata.tsv](../../data/metadata.tsv)`

**Cause**: `SampleID` values in
[data/metadata.tsv](../../data/metadata.tsv) do not exactly match the
sample identifiers used as column names in the Step 1 BIOM table.

**Actions**:

- Compare `sample_metadata$SampleID` against
  `colnames(asv_count_matrix)` directly in an R console to spot the
  exact mismatch (case, whitespace, or naming convention differences are
  the usual culprits).
- Re-export `data/metadata.tsv` using the exact sample identifiers from
  [data/asv_count_table.csv](../../data/asv_count_table.csv).

### [combined_marker_predicted_and_nsti.tsv](../../results/2_picrust2_pipeline/picrust2_out_pipeline/combined_marker_predicted_and_nsti.tsv) or `weighted_nsti.tsv` Not Found

**Error**: `Could not find the expected file: ...`

**Possible causes**:

- [Step 2](2_picrust2_pipeline.md) was run with
  `decompress_output_files <- FALSE`, so the file exists only as
  `weighted_nsti.tsv.gz` (this notebook does not look for `.gz` files,
  since the default configuration always decompresses Step 2’s output
  first).
- [Step 2](2_picrust2_pipeline.md) has not been run yet, or was run with
  `--no_pathways`/older PICRUSt2 flags that skip NSTI computation
  entirely (uncommon; NSTI is included by default).

**Actions**:

- Re-run Step 2 with `decompress_output_files <- TRUE` (the default), or
  manually decompress the specific `.gz` file this error names.
- Confirm
  [results/2_picrust2_pipeline/picrust2_out_pipeline/](../../results/2_picrust2_pipeline/picrust2_out_pipeline/)
  exists and contains `KO_metagenome_out/`.

### Log Cross-Check Count Mismatch

**Symptom**: A warning that this notebook’s own count of flagged ASVs
does not match the count reported in [Step 2](2_picrust2_pipeline.md)’s
log.

**Cause**: `max_nsti_threshold` in this notebook differs from the
`--max_nsti` value used for the completed Step 2 run. Steps 2, 3, and 5
read the same `PICRUST2_MAX_NSTI` environment variable, but an existing
result folder may have been produced earlier with another value.

**Actions**:

- Set `PICRUST2_MAX_NSTI` to the value used for the completed Step 2 run
  before rendering or executing this notebook, or regenerate Step 2 with
  the intended shared value.

### Suspiciously High NSTI Values Across the Entire Dataset

**Symptom**: Nearly all ASVs/samples are flagged, even ones that should
be well-characterized (e.g. common human gut taxa).

**Possible causes**:

- The input sequences are not actually 16S rRNA sequences, or are from a
  very different marker gene/region than PICRUSt2’s default reference
  tree expects.
- A sequence-orientation or primer-trimming issue in the companion DADA2
  workflow that alters placement accuracy.

**Actions**:

- Spot-check a few high-NSTI, high-abundance sequences (see the
  [scatterplot](#plot-scatter)) with a BLAST search against a general
  reference database to confirm they are the expected marker gene and
  organism group.
- Revisit the amplicon primer/region used in the original sequencing run
  against PICRUSt2’s expectations (see the [PICRUSt2 Key
  Limitations](https://github.com/picrust/picrust2/wiki/Key-Limitations)
  reference).

### `mean_relative_abundance` Shows `#NUM!` for (Nearly) Every ASV in Excel

**Symptom**: In the exported `Per_ASV_NSTI` sheet, the
`mean_relative_abundance` column shows `#NUM!` in almost every row
rather than a percentage.

**Cause**: One (or more) samples in [Step
1](1_prepare_picrust2_inputs.md)’s
[study_seqs.biom](../../results/1_prepare_picrust2_inputs/study_seqs.biom)
has zero total ASV counts (commonly a negative/blank control that
yielded no amplicon reads). Converting that column to relative abundance
divides by a column sum of `0`, producing `0/0 = NaN` for every ASV in
that one sample. `rowMeans()` returns `NaN` for a row if *any* value in
that row is `NaN` (its default is `na.rm = FALSE`), so a single
all-`NaN` sample column silently turns essentially every ASV’s
`mean_relative_abundance` into `NaN` — and `openxlsx` writes R’s `NaN`
to a `.xlsx` cell as the Excel error `#NUM!`, not as a blank cell, which
is why it appears as a spreadsheet error rather than an obviously
missing value.

**What this notebook does about it**: [Data
Import](#import-asv-abundance) now detects any zero-total-count sample
column immediately after importing the BIOM table and prints a
`warning()` naming it. [Compute QC Metrics](#compute-per-asv) then
computes `mean_relative_abundance` with `rowMeans(..., na.rm = TRUE)`,
so that one sample is excluded from the average instead of corrupting
every ASV’s value, and the per-sample
`pct_reads_excluded_by_nsti_filter` metric reports `NA` (not a
misleading `0%`) for that specific sample.

**Actions**:

- If the printed warning names a sample you recognize as an
  intentionally empty negative/blank control, no action is needed — this
  is expected, and the fix above already handles it correctly.
- If the printed warning names a sample that was supposed to contain
  real biological data, treat this as a Step 1 upstream data problem
  (e.g. a failed sequencing run, a sample dropped during DADA2
  filtering, or a sample ID typo in
  [data/asv_count_table.csv](../../data/asv_count_table.csv)) rather
  than something to fix in this notebook.

### BIOM Table Fails to Import as a Numeric Matrix

**Error**:
`Error in colSums(asv_count_matrix) : 'x' must be numeric or complex`

**Cause**: `as(biom_object, "matrix")` — the S4 coercion method
previously used to extract the count matrix from a `biom`-class object —
is not reliably registered across all `biomformat` versions, and can
silently return a non-numeric object instead of raising an error itself.
The failure then only surfaces later, at `colSums()`.

**What this notebook does about it**: The import chunk now uses
`as.matrix(biom_data(read_biom(...)))` instead — `biom_data()` is
`biomformat`’s own documented accessor for the count matrix, and
wrapping it in `as.matrix()` guarantees a plain, dense, numeric matrix
regardless of whether the file stores data densely or sparsely. A
defensive `is.numeric()` check immediately follows and will now stop
with a clear message (reporting the actual storage mode) if this ever
occurs again, rather than failing several lines later inside
`colSums()`.

**Actions**:

- Update to the current version of this notebook if you see the
  `colSums` error above — this is a code fix, not a data problem.
- If the new `is.numeric()` check itself fails, confirm
  `study_seqs.biom` from Step 1 was written correctly (re-running Step 1
  is usually sufficient) and check your installed `biomformat` version
  (`packageVersion("biomformat")`).
