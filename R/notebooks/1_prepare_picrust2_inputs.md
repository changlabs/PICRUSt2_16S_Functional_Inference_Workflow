Step 1: Prepare PICRUSt2 Input Files from DADA2 Results
================

- [Introduction](#introduction)
  - [Purpose](#purpose)
  - [Prerequisites](#prerequisites)
  - [What This Notebook Does](#what-this-notebook-does)
  - [Statistics and Calculations](#statistics-and-calculations)
  - [Expected Input](#expected-input)
  - [Expected Output](#expected-output)
- [Environment Setup](#environment-setup)
  - [Load Required Packages](#load-packages)
- [Configuration](#configuration)
  - [Define Path Parameters](#define-paths)
- [Data Import](#data-import)
  - [Import ASV Count Table](#import-count-table)
  - [Import ASV Sequence Table](#import-sequence-table)
- [BIOM Table Construction](#biom-table-construction)
  - [Transpose Count Table](#transpose-counts)
  - [Validate ASV/Sequence Consistency](#validate-asv-consistency)
  - [Build and Export the BIOM Table](#build-export-biom)
- [FASTA File Construction](#fasta-file-construction)
  - [Validate Sequences](#validate-sequences)
  - [Construct and Export the FASTA File](#build-export-fasta)
- [Results Summary](#results-summary)
  - [Compute Summary Statistics](#summary-stats)
  - [Export Summary to Excel](#export-summary)
  - [Document and Export Column Dictionary](#column-dictionary)
- [Output File Summary](#output-file-summary)
- [Output Interpretation](#output-interpretation)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [Methods](#methods)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [Input File Not Found](#input-file-not-found)
    - [ASVs Missing a Corresponding
      Sequence](#asvs-missing-a-corresponding-sequence)
    - [Invalid Sequence Characters](#invalid-sequence-characters)
    - [BIOM Table Shape Looks Wrong](#biom-table-shape-looks-wrong)

<!-- The hidden setup chunk below controls whether analysis code is executed. HTML reports run the workflow; GitHub Markdown remains a non-executing tutorial and code reference. -->
<!-- The hidden CSS chunk below affects only the rendered HTML report and is deliberately excluded from GitHub Markdown. -->

# Introduction

## Purpose

This notebook is **Step 1** of the PICRUSt2 functional prediction
workflow. It prepares the two input files
[PICRUSt2](https://github.com/picrust/picrust2/wiki) requires — a
[BIOM](http://biom-format.org/)-format feature (OTU/ASV) table and a
FASTA file of representative marker-gene sequences — by reformatting the
ASV abundance table and representative sequences produced by the DADA2
pipeline.

## Prerequisites

Before running this notebook, ensure that:

1.  Required R packages are installed — run
    [setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
    once per R environment (see the repository [README](../../README.md)
    for details). The [Load Required Packages](#load-packages) chunk
    below will only *use* these packages; it does not install them.
2.  [asv_count_table.csv](../../data/asv_count_table.csv) is present in
    the [data/](../../data/) directory — a sample x ASV read-count
    matrix. This is a direct output of the DADA2 pipeline.
3.  [asv_sequences.csv](../../data/asv_sequences.csv) is present in the
    [data/](../../data/) directory — an ASV ID-to-sequence mapping. This
    corresponds to [asv_sequences.csv](../../data/asv_sequences.csv)
    from the same DADA2 pipeline step (copy or rename it into `data/` as
    [asv_sequences.csv](../../data/asv_sequences.csv), or adjust the
    input filename in the [Configuration](#define-paths) section below).

## What This Notebook Does

The workflow accomplishes the following tasks:

1.  **Count Table Import**: Reads the sample x ASV abundance matrix
    exported by DADA2.
2.  **Matrix Transposition**: Transposes the table to the feature (ASV)
    x sample orientation required by the [BIOM
    specification](http://biom-format.org/).
3.  **ASV/Sequence Consistency Check**: Verifies that every ASV in the
    count table has a corresponding sequence, and flags any mismatches
    before writing output.
4.  **BIOM Table Construction**: Converts the transposed count matrix
    into a [BIOM](http://biom-format.org/) object and writes it to disk
    in BIOM 1.0 (JSON) format.
5.  **FASTA Construction**: Validates each ASV sequence as a proper DNA
    string using
    [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
    and writes a standard-compliant FASTA file of representative
    sequences.
6.  **Output Read-Back Validation**: Reopens both files after writing
    and verifies their identifiers and values against the in-memory
    source objects, so a truncated or malformed export cannot be passed
    silently to Step 2.
7.  **Report Generation**: Exports a summary table (ASV counts, sample
    counts, sequence length distribution, output file sizes) to an Excel
    file for documentation and review.

## Statistics and Calculations

This preparation step performs **no hypothesis test**. Its reported
statistics are descriptive validation summaries: numbers of samples and
ASVs, total reads, minimum/maximum/mean sequence length, and output file
sizes. Transposition changes only table orientation; it does not
normalize, transform, filter, or otherwise alter the read counts.
Sequence and identifier checks are deterministic data- integrity checks
rather than statistical tests.

## Expected Input

- **Location**: Both input files should be placed in the
  [data/](../../data/) directory relative to your project root.
- [data/asv_count_table.csv](../../data/asv_count_table.csv) — first
  column contains sample identifiers (read in as row names); remaining
  columns are ASV IDs (e.g. `ASV1`, `ASV2`, …) with integer read counts.
- [data/asv_sequences.csv](../../data/asv_sequences.csv) — two columns:
  `ASV_ID` and `sequence`, with one row per ASV.

## Expected Output

- [results/1_prepare_picrust2_inputs/study_seqs.biom](../../results/1_prepare_picrust2_inputs/study_seqs.biom)
  — BIOM-format feature table (ASV x sample counts), used as the `-i`
  input to `picrust2_pipeline.py`.
- [results/1_prepare_picrust2_inputs/study_seqs.fna](../../results/1_prepare_picrust2_inputs/study_seqs.fna)
  — FASTA file of representative ASV sequences, used as the `-s` input
  to `picrust2_pipeline.py`.
- [results/1_prepare_picrust2_inputs/picrust2_input_summary.xlsx](../../results/1_prepare_picrust2_inputs/picrust2_input_summary.xlsx)
  — Summary statistics for the reformatting step, plus a trailing
  `Column_Dictionary` sheet documenting every column in the workbook.

------------------------------------------------------------------------

# Environment Setup

<div class="alert alert-info">

**Before running this notebook**: install the required R packages by
running
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
once per R environment. This only needs to be done once for this
workflow.

</div>

## Load Required Packages

The chunk below loads every R package this notebook depends on —
[biomformat](https://bioconductor.org/packages/biomformat/) for
reading/writing the BIOM feature table and
[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
for sequence validation and FASTA writing — plus three helper functions
shared across the whole pipeline: Excel writing, column-dictionary
documentation, and clickable output links (each sourced from
[R/functions/](../functions/)).

``` r
# NOTE: This chunk assumes the packages below are already installed. If any
# library() call fails with "there is no package called ...", run
# setup/install_R_dependencies.R once to install every package this
# workflow needs, then re-run this notebook.

# biomformat: Bioconductor package for reading and writing BIOM-format files
# Provides make_biom() and write_biom() used to construct the PICRUSt2 feature table
library(biomformat)

# Biostrings: Bioconductor package for biological sequence handling
# Used here to validate that every ASV sequence is a well-formed DNA string
# and to write a standards-compliant FASTA file with a battle-tested writer
# rather than hand-rolled string concatenation
library(Biostrings)

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# openxlsx: Excel file creation and manipulation
# Allows reading and writing Excel files without requiring Java dependencies
library(openxlsx)

# Source the custom Excel utility function from the project's function library
# This function handles creating or appending sheets to Excel workbooks
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source the custom column-dictionary utility function from the project's
# function library. This function documents every column of a sheet being
# exported, so a trailing "Column_Dictionary" sheet can be added to the
# workbook (see the Document and Export Column Dictionary section below).
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source the custom output-links utility function from the project's function
# library. This function renders a clickable Markdown link for each output
# file written by a later chunk, for reproducibility across machines.
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

## Define Path Parameters

Set up the input and output directory paths. Using `here()` ensures
paths are relative to the project root, making the script portable
across different systems.

``` r
# Base data folder containing the DADA2-derived ASV table and sequence mapping
data_folder <- here("data")

# Base results folder for all pipeline outputs
results_folder <- here("results")

# Output folder for this step's results, numbered for pipeline step tracking
output_folder <- here(results_folder, "1_prepare_picrust2_inputs")

# Create the output directory and any necessary parent directories
# If the directory already exists, this operation has no effect (safe to re-run)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Path to the sample x ASV read-count table exported by the companion DADA2 workflow
input_asv_count_table_path <- here(data_folder, "asv_count_table.csv")

# Path to the input ASV ID-to-sequence mapping table
input_asv_sequences_path <- here(data_folder, "asv_sequences.csv")

# Path to the BIOM-format feature table that will be written for PICRUSt2 (-i argument)
output_biom_path <- here(output_folder, "study_seqs.biom")

# Path to the FASTA file of representative sequences that will be written for PICRUSt2 (-s argument)
output_fasta_path <- here(output_folder, "study_seqs.fna")

# Path to the Excel workbook summarizing this reformatting step
output_excel_path <- here(output_folder, "picrust2_input_summary.xlsx")

# Name of the worksheet within the summary Excel file
output_excel_sheet <- "PICRUSt2_Input_Summary"
```

------------------------------------------------------------------------

# Data Import

## Import ASV Count Table

Read the sample x ASV abundance matrix. The first column (sample
identifiers) is used as row names.

``` r
# Read the count table with the first column (sample IDs) as row names
# header = TRUE: the first row contains ASV ID column names
asv_count_table <- read.csv(
  input_asv_count_table_path,
  header    = TRUE,
  row.names = 1,
  check.names = FALSE  # preserve ASV IDs exactly as written (e.g. "ASV1"), no syntactic renaming
)

# Validate the count-table contract before transposition or BIOM conversion.
# PICRUSt2 requires a non-empty, finite, non-negative integer count matrix;
# coercing a malformed character column later would otherwise produce a much
# less informative make_biom() error.
if (nrow(asv_count_table) == 0L || ncol(asv_count_table) == 0L) {
  stop("data/asv_count_table.csv must contain at least one sample and one ASV.")
}
if (anyDuplicated(colnames(asv_count_table))) {
  stop("data/asv_count_table.csv contains duplicate ASV column names.")
}
count_values <- as.matrix(asv_count_table)
if (!is.numeric(count_values)) {
  stop("Every ASV abundance column in data/asv_count_table.csv must be numeric.")
}
if (any(!is.finite(count_values)) || any(count_values < 0)) {
  stop("ASV counts must be finite and non-negative; NA/NaN/Inf/negative values are not allowed.")
}
if (any(count_values != floor(count_values))) {
  stop("ASV counts must be whole read counts; fractional values were detected.")
}

# Report the dimensions read in: rows = samples, columns = ASVs
cat(
  "Imported ASV count table:\n",
  "- Samples:", nrow(asv_count_table), "\n",
  "- ASVs:", ncol(asv_count_table), "\n"
)
```

## Import ASV Sequence Table

Read the ASV ID-to-sequence mapping table.

``` r
# Read the ASV ID / sequence mapping table
asv_sequences <- read.csv(
  input_asv_sequences_path,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Validate that the expected columns are present before proceeding
required_columns <- c("ASV_ID", "sequence")
missing_columns <- setdiff(required_columns, colnames(asv_sequences))
if (length(missing_columns) > 0) {
  stop(
    "The ASV sequence table is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    ". Expected columns: ", paste(required_columns, collapse = ", ")
  )
}
if (nrow(asv_sequences) == 0L) {
  stop("data/asv_sequences.csv must contain at least one ASV sequence.")
}
if (anyDuplicated(asv_sequences$ASV_ID)) {
  stop("data/asv_sequences.csv contains duplicate ASV_ID values.")
}
if (any(is.na(asv_sequences$ASV_ID) | trimws(asv_sequences$ASV_ID) == "")) {
  stop("data/asv_sequences.csv contains an empty ASV_ID value.")
}
if (any(is.na(asv_sequences$sequence) | trimws(asv_sequences$sequence) == "")) {
  stop("data/asv_sequences.csv contains an empty sequence value.")
}

# Report the number of ASV sequences read in
cat("Imported", nrow(asv_sequences), "ASV sequences.\n")
```

------------------------------------------------------------------------

# BIOM Table Construction

## Transpose Count Table

[BIOM](http://biom-format.org/) tables are stored as observations
(features/ASVs) x samples, which is the transpose of the sample x ASV
layout produced by DADA2.

``` r
# Transpose so that rows = ASVs (features) and columns = samples
# t() on a data.frame returns a matrix, which is the format make_biom() expects
asv_count_matrix <- t(as.matrix(asv_count_table))

cat(
  "Transposed count matrix:\n",
  "- Features (ASVs):", nrow(asv_count_matrix), "\n",
  "- Samples:", ncol(asv_count_matrix), "\n"
)
```

## Validate ASV/Sequence Consistency

Before writing any output, confirm that every ASV in the count table has
a matching sequence, and vice versa. Silent mismatches here would
otherwise surface as a cryptic error deep inside PICRUSt2.

``` r
# ASV IDs present in the (transposed) count matrix
asv_ids_in_counts <- rownames(asv_count_matrix)

# ASV IDs present in the sequence table
asv_ids_in_sequences <- asv_sequences$ASV_ID

# ASVs with counts but no corresponding sequence — these cannot be placed by PICRUSt2
asvs_missing_sequence <- setdiff(asv_ids_in_counts, asv_ids_in_sequences)

# ASVs with a sequence but absent from the count table — harmless, but reported for transparency
asvs_missing_counts <- setdiff(asv_ids_in_sequences, asv_ids_in_counts)

if (length(asvs_missing_sequence) > 0) {
  stop(
    "Found ", length(asvs_missing_sequence), " ASV(s) in the count table with no matching sequence:\n",
    paste("  -", head(asvs_missing_sequence, 10), collapse = "\n"),
    if (length(asvs_missing_sequence) > 10) "\n  - ... (truncated)" else ""
  )
}

if (length(asvs_missing_counts) > 0) {
  warning(
    "Found ", length(asvs_missing_counts), " ASV sequence(s) with no corresponding entry in the count table. ",
    "These will be excluded from the FASTA file to keep it consistent with the BIOM table."
  )
  # Restrict the sequence table to only the ASVs actually present in the count table,
  # so the FASTA file and the BIOM table describe exactly the same feature set.
  asv_sequences <- asv_sequences[asv_sequences$ASV_ID %in% asv_ids_in_counts, ]
}

cat("ASV/sequence consistency check passed:", nrow(asv_sequences), "ASVs shared between both inputs.\n")
```

## Build and Export the BIOM Table

Convert the validated count matrix into a
[BIOM](http://biom-format.org/) object and write it to disk in BIOM 1.0
(JSON) format, as expected by picrust2.

``` r
# Convert the ASV x sample count matrix into a BIOM object
biom_object <- make_biom(asv_count_matrix)

# Write the BIOM object to disk
write_biom(biom_object, output_biom_path)

# Reopen the file from disk and compare it with the matrix used to create it.
# This validates the serialized artifact itself—not only the in-memory object—
# before Step 2 is allowed to consume it.
written_biom_matrix <- as.matrix(biom_data(read_biom(output_biom_path)))
biom_ids_match <- identical(
  dimnames(written_biom_matrix),
  dimnames(asv_count_matrix)
)
biom_counts_match <- isTRUE(all.equal(
  unname(written_biom_matrix),
  unname(asv_count_matrix),
  tolerance = 0,
  check.attributes = FALSE
))
if (!biom_ids_match || !biom_counts_match) {
  stop(
    "The BIOM file failed read-back validation: its identifiers or counts ",
    "do not match the validated ASV matrix."
  )
}

cat("BIOM table written to:\n", output_biom_path, "\n")
```

The reporting chunk below adds a portable link to the verified BIOM file
without repeating the full machine-specific path in the narrative.

------------------------------------------------------------------------

# FASTA File Construction

## Validate Sequences

Construct a
[Biostrings::DNAStringSet](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
from the sequence table. This both validates that every sequence
contains only legal IUPAC nucleotide codes and provides a robust,
standards-compliant object for writing FASTA output — avoiding the
pitfalls of manual string concatenation (e.g. stray whitespace,
inconsistent line endings, or malformed headers).

``` r
# Build a named DNAStringSet: names become FASTA headers, values become sequences
# This step will raise an informative error if any sequence contains characters
# outside the standard IUPAC nucleotide alphabet
asv_sequences_stringset <- DNAStringSet(setNames(asv_sequences$sequence, asv_sequences$ASV_ID))

# Report basic sequence length statistics as a sanity check
sequence_lengths <- width(asv_sequences_stringset)
cat(
  "Sequence validation passed:\n",
  "- Total sequences:", length(asv_sequences_stringset), "\n",
  "- Length range:", min(sequence_lengths), "-", max(sequence_lengths), "bp\n",
  "- Mean length:", round(mean(sequence_lengths), 1), "bp\n"
)
```

## Construct and Export the FASTA File

Write the validated sequences to
[study_seqs.fna](../../results/1_prepare_picrust2_inputs/study_seqs.fna),
as expected by PICRUSt2.

``` r
# writeXStringSet() writes a standards-compliant FASTA file directly from the
# validated DNAStringSet object
writeXStringSet(asv_sequences_stringset, filepath = output_fasta_path, format = "fasta")

# Reopen the FASTA file and verify both header order and sequence content.
# Matching the serialized file protects against malformed headers, incomplete
# writes, or an accidental mismatch between the FASTA and BIOM feature sets.
written_fasta <- readDNAStringSet(output_fasta_path, format = "fasta")
fasta_ids_match <- identical(names(written_fasta), names(asv_sequences_stringset))
fasta_sequences_match <- identical(
  as.character(written_fasta),
  as.character(asv_sequences_stringset)
)
if (!fasta_ids_match || !fasta_sequences_match) {
  stop(
    "The FASTA file failed read-back validation: its ASV headers or ",
    "sequences do not match the validated sequence set."
  )
}

cat("FASTA file written to:\n", output_fasta_path, "\n")
```

The next reporting chunk links the verified FASTA artifact that will be
passed to the `-s` argument of `picrust2_pipeline.py`.

------------------------------------------------------------------------

# Results Summary

## Compute Summary Statistics

This chunk assembles a one-row summary of the reformatting step
(input/output file names, sample and ASV counts, sequence length
statistics, and output file sizes), then previews it below as a tidy
Field/Value table before it is written to the Excel workbook.

``` r
# Assemble a one-row summary table capturing the key facts about this reformatting step
picrust2_input_summary <- data.frame(
  `Input ASV Count Table` = basename(input_asv_count_table_path),
  `Input ASV Sequence Table` = basename(input_asv_sequences_path),
  `Samples` = ncol(asv_count_matrix),
  `ASVs (Features)` = nrow(asv_count_matrix),
  `Sequence Length Min (bp)` = min(sequence_lengths),
  `Sequence Length Max (bp)` = max(sequence_lengths),
  `Sequence Length Mean (bp)` = round(mean(sequence_lengths), 1),
  `Total Reads` = sum(asv_count_matrix),
  `Output BIOM File` = basename(output_biom_path),
  `Output BIOM Size (KB)` = round(file.info(output_biom_path)$size / 1e3, 1),
  `Output FASTA File` = basename(output_fasta_path),
  `Output FASTA Size (KB)` = round(file.info(output_fasta_path)$size / 1e3, 1),
  check.names = FALSE
)

# Display the summary in the notebook. DT is namespace-qualified (not
# library()-loaded) since this is its only use in this notebook. The one-row
# summary above is reshaped into a tidy two-column Field/Value table first,
# since DT::datatable() (unlike kable()) is not designed to transpose a
# single-row data frame directly; the fixed-height scroll box matches the
# interactive DT::datatable() tables used elsewhere in this pipeline.
picrust2_input_summary_long <- data.frame(
  Field = names(picrust2_input_summary),
  Value = unlist(picrust2_input_summary, use.names = FALSE),
  check.names = FALSE
)
DT::datatable(
  picrust2_input_summary_long,
  options = list(dom = 't', pageLength = nrow(picrust2_input_summary_long),
                 scrollY = "400px", scrollCollapse = TRUE, paging = FALSE,
                 ordering = FALSE, searching = FALSE),
  rownames = FALSE,
  caption = "PICRUSt2 input preparation summary"
)
```

## Export Summary to Excel

This chunk writes the summary table above to the output workbook using
the project’s shared Excel-writing helper,
[add_sheet_to_excel_function.R](../functions/add_sheet_to_excel_function.R).

``` r
# Start each full notebook run with a clean workbook so stale sheets from an
# earlier run cannot survive. overwrite = TRUE also makes this chunk safe to
# rerun interactively after the workbook has been created.
if (file.exists(output_excel_path)) {
  unlink(output_excel_path)
}
add_sheet_to_excel(output_excel_path, output_excel_sheet, picrust2_input_summary,
                   rownames = FALSE, overwrite = TRUE)
```

This short reporting chunk confirms where the first workbook sheet was
written; the workbook is not complete until the column dictionary is
appended below.

## Document and Export Column Dictionary

Every column written to this workbook is documented in a trailing
`Column_Dictionary` sheet (Sheet / Column / Explanation), so the file is
self-explanatory to anyone opening it outside of this notebook.

``` r
# Human-readable explanation for every column written to the
# PICRUSt2_Input_Summary sheet, keyed by exact column name. build_column_dictionary()
# (sourced in Load Required Packages above) cross-checks this list against
# the actual columns of picrust2_input_summary and stops with an informative
# error if any column would otherwise be exported without documentation.
picrust2_input_summary_column_descriptions <- c(
  `Input ASV Count Table` =
    "File name of the DADA2-derived sample x ASV read-count matrix supplied as input (data/asv_count_table.csv).",
  `Input ASV Sequence Table` =
    "File name of the ASV ID-to-sequence mapping table supplied as input (data/asv_sequences.csv).",
  `Samples` =
    "Number of samples (columns) in the transposed ASV x sample count matrix written to the BIOM table.",
  `ASVs (Features)` =
    "Number of ASVs (features/rows) in the transposed ASV x sample count matrix written to the BIOM table.",
  `Sequence Length Min (bp)` =
    "Shortest representative ASV sequence length, in base pairs, among all sequences written to the FASTA file.",
  `Sequence Length Max (bp)` =
    "Longest representative ASV sequence length, in base pairs, among all sequences written to the FASTA file.",
  `Sequence Length Mean (bp)` =
    "Mean representative ASV sequence length, in base pairs, rounded to 1 decimal place.",
  `Total Reads` =
    "Sum of every read count in the ASV x sample count matrix, i.e. total reads across all samples and ASVs combined.",
  `Output BIOM File` =
    "File name of the BIOM-format feature table written for PICRUSt2 (study_seqs.biom).",
  `Output BIOM Size (KB)` =
    "File size of study_seqs.biom on disk, in kilobytes, rounded to 1 decimal place.",
  `Output FASTA File` =
    "File name of the FASTA file of representative ASV sequences written for PICRUSt2 (study_seqs.fna).",
  `Output FASTA Size (KB)` =
    "File size of study_seqs.fna on disk, in kilobytes, rounded to 1 decimal place."
)

# Build the Sheet / Column / Explanation table for this workbook's one data
# sheet. Passing the actual output_excel_sheet variable and picrust2_input_summary
# object (rather than re-typing the sheet name or column list) keeps this
# dictionary from drifting out of sync with what Export Summary to Excel above
# actually wrote.
column_dictionary <- build_column_dictionary(
  sheet_name    = output_excel_sheet,
  data          = picrust2_input_summary,
  descriptions  = picrust2_input_summary_column_descriptions,
  workbook_path = output_excel_path,
  rownames      = FALSE
)
# Deliberately not printed here (knitr::kable(column_dictionary)) -- the
# dictionary is exported to the workbook below and is meant to be read
# there, not duplicated in this notebook's own rendered output.
```

This chunk appends the completed column dictionary as the workbook’s
final sheet, then links the finished workbook in the report.

``` r
# Append the dictionary as the final sheet in the workbook. Because
# add_sheet_to_excel() appends worksheets in call order and this is the last
# call made against output_excel_path in this notebook, Column_Dictionary
# ends up as the trailing sheet in the saved .xlsx file.
add_sheet_to_excel(output_excel_path, "Column_Dictionary", column_dictionary,
                   rownames = FALSE, overwrite = TRUE)
```

The workbook is now complete. This reporting chunk confirms the final
path and provides a clickable link for users reading the rendered
report.

------------------------------------------------------------------------

# Output File Summary

The tree below lists every file this notebook has written to its own
output folder,
[results/1_prepare_picrust2_inputs/](../../results/1_prepare_picrust2_inputs/),
as a clickable, portable link (relative to this notebook’s own location)
with a short description – built live from what is actually on disk at
knit time via the project’s shared
[render_output_tree_function.R](../functions/render_output_tree_function.R)
helper, so it always matches this run’s real output rather than a
hand-maintained list.

------------------------------------------------------------------------

# Output Interpretation

- **`study_seqs.biom` and `study_seqs.fna` form one matched input
  pair.** Do not replace or edit only one of them: both files must
  describe the same ASV identifiers for PICRUSt2 sequence placement and
  abundance prediction to remain aligned.
- **Counts are unchanged.** This notebook only validates and transposes
  the DADA2 count table; it does not normalize, rarefy, transform, or
  filter the abundances.
- **Successful export confirms file integrity, not biological
  suitability.** Primer choice, sequence quality, taxonomic scope, and
  the distance to PICRUSt2 reference genomes still affect the
  predictions. Step 3 evaluates that last point through NSTI.
- **The Excel workbook is an audit summary.** Compare its sample count,
  ASV count, total reads, and sequence-length range with the upstream
  DADA2 results before continuing.

------------------------------------------------------------------------

# Recommended Next Step

The BIOM feature table and FASTA file exported above are the two
required inputs to PICRUSt2 itself. Proceed to [Step
2](2_picrust2_pipeline.md) — PICRUSt2 Pipeline, which reads
[study_seqs.biom](../../results/1_prepare_picrust2_inputs/study_seqs.biom)
and
[study_seqs.fna](../../results/1_prepare_picrust2_inputs/study_seqs.fna)
from
[results/1_prepare_picrust2_inputs/](../../results/1_prepare_picrust2_inputs/)
(this notebook’s output folder) and runs the PICRUSt2 functional
prediction pipeline to produce per-sample KO, EC, and MetaCyc pathway
abundance tables.

------------------------------------------------------------------------

# Session Information

Record the R environment for reproducibility.

------------------------------------------------------------------------

# References

## Methods

- McDonald D, Clemente JC, Kuczynski J, et al. (2012). The Biological
  Observation Matrix (BIOM) format or: how I learned to stop worrying
  and love the ome-ome. *GigaScience* 1, 7.
  <https://doi.org/10.1186/2047-217X-1-7> (defines the
  [BIOM](http://biom-format.org/) format written by this notebook.)
- Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for
  prediction of metagenome functions. *Nat Biotechnol* 38, 685-688.
  <https://doi.org/10.1038/s41587-020-0548-6>
- Pagès H, Aboyoun P, Gentleman R, DebRoy S. Biostrings: Efficient
  manipulation of biological strings. Bioconductor.
  <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>
  (provides `DNAStringSet()` and `writeXStringSet()`, used to validate
  sequences and write the FASTA output.)
- [biomformat Bioconductor
  package](https://bioconductor.org/packages/release/bioc/html/biomformat.html)
  — provides `make_biom()` and `write_biom()`, used to construct and
  write the BIOM table.

## Related

- [Step 2 — PICRUSt2 Pipeline](2_picrust2_pipeline.md) — this notebook’s
  recommended next step.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### Input File Not Found

**Error**:
`cannot open file '.../data/asv_count_table.csv': No such file or directory`

**Solutions**:

- Confirm [asv_count_table.csv](../../data/asv_count_table.csv) and
  [asv_sequences.csv](../../data/asv_sequences.csv) are present in
  [data/](../../data/).
- These files are generated by the DADA2 pipeline - copy them from in
  that project.

### ASVs Missing a Corresponding Sequence

**Error**: `Found N ASV(s) in the count table with no matching sequence`

**Possible causes**:

- The count table and sequence table were exported from different DADA2
  runs.
- The sequence table was truncated or filtered independently of the
  count table.

**Actions**:

- Re-export both files together from the same DADA2 pipeline run.
- Confirm `ASV_ID` values in
  [asv_sequences.csv](../../data/asv_sequences.csv) exactly match the
  column names of [asv_count_table.csv](../../data/asv_count_table.csv)
  (case-sensitive, no extra whitespace).

### Invalid Sequence Characters

**Error**: `DNAStringSet` construction fails with a message about
invalid characters

**Possible causes**:

- The `sequence` column contains non-IUPAC characters (e.g. stray commas
  from a malformed CSV export, or amino acid codes instead of
  nucleotides).

**Actions**:

- Inspect the offending rows:
  `asv_sequences[!grepl("^[ACGTURYSWKMBDHVN]+$", asv_sequences$sequence, ignore.case = TRUE), ]`
- Re-export `asv_sequences.csv` from DADA2, ensuring the sequence column
  contains only raw nucleotide strings.

### BIOM Table Shape Looks Wrong

**Symptom**: `study_seqs.biom` shape does not match the expected
`(n_ASVs, n_samples)`.

**Possible causes**:

- The count table was not transposed before calling `make_biom()`.
- `read.csv()` picked up an unexpected column as row names (e.g. if the
  first column header was not the sample ID column).

**Actions**:

- Confirm `dim(asv_count_matrix)` reports `(ASVs, samples)`, not
  `(samples, ASVs)`.
- Open [asv_count_table.csv](../../data/asv_count_table.csv) and verify
  the first column contains sample identifiers.
