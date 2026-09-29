#!/usr/bin/env Rscript

# Run all five workflow steps against the bundled example. Inputs, generated
# outputs, and rendered reports remain below example/ so this test never reads
# from data/ or writes to the normal results/ directory.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Run this file with Rscript: Rscript example/run_example.R", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

if (length(commandArgs(trailingOnly = TRUE))) {
  stop("This runner accepts no arguments. Use: Rscript example/run_example.R", call. = FALSE)
}

data_root <- file.path(project_root, "example", "data")
run_results <- file.path(project_root, "example", "run_results")
reference_results <- file.path(project_root, "example", "reference_results")
report_directory <- file.path(run_results, "reports")

Sys.setenv(
  PICRUST2_DATA_DIR = data_root,
  PICRUST2_RESULTS_DIR = run_results,
  PICRUST2_REUSE_EXISTING_OUTPUTS = "false"
)

required_packages <- c(
  "rmarkdown", "here", "fs", "biomformat", "Biostrings", "readr",
  "dplyr", "tidyr", "tibble", "ggplot2", "plotly", "htmlwidgets",
  "DT", "openxlsx", "ggpicrust2", "MicrobiomeStat", "Maaslin2",
  "data.table", "htmltools", "KEGGREST", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the workflow R dependencies first; missing: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

required_inputs <- file.path(data_root, c(
  "asv_count_table.csv", "asv_sequences.csv", "asv_taxonomy.csv", "metadata.tsv"
))
if (!all(file.exists(required_inputs))) {
  stop(
    "The bundled example is incomplete. Missing:\n  - ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "\n  - "),
    call. = FALSE
  )
}

count_table <- utils::read.csv(required_inputs[[1]], check.names = FALSE)
sequence_table <- utils::read.csv(required_inputs[[2]], check.names = FALSE)
taxonomy_table <- utils::read.csv(required_inputs[[3]], check.names = FALSE)
metadata <- utils::read.delim(required_inputs[[4]], check.names = FALSE)

if (!identical(names(count_table)[[1]], "SampleID") ||
    !all(c("ASV_ID", "sequence") %in% names(sequence_table)) ||
    !identical(names(taxonomy_table)[[1]], "ASV_ID") ||
    !all(c("SampleID", "Condition", "Reference") %in% names(metadata))) {
  stop("The example files do not follow the documented input schemas.", call. = FALSE)
}
if (!identical(sort(count_table$SampleID), sort(metadata$SampleID))) {
  stop("Example sample identifiers differ between counts and metadata.", call. = FALSE)
}
asv_ids <- names(count_table)[-1]
if (!identical(sort(asv_ids), sort(sequence_table$ASV_ID)) ||
    !all(asv_ids %in% taxonomy_table$ASV_ID)) {
  stop("Example ASV identifiers are inconsistent across input files.", call. = FALSE)
}

reference_by_condition <- tapply(
  as.logical(metadata$Reference), metadata$Condition, function(x) unique(x)
)
if (any(lengths(reference_by_condition) != 1L) ||
    sum(vapply(reference_by_condition, isTRUE, logical(1))) != 1L) {
  stop("Exactly one complete example condition must be the reference.", call. = FALSE)
}

if (!rmarkdown::pandoc_available()) {
  quarto_executable <- Sys.which("quarto")
  if (nzchar(quarto_executable)) {
    pandoc_candidates <- list.files(
      file.path(dirname(normalizePath(quarto_executable)), "tools"),
      pattern = "^pandoc$", recursive = TRUE, full.names = TRUE
    )
    pandoc_candidates <- pandoc_candidates[file.access(pandoc_candidates, 1L) == 0L]
    if (length(pandoc_candidates)) Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_candidates[[1]]))
  }
}
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required. Run from RStudio or install Quarto/Pandoc.", call. = FALSE)
}

conda_candidates <- unique(c(
  Sys.which("conda"), "/opt/anaconda3/bin/conda", "/opt/homebrew/bin/conda"
))
conda_candidates <- conda_candidates[nzchar(conda_candidates)]
if (!any(file.access(conda_candidates, 1L) == 0L)) {
  stop(
    "Conda with a PICRUSt2 environment is required for Steps 2 and 5. Run ",
    "setup/install_picrust2.sh first.", call. = FALSE
  )
}

if (dir.exists(run_results)) unlink(run_results, recursive = TRUE, force = TRUE)
dir.create(report_directory, recursive = TRUE, showWarnings = FALSE)

render_step <- function(filename) {
  message("\nRendering ", filename, " ...")
  rmarkdown::render(
    input = file.path(project_root, "R", "notebooks", filename),
    output_format = "html_document",
    output_dir = report_directory,
    envir = new.env(parent = globalenv()),
    clean = TRUE,
    quiet = FALSE
  )
}

render_step("1_prepare_picrust2_inputs.Rmd")
render_step("2_picrust2_pipeline.Rmd")
render_step("3_nsti_quality_assessment.Rmd")
render_step("4_functional_differential_abundance_analysis.Rmd")
render_step("5_taxon_contribution.Rmd")

# Output links are generated relative to R/notebooks/, where reports normally
# live. Rewrite them for the isolated report directory and remove the local
# checkout path from displayed console output before reports are shared.
report_files <- list.files(report_directory, pattern = "[.]html$", full.names = TRUE)
github_blob_root <- paste0(
  "https://github.com/changlabs/",
  "PICRUSt2_16S_Functional_Inference_Workflow/blob/main/"
)
github_tree_root <- paste0(
  "https://github.com/changlabs/",
  "PICRUSt2_16S_Functional_Inference_Workflow/tree/main/"
)
for (report_file in report_files) {
  report_html <- readLines(report_file, warn = FALSE, encoding = "UTF-8")
  report_html <- gsub(
    'href="../../example/run_results/', 'href="../', report_html, fixed = TRUE
  )
  report_html <- gsub(
    'href="../../results/', 'href="../', report_html, fixed = TRUE
  )
  # Reports are rendered outside R/notebooks/, so links that were originally
  # relative to that source directory must become stable repository links.
  report_html <- gsub(
    'href="../../README.md"',
    paste0('href="', github_blob_root, 'README.md"'),
    report_html, fixed = TRUE
  )
  report_html <- gsub(
    'href="../../setup/([^"#]+)(#[^"]*)?"',
    paste0('href="', github_blob_root, 'setup/\\1\\2"'),
    report_html, perl = TRUE
  )
  report_html <- gsub(
    'href="../../data/README.md"',
    paste0('href="', github_blob_root, 'data/README.md"'),
    report_html, fixed = TRUE
  )
  report_html <- gsub(
    'href="../../data/([^"#]+)(#[^"]*)?"',
    paste0('href="', github_blob_root, 'example/data/\\1\\2"'),
    report_html, perl = TRUE
  )
  report_html <- gsub(
    'href="../../data/"',
    paste0('href="', github_tree_root, 'example/data/"'),
    report_html, fixed = TRUE
  )
  report_html <- gsub(
    'href="../functions/([^"#]+)(#[^"]*)?"',
    paste0('href="', github_blob_root, 'R/functions/\\1\\2"'),
    report_html, perl = TRUE
  )
  report_html <- gsub(
    'href="../functions/"',
    paste0('href="', github_tree_root, 'R/functions/"'),
    report_html, fixed = TRUE
  )
  report_html <- gsub(
    'href="([1-5]_[^"#]+[.](?:Rmd|md))(#[^"]*)?"',
    paste0('href="', github_blob_root, 'R/notebooks/\\1\\2"'),
    report_html, perl = TRUE
  )
  report_html <- gsub(project_root, "PROJECT_ROOT", report_html, fixed = TRUE)
  report_html <- gsub(
    utils::URLencode(project_root, reserved = FALSE),
    "PROJECT_ROOT", report_html, fixed = TRUE
  )
  writeLines(report_html, report_file, useBytes = TRUE)
}

required_outputs <- c(
  file.path(run_results, "1_prepare_picrust2_inputs", c(
    "study_seqs.biom", "study_seqs.fna", "picrust2_input_summary.xlsx"
  )),
  file.path(run_results, "2_picrust2_pipeline", "picrust2_out_pipeline", c(
    "combined_marker_predicted_and_nsti.tsv",
    "KO_metagenome_out/pred_metagenome_unstrat.tsv",
    "KO_metagenome_out/weighted_nsti.tsv",
    "EC_metagenome_out/pred_metagenome_unstrat.tsv",
    "pathways_out/path_abun_unstrat.tsv"
  )),
  file.path(run_results, "3_nsti_quality_assessment", "nsti_quality_summary.xlsx"),
  file.path(run_results, "4_functional_differential_abundance_analysis", c(
    "functional_daa_summary_MetaCyc.xlsx", "functional_daa_summary_KEGG.xlsx"
  )),
  file.path(run_results, "5_taxon_contribution", c(
    "MetaCyc/taxon_contribution_summary.xlsx",
    "KEGG/taxon_contribution_summary.xlsx"
  )),
  file.path(report_directory, paste0(c(
    "1_prepare_picrust2_inputs", "2_picrust2_pipeline",
    "3_nsti_quality_assessment", "4_functional_differential_abundance_analysis",
    "5_taxon_contribution"
  ), ".html"))
)
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs)) {
  stop(
    "The example run did not produce the complete required deliverable set:\n  - ",
    paste(missing_outputs, collapse = "\n  - "),
    call. = FALSE
  )
}

message(
  "\nExample run complete.\n",
  "Completed: all Steps 1-5, including both MetaCyc and KEGG analyses.\n",
  "Generated results: ", run_results, "\n",
  "Bundled reference results: ", reference_results, "\n",
  "The normal data/ and results/ directories were not touched."
)
