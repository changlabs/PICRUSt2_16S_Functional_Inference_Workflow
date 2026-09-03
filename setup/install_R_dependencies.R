################################################################################
# Script: install_R_dependencies.R
# Purpose:
#   - Install all R packages required by the workflow's five R notebooks,
#     covering data preparation, reporting, visualization, differential-
#     abundance analysis, and taxon-contribution analysis.
#   - Provide clear progress feedback at every installation step so the user
#     can follow along and quickly identify failures.
#   - Verify that every package loaded correctly after installation and print a
#     human-readable summary showing how many packages succeeded or failed.
#   - Report the installed version of each key package and the active R /
#     Bioconductor version for reproducibility documentation.
#
# When to run this script:
#   Run this script ONCE after cloning or downloading the repository for the
#   first time, or whenever you set up a new R environment (e.g., a new conda
#   environment, a new RStudio project on a different machine, or after a major
#   R version upgrade). You do NOT need to re-run it on every analysis session.
#
#   Note: this installs the R packages used by this repository's R notebooks
#   (e.g. R/notebooks/1_prepare_picrust2_inputs.Rmd). PICRUSt2 itself is a separate,
#   Python-based tool installed via conda — see setup/install_picrust2.sh in
#   the repository root, not this script.
#
# Prerequisites:
#   - R >= 4.1.0
#   - An active internet connection (CRAN and Bioconductor servers must be reachable)
#
# Output:
#   - Packages are installed into your active R library path (.libPaths()[1]).
#   - No files are written to disk by this script; it only modifies the R library.
#   - A final summary table is printed to the console showing success/failure
#     status and version strings for the key packages.
#
# Troubleshooting:
#   - If a CRAN package fails, try installing it manually:
#       install.packages("<package_name>", dependencies = TRUE)
#   - If a Bioconductor package fails, ensure BiocManager is up to date:
#       install.packages("BiocManager")
#       BiocManager::install("<package_name>", update = FALSE, ask = FALSE)
#   - On Linux, some packages (e.g., xml2, curl, openssl) require system
#     libraries. Install them with your package manager first, for example:
#       sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev
#
################################################################################


# ==============================================================================
# Print banner
# ==============================================================================
# A visual banner is printed at startup so the user can immediately confirm that
# the correct script is running. This is especially helpful when sourcing from a
# larger workflow or calling via Rscript in a batch job.
cat("
╔══════════════════════════════════════════════════════════════════════════════╗
║           PICRUSt2 WORKFLOW - R DEPENDENCY INSTALLATION                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")

minimum_r_version <- package_version("4.1.0")
if (getRversion() < minimum_r_version) {
    stop(
        "This workflow requires R >= ", minimum_r_version,
        "; current version: ", getRversion()
    )
}

# Use a reproducible non-interactive CRAN mirror unless the caller already set
# one. PICRUST2_CRAN_REPO provides a configuration override without editing the
# installer.
configured_cran_repository <- Sys.getenv(
    "PICRUST2_CRAN_REPO", unset = "https://cloud.r-project.org"
)
if (!nzchar(trimws(configured_cran_repository))) {
    stop("PICRUST2_CRAN_REPO must not be empty.")
}
options(repos = c(CRAN = configured_cran_repository))


# ==============================================================================
# Helper function: install_if_missing
# ==============================================================================
# Rather than calling install.packages() or BiocManager::install() directly
# everywhere in the script, we centralise the install-and-verify logic in this
# single helper. It checks whether the package is already present (skipping the
# download if it is), attempts the appropriate install command, and immediately
# confirms whether the install succeeded.
#
# Arguments:
#   packages  : Character vector of package names to check and install.
#   source    : Either "CRAN" (default) or "Bioconductor". Controls which
#               backend is used for packages not already present.
#
# Behaviour:
#   - Packages that are already available produce a "✓ already installed" line.
#   - Packages that need installing produce "Installing: <name> ..." followed
#     by either "✓ installed successfully" or "✗ Failed to install <name>".
#   - The function does NOT stop the script on failure; it merely reports the
#     failure so all packages are attempted before the user reviews the summary.
install_if_missing <- function(packages, source = "CRAN") {
    if (!(source %in% c("CRAN", "Bioconductor"))) {
        stop("install_if_missing(): source must be 'CRAN' or 'Bioconductor'.")
    }
    for (pkg in packages) {

        # -------------------------------------------------------------------
        # Check availability before attempting installation
        # -------------------------------------------------------------------
        # requireNamespace() loads the package namespace without attaching it,
        # making it the lightweight way to test whether a package is present.
        # quietly = TRUE suppresses the "there is no package called ..." message
        # so our own formatted output remains clean.
        if (!requireNamespace(pkg, quietly = TRUE)) {

            cat("Installing:", pkg, "...\n")

            # ---------------------------------------------------------------
            # Install from the appropriate source
            # ---------------------------------------------------------------
            tryCatch(
                {
                    if (source == "CRAN") {
                        # dependencies = TRUE installs optional packages needed by
                        # package features as well as mandatory dependencies.
                        install.packages(pkg, dependencies = TRUE, quiet = TRUE)

                    } else {
                        # update = FALSE avoids unrelated upgrades; ask = FALSE
                        # permits unattended installation.
                        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
                    }
                },
                error = function(error) {
                    warning("Installation attempt failed for ", pkg, ": ",
                            conditionMessage(error))
                }
            )

            # ---------------------------------------------------------------
            # Verify that the install succeeded
            # ---------------------------------------------------------------
            # We call requireNamespace() again after installation. If it still
            # returns FALSE, the install silently failed (e.g., a missing system
            # library, network error, or compilation error).
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat("  ✓", pkg, "installed successfully\n")
            } else {
                cat("  ✗ Failed to install", pkg, "\n")
            }

        } else {
            # Package was already present; skip the download.
            cat("  ✓", pkg, "already installed\n")
        }
    }
}


# ==============================================================================
# Step 1 of 4 — Install BiocManager
# ==============================================================================
# BiocManager is the official gateway to Bioconductor packages. It must be
# present before we can call BiocManager::install() for biomformat, Biostrings,
# etc. Installing it from CRAN first is safe even on systems that already have
# it; install.packages() will silently do nothing if the package is up to date.
cat("\n[1/4] Checking BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    install.packages("BiocManager", quiet = TRUE)
}
cat("  ✓ BiocManager available\n")


# ==============================================================================
# Step 2 of 4 — CRAN packages
# ==============================================================================
# These packages are all available on the Comprehensive R Archive Network
# (CRAN) and support project-relative paths, reporting, data manipulation,
# visualization, and the functional-analysis notebooks:
#
#   File system / paths: here
#   Excel & reporting  : openxlsx, knitr, DT, fs, htmltools, rmarkdown
#
# All packages are passed as a single character vector to `install_if_missing`
# so that the loop handles each one consistently.
cat("\n[2/4] Installing CRAN packages...\n")

cran_packages <- c(
    # ------------------------------------------------------------------
    # File system and paths
    # ------------------------------------------------------------------
    # here        : Constructs paths relative to the project root using
    #               .here or .Rproj anchors; prevents hardcoded absolute paths.
    "here",

    # ------------------------------------------------------------------
    # Excel and reporting
    # ------------------------------------------------------------------
    # openxlsx    : Read/write .xlsx files without Java; used by
    #               add_sheet_to_excel_function.R for per-step summary reports.
    "openxlsx",

    # knitr       : Dynamic report generation that embeds R code inside
    #               R Markdown notebooks, and provides kable() for
    #               formatted summary tables.
    "knitr",

    # DT          : Interactive, sortable/searchable HTML tables, used for
    #               summary, per-sample, and per-ASV tables in Steps 1-3 that can grow
    #               to hundreds of rows in a real run. Was already used via
    #               library(DT) in those notebooks but was missing from this
    #               installer -- added here so a fresh environment set up
    #               from this script alone has everything those notebooks need.
    "DT",

    # fs          : Portable file/directory inspection and relative-path
    #               construction used by the shared output-link/tree helpers.
    "fs",

    # htmltools   : Escapes output-tree link text safely and combines multiple
    #               Step 5 heatmap widgets in the rendered notebook.
    "htmltools",

    # ------------------------------------------------------------------
    # Data import and wrangling (used from Step 3 onward for NSTI quality
    # assessment and differential-abundance analysis)
    # ------------------------------------------------------------------
    # readr       : Fast, consistent delimited-file import (metadata.tsv,
    #               PICRUSt2's .tsv/.tsv.gz output tables); transparently
    #               reads gzip-compressed files without manual decompression.
    "readr",

    # dplyr       : Data manipulation grammar (filter/mutate/summarise/join)
    #               used throughout condition-wise QC and statistics steps.
    "dplyr",

    # tidyr       : Reshaping data between wide and long formats, needed for
    #               per-sample/per-condition summary tables and plotting inputs.
    "tidyr",

    # tibble      : Explicitly loaded by Step 4 for row-name/column
    #               conversions.
    "tibble",

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    # ggplot2     : Grammar-of-graphics plotting; the static plot object that
    #               plotly::ggplotly() converts into an interactive widget.
    "ggplot2",

    # plotly      : Converts ggplot2 plots into interactive (zoomable,
    #               hoverable) HTML widgets, embedded directly in the
    #               notebook's HTML output.
    "plotly",

    # htmlwidgets : Saves plotly (and other JavaScript-based) widgets as
    #               standalone, self-contained HTML files that can be shared
    #               or opened independently of the notebook.
    "htmlwidgets",

    # scales      : Percentage-axis formatting in Step 3's NSTI plot.
    "scales",

    # rmarkdown   : Detects Pandoc when Step 5 saves self-contained interactive
    #               heatmaps and renders HTML/GitHub notebook documents.
    "rmarkdown"
)

install_if_missing(cran_packages, source = "CRAN")


# ==============================================================================
# Step 3 of 4 — Functional-analysis packages
# ==============================================================================
# These packages support the Step 4 differential-abundance analysis and the
# Step 5 taxon-contribution analysis.
cat("\n[3/4] Installing functional-analysis CRAN packages...\n")

functional_analysis_packages <- c(
    # ------------------------------------------------------------------
    # Core functional analysis
    # ------------------------------------------------------------------
    # MicrobiomeStat: Provides the LinDA backend used by Step 4.
    "MicrobiomeStat",

    "data.table"   # Memory-efficient streaming import used by Step 5
)

install_if_missing(functional_analysis_packages, source = "CRAN")

# Step 4 depends on corrected KEGG upper-half aggregation and effective
# prokaryote filtering introduced in ggpicrust2 2.5.16. install_if_missing()
# intentionally skips installed packages, so enforce this minimum separately
# and upgrade an older installation before Step 4 can be run.
minimum_ggpicrust2_version <- package_version("2.5.16")
if (!requireNamespace("ggpicrust2", quietly = TRUE) ||
    packageVersion("ggpicrust2") < minimum_ggpicrust2_version) {
    cat("Updating ggpicrust2 to version", as.character(minimum_ggpicrust2_version), "or newer...\n")
    # The corrected release may lead CRAN's macOS binary, so use the
    # maintainer's R-universe repository first and CRAN for dependencies.
    ggpicrust2_repositories <- c(
        ggpicrust2 = "https://cafferychen777.r-universe.dev",
        CRAN       = configured_cran_repository
    )
    install.packages(
        "ggpicrust2",
        repos = ggpicrust2_repositories,
        dependencies = NA,
        quiet = TRUE
    )
}
if (!requireNamespace("ggpicrust2", quietly = TRUE) ||
    packageVersion("ggpicrust2") < minimum_ggpicrust2_version) {
    stop(
        "ggpicrust2 >= ", minimum_ggpicrust2_version,
        " is required for reliable Step 4 KEGG aggregation/filtering. Installed: ",
        if (requireNamespace("ggpicrust2", quietly = TRUE)) as.character(packageVersion("ggpicrust2")) else "not installed"
    )
}

# ==============================================================================
# Step 4 of 4 — Bioconductor packages
# ==============================================================================
# These packages handle microbiome-specific file formats, sequence operations,
# pathway annotation, and the selected statistical backend used by Step 4.
#
#   biomformat  : Reads and writes BIOM-format feature tables. Used to convert
#                 the ASV x sample count matrix into study_seqs.biom, the
#                 -i input expected by picrust2_pipeline.py.
#                 Reference: https://github.com/joey711/biomformat
#
#   Biostrings  : Efficient manipulation and validation of biological sequences
#                 (DNA, RNA, amino acid). Used to validate ASV sequences against
#                 the IUPAC nucleotide alphabet and to write the standards-
#                 compliant study_seqs.fna FASTA file, the -s input expected by
#                 picrust2_pipeline.py.
#
#   KEGGREST    : Completes human-readable names for selected KEGG pathways
#                 when a name is absent from the bundled local reference data.
cat("\n[4/4] Installing Bioconductor packages...\n")

bioc_packages <- c(
    "biomformat",      # BIOM-format feature table I/O
    "Biostrings",      # Biological sequence validation and FASTA I/O
    "Maaslin2",        # Complementary transformed linear-model DAA backend
    "KEGGREST"         # Optional KEGG pathway-name completion in Step 5
)

install_if_missing(bioc_packages, source = "Bioconductor")


# ==============================================================================
# Installation summary
# ==============================================================================
# After all installation attempts have completed, we compile the full list of
# packages and use requireNamespace() to test each one. The results are counted
# and displayed in a clear pass/fail format so the user can immediately see
# whether any packages need manual attention.
cat("\n
╔══════════════════════════════════════════════════════════════════════════════╗
║                         INSTALLATION SUMMARY                                 ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")

# Combine all package vectors into one master list for the summary check.
all_packages <- unique(c(
    cran_packages, functional_analysis_packages, "ggpicrust2", bioc_packages
))

# sapply returns a named logical vector: TRUE if the package loaded, FALSE if not.
installed <- sapply(all_packages, requireNamespace, quietly = TRUE)

cat("Successfully installed:", sum(installed), "/", length(all_packages), "packages\n\n")

if (all(installed)) {
    # All packages are present; the environment is ready for analysis.
    cat("✓ All packages installed successfully!\n")
    cat("\nYou can now run the pipeline notebooks.\n")
} else {
    # At least one package is missing. List each failed package so the user
    # knows exactly what to investigate or install manually.
    cat("✗ Some packages failed to install:\n")
    cat(paste("  -", all_packages[!installed], collapse = "\n"), "\n")
    cat("\nTry installing failed packages manually.\n")
}


# ==============================================================================
# Installed version report
# ==============================================================================
# For reproducibility and collaborative troubleshooting it is important to know
# exactly which versions are active. We print version strings for the key
# packages, plus the R and Bioconductor release versions. This output can be
# copied directly into a Methods section or a lab notebook.
cat("\n
══════════════════════════════════════════════════════════════════════════════
                          INSTALLED VERSIONS
══════════════════════════════════════════════════════════════════════════════
\n")

# Check and report version for each key package.
# packageVersion() returns a package_version object; we coerce it to character
# for clean printing.
key_packages <- c("biomformat", "Biostrings", "here", "openxlsx", "ggplot2",
                  "plotly", "dplyr", "ggpicrust2", "htmltools",
                  "MicrobiomeStat", "Maaslin2", "KEGGREST", "data.table",
                  "rmarkdown")
for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        ver <- as.character(packageVersion(pkg))
        # sprintf with %-15s left-aligns the package name in a 15-character field
        # so all version numbers line up in a tidy column.
        cat(sprintf("%-15s %s\n", pkg, ver))
    }
}

# Report the R and Bioconductor release versions.
# These are essential context for anyone trying to reproduce the analysis on a
# different machine or at a later date.
cat("\nR version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
