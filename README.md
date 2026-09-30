# PICRUSt2 16S Functional Inference Workflow

[![R Version](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://www.r-project.org/) [![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![PICRUSt2](https://img.shields.io/badge/PICRUSt2-functional%20inference-087f86)](https://github.com/picrust/picrust2) [![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)](#setup)

A reproducible R-based workflow for converting DADA2 ASV outputs into [PICRUSt2](https://github.com/picrust/picrust2) functional predictions, evaluating prediction quality with NSTI (Nearest Sequenced Taxon Index), testing MetaCyc and KEGG pathway differential abundance, and tracing selected pathway predictions back to contributing taxa.

PICRUSt2 results represent **predicted functional potential** inferred from marker-gene data. They are not direct measurements of genes, transcripts, proteins, metabolites, or pathway activity.

------------------------------------------------------------------------

## Table of Contents

- [Key Features](#key-features)
- [Pipeline Overview](#pipeline-overview)
- [Bundled Example Dataset](#bundled-example-dataset)
- [Setup](#setup)
  - [1. Clone the Repository](#1-clone-the-repository)
  - [2. Open the R Project in RStudio](#2-open-the-r-project-in-rstudio)
  - [3. Install R Dependencies](#3-install-r-dependencies)
  - [4. Install PICRUSt2](#4-install-picrust2)
  - [5. Add Your Input Data](#5-add-your-input-data)
- [Running the Pipeline](#running-the-pipeline)
  - [Step 1 — Prepare PICRUSt2 Inputs](#step-1--prepare-picrust2-inputs-required)
  - [Step 2 — Run the PICRUSt2 Pipeline](#step-2--run-the-picrust2-pipeline-required)
  - [Step 3 — NSTI Quality Assessment](#step-3--nsti-quality-assessment-recommended)
  - [Step 4 — Functional Differential Abundance](#step-4--functional-differential-abundance-optional)
  - [Step 5 — Taxon Contribution Analysis](#step-5--taxon-contribution-analysis-optional)
- [Column Dictionaries](#column-dictionaries)
- [Project Structure](#project-structure)
- [References](#references)
- [License](#license)
- [Acknowledgments](#acknowledgments)

------------------------------------------------------------------------

## Key Features

- **DADA2-compatible inputs** — Validates sample-by-ASV count data and representative sequences before creating matched BIOM and FASTA files
- **Complete PICRUSt2 prediction** — Produces marker-copy-number, KO, EC, MetaCyc pathway, and NSTI outputs using the standard PICRUSt2 pipeline
- **Prediction-quality assessment** — Summarizes per-ASV NSTI, abundance-weighted per-sample NSTI, and the fraction of reads excluded by the configured NSTI threshold
- **MetaCyc and KEGG analysis** — Runs both functional levels throughout differential-abundance and taxon-contribution analyses
- **Complementary differential-abundance methods** — Reports LinDA and MaAsLin2 separately, with Benjamini-Hochberg adjusted p-values and descriptive method agreement
- **Stratified taxon contributions** — Links selected MetaCyc and KEGG pathway predictions to contributing ASVs, genera, and phyla
- **Interactive visualizations** — Generates shareable Plotly NSTI plots and condition-split genus/phylum contribution heatmaps
- **Documented Excel outputs** — Every workbook ends with a `Column_Dictionary` sheet explaining its exported columns
- **Reproducible paths and reports** — All code uses project-relative paths, and every notebook supports an executable HTML report plus a non-executing GitHub Markdown document

------------------------------------------------------------------------

## Pipeline Overview

``` text
        DADA2 ASV count, sequence,
        taxonomy, and metadata files
                     │
                     ▼
        ┌──────────────────────────┐
        │     Step 1 (Required)    │
        │  Prepare PICRUSt2 Inputs │
        │       BIOM + FASTA       │
        └────────────┬─────────────┘
                     ▼
        ┌──────────────────────────┐
        │     Step 2 (Required)    │
        │   Run PICRUSt2 Pipeline  │
        │ KO, EC, MetaCyc, & NSTI  │
        └────────────┬─────────────┘
                     ▼
        ┌──────────────────────────┐
        │   Step 3 (Recommended)   │
        │ NSTI Quality Assessment  │
        └────────────┬─────────────┘
                     ▼
        ┌──────────────────────────┐
        │     Step 4 (Optional)    │
        │ Functional Differential  │
        │   Abundance: MetaCyc +   │
        │           KEGG           │
        └────────────┬─────────────┘
                     ▼
        ┌──────────────────────────┐
        │     Step 5 (Optional)    │
        │    Taxon Contributions   │
        │   MetaCyc + KEGG paths   │
        └────────────┬─────────────┘
                     ▼
          Interpretable predicted
          functional profiles and
             contribution reports
```

[Steps 1](R/notebooks/1_prepare_picrust2_inputs.md) and [2](R/notebooks/2_picrust2_pipeline.md) are required to generate the functional predictions. [Step 3](R/notebooks/3_nsti_quality_assessment.md) is a recommended quality-assessment checkpoint before biological interpretation. [Step 4](R/notebooks/4_functional_differential_abundance_analysis.md) is optional and performs the workflow's inferential analysis. [Step 5](R/notebooks/5_taxon_contribution.md) is optional for normal studies and generates its own stratified PICRUSt2 run to explain which taxa contribute to selected pathway predictions. The bundled example runs all five steps.

------------------------------------------------------------------------

## Bundled Example Dataset

The repository includes an isolated, clone-ready example under [`example/data/`](example/data/). It contains 10 samples, 1,468 representative 16S ASVs, a sample-by-ASV count table, metadata in the workflow's normal TSV layout, and a DADA2-compatible taxonomy table. The `Control`/`Treatment` design and all taxonomy labels are synthetic demonstration annotations intended only to exercise the workflow. They are not suitable for biological interpretation.

Run the complete example from the repository root:

``` bash
Rscript example/run_example.R
```

The runner executes **all five steps**, including NSTI assessment, both MetaCyc and KEGG differential-abundance analyses, and both taxon-contribution branches. It reads only [`example/data/`](example/data/) and writes only to the ignored `example/run_results/` directory, so it cannot mix with files placed in the normal `data/` directory or overwrite normal `results/`.

To inspect the example without installing PICRUSt2 or running the workflow, browse the committed reports under [`example/reference_results/reports/`](example/reference_results/reports/), the compact result files under [`example/reference_results/`](example/reference_results/), or the GitHub Pages site.

The abundance and sequence tables were derived from human stool samples generated in our laboratory, then randomly reduced and anonymized for demonstration purposes. The analysis methods should be cited using the PICRUSt2 and ggpicrust2 papers listed under [References](#references).

------------------------------------------------------------------------

## Setup

### 1. Clone the Repository

Download the repository from GitHub or clone it from the terminal:

``` bash
git clone https://github.com/changlabs/PICRUSt2_16S_Functional_Inference_Workflow.git
```

### 2. Open the R Project in RStudio

Open [PICRUSt2_16S_Functional_Inference_Workflow.Rproj](PICRUSt2_16S_Functional_Inference_Workflow.Rproj) by double-clicking it, or from inside [RStudio](https://posit.co/products/open-source/rstudio/) via **File → Open Project**. The notebooks use `here::here()` to resolve paths relative to the project root, so run them from within this project.

### 3. Install R Dependencies

Open [setup/install_R_dependencies.R](setup/install_R_dependencies.R) in RStudio and run it with **Source**. It installs the CRAN and [Bioconductor](https://bioconductor.org/) packages used across the workflow, including `biomformat`, `Biostrings`, `ggpicrust2`, `MicrobiomeStat`, `Maaslin2`, `KEGGREST`, `data.table`, `openxlsx`, `DT`, `ggplot2`, and `plotly`.

Run the installer once for each new R environment. On Linux, some packages may first require system libraries:

``` bash
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev
```

### 4. Install PICRUSt2

[PICRUSt2](https://github.com/picrust/picrust2) is installed separately because it is distributed through Bioconda and depends on compiled phylogenetic-placement tools. A working conda or Miniconda installation is required.

On Linux, run:

``` bash
bash setup/install_picrust2.sh
```

The installer creates a conda environment named `picrust2` by default. Set `PICRUST2_CONDA_ENV_NAME` before running the installer and notebooks if you need a different environment name. The supplied installer is intentionally Linux-only; on macOS, provide an existing compatible PICRUSt2 conda environment and configure the notebooks to use it.

### 5. Add Your Input Data

Copy the following study files into `data/`:

| File | Used by | Required content |
|:-----------------------|:-----------------------|:-----------------------|
| `asv_count_table.csv` | Step 1 | Sample-by-ASV integer read-count matrix |
| `asv_sequences.csv` | Step 1 | `ASV_ID` and representative nucleotide sequence |
| `metadata.tsv` | Steps 3–5 | `SampleID`, `Condition`, and `Reference` |
| `asv_taxonomy.csv` | Step 5 | `ASV_ID` plus standard taxonomic-rank columns |

Normal study files are not included in the repository and are ignored by Git. The bundled example lives separately under `example/data/`. See [data/README.md](data/README.md) for the complete schemas and consistency requirements.

------------------------------------------------------------------------

## Running the Pipeline

Run the applicable numbered notebooks in order from the project root. In RStudio, use **Knit** with the HTML output to execute the analysis and create a complete interactive report. Select **Knit to github_document** when you want to refresh the repository Markdown without executing the analysis. The bundled example instead uses `Rscript example/run_example.R` so its inputs and outputs remain isolated.

### Step 1 — [Prepare PICRUSt2 Inputs](R/notebooks/1_prepare_picrust2_inputs.md) *(required)*

Validates the ASV count and sequence tables, confirms their identifiers agree, transposes the abundance matrix into PICRUSt2 orientation, and creates `study_seqs.biom` and `study_seqs.fna`. It also exports an Excel validation summary.

### Step 2 — [Run the PICRUSt2 Pipeline](R/notebooks/2_picrust2_pipeline.md) *(required)*

Runs `picrust2_pipeline.py` on the Step 1 BIOM and FASTA files. It produces predicted marker-copy-number, KO, EC, MetaCyc pathway, and NSTI outputs, records command provenance and software versions, and optionally decompresses downstream tables.

### Step 3 — [NSTI Quality Assessment](R/notebooks/3_nsti_quality_assessment.md) *(recommended)*

Summarizes PICRUSt2's per-ASV and abundance-weighted per-sample NSTI values, cross-checks the Step 2 filtering log, and calculates the percentage of reads excluded by the configured NSTI threshold. Its Tukey upper-fence flag (`Q3 + 1.5 × IQR` by default) is a descriptive screening rule, not a hypothesis test or p-value.

### Step 4 — [Functional Differential Abundance](R/notebooks/4_functional_differential_abundance_analysis.md) *(optional)*

Analyzes both MetaCyc and KEGG pathway predictions with LinDA and MaAsLin2. Each non-reference condition is compared independently with the condition marked `Reference = TRUE` in `metadata.tsv`. Benjamini-Hochberg adjustment is applied separately within each method, functional level, and comparison; `Method_Agreement` is descriptive and is not a combined p-value or additional test.

### Step 5 — [Taxon Contribution Analysis](R/notebooks/5_taxon_contribution.md) *(optional)*

Runs the stratified KO, EC, and MetaCyc prediction stages and analyzes both MetaCyc and KEGG pathway contributions. It produces documented workbooks and interactive condition-split heatmaps at genus and phylum levels. Taxa and pathways are ordered deterministically by contribution; no clustering or inferential test is performed in this step.

------------------------------------------------------------------------

## Column Dictionaries

Every Excel workbook ends with a `Column_Dictionary` sheet documenting each exported column in plain language. The dictionaries are generated from the actual output tables with [R/functions/build_column_dictionary_function.R](R/functions/build_column_dictionary_function.R), reducing the risk that workbook documentation drifts from the analysis code.

------------------------------------------------------------------------

## Project Structure

``` text
PICRUSt2_16S_Functional_Inference_Workflow/
├── PICRUSt2_16S_Functional_Inference_Workflow.Rproj
├── README.md
├── LICENSE
├── pages/
│   └── index.html
├── .github/workflows/
│   └── pages.yml
├── setup/
│   ├── install_R_dependencies.R
│   └── install_picrust2.sh
├── R/
│   ├── notebooks/
│   │   ├── 1_prepare_picrust2_inputs.Rmd
│   │   ├── 2_picrust2_pipeline.Rmd
│   │   ├── 3_nsti_quality_assessment.Rmd
│   │   ├── 4_functional_differential_abundance_analysis.Rmd
│   │   ├── 5_taxon_contribution.Rmd
│   │   └── *.md
│   ├── functions/
│   │   ├── add_sheet_to_excel_function.R
│   │   ├── build_column_dictionary_function.R
│   │   ├── decompress_output_files_function.R
│   │   ├── render_output_links_function.R
│   │   ├── render_output_tree_function.R
│   │   └── run_kegg_taxon_contribution_function.R
│   └── images/
│       └── 2_picrust2_flowchart.png
├── data/
│   └── README.md
├── example/
│   ├── README.md
│   ├── run_example.R
│   ├── data/
│   │   ├── asv_count_table.csv
│   │   ├── asv_sequences.csv
│   │   ├── asv_taxonomy.csv
│   │   └── metadata.tsv
│   └── reference_results/
│       ├── README.md
│       ├── reports/
│       │   └── *.html
│       └── <compact outputs from Steps 1–5>
└── results/
    └── <regenerable outputs from Steps 1–5>
```

------------------------------------------------------------------------

## References

### Core Functional-Inference Methods

- Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology* 38:685–688. [DOI:10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- Langille MGI, Zaneveld J, Caporaso JG, et al. (2013). Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. *Nature Biotechnology* 31:814–821. [DOI:10.1038/nbt.2676](https://doi.org/10.1038/nbt.2676)
- Yang C, Mai J, Cao X, Burberry A, Cominelli F, Zhang L (2023). ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics* 39(8):btad470. [DOI:10.1093/bioinformatics/btad470](https://doi.org/10.1093/bioinformatics/btad470)

### Differential-Abundance Methods

- Zhou H, He K, Chen J, Zhang X (2022). LinDA: linear models for differential abundance analysis of microbiome compositional data. *Genome Biology* 23:95. [DOI:10.1186/s13059-022-02655-5](https://doi.org/10.1186/s13059-022-02655-5)
- Mallick H, Rahnavard A, McIver LJ, et al. (2021). Multivariable association discovery in population-scale meta-omics studies. *PLoS Computational Biology* 17(11):e1009442. [DOI:10.1371/journal.pcbi.1009442](https://doi.org/10.1371/journal.pcbi.1009442)

### Databases and Formats

- Caspi R, Billington R, Keseler IM, et al. (2020). The MetaCyc database of metabolic pathways and enzymes — a 2019 update. *Nucleic Acids Research* 48(D1):D445–D453. [DOI:10.1093/nar/gkz862](https://doi.org/10.1093/nar/gkz862)
- Kanehisa M, Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28(1):27–30. [DOI:10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- McDonald D, Clemente JC, Kuczynski J, et al. (2012). The Biological Observation Matrix (BIOM) format. *GigaScience* 1:7. [DOI:10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7)

------------------------------------------------------------------------

## License

This project is released under the [MIT License](LICENSE).

------------------------------------------------------------------------

## Acknowledgments

- [PICRUSt2](https://github.com/picrust/picrust2) developers and the Huttenhower Lab for the functional-inference framework and documentation
- [ggpicrust2](https://github.com/cafferychen777/ggpicrust2), [MicrobiomeStat](https://github.com/cafferychen777/MicrobiomeStat), and [MaAsLin2](https://bioconductor.org/packages/release/bioc/html/Maaslin2.html) developers for the downstream functional-analysis methods and implementations
- [DADA2](https://benjjneb.github.io/dada2/) developers for the ASV inference framework that produces this workflow's expected inputs
- [Bioconductor](https://bioconductor.org/) and [CRAN](https://cran.r-project.org/) communities for the R packages used throughout the workflow
- [Claude AI](https://claude.ai/) for assistance with code and documentation

------------------------------------------------------------------------

**Author**: Amro Abbas - Generated with Claude AI assistance\
**Last Updated**: September 2026
