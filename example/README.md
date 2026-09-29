# Bundled PICRUSt2 example

This directory contains a self-contained test profile for the complete workflow. It is deliberately separate from [`data/`](../data/) and the generated `results/` directory: running it cannot mix the example files with a user's study inputs, and adding another dataset to `data/` cannot alter the example run.

The example contains 10 samples and 1,468 representative 16S ASVs. It exercises input validation, BIOM/FASTA creation, the community-level PICRUSt2 pipeline, NSTI assessment, MetaCyc and KEGG differential-abundance analyses, and MetaCyc and KEGG taxon-contribution analyses.

## Contents

- `data/asv_count_table.csv`: sample-by-ASV integer read-count matrix for samples `S01`–`S10`.
- `data/asv_sequences.csv`: `ASV_ID` and representative nucleotide sequence for all 1,468 ASVs.
- `data/asv_taxonomy.csv`: DADA2-compatible taxonomy layout used to exercise genus- and phylum-level contribution summaries.
- `data/metadata.tsv`: five synthetic `Control` and five synthetic `Treatment` samples; `Control` is the reference condition.
- `reference_results/`: compact outputs and executed HTML reports from all five steps for viewing without running the workflow.
- `run_example.R`: executes the complete example without reading or modifying the normal `data/` or `results/` trees.

The filenames already match the workflow's documented input names, so no special example-only filename mapping is needed.

## Run the complete example

First install the repository's R dependencies and PICRUSt2 conda environment as described in the root [`README.md`](../README.md). Then run this command from the repository root:

``` bash
Rscript example/run_example.R
```

All five steps are required for this example. The runner validates the four input files before starting, sets isolated input and output paths, refuses extra arguments, removes only a previous `example/run_results/` run, and then renders the five notebooks in order. Both MetaCyc and KEGG analyses run in Steps 4 and 5. Generated files go only to the ignored `example/run_results/` directory.

The run uses the same PICRUSt2 and workflow settings as the normal notebooks. You may set the documented environment variables such as `PICRUST2_NUM_THREADS`, `PICRUST2_MAX_NSTI`, or `PICRUST2_CONDA_ENV_NAME` before starting. The runner always disables reuse so Steps 2 and 5 are executed from the bundled inputs rather than silently adopting an older run.

To inspect the outputs without installing or running anything, open the HTML files in [`reference_results/reports/`](reference_results/reports/) or use the repository's GitHub Pages site. Compact BIOM, FASTA, prediction, quality-assessment, differential-abundance, and taxon-contribution files are provided beside the reports.

## Interpretation limits

The `Control`/`Treatment` assignments are synthetic test metadata. Every row in `asv_taxonomy.csv` is explicitly labelled `Synthetic demo taxonomy in DADA2-compatible format; replace before biological interpretation`. These labels are used only to test joins, summaries, workbooks, and heatmaps. They must not be treated as sequence-based taxonomic assignments.

PICRUSt2 output is predicted functional potential inferred from marker-gene placement. It is not a direct measurement of genes, transcripts, proteins, metabolites, or pathway activity. The synthetic group labels also mean that differential-abundance results demonstrate execution and output structure, not a biological treatment effect.

## Data provenance and citation

The abundance and sequence tables were derived from human stool samples generated in our laboratory, then randomly reduced and anonymized for demonstration purposes.

For the analysis itself, cite:

> Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology*, 38, 685–688. <https://doi.org/10.1038/s41587-020-0548-6>

> Yang C, Mai J, Cao X, Burberry A, Cominelli F, Zhang L. (2023). ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics*, 39(8), btad470. <https://doi.org/10.1093/bioinformatics/btad470>

The root README lists the additional statistical-method, database, and BIOM-format references used by the workflow.
