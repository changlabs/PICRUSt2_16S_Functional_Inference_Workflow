# Reference results

These are compact outputs from the completed local example run rendered on 3 September 2026 with R 4.5.2 and PICRUSt2 2.6.3. The source example contains 10 samples and 1,468 ASVs. PICRUSt2 excluded two poorly aligned sequences (`ASV286` and `ASV13801`) during placement; consult the executed Step 2 report for the recorded run details.

The `reports/` directory contains executed HTML reports for every step. Together, the compact step folders contain:

- the validated BIOM table, FASTA sequences, and Step 1 summary workbook;
- community-level KO, EC, MetaCyc, marker/NSTI tables and the Step 2 run-summary workbook;
- the NSTI quality-assessment workbook;
- separate MetaCyc and KEGG differential-abundance workbooks; and
- separate MetaCyc and KEGG taxon-contribution workbooks.

Large placement intermediates, combined per-ASV KO/EC predictions, stratified contribution tables, logs containing machine-specific paths, and duplicate standalone plot widgets are intentionally omitted. Links to bundled outputs remain clickable in the executed reports; omitted targets are displayed as annotated plain text instead of broken links. The reports embed the visualizations needed for inspection. Rerunning the example recreates the full output tree, including the omitted targets, under the ignored `example/run_results/` directory and never overwrites these committed reference results automatically.

The reports were generated before the inputs were moved from the normal ignored `data/` directory into `example/data/`. Their displayed paths were sanitized and their result links were adjusted for the committed reference layout; their numerical results were not changed.

The metadata groups and taxonomy are synthetic demonstration annotations, so the reports are execution examples rather than biological findings. The abundance and sequence tables were derived from human stool samples generated in our laboratory, then randomly reduced and anonymized for demonstration purposes.
