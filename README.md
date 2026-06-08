# Metabarcoding-Analysis

## Introduction


## Requirements


## How to use

### Create a new project
```
EXAMPLE_PROJECT
|----01_blast
|    |---results
|        |--blast.out
|    |---tmp
|
|----02_analysis
|    |---results
|        |--best_blast_hits.xlsx
|        |--filtered_blast_results.xlsx
|        |--results.xlsx
|
|----03_gbif
|    |---results
|        |--gbif_occurrence_counts.xlsx
|        |--taxonomy_table.xlsx
|        |--results_with_gbif_presence
|
|----input
|
|----log.txt
|----settings.xlsx
```
### Settings file
#### blast-tab
|parameter|value|description|
|---|---|---|
|reference fasta|_filename_|FASTA file used to build the BLAST databease (see above for requirements). File in input folder or give path to file.|
|esv_fasta|_filename_|Fasta file with the ESVs to be analysed (from the apscale output). File in input folder or give path to file|
|max_target_seqs|_int_|maximum number of target sequqences to be kept for further analysis from the blast|

### analysis-tab
|parameter|value|description|
|---|---|---|
|esv_table|_filename_|Excel ESV table from the Apscale output. File in input folder or give path to file.|
|score_file|_filename_|Scoring file. Optional (if left empty analysis is run without considering the quality of the barcodes in the reference database. File in input or give path to file.|
|score_type|_str_|can be either _local_ (for species from Germany, _all_ (gives no extra scoring points for local species), or _column_name_ (name of a column in the scoring file with custom scoring).|
|length_threshold|_int_|Minimum alignment length to keep from the BLAST results|
|blast_results|_filename_|Optional override path to a BLAST output file. If empty, uses the one from _01_blast/results_.|

### gbif-tab
|parameter|value|description|
|---|---|---|
|input_results_xlsx|_filename_|Optional override. If blank, uses the one from _02_analysis_results_.|
|year_from|_int_|Count occurrences from this year onwards (inclusive).|
|lon|_float_|Longitude of circle center (WGS84)|
|lat|_float_|Latitude of circle center (WGS84)|
|radius_m|_int_|Radius for the occurrence search circle in meters|
|min_occurences|_int_|Minimum number of occurrences required for _occurs_in_area = Ture_|
|n_points|_int_|Number of polygon points to approximate the circle|
|workers_unauth|_int_|Worker threads when running the GBIF search WITHOUT credentials|
|workers_auth|_int_|Worker threads when running the GBIF search WITH credentials.|
|rps|_int_|number of requests per second for the GBIF API|
|gbif_user|_str_|GBIF username|
|gbif_pwd|_str_|GBIF password|


## Run the analysis


## Resulting files
