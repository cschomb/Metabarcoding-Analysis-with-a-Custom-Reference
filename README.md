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
|reference fasta||FASTA file used to build the BLAST databease (see above for requirements)|
|esv_fasta||Fasta file with the ESVs to be analysed (from the apscale output)|
|max_target_seqs||maximum number of target sequqences to be kept for further analysis from the blast|

### Run the analysis


## Resulting files
