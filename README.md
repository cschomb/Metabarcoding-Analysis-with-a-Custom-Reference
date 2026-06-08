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
|reference fasta|**filename**|FASTA file used to build the BLAST databease (see above for requirements). File in input folder or give path to file.|
|esv_fasta|**filename**|Fasta file with the ESVs to be analysed (from the apscale output). File in input folder or give path to file|
|max_target_seqs|**int**|maximum number of target sequqences to be kept for further analysis from the blast|

### analysis-tab
|parameter|value|description|
|---|---|---|
|esv_table||Excel ESV table from the Apscale output. File in input folder or give path to file.|
||||
||||
||||
||||

### Run the analysis


## Resulting files
