# MEtabarcoding-ANalysis against a CUstom REference (MEAN-CURE)

## Introduction
This script is used for the Analysis of Metabarcoding results from the Apscale Pipeline from a custom reference database. It assignes taxonomic information to the ESV sequences and can check the plausibility of these assignments by querying GBIF occurrence records.

<img width="966" height="504" alt="analysis_overview" src="https://github.com/user-attachments/assets/40d953df-40e7-47f0-b07c-8284c0cfa35f" />

### BLAST search
Performs a BLAST search of ESV sequences against a custom reference database (see Requirements for details). The maximum number of BLAST hits to keep can be set in the settings file. It also incorporates information about the length of the query coverage, which is used in the analysis step.

### Analysis
Filters the BLAST hits, so that the top 3% of the most identical records remain as candidates. These results are then calculated into a BLAST score:

$$BLAST score = (percent identity + query coverage)/2$$

The script will apply information thresholds according to the percent identity bracket of the BLAST hit:

```
100-97%  - all taxonomic information is kept
 97-94%  - taxonomic information about species is deleted
 94-90%  - taxonomic information about genus and species is deleted
  < 90%  - no taxonomic information is kept
```




### GBIF occurrence search

## Requirements
This script needs the output from a Apscale run (v4 or higher), and a reference FASTA file with structured headers, similar to the SINTAX format.

>[!IMPORTANT]
>```[Accession number];tax=p:[PHYLUM],c:[CLASS],o:[ORDER],f:[FAMILY],g:[GENUS],s:[SPECIES],gbiftax:[GBIF TAXONOMY],ncbitax:[NCBI TAXONOMY]```
>
>Example:
>
>```OR234967.1;tax=p:Tracheophyta,c:Magnoliopsida,o:Boraginales,f:Boraginaceae,g:Borago,s:Borago_officinalis,gbiftax:2926110,ncbitax:13363```


The taxonomy IDs for GBIF and NCBI are optional and can help to resolve taxonomic issues when querying sequences that come from NCBI originally against the GBIF taxonomy backbone. When left out, this will result in the corresponding columns to be empty in the results files. Prepepared reference fasta file can be found [here](https://github.com/cschomb/Primer-specific-Reference-Databases-for-Metabarcoding).

### Optional scoring file
The Primer specific reference databases also come with a scoring file, which reflect the quality of the reference sequences and the possibility to use them as a country-specific local database. Sequences score points according how well they can be distinguished from similar sequences and are optionally awarded points for local occurrence.



### Dependencies
This script was written for a UNIX Command Line Interface and requires `Python v3.12` and `ncbi-blast+ v2.12`

Further Python dependencies are:
`
numpy
pandas
requests
openpyxl
pyproj
tqdm
urllib3
`



## How to use
### Installation
```
pip install meancure
```

Upgrade:
```
pip install --upgrade meancure
```

### Create a new project
Initialise the folder structure with the `init` command
```
python3 analysis.py init [NAME]
```
This will results in the following folder structure:
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
|        |--filtered_blast_results.csv
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
Fill out the three tabs of the settings file with the required information.
#### blast-tab
|parameter|value|description|
|---|---|---|
|reference fasta|_filename_|FASTA file used to build the BLAST databease (see above for requirements). File in input folder or give path to file.|
|esv_fasta|_filename_|Fasta file with the ESVs to be analysed (from the apscale output). File in input folder or give path to file|
|max_target_seqs|_int_|maximum number of target sequqences to be kept for further analysis from the blast|

#### analysis-tab
|parameter|value|description|
|---|---|---|
|esv_table|_filename_|Excel ESV table from the Apscale output. File in input folder or give path to file.|
|score_file|_filename_|Scoring file. Optional (if left empty analysis is run without considering the quality of the barcodes in the reference database. File in input or give path to file.|
|score_type|_str_|can be either _local_ (for species from Germany, _all_ (gives no extra scoring points for local species), or _column_name_ (name of a column in the scoring file with custom scoring).|
|length_threshold|_int_|Minimum alignment length to keep from the BLAST results|
|blast_results|_filename_|Optional override path to a BLAST output file. If empty, uses the one from _01_blast/results_.|

#### gbif-tab
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


### Run the analysis
The script can be started with the following command
```
python3 analysis.py run --project [NAME] --steps [all]
```
When chosing `--steps all` the complete pipeline is run from BLAST to GBIF query.
It is also possible to run a subset of the steps. For instance if the BLAST is already done and the analysis and/or GBIF step should be repeated with different settings, the steps flag would be `--steps analysis,gbif`.
The new analysis is appended to the `log.txt` file and every run of the pipeline is distinguishable with the timestamp appended to the resulting files.



## Resulting files
### BLAST results
The BLAST produces a custom output tab-separated text with the following information in the columns:
|query id (from the ESV fasta file)|sequence id (Fasta header from the reference file)|percent identity|alignment length|number of mismatches|number of gaps|query start|query end|sequence start|sequence end|E-value|bitscore|query coverage|
|---|---|---|---|---|---|---|---|---|---|---|---|---|

### Analysis
#### filtered BLAST results
Retains the top 3% percent identity of each BLAST query, splits the information from the FASTA header of the reference sequences and adds information about the percent identity intervall of the remaining candidates.

#### best BLAST hits
Contains information of the best BLAST hit after the analysis. Includes information about the score of the reference sequence, the quality of the reference sequence, as well as a comment if multiple candidates have the same BLAST score and one of them is chosen because of its quality. Also adds a list of the other candidates with their corresponding percent identity.

#### results
Puts together the information from the ESV table (Apscale output) and the best BLAST hits. 
|hash|sequence|s_acc|s_acc ... s_species|gbif_taxid|ncbi_taxid|pident|evlaue|total_score|comment|quality|other_candidates|sum_of reads|samples|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|unique identifier from the ESV fasta file|Nucleotide sequence of the ESV|Accession number of the reference sequence|Taxonomic information of the reference sequence|GBIF taxid of the reference sequence|NCBI taxid of the reference sequence|Percent identity of the reference sequence to the ESV sequence|E-Value from the BLAST analysis|combined score of the reference sequence from blast score and scoring file|More than one reference sequence with the same BLAST score|quality of the reference sequence (from the scoring file)|List of other candidates from the BLAST, that were within 3% percent identity of the best hit|Read count of the ESV over all samples|Read counts for the individual samples|

### GBIF
#### taxonomy table
Table of taxonomies that are queried at GBIF. The taxonomies of the best BLAST hits are split up into _Species_, _Genus_, and _Family_; all of these levels are queried individually, while the list of other candidates are queried as _Species_.

#### gbif occurence counts
Results from the GBIF query with information about the original query that was send, the _gbif_match_type_ for the query, confidence for the match, number of occurences in the area queried, and a True / False column whether the query is accepted to occur in the area based on the threshold set in the setting file.

#### results with gbif presence
Combines the results file with the GBIF query results, by adding in new columns.
|...|quality|GBIF_presence_family|GBIF_presence_genus|GBIF_presence_species|remaining_other_candidates|sum_of_reads|...|
|---|---|---|---|---|---|---|---|
|...|...|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|Remaining other candidates, that occur in the queried area according to GBIF|...|...|

