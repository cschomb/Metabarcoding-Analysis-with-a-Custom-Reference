# MEtabarcoding-ANalysis against a CUstom REference (MEAN-CURE)

## Introduction
This script is used for the Analysis of Metabarcoding results from the Apscale Pipeline from a custom reference database. It assignes taxonomic information to the ESV sequences and can check the plausibility of these assignments by querying GBIF occurrence records.

<img width="966" height="504" alt="analysis_overview- Wiederhergestellt" src="https://github.com/user-attachments/assets/86b9c0f9-bdc1-431b-be2d-85e56df3525a" />


### BLAST search
Performs a BLAST search of ESV sequences against a custom reference database (see Requirements for details). The maximum number of BLAST hits to keep can be set in the settings file. It also incorporates information about the length of the query coverage, which is used in the analysis step.

### Analysis
Filters the BLAST hits, so that the top 3% of the most identical records remain as candidates. These results are then calculated into a BLAST score:

$$BLAST score = percent identity + (query coverage/5)$$

The script will apply information thresholds according to the percent identity bracket of the BLAST hit:

```
100-97%  - all taxonomic information is kept
 97-94%  - taxonomic information about species is deleted
 94-90%  - taxonomic information about genus and species is deleted
  < 90%  - no taxonomic information is kept
```
#### Optional scoring file
If there is additional information on the quality and reliability and / or the local occurrence of species in the reference sequences (as outlined in the example [here](https://github.com/cschomb/Primer-specific-Reference-Databases-for-Metabarcoding)). This information can also be incorporated into the analysis, in order to prefer more reliable and local references.


### GBIF occurrence search
The script also provides the possibility to perform a plausibility check against occurrences in GBIF. Taxa from the analysis can be compared against GBIF in an area and time period, which can be set and a threshold for the number of occurrences within the search parameters. The species, genera and families are then assigned either 'plausible' or 'implausible' status, depending on the outcome of the GBIF search.


## Dependencies and requirements

This script requires `Python v3.9x` and `ncbi-blast+ v2.12`

Further Python dependencies are:
`
numpy
pandas
requests
openpyxl
pyproj
tqdm
urllib3
streamlit>=1.37 (for the GUI version)
folium (for the GUI version)
streamlit-folium (for the GUI version)
`
### Reference FASTA
The reference FASTA is used to BLAST the ESVs against it. In order to apply the information thresholds during the analysis (as described above), the FASTA headers need to be structured similar to the SINTAX format:
>[!IMPORTANT]
>```[Accession number];tax=p:[PHYLUM],c:[CLASS],o:[ORDER],f:[FAMILY],g:[GENUS],s:[SPECIES],ncbitax:[NCBI TAXONOMY],gbiftax:[GBIF TAXONOMY],coltax[CATALOGUE OF LIFE TAXONOMY]```
>
>Example:
>
>```>AF091151.1;tax=p:Tracheophyta,c:Magnoliopsida,o:Boraginales,f:Boraginaceae,g:Borago,s:Borago_officinalis,ncbitax:13363,gbiftax:2926110,coltax:MJ92```

The taxonomy IDs for GBIF and NCBI are optional and can help to resolve taxonomic issues when querying sequences that originally come from NCBI against the GBIF taxonomy backbone. When left out, this will result in the corresponding columns to be empty in the results files. An example of a reference fasta file can be found [here](https://github.com/cschomb/Primer-specific-Reference-Databases-for-Metabarcoding).


### ESV table
The ESV table must contain a unique identifier (hash) for each ESV sequence, its actual sequence, and columns denoting the amount of reads per ESV and sample. Easiest to use is the output from the [Apscale](https://github.com/DominikBuchner/apscale) pipeline.

Example:
| hash                                    | sequence                                     | sample_1 | ... | sample_n |
|-----------------------------------------|----------------------------------------------|----------|-----|----------|
| d03850702e6fe2433d53fadddfce4f378add... | TTAAACTCAGCGGGTAGTCCCGCCTGACCTGGGGTCACATT... | 216      | 70  | 400      |
| 7ddcf83c8cfa7ed3a6c703a1b4b06cb77fa0... | TTAAACTCAGCGGGTAGTCCCGCCTGACCTGGGGTCACAAT... | 503      | 211 | 734      |
| ...                                     | ...                                          | ...      | ... | ...      |

### ESV FASTA
A FASTA file with the sequences of the ESVs. The FASTA headers must correspond to the unique identifiers (hash) in the ESV table. Here, as well the easiest to use is the ouptput from Apscale.


### Scoring file
To incorporate more information about the sequences in the reference database, a scoring file can be used (as described [here](https://github.com/cschomb/Primer-specific-Reference-Databases-for-Metabarcoding)). The file must have at least two columns: one named 'Accession' containing the accession numbers of the sequences in the reference database and a (named) column containing their score that is used in the analysis. 

Example (for a scoring file containing the mandatory columns, but also taxonomic information, a category how well the sequences can be distinguished from others in the database, whether the species occurs naturally in Germany (occurence_score) and the combined score (local_score) for use as a local database, as well as a general quality indicator of the reference sequence):

| Accession | Phylum | Class | Order | Family | Genus | Species | coltax | gbiftax | ncbitax | category | Score | occurence_score | local_score | quality |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AF046950.1 | Tracheophyta | Magnoliopsida | Asterales | Asteraceae | Bellis | Bellis_perennis | 41492 | 3117424 | 5WH44 | accepted | 4 | 5 | 9 | gold |
| AF069464.2 | Tracheophyta | Magnoliopsida | Lamiales | Plantaginaceae | Veronica | Veronica_annulata | 86152 | 7923758 | B443S | low_coverage | 2 | 0 | 2 | silver |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |




## How to use
### Command line version
#### Installation
```
pip install meancure
```

Upgrade:
```
pip install --upgrade meancure
```

#### Create a new project
Initialise the folder structure with the `init` command
```
meancure init [NAME]
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
|        |--occurrence_cache.json
|        |--taxonomy_table.xlsx
|        |--results_with_gbif_presence
|
|----input
|
|----log.txt
|----settings.xlsx
```
#### Settings file
Fill out the three sheets of the settings file with the required information.
##### Blast settings
|parameter|value|description|
|---|---|---|
|reference fasta|_filename_|FASTA file used to build the BLAST databease (see above for requirements). File in input folder or give path to file.|
|esv_fasta|_filename_|Fasta file with the ESVs to be analysed (from the apscale output). File in input folder or give path to file|
|max_target_seqs|_int_|maximum number of target sequqences to be kept for further analysis from the blast|

##### Analysis settings
|parameter|value|description|
|---|---|---|
|esv_table|_filename_|Excel ESV table from the Apscale output. File in input folder or give path to file.|
|score_file|_filename_|Scoring file. Optional (if left empty analysis is run without considering the quality of the barcodes in the reference database. File in input or give path to file.|
|score_type|_str_|can be either _local_ (for species from Germany, _all_ (gives no extra scoring points for local species), or _column_name_ (name of a column in the scoring file with custom scoring).|
|length_threshold|_int_|Minimum alignment length to keep from the BLAST results|
|blast_results|_filename_|Optional override path to a BLAST output file. If empty, uses the one from _01_blast/results_.|

##### GBIF settings
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


#### Run the analysis
The script can be started with the following command
```
meancure run --project [NAME] --steps [all]
```
When chosing `--steps all` the complete pipeline is run from BLAST to GBIF query.
It is also possible to run a subset of the steps. For instance if the BLAST is already done and the analysis and/or GBIF step should be repeated with different settings, the steps flag would be `--steps analysis,gbif`.
The new analysis is appended to the `log.txt` file and every run of the pipeline is distinguishable with the timestamp appended to the resulting files.

### GUI version
#### Installation
A version of this script with a graphical user interface is also available. It can either be started via command line:
```
pip install meancure-gui
```
The program can be started simply by entering its name in the promt:
```
meancure-gui
```
A browser window will open automatically and there new projects can be created by chosing a folder, Files needed for the analysis can be loaded and the analysis settings can be entered. (Info: when using WSL, the browser will not automatically open, simply copy the IP adress from the command line into your browser)

Upgrade:
```
pip install --upgrade meancure-gui
```

>[!NOTE]
> Windows users can download the zipped archive, extract it and start the App via the 'Meancure.vbs' file. Users need only provide a Python3.9x version and the BLAST executables (information on how to install are given as online help in the app, in case the executables are not found).


## Resulting files
### BLAST results
The BLAST produces a custom output tab-separated text with the following information in the columns:
|query id (from the ESV fasta file)|sequence id (Fasta header from the reference file)|percent identity|alignment length|number of mismatches|number of gaps|query start|query end|sequence start|sequence end|E-value|bitscore|query coverage|
|---|---|---|---|---|---|---|---|---|---|---|---|---|

### Analysis
#### Filtered BLAST results
Retains the BLAST hits within 3% percent identity of the top hit of each BLAST query, splits the information from the FASTA header of the reference sequences and adds information about the percent identity intervall of the remaining candidates.

#### Best BLAST hits
Contains information of the best BLAST hit after the analysis. Includes information about the score of the reference sequence, the quality of the reference sequence, as well as a comment if multiple candidates have the same BLAST score and one of them is chosen because of its quality from the optional scoring file. Also adds a list of the other candidates with their corresponding percent identities.

#### Results
Puts together the information from the ESV table (Apscale output) and the best BLAST hits. 
|hash|sequence|s_acc|s_phylum ... s_species|ncbitax|gbiftax|coltax|pident|evalue|total_score|comment|quality|other_candidates|sum_of reads|samples|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|unique identifier from the ESV fasta file|Nucleotide sequence of the ESV|Accession number of the reference sequence|Taxonomic information of the reference sequence|NCBI taxid of the reference sequence|GBIF taxid of the reference sequence|Catalogue of Life taxid of the reference sequence|Percent identity of the reference sequence to the ESV sequence|E-Value from the BLAST analysis|combined score of the reference sequence from blast score and scoring file|More than one reference sequence with the same BLAST score|quality of the reference sequence (from the scoring file)|List of other candidates from the BLAST, that were within 3% percent identity of the best hit|Read count of the ESV over all samples|Read counts for the individual samples|

### GBIF
#### Taxonomy table
Table of taxonomies that are queried at GBIF. The taxonomies of the best BLAST hits are split up into _Species_, _Genus_, and _Family_; all of these levels are queried individually, while the list of other candidates are queried as _Species_.

#### GBIF occurence counts
Results from the GBIF query, number of occurences in the area queried, and a True / False statement whether the query is accepted to occur in the area based on the threshold set in the setting file, written in json format. When repeating GBIF occurrence searches, already queried taxon keys are take from this cache file.

#### Results with gbif presence
Combines the results file with the GBIF query results, by adding in new columns.
|...|quality|GBIF_presence_family|GBIF_presence_genus|GBIF_presence_species|remaining_other_candidates|sum_of_reads|...|
|---|---|---|---|---|---|---|---|
|...|...|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|$${\color{green}plausible}$$ / $${\color{red}implausible}$$ depending on presence in GBIF query|Remaining other candidates, that occur in the queried area according to GBIF|...|...|

