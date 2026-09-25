# ETEC Plasmidome

This repository contains the R scripts and supporting data used for the analysis and visualization of the plasmidome of enterotoxigenic *Escherichia coli* (ETEC), with a particular focus on plasmids carrying colonisation factors (CFs) and other virulence-associated determinants.

The repository accompanies the associated manuscript and is intended to support transparency and reproducibility of the analyses and figures presented in the study.

## Repository contents

The repository contains R scripts used to generate the main and Extended Data figures, together with the corresponding input data.

### Main figures

* `Figure1.R` — Figure 1
* `Fig2A.R` — Figure 2A
* `Figure3BC.R` — Figure 3B–C
* `Figure4.R` — Figure 4
* `Figure5.R` — Figure 5

### Extended Data figures

* `Extended_Data_1.R` — Extended Data Figure 1
* `Extended_Data_2.R` — Extended Data Figure 2
* `Extended_Data_3.R` — Extended Data Figure 3
* `Extended_Data_4.R` — Extended Data Figure 4
* `Extended_Data_5.R` — Extended Data Figure 5
* `Extended_Data_6.R` — Extended Data Figure 6
* `Extended_Data_7.R` — Extended Data Figure 7
* `Extended_Data_8.R` — Extended Data Figure 8
* `Extended_Data_9.R` — Extended Data Figure 9

## Data

Input files used by the R scripts are provided in the repository where redistribution is permitted.

These include processed plasmid and isolate metadata, plasmid similarity/distance data, phylogenetic trees, colonisation factor datasets, and intermediate datasets used to generate the figures.

The repository also includes FASTA files containing sequences associated with the newly identified putative colonisation factors (NEW_CFs) described in the study. These sequence files are provided to facilitate further investigation, comparison, and reproducibility of the NEW_CF analyses.

The analyses make use of files in CSV, TSV, XLSX, Newick (`.nwk`), and FASTA formats.

## Running the scripts

The analyses were performed in R.

Required R packages are specified at the beginning of each script.

### File paths

The scripts retain the local file paths used during the original analysis. Before running a script on another system, these paths should be changed to point to the corresponding input files downloaded from this repository.

Apart from updating file paths as necessary, the scripts document the analysis and visualization procedures used to generate the corresponding figures.

## Reproducibility

Each figure script is intended to be run independently with its required input files.

Some final figures were subsequently adjusted for presentation after export (for example, where indicated in the corresponding R script). Therefore, minor graphical differences may exist between the direct R output and the final version presented in the manuscript.

## Citation

If you use the code or data from this repository, please cite the associated manuscript.

**Citation:**
To be added upon publication.


## Contact

For questions regarding the analyses or data in this repository, please contact the authors through the von Mentzer Lab.
