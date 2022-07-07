# SARS-CoV2_ONT_data
This repository contains the scripts that were used in Tombácz et al. 2022. High Temporal-Resolution Nanopore Sequencing Dataset of  SARS-CoV-2 and Host Cell RNAs.

The sequences can be found at:
https://www.ebi.ac.uk/ena/browser/view/PRJEB51064

The .fastq files should be dowlnoaded and renamed to match the sample_name column in the "metadata.tsv" file. Alternatively, the viral .bam files can be downloaded as well and converted into .fastq files (this is possibly faster). Either way, the genereated .fastq files should be mapped using the minimap.sh executable.




