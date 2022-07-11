# SARS-CoV2_ONT_data
This repository contains the scripts that were used in Tombácz et al. 2022. High Temporal-Resolution Nanopore Sequencing Dataset of  SARS-CoV-2 and Host Cell RNAs.

The sequences can be found at:
https://www.ebi.ac.uk/ena/browser/view/PRJEB51064

The .fastq files should be dowlnoaded with *PRJEB51064_SARS-CoV2_pass.download.sh* and then mapped using the *minimap.sh* executable.

The *gigasci.worklfow.R* contains the analysis of the reads, carried out in four steps: \
1.) Import the *.bam* files (itt will rename them to match the *sample_name* column in the "metadata.tsv" file) \
2.) Clustering the alignments into unique alignments, i.e. *Transcripts* \
3.) Detection of "Leaders" and "Trailers"; Differentiating betweem *genomic* and *sub-genomic RNAs* \
4.) Summing and statistics

Finally, the script also contains: \
5.) Generating plots

The *gigasci.worklfow.R* file can be sourced, which will carry out the whole analyis. For each step and parameter there are comments to make the code understandable.
The required R packages must be installed before the scirpts can be run, obviously.


