# SARS-CoV2_ONT_data
This repository contains the scripts that were used in Tombácz et al. 2022. High Temporal-Resolution Nanopore Sequencing Dataset of  SARS-CoV-2 and Host Cell RNAs.

The sequences can be found at:
https://www.ebi.ac.uk/ena/browser/view/PRJEB51064

The original *.bam* files (inital mapping) can be downloaded using the *PRJEB51064_SARS-CoV2_pass.download.sh* script, \
then the *fastqfrombam.sh* script can be used to generate the *.fastq* files (containing the viral reads only), which, in turn can mapped to the reference genome using the *minimap.sh* executable. \
Alternatively, *.fastq* files (containing every *passed* reads) can be downloaded from ENA and mapped using *minimap.sh*

The *gigasci.worklfow.R* contains the analysis of the reads, carried out in four steps: \
1.) Import the *.bam* files (itt will rename them to match the *sample_name* column in the "metadata.tsv" file) \
2.) Clustering the alignments into unique alignments, i.e. *Transcripts* \
3.) Detection of "Leaders" and "Trailers"; Differentiating betweem *genomic* and *sub-genomic RNAs* \
4.) Summing and statistics

The *gigasci.worklfow.R* file can be sourced, which will carry out the whole analyis. \
For each step and parameter there are comments to make the code understandable, it is suggested to read these before running the workflow.
The worklofw (by default) will save the anaylsis into an .RData image ("all.RData"). Finally, the *giga.plots.R* script should be sourced, which will load the data and then ultimately generate the plots.

For testing pourposes not every *.bam* file is required; the *PRJEB51064_SARS-CoV2_pass.download.sh* script can be modified (for example using comments (#) at the start of the lines) to download only the intented files to be analyzed (it is possible tough that this will result in errors, especially at the plotting part).


The scripts can be used with other *.bam* files, reference genomes and/or parameters as well to import, filter and anyalze of alignments in *R*. For example to dereplicate alignments into *Transcripts*. The settings and the "metadata.tsv" file has to be modified accordingly.

minimap2 has to be in the path for the mapping part ("minimap.sh"), the genome analysis toolkit (GATK) for the conversion of *.bam* files to *.fastq* ("fastqfrombam.sh"); and the required R packages must be installed before the scirpts can be run.


