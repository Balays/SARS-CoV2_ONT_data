#!/bin/bash

NPROC=$(nproc)
input=.bam.ori
outdir=.fastq

mkdir $outdir


for bamfile in $(ls $input); do
	bam=$(basename $bamfile .bam)
	fastq=${bam}.fastq
	echo 'making ' $fastq ' from ' $bamfile '...'
	gatk SamToFastq -I ${input}/${bamfile} -F ${outdir}/${fastq}
done


