#!/bin/bash

NPROC=8
input='.fastq'
index='NC_045512.2.fasta'
outdir='.bam'

mkdir $outdir

for fastq in $(ls $input); do
	base=$(basename $input/$fastq .fastq)

	echo 'mapping ' $base ' to ' $index '...'
	minimap2 -ax splice -Y -C5 --MD -un \
		--secondary=no -g 30000 -G 30000 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 --splice-flank=no --max-chain-skip=40 --for-only \
		-t ${NPROC} ${index} ${input}/${fastq} > ${outdir}/${base}.sam
	
	samtools view ${outdir}/${base}.sam -b -@ {NPROC} -o ${outdir}/${base}.bam
	samtools sort  -@ ${NPROC} -o ${outdir}/${base}.bam ${outdir}/${base}'.bam'
	samtools index -@ ${NPROC} ${outdir}/${base}.bam
	
done
