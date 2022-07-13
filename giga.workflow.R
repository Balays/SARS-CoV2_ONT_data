##### Import libraries and functions
library(seqinr)
library(Rsamtools)
library(tidyverse)
library(stringi)
library(ggpubr)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)

source('functions/ov.from.bam2.R')
source('functions/get.best.aln.R')
source('functions/filter.and.import.bams.R')
bam.flags <- read.delim('functions/bam.flags.tsv')

##### Metadata
metadata <- read.delim('metadata.tsv', as.is = 12)
metafilt <- metadata#[grep('dRNA', metadata$sample, invert = T), ]
metafilt$hpi <- factor(metafilt$hpi, levels=c('Hpi1', 'Hpi2' ,'Hpi4', 'Hpi6', 'Hpi8', 'Hpi10',
                                              'Hpi12', 'Hpi14', 'Hpi16', 'Hpi18', 'Hpi20',
                                              'Hpi24', 'Hpi36', 'Hpi48' ,'Hpi72', 'Hpi96', 'dRNA'))

## which columns of the metadata table to be included in the downstream analysis?
## The first element of this vector should be sample ID, which should be the same as the .bam file's names.
metacols <- c('sample_name', 'hpi', 'Time')

#### Settings
## Reference genome
genome <- 'NC_045512.2' 
fasta  <- seqinr::read.fasta('NC_045512.2.fasta')
l_genome <- length(fasta[[1]])
## Alignment import
seqnames.tofilt <- genome ## Filter for alignments that mapped to this contig (viral genome)
flag.tokeep <- NA # This will not drop alignments in .bam files based on bamflag only
#c(0,16) ## this will leave primary alignments only and filter out supplementary and secondary alignments as well
mapq.filt   <- 1 ## Mapping quality threshold (in minimap2 60 is the highest)
write.filtered <- T ## Write out filtered .bam files
force.create <- T ## Overwrite?
filter.bams <- T ## filter that one read only has the best alignment
by <- c("seqnames", "start",  "end")

rename.files <- F ## if bam file's names are what the original (metadata "sample" column; downloaded from ENA), rename it to metadata 'sample_name' column

## Miscallenaous
luniq <- function(x) length(unique(x))
save.data <- NA #'all.RData'
make.plots <- T ## carry out plotting (run 'gigasci.plots.R')
fig.dir <- 'figures' ## save plots to this directory
by  <- c('seqnames', 'start', 'end')

## file name pattern for .bam files (this will be cropped from the file name to get sample name)
pattern  <- '.bam'

## load already saved data?
load <- F

if (load) {load(save.data)} else {
  
  
  ######  1. Import .bam files      ####
  bamdir   <- '.bam'
  bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
  bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]
  
  nbam <- length(bamfiles)
  bamnames <- gsub('.*\\/', '', bamfiles)
  bamnames <- gsub(pattern, '', bamnames)
  
  if(rename.files) {
    for (i in 1:nbam) {
      file   <- bamnames[i]
      sample <- metadata$sample_name[metadata$sample == file]
      source <- paste0(bamdir, '/', file,   '.bam')
      target <- paste0(bamdir, '/', sample, '.bam')
      file.rename(source, target)
    }
    bamdir   <- '.bam'
    bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
    bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]
    
    nbam <- length(bamfiles)
    bamnames <- gsub('.*\\/', '', bamfiles)
    bamnames <- gsub(pattern, '', bamnames)
    
  }
  
  
  
  bam.all <- import.bams(bamfiles, bamnames, write.filtered=F, 
                         rm.gaps.in.aln = T, mapq.filt = mapq.filt,
                         flag.tokeep = NA, flag.tocrop = NA, seqnames.tofilt = seqnames.tofilt)
  
  
  #### ####
}

######  2. Cluster reads   ####
## Cluster the reads on exact matching in the dataframe of alignments (bam.all)
source('cluster.reads.R')
 
## Merge 'Transcripts' with 'Exons'
all.data <- merge(tr.gt, ex.sp[,c("EX_ID", by, "strand")], by.x=c("start",  "end", "strand"), by.y=c("start",  "end", "strand"), all=T)
## Merge with reads
all.data <- merge(all.data, unique.data.frame(bam.filt[,c('qname', 'sample')]), by='qname', all=T)
all.data <- all.data[!is.na(all.data$qname), ]
## Merge with metadata
all.data <- merge(all.data, metafilt[,metacols], by.x='sample', by.y=metacols[1])

all.data$tr.ORF <- NA
TR.EX    <- unique.data.frame(all.data[,c("EX_ID", "TR_ID",by,"strand")]) 
#### ####


## save image
if(!is.na(save.data)) {
    save.image(save.data)
}



######  3. Detection of Leaders and Trailers, Differentiating of sub-genomic RNAs and genomic RNAs  ####
## Find part of ORF1ab for genomic transcripts ####
## Making a data.frame of the whole ORF1ab gene
CDS.df <-read.delim('NC_045512.2.ORFs.tsv')
orf1ab.frag.df        <- CDS.df[CDS.df$protein == 'ORF1ab', ]
orf1ab.frag.df$end    <- max(orf1ab.frag.df$end)
orf1ab.frag.df$start  <- min(orf1ab.frag.df$start)
orf1ab.frag.df$ORF    <- 'ORF1ab'
orf1ab.frag.df        <- unique.data.frame(orf1ab.frag.df)
## Find the overlaps between the ORF1ab and each read
orf1ab.frag.thresh    <- 10 ## This is the minimum overlap for each read with ORF1ab to be considered as of genomic origin.
ex.ov <- genome_join(ex.sp, orf1ab.frag.df, type='any', minoverlap=orf1ab.frag.thresh, by=by)
orf1ab.frag.data <- TR.EX[is.element(TR.EX$EX_ID, ex.ov$EX_ID), ]

## Parameters for leader and trailer identification
leader.thresh  <- c(85,55) ## Leaders will be considered as mapped regions ('exons') that end in this range
trailer.thresh <- 29749 ## Trailers will be considered if the mapped region ends after this position.

source('leaders.and.trailers.R')
#### ####
##


######  4. Sum data ####
##
genomic.tr.df <- as.data.frame(all.data %>% group_by(sample, hpi, Time, is.genomic) %>% summarise(freq=n()) %>% spread(is.genomic, freq, fill=0))
genomic.tr.df$ratio <-  genomic.tr.df$'TRUE' / genomic.tr.df$'FALSE' 
##
##### Data frame of reads with TR, exon and ORF info
tr.sp <- spread(all.data[,c("TR_ID","sample",metacols[2:length(metacols)], "qname","strand","exon.composition", 
                            "num_switches", "pos_nr","EX_ID","tr.ORF","TR.leader","TR.width","is.genomic", "is.subgenomic", "cluster")],
                pos_nr, EX_ID)   #c(1,2,3,6,8,9,10,17,18,19,21,24,25,26,27)
tr.sp <- tr.sp[,c(colnames(tr.sp)[1:14], cols)]

#### ####
##


## save image
if(!is.na(save.data)) {
  save.image(save.data)
}


if(make.plots) {
  source('giga.plots.R')
}

