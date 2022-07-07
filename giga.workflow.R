##### Import libraries and functions
library(seqinr)
library(misc)
library(Rsamtools)
library(tidyverse)
library(stringi)
library(ggpubr)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)

setwd('D:/ubuntu/my.repos/SARS-CoV2_ONT_data')

source('functions/ov.from.bam2.R')
source('functions/get.best.aln.R')
source('functions/feature.OV.from.polyC.TR.R')
source('functions/filter.and.import.bams.R')
#source('D:/odrive/Google Drive/my.R.packages/misc/R/ov.from.bam2.R')

by <- c("seqnames", "start",  "end")
#palette <-  c(pal_npg()(10), pal_aaas()(10))[c(1:10,14,18,15,13:11,16,17,20)]

##### Metadata
metadata <- read.delim('metadata.tsv', as.is = 12)
metafilt <- metadata#[grep('dRNA', metadata$sample, invert = T), ]
metafilt$hpi <- factor(metafilt$hpi, levels=c('Hpi1', 'Hpi2' ,'Hpi4', 'Hpi6', 'Hpi8', 'Hpi10',
                                              'Hpi12', 'Hpi14', 'Hpi16', 'Hpi18', 'Hpi20',
                                              'Hpi24', 'Hpi36', 'Hpi48' ,'Hpi72', 'Hpi96', 'dRNA'))

#### Settings
## genome
genome <- 'NC_045512.2' 
fasta  <- seqinr::read.fasta('NC_045512.2.fasta')
l_genome <- length(fasta[[1]])
## alignement import
seqnames.tofilt <- genome ## Filter for alignments that mapped to the viral genome
flag.tokeep <- NA #
#c(0,16) ## this will leave primary alignments only and filter out supplementary and secondary alignmnents as well
mapq.filt   <- 1 ## Mapping quality treshold (in minimap2 60 is the highest)
write.filtered <- T
force.create <- T
filter.bams <- T ## filter that one read only has the best alignment

## misc
save.data <- 'all.RData'
#make.plots <- T
fig.dir <- 'figures'
by  <- c('seqnames', 'start', 'end')

## file name pattern for .bam files (this will be cropped from the filename to get sample name)
pattern  <- '.bam'

## load already saved data?
load <- F

if (load) {load(save.data)} else {
  
  ######  1. Import .bam files      ####
  bamdir   <- '.bam'
  bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
  bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]
  
  filt <- F
  tokeep <- NULL
  if(filt) { bamfiles <- bamfiles[tokeep] }
  
  nbam <- length(bamfiles)
  bamnames <- gsub('.*\\/', '', bamfiles)
  bamnames <- gsub(pattern, '', bamnames)
  
  bam.all <- import.bams(bamfiles, bamnames, write.filtered=F, 
                         rm.gaps.in.aln = T, mapq.filt = mapq.filt,
                         flag.tokeep = NA, flag.tocrop = NA, seqnames.tofilt = seqnames.tofilt)
  
  
  bam.sum <- bam.all %>% group_by(sample, strand) %>% summarise(count=n())
  
  #### ####


  ######  2. Cluster reads      ####
  ## Cluster the reads on exact matching in the dataframe of alignments (bam.all)
  source('cluster.reads.R')
 
  #### ####
}

if(!is.na(save.data)) {
    save.image(save.data)
}


######  3. Feature annotation ####
## generate feature tables
CDS.df <-read.delim('NC_045512.2.ORFs.tsv')
feature.df <- CDS.df 
colnames(feature.df)[5:6] <- c('ORF', 'parent')
feature.colname <- 'ORF' 
feature.df <- feature.df[, c(feature.colname, "strand", by, 'parent')]
colnames(feature.df)[1] <- 'ORF'

## carry out feature annotation
source('feature.annotation.R')

all.data <- merge(all.data, metafilt[,c("sample_name", "hpi", "Time")], by.x='sample', by.y='sample_name')

#### ####

######  4. Detection of Leaders and Trailers, Differentiating of sub-genomic RNAs and genomic RNAs
leader.thresh  <- c(85,55)
trailer.thresh <- 29749
source('leaders.and.trailers.R')
## "ALL.DATA" DATAFRAME READY!
genomic.tr.df <- as.data.frame(all.data %>% group_by(sample, hpi, Time, is.genomic) %>% summarise(freq=n()) %>% spread(is.genomic, freq, fill=0))
genomic.tr.df$ratio <-  genomic.tr.df$'TRUE' / genomic.tr.df$'FALSE' 
#### ####

######  5. Sum data ####
source('sum.data.R')
reads <- unique.data.frame(all.data[,c("qname", "sample", 'hpi', 'Time')])
read.counts <- as.data.frame(reads %>% group_by(sample, hpi, Time) %>% summarise(read_count=n()))

#### ####

if(!is.na(save.data)) {
  save.image(save.data)
}


######  6. Generatin plots ####

## Figure 2. Ratio of sub-genomic / genomic transcripts

all.data$category <- NA
all.data$category[ all.data$is.genomic == T & all.data$is.subgenomic == F ] <- 'genomic'
all.data$category[ all.data$is.genomic == F & all.data$is.subgenomic == T ] <- 'sub-genomic'
all.data$category[ all.data$is.genomic == T & all.data$is.subgenomic == T ] <- 'sub-genomic'
all.data$category[is.na(all.data$category)] <- 'unclassified'

plot.data <- all.data %>% group_by(sample, category) %>% summarise(count=n()) %>% spread(category, count, fill=0)
#colnames(plot.data)[-1] <- c('non-genomic', 'genomic')
plot.data$ratio <- plot.data$'sub-genomic' / plot.data$genomic
plot.data <- merge(plot.data, metafilt[,c("sample_name", "h", "hpi", "Time", 'rep')], by.x='sample', by.y='sample_name')

plot.sum  <- plot.data[plot.data$hpi != 'dRNA', ] %>% 
  group_by(h, hpi, Time) %>% summarise(mean=mean(ratio), sd=sd(ratio), ymin=mean-sd, ymax=mean+sd )

gg.clust <- ggplot(plot.data[plot.data$hpi != 'dRNA', ]) + 
  geom_point(aes(x=Time, y=ratio)) + 
  geom_smooth(aes(x=Time, y=ratio)) +
  #geom_pointrange(aes(x=Time, y=mean, ymin=ymin, ymax=ymax, colour='mean'), size=0.25, data = plot.sum) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(name = 'Time (hours past infection)') +
  scale_y_continuous(name = 'ratio (sgRNA/gRNA)') +
  theme_bw() + theme(legend.position = 'none') 
gg.clust

ggsave('giga.Fig2.jpg', width = 20, height = 14)


## SuppFig S1 Violinplot of sub-genomic and genomic RNA lengths

tr.stats <- tr.sp[,1:14] #merge(tr.sp[,1:7], tr.uni[,1:4], by='TR_ID')
colnames(tr.stats)[11] <- "Mapped Length (nt)"
tr.stats$method <- 'cDNA'
tr.stats$method[tr.stats$sample == 'dRNA' ] <- 'dRNA'

tr.stats$category <- NA
tr.stats$category[ tr.stats$is.genomic == T & tr.stats$is.subgenomic == F ] <- 'genomic'
tr.stats$category[ tr.stats$is.genomic == F & tr.stats$is.subgenomic == T ] <- 'sub-genomic'
tr.stats$category[ tr.stats$is.genomic == T & tr.stats$is.subgenomic == T ] <- 'sub-genomic'
tr.stats$category[is.na(tr.stats$category)] <- 'unclassified'

tr.stats <- merge(tr.stats, metafilt[c('sample_name', 'h')], by.y='sample_name', by.x='sample')

ggv.cDNA <- ggviolin(tr.stats[tr.stats$method == 'cDNA', ], size = 0.25,
                     'h', "Mapped Length (nt)", color = 'h', #palette = 'npg', 
                     add=c('boxplot', 'median', 'mean_sd'), fill='lightgrey') + 
  scale_color_viridis_d() +
  coord_cartesian(ylim=c(0, 32000)) +
  theme_bw() + 
  theme(legend.position='none',
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()
  ) + 
  facet_grid(rows = vars(category), 
             cols = vars(method), scales = 'free_x')

ggv.dRNA <- ggviolin(tr.stats[tr.stats$method != 'cDNA', ], size = 0.25,
                     'method', "Mapped Length (nt)", color = 'method', palette = 'npg', 
                     add=c('boxplot', 'median', 'mean_sd'), fill='lightgrey') + 
  #scale_color_viridis_d() +
  coord_cartesian(ylim=c(0, 32000)) +
  theme_bw() + 
  theme(legend.position='none',
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
  ) + 
  facet_grid(rows = vars(category), 
             cols = vars(method), scales = 'free_x')
ggv <- cowplot::plot_grid(ggv.cDNA, ggv.dRNA, rel_widths = c(6,1), align = 'vh', axis = 'tblr')

ggsave('giga.SuppFig_S1.jpg', ggv, width = 24, height = 12)

#### ####
