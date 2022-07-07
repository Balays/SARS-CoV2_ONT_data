####  Mapped v9
mapv     <- 'mapv9'

#### RAW Mapped v9
pattern  <- '.merged.bam'

bamdir   <- 'I:/data/SARS-CoV2/mapped_v9/.bam'
bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]

filt <- F
tokeep <- 16:24
if(filt) { bamfiles <- bamfiles[tokeep] }

nbam <- length(bamfiles)
bamnames <- gsub('.*\\/', '', bamfiles)
bamnames <- gsub(pattern, '', bamnames)

pattern  <- '.raw'
bam.raw.mapv9 <- import.bams(bamfiles, bamnames, write.filtered=F, 
                             rm.gaps.in.aln = T, mapq.filt = mapq.filt,
                             flag.tokeep = NA, flag.tocrop = NA, seqnames.tofilt = seqnames.tofilt)
bam.raw.mapv9$pattern <- paste0(mapv, pattern)
bam.raw.mapv9.sum <- bam.raw.mapv9 %>% group_by(pattern, sample, strand) %>% summarise(count=n())


#### Mapped v8 dRNA
mapv     <- 'mapv8'

bam.raw.mapv8.dRNA <- import.bams(bamfiles = 'I:/data/SARS-CoV2/mapped_v8/dRNA.bam/dRNA.bam', bamnames = 'dRNA', 
                                  write.filtered=F, is.lortia = F, 
                                  rm.gaps.in.aln = T, mapq.filt = mapq.filt, 
                                  flag.tokeep = NA, flag.tocrop = NA, seqnames.tofilt = seqnames.tofilt)
bam.raw.mapv8.dRNA$pattern <- paste0(mapv, '.dRNA')
bam.raw.mapv8.dRNA.sum     <- bam.raw.mapv8.dRNA %>% group_by(pattern, sample, strand) %>% summarise(count=n())


#### Resequenced mapped v9 RAW
mapv     <- 'mapv9'
pattern  <- '.bam'

bamdir   <- 'I:/data/SARS-CoV2/rebasecall_reseq/.bam'
bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]

filt <- F
tokeep <- 16:24
if(filt) { bamfiles <- bamfiles[tokeep] }

nbam <- length(bamfiles)
bamnames <- gsub('.*\\/', '', bamfiles)
bamnames <- gsub(pattern, '', bamnames)

bam.reseq.raw.mapv9 <- import.bams(bamfiles, bamnames, write.filtered=F, 
                             rm.gaps.in.aln = T, mapq.filt = mapq.filt,
                             flag.tokeep = NA, flag.tocrop = NA, seqnames.tofilt = seqnames.tofilt)
pattern  <- '.reseq.raw'
bam.reseq.raw.mapv9$pattern <- paste0(mapv, pattern)
bam.reseq.raw.mapv9.sum <- bam.reseq.raw.mapv9 %>% group_by(pattern, sample, strand) %>% summarise(count=n())


#### SUM
bam.sum <- plyr::rbind.fill( bam.raw.mapv9.sum, bam.raw.mapv8.dRNA.sum, bam.reseq.raw.mapv9.sum )
bam.sum$pattern <- factor(bam.sum$pattern, levels = c('mapv9.raw', 'mapv9.reseq.raw', 'mapv8.dRNA'))
bam.sum$hpi  <- gsub('_.*', '', bam.sum$sample)
bam.sum$time <- as.integer(gsub('Hpi', '', bam.sum$hpi))
bam.sum$hpi  <- factor(bam.sum$hpi, levels = levels(metafilt$hpi))


if (make.plots) {
  gg1 <- ggplot(bam.sum[bam.sum$time <= 10
              &!is.na(bam.sum$time), ]) + 
    geom_col(aes(hpi, count, fill=pattern), position = position_dodge()) +
    scale_y_continuous("Read count", labels=function(x) format(x, big.mark = " ", scientific = FALSE) ) +
    scale_fill_manual(values = pal_aaas()(10)[c(1,2)]) +
    theme_ipsum() + 
    theme(plot.margin = unit(c(1,1,1,1), units = 'cm')) +
    facet_wrap(~hpi, scales = 'free_x', nrow = 1)
  
  gg2 <- ggplot(bam.sum[bam.sum$time >  10
                      & bam.sum$time <= 24
                      & !is.na(bam.sum$time), ]) + 
    geom_col(aes(hpi, count, fill=pattern), position = position_dodge()) +
    scale_y_continuous("Read count", labels=function(x) format(x, big.mark = " ", scientific = FALSE) ) +
    scale_fill_manual(values = pal_aaas()(10)[c(1,2)]) +
    theme_ipsum() + 
    theme(plot.margin = unit(c(1,1,1,1), units = 'cm')) +
    facet_wrap(~hpi, scales = 'free_x', nrow = 1)
  
  gg3 <- ggplot(bam.sum[bam.sum$time >  24
                      | is.na(bam.sum$time), ]) + 
    geom_col(aes(hpi, count, fill=pattern), position = position_dodge()) + 
    scale_y_continuous("Read count", labels=function(x) format(x, big.mark = " ", scientific = FALSE) ) +
    scale_fill_manual(values = pal_aaas()(10)[c(1,3)]) +
    theme_ipsum() +
    theme(plot.margin = unit(c(1,5,1,1), units = 'cm')) +
    facet_wrap(~hpi, scales = 'free_x', nrow = 1)
  
  #gg3 <- cowplot::plot_grid(gg3,NULL,nrow = 1,rel_widths = c(60,1))
  ggc <- cowplot::plot_grid(gg1,gg2,gg3, ncol = 1, align = 'v', axis = 'br')
  ggsave(paste0(fig.dir, '/libsizes.rebasecalled.jpg'), plot = ggc, height = 12, width = 10)
  
  
  
  gg <- ggplot(bam.sum) + 
    geom_col(aes(hpi, count, fill=pattern), position = position_dodge()) + 
    scale_fill_manual(values = pal_aaas()(10)[c(1:3)]) +
    scale_y_continuous("Read count", labels=function(x) format(x, big.mark = " ", scientific = FALSE) ) +
    theme_ipsum() + facet_wrap(~hpi, scales = 'free_x', nrow = 1)
  ggsave(paste0(fig.dir, '/libsizes.rebasecalled.2.jpg'), plot = gg, height = 6, width = 18)
  
}

