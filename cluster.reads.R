


##### Filter .bam files?
if (filter.bams) {
  bam.flags <- read.delim('D:/odrive/Google Drive/my.R.packages/bam.flags.tsv')
  bam.filt  <- get.best.aln(bam.all, bam.flags)
} else {
  bam.filt <- bam.all
}

plyr::count(bam.filt$strand)
plyr::count(bam.all$strand)


if(!is.na(save.data)) {
  save.image(save.data)
}
############################

#### Cluster reads into Transcripts ##### 
mtable <- bam.filt[,c("seqnames", "qname", 'start', 'end', 'strand')]
mtable$position <- paste0(mtable$start, '-', mtable$end, '::', mtable$strand)
tr.df <- tidyr::pivot_wider(mtable[,c("seqnames", "qname", "position")], names_from = seqnames, values_from = position)

#### omit full NA column
tr.df <- tr.df[,apply(tr.df, 2, function(x) {length(x) != sum(is.na(x))}) ]

ncol <- max(stri_count_regex(tr.df$NC_045512.2[], ',')) + 1
cols <- paste0('position_', c(1:ncol))
tr.df$NC_045512.2 <- gsub('c\\(', '', tr.df$NC_045512.2)
tr.df$NC_045512.2 <- gsub('\\)',  '', tr.df$NC_045512.2)
tr.df$NC_045512.2 <- gsub('\\"',  '', tr.df$NC_045512.2)

tr.df  <- separate(tr.df, NC_045512.2, sep=', ', into = cols, fill = "right", extra='merge')
tr.df  <- unite(tr.df, 'exon.composition', sep=';', cols, remove = F, na.rm = T)

tr.uni <- unique.data.frame(tr.df[,-1])
tr.uni <- data.frame(TR_ID=paste0('TR_', c(1:nrow(tr.uni))), tr.uni)[,1:2]

tr.df  <- merge(tr.uni, tr.df, by='exon.composition')

tr.gt  <- gather(tr.df, pos_nr, pos, -c(1:3))
tr.gt  <- tr.gt[!is.na(tr.gt$pos),]
nrow(tr.gt) == nrow(bam.filt)
## OK
#tr.gt <- merge(mtable, tr.gt, by='qname')

tr.gt$start  <- as.integer(gsub('-.*', '',  tr.gt$pos))
tr.gt$end    <- gsub('::.*', '', tr.gt$pos)
tr.gt$end    <- as.integer(gsub('.*-', '',  tr.gt$end))
tr.gt$strand <- gsub('.*::', '', tr.gt$pos)

### Exons
ex.df  <- bam.filt[,c("seqnames", "start", "end", "strand","sample")]
ex.df$pos <- paste0(ex.df$start, '-', ex.df$end)
ex.sum <- ex.df %>% group_by(seqnames, start, end, strand, sample) %>% summarise(count=n())
ex.sum <- ex.sum[order(ex.sum$strand, ex.sum$start), ]
ex.sp  <- spread(ex.sum, sample, count, fill = 0)
ex.sp  <- data.frame(EX_ID=paste0('ex_', seq(1, nrow(ex.sp))), ex.sp)
ex.sum <- gather(ex.sp, sample, count, -c(1:5))

