########## Feature annotation ##### 
#### Find features in exons

all.exdf <- feature.OV.from.polyC.TR(feature.df, feature.colname=feature.colname, ex.sp, EX.only = T)[[2]]

##
all.exdf <- as.data.frame(all.exdf[,-1])
#all.exdf$all.ORF <- apply(all.exdf[,grep('ORF', colnames(all.exdf))], 1, function(x) paste(unique(na.omit(as.character(x))), collapse = ';') )
#all.exdf <- paste.cols(all.exdf, grep('ORF', colnames(all.exdf),   value = T),  sep = ';', name = 'all.ORF')
all.exdf <- unite(all.exdf, col='all.ORF', starts_with('ORF'), na.rm = T, sep = ';', remove = F)
all.exdf$n_ORF  <- apply(all.exdf[,grep('ORF', colnames(all.exdf), value = T)], 1, function(x) sum(!is.na(x)))
all.exdf$EX.ORF <- all.exdf$ORF_1

##### merge exon data and TR data

#all.data <- merge(tr.gt, all.exdf, by.x=c("start",  "end", "strand"), by.y=c("start.TR",  "end.TR", "strand.TR"), all=T)
all.data <- merge(tr.gt, ex.sp[,c("EX_ID", by, "strand")], by.x=c("start",  "end", "strand"), by.y=c("start",  "end", "strand"), all=T)
all.data <- merge(all.data, all.exdf[,c("EX_ID", "EX.ORF", "all.ORF", "n_ORF")], by='EX_ID', all=T)

all.data <- merge(all.data, unique.data.frame(bam.filt[,c('qname', 'sample')]), by='qname', all=T)
all.data <- all.data[!is.na(all.data$qname), ]
#all.data$seqnames <- genome

#### Find features sequences in TRs
TR.EX    <- unique.data.frame(all.data[,c("EX_ID", "TR_ID",by,"strand")]) 

tr.data  <- feature.OV.from.polyC.TR(feature.df, feature.colname='ORF', TR.EX, EX.only = F)[[4]]
tr.data$tr.ORF <- tr.data$ORF_1

##### Data frame of exons with TR, read alignments and ORF info
all.data <- merge(all.data, unique.data.frame(tr.data[,c("TR_ID", "tr.ORF")]), by='TR_ID', all.x=T)



###### Find part of ORF1ab for genomic transcripts 

orf1ab.frag.thresh    <- 10
orf1ab.frag.df        <- CDS.df[CDS.df$protein == 'ORF1ab', ]
#orf1ab.frag.df       <- orf1ab.frag.df[orf1ab.frag.df$end == max(orf1ab.frag.df$end), ]
#orf1ab.frag.df$start <- orf1ab.frag.df$end - orf1ab.frag.thresh
#orf1ab.frag.df$ORF   <- 'ORF1ab'

orf1ab.frag.data  <- feature.OV.from.polyC.TR(feature.df = orf1ab.frag.df, feature.colname='ORF', TR.EX, EX.only = F,
                                              type='any', minoverlap=orf1ab.frag.thresh)[[4]]


###### Find full ORF1ab for genomic transcripts 

orf1ab.df        <- CDS.df[CDS.df$protein == 'ORF1ab', ]
orf1ab.df$end <- max(orf1ab.frag.df$end); orf1ab.df$start <- min(orf1ab.frag.df$start)
orf1ab.df <- unique.data.frame(select(orf1ab.df, -ORF))
orf1ab.df$ORF   <- 'ORF1ab'

orf1ab.full.data  <- feature.OV.from.polyC.TR(feature.df = orf1ab.df, feature.colname='ORF', TR.EX, EX.only = F, type='within')[[4]]
