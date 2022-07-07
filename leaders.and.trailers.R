##### Leaders  ##### 
### Which exons are leaders?
#leader.ex <- ex.sp$EX_ID[ex.sp$start <= leader.min ]
#leader.ex <- genome_join(TRS.L[,c('ORF', by, "strand")], 
#                         unique.data.frame(all.data[,c("EX_ID", by, "strand")]),
#                         by=by, type='within')
leader.ex <- ex.sp[ex.sp$end <= leader.thresh[1] & ex.sp$end >= leader.thresh[2], ]
leader.ex <- unique(leader.ex$EX_ID)

### Exons with leaders
all.data$is.leader <- F
all.data$is.leader[is.element(all.data$EX_ID, leader.ex)] <- T
plyr::count(all.data$is.leader)
### Transcripts with leaders
leader.tr <- unique.data.frame(all.data[, c("is.leader", "TR_ID")])
leader.tr <- leader.tr[leader.tr$is.leader == T, ]
all.data$TR.leader <- F
all.data$TR.leader[is.element(all.data$TR_ID, leader.tr$TR_ID)] <- T

plyr::count(all.data$TR.leader)

#tr.orf$TR.leader[is.element(tr.orf$TR_ID, leader.tr)] <- T

##### Data frame of reads with TR, exon and ORF info ##### 

tr.sp <- spread(all.data[,c("TR_ID","sample","qname","strand","exon.composition", "pos_nr","EX_ID","tr.ORF","hpi","Time","TR.leader")], 
                pos_nr, EX_ID)
cols <- colnames(tr.sp %>% select(starts_with('position')))
cols <- cols[order(cols)]
tr.sp <- tr.sp[,c(colnames(tr.sp)[1:9], cols)]

tr.sp$num_exon <- apply(tr.sp %>% select(starts_with('position')), 1, function(x) sum(!is.na(x)) )
tr.sp$num_switches <- tr.sp$num_exon - 1


### Transcirpt non-redundant info
tr.ex  <- unique.data.frame(dplyr::select(all.data, -c(sample, qname, hpi, Time)))
tr.ex$width <- abs(tr.ex$end - tr.ex$start)
##
tr.pos <- tr.ex %>% group_by(TR_ID) %>% summarise(TR.start=min(start), TR.end=max(end), TR.width=sum(width))
tr.uni <- unique.data.frame(dplyr::select(tr.sp, -c(sample, qname, hpi, Time)))
tr.uni <- merge(tr.pos, tr.uni, by='TR_ID')

### Trailers
tr.uni$TR.trailer <- F
tr.uni$TR.trailer[# tr.uni$TR.start >  85 &
  tr.uni$TR.end   >= trailer.thresh 
  # is.na(tr.uni$position_2)
  #& tr.uni$TR.width >= 5000
] <- T
plyr::count(tr.uni$TR.trailer)

### Cluster
tr.uni$cluster <- NA
tr.uni$cluster[tr.uni$TR.leader == T & tr.uni$TR.trailer == T] <- '5_3'
tr.uni$cluster[tr.uni$TR.leader == T & tr.uni$TR.trailer == F] <- '5_non3'
tr.uni$cluster[tr.uni$TR.leader == F & tr.uni$TR.trailer == T] <- 'non5_3'
tr.uni$cluster[tr.uni$TR.leader == F & tr.uni$TR.trailer == F] <- 'non5_non3'

tr.uni <- data.frame(dplyr::select(tr.uni, !cols), dplyr::select(tr.uni, cols))
plyr::count(tr.uni$cluster)

### Genomic reads
tr.uni$is.genomic <- F

## 1.) Very long TRs
tr.uni$is.genomic[ tr.uni$TR.width > (l_genome - 1000) ] <- T

## 2.) Full ORF1ab 
tr.uni$is.genomic[
       is.element(tr.uni$TR_ID,   orf1ab.full.data$TR_ID)
   #&  is.element(tr.uni$cluster, c('non5_non3', 'non5_3')) ->> same result when we filter out TRs with a JS after the TRS-L
       ] <- T

## 3.) 3-end of ORF1ab (10 nt) and no leader
#tr.uni$is.genomic[tr.uni$TR.leader == F & is.element(tr.uni$TR_ID, orf1ab.frag.data$TR_ID)] <- T
genomic.tr <- unique(all.data$TR_ID[is.element(all.data$EX_ID, ex.sp$EX_ID[ex.sp$start <= 21542 & ex.sp$end >= 29000 ])])
tr.uni$is.genomic[is.element(tr.uni$TR_ID, genomic.tr)] <- T
  
## 4.) Any 10-nt fragment of ORF1ab
tr.uni$is.genomic[is.element(tr.uni$TR_ID, orf1ab.frag.data$TR_ID)] <- T

### Subgenomic reads
tr.uni$is.subgenomic <- F
subgenomic.tr <- unique(all.data$TR_ID[is.element(all.data$EX_ID, ex.sp$EX_ID[ex.sp$start >= 21521 ])])
tr.uni$is.subgenomic[is.element(tr.uni$cluster, c('5_3', '5_non3')) &
                     is.element(tr.uni$TR_ID, subgenomic.tr)] <- T

table(tr.uni[,c("is.genomic", "is.subgenomic")] )


#tr.pos <- dplyr::select(tr.uni,c(TR_ID, TR.start, TR.end, is.genomic))
#tr.pos <- spread(all.data[,c(1,6,7,8,9,10,17,18,19,21)], pos_nr, pos)
tr.pos <- unique.data.frame(all.data[,c("TR_ID", "pos", "pos_nr")])
tr.pos <- spread(tr.pos, pos_nr, pos)
tr.pos <- tr.pos[,c(colnames(tr.pos)[1], cols)]
tr.pos <- merge(dplyr::select(tr.uni, !cols), tr.pos, by='TR_ID')

all.data <- dplyr::select(all.data, !any_of(c( 'TR.start', 'TR.end', 'TR.width','is.genomic', 'is.subgenomic', 'num_switches', 'TR.trailer', 'cluster')))
all.data <- merge(all.data, dplyr::select(tr.uni, c(TR_ID, TR.start, TR.end, TR.width, is.genomic, is.subgenomic, num_switches, TR.trailer, cluster)), by='TR_ID')

### Number of exons in TRs


######### ALL DATA DATAFRAME READY!