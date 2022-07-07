#' This function finds the most upstream overlapping feature (usually CDS) from a genome annotation (feature.df) for
#' each exon and subsequently each transcript. The default is fully overlapping features, but 
#' the 'type', 'maxgap' and 'minoverlap' arguments can be used to change that.
#' 
#' @export
#'
#'

require(tidyr)
require(fuzzyjoin)

feature.OV.from.polyC.TR <- function(feature.df, TR.EX, EX.only=F, feature.colname='gene', antisense=F,
                                     type='within', by = c('seqnames', 'start', 'end'), maxgap=-1, minoverlap=1) {
  
  ### Giving dummy TR_IDs if there are none
  if (EX.only) { TR.EX$TR_ID <- TR.EX$EX_ID }
  
  ### Dropping non-mandatory columns
  TR.EX <- TR.EX[,c('TR_ID', 'EX_ID', "strand", by)]
  
  feature.df <- feature.df[, c(feature.colname, "strand", by)]
  # replace feature colname for convenience ->> have to re-rename at the end!
  colnames(feature.df)[1] <- 'feature'
  
  #### find overlaps of features and Transcripts
  tr.ov <- genome_join(feature.df[,],
                       TR.EX[,], 
                       type=type, by=by, maxgap=maxgap, minoverlap=minoverlap)
  colnames(tr.ov)[grep('.x', colnames(tr.ov))] <- gsub('.x', '.feature', colnames(tr.ov)[grep('.x', colnames(tr.ov))])
  colnames(tr.ov)[grep('.y', colnames(tr.ov))] <- gsub('.y', '.TR',      colnames(tr.ov)[grep('.y', colnames(tr.ov))])
  
  ### Sense overlaps
  tr.ov.s <- data.frame(NULL)
  try({
    tr.ov.s   <- tr.ov[tr.ov$strand.feature == tr.ov$strand.TR, ]
  })
  ## positive strand
  tr.ov.s.pos <- data.frame(NULL)
  try({
    tr.ov.s.pos <- tr.ov.s[tr.ov.s$strand.TR == '+', ]
    tr.ov.s.pos <- tr.ov.s.pos[order(tr.ov.s.pos$TR_ID, tr.ov.s.pos$start.feature), ]
  })
  ## negative strand
  tr.ov.s.neg <- data.frame(NULL)
  try({
    tr.ov.s.neg <- tr.ov.s[tr.ov.s$strand.TR == '-', ]
    tr.ov.s.neg <- tr.ov.s.neg[order(tr.ov.s.neg$TR_ID, tr.ov.s.neg$end.feature, decreasing = T), ]
  })
  tr.ov.s <- plyr::rbind.fill(tr.ov.s.pos, tr.ov.s.neg)
  
  ### Antisense overlaps
  if (antisense) {
    tr.ov.as <- data.frame(NULL)
    try({
      tr.ov.as   <- tr.ov[tr.ov$strand.TR != tr.ov$strand.feature, ]
    })
    ## positive strand
    tr.ov.as.pos <- data.frame(NULL)
    try({
      tr.ov.as.pos <- tr.ov.as[tr.ov.as$strand.TR == '+', ]
      tr.ov.as.pos <- tr.ov.as.pos[order(tr.ov.as.pos$TR_ID, tr.ov.as.pos$start.feature), ]
    })
    ## negative strand
    tr.ov.as.neg <- data.frame(NULL)
    try({
      tr.ov.as.neg <- tr.ov.as[tr.ov.as$strand.TR == '-',  ]
      tr.ov.as.neg <- tr.ov.as.neg[order(tr.ov.as.neg$TR_ID, tr.ov.as.neg$end.feature, decreasing = T), ]
    })
    tr.ov.as <- plyr::rbind.fill(tr.ov.as.pos, tr.ov.as.neg)
  
    ## Combine Sense and Antisense overlaps
    tr.ov <- plyr::rbind.fill(tr.ov.s, tr.ov.as) 
    
  } else {
    ## Keep only Sense overlaps
    tr.ov <- tr.ov.s
  }
  
  ## EXON-level features
  ex.ov.sum <- tr.ov %>% group_by(TR_ID, EX_ID, strand.TR) %>% 
    summarise(TR_ID, EX_ID, start.TR, end.TR, strand.TR, feature,  start.feature, end.feature, strand.feature, feature_n=seq_along(feature))
  levels <- paste0(feature.colname, '_', unique(ex.ov.sum$feature_n))
  ex.ov.sum$feature_n <- paste0(feature.colname, '_', ex.ov.sum$feature_n )
  ex.ov.sum$feature_n <- factor(ex.ov.sum$feature_n, levels = levels)
  ex.ov.sp  <- spread(ex.ov.sum[, c('TR_ID', 'EX_ID', 'start.TR', 'end.TR', 'strand.TR','feature','feature_n')], feature_n, feature)
  
  
  ## Transcript-level features
  tr.ov.sum <- tr.ov %>% group_by(TR_ID, strand.TR) %>% 
    summarise(TR_ID, feature,  start.feature, end.feature, strand.feature, feature_n=seq_along(feature))
  levels <- paste0(feature.colname, '_', unique(tr.ov.sum$feature_n))
  tr.ov.sum$feature_n <- paste0(feature.colname, '_', tr.ov.sum$feature_n )
  tr.ov.sum$feature_n <- factor(tr.ov.sum$feature_n, levels = levels)
  tr.ov.sp  <- spread(tr.ov.sum[, c('TR_ID', 'strand.TR','feature', 'feature_n')], feature_n, feature)
  
  return(list(ex.ov.sum, ex.ov.sp, tr.ov.sum, tr.ov.sp ))
}
