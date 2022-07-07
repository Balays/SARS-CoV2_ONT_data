#' Filter reads for the highest mapq, and/or filter out supplementary and/or secondary alignments.
#' Requires data.frame containing bam flag information
#' @export
#'
get.best.aln <- function(bam, bam.flags, best.mapq=T, rm.supplementary=T, rm.secondary=T, keep.chim=F) {
  
  require(tidyr)
  require(dplyr)
  require(misc)
  require(stringr)
  
  bam.seq <- bam
  luniq.qname <- length(unique(bam.seq$qname))
  ### get best mapq
  if (best.mapq) {
    bam.mapq <- bam.seq %>% group_by(qname) %>% summarise(qname, mapq.max=max(mapq))
    bam.mapq <- unique.data.frame(bam.mapq)
    bam.seq  <- merge(bam.seq, bam.mapq, by.x = c('qname', 'mapq'), by.y = c('qname', 'mapq.max'))
    
    if (luniq.qname != length(unique(bam.seq$qname))) {
      message('When filtering for best mapq, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {
      message('Reads were filtered, based on mapq for the best alignment only!')
    }
    
  }
  
  if (rm.secondary) {
    sec.flag   <- 256
    sec.flags  <- c(sec.flag, bam.flags$Decimal + sec.flag)
    bam.nosec  <- bam.seq[!is.element(bam.seq$flag, c(sec.flags)), ]
    bam.seq    <- bam.nosec
    if (luniq.qname != length(unique(bam.seq$qname))) {
      noprim <- setdiff(bam$qname, bam.nosupp$qname)
      message('When filtering out secondary alignments, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {
      noprim <- NULL
      message('Secondary alignments were fitered out!')
    }
    
  }
  
  if (rm.supplementary) {
    supp.flag   <- 2048
    supp.flags  <- c(supp.flag, bam.flags$Decimal + supp.flag)
    bam.nosupp  <- bam.seq[!is.element(bam.seq$flag, c(supp.flags)), ]
    bam.seq     <- bam.nosupp
    if (luniq.qname != length(unique(bam.seq$qname))) {
      nolinear <- setdiff(bam$qname, bam.nosupp$qname)
      message('When filtering out supplementary alignments, ', luniq.qname - length(unique(bam.seq$qname)), ' reads were dropped!')
    } else {
      nolinear <- NULL
      message('Supplementary alignments were fitered out!')
    } 
  }
  
  chim <- c(noprim, nolinear)
  if(!is.null(chim)) {
    bam.chim <- bam[is.element(bam$qname, chim), ]
    
    if (keep.chim) {
      message('Keeping the alignment, which has highest mapq and and greatest width')
      #qname <- '25bb9865-f2e5-44f8-b88e-8956af95107f'
      best.reads <- data.frame(NULL)
      for (qname in unique(bam.chim$qname)) {
        read.df <- bam[bam$qname == qname, ]
        read.df <- read.df[order(read.df$mapq, read.df$width, decreasing = T), ]
        read.df <- read.df[1,]
        best.reads <- plyr::rbind.fill(read.df, best.reads)
      } 
      bam.chim <- best.reads
    } else {
      message('Keeping primary alignments in those cases where the secondary or supplementary had a higher mapq.')
      bam.chim <- bam.chim[bam.chim$flag == 0, ]
    }
    bam.seq <- plyr::rbind.fill(bam.chim, bam.seq)
  }  
  return(bam.seq)
}


