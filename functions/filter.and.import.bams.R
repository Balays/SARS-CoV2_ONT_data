
import.bams <- function( bamfiles, bamnames, write.filtered, bam.all=data.frame(NULL), 
                        filtering=c('all.reads.w.supp.ali', 'only.supp.ali.of.reads'), ...) {
  
  nbam <- length(bamfiles)
  
  for (i in seq(1, nbam) ) { 
    try({
    bamfile <- bamfiles[i]
    message('Start analyzing ', bamfile, '...')
    
    ## Filter alignments and write out filtered alignments
    if(write.filtered) {
      destination <- paste0(stringi::stri_replace_last_regex(bamfile, '.bam', ''), '.filtered.bam')
      if(any(!file.exists(destination), force.create)) {
        params  <- ScanBamParam(what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand"))
        bam.ori <- as.data.frame(scanBam(bamfile, param = params))
        if(filtering == 'all.reads.w.supp.ali') {
          message('This filters out every read that has (a) supplementary alignment(s), since these are most likely chimeric reads.')
          tokeep <- 
            !is.na(bam.ori$rname) &
            bam.ori$rname == seqnames.tofilt &
            bam.ori$mapq >= mapq.filt &
            !is.element(bam.ori$qname, dup(bam.ori$qname))
            
        } else if(filtering == 'only.supp.ali.of.reads') {
          message(' This filters out supplementary alignments, but leaves the other parts of the read.')
          tokeep <- 
            !is.na(bam.ori$rname) &
            bam.ori$rname == seqnames.tofilt &
            bam.ori$mapq >= mapq.filt &
            is.element(bam.ori$flag, flag.tokeep)
          
        } else {
          message('Do not filter reads, based on supplementary alignments.')
          tokeep <- 
            !is.na(bam.ori$rname) &
            bam.ori$rname == seqnames.tofilt &
            bam.ori$mapq >= mapq.filt
        }
          
        message(sum(!tokeep), ' reads were filtered out in ', bamnames[i], '!')
        message('Generating ', destination, ' ...')
        filterBam(BamFile(bamfile, yieldSize = length(tokeep)),  filter=tokeep, destination=destination)
        rm(bam.ori)
      }
      bamfile <- destination
    }
    ## Import filtered alignments into a data.frame
    bam <- ov.from.bam2(bamfile, ...
                        )
    bam$sample <- bamnames[i]
    bam.all <- plyr::rbind.fill(bam, bam.all)
  })
  }
  return(bam.all)
}  
  
  
  
