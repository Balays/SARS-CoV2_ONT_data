#' Import alignments from a .BAM file.
#' Optional filtering can be applied for the alignments, based on bam flags, mapping quality, and reference sequences.
#' LoRTIA ouptuts can be imported, along with LoRTIA tags for further filtering.
#' Gaps (Ns) in the alignments can be removed into distinct alignments for the same read using the CIGAR values.
#'
#' @export

ov.from.bam2 <- function (bamfile, flag.tokeep = c(0, 16), flag.tocrop = c(4),
          seqnames.tofilt = NA, mapq.filt = NA, rm.gaps.in.aln = F, primes = T,
          lortia.tags = c("l3", "l5", "r3", "r5"), is.lortia = F, rm.false.exons = F, rm.non.correct = F) {
  require(GenomicAlignments)
  require(Rsamtools)
  #bamfile <- 'I:/data/SARS-CoV2/mapped_v8/.bam.fastq.pych.bam/Hpi10_A.merged_pychopped.bam'

  if (is.lortia) {
    params <- ScanBamParam(what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand"), tag = lortia.tags)
    bam <- as.data.frame(scanBam(bamfile, param = params))
    bam$tags <- paste(bam[1, paste0("tag.", lortia.tags)], collapse = ",")
    bam$tags <- apply(as.data.frame(
      apply(bam[, paste0("tag.", lortia.tags)], 2, function(x) gsub(".*,", "", x))),
                      1, function(x) paste0(unique(x), collapse = ","))
    if (rm.non.correct) {
      bam <- bam[grepl("correct", bam$tags), ]
    }
    if (rm.false.exons) {
      bam <- bam[!grepl("false", bam$tags), ]
    }
  } else {
    params <- ScanBamParam(what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand"))
    bam <- as.data.frame(scanBam(bamfile, param = params))
  }

  if (!is.na(seqnames.tofilt)) {
    bam <- bam[!is.na(bam$rname), ]
    bam <- bam[bam$rname == seqnames.tofilt, ]
  }
  if (!is.na(mapq.filt)) {
    bam <- bam[bam$mapq >= mapq.filt, ]
  }
  if (!is.na(flag.tokeep)) {
    bam <- bam[is.element(bam$flag, flag.tokeep), ]
  }
  if (!is.na(flag.tocrop)) {
    bam <- bam[!is.element(bam$flag, flag.tocrop), ]
  }
  if (rm.gaps.in.aln) {
    bam$group <- seq(1, nrow(bam))
    aln  <- as.data.frame(extractAlignmentRangesOnReference(bam$cigar, pos = bam$pos, f = NULL))
    stopifnot(luniq(aln$group) == nrow(bam))
  } else {
    bam$group <- seq(1, nrow(bam))
    aln <- as.data.frame(cigarRangesAlongReferenceSpace(bam$cigar, pos = bam$pos, f = NULL, reduce.ranges=T))
    stopifnot(nrow(aln) == nrow(bam))
  }
  aln <- dplyr::select(aln, -group_name)
  aln <- merge(aln, bam, by = "group")[, -1]
  bam <- aln
  if (is.lortia) {
    bam <- bam[, c("rname", "strand", "cigar", "qwidth", "start", "end", "width", "qname", "flag", "pos", "mapq",
                   paste0("tag.", lortia.tags), "tags")]
  } else {
    bam <- bam[, c("rname", "strand", "cigar", "qwidth", "start", "end", "width", "qname", "flag", "pos", "mapq")]
  }

  colnames(bam)[1] <- "seqnames"
  if (primes) {
    bam$prime5 <- bam$start
    bam$prime3 <- bam$end
    bam$prime5[bam$strand == "-"] <- bam$end[bam$strand == "-"]
    bam$prime3[bam$strand == "-"] <- bam$start[bam$strand == "-"]
  }
  return(bam)
}
