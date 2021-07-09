# BAM file operations

#' Load a BAM file
#'
#' @description Load a BAM file. This function is essentially a wrapper of
#' \code{\link[GenomicAlignments]{readGAlignments}}. Note that paried-end reads
#' will be treated as single-end.
#'
#' @param bamFile The path of a BAM file. Note that the BAM index file must also exist in the
#' same folder and have the same prefix. For example, if the BAM file is aln.bam, then the BAM
#' index file should be aln.bam.bai. (Required).
#' @param param A \code{\link[Rsamtools]{ScanBamParam}} object. See \code{Rsamtools} package
#' for more details. (Default: ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
#' isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE), mapqFilter=0)).
#'
#' @return A \code{\link[GenomicAlignments]{GAlignments}} object of the aligned reads.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom Rsamtools ScanBamParam scanBamFlag testPairedEndBam
#' @importFrom GenomicAlignments readGAlignments
#'
loadBam = function(bamFile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
  isSecondaryAlignment=FALSE, isNotPassingQualityControls=FALSE), mapqFilter=0)) {
  # sanity check
  if(!is(bamFile, 'character')) {
    stop('bamFile must be a character.')
  }
  if(!is(param, 'ScanBamParam') & !is.null(param)) {
    stop('param must be a GAlignments object or set to NULL.')
  }

  # check if bam file is single or paired end
  peCheck = testPairedEndBam(bamFile)
  if(peCheck) {
    message(sprintf('%s paired-end bam file will be treated as single-end', .now()))
  }

  # load bam file
  bam = readGAlignments(bamFile, param=param)

  return(bam)
}


#' Shift the 5'-end of reads of towards downstream position
#'
#' @param bam A GAlignments object of aligned reads of the same length.
#' @param readLen The read length to use.
#' @param shiftLen The shift size for this read length.
#'
#' @return A GAlignments object of shifted reads.
#'
#' @importFrom GenomicAlignments cigar cigarNarrow cigarQNarrow
#' @importFrom BiocGenerics strand
#'
.shift = function(bam, readLen, shiftLen) {
  cigars = cigar(bam)
  negStrands = as.character(strand(bam)) == '-'

  cigars = as.character(cigarNarrow(cigars)) # fix soft/hard clipping
  cigars = cigarQNarrow(cigars, start=ifelse(negStrands, 1, shiftLen+1),
    end=ifelse(negStrands, readLen-shiftLen, -1))
  bam@cigar = as.character(cigars)
  bam@start = bam@start + attributes(cigars)$rshift

  return(bam)
}

#' Shift the 5'-end of reads of towards downstream direction
#'
#' @description Shift the 5'-end of reads of towards downstream direction. Ribosome
#' profiling reads are known to have P-site offset. Reads on the positive strand
#' will be shifted to the right, and reads on the negative strand will be shifted
#' to the left.
#'
#' @param bam A \code{\link[GenomicAlignments]{GAlignments}} object of aligned reads.
#' (Required).
#' @param readLens A vector of read lengths to use (positive). If \code{NULL}, all
#' lengths will be kept. (Default: NULL).
#' @param shiftLens A vector of shift size (non-negative) for each read length.
#' If \code{NULL}, reads will not be shifted. If both \code{shiftLens} and \code{readLens}
#' are not \code{NULL}, they must have the same length. If \code{readlens} is \code{NULL}
#' and \code{shiftLens} is a single value, all reads will be shifted by this length.
#' (Default: NULL).
#' @param fiveEndOnly A logical variable indicating if only keeping the shifted 5'-end.
#' (Default: TRUE).
#'
#' @return If \code{fiveEndOnly} is \code{TRUE}, return a \code{\link[GenomicRanges]{GRanges}}
#' object of shifted reads. Otherwise, return a \code{\link[GenomicAlignments]{GAlignments}}
#' object of shifted reads.
#'
#' @export
#'
#' @importFrom methods is as
#' @importFrom GenomicAlignments qwidth
#' @importFrom IRanges resize
#' @importFrom GenomicRanges GRanges
#'
shiftReads = function(bam, readLens=NULL, shiftLens=NULL, fiveEndOnly=TRUE) {
  # sanity check
  if(!is(bam, 'GAlignments')) {
    stop('bam must be a GAlignments object.')
  }
  # resolve the cases when readLens or shiftLens is NULL
  if(is.null(readLens)) {
    readLens = sort(unique(qwidth(bam)))

    # sanity check
    if(length(shiftLens) > 1 & !is.null(shiftLens)) {
      stop('When readLens is NULL, shiftLens can only be a single value or NULL.')
    }

    if(!is.null(shiftLens)) {
      shiftLens = rep(shiftLens, length(readLens))
    }
  }
  if(is.null(shiftLens)) {
    shiftLens = rep(0, length(readLens))
  }
  if(!all(readLens > 0)) {
    stop('All values of readLens must be positive.')
  }
  if(!all(shiftLens >= 0)) {
    stop('All values of shiftLens must be non-negative.')
  }
  if(length(readLens) != length(shiftLens)) {
    stop('readLens and shiftLens do not match.')
  }
  if(!is.logical(fiveEndOnly)) {
    stop('fiveEndOnly must be a logical variable.')
  }

  # shift reads
  bamShifted = NULL
  qwidths = NULL
  for(i in seq(1, length(readLens))) {
    readLen = readLens[i]
    shiftLen = shiftLens[i]

    subBam = bam[qwidth(bam) == readLen] # select by read length
    subBamShifted = subBam

    if(length(subBam) > 0) {
      if(shiftLen != 0) {
        subBamShifted = .shift(subBam, readLen, shiftLen)
      }
    } else {
      message(sprintf('%s No reads have length of %d. Skip', .now(), readLen))
    }

    if(i == 1) {
      bamShifted = subBamShifted
      qwidths = qwidth(subBam)
    } else {
      bamShifted = c(bamShifted, subBamShifted)
      qwidths = c(qwidths, qwidth(subBam))
    }
  }

  # only keep 5-end position
  if(fiveEndOnly) {
    bamShifted = as(bamShifted, 'GRanges')
    bamShifted = resize(bamShifted, width=1)
    bamShifted$qwidth = qwidths
  }

  return(bamShifted)
}
