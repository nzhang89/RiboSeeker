# metagene

#' Filter seqlevels
#'
#' @param regionGR A GRanges object of the target regions to calculate metagene. Note that
#' all regions must have equal width.
#' @param seqInfoKeep A Seqinfo object specifying which seqlevels in regionGR to keep.
#'
#' @return
#'
#' @importFrom BiocGenerics sort
#' @importFrom GenomeInfoDb sortSeqlevels seqlevels keepSeqlevels
#'
.filterSeqlevels = function(regionGR, seqInfoKeep) {
  # sort by seqlevels and filter seqlevels
  seqInfoKeep = sortSeqlevels(seqInfoKeep)

  regionGR = sortSeqlevels(regionGR)
  regionGR = sort(regionGR)

  if(!all(seqlevels(regionGR) %in% seqlevels(seqInfoKeep))) {
    message(sprintf('%s seqlevels in bamGR and regionGR not identical. Use common seqlevels.',
      .now()))

    seqlevelsKeep = seqlevels(seqInfoKeep)
    seqlevelsRegionGR = seqlevels(regionGR)
    seqlevelsKeep = intersect(seqlevelsKeep, seqlevelsRegionGR)

    if(length(seqlevelsKeep) == 0) {
      stop(paste('No common seqlevels between bamGR and regionGR.',
        'Make sure they are from the same assembly.'))
    }

    # prune seqlevels
    regionGR = keepSeqlevels(regionGR, seqlevelsKeep, pruning.mode='coarse')
    message(sprintf('%s The following seqlevels are dropped from regionGR: %s.', .now(),
      paste(seqlevelsRegionGR[!seqlevelsRegionGR %in% seqlevelsKeep], collapse=',')))
  }

  return(regionGR)
}

#' Get CDS start and end regions
#'
#' @param txdb A TxDb object of genome annotation. See GenomicFeatures package for more
#' details.
#' @param seqInfoKeep A Seqinfo object specifying which seqlevels in txdb to keep.
#' @param cdsStartUpstream A numeric variable indicating the width to use for the upstream
#' region of CDS start site (not including CDS start site).
#' @param cdsStartDownstream A numeric variable indicating the width to use for the downstream
#' region of CDS start site (including CDS start site).
#' @param cdsEndUpstream A numeric variable indicating the width to use for the upstream
#' region of CDS end site (including CDS end site).
#' @param cdsEndDownstream A numeric variable indicating the width to use for the downstream
#' region of CDS end site (not including CDS end site).
#'
#' @return a list of 3 elements, CDS start region, CDS end region, and transcript id (names of
#' the list are "start", "end", and "tx"). CDS start and end regions are both GRanges objects
#' with a mcol of "tx".
#'
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BiocGenerics sort as.data.frame
#' @importFrom IRanges trim
#' @importFrom GenomicFeatures cdsBy
#' @importFrom GenomicRanges makeGRangesFromDataFrame start end width
#' @importFrom dplyr %>% select group_by filter summarize
#'
.getCDSRegions = function(txdb, seqInfoKeep, cdsStartUpstream, cdsStartDownstream,
  cdsEndUpstream, cdsEndDownstream) {
  # CDS by transcript
  cdsGR = unlist(cdsBy(txdb, by='tx', use.names=TRUE))
  cdsGR$tx = names(cdsGR)
  names(cdsGR) = NULL

  # sort by seqlevels and filter seqlevels
  cdsGR = sortSeqlevels(cdsGR)
  cdsGR = sort(cdsGR)
  cdsGR = .filterSeqlevels(cdsGR, seqInfoKeep)
  cdsGRSeqInfo = seqinfo(cdsGR)

  # get cds start and end position
  cdsDF = as.data.frame(cdsGR)
  cdsDF = cdsDF %>% select(seqnames, start, end, strand, tx)
  cdsDF = cdsDF %>% group_by(tx) %>% filter(
    length(unique(seqnames)) == 1 &
    length(unique(strand)) == 1
  )
  cdsDF = cdsDF %>% group_by(tx) %>% summarize(
    seqnames = unique(as.character(seqnames)),
    start = min(start),
    end = max(end),
    strand = unique(strand)
  )

  # convert back to GRanges
  cdsGRUse = makeGRangesFromDataFrame(as.data.frame(cdsDF), seqinfo=cdsGRSeqInfo,
    keep.extra.columns=TRUE, starts.in.df.are.0based=FALSE)
  cdsGRUsePos = cdsGRUse[strand(cdsGRUse) == '+']
  cdsGRUseNeg = cdsGRUse[strand(cdsGRUse) == '-']

  # extract cds start and end regions
  cdsStartPos = cdsGRUsePos
  cdsStartNeg = cdsGRUseNeg
  start(cdsStartPos) = start(cdsGRUsePos) - cdsStartUpstream
  end(cdsStartPos) = start(cdsGRUsePos) + cdsStartDownstream - 1
  start(cdsStartNeg) = end(cdsGRUseNeg) - cdsStartDownstream + 1
  end(cdsStartNeg) = end(cdsGRUseNeg) + cdsStartUpstream
  cdsStart = trim(c(cdsStartPos, cdsStartNeg))

  cdsEndPos = cdsGRUsePos
  cdsEndNeg = cdsGRUseNeg
  start(cdsEndPos) = end(cdsGRUsePos) - cdsEndUpstream + 1
  end(cdsEndPos) = end(cdsGRUsePos) + cdsEndDownstream
  start(cdsEndNeg) = start(cdsGRUseNeg) - cdsEndDownstream
  end(cdsEndNeg) = start(cdsGRUseNeg) + cdsEndUpstream - 1
  cdsEnd = trim(c(cdsEndPos, cdsEndNeg))

  # keep regions not out of bounds
  idxKeep = ((width(cdsStart) == cdsStartUpstream + cdsStartDownstream) &
    (width(cdsEnd) == cdsEndUpstream + cdsEndDownstream)
  )
  cdsStart = cdsStart[idxKeep]
  cdsEnd = cdsEnd[idxKeep]
  cdsStartTx = cdsStart$tx
  cdsEndTx = cdsEnd$tx
  if(!all(cdsStartTx == cdsEndTx)) {
    stop('Inconsistent CDS start and end transcript order.')
  }

  return(list(
    start=cdsStart, end=cdsEnd, tx=cdsStartTx
  ))
}

#' Calculate metagene matrix for fixed length regions
#'
#' @param bamGR A GRanges object of aligned reads.
#' @param regionGR A GRanges object of the target regions to calculate metagene. Note that
#' all regions must have equal width.
#'
#' @return A matrix of metagene.
#'
#' @importFrom GenomeInfoDb seqnames seqlevels sortSeqlevels
#' @importFrom IRanges coverage Views viewApply
#' @importFrom GenomicRanges GRangesList
#' @importFrom BiocGenerics do.call cbind sort
#'
.calcMetagene = function(bamGR, regionGR) {
  # sort by seqlevels and filter seqlevels
  bamGR = sortSeqlevels(bamGR)
  bamGR = sort(bamGR)

  regionGR = sortSeqlevels(regionGR)
  regionGR = sort(regionGR)

  # calculate coverage
  cvg = coverage(bamGR)

  # split metagene regions by strand
  names(regionGR) = as.character(seq(1, length(regionGR)))
  regionGRPos = regionGR[strand(regionGR) == '+' | strand(regionGR) == '*']
  regionGRNeg = regionGR[strand(regionGR) == '-']
  # reduce seqnames
  seqlevels(regionGRPos) = unique(as.character(seqnames(regionGRPos)))
  seqlevels(regionGRNeg) = unique(as.character(seqnames(regionGRNeg)))
  regionGRLPos = GRangesList(split(regionGRPos, seqnames(regionGRPos)), compress=TRUE)
  regionGRLNeg = GRangesList(split(regionGRNeg, seqnames(regionGRNeg)), compress=TRUE)

  # metagene
  metagenePos = viewApply(Views(cvg[names(cvg) %in% names(regionGRLPos)], regionGRLPos),
    function(x) as.data.frame(x)$value)
  metageneNeg = viewApply(Views(cvg[names(cvg) %in% names(regionGRLNeg)], regionGRLNeg),
    function(x) as.data.frame(x)$value)

  # to matrix
  metagenePos = t(do.call(cbind, metagenePos))
  metageneNeg = t(do.call(cbind, metageneNeg))
  metageneNeg = metageneNeg[, ncol(metageneNeg):1] # reverse columns
  metagene = rbind(metagenePos, metageneNeg)

  rownames(metagene) = c(names(regionGRPos), names(regionGRNeg))
  metagene = metagene[names(regionGR), ]
  rownames(metagene) = NULL

  return(metagene)
}

#' Subset metagene matrices
#'
#' @description Keep transcripts with the highest expression in coding start and end regions.
#'
#' @param metagene A list of two matrices representing read counts in CDS start and end
#' regions. Each row is a transcript, each column is a position, and value represents read
#' counts in a transcript at a position.
#' @param nTx A numeric variable of the number of transcripts to keep.
#'
#' @return A filtered list of two matrices where the highest expressed transcripts (in coding
#' start and end regions) are kept.
#'
.subsetMetagene = function(metagene, nTx) {
  readsSum = rowSums(metagene[[1]]) + rowSums(metagene[[2]])
  names(readsSum) = rownames(metagene[[1]])
  readsSum = sort(readsSum, decreasing=TRUE) # sort in decreasing order

  # transcripts to keep
  txKeep = names(readsSum)[1:nTx]
  metageneCDSStart = metagene[[1]][which(rownames(metagene[[1]]) %in% txKeep), ]
  metageneCDSEnd = metagene[[2]][which(rownames(metagene[[2]]) %in% txKeep), ]

  # subset metagene
  metagene = list(start=metageneCDSStart, end=metageneCDSEnd)

  return(metagene)
}

#' Convert to row-wise relative frequency
#'
#' @param metagene A matrix of read counts or a list of two matrices representing read counts
#' in CDS start and end regions. Each row is a transcript, each column is a position, and
#' value represents read counts in a transcript at a position.
#' @param mode A numeric variable (1 or 2) indicating if metagene is a single matrix (mode
#' is 1) or a list of two matrices (mode is 2).
#' @param naVal A numeric variable to replace na with. (Default: 0).
#'
#' @return A matrix or a list of two matrices with values normalized as the row-wise relative
#' frequency. When mode is 2, the row-wise total counts are calculated from both matrices.
#'
.calcRelFreq = function(metagene, mode, naVal=0) {
  if(mode == 1) { # metagene is a single matrix
    totalCounts = rowSums(metagene)
    metageneRelFreq = metagene / totalCounts
    metageneRelFreq[is.na(metageneRelFreq)] = naVal
  } else if(mode == 2) { # metagene is a list of two matrices
    # total counts from both matrices
    totalCounts = rowSums(metagene[[1]]) + rowSums(metagene[[2]])
    metageneRelFreq1 = metagene[[1]] / totalCounts
    metageneRelFreq2 = metagene[[2]] / totalCounts
    metageneRelFreq1[is.na(metageneRelFreq1)] = naVal
    metageneRelFreq2[is.na(metageneRelFreq2)] = naVal
    metageneRelFreq = list(metageneRelFreq1, metageneRelFreq2)
    names(metageneRelFreq) = names(metagene)
  }

  return(metageneRelFreq)
}

#' Metagene matrix calculation
#'
#' @description Two modes are provided for metagene matrix calculation. If \code{regionGR} is
#' specified, metagene will be calculated for the regions given in \code{regionGR}. Otherwise,
#' \code{txdb} must be specified. If \code{txList} is specified, metagene will be calculated in
#' the CDS start and end regions for the transcripts specified. If \code{txList} is not
#' specified, the highest-expressed transcripts will be selected and metagene matrix will be
#' calculated in the CDS start and end regions. Also, the 5'-end position will be kept for the
#' reads.
#'
#' @param bam A \code{GAlignments} object of aligned reads.
#' (Required).
#' @param regionGR A \code{GRanges} object of the target regions to calculate metagene. Note
#' that all regions must have equal width. If \code{regionGR} is set, \code{txdb} will be
#' ignored. If \code{NULL}, \code{txdb} must be set. (Default: NULL).
#' @param txdb A \code{TxDb} object of genome annotation. See \code{GenomicFeatures} package
#' for more details. If \code{NULL}, \code{regionGR} must be set. (Default: NULL).
#' @param txList A character vector of transcript IDs. Note that the transcript IDs set here
#' should also be found in the \code{txdb}. (Default: NULL).
#' @param readLen A vector of read lengths to use (positive). If \code{NULL}, all
#' lengths will be kept. (Default: NULL).
#' @param relFreq A logical variable indicating if transcript wise relative frequency (
#' normalized by the total read counts in each metagene region) should be returned instead of
#' raw read counts. Note that if \code{txdb} is set but \code{regionGR} is \code{NULL},
#' metagene for both CDS start and end regions will be calculated. In this case, if setting
#' \code{relFreq} to \code{TRUE}, metagene will be normalized by the total read counts in both
#' CDS start and end region for each transcript. (Default: TRUE).
#' @param nTx A numeric variable of the number of transcripts to keep when \code{txdb} is
#' set but \code{txList} is \code{NULL}. The top \code{nTx} transcripts with the most reads in
#' coding start and end regions will be kept. (Default: 2000).
#' @param cdsStartUpstream A numeric variable indicating the width to use for the upstream
#' region of CDS start site (not including CDS start site) if \code{txdb} is set but
#' \code{regionGR} is \code{NULL}. (Default: 50).
#' @param cdsStartDownstream A numeric variable indicating the width to use for the downstream
#' region of CDS start site (including CDS start site) if \code{txdb} is set but
#' \code{regionGR} is \code{NULL}. (Default: 50).
#' @param cdsEndUpstream A numeric variable indicating the width to use for the upstream
#' region of CDS end site (including CDS end site) if \code{txdb} is set but
#' \code{regionGR} is \code{NULL}. (Default: 50).
#' @param cdsEndDownstream A numeric variable indicating the width to use for the downstream
#' region of CDS end site (not including CDS end site) if \code{txdb} is set but
#' \code{regionGR} is \code{NULL}. (Default: 50).
#'
#' @return A \code{list} containing three elements. The first element is metagne, either one
#' \code{matrix} if \code{regionGR} is set or a list of two \code{matrices} in CDS start and
#' end regions (with names of start and end). Each row is a transcript, each column is a
#' position, and a value represents read counts or relative frequency for a transcript at a
#' position. The second element is the \code{GRanges} for the metagene, either one range if
#' \code{regionGR} is set or a list of two ranges representing the CDS start and end regions.
#' This list also contains information of the transcripts selected and flanking sequence
#' lengths for CDS start and end. The third element is an internal variable indicating if
#' \code{regionGR} is specified or not (1 means \code{regionGR} is set and 2 means not set).
#'
#' @importFrom methods is as
#' @importFrom GenomeInfoDb  seqinfo
#' @importFrom GenomicAlignments qwidth
#' @importFrom IRanges resize
#' @importFrom GenomicRanges GRanges width
#'
#' @export
#'
calcMetagene = function(bam, regionGR=NULL, txdb=NULL, txList=NULL, readLen=NULL,
  relFreq=TRUE, nTx=2000, cdsStartUpstream=50, cdsStartDownstream=50,
  cdsEndUpstream=50, cdsEndDownstream=50) {
  # sanity check
  if(!is(bam, 'GAlignments')) {
    stop('bam must be a GAlignments object.')
  }
  if(!is.null(readLen)) {
    if(length(readLen) == 0) {
      stop('When readLen is not NULL, it must contain at least one element.')
    }
    readLen = round(unique(as.numeric(readLen)))
  }
  if(!is.logical(relFreq)) {
    stop('relFreq must be a logical variable.')
  }

  # metagene mode
  mode = 0
  if(!is.null(regionGR)) {
    if(!is(regionGR, 'GRanges')) {
      stop('If regionGR is set, it must be a GRanges object.')
    }
    if(length(unique(width(regionGR))) != 1) {
      stop('Not all regions of regionGR have equal widths.')
    }
    message(sprintf('%s regionGR is set and all regions have equal widths.', .now()))
    mode = 1
  } else {
    if(is.null(txdb)) {
      stop('When regionGR is NULL, txdb must be set.')
    }
    if(!is(txdb, 'TxDb')) {
      stop('txdb must be a TxDb object.')
    }

    txListStatus = 'not set'
    if(!is.null(txList)) {
      if(length(txList) == 0) {
        stop('When txList is not NULL, it must contain at least one element.')
      }
      txList = as.character(txList)
      txListStatus = 'set'
    }
    message(sprintf('%s regionGR is not set, txdb is set and txList is %s.', .now(),
      txListStatus))
    mode = 2
  }

  if(mode == 2) {
    if(!(length(nTx) == 1 & all(nTx > 0))) {
      stop('nTx must be a positive and numeric value.')
    }
    if(!(length(cdsStartUpstream) == 1 & all(cdsStartUpstream > 0))) {
      stop('cdsStartUpstream must be a positive and numeric value.')
    }
    if(!(length(cdsStartDownstream) == 1 & all(cdsStartDownstream > 0))) {
      stop('cdsStartDownstream must be a positive and numeric value.')
    }
    if(!(length(cdsEndUpstream) == 1 & all(cdsEndUpstream > 0))) {
      stop('cdsEndUpstream must be a positive and numeric value.')
    }
    if(!(length(cdsEndDownstream) == 1 & all(cdsEndDownstream > 0))) {
      stop('cdsEndDownstream must be a positive and numeric value.')
    }
    nTx = round(nTx)
    cdsStartUpstream = round(cdsStartUpstream)
    cdsStartDownstream = round(cdsStartDownstream)
    cdsEndUpstream = round(cdsEndUpstream)
    cdsEndDownstream = round(cdsEndDownstream)
  }

  # subset reads
  if(!is.null(readLen)) {
    bam = bam[qwidth(bam) %in% readLen] # select by read length
  }
  bamGR = as(bam, 'GRanges')
  bamGR = resize(bamGR, width=1)

  # metagene
  metageneObj = NULL
  if(mode == 1) {
    regionGR = .filterSeqlevels(regionGR, seqinfo(bamGR))
    metagene = .calcMetagene(bamGR, regionGR)
    metageneObj = list(metagene=metagene, region=regionGR, mode=mode)
  } else if(mode == 2) {
    # get CDS start and end regions
    cdsRegions = .getCDSRegions(txdb, bamGR, cdsStartUpstream, cdsStartDownstream,
      cdsEndUpstream, cdsEndDownstream)

    # filter or random sample
    txKeep = NULL
    subsetTx = FALSE
    if(!is.null(txList)) { # external tx list is set
      txKeep = which(cdsRegions$tx %in% txList)

      if(length(txKeep) == 0) {
        stop('No transcript can be found in txdb.')
      }
    } else { # extternal tx list is not set, use all transcripts for now
      if(length(cdsRegions$tx) < nTx) {
        message(sprintf('%s Number of CDS regions (n = %d) < nTx (%d) specified. Use all
          regions', .now(), length(cdsRegions$tx), nTx))
      }
      # txKeep = sample(cdsRegions$tx, min(length(cdsRegions$tx), nTx), replace=FALSE)
      nTx = min(length(cdsRegions$tx), nTx)
      txKeep = cdsRegions$tx
      subsetTx = TRUE
    }

    message(sprintf('%s %d transcripts are used to calculate metagene.', .now(),
      length(txKeep)))
    cdsRegionsStart = cdsRegions$start[which(cdsRegions$start$tx %in% txKeep)]
    cdsRegionsEnd = cdsRegions$end[which(cdsRegions$end$tx %in% txKeep)]
    cdsRegions = list(
      start=cdsRegionsStart,
      end=cdsRegionsEnd,
      tx=txKeep,
      startFlankLength=c(cdsStartUpstream, cdsStartDownstream),
      endFlankLength=c(cdsEndUpstream, cdsEndDownstream)
    )

    metageneCDSStart = .calcMetagene(bamGR, cdsRegions$start)
    metageneCDSEnd = .calcMetagene(bamGR, cdsRegions$end)
    rownames(metageneCDSStart) = cdsRegions$tx
    rownames(metageneCDSEnd) = cdsRegions$tx
    metagene = list(start=metageneCDSStart, end=metageneCDSEnd)

    if(subsetTx) {
      message(paste(sprintf('%s Keeping %d transcripts with the highest expression', .now(),
        nTx), 'in CDS start and end regions.'))
      metagene = .subsetMetagene(metagene, nTx)
      cdsRegions$tx = rownames(metagene[[1]])
      cdsRegions$start = cdsRegions$start[which(cdsRegions$start$tx %in% cdsRegions$tx)]
      cdsRegions$end = cdsRegions$end[which(cdsRegions$end$tx %in% cdsRegions$tx)]
    }

    metageneObj = list(metagene=metagene, region=cdsRegions, mode=mode)
  }

  # row-wise relative frequency
  if(relFreq) {
    metageneObj$metagene = .calcRelFreq(metageneObj$metagene, mode=metageneObj$mode)
  }

  return(metageneObj)
}
