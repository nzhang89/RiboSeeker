# ORF expression

#' Trim ORF
#'
#' @description Trim ORF from both ends. If the trimmed length is longer than one CDS, that
#' whole CDS will be removed and the trimming will continue on the next CDS using the
#' remaining trimming length. For example, given an ORF and the first CDS is 3bp, if trimStart
#' is 6bp, then the first CDS will be removed and the second CDS will be trimmed 3bp.
#'
#' @param orfGRL A GRangesList object of the ORFs.
#' @param trimStart A numeric variable indicating how many bases to trim for ORF start.
#' @param trimEnd A numeric variable indicating how many bases to trim for ORF end.
#'
#' @return a GRangesList object of trimmed ORFs.
#'
#' @importFrom BiocGenerics as.data.frame
#' @importFrom dplyr %>% filter group_by summarize inner_join mutate select
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom GenomeInfoDb seqinfo
#'
.trimORF = function(orfGRL, trimStart, trimEnd) {
  # standardize ORF
  orfGRL = .stdGRL(orfGRL, sameChr=TRUE, sameStrand=TRUE,
    orfLengthDivisibleBy3=FALSE)

  # work with a dataframe version of orfGRL as a personal preference
  orfDF = as.data.frame(orfGRL)
  colnames(orfDF) = make.names(colnames(orfDF), unique=TRUE)

  # split into positive and negative strand
  orfDFPos = orfDF %>% filter(strand == '+')
  orfDFNeg = orfDF %>% filter(strand == '-')

  # trim ORF boundary. We also check if the trimmed boundary is valid
  # for positive strand, ORF start and end are normal start and end columns
  # for negative strand, ORF start and end are actually end and start columns
  orfRangeDFPos = orfDFPos %>% group_by(group) %>% summarize(
    start = min(start) + trimStart,
    end = max(end) - trimEnd,
    valid1 = end > 0 & end >= start
  ) %>% as.data.frame()
  colnames(orfRangeDFPos) = c('group', 'orfStart', 'orfEnd', 'valid1')

  orfRangeDFNeg = orfDFNeg %>% group_by(group) %>% summarize(
    start = min(start) + trimEnd,
    end = max(end) - trimStart,
    valid1 = end > 0 & end >= start
  ) %>% as.data.frame()
  colnames(orfRangeDFNeg) = c('group', 'orfStart', 'orfEnd', 'valid1')

  orfRangeDF = rbind(orfRangeDFPos, orfRangeDFNeg) # row bind

  # add the trimmed boundary info to ORF dataframe
  orfDF = inner_join(orfDF, orfRangeDF, by='group')

  # trim each range if applicable
  orfDF = orfDF %>% mutate(
    trimmedStart = pmax(start, orfStart),
    trimmedEnd = pmin(end, orfEnd),
    valid2 = trimmedEnd >= trimmedStart
  )
  orfDF = orfDF %>% filter(valid1 & valid2) %>% select(group, group_name, seqnames,
    trimmedStart, trimmedEnd, strand)
  colnames(orfDF) = c('group', 'group_name', 'seqnames', 'start', 'end', 'strand')
  message(sprintf('%s %d ORFs are kept after trimming both ends.', .now(),
    length(unique(orfDF$group))))

  # rebuild GrangesList
  orfGRLTrimmedNames = unique(orfDF$group_name)
  names(orfGRLTrimmedNames) = unique(orfDF$group)

  orfGRLTrimmed = makeGRangesListFromDataFrame(orfDF, split.field='group',
    seqinfo=seqinfo(orfGRL), keep.extra.columns=TRUE)
  names(orfGRLTrimmed) = orfGRLTrimmedNames[names(orfGRLTrimmed)]

  return(orfGRLTrimmed)
}

#' Calculate RPKM
#'
#' @description Calculate read counts per kilo base per million reads (RPKM).
#'
#' @param bam A \code{GRanges} or \code{GAlignments} object of reads. Note that for Ribo-seq
#' data, the reads should be already size selected and shifted. Check function
#' \code{shiftReads} on how to shift reads. For RNA-seq data, there is no need to shift or size
#' select reads. Also, for each read, only the 5'-most position is used. (Required).
#' @param orfGRL A \code{GRangesList} object of ORFs. We recommend assigning a unique name to
#' each ORF using \code{names(orfGRL)}. In addition, the following modifications are also
#' applied: 1. If the names of orfGRL are NULL, rename each element as "orf_1", "orf_2", etc;
#' 2. Strands marked as "*" are replaced with "+"; 3. Remove elements with multiple chromosomes
#' or strands (one ORF is on multiple chromosomes or different strands); 4. Remove elements
#' where the ORF length is not divisible by 3; and 5. MOST IMPORTANTLY, if an ORF is on
#' positive strand, sort by coordinates (seqnames, start, end) in ascending order. Otherwise,
#' sort by coordinates (seqnames, end, start) in descending order. The purpose is to achieve
#' the same behavior as \code{cdsBy} function in \code{GenomicFeatures} package. (Required).
#' @param libSize A positive numeric variable indicating the library size of the reads. By
#' default, we use the number of reads in \code{bam} object specified. (Default: length(bam)).
#' @param trimStart A non-negative numeric variable indicating how many bases to trim for
#' ORF start. (Default: 6).
#' @param trimEnd A non-negative numeric variable indicating how many bases to trim for
#' ORF end. (Default: 6).
#' @param ignoreStrand A logical variable indicating if ignoring that reads and ORFs must
#' be on the same strand. (Default: TRUE).
#'
#' @return A \code{data.frame} with 4 columns, specified below: 1. Column 1 is
#' ORF ID (\code{orfId}, either user specified in \code{orfGRL} or internally generated);  2.
#' Column 2 is trimmed ORF length (\code{orfLenTrimmed}); 3. Column 3 is the read counts
#' (\code{countORF}) in the trimmed ORF region; Column 4 is the RPKM value (\code{rpkmORF}).
#'
#' @importFrom methods is as
#' @importFrom BiocGenerics sapply
#' @importFrom IRanges resize IntegerList
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicFeatures coverageByTranscript
#' @importFrom dplyr %>% mutate
#'
#' @export
#'
calcRPKM = function(bam, orfGRL, libSize=length(bam), trimStart=6, trimEnd=6,
  ignoreStrand=TRUE) {
  # sanity check
  if(!is(bam, 'GRanges') & !is(bam, 'GAlignments')) {
    stop('bam must be a GRanges or GAlignments object.')
  }
  if(!is(orfGRL, 'GRangesList')) {
    stop('orfGRL must be a GRangesList object.')
  }
  if(!(length(libSize) == 1 & all(is.numeric(libSize)) & all(libSize > 0))) {
    stop('libSize must be a positive numeric variable.')
  }
  if(!(length(trimStart) == 1 & all(is.numeric(trimStart)) & all(trimStart >= 0))) {
    stop('trimStart must be a non-negative numeric variable.')
  }
  if(!(length(trimEnd) == 1 & all(is.numeric(trimEnd)) & all(trimEnd >= 0))) {
    stop('trimEnd must be a non-negative numeric variable.')
  }
  if(!is.logical(ignoreStrand)) {
    stop('ignoreStrand must be a logical variable.')
  }

  # convert bam to GRanges and keep only the 5'-most position
  bamGR = as(bam, 'GRanges')
  bamGR = resize(bamGR, width=1)

  # trim ORF
  orfGRLTrimmed = .trimORF(orfGRL, trimStart=trimStart, trimEnd=trimEnd)

  # counts for each position for each ORF
  counts = coverageByTranscript(bamGR, orfGRLTrimmed, ignore.strand=ignoreStrand)
  counts = IntegerList(counts)

  # RPKM dataframe
  rpkmDF = data.frame(
    orfId = names(counts),
    orfLenTrimmed = sapply(counts, length),
    countORF = sapply(counts, sum)
  )

  rpkmDF = rpkmDF %>% mutate(
    rpkmORF = (1e+06 / libSize) * (1e+03 / orfLenTrimmed) * countORF
  )
  rownames(rpkmDF) = NULL

  return(rpkmDF)
}

#' Calculate translation efficiency
#'
#' @description Calculate translation efficiency given a Ribo- and RNA-seq sample, and a list
#' of ORF ranges. Basically, RPKM values for Ribo- and RNA-seq samples are calculated first,
#' and then the translation efficiency is calculated as the log2(riboRPKM + pseudoCount) -
#' log2(rnaRPKM + pseudoCount) with pseudoCount being a small value to prevent producing Inf.
#'
#' @param riboBam A \code{GRanges} or \code{GAlignments} object of reads. For Ribo-seq data,
#' the reads should be already size selected and shifted. Check function \code{shiftReads}
#' on how to shift reads. Also, for each read, only the 5'-most position is used. (Required).
#' @param rnaBam A \code{GRanges} or \code{GAlignments} object of reads. Note that for
#' RNA-seq data, there is no need to shift or size select reads. Also, for each read, only the
#' 5'-most position is used. (Required).
#' @param orfGRL A \code{GRangesList} object of ORFs. We recommend assigning a unique name to
#' each ORF using \code{names(orfGRL)}. In addition, the following modifications are also
#' applied: 1. If the names of orfGRL are NULL, rename each element as "orf_1", "orf_2", etc;
#' 2. Strands marked as "*" are replaced with "+"; 3. Remove elements with multiple chromosomes
#' or strands (one ORF is on multiple chromosomes or different strands); 4. Remove elements
#' where the ORF length is not divisible by 3; and 5. MOST IMPORTANTLY, if an ORF is on
#' positive strand, sort by coordinates (seqnames, start, end) in ascending order. Otherwise,
#' sort by coordinates (seqnames, end, start) in descending order. The purpose is to achieve
#' the same behavior as \code{cdsBy} function in \code{GenomicFeatures} package. (Required).
#' @param riboLibSize A positive numeric variable indicating the library size of the Ribo-seq
#' reads. By default, we use the number of reads in \code{riboBam} object specified.
#' (Default: length(riboBam)).
#' @param rnaLibSize A positive numeric variable indicating the library size of the RNA-seq
#' reads. By default, we use the number of reads in \code{rnaBam} object specified.
#' (Default: length(rnaBam)).
#' @param trimStart A non-negative numeric variable indicating how many bases to trim for
#' ORF start. (Default: 6).
#' @param trimEnd A non-negative numeric variable indicating how many bases to trim for
#' ORF end. (Default: 6).
#' @param ignoreStrand A logical variable indicating if ignoring that reads and ORFs must
#' be on the same strand. (Default: TRUE).
#' @param pseudoCount A non-negative numeric variable specifying a small value added to the
#' RPKMs calculated before taking log2 transformation in order to prevent producing Inf.
#' (Default: 1e-03).
#'
#' @return A \code{data.frame} with 7 columns, specified below: 1. Column 1 is ORF ID
#' (\code{orfId}, either user specified in \code{orfGRL} or internally generated);  2. Column
#' 2 is trimmed ORF length (\code{orfLenTrimmed}); 3. Column 3 and 4 are the read counts for
#' Ribo- and RNA-seq samples (\code{countRibo} and \code{countRNA}) in the trimmed ORF region,
#' respectively; Column 5 and 6 are the RPKM values for Ribo- and RNA-seq samples
#' (\code{rpkmRibo} and \code{rpkmRNA}), respectively; Column 7 is the calculated translation
#' efficiency (log2 transformed).
#'
#' @export
#'
#' @importFrom methods is as
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr %>% inner_join select mutate
#'
calcTransEff = function(riboBam, rnaBam, orfGRL, riboLibSize=length(riboBam),
  rnaLibSize=length(rnaBam), trimStart=6, trimEnd=6, ignoreStrand=TRUE,
  pseudoCount=1e-03) {
  # sanity check
  if(!is(riboBam, 'GRanges') & !is(riboBam, 'GAlignments')) {
    stop('riboBam must be a GRanges or GAlignments object.')
  }
  if(!is(rnaBam, 'GRanges') & !is(rnaBam, 'GAlignments')) {
    stop('rnaBam must be a GRanges or GAlignments object.')
  }
  if(!is(orfGRL, 'GRangesList')) {
    stop('orfGRL must be a GRangesList object.')
  }
  if(!(length(riboLibSize) == 1 & all(is.numeric(riboLibSize)) & all(riboLibSize > 0))) {
    stop('riboLibSize must be a positive numeric variable.')
  }
  if(!(length(rnaLibSize) == 1 & all(is.numeric(rnaLibSize)) & all(rnaLibSize > 0))) {
    stop('rnaLibSize must be a positive numeric variable.')
  }
  if(!(length(trimStart) == 1 & all(is.numeric(trimStart)) & all(trimStart >= 0))) {
    stop('trimStart must be a non-negative numeric variable.')
  }
  if(!(length(trimEnd) == 1 & all(is.numeric(trimEnd)) & all(trimEnd >= 0))) {
    stop('trimEnd must be a non-negative numeric variable.')
  }
  if(!is.logical(ignoreStrand)) {
    stop('ignoreStrand must be a logical variable.')
  }
  if(!(length(pseudoCount) == 1 & all(is.numeric(pseudoCount)) & all(pseudoCount >= 0))) {
    stop('pseudoCount must be a non-negative numeric variable.')
  }

  # calculate RPKM
  riboRPKMDF = calcRPKM(riboBam, orfGRL, libSize=riboLibSize, trimStart=trimStart,
    trimEnd=trimEnd, ignoreStrand=ignoreStrand)
  rnaRPKMDF = calcRPKM(rnaBam, orfGRL, libSize=rnaLibSize, trimStart=trimStart,
    trimEnd=trimEnd, ignoreStrand=ignoreStrand)
  colnames(riboRPKMDF) = c('orfId', 'orfLenTrimmed', 'countRibo', 'rpkmRibo')
  colnames(rnaRPKMDF) = c('orfId', 'orfLenTrimmed', 'countRNA', 'rpkmRNA')

  # calculate translation efficiency
  teDF = inner_join(riboRPKMDF, rnaRPKMDF, by=c('orfId', 'orfLenTrimmed'))
  teDF = teDF %>% select(orfId, orfLenTrimmed, countRibo, countRNA, rpkmRibo, rpkmRNA)
  teDF = teDF %>% mutate(
    transEff = log2(pseudoCount + rpkmRibo) - log2(pseudoCount + rpkmRNA)
  )

  return(teDF)
}
