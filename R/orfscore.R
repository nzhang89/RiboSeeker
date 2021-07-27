# ORFScore

#' Standadize ORF GrangesList
#'
#' @description We modify the ORF GrangesList in the following ways: 1. If the names of the
#' orfGRL are NULL, rename each element as "orf_1", "orf_2", etc; 2. Strands marked as "*" are
#' replaced with "+"; 3. Remove elements with multiple chromosomes or strands (one ORF is on
#' multiple chromosomes or different strands); 4. Remove elements where the ORF length is not
#' divisible by 3; and 5. MOST IMPORTANTLY, if an ORF is on positive strand, sort by
#' coordinates (seqnames, start, end) in ascending order. Otherwise, sort by coordinates
#' (seqnames, end, start) in descending order. The purpose is to achieve the same behavior as
#' cdsBy function in GenomicFeatures package.
#'
#' @param orfGRL A GrangesList object of the ORFs.
#' @param sameChr A logical variable indicating if or not filtering ORFs on multiple
#' chromosomes. (Default: TRUE).
#' @param sameStrand A logical variable indicating if or not filtering ORFs on different
#' strands. (Default: TRUE).
#' @param orfLengthDivisibleBy3 A logical variable indicating if or not filtering ORFs where
#' the ORF length is not divisible by 3. (Default: TRUE).
#'
#' @return A modified GrangesList object of the ORFs with changes stated above.
#'
#' @importFrom BiocGenerics as.data.frame
#' @importFrom dplyr %>% mutate if_else group_by summarize filter arrange desc
#' @importFrom GenomicRanges makeGRangesListFromDataFrame
#' @importFrom GenomeInfoDb seqinfo
.stdGRL = function(orfGRL, sameChr=TRUE, sameStrand=TRUE, orfLengthDivisibleBy3=TRUE) {
  # add names
  if(is.null(names(orfGRL))) {
    names(orfGRL) = sprintf('orf_%d', seq(1, length(orfGRL)))
  }

  # work with a dataframe version of orfGRL as a personal preference
  orfDF = as.data.frame(orfGRL)
  # drop meta columns. The purpose is to minimize the possibility of repeated names with
  # the first 7 columns (group, group_name, seqnames, start, end, width, strand).
  orfDF = orfDF[, 1:7]
  orfDF$seqnames = as.character(orfDF$seqnames)
  orfDF$strand = as.character(orfDF$strand)

  # replace '*' with '+' for strand
  orfDF = orfDF %>% mutate(strand = if_else(strand == '*', '+', strand))

  # check if one ORF has multiple chromosomes or strands, and
  # calculate ORF length
  orfDFSummary = orfDF %>% group_by(group) %>% summarize(
    groupName = unique(group_name),
    chrCount = length(unique(seqnames)),
    strandCount = length(unique(strand)),
    orfLength = sum(width)
  ) %>% as.data.frame()

  # filter
  if(sameChr) {
    lengthBeforeFilter = nrow(orfDFSummary)
    orfDFSummary = orfDFSummary %>% filter(chrCount == 1)
    lengthAfterFilter = nrow(orfDFSummary)
    message(sprintf('%s filtered %d ORFs with multiple chromosomes.', .now(),
      lengthBeforeFilter - lengthAfterFilter))
  }
  if(sameStrand) {
    lengthBeforeFilter = nrow(orfDFSummary)
    orfDFSummary = orfDFSummary %>% filter(strandCount == 1)
    lengthAfterFilter = nrow(orfDFSummary)
    message(sprintf('%s filtered %d ORFs with multiple strands.', .now(),
      lengthBeforeFilter - lengthAfterFilter))
  }
  if(orfLengthDivisibleBy3) {
    lengthBeforeFilter = nrow(orfDFSummary)
    orfDFSummary = orfDFSummary %>% filter(orfLength %% 3 == 0)
    lengthAfterFilter = nrow(orfDFSummary)
    message(sprintf('%s filtered %d ORFs with length not divisible by 3.', .now(),
      lengthBeforeFilter - lengthAfterFilter))
  }
  if(nrow(orfDFSummary) == 0) {
    stop('No ORF passes filtering.')
  }
  message(sprintf('%s %d out of %d ORFs are kept.', .now(), nrow(orfDFSummary),
    length(orfGRL)))

  orfDF = orfDF %>% filter(group %in% unique(orfDFSummary$group))

  # split into positive and negative strand and sort differently
  orfDFPos = orfDF %>% filter(strand == '+')
  orfDFNeg = orfDF %>% filter(strand == '-')

  orfDFPos = orfDFPos %>% arrange(group, seqnames, start, end)
  orfDFNeg = orfDFNeg %>% arrange(group, seqnames, dplyr::desc(end), start)

  orfDF = rbind(orfDFPos, orfDFNeg)

  # rebuild GrangesList
  orfGRLStd = makeGRangesListFromDataFrame(orfDF, split.field='group',
    seqinfo=seqinfo(orfGRL), keep.extra.columns=TRUE)
  names(orfGRLStd) = unique(orfDF$group_name)

  return(orfGRLStd)
}

#' Calculate the percentage of positive elements for each row in a matrix
#'
#' @param counts A matrix of read counts with three rows. Rows represent three frames, and
#' columns represent positions.
#'
#' @return A vector showing the percentage of positive elements for each row.
#'
#'
.getFramePosPct = function(counts) {
  return(rowSums(counts > 0) / ncol(counts))
}

#' Calculate ORFScore sign
#'
#' @param targetFrame A numeric variable of the target frame. This frame is expected to have
#' larger counts than the other two.
#' @param frame1 A numeric variable of frame 1 counts.
#' @param frame2 A numeric variable of frame 2 counts.
#' @param frame3 A numeric variable of frame 3 counts.
#'
#' @return A numeric variable (-1 or +1), depending on:
#' 1. targetFrame == 1: if frame1 > frame2 and frame3, then +1, else -1;
#  2. targetFrame == 2: if frame2 > frame1 and frame3, then +1, else -1;
#  and 3. targetFrame == 3: if frame3 > frame1 and frame2, then +1, else -1.
#'
.getORFScoreSign = function(targetFrame, frame1, frame2, frame3) {
  sign = 0
  if(targetFrame == 1) {
    sign = ifelse(frame1 > max(frame2, frame3), 1, -1)
  } else if(targetFrame == 2) {
    sign = ifelse(frame2 > max(frame1, frame3), 1, -1)
  } else if(targetFrame == 3) {
    sign = ifelse(frame3 > max(frame1, frame2), 1, -1)
  }

  return(sign)
}

#' Calculate raw ORFScore
#'
#' @param frameCounts A numeric vector of length 3 showing the read counts for the three
#' frames.
#' @param p A numeric vector of length 3 indicating the null ditribution for chi-squared test.
#'
#' @return A numeric variable of the test statistic. If the counts for all three frames are
#' zero, return NA.
#'
#' @importFrom stats chisq.test
#'
.getRawORFScore = function(frameCounts, p) {
  rawORFScore = NA
  if(max(frameCounts) > 0) {
    rawORFScore = suppressWarnings(chisq.test(frameCounts, p=p)$statistic)
  }
  names(rawORFScore) = NULL

  return(rawORFScore)
}

#' ORFScore calculation
#'
#' @description ORFScore is firstly defined in Bazzini et al., 2014 (PMID: 24705786), and is
#' used to dscover novel open reading frames (ORF) or rank ORFs showing active translation.
#' Basically, given an ORF, read counts for the three frames are calculated. Then a Chi-squared
#' test statistic is computed by comparing the read counts with an equal null distribution p =
#' c(1/3, 1/3, 1/3). The log2(1 + test statistic) is called ORFScore. In addition, the sign of
#' the ORFScore is positive if the target frame (by default is frame 1) counts are larger than
#' the counts of the other two frames, and negative otherwise.
#'
#' @param bam A \code{GRanges} or \code{GAlignments} object of reads. Note that the reads
#' should be already size selected and shifted. Check function \code{shiftReads} on how to
#' shift reads. Also, for each read, only the 5'-most position is used. (Required).
#' @param orfGRL A \code{GRangesList} object of ORFs. We recommend assigning a unique name to
#' each ORF using \code{names(orfGRL)}. In addition, the following modifications are also
#' applied: 1. If the names of orfGRL are NULL, rename each element as "orf_1", "orf_2", etc;
#' 2. Strands marked as "*" are replaced with "+"; 3. Remove elements with multiple chromosomes
#' or strands (one ORF is on multiple chromosomes or different strands); 4. Remove elements
#' where the ORF length is not divisible by 3; and 5. MOST IMPORTANTLY, if an ORF is on
#' positive strand, sort by coordinates (seqnames, start, end) in ascending order. Otherwise,
#' sort by coordinates (seqnames, end, start) in descending order. The purpose is to achieve
#' the same behavior as \code{cdsBy} function in \code{GenomicFeatures} package. (Required).
#' @param frameOrder A numeric vector of length 3 showing the frames for each position in each
#' ORF. By default, the first position in each ORF is frame 1, the second position is frame 2,
#' and the third position is frame 3. Repeat this pattern afterwards (e.g. 4th position is
#' frame 1, 5th is frame 2, and 6th is frame 3. So on and so forth). (Default: c(1, 2, 3)).
#' @param targetFrame A numeric variable indicating which frame is expected to have higher
#' read counts. By default, frame 1 is expected to have higher read counts than frame 2 and 3.
#' (Default: 1).
#' @param ignoreStrand A logical variable indicating if ignoring that reads and ORFs must
#' be on the same strand. (Default: TRUE).
#' @param probNULL A numeric vector of length 3 showing the null distribution of the read
#' counts in the three frames of an ORF. Must be non-negative and sum up to 1. By default, an
#' equal null distribution is used in the chi-squared test. (Default: c(1/3, 1/3, 1/3)).
#'
#' @return A \code{data.frame} with 9 columns, specified below: 1. Column 1 is ORF ID
#' (\code{orfId}, either user specified in \code{orfGRL} or internally generated);  2. Columns
#' 2 to 4 are the read counts for the three frames where the order is specified by
#' \code{frameOrder} (e.g. \code{frame1Count}, \code{frame2Count}, and \code{frame3Count});
#' 3. Columns 5 to 7 are the percentages of positive counts for the three frames where the
#' order is specified by \code{frameOrder} (e.g. \code{frame1PosPct}, \code{frame2PosPct},
#' and \code{frame3PosPct}). For example, if an ORF has 30 positions (10 positions for each
#' frame), 8 positions for frame 1 are positive, 1 position for frame 2 is positive, and 0
#' position for frame 3 is positive, then column 5 to 7 are 0.8, 0.1, and 0. The purpose of
#' these three columns is to help filtering ORFs with high ORFScore, but the reads only show
#' up in very few positions in the target frame. An example would be an ORF has 300 positions.
#' Frame 1 has 100 read counts, and frame 2 and 3 has 0 read counts. But all the 100 read
#' counts for frame 1 are located in the same position. In this case, the ORFScore will be
#' large (if frame 1 is the target frame), but \code{frame1PosPct} is small (only 0.01). This
#' ORF might be more likely to be a false positive; and 4. Columns 8 and 9 are raw ORFScores
#' (\code{rawORFScore}, test statistics with signs) and final ORFScores (\code{ORFScore},
#' log2(1 + rawORFScore) with signs). If the read counts for all three frames are zero, the
#' raw and final ORFScore is set to \code{NA}. The dataframe is sorted by ORFScore in
#' descending order.
#'
#' @importFrom methods is as
#' @importFrom IRanges resize IntegerList
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicFeatures coverageByTranscript
#' @importFrom dplyr %>% rowwise mutate arrange desc
#'
#' @export
#'
calcORFScore = function(bam, orfGRL, frameOrder=c(1, 2, 3), targetFrame=1, ignoreStrand=TRUE,
  probNULL=c(1/3, 1/3, 1/3)) {
  # sanity check
  if(!is(bam, 'GRanges') & !is(bam, 'GAlignments')) {
    stop('bam must be a GRanges or GAlignments object.')
  }
  if(!is(orfGRL, 'GRangesList')) {
    stop('orfGRL must be a GRangesList object.')
  }
  if(!(all(is.numeric(frameOrder)) & length(frameOrder) == 3 &
    all(sort(frameOrder) == c(1, 2, 3)))) {
    stop(paste('frameOrder must be a numeric vector of length 3 containing the combination',
    'of c(1, 2, 3).'))
  }
  if(!(all(is.numeric(targetFrame)) & length(targetFrame) == 1 &
    all(targetFrame %in% c(1, 2, 3)))) {
    stop('targetFrame must be a numeric variable selecting from c(1, 2, 3).')
  }
  if(!is.logical(ignoreStrand)) {
    stop('ignoreStrand must be a logical variable.')
  }
  if(!(all(is.numeric(probNULL)) & all(probNULL >= 0) & sum(probNULL) == 1)) {
    stop(paste('probNULL must be a numeric vector of length 3, and all elements are',
      'non-negative and sum up to 1.'))
  }

  # convert bam to granges and keep only the 5'-most position
  bamGR = as(bam, 'GRanges')
  bamGR = resize(bamGR, width=1)

  # standadize orf grangeslist
  orfGRLStd = .stdGRL(orfGRL)

  # counts for each position for each orf
  counts = coverageByTranscript(bamGR, orfGRLStd, ignore.strand=ignoreStrand)
  counts = IntegerList(counts)

  # read counts and positive percentage for each frame
  # each orf length in orfGRLStd is divisible by 3
  frameCounts = lapply(lapply(counts, matrix, nrow=3), rowSums)
  frameCounts = as.data.frame(do.call(rbind, frameCounts))
  framePosPct = lapply(lapply(counts, matrix, nrow=3), .getFramePosPct)
  framePosPct = as.data.frame(do.call(rbind, framePosPct))

  # construct orfScore dataframe
  orfScoreDF = cbind(frameCounts, framePosPct)
  colnames(orfScoreDF) = c(sprintf('frame%dCount', frameOrder),
    sprintf('frame%dPosPct', frameOrder))
  orfScoreDF$orfId = rownames(orfScoreDF)
  rownames(orfScoreDF) = NULL
  orfScoreDF = orfScoreDF[, c('orfId', sprintf('frame%dCount', frameOrder),
    sprintf('frame%dPosPct', frameOrder))]

  # ORFScore
  orfScoreDF = orfScoreDF %>% rowwise() %>% mutate(
    orfScoreSign = .getORFScoreSign(targetFrame, frame1Count, frame2Count, frame3Count),
    rawORFScore = orfScoreSign * .getRawORFScore(c(frame1Count, frame2Count, frame3Count),
      p=probNULL),
    orfScore = orfScoreSign * log2(1 + abs(rawORFScore))
  )
  orfScoreDF = orfScoreDF[, c(1, 2:4, 5:7, 9, 10)] # remove ORFScore sign column
  orfScoreDF = orfScoreDF %>% arrange(dplyr::desc(orfScore))
  orfScoreDF = as.data.frame(orfScoreDF)

  return(orfScoreDF)
}
