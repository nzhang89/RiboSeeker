# read distribution

#' Read length distribution
#'
#' @description Read length distribution.
#'
#' @param bam A \code{GAlignments} object of aligned reads. (Required).
#'
#' @return A \code{data.frame} with 3 columns (\code{readLen}, \code{count}, and \code{pct})
#' of read length, number of reads with this length, and percentage of reads with this length.
#'
#' @export
#'
#' @importFrom methods is
#' @importFrom dplyr %>% group_by summarize n mutate
#'
readLengthDist = function(bam) {
  # sanity check
  if(!is(bam, 'GAlignments')) {
    stop('bam must be a GAlignments object.')
  }

  # read length distribution
  dist = data.frame(
    readLen = qwidth(bam)
  )
  dist = dist %>% group_by(readLen) %>% summarize(count = n()) %>%
    mutate(pct = 100 * count / sum(count)) %>% as.data.frame()

  return(dist)
}

#' Genomic feature distribution of aligned reads
#'
#' @description Genomic feature distribution of aligned reads. Reads are assigned to a feature
#' based on the following priority: CDS, UTR, Intron, and Intergenic.
#'
#' @param bam A \code{GAlignments} or \code{GRanges} object of aligned reads. (Required).
#' @param txdb A \code{TxDb} object of genome annotation. See \code{GenomicFeatures} package
#' for more details. (Required).
#' @param category A vector of characters. By default, distribution of reads on CDS, UTR,
#' intron, and intergenic regions is calculated. Must be selected from
#' \code{c('CDS', 'UTR', 'Intron', 'Intergenic')}. (Default: c('CDS', 'UTR', 'Intron',
#' 'Intergenic')).
#' @param fiveEndOnly A logical variable indicating if only keeping the 5'-ends of reads.
#' This option can be used to decrease the ambiguity for reads spanning multiple features.
#' (Default: TRUE).
#' @param ignoreStrand A logical variable indicating if ignoring that reads and features must
#' be on the same strand (Default: TRUE).
#'
#' @return A \code{data.frame} with 3 columns (\code{feature}, \code{count}, and \code{pct})
#' of genomic feature, number of reads falling into this feature, and percentage of reads
#' falling into this feature.
#'
#' @export
#'
#' @importFrom methods is as
#' @importFrom IRanges resize reduce countOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicFeatures genes exons cds fiveUTRsByTranscript threeUTRsByTranscript
#' intronsByTranscript
#' @importFrom dplyr %>% group_by summarize n mutate
#'
readGenomeDist = function(bam, txdb, category=c('CDS', 'UTR', 'Intron', 'Intergenic'),
  fiveEndOnly=TRUE, ignoreStrand=TRUE) {
  # sanity check
  if(!is(bam, 'GAlignments') & !is(bam, 'GRanges')) {
    stop('bam must be a GAlignments or GRanges object.')
  }
  if(!is(txdb, 'TxDb')) {
    stop('txdb must be a TxDb object.')
  }
  if(!all(category %in% c('CDS', 'UTR', 'Intron', 'Intergenic'))) {
    stop('category must be selected from c(\'CDS\', \'UTR\', \'Intron\', \'Intergenic\')')
  }
  if(!is.logical(fiveEndOnly)) {
    stop('fiveEndOnly must be a logical variable.')
  }
  if(!is.logical(ignoreStrand)) {
    stop('ignoreStrand must be a logical variable.')
  }

  # keep only 5'-end to minimize a read overlapping multiple features
  bamGR = as(bam, 'GRanges')
  if(fiveEndOnly) {
    bamGR = resize(as(bam, 'GRanges'), width=1)
  }

  # extract each feature
  geneGR = suppressMessages(reduce(genes(txdb))) # gene
  cdsGR = suppressMessages(reduce(cds(txdb))) # CDS
  fiveUTRGR = suppressMessages(reduce(unlist(fiveUTRsByTranscript(txdb)))) # 5'-UTR
  threeUTRGR = suppressMessages(reduce(unlist(threeUTRsByTranscript(txdb)))) # 3'-UTR
  utrGR = reduce(c(fiveUTRGR, threeUTRGR)) # UTR
  intronGR = suppressMessages(reduce(unlist(intronsByTranscript(txdb)))) # introns

  # internal category counter
  # UA (unassigned) means read has not been assigned to a feature
  categoryIndex = rep('UA', length(bamGR))

  # overlap with features
  # reads are first matched to CDS, then
  # unassigned reads are matched to UTR, then
  # unassigned reads are matched to Intron, then
  # unassigned reads are matched to Intergenic (because reads do not match genes).
  categoryIndex[categoryIndex == 'UA'][countOverlaps(bamGR[categoryIndex == 'UA'], cdsGR,
    ignore.strand=ignoreStrand) > 0] = 'CDS'
  categoryIndex[categoryIndex == 'UA'][countOverlaps(bamGR[categoryIndex == 'UA'], utrGR,
    ignore.strand=ignoreStrand) > 0] = 'UTR'
  categoryIndex[categoryIndex == 'UA'][countOverlaps(bamGR[categoryIndex == 'UA'], intronGR,
    ignore.strand=ignoreStrand) > 0] = 'Intron'
  categoryIndex[categoryIndex == 'UA'][countOverlaps(bamGR[categoryIndex == 'UA'], geneGR,
    ignore.strand=ignoreStrand) > 0] = 'GeneOther'
  categoryIndex[categoryIndex == 'UA'] = 'Intergenic'

  # to summary dataframe
  features = categoryIndex[categoryIndex %in% category]
  features = as.data.frame(table(features))
  colnames(features) = c('feature', 'count')
  features = features %>% mutate(pct = 100 * count / sum(count))

  return(features)
}
