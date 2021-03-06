% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metagene.R
\name{.calcRelFreq}
\alias{.calcRelFreq}
\title{Convert to row-wise relative frequency}
\usage{
.calcRelFreq(metagene, mode, naVal = 0)
}
\arguments{
\item{metagene}{A matrix of read counts or a list of two matrices representing read counts
in CDS start and end regions. Each row is a transcript, each column is a position, and
value represents read counts in a transcript at a position.}

\item{mode}{A numeric variable (1 or 2) indicating if metagene is a single matrix (mode
is 1) or a list of two matrices (mode is 2).}

\item{naVal}{A numeric variable to replace na with. (Default: 0).}
}
\value{
A matrix or a list of two matrices with values normalized as the row-wise relative
frequency. When mode is 2, the row-wise total counts are calculated from both matrices.
}
\description{
Convert to row-wise relative frequency
}
