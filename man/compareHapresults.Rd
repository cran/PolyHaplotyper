% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyHaplotyper.R
\name{compareHapresults}
\alias{compareHapresults}
\title{compare two haplotyping results}
\usage{
compareHapresults(haploblock, hapresultsA, hapresultsB)
}
\arguments{
\item{haploblock}{a list of character vectors. The names are the names of the
haploblocks, the character vectors have the names of the markers in each
haploblock.}

\item{hapresultsA}{and}

\item{hapresultsB}{two list as returned by inferHaplotypes,
with one item (itself a list) per haploblock with at least a matrix hapdos
and a character vector markers. All
haploblocks in param haploblock must occur in hapresultsA and in hapresultsB.
The individual names (colnames of the hapdos items for each haploblock)
must be identical and in the same order in hapresultsA and hapresultsB}
}
\value{
a list with one element per haploblock in param haploblock.
Each element is itself a list with elements:
$identical: TRUE or FALSE\cr
$message: a single string, "" if the comparison is possible, else the reason
why not (if $message is not "", $identical is always FALSE). The next
elements are only present if $message is "":\cr
$compindiv: a matrix comparing the two hapdos, with one column per individual
and 5 rows: Both_NA, A_NA, B_NA, Equal, Uneq. The last 2 have NA values if A
and/or B is NA
$haplofreq: a matrix with one row per haplotype occurring in A and/or B,
and columns A and B, with the total frequency of each haplotype in
hapdos A or hapdos B
}
\description{
compare two haplotyping results, e.g. PolyHaplotyper and
SATlotyper
}
