% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyHaplotyper.R
\name{totHapcombCount}
\alias{totHapcombCount}
\title{calculate the total nr of possible haplotype combinations}
\usage{
totHapcombCount(ploidy, nmrk)
}
\arguments{
\item{ploidy}{a vector of 1 or more ploidy levels}

\item{nmrk}{a vector of 1 or more numbers of markers per haploblock}
}
\value{
a vector with the number of possible haplotype combinations for
each pair of ploidy and nmrk values
}
\description{
calculate the total nr of possible haplotype combinations
}
\details{
nmrk is used to calculate the total number of possible haplotypes
= 2^nmrk. The shorter vector of ploidy and nmrk is recycled.
}
\examples{
totHapcombCount(ploidy=4, nmrk=c(1:8))
}
