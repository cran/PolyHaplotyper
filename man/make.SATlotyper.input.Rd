% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyHaplotyper.R
\name{make.SATlotyper.input}
\alias{make.SATlotyper.input}
\title{convert PolyHaplotyper marker data to SATlotyper format}
\usage{
make.SATlotyper.input(mrkDosage, indiv=NULL, markers,
ploidy, phenotype=0, fname)
}
\arguments{
\item{mrkDosage}{matrix or data.frame of allele dosages; same as input for
inferHaplotypes. Markers are in rows, individuals in columns, each cell has
a marker dosage. All marker dosages must be in 0:ploidy or NA}

\item{indiv}{the names of the individuals to include in the ShesisPlus
input data. Default NULL includes all individuals}

\item{markers}{character vector with the names of the markers in the
haploblock; all must occur in mrkDosage}

\item{ploidy}{single integer: the ploidy level}

\item{phenotype}{vector with the phenotypes of all individuals, in order
of the columns of dosmat; default 0}

\item{fname}{filename of a *.csv output file: this will contain the data
in the SATlotyper format (the saved data.frame is also the return value).
If "" no file is written}
}
\value{
a data.frame in the SATlotyper input format: a header row with
"Genotype" and the marker names, and one row per individual with the
individual name plus for each marker the genotype as a sorted
string of <ploidy> A's and B's, or <ploidy> N's
}
\description{
convert PolyHaplotyper marker data (input) to SATlotyper
format for a single haploblock.
}
\examples{
data(PolyHaplotyper_small)
SATin <- make.SATlotyper.input(mrkDosage=phdos, markers=phblocks[[1]],
                              ploidy=6, fname="")
head(SATin)
}
