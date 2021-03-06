% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyHaplotyper.R
\name{mrkdos2mrkdid}
\alias{mrkdos2mrkdid}
\title{get marker dosage IDs from marker dosages}
\usage{
mrkdos2mrkdid(mrkDosage, indiv=NULL, ploidy, check=TRUE)
}
\arguments{
\item{mrkDosage}{matrix or data.frame. Markers are in rows, individuals in
columns, each cell has a marker dosage. Names of individuals are the column
names, marker names are the row names or (if a data.frame) in a column named
MarkerNames. All marker dosages must be in 0:ploidy or NA. If a data.frame,
additional columns may be present.}

\item{indiv}{NULL (default) or a character vector with names of individuals
to be selected. If NULL, all columns are selected;
if mrkDosage is a data.frame, that is probably not what is intended.}

\item{ploidy}{all marker dosages are checked to be in 0:ploidy or NA}

\item{check}{if TRUE (default) checkmrkDosage is called. If FALSE it is
assumed that mrkDosage is a matrix (not a data.frame) and it is not checked.}
}
\value{
a vector of marker dosage IDs, one for each column of mrkDosage:
each a number in 1:((ploidy+1)^nrow(mrkDosage)), NA for each column in
dosages where any of the dosages are NA
}
\description{
get marker dosage IDs (mrkdid) from marker dosages
}
\details{
with ploidy==1 and (of course) all dosages 0 or 1 this function
returns the haplotype numbers for the haplotype specified by each column
}
\examples{
# dosages of 3 markers in 3 tetraploid individuals:
mrkdosg <-
 matrix(c(1,2,2, 4,0,0, 3,0,2), nrow=3,
        dimnames=list(c("mrkA", "mrkB", "mrkC"), c("indiv1", "indiv2", "indiv3")))
# get the "marker dosage IDs":
dids <- mrkdos2mrkdid(mrkDosage=mrkdosg, ploidy=4)
# convert dids back to marker dosages:
mrkdid2mrkdos(dosageIDs=dids, nmrk=3, ploidy=4, mrknames=c("mrkA", "mrkB", "mrkC"))
}
