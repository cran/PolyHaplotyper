% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolyHaplotyper.R
\name{padded}
\alias{padded}
\title{pad an integer (prefix with zeroes to a fixed length)}
\usage{
padded(x, maxx=0)
}
\arguments{
\item{x}{vector of non-negative integers}

\item{maxx}{a single integer to whose nchar all x elements will be padded;
if 0 (default) the largest value in x will be used}
}
\value{
a character vector representing the values of x left-padded with 0's
to the length of integer maxx or of max(x)
}
\description{
pad an integer (prefix with zeroes to a fixed length)
}
\examples{
padded(c(21, 1, 121, NA, 0))
padded(c(21, 1, 121, NA, 0), maxx=1000)
}
