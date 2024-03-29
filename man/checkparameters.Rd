% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{checkparameters}
\alias{checkparameters}
\title{check S and E steps based on isotopologue distribution.}
\usage{
checkparameters(resp, M, n, SE, imported, D2)
}
\arguments{
\item{resp}{isotopologue distribution}

\item{M}{total number of carbons of the fatty acid.}

\item{n}{maximum number of elongation steps.}

\item{SE}{list with S and E parameters. NA indicates they must be estimated 
while 0 indicates it does not occur.}

\item{imported}{logical. TRUE if S16 is predefined as 0 (n3 or n6 series).}

\item{D2}{numeric between 0 and 1. Only if D2 >= 0.4 parameters are checked 
based on distribution.}
}
\value{
SE list corrected.
}
\description{
check S and E steps based on isotopologue distribution.
}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
\keyword{internal}
