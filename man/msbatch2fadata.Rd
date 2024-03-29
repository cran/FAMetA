% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{msbatch2fadata}
\alias{msbatch2fadata}
\title{Extract FA data from an annotated msbatch.}
\usage{
msbatch2fadata(msbatch, faid)
}
\arguments{
\item{msbatch}{annotated msbatch.}

\item{faid}{data frame with two columns (ID and Compound) specifying FA ids
and FA names. FA names must be unique and omega series must be indicated
(i.e. FA(20:4)n3, FA(24:1)n9, FA(16:0)). Unknown FA series can be named as nx,
ny, nz to differentiate between isomers.}
}
\value{
fadata.
}
\description{
Extract FA data from an annotated msbatch.
}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
\keyword{internal}
