% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HPDataFrame.R
\name{nside.HPDataFrame}
\alias{nside.HPDataFrame}
\title{HEALPix Nside parameter from a \code{\link{HPDataFrame}}}
\usage{
\method{nside}{HPDataFrame}(x)
}
\arguments{
\item{x}{A \code{\link{HPDataFrame}}.}
}
\value{
The HEALPix Nside parameter.
}
\description{
This function returns the HEALPix Nside parameter
of a \code{\link{HPDataFrame}}
}
\examples{
df <- HPDataFrame(I = rep(0,12), nside = 1)
nside(df)

}
