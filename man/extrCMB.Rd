% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Statistics.R
\name{extrCMB}
\alias{extrCMB}
\title{Extreme values}
\usage{
extrCMB(cmbdf, win, n, intensities = "I")
}
\arguments{
\item{cmbdf}{A \code{\link{CMBDataFrame}}.}

\item{win}{A \code{\link{CMBWindow}}}

\item{n}{An integer value.}

\item{intensities}{A \code{\link{CMBDataFrame}} column with measured values.}
}
\value{
A \code{\link{CMBDataFrame}} with \code{n} largest extreme values
}
\description{
This function returns \code{n} largest extreme values for the specified
\code{\link{CMBDataFrame}} column  \code{intensities} and
\code{\link{CMBWindow}} region.
}
\examples{
## Download the map first
# downloadCMBMap(foreground = "smica", nside = 1024)
# df <- CMBDataFrame("CMB_map_smica1024.fits")
# cmbdf <- sampleCMB(df, sample.size = 1000)
#
# win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(0,0,pi/2))
# extrCMB(cmbdf, win1,5)
#
## Ploting the window and 5 top extreme values
# plot(win1)
# plot(extrCMB(cmbdf, win1,5), col ="blue", size = 4,add = TRUE)

}
