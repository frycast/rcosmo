% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMBDataFrameGenerics.R
\name{plot.CMBDataFrame}
\alias{plot.CMBDataFrame}
\title{Plot CMB Data}
\usage{
\method{plot}{CMBDataFrame}(
  x,
  intensities = "I",
  add = FALSE,
  sample.size,
  type = "p",
  size = 1,
  box = FALSE,
  axes = FALSE,
  aspect = FALSE,
  col,
  back.col = "black",
  labels,
  hp.boundaries = 0,
  hpb.col = "gray",
  depth_test = "less",
  lab_depth_test = "always",
  ...
)
}
\arguments{
\item{x}{A \code{\link{CMBDataFrame}}.}

\item{intensities}{The name of a column that specifies CMB intensities.
This is only used if \code{col} is unspecified.}

\item{add}{If TRUE then this plot will be added to any existing plot.
Note that if \code{back.col} (see below) is specified then a new plot
window will be opened and \code{add = TRUE} will have no effect.}

\item{sample.size}{Optionally specifies the size of a simple random
sample to take before plotting. This can make the plot less
computationally intensive.}

\item{type}{A single character indicating the type of item to plot.
Supported types are: 'p' for points, 's' for spheres, 'l' for lines,
'h' for line segments from z = 0, and 'n' for nothing.}

\item{size}{The size of plotted points.}

\item{box}{Whether to draw a box.}

\item{axes}{Whether to draw axes.}

\item{aspect}{Either a logical indicating whether to adjust the
aspect ratio, or a new ratio.}

\item{col}{Specify the colour(s) of the plotted points.}

\item{back.col}{Optionally specifies the background colour of
the plot. This argument is passed to rgl::bg3d.}

\item{labels}{Optionally specify a vector of labels to plot,
such as words or vertex indices. If this is specified then
\code{rgl::text3d} is used instead of \code{rgl::plot3d}. Then
\code{length(labels)} must equal \code{nrow(x)}.}

\item{hp.boundaries}{Integer. If greater than 0 then HEALPix
pixel boundaries at \code{nside = hp.boundaries} will be
added to the plot.}

\item{hpb.col}{Colour for the \code{hp.boundaries}.}

\item{depth_test}{The depth test to be applied to the
plotted points. This controls how resistant the plotted
object is to being obscured. See \code{\link[rgl]{rgl.material}}}

\item{lab_depth_test}{The \code{\link{rgl}} depth test
to be applied to the labels and pixel boundaries if present.
See \code{\link[rgl]{rgl.material}}}

\item{...}{Arguments passed to rgl::plot3d.}
}
\value{
A plot of the CMB data
}
\description{
This function produces a plot from a \code{\link{CMBDataFrame}}.
}
\examples{
## First download the map
# downloadCMBMap(foreground = "smica", nside = 1024)
# sky <- CMBDataFrame("CMB_map_smica1024.fits")
# plot(sky, sample.size = 800000)

}
