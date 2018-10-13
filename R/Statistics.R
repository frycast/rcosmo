#' Sample covariance function
#'
#' This function provides an empirical covariance function for data
#' in a \code{\link{CMBDataFrame}} or data.frame. It assumes that data are from a stationary spherical
#' random field and the covariance depends only on a geodesic distance between locations.
#' Output is a binned covariance.
#'
#' @param cmbdf is a \code{\link{CMBDataFrame}} or \code{\link{data.frame}}
#' @param num.bins specifies the number of bins
#' @param sample.size optionally specify the size of a simple random
#' sample to take before calculating covariance. This may be useful if
#' the full covariance computation is too slow.
#' @param max.dist an optional number between 0 and pi specifying the
#' maximum geodesic distance to use for calculating covariance. Only
#' used if \code{breaks} are unspecified.
#' @param breaks optionally specify the breaks manually using a
#' vector giving the break points between cells. This vector
#' has length \code{num.bins} since the last break point is taken
#' as \code{max.dist}. If \code{equiareal = TRUE} then
#' these breaks should be \eqn{cos(r_i)} where \eqn{r_i} are radii.
#' If \code{equiareal = FALSE} then these breaks should be \eqn{r_i}.
#' @param equiareal if TRUE then the bins have equal spherical
#' area. If false then the bins have equal annular widths.
#' Default is TRUE.
#' @param calc.max.dist if TRUE then the \code{max.dist} will be
#' calculated from the locations in cmbdf. Otherwise
#' either \code{max.dist} must be specified or max.dist will
#' default to pi.
#'
#' @return
#'
#' An object of the class CMBCovariance that is a modification of \code{\link[geoR]{variog}}
#' from the package \strong{geoR}  with variogram values replaced by covariances.
#'
#' The attribute "breaks" contains the break points used to create bins.
#' The result has \code{num.bins + 1} values since the first value, the sample
#' variance,  is not counted as a bin.
#'
#'
#' \describe{
#' \item{u}{a vector with distances.}
#' \item{v}{a vector with estimated covariance values at distances given in u.}
#' \item{n}{number of pairs in each bin}
#' \item{sd}{standard deviation of the values in each bin}
#' \item{bins.lim}{ limits defining the interval spanned by each bin}
#' \item{ind.bin}{a logical vector indicating whether the number of pairs
#' in each bin is greater or equal to the value in the argument pairs.min}
#' \item{var.mark}{variance of the data}
#' \item{beta.ols}{parameters of the mean part of the model fitted by ordinary
#' least squares}
#' \item{output.type}{echoes the option argument}
#' \item{max.dist}{maximum distance between pairs allowed in the covariance calculations}
#' \item{n.data}{number of data}
#' \item{direction}{direction for which the covariance was computed}
#' \item{call}{the function call}
#' }
#'
#'
#'@references \strong{geoR} package, \code{\link[geoR]{variog}}, \code{\link{variogramCMB}}, \code{\link{corrCMB}}
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # Cov <- covCMB(cmbdf, max.dist = 0.03, num.bins = 10)
#' # Cov
#'
#' @export
covCMB <- function(cmbdf,
                   num.bins = 10,
                   sample.size,
                   max.dist = pi,
                   breaks,
                   equiareal = TRUE,
                   calc.max.dist = FALSE)
{
  if (!is.CMBDataFrame(cmbdf)) {
    stop("cmbdf must be a CMBDataFrame")

  }

  if (!missing(sample.size)) {
    cmbdf <- rcosmo::sampleCMB(cmbdf, sample.size = sample.size)

  }

  if (is.null(coords(cmbdf)) || coords(cmbdf) != "cartesian") {
    coords(cmbdf) <- "cartesian"

  }

  if (!all(c("x", "y", "z", "I") %in% names(cmbdf)))
  {
    stop("cmbdf must have columns named 'x', 'y', 'z', 'I' in that order")
  }

  if (missing(breaks))
  {
    if (calc.max.dist)
    {
      max.dist <- maxDist_internal2(cmbdf)
    }


    if (equiareal)
    {
      breaks <-
        seq(cos(max.dist), 1, length.out = num.bins + 1)[-(num.bins + 1)]
    }
    else
    {
      breaks <- cos(rev(seq(0, max.dist, length.out = num.bins + 1)[-1]))
    }

  }

  covs <- covCMB_internal_var(cmbdf[, c("x", "y", "z", "I")], breaks)

  # Reverse order since cosine is decreasing
  covs[2:nrow(covs), 1] <- rev(covs[2:nrow(covs), 1])
  covs[2:nrow(covs), 2] <- rev(covs[2:nrow(covs), 2])
  covs[2:nrow(covs), 3] <- rev(covs[2:nrow(covs), 3])
  v <- c(0, rev(acos(breaks)))
  # Drop the throw-away bin (distances greater than max.dist)
  covs <- covs[-nrow(covs), ]

  # Find the bin centers (break_{i+1} - break_i)/2 with break_0 = 0.
  centers <-
    c(0, v[1:(length(v) - 1)] + (v[2:length(v)] - v[1:(length(v) - 1)]) / 2)

  pairs.min <- 2
  n <- as.integer(c(covs[, 2][1], covs[, 2][-1] / 2))
  indp <- (n >= pairs.min)
  sd1 <- sqrt(covs[,3])
  result <- list(
    u = centers,
    v = covs[, 1],
    n = n,
    sd1 = sd1,
    bins.lim = v,
    ind.bin = indp
  )

  call.fc <- match.call()
  option <- "bin"
  estimator.type <- "classical"
  trend <- "cte"
  data.var <-  stats::var(cmbdf$I)
  n.data <- length(cmbdf$I)
  ##
  beta.ols <- mean(cmbdf$I)
  umax <- max(centers[centers < max.dist])

  result <- c(
    result,
    list(
      var.mark = data.var,
      beta.ols = beta.ols,
      output.type = option,
      max.dist = max.dist,
      estimator.type = estimator.type,
      n.data = n.data,
      lambda = 1,
      trend = trend,
      pairs.min = pairs.min
    )
  )
  result$nugget.tolerance <- 1e-12
  result$direction <- "omnidirectional"
  result$tolerance <- "none"
  result$uvec <- centers
  result$call <- call.fc
  oldClass(result) <- "CMBCovariance"
  return(result)
}

#' Sample correlation function
#'
#' This function provides an empirical correlation function for data
#' in a \code{\link{CMBDataFrame}} or data.frame. It assumes that data are from a stationary spherical
#' random field and the correlation depends only on a geodesic distance between locations.
#' Output is a binned correlation.
#'
#' @param cmbdf is a \code{\link{CMBDataFrame}} or \code{\link{data.frame}}
#' @param num.bins specifies the number of bins
#' @param sample.size optionally specify the size of a simple random
#' sample to take before calculating correlation. This may be useful if
#' the full correlation computation is too slow.
#' @param max.dist an optional number between 0 and pi specifying the
#' maximum geodesic distance to use for calculating correlation. Only
#' used if \code{breaks} are  unspecified.
#' @param breaks optionally specify the breaks manually using a
#' vector giving the break points between cells. This vector
#' has length \code{num.bins} since the last break point is taken
#' as \code{max.dist}. If \code{equiareal = TRUE} then
#' these breaks should be \eqn{cos(r_i)} where \eqn{r_i} are radii.
#' If \code{equiareal = FALSE} then these breaks should be \eqn{r_i}.
#' @param equiareal if TRUE then the bins have equal spherical
#' area. If false then the bins have equal annular widths.
#' Default is TRUE.
#' @param calc.max.dist if TRUE then the \code{max.dist} will be
#' calculated from the locations in cmbdf. Otherwise
#' either \code{max.dist} must be specified or max.dist will
#' default to pi.
#'
#' @return
#'
#' An object of the class CMBCorrelation that is a modification of \code{\link[geoR]{variog}}
#' from the package \strong{geoR} with variogram values replaced by correlation.
#'
#' The attribute "breaks" contains the break points used to create bins.
#' The result has \code{num.bins + 1} values since the first value at distance 0 is not
#' counted as a bin.
#'
#' \describe{
#' \item{u}{a vector with distances.}
#' \item{v}{a vector with estimated correlation values at distances given in u.}
#' \item{n}{number of pairs in each bin}
#' \item{sd}{standard deviation of the values in each bin}
#' \item{bins.lim}{ limits defining the interval spanned by each bin}
#' \item{ind.bin}{a logical vector indicating whether the number of pairs
#' in each bin is greater or equal to the value in the argument pairs.min}
#' \item{var.mark}{variance of the data}
#' \item{beta.ols}{parameters of the mean part of the model fitted by ordinary
#' least squares}
#' \item{output.type}{echoes the option argument}
#' \item{max.dist}{maximum distance between pairs allowed in the correlation calculations}
#' \item{n.data}{number of data}
#' \item{direction}{direction for which the correlation was computed}
#' \item{call}{the function call}
#' }
#'
#'
#'@references \strong{geoR} package,\code{\link[geoR]{variog}},  \code{\link{variogramCMB}}, \code{\link{covCMB}}
#'
#'@examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # corcmb <- corrCMB(cmbdf, max.dist = 0.03, num.bins = 10, sample.size=1000)
#' # corcmb
#'
#' @export
corrCMB <- function(cmbdf,
                         num.bins = 10,
                         sample.size,
                         max.dist = pi,
                         breaks,
                         equiareal = TRUE,
                         calc.max.dist = FALSE) {

  corrCMB<- covCMB(cmbdf, num.bins,
                 sample.size,
                 max.dist,
                 breaks,
                 equiareal ,
                 calc.max.dist)
  corrCMB$v <- corrCMB$v/corrCMB$v[1]
  oldClass(corrCMB) <- "CMBCorrelation"
  return(corrCMB)
}

#' Sample variogram
#'
#' This function provides an empirical variogram for data in a
#' \code{\link{CMBDataFrame}} or data.frame. It assumes that data are from a
#' stationary spherical random field and the covariance depends only on a
#' geodesic distance between locations. Output is a binned variogram.
#'
#'
#' @param cmbdf is a \code{\link{CMBDataFrame}} or \code{\link{data.frame}}
#' @param num.bins specifies the number of bins
#' @param sample.size optionally specify the size of a simple random
#' sample to take before calculating variogram. This may be useful if
#' the full covariance computation is too slow.
#' @param max.dist an optional number between 0 and pi specifying the
#' maximum geodesic distance to use for calculating covariance. Only
#' used if \code{breaks} are  unspecified.
#' @param breaks optionally specify the breaks manually using a
#' vector giving the break points between cells. This vector
#' has length \code{num.bins} since the last break point is taken
#' as \code{max.dist}. If \code{equiareal = TRUE} then
#' these breaks should be \eqn{cos(r_i)} where \eqn{r_i} are radii.
#' If \code{equiareal = FALSE} then these breaks should be \eqn{r_i}.
#' @param equiareal if TRUE then the bins have equal spherical
#' area. If false then the bins have equal annular widths.
#' Default is TRUE.
#' @param calc.max.dist if TRUE then the \code{max.dist} will be
#' calculated from the locations in cmbdf. Otherwise
#' either \code{max.dist} must be specified or max.dist will
#' default to pi.
#'
#' @return
#' An object of class \code{\link[geoR]{variog}} specified in the package \strong{geoR}.
#'
#' The attribute "breaks" contains the break points used to create bins.
#' The result has \code{num.bins + 1} values since the first value at distance 0 is not
#' counted as a bin.
#'
#'\describe{
#' \item{u}{a vector with distances.}
#' \item{v}{a vector with estimated variogram values at distances given in u.}
#' \item{n}{number of pairs in each bin}
#' \item{sd}{standard deviation of the values in each bin}
#' \item{bins.lim}{ limits defining the interval spanned by each bin}
#' \item{ind.bin}{a logical vector indicating whether the number of pairs
#' in each bin is greater or equal to the value in the argument pairs.min}
#' \item{var.mark}{variance of the data}
#' \item{beta.ols}{parameters of the mean part of the model fitted by ordinary
#' least squares}
#' \item{output.type}{echoes the option argument}
#' \item{max.dist}{maximum distance between pairs allowed in the variogram calculations}
#' \item{n.data}{number of data}
#' \item{direction}{direction for which the variogram was computed}
#' \item{call}{the function call}
#' }
#'
#'
#'@references \strong{geoR} package, \code{\link[geoR]{variog}},  \code{\link{covCMB}}, \code{\link{corrCMB}}
#'
#'@examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # varcmb <- variogramCMB(cmbdf, max.dist = 0.1, num.bins = 30, sample.size=100)
#' # varcmb
#'
#' @export
variogramCMB <- function(cmbdf,
                   num.bins = 10,
                   sample.size,
                   max.dist = pi,
                   breaks,
                   equiareal = TRUE,
                   calc.max.dist = FALSE) {

  varCMB<- covCMB(cmbdf, num.bins ,
                     sample.size,
                     max.dist ,
                     breaks,
                     equiareal ,
                     calc.max.dist )
  varCMB$v <- varCMB$v[1]-varCMB$v
  oldClass(varCMB) <- "variogram"
  return(varCMB)
}

#'Plot sample variogram
#'
#'Plots sample (empirical) variogram. Uses \code{\link[geoR]{plot.variogram}} from
#'\strong{geoR} package.
#'
#'@param x An object of class variogram.
#'@param ... Extra arguments as in \code{\link[geoR]{plot.variogram}} passed to plot.default.
#'
#'
#'@return Produces a plot with the sample variogram.
#'
#'@references \strong{geoR} package, \code{\link{variogramCMB}}, \code{\link[geoR]{variog}},
#'\code{\link[geoR]{plot.variogram}}
#'
#'@examples
#'
#' ## Download the map first and call library(geoR)
#' # library(geoR)
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # varcmb <- variogramCMB(cmbdf, max.dist = 0.1, num.bins = 30, sample.size=1000)
#' # plot(varcmb)
#'
#'@name plot.variogram
#'
NULL

#'Plot sample CMBCovariance
#'
#'Plots sample (empirical) covariance function. Uses \code{\link[geoR]{plot.variogram}} from
#'\strong{geoR} package.
#'
#'@param x  An object of class CMBCovariance.
#'@param ...  Extra arguments as in \code{\link[geoR]{plot.variogram}} passed to plot.default.
#'
#'@return Produces a plot with the sample covariance function.
#'
#'@references \strong{geoR} package, \code{\link{covCMB}}, \code{\link[geoR]{variog}},
#'\code{\link[geoR]{plot.variogram}}
#'
#'@examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # Cov <- covCMB(cmbdf, max.dist = 0.03, num.bins = 10)
#' # plot(Cov)
#'
#' @export
plot.CMBCovariance <-  function (x, ...) {
    x0 <- x
    attributes(x0)$class <- "variogram"
    graphics::plot(x0, ylab = "sample covariance", ...)
}

#'Plot sample CMBCorrelation
#'
#'Plots sample (empirical) correlation function. Uses \code{\link[geoR]{plot.variogram}} from
#'\strong{geoR} package.
#'
#'@param x  An object of class CMBCorrelation.
#'@param ...  Extra arguments as in \code{\link{plot.variogram}} passed to plot.default.
#'
#'@return Produces a plot with the sample correlation function.
#'
#'@references \strong{geoR} package, \code{\link{corrCMB}}, \code{\link[geoR]{variog}},
#'\code{\link{plot.variogram}}
#'
#'@examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # corcmb <- corrCMB(cmbdf, max.dist = 0.03, num.bins = 10, sample.size=1000)
#' # plot(corcmb)
#'
#' @export
plot.CMBCorrelation <-  function (x, ...) {
    x0 <- x
    attributes(x0)$class <- "variogram"
    graphics::plot(x0, ylab= "sample correlation", ...)
}


#' Covariance estimate via power spectra
#'
#'This function provides a covariance estimate
#'using the values of the estimated
#'power spectra.
#'
#'
#'@param PowerSpectra a data frame which first
#'column lists values of multipole
#'moments and the second column gives
#'the corresponding values of CMB power
#'spectra.
#'@param Ns a number of points in which
#'the covariance estimate is computed on
#'the interval [-1,1]
#'
#'
#'@return
#' A data frame which first column is 1-d grid starting at -1+1/Ns and
#' finishing at 1 with the step 2/Ns. The second column is the values of
#' estimated covariances on this grid.
#'
#'@references Formula (2.1) in Baran A., Terdik G. Power spectrum estimation
#' of spherical random fields based on covariances. Annales Mathematicae et
#' Informaticae 44 (2015) pp. 15–22.
#'
#' Power Spectra data are from Planck Legacy Archive
#' \url{http://pla.esac.esa.int/pla/#cosmology}
#'
#'
#'@examples
#'
#' ## Download the power spectrum first
#' # N <- 20000
#' # COM_PowerSpectra <- downloadCMBPS(link=1)
#' #
#' # Cov_est <- covPwSp(COM_PowerSpectra[,1:2], N)
#' # plot(Cov_est, type="l")
#'
#' ## Plot the covariance estimate as a function of angular distances
#' # plot(acos(Cov_est[,1]), Cov_est[,2], type ="l",
#' #      xlab ="angular distance", ylab ="Estimated Covariance")
#'
#'@export
covPwSp <- function(PowerSpectra, Ns)
{
  if (requireNamespace("gsl", quietly = TRUE)) {

    Nl <- length(PowerSpectra[, 1])
    Pls <- gsl::legendre_Pl_array(lmax = PowerSpectra[Nl, 1],
                                  x = seq(-1 + 1 / Ns, 1, 2 / Ns))

    Cov_est <- Cov_func(Pls, PowerSpectra[, 2], PowerSpectra[, 1])
    df <- data.frame(t = seq(-1 + 1 / Ns, 1, 2 / Ns), Estimated_Cov = Cov_est)
    return(df)
  } else {

    stop("Package \"gsl\" needed for this function. Please install it.")
  }
}

# Helper function for covPwSp
Cov_func <- function(mat, Dfl , l)  {

    apply(mat[l + 1, ], function(col) {
      sum(Dfl * (2 * l + 1) / (4 * pi * l * (l + 1)) * col)
    }, MARGIN = 2)
}



#'Plot angular scatterplots and means
#'
#'For specified measurements from \code{\link{CMBDataFrame}} this function
#'produces scatterplots and binned means versus theta and phi angles.
#'
#'@param cmbdf A  full \code{\link{CMBDataFrame}} or a windowed
#'\code{\link{CMBDataFrame}}
#'
#'@param intensities  A CMBDataFrame column with measured values
#'
#'@return
#' 2x2 plot. The first row shows scatterplots. The second row gives
#' barplots of the corresponding means computed over bins. The first column
#' corresponds to the values of theta  and the second one is for psi.
#'
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # df.sample <- CMBDataFrame(df, sample.size = 80000)
#' # win <- CMBWindow(theta = c(pi/4,pi/2,pi/2,pi/4), phi = c(0,0,pi/2,pi/2))
#' # cmbdf.win <- window(df.sample, new.window = win)
#' #
#' # intensities <- "I"
#' # plotAngDis(cmbdf.win, intensities)
#'
#'@export
plotAngDis <- function(cmbdf, intensities = "I")
{
  coords(cmbdf) <- "spherical"

  thetabreaks <-
    cut(
      cmbdf$theta,
      breaks = graphics::hist(cmbdf$theta, plot = FALSE)$breaks,
      right = FALSE
    )
  theta.mean <-
    t(tapply(as.data.frame(cmbdf)[, intensities], thetabreaks, mean,
             na.rm = TRUE))

  phibreaks <-
    cut(cmbdf$phi,
        breaks = graphics::hist(cmbdf$phi, plot = FALSE)$breaks,
        right = FALSE)
  phi.mean <-
    t(tapply(as.data.frame(cmbdf)[, intensities], phibreaks, mean,
             na.rm = TRUE))

  graphics::par(mfrow = c(2, 2), mar = c(1, 4, 1, 1) + 0.5)
  graphics::plot(as.data.frame(cmbdf[, c("theta", intensities)]), type = "p", col = "red")
  graphics::plot(
    as.data.frame(cmbdf[, c("phi", intensities)]),
    ylab = "",
    type = "p",
    col = "blue"
  )

  graphics::barplot(
    theta.mean,
    legend = rownames(theta.mean),
    col = "red",
    ylim = 1.1 * range(theta.mean),
    xlab = " ",
    ylab = "Mean Values",
    xaxt = 'n'
  )
  graphics::title(xlab = "theta",
        line = 0,
        cex.lab = 1)
  graphics::barplot(
    phi.mean,
    legend = rownames(phi.mean),
    col = "blue",
    ylim = 1.1 * range(phi.mean),
    xlab = "",
    ylab = "",
    xaxt = 'n'
  )
  graphics::title(xlab = "phi",
        line = 0,
        cex.lab = 1)
}



#' Take a simple random sample from a \code{\link{CMBDataFrame}}
#'
#' This function returns a \code{\link{CMBDataFrame}} with size sample.size,
#' whose rows comprise a simple random sample of the rows
#' from the input CMBDataFrame.
#'
#'@param cmbdf a \code{\link{CMBDataFrame}}.
#'@param sample.size the desired sample size.
#'
#'@return
#' A \code{\link{CMBDataFrame}} with size sample.size,
#' whose rows comprise a simple random sample of the rows
#' from the input CMBDataFrame.
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # plot(sampleCMB(df, sample.size = 800000))
#'
#' df <- CMBDataFrame(nside = 16, I = rnorm(12 * 16 ^ 2), ordering = "nested")
#' df.sample <- sampleCMB(df, sample.size = 100)
#' df.sample
#'
#'@export
sampleCMB <- function(cmbdf, sample.size)
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  srows <- sample(1:nrow(cmbdf), sample.size)
  cmbdf[srows, ]
}




#' First Minkowski functional
#'
#' This function returns an area of the spherical region
#' where measured values
#' are above of the specified threshold level \eqn{alpha}.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param alpha A numeric threshold level.
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'@return
#' The area of the exceedance region
#'
#' @references  Leonenko N., Olenko A. (2014) Sojourn measures of Student
#' and Fisher-Snedecor random fields.  Bernoulli, 20:1454-1483.
#'
#'@examples
#'
#' n <- 64
#' cmbdf <- CMBDataFrame(nside=n, I = rnorm(12*n^2),
#'                       coords = "cartesian",
#'                       ordering = "nested")
#' fmf(cmbdf, 0, 4)
#' fmf(cmbdf, 2, 4)
#'
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' cmbdf.win <- window(cmbdf, new.window = win)
#' fmf(cmbdf.win, 0, 4)
#'
#'@export
fmf <- function(cmbdf, alpha, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  pixelArea(cmbdf)*sum(cmbdf[,intensities, drop = TRUE] > alpha)
}



#' Threshold exceedance probability
#'
#' This function returns an estimated exceedance probability for the specified
#' \code{\link{CMBDataFrame}} column  \code{intensities}, threshold level
#' \eqn{alpha} and \code{\link{CMBWindow}} region.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win A \code{\link{CMBWindow}}
#'@param alpha A numeric threshold level.
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@return
#'
#'Estimated threshold exceedance probability
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 1000)
#'
#' # intensities <- "I"
#' # alpha <- mean(cmbdf[,intensities, drop = TRUE])
#' # alpha
#'
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # exprob(cmbdf, win1, alpha)
#'
#'@export
exprob <- function(cmbdf, win, alpha, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win <- window(cmbdf, new.window = win)
  sum(cmbdf.win[,intensities, drop = TRUE] > alpha)/sum(1-is.na(cmbdf.win[,intensities, drop = TRUE]))
}



#' Quantile-Quantile plots for \code{\link{CMBWindow}}s
#'
#' This function is a modification of standard \link{qqplot} functions to work
#' with \code{\link{CMBWindow}} regions.
#'
#' \code{\link{qqplotWin}} produces a QQ plot of observations in two
#' \code{\link{CMBWindow}}s for the specified \code{\link{CMBDataFrame}} column
#' \code{intensities}. The function automatically adds a diagonal line.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win1 A \code{\link{CMBWindow}}
#'@param win2 A \code{\link{CMBWindow}}
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@return
#'
#' A list with quantile components x and y and a QQ plot with a diagonal line
#'
#'@references \link{qqnormWin}, \link{qqnorm}, \link{qqplot}
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 10000)
#'
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # win2 <- CMBWindow(theta = c(2*pi/3,3*pi/4,3*pi/4, 2*pi/3),
#' #                   phi = c(pi/4,pi/4,pi/3,pi/3))
#'
#' # qqplotWin(cmbdf, win1, win2)
#'
#'@export
qqplotWin <- function(cmbdf, win1, win2, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win1) | !is.CMBWindow(win2) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win1 <- window(cmbdf, new.window = win1)
  cmbdf.win2 <- window(cmbdf, new.window = win2)
  stats::qqplot(cmbdf.win1[,intensities, drop = TRUE], cmbdf.win2[,intensities, drop = TRUE])
  graphics::abline(c(0,1))
  stats::qqplot(cmbdf.win1[,intensities, drop = TRUE], cmbdf.win2[,intensities, drop = TRUE],plot.it = FALSE)
}

#' Normal QQ plot for \code{\link{CMBWindow}}
#'
#' This function is a modification of standard \link{qqnorm} functions to work
#' with \code{\link{CMBWindow}} regions.
#'
#'\code{\link{qqnormWin}} returns a normal QQ plot of for the specified
#'\code{\link{CMBDataFrame}} column \code{intensities} and \code{\link{CMBWindow}}
#'region. The function automatically adds a line of a “theoretical” normal
#'quantile-quantile plot.
#'
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win A \code{\link{CMBWindow}}
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@return
#'
#' A list with quantile components x and y and a normal QQ plot with QQ line
#'
#'@references \link{qqnorm}, \link{qqplot}, \link{qqplotWin}
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 1000)
#'
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # qqnormWin(cmbdf, win1)
#'
#'@export
qqnormWin <- function(cmbdf, win, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win <- window(cmbdf, new.window = win)
  stats::qqnorm(cmbdf.win[,intensities, drop = TRUE])
  stats::qqline(cmbdf.win[,intensities, drop = TRUE])
  stats::qqnorm(cmbdf.win[,intensities, drop = TRUE],plot.it = FALSE)
}


#' CMB Entropy
#'
#'This function returns an estimated entropy for the specified
#'\code{\link{CMBDataFrame}} column  \code{intensities} and \code{\link{CMBWindow}}
#'region. The functions employs the function \link{entropy} and uses histogram
#'counts of \code{intensities} for computations. All arguments of the standard
#'\link{entropy} can be used.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win A \code{\link{CMBWindow}}
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'@param method	 A method to estimate entropy, see \link{entropy}
#'
#'@return
#'
#'Estimated Shannon entropy for observations in \code{\link{CMBWindow}}
#'
#'@references \link{entropy}
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 10000)
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # entropyCMB(cmbdf, win1)
#'
#'@export
entropyCMB <- function(cmbdf, win, intensities = "I", method)
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win <- window(cmbdf, new.window = win)
  y <- graphics::hist(cmbdf.win[,intensities, drop = TRUE], plot = FALSE)$counts

  if (missing(method)) return(entropy::entropy(y))
  return(entropy::entropy(y, method = method))
}

#'Chi-squared statistic for two \code{\link{CMBWindow}}s
#'
#'This function returns the empirical chi-squared statistic for \code{intensities}
#'observations from two \code{\link{CMBWindow}}s of the specified
#'\code{\link{CMBDataFrame}}. The functions employs the function \link{chi2.empirical} and uses histogram
#'counts of \code{intensities} for computations.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win1 A \code{\link{CMBWindow}}
#'@param win2 A \code{\link{CMBWindow}}
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@return
#'
#'Estimated Chi-squared statistic for observations in two
#'\code{\link{CMBWindow}}s.  For small sample sizes and many zero counts
#'Chi-squared statistic is inefficient.
#'
#'@references \link{chi2.empirical}
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 1000)
#' #
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # win2 <- CMBWindow(theta = c(pi/2,pi,pi/2),  phi = c(0,0,pi/2))
#' # plot(win1)
#' # plot(win2,col="green")
#' #
#' # chi2CMB(cmbdf, win1, win2)
#'
#'@export
chi2CMB <- function(cmbdf, win1, win2, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win1) | !is.CMBWindow(win2) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win1 <- window(cmbdf, new.window = win1)
  cmbdf.win2 <- window(cmbdf, new.window = win2)
  y10 <- graphics::hist(cmbdf.win1[,intensities, drop = TRUE], plot = FALSE)$counts
  y20 <- graphics::hist(cmbdf.win2[,intensities, drop = TRUE], plot = FALSE)$counts
  nb <- min(length(y10),length(y20))

  combI <- c(cmbdf.win1[,intensities, drop = TRUE],cmbdf.win2[,intensities, drop = TRUE])
  yb <- seq(min(combI), max(combI), length.out=min(nb))

  y1 <- graphics::hist(cmbdf.win1[,intensities, drop = TRUE], breaks=yb, plot = FALSE)$counts
  y2 <- graphics::hist(cmbdf.win2[,intensities, drop = TRUE], breaks=yb, plot = FALSE)$counts

  entropy::chi2.empirical(y1, y2)
}

#' Sample Renyi function
#'
#' This function computes values of the sample Renyi function. Returns
#' the estimated values of \eqn{T(q)}  for \eqn{q} taking values on a grid.
#' For large data sets could be rather time consuming.
#'
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param q.min Left endpoint of the interval to compute the Renyi function. The default
#'value is 1.01,
#'@param q.max Right endpoint of the interval to compute the Renyi function. The default
#'value is 10
#'@param N Number of points to compute the Renyi function. The default value is 20.
#'@param k.box  A  dyadic decomposition level in computing the Renyi function,
#'see the references in Details. The default value is \eqn{log2(nside(cmbdf)) - 3}
#'@param intensities  A CMBDataFrame column with measured values
#'
#'
#'@return Data frame which first column is the sampling grid
#'\eqn{seq(q.min, q.max, length.out = N)} of \eqn{q} values. Another column
#'consists of values of the sample Renyi function \eqn{T(q)} computed
#'on the grid using the  \eqn{k.box}th level dyadic decomposition
#' of the unit ball.
#'
#'
#'@references
#'(1) Leonenko, N., and Shieh, N. 2013. Rényi function for multifractal
#'random fields. Fractals 21, Article No. 1350009.
#'
#'(2) http://mathworld.wolfram.com/RenyiEntropy.html
#'
#' @examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' #
#' # cmbdf <- CMBDataFrame("CMB_map_smica1024.fits")
#' # win <- CMBWindow(theta = c(pi/4,pi/2,pi/2), phi = c(0,0,pi/2))
#' # cmbdf<- window(cmbdf, new.window = win)
#' # Tq <- fRen(cmbdf)
#' #
#' # plot(Tq[,1], Tq[,2], ylab =expression(D[q]), xlab = "q",
#' # main = "Sample Renyi function", pch = 20, col = "blue")
#'
#'
#' @export
fRen <- function(cmbdf,
                 q.min = 1.01,
                 q.max = 10,
                 N = 20,
                 k.box = log2(nside(cmbdf)) - 3,
                 intensities = "I") {
  if (!is.CMBDataFrame(cmbdf))
  {
    stop("Argument must be a CMBDataFrame")
  }
  ns1 <- nside(cmbdf)
  pixind <- pix.CMBDataFrame(cmbdf)
  nagrpix <- setdiff(1:(12 * ns1 ^ 2), pixind)
  field.comp <- rep(0, 12 * ns1 ^ 2)
    field.in <- cmbdf[, intensities, drop = T]
     minint <- min(field.in)
    field.in <- field.in - minint
  field.final <- replace(field.comp, pixind, field.in)

  res.max <- log2(ns1)
  npix <- 12 * 4 ^ k.box
  delta <- sqrt(4 * pi / npix)
  lev.diff <- 4 ^ (res.max - k.box)

  if (res.max - k.box > 0) {
    nagrpix <- unique(ancestor(nagrpix, res.max - k.box))
  }

  agrpix <- setdiff((1:npix), nagrpix)
  mu <-  vector(mode = "numeric", length = length(agrpix))
  field.total <- 0
  i <- 1
  for (j in agrpix) {
    pixd <- (lev.diff * (j - 1) + 1):(lev.diff * j)
    field.total <- field.total + sum(field.final[pixd])
    mu[i] <- sum(field.final[pixd])
    i <- i + 1
  }
  mu <- mu / field.total
  Q <- seq(q.min, q.max, length.out = N)
  Tq <-  vector(mode = "numeric", length = N)
  ri <- 1
  for (q in Q) {
    Tq[ri] <- 1 / (q - 1) * log2(sum(mu ^ q)) / log2(delta)
    ri <- ri + 1
  }
  Tqf <- data.frame(q = Q, tq = Tq)
  return(Tqf)
}



#' Extreme values
#'
#' This function returns \code{n} largest extreme values for the specified
#' \code{\link{CMBDataFrame}} column  \code{intensities} and
#' \code{\link{CMBWindow}} region.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param win A \code{\link{CMBWindow}}
#'@param n An integer value.
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@return
#'
#'A \code{\link{CMBDataFrame}} with \code{n} largest extreme values
#'
#'@examples
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 1000)
#' #
#' # win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(0,0,pi/2))
#' # extrCMB(cmbdf, win1,5)
#' #
#' ## Ploting the window and 5 top extreme values
#' # plot(win1)
#' # plot(extrCMB(cmbdf, win1,5), col ="blue", size = 4,add = TRUE)
#'
#'@export
extrCMB <- function(cmbdf, win, n, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  if ( !is.CMBWindow(win) )
  {
    stop("Argument must be a CMBWindow")
  }
  cmbdf.win <- window(cmbdf, new.window = win)
  cmbdf.win[order(cmbdf.win[,intensities, drop = TRUE]),][1:n,]
}


#' q-statistic
#'
#'
#'This function returns an estimated q-statistic for the specified  column
#'\code{intensities} in a \code{\link{CMBDataFrame}} and the list of
#'\code{\link{CMBWindow}}s.
#'
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
#'@param listwin A list of \code{\link{CMBWindow}}s
#'@param intensities A \code{\link{CMBDataFrame}} column with measured values.
#'
#'@details The q-statistics is used to measure spatial stratified heterogeneity
#'and takes values in [0, 1]. It gives the percent of the variance of
#'\code{intensities} explained by the stratification. 0 corresponds to no spatial
#'stratified heterogeneity, 1 to perfect spatial stratified heterogeneity.
#'
#'@return
#'
#'Estimated q-statistics for observations in a list of \code{\link{CMBWindow}}s
#'
#'@references
#'Wang, J.F, Zhang, T.L, Fu, B.J. (2016). A measure of spatial
#'stratified heterogeneity. Ecological Indicators. 67: 250–256.
#'
#'@examples
#'
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 1000)
#' # win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' # win2 <- CMBWindow(theta = c(pi/2,pi,pi/2),  phi = c(0,0,pi/2))
#' #
#' # lw <- list(win1, win2)
#' # qstat(cmbdf, lw)
#'
#'@export
qstat <- function(cmbdf, listwin, intensities = "I")
{
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  cmbint <- sapply(listwin,function(v1,v2){
    window(v1, new.window = v2)},v1=cmbdf[,intensities, drop = TRUE])

  ni <- sapply(cmbint, length)
  sigmai <- sapply(cmbint, stats::var)
  sigma <- stats::var(unlist(cmbint, recursive = FALSE))

  1- sum((ni-1)*sigmai)/(sigma*(sum(ni)-1))
}


#' Computes values of covariance functions
#'
#'
#'This function computes the covariances given the separation distance of  locations.
#'Options for different covariance functions on spheres are available. The function uses
#'the function \code{\link[geoR]{cov.spatial}} for covariance models from the package
#'\strong{geoR} and modifies it for additional new models on spheres.
#'
#'
#'@param obj Vector of distances between pairs of spatial locations.
#'@param cov.model A type of the correlation function. Available choices are: "matern",
#'"exponential","spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric". Default is "matern"
#'@param cov.pars A vector with two covariance parameters. The first parameter
#'corresponds to the variance sigma^2. The second parameter corresponds
#'to the range phi of the correlation function.
#'@param kappa A smoothness parameter of the correlation function.
#'
#'@details
#'The function returns the value of the covariance \code{C(h)} at the distance \code{h}.
#'The covariance function has the form
#'
#'\deqn{C(h) = sigma^2 * rho(h/phi).}
#'
#'The parameters of the covariance are positive numbers \code{sigma^2}, \code{phi}
#' and \code{kappa}.
#'
#'Expressions for the correlation functions which are not included in the package
#'\strong{geoR}:
#'
#'\describe{
#' \item{\strong{askey}}{
#' \deqn{rho(h/phi) = (1 - h/phi)^kappa, if h < phi;}
#' \deqn{0, otherwise.}}
#' \item{\strong{c2wendland}}{
#' \deqn{rho(h/phi) =  (1 + kappa * h/phi) * (1 - h/phi)^kappa, if h < phi;}
#' \deqn{0, otherwise.}}
#' \item{\strong{c4wendland}}{
#' \deqn{rho(h/phi) =  (1 + kappa * h/phi + (kappa^2 - 1) * (h/phi)^2 / 3) * (1 - h/phi)^kappa, if h < phi;}
#' \deqn{0, otherwise.}}
#' \item{\strong{sinepower}}{
#' \deqn{rho(h/phi) = 1 - (sin(h/(2 phi))) ^ kappa}}
#'  \item{\strong{multiquadric}}{
#'  \deqn{C(h) =   (1 - phi) ^ (2 * kappa) / (1 + phi^2 - 2 * phi * cos(h))^kappa,
#'  0<phi<1}}
#'  }
#'
#'Additional information can be found in the section Details in
#'\code{\link[geoR]{cov.spatial}}.
#'
#'@return
#'
#'Values of a covariance function for the given distances.
#'
#'@references
#'\strong{geoR} package, \code{\link[geoR]{cov.spatial}}
#'
#'T. Gneiting. Strictly and non-strictly positive definite functions on spheres.
#'Bernoulli 19 (2013), no. 4, 1327-1349.
#'
#'@examples
#'
#'## Compute Askey variogram at x = pi/4
#'
#' 1 - covmodelCMB(pi/4, cov.pars = c(1, pi), kappa = 3, cov.model = "askey" )
#'
#'## Plot of the Askey covariance function
#'
#' v1.f <- function(x, ...) {covmodelCMB(x, ...)}
#' curve( v1.f(x, cov.pars = c(1, 0.03), kappa = 3, cov.model = "askey"),
#' from = 0, to = 0.1, xlab = "distance", ylab = expression(cov(h)), lty = 2,
#' main = "covariance")
#'

#'@export
covmodelCMB  <-  function (obj,
                           cov.model = "matern",
                           cov.pars = stop("no cov.pars argument provided"),
                           kappa = 0.5) {

  fn.env <- sys.frame(sys.nframe())
  .checkCMB.cov.model(
    cov.model = cov.model,
    cov.pars = cov.pars,
    kappa = kappa,
    env = fn.env,
    output = FALSE
  )
  phi <- get("phi", envir = fn.env)
  sigmasq <- get("sigmasq", envir = fn.env)
  covs <- obj
  covs[] <- 0
  for (i in 1:get("ns", envir = fn.env)) {
    if (phi[i] < 1e-16)
      cov.model[i] <- "pure.nugget"
    obj.sc <- obj / phi[i]
    cov.values <- switch(
      cov.model[i],
      askey = ifelse(obj < phi[i], (1 - (obj.sc)) ^ (kappa[i]), 0),
      c2wendland = ifelse(obj < phi[i], (1 + kappa[i] * (obj.sc)) * (1 - (obj.sc)) ^
                            (kappa[i]), 0),
      c4wendland = ifelse(obj < phi[i], (1 + kappa[i] * (obj.sc) + (((kappa[i]) ^
                                                                       2 - 1) * (obj.sc) ^ 2
      ) / 3) * (1 - (obj.sc)) ^ (kappa[i]), 0),
      sinepower = (1 - (sin(obj.sc/ 2)) ^ (kappa[i])),
      multiquadric =  (1 - phi[i]) ^ (2 * kappa[i]) / (1 + (phi[i]) ^ 2 - 2 *
                                                         (phi[i]) * cos(obj)) ^ (kappa[i]),
      pure.nugget = rep(0, length(obj)),
      wave = (1 / obj) * (phi[i] * sin(obj.sc)),
      exponential = exp(-(obj.sc)),
      matern = {
        if (kappa[i] == 0.5)
          exp(-(obj.sc))
        else
          geoR::matern(u = obj,
                 phi = phi[i],
                 kappa = kappa[i])
      },
      gaussian = exp(-((obj.sc) ^ 2)),
      spherical = ifelse(obj < phi[i], (1 - 1.5 * (obj.sc) +
                                          0.5 * (obj.sc) ^ 3), 0),
      circular = {
        obj.sc[obj.sc > 1] <- 1

        ifelse(obj < phi[i], (1 - (2 * ((obj.sc) *
                                          sqrt(1 - ((obj.sc) ^ 2)) +
                                          asin(obj.sc)
        )) / pi), 0)
      },
      cubic = {
        ifelse(obj < phi[i], (1 - (
          7 * (obj.sc ^ 2) -
            8.75 * (obj.sc ^ 3) +
            3.5 * (obj.sc ^ 5) -
            0.75 * (obj.sc ^ 7)
        )), 0)
      },
      power = (obj) ^ phi,
      powered.exponential = exp(-((obj.sc) ^ kappa[i])),
      cauchy = (1 + (obj.sc) ^ 2) ^ (-kappa[i]),
      gneiting = {
        obj.sc <- 0.301187465825 * obj.sc

        t2 <- (1 - obj.sc)

        t2 <- ifelse(t2 > 0, (t2 ^ 8), 0)

        (1 + 8 * obj.sc + 25 * (obj.sc ^ 2) + 32 * (obj.sc ^ 3)) * t2
      },
      gencauchy = (1 + (obj.sc) ^ kappa[2]) ^ (-kappa[1] / kappa[2]),
      gneiting.matern = {
        obj.sc <- 0.301187465825 * obj.sc / kappa[2]

        t2 <- (1 - obj.sc)

        t2 <- ifelse(t2 > 0, (t2 ^ 8), 0)

        cov.values <-
          (1 + 8 * obj.sc + 25 * (obj.sc ^ 2) + 32 * (obj.sc ^ 3)) * t2

        cov.values * geoR::matern(u = obj,
                            phi = phi[i],
                            kappa = kappa[1])

      },
      stop("wrong or no specification of cov.model")
    )
    if (cov.model[i] == "power") {
      A <- (max(cov.values) / sqrt(pi)) * gamma((1 + phi[i]) / 2) *
        gamma(1 - (phi[i] / 2))
      A <- A * max(gamma(phi[i] + (1 + (1:2)) / 2) / (gamma(1 +
                                                              phi[i]) * gamma((1 + (
                                                                1:2
                                                              )) / 2)))
      cov.values <- A - cov.values
      cov.values <- cov.values / max(cov.values)
    }
    cov.values <- ifelse(obj < 1e-16, sigmasq[i], sigmasq[i] *
                           cov.values)
    covs <- covs + cov.values
  }
  if (sum(sigmasq) < 1e-16)
    covs[obj < 1e-16] <- 1
  if (any(!is.finite(covs)))
    warning("Infinity value in covmodelCMB ")
  if (any(is.na(covs)))
    warning("NA value in covmodelCMB ")
  if (any(is.nan(covs)))
    warning("NaN value in covmodelCMB ")
  return(covs)
}

#' Estimates parameters of variograms
#'
#'
#'This function estimates variogram parameters by fitting a parametric model
#'from \code{\link{covmodelCMB}} to a sample variogram. The function modifies
#'\code{\link[geoR]{variofit}} from the package \strong{geoR}
#'for additional covariance models on spheres.
#'
#'@param vario An object of the class \code{variogram} obtained as an output of
#'the function \code{\link{variogramCMB}}.
#'@param ini.cov.pars A vector with initial values for the variogram parameters.
#'The first parameter corresponds to the variance sigma^2. The second parameter
#'corresponds to the range phi of the correlation function.
#'@param cov.model A type of the variogram function. Available choices are: "matern",
#'"exponential","spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric". The default is "matern"
#'@param fix.nugget logical. Indicates whether the nugget variance should be regarded
#'as fixed or be estimated. The default is FALSE.
#'@param nugget A value for the nugget parameter. Regarded as a fixed values if
#'\code{fix.nugget = TRUE} or as a initial value for the minimization algorithm if
#'\code{fix.nugget = FALSE}. The default is zero.
#'@param fix.kappa logical. Indicates whether the parameter kappa should be regarded
#'as fixed or be estimated. The default is TRUE.
#'@param kappa A value for the smoothness parameter. Regarded as a fixed values if
#'\code{fix.kappa = TRUE} or as a initial value for the minimization algorithm if
#'\code{fix.kappa = FALSE}. Required not in all covariance models, see
#'\code{\link{covmodelCMB}}. The default is 0.5.
#'@param simul.number number of simulation. Used if \code{vario} has empirical variograms
#'for more than one data-set (simulations). The default is NULL
#'@param max.dist A maximum distance to fit a variogram model. The default is
#'\code{x$max.dist}.
#'@param weights Weights used in the loss function in the minimization algorithm.
#'@param limits Lower and upper limits for the model parameters used
#'in the numerical minimisation by \code{minimisation.function = "optim"}.
#'@param minimisation.function Minimization function ("optim", "nlm", "nls") to estimate
#'the parameters.
#'@param messages logical. Indicates whether or not status messages are printed on
#'the screen.
#'@param ... other minimisation parameters
#'
#'
#'@details
#'The parameter values of a variogram function from \code{\link{covmodelCMB}} are
#'found by numerical optimization using one of the functions: \code{\link{optim}},
#'\code{\link{nlm}} and \code{\link{nls}}.
#'
#'The function modifies \code{\link[geoR]{variofit}} from the package \strong{geoR}
#'for additional variogram models on spheres.  Available models are: "matern",
#'"exponential", "spherical", "powered.exponential", "cauchy", "gencauchy",
#'"pure.nugget", "askey", "c2wendland", "c4wendland", "sinepower", "multiquadric".
#'
#'Additionally it rescales an empirical variogram to the range \code{[0,1]} before
#'numerical optimisation and then transforms all obtained results to the original
#'scale. If \code{ini.cov.pars} are not provided then the 5x5 grid
#'\code{(seq(0,max(vario$v),l=5), seq(0,vario$max.dist,l=5))}
#'of initial values of sigma^2  and phi is used.
#'
#'
#'@return
#'An object of the class \code{variomodel} and \code{variofit}, see
#'\code{\link[geoR]{variofit}}
#'
#'@references
#'\strong{geoR} package, \code{\link[geoR]{variofit}}, \code{\link{covmodelCMB}}
#'
#'@examples
#' #
#' # df <- CMBDataFrame("../CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 10000)
#' # varcmb <- variogramCMB(cmbdf, max.dist = 0.1, num.bins = 30)
#' # varcmb
#' #
#' # ols <- variofitCMB(varcmb,  fix.nug=FALSE, wei="equal", cov.model= "matern")
#' # plot(varcmb)
#' # lines(ols, lty=2)
#' # str(ols)
#' #
#' # ols <- variofitCMB(vario1, fix.nug = TRUE, kappa = 3, wei = "equal",
#' # cov.model = "askey")
#' # plot(varcmb, main = ols$cov.model)
#' # linesCMB(ols, lty = 2)
#' # str(ols)
#'
#'@export
variofitCMB <- function (vario, ini.cov.pars, cov.model, fix.nugget = FALSE,
                         nugget = 0, fix.kappa = TRUE, kappa = 0.5, simul.number = NULL,
                         max.dist = vario$max.dist, weights, minimisation.function,
                         limits = geoR::pars.limits(), messages, ...) {
  cov.model <- match.arg(cov.model, choices = CMB.cov.models)
  vario1 <- vario
  vario1$v <- vario1$v/max(vario1$v)
  if (missing(ini.cov.pars)){
    ini.cov.pars <- expand.grid(seq(0,max(vario1$v),l=5), seq(0,vario1$max.dist,l=5))
  }
  if (cov.model %in% c("matern",
                       "exponential",
                       "spherical",
                       "powered.exponential",
                       "cauchy",
                       "gencauchy",
                       "pure.nugget")){
    variofitCMB <- geoR::variofit(vario1, ini.cov.pars=ini.cov.pars, cov.model=cov.model, fix.nugget=fix.nugget,
                                  nugget=nugget, fix.kappa=fix.kappa, kappa=kappa, simul.number=simul.number,
                                  max.dist=max.dist, weights=weights, minimisation.function=minimisation.function,
                                  limits = limits, messages=messages, ...)
  }
  if (cov.model %in% c( "askey",
                        "c2wendland",
                        "c4wendland",
                        "sinepower",
                        "multiquadric")){
    variofitCMB <- variofit1(vario1, ini.cov.pars=ini.cov.pars, cov.model=cov.model, fix.nugget=fix.nugget,
                             nugget=nugget, fix.kappa=fix.kappa, kappa=kappa, simul.number=simul.number,
                             max.dist=max.dist, weights=weights, minimisation.function=minimisation.function,
                             limits = limits, messages=messages, ...)
  }
  variofitCMB$cov.pars[1] <- variofitCMB$cov.pars[1]*max(vario$v)
  variofitCMB$nugget <- variofitCMB$nugget*max(vario$v)
  return(variofitCMB)
}

# Helper function for variofitCMB
variofit1 <-   function (vario,
            ini.cov.pars,
            cov.model,
            fix.nugget,
            nugget,
            fix.kappa,
            kappa,
            simul.number,
            max.dist,
            weights,
            minimisation.function,
            limits,
            messages,
            ...) {
    call.fc <- match.call()
    if (missing(messages))
      messages.screen <-
        as.logical(ifelse(is.null(getOption("geoR.messages")),
                          TRUE, getOption("geoR.messages")))
    else
      messages.screen <- messages
    if (length(class(vario)) == 0 ||
        all(class(vario) != "variogram"))
      warning("object vario should preferably be of the geoR's class \"variogram\"")
    if (!missing(ini.cov.pars)) {
      if (any(class(ini.cov.pars) == "eyefit"))
        cov.model <- ini.cov.pars[[1]]$cov.model
      if (any(class(ini.cov.pars) == "variomodel"))
        cov.model <- ini.cov.pars$cov.model
    }
    if (missing(cov.model))
      cov.model <- "matern"
    cov.model <- match.arg(cov.model, choices = CMB.cov.models)
    if (cov.model == "stable")
      cov.model <- "powered.exponential"
    if (cov.model == "powered.exponential")
      if (limits$kappa["upper"] > 2)
        limits$kappa["upper"] <- 2
    if (missing(weights)) {
      if (vario$output.type == "cloud")
        weights <- "equal"
      else
        weights <- "npairs"
    }
    else
      weights <- match.arg(weights, choices = c("npairs",
                                                "equal", "cressie"))
    if (messages.screen) {
      cat(paste("variofit: covariance model used is", cov.model,
                "\n"))
      cat(paste("variofit: weights used:", weights, "\n"))
    }
    if (missing(minimisation.function))
      minimisation.function <- "optim"
    if (any(cov.model == c("linear", "power")) &
        minimisation.function ==
        "nls") {
      cat(
        "warning: minimisation function nls can not be used with given cov.model.\n          changing for \"optim\".\n"
      )
      minimisation.function <- "optim"
    }
    if (minimisation.function == "nls" & weights != "equal") {
      warning(
        "variofit: minimisation function nls can only be used with weights=\"equal\".\n          changing for \"optim\".\n"
      )
      minimisation.function <- "optim"
    }
    if (is.matrix(vario$v) & is.null(simul.number))
      stop(
        "object in vario$v is a matrix. This function works for only 1 empirical variogram at once\n"
      )
    if (!is.null(simul.number))
      vario$v <- vario$v[, simul.number]
    if (mode(max.dist) != "numeric" || length(max.dist) > 1)
      stop("a single numerical value must be provided in the argument max.dist")
    if (max.dist == vario$max.dist)
      XY <- list(u = vario$u, v = vario$v, n = vario$n)
    else
      XY <- list(u = vario$u[vario$u <= max.dist],
                 v = vario$v[vario$u <=
                               max.dist],
                 n = vario$n[vario$u <= max.dist])
    if (cov.model == "pure.nugget") {
      minimisation.function <- "not used"
      message <-
        "correlation function does not require numerical minimisation"
      if (weights == "equal")
        lm.wei <- rep(1, length(XY$u))
      else
        lm.wei <- XY$n
      if (cov.model == "pure.nugget") {
        if (fix.nugget) {
          temp <- stats::lm((XY$v - nugget) ~ 1, weights = lm.wei)
          cov.pars <- c(temp$coef, 0)
        }
        else {
          temp <- stats::lm(XY$v ~ 1, weights = lm.wei)
          nugget <- temp$coef
          cov.pars <- c(0, 0)
        }
      }
      value <- sum((temp$residuals) ^ 2)
    }
    else {
      if (messages.screen)
        cat(paste(
          "variofit: minimisation function used:",
          minimisation.function,
          "\n"
        ))
      umax <- max(vario$u)
      vmax <- max(vario$v)
      if (missing(ini.cov.pars)) {
        ini.cov.pars <- as.matrix(expand.grid(c(vmax / 2, 3 *
                                                  vmax / 4, vmax), seq(0, 0.8 * umax, len = 6)))
        if (!fix.nugget)
          nugget <- unique(c(nugget, vmax / 10, vmax / 4, vmax / 2))
        if (!fix.kappa)
          kappa <- unique(c(kappa, 0.25, 0.5, 1, 1.5, 2))
        if (messages.screen)
          warning("initial values not provided - running the default search")
      }
      else {
        if (any(class(ini.cov.pars) == "eyefit")) {
          init <- nugget <- kappa <- NULL
          for (i in 1:length(ini.cov.pars)) {
            init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
            nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
            if (cov.model == "gneiting.matern")
              kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
            else
              kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
          }
          ini.cov.pars <- init
        }
        if (any(class(ini.cov.pars) == "variomodel")) {
          nugget <- ini.cov.pars$nugget
          kappa <- ini.cov.pars$kappa
          ini.cov.pars <- ini.cov.pars$cov.pars
        }
      }
      if (is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)) {
        ini.cov.pars <- as.matrix(ini.cov.pars)
        if (nrow(ini.cov.pars) == 1)
          ini.cov.pars <- as.vector(ini.cov.pars)
        else {
          if (ncol(ini.cov.pars) != 2)
            stop(
              "\nini.cov.pars must be a matrix or data.frame with 2 components:
              \ninitial values for sigmasq (partial sill) and phi (range parameter)\n"
            )
        }
      }
      else if (length(ini.cov.pars) > 2)
        stop("\nini.cov.pars must provide initial values for sigmasq and phi\n")
      if (is.matrix(ini.cov.pars) | (length(nugget) > 1) |
          (length(kappa) > 1)) {
        if (messages.screen)
          cat("variofit: searching for best initial value ...")
        ini.temp <- matrix(ini.cov.pars, ncol = 2)
        grid.ini <- as.matrix(expand.grid(
          sigmasq = unique(ini.temp[,
                                    1]),
          phi = unique(ini.temp[, 2]),
          tausq = unique(nugget),
          kappa = unique(kappa)
        ))
        v.loss <- function(parms, u, v, n, cov.model, weights) {
          sigmasq <- parms[1]
          phi <- parms[2]
          if (cov.model == "power")
            phi <- 2 * exp(phi) / (1 + exp(phi))
          tausq <- parms[3]
          kappa <- parms[4]
          if (cov.model == "power")
            v.mod <- tausq + covmodelCMB(u,
              cov.pars = c(sigmasq,
                           phi),
              cov.model = "power",
              kappa = kappa
            )
          else
            v.mod <- (sigmasq + tausq) - covmodelCMB(u,
              cov.pars = c(sigmasq, phi),
              cov.model = cov.model,
              kappa = kappa
            )
          if (weights == "equal")
            loss <- sum((v - v.mod) ^ 2)
          if (weights == "npairs")
            loss <- sum(n * (v - v.mod) ^ 2)
          if (weights == "cressie")
            loss <- sum((n / (v.mod ^ 2)) * (v - v.mod) ^ 2)
          return(loss)
        }
        grid.loss <- apply(
          grid.ini,
          1,
          v.loss,
          u = XY$u,
          v = XY$v,
          n = XY$n,
          cov.model = cov.model,
          weights = weights
        )
        ini.temp <- grid.ini[which(grid.loss == min(grid.loss))[1],
                             , drop = FALSE]
        if (is.R())
          rownames(ini.temp) <- "initial.value"
        if (messages.screen) {
          cat(" selected values:\n")
          print(rbind(round(ini.temp, digits = 2), status = ifelse(
            c(FALSE,
              FALSE, fix.nugget, fix.kappa), "fix", "est"
          )))
          cat(paste("loss value:", min(grid.loss), "\n"))
        }
        names(ini.temp) <- NULL
        ini.cov.pars <- ini.temp[1:2]
        nugget <- ini.temp[3]
        kappa <- ini.temp[4]
        grid.ini <- NULL
      }
      if (ini.cov.pars[1] > 2 * vmax)
        warning("unreasonable initial value for sigmasq (too high)")
      if (ini.cov.pars[1] + nugget > 3 * vmax)
        warning("unreasonable initial value for sigmasq + nugget (too high)")
      if (vario$output.type != "cloud") {
        if (ini.cov.pars[1] + nugget < 0.3 * vmax)
          warning("unreasonable initial value for sigmasq + nugget (too low)")
      }
      if (nugget > 2 * vmax)
        warning("unreasonable initial value for nugget (too high)")
      if (ini.cov.pars[2] > 1.5 * umax)
        warning("unreasonable initial value for phi (too high)")
      if (!fix.kappa) {
        if (cov.model == "powered.exponential")
          Tkappa.ini <- log(kappa / (2 - kappa))
        else
          Tkappa.ini <- log(kappa)
      }
      if (minimisation.function == "nls") {
        if (ini.cov.pars[2] == 0)
          ini.cov.pars <- max(XY$u) / 10
        if (kappa == 0)
          kappa <- 0.5
        if (cov.model == "power")
          Tphi.ini <- log(ini.cov.pars[2] / (2 - ini.cov.pars[2]))
        else
          Tphi.ini <- log(ini.cov.pars[2])
        XY$cov.model <- cov.model
        if (fix.nugget) {
          XY$nugget <- as.vector(nugget)
          if (fix.kappa) {
            XY$kappa <- as.vector(kappa)
            res <- stats::nls((v - nugget) ~ matrix((
              1 - covmodelCMB (u,
                cov.pars = c(1, exp(Tphi)),
                cov.model = cov.model,
                kappa = kappa
              )
            ), ncol = 1),
            start = list(Tphi = Tphi.ini),
            data = XY,
            algorithm = "plinear",
            ...
            )
          }
          else {
            if (cov.model == "powered.exponential")
              res <- stats::nls((v - nugget) ~ matrix((
                1 - covmodelCMB(u,
                  cov.pars = c(1, exp(Tphi)),
                  cov.model = cov.model,
                  kappa = (2 * exp(Tkappa) /
                             (1 + exp(Tkappa)))
                )
              ),
              ncol = 1),
              start = list(Tphi = Tphi.ini,
                           Tkappa = Tkappa.ini),
              data = XY,
              algorithm = "plinear",
              ...
              )
            else
              res <- stats::nls((v - nugget) ~ matrix((
                1 -
                  covmodelCMB(u,
                    cov.pars = c(1, exp(Tphi)),
                    cov.model = cov.model,
                    kappa = exp(Tkappa)
                  )
              ),
              ncol = 1),
              start = list(Tphi = Tphi.ini,
                           Tkappa = Tkappa.ini),
              data = XY,
              algorithm = "plinear",
              ...
              )
            kappa <- exp(stats::coef(res)["Tkappa"])
            names(kappa) <- NULL
          }
          cov.pars <- stats::coef(res)[c(".lin", "Tphi")]
          names(cov.pars) <- NULL
        }
        else {
          if (fix.kappa) {
            XY$kappa <- kappa
            res <- stats::nls(
              v ~ cbind(1, (
                1 - covmodelCMB (
                  u,
                  cov.pars = c(1, exp(Tphi)),
                  cov.model = cov.model,
                  kappa = kappa
                )
              )),
              start = list(Tphi = Tphi.ini),
              algorithm = "plinear",
              data = XY,
              ...
            )
          }
          else {
            if (cov.model == "powered.exponential")
              res <- stats::nls(
                v ~ cbind(1, (
                  1 - covmodelCMB (
                    u,
                    cov.pars = c(1, exp(Tphi)),
                    cov.model = cov.model,
                    kappa = (2 * exp(Tkappa) /
                               (1 + exp(Tkappa)))
                  )
                )),
                start = list(Tphi = Tphi.ini, Tkappa = Tkappa.ini),
                algorithm = "plinear",
                data = XY,
                ...
              )
            else
              res <- stats::nls(
                v ~ cbind(1, (
                  1 - covmodelCMB (
                    u,
                    cov.pars = c(1, exp(Tphi)),
                    cov.model = cov.model,
                    kappa = exp(Tkappa)
                  )
                )),
                start = list(Tphi = Tphi.ini,
                             Tkappa = Tkappa.ini),
                algorithm = "plinear",
                data = XY,
                ...
              )
            kappa <- exp(stats::coef(res)["Tkappa"])
            names(kappa) <- NULL
          }
          nugget <- stats::coef(res)[".lin1"]
          names(nugget) <- NULL
          cov.pars <- stats::coef(res)[c(".lin2", "Tphi")]
          names(cov.pars) <- NULL
        }
        if (cov.model == "power")
          cov.pars[2] <-
          2 * exp(cov.pars[2]) / (1 + exp(cov.pars[2]))
        else
          cov.pars[2] <- exp(cov.pars[2])
        if (nugget < 0 | cov.pars[1] < 0) {
          warning(
            "\nvariofit: negative variance parameter found using the default option \"nls\".\n        Try another minimisation function and/or fix some of the parameters.\n"
          )
          temp <- c(
            sigmasq = cov.pars[1],
            phi = cov.pars[2],
            tausq = nugget,
            kappa = kappa
          )
          print(rbind(round(temp, digits = 4), status = ifelse(
            c(FALSE,
              FALSE, fix.nugget, fix.kappa), "fix", "est"
          )))
          return(invisible())
        }
        value <- sum(stats::resid(res) ^ 2)
        message <- "nls does not provides convergence message"
      }
      if (minimisation.function == "nlm" | minimisation.function ==
          "optim") {
        .global.list <- list(
          u = XY$u,
          v = XY$v,
          n = XY$n,
          fix.nugget = fix.nugget,
          nugget = nugget,
          fix.kappa = fix.kappa,
          kappa = kappa,
          cov.model = cov.model,
          m.f = minimisation.function,
          weights = weights
        )
        ini <- ini.cov.pars
        if (cov.model == "power")
          ini[2] <- log(ini[2] / (2 - ini[2]))
        if (cov.model == "linear")
          ini <- ini[1]
        if (fix.nugget == FALSE)
          ini <- c(ini, nugget)
        if (!fix.kappa)
          ini <- c(ini, Tkappa.ini)
        names(ini) <- NULL
        if (minimisation.function == "nlm") {
          result <- stats::nlm(.loss1.vario, ini, g.l = .global.list,
                        ...)
          result$par <- result$estimate
          result$value <- result$minimum
          result$convergence <- result$code
          if (!is.null(get(".temp.theta", pos = 1)))
            result$par <- get(".temp.theta", pos = 1)
        }
        else {
          lower.l <- sapply(limits, function(x)
            x[1])
          upper.l <- sapply(limits, function(x)
            x[2])
          if (fix.kappa == FALSE) {
            if (fix.nugget) {
              lower <- lower.l[c("sigmasq.lower", "phi.lower",
                                 "kappa.lower")]
              upper <- upper.l[c("sigmasq.upper", "phi.upper",
                                 "kappa.upper")]
            }
            else {
              lower <- lower.l[c("sigmasq.lower",
                                 "phi.lower",
                                 "tausq.rel.lower",
                                 "kappa.lower")]
              upper <- upper.l[c("sigmasq.upper",
                                 "phi.upper",
                                 "tausq.rel.upper",
                                 "kappa.upper")]
            }
          }
          else {
            if (cov.model == "power") {
              if (fix.nugget) {
                lower <- lower.l[c("sigmasq.lower", "phi.lower")]
                upper <- upper.l[c("sigmasq.upper", "phi.upper")]
              }
              else {
                lower <- lower.l[c("sigmasq.lower",
                                   "phi.lower",
                                   "tausq.rel.lower")]
                upper <- upper.l[c("sigmasq.upper",
                                   "phi.upper",
                                   "tausq.rel.upper")]
              }
            }
            else {
              lower <- lower.l["phi.lower"]
              upper <- upper.l["phi.upper"]
            }
          }
          result <- stats::optim(
            ini,
            .loss1.vario,
            method = "L-BFGS-B",
            hessian = TRUE,
            lower = lower,
            upper = upper,
            g.l = .global.list,
            ...
          )
        }
        value <- result$value
        message <- paste(minimisation.function,
                         "convergence code:",
                         result$convergence)
        if (cov.model == "linear")
          result$par <- c(result$par[1], 1, result$par[-1])
        cov.pars <- as.vector(result$par[1:2])
        if (cov.model == "power")
          cov.pars[2] <-
          2 * exp(cov.pars[2]) / (1 + exp(cov.pars[2]))
        if (!fix.kappa) {
          if (fix.nugget)
            kappa <- result$par[3]
          else {
            nugget <- result$par[3]
            kappa <- result$par[4]
          }
          if (.global.list$cov.model == "powered.exponential")
            kappa <- 2 * (exp(kappa)) / (1 + exp(kappa))
          else
            kappa <- exp(kappa)
        }
        else if (!fix.nugget)
          nugget <- result$par[3]
      }
    }
    if (cov.model == "matern")
      #     kappa <- min(0.5,kappa)
      if (cov.model == "powered.exponential")
        kappa <- min(1,kappa)
    if (cov.model == "askey")
      kappa <- max(2,kappa)
    if (cov.model == "c2wendland")
      kappa <- max(4,kappa)
    if (cov.model == "c4wendland")
      kappa <- max(6,kappa)
    if (cov.model == "sinepower")
      kappa <- min(2,kappa)
    estimation <- list(
      nugget = nugget,
      cov.pars = cov.pars,
      cov.model = cov.model,
      kappa = kappa,
      value = value,
      trend = vario$trend,
      beta.ols = vario$beta.ols,
      practicalRange = practicalRangeCMB(
        cov.model = cov.model,
        phi = cov.pars[2],
        kappa = kappa
      ),
      max.dist = max.dist,
      minimisation.function = minimisation.function
    )
    estimation$weights <- weights
    if (weights == "equal")
      estimation$method <- "OLS"
    else
      estimation$method <- "WLS"
    estimation$fix.nugget <- fix.nugget
    estimation$fix.kappa <- fix.kappa
    estimation$lambda <- vario$lambda
    estimation$message <- message
    estimation$call <- call.fc
    oldClass(estimation) <- c("variomodel", "variofit")
    return(estimation)
  }

#'Plot theoretical CMBCovariance
#'
#'Plots theoretical covariance functions from the list defined in \code{\link{covmodelCMB}}
#'
#'@param cov.model A type of the correlation function. Available choices are: "matern",
#'"exponential","spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric". The default is "matern"
#'@param sigmasq The variance parameter as documented in \code{\link{covmodelCMB}}.
#'The default is 1.
#'@param phi The range parameter as documented in \code{\link{covmodelCMB}}. The default is
#'\code{pi}.
#'@param kappa A smoothness parameter of the correlation function. The default is 0.5.
#'@param from A lower range of the plotting region. The default is \code{lb =0}
#'@param to An upper range of the plotting region. The default is \code{ub= pi}.
#'@param ...  optional plotting parameters.
#'
#'
#'@return Produces a plot with the theoretical covariance function.
#'
#'@references  \code{\link{covmodelCMB}}
#'
#'@examples
#'
#' plotcovmodelCMB("matern", sigmasq = 5)
#' plotcovmodelCMB("askey", phi = pi/4, to  = pi/2, kappa = 4)
#'
#'
#'@export
plotcovmodelCMB <- function (cov.model = "matern", sigmasq=1,
                             phi = pi,
                             kappa = 0.5,
                             from =0 , to = pi, ...){
  cov.model <- match.arg(cov.model, choices = CMB.cov.models)
  .checkCMB.cov.model(cov.model = cov.model,
                      cov.pars = c(sigmasq, phi),
                      kappa = kappa,
                      output = FALSE)
    x <- NULL; rm(x) # Trick R CMD Check
    graphics::curve(covmodelCMB(x, cov.pars = c(sigmasq, phi),
                                kappa = kappa, cov.model = cov.model),
                  from = from,
                  to = to,
                  xlab = "distance",
                  ylab = expression(C(h)),
                  lty = 2,
                  main = bquote(paste(.(cov.model)~"covariance function"))
  )
}

#'Plot theoretical variogram
#'
#'Plots theoretical variogram functions from the list defined in \code{\link{covmodelCMB}}
#'
#'@param cov.model A type of the variogram function. Available choices are: "matern",
#'"exponential","spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric". The default is "matern"
#'@param sigmasq The variance parameter as documented in \code{\link{covmodelCMB}}.
#'The default is 1.
#'@param phi The range parameter as documented in \code{\link{covmodelCMB}}. The default is
#'\code{pi}.
#'@param kappa A smoothness parameter of the variogram function. The default is 0.5.
#'@param from A lower range of the plotting region. The default is \code{lb =0}
#'@param to An upper range of the plotting region. The default is \code{ub= pi}.
#'@param ...  optional plotting parameters.
#'
#'
#'@return Produces a plot with the theoretical variogram.
#'
#'@references  \code{\link{covmodelCMB}}
#'
#'@examples
#'
#' plotvariogram("matern", sigmasq = 5)
#' plotvariogram("askey", phi = pi/4, to  = pi/2, kappa = 4)
#'
#'
#'@export
plotvariogram <- function (cov.model = "matern", sigmasq=1,
                             phi = pi,
                             kappa = 0.5,
                             from =0 , to = pi, ...){
  cov.model <- match.arg(cov.model, choices = CMB.cov.models)
  .checkCMB.cov.model(cov.model = cov.model,
                      cov.pars = c(sigmasq, phi),
                      kappa = kappa,
                      output = FALSE)
  x <- NULL; rm(x) # Trick R CMD Check
  graphics::curve(covmodelCMB(0, cov.pars = c(sigmasq, phi),
                              kappa = kappa,
                              cov.model = cov.model)
                  - covmodelCMB(x, cov.pars = c(sigmasq, phi),
                                kappa = kappa,
                                cov.model = cov.model),
                  from = from,
                  to = to,
                  xlab = "distance",
                  ylab = expression(gamma(h)),
                  lty = 2,
                  main = bquote(paste(.(cov.model)~"variogram"))
  )
}

#' Adds lines of fitted variograms to variogram plots
#'
#'
#'This function adds a line with the variogram model fitted by the function
#'\code{\link{variofitCMB}} to a current variogram plot. The function modifies
#'\code{\link[geoR]{lines.variomodel.variofit}} from the package \strong{geoR}
#'for additional covariance models on spheres.
#'
#'
#'@param x An object of the class \code{variofit} containing information about
#'the fitted model obtained as an output of the function \code{\link{variofitCMB}}.
#'@param max.dist A maximum distance to draw the variogram line. The default is
#'\code{x$max.dist}.
#'@param scaled logical. If TRUE the sill in the plot is 1.
#'@param ... other plotting parameters passed to \code{\link[graphics]{curve}}
#'
#'@details
#'The function adds a line with fitted variogram model to a plot. It is used
#'to compare empirical variograms against fitted models returned by
#'\code{\link{variofitCMB}}. #' Available models are: "matern", "exponential",
#'"spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric".
#'
#'@return
#'A line with a fitted variogram model is added to a plot.
#'
#'@references
#'\strong{geoR} package, \code{\link[geoR]{lines.variomodel.variofit}},
#'\code{\link{covmodelCMB}}, \code{\link{variofitCMB}}
#'
#'@examples
#' ## Plot the fitted Matern variogram versus its empirical variogram
#' #
#' # df <- CMBDataFrame("../CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 10000)
#' # varcmb <- variogramCMB(cmbdf, max.dist = 0.1, num.bins = 30)
#' # varcmb
#' # ols <- variofitCMB(varcmb,  fix.nug=FALSE, wei="equal", cov.model= "matern")
#' # plot(varcmb)
#' # lines(ols, lty=2)
#' #
#' ## Plot the fitted Askey variogram versus its empirical variogram
#' #
#' # ols <- variofitCMB(vario1, ini.cov.pars = c(1, 0.03), fix.nug = TRUE,
#' #     kappa = 3, wei = "equal", cov.model = "askey")
#' # plot(varcmb, main = ols$cov.model)
#' # linesCMB(ols, lty = 2)
#'
#'@export
linesCMB <-  function (x, max.dist, scaled = FALSE, ...) {
  my.l <- list()
  if (missing(max.dist)) {
    my.l$max.dist <- x$max.dist
    if (is.null(my.l$max.dist))
      stop("argument max.dist needed for this object")
  }
  else
    my.l$max.dist <- max.dist
  if (any(
    x$cov.model == c(
      "matern",
      "powered.exponential",
      "cauchy",
      "gencauchy",
      "gneiting.matern",
      "askey",
      "c2wendland",
      "c4wendland",
      "sinepower",
      "multiquadric"
    )
  ))
  my.l$kappa <- x$kappa
  else
    kappa <- NULL
  if (is.vector(x$cov.pars))
    my.l$sill.total <- x$nugget + x$cov.pars[1]
  else
    my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
  my.l$nugget <- x$nugget
  my.l$cov.pars <- x$cov.pars
  my.l$cov.model <- x$cov.model
  if (scaled) {
    if (is.vector(x$cov.model))
      my.l$cov.pars[1] <-  my.l$cov.pars[1] / my.l$sill.total
    else
      my.l$cov.pars[, 1] <-
        my.l$cov.cov.pars[, 1] / my.l$sill.total
    my.l$sill.total <- 1
  }
  gamma.f <- function(x, my.l) {
    if (any(my.l$cov.model == c("linear", "power")))
      return(my.l$nugget + my.l$cov.pars[1] * (x ^ my.l$cov.pars[2]))
    else
      return(
        my.l$sill.total -
          covmodelCMB(x,
            cov.model = my.l$cov.model,
            kappa = my.l$kappa,
            cov.pars = my.l$cov.pars
          )
      )
  }
  x <- NULL; rm(x) #Trick R CMD Check
  graphics::curve(gamma.f(x, my.l = my.l),
    from = 0,
    to = my.l$max.dist,
    add = TRUE,
    ...
  )
  return(invisible())
}

#' Practical range for covariance  function
#'
#' This function computes the practical range for covariance  functions on spheres.
#' The function modifies \code{\link[geoR]{practicalRange}} from the
#' package \strong{geoR} for additional covariance models on spheres.
#'
#'
#'@param cov.model A type of the correlation function. Available choices are: "matern",
#'"exponential","spherical", "powered.exponential", "cauchy", "gencauchy", "pure.nugget",
#'"askey", "c2wendland", "c4wendland", "sinepower", "multiquadric".
#'@param phi The range parameter as documented in \code{\link{covmodelCMB}}
#'@param kappa A smoothness parameter of the correlation function.
#'@param correlation A correlation threshold (default is 0.05)
#'@param ... other optimisation parameters
#'
#'@details
#' The practical(effective)  range for a covariance  function is the distance at which
#' a covariance function first time reaches the specified value \code{correlation}.  For
#' covariance functions on  spheres the practical range does not exceed \eqn{pi}, the
#' distance beyond which a covariance function is not defined. For the covariance
#' functions "spherical", "askey", "c2wendland", "c4wendland" their practical ranges
#' are equal to lengths of their support.
#'
#'@return
#'Value of the practical range for the covariance  function specified in \code{\link{covmodelCMB}}
#'
#'@references
#'\strong{geoR} package, \code{\link[geoR]{practicalRange}}, \code{\link{covmodelCMB}}
#'
#'@examples
#'
#'practicalRangeCMB(cov.model = "sinepower", phi = 0.1,  kappa = 0.5)
#'practicalRangeCMB(cov.model = "askey", phi = 0.1,  kappa = 0.5)
#'
#'@export
practicalRangeCMB <- function (cov.model,
            phi,
            kappa = 0.5,
            correlation = 0.05,...){
    cov.model <- match.arg(cov.model, choices = CMB.cov.models)
    .checkCMB.cov.model(
      cov.model = cov.model,
      cov.pars = c(1, phi),
      kappa = kappa,
      output = FALSE
    )
    if (cov.model %in% c("circular",
                         "cubic",
                         "spherical",
                         "askey",
                         "c2wendland",
                         "c4wendland"))
      return(min(phi, pi))
    if (any(cov.model %in% c("pure.nugget")))
      return(0)
    if (any(cov.model %in% c("linear", "sinepower", "multiquadric")))
      return(pi)
    if (any(cov.model %in% c("power")))
      return(Inf)
    findRange <- function(range, cm, p, k, cor){
      covmodelCMB (
        range,
        cov.model = cm,
        kappa = k,
        cov.pars = c(1, p)
      ) - cor
    }
    pr <- stats::uniroot(findRange,
        interval = c(0, 50 * phi + 1),
        cm = cov.model,
        p = phi,
        k = kappa,
        cor = correlation,...)$root
    return(min(pr, pi))
  }

CMB.cov.models <- c(
    "matern",
    "exponential",
    "spherical",
    "powered.exponential",
    "cauchy",
    "gencauchy",
    "pure.nugget",
    "askey",
    "c2wendland",
    "c4wendland",
    "sinepower",
    "multiquadric"
  )

# Helper function for variofitCMB
.loss1.vario <-  function (theta, g.l)
{
  if (g.l$cov.model == "linear")
    theta <- c(theta[1], 1, theta[-1])
  ##
  ## Imposing constraints for nlm
  ##
  if (g.l$m.f == "nlm") {
    assign(".temp.theta",  NULL, pos = 1)
    if (!g.l$fix.kappa) {
      if (g.l$fix.nugget) {
        if (g.l$cov.model == "power")
          theta.minimiser <- theta[1]
        else
          theta.minimiser <- theta[1:2]
        Tkappa <- theta[3]
      }
      else{
        if (g.l$cov.model == "power")
          theta.minimiser <- theta[c(1:3)]
        else
          theta.minimiser <- theta[1:3]
        Tkappa <- theta[4]
      }
    }
    else
      theta.minimiser <- theta
    penalty <- 10000 * sum(0 - pmin(theta.minimiser, 0))
    theta <- pmax(theta.minimiser, 0)
    if (!g.l$fix.kappa)
      theta <- c(theta.minimiser, Tkappa)
    if (any(theta.minimiser < 0))
      assign(".temp.theta", theta, pos = 1)
    else
      penalty <- 0
  }
  else
    penalty <- 0
  ##
  ## reading parameters
  ##
  if (!g.l$fix.kappa) {
    if (g.l$fix.nugget) {
      tausq <- g.l$nugget
      Tkappa <- theta[3]
    }
    else{
      tausq <- theta[3]
      Tkappa <- theta[4]
    }
    ## kappa now > 0 for both nlm() and optim()
    ##if(g.l$m.f == "nlm"){
    if (g.l$cov.model == "powered.exponential")
      kappa <-  2 * (exp(Tkappa)) / (1 + exp(Tkappa))
    else
      kappa <- exp(Tkappa)
    ##}
    ##else kappa <- Tkappa
  }
  else{
    kappa <- g.l$kappa
    if (g.l$fix.nugget)
      tausq <- g.l$nugget
    else
      tausq <- theta[3]
  }
  ##
  sigmasq <- theta[1]
  phi <- theta[2]
  if (g.l$cov.model == "power")
    phi <- 2 * exp(phi) / (1 + exp(phi))
  sill.total <- sigmasq + tausq
  ##
  ## Computing values for the theoretical variogram
  ##
  if (any(g.l$cov.model == c("linear", "power")))
    gammaU <- tausq + sigmasq * (g.l$u ^ phi)
  else
    gammaU <-
    sill.total - covmodelCMB(g.l$u,
      cov.model = g.l$cov.model,
      kappa = kappa,
      cov.pars = c(sigmasq, phi)
    )
  ##
  ## Computing loss function
  ##
  if (g.l$weight == "equal")
    loss <- sum((g.l$v - gammaU) ^ 2)
  if (g.l$weights == "npairs")
    loss <- sum(g.l$n * (g.l$v - gammaU) ^ 2)
  if (g.l$weights == "cressie")
    loss <- sum((g.l$n / (gammaU ^ 2)) * (g.l$v - gammaU) ^ 2)
  if (loss > (.Machine$double.xmax ^ 0.5) |
      loss == Inf | loss == -Inf | is.nan(loss))
    loss <- .Machine$double.xmax ^ 0.5
  return(loss + penalty)
}

# Helper function to check validity of cov.model
.checkCMB.cov.model <-
  function(cov.model,
           cov.pars,
           kappa,
           env = NULL,
           output = TRUE)
  {
    ## extracting covariance parameters
    if (is.vector(cov.pars))
      sigmasq <- cov.pars[1]
    else
      sigmasq <- cov.pars[, 1]
    if (is.vector(cov.pars))
      phi <-  cov.pars[2]
    else
      phi <- cov.pars[, 2]
    if (missing(kappa) || is.null(kappa))
      kappa <- NA
    ## checking for nested models
    cov.pars <- drop(cov.pars)
    if (is.vector(cov.pars))
      ns <- 1
    else{
      ns <- nrow(cov.pars)
      if (length(cov.model) == 1)
        cov.model <- rep(cov.model, ns)
      if (length(kappa) == 1)
        kappa <- rep(kappa, ns)
    }
    if (length(cov.model) != ns)
      stop("wrong length for cov.model")
    ##
    cov.model <- sapply(cov.model, match.arg, CMB.cov.models)
    cov.model[cov.model == "stable"] <- "powered.exponential"
    if (any(cov.model == c("gneiting.matern", "gencauchy"))) {
      if (length(kappa) != 2 * ns)
        stop(
          paste(
            "wrong length for kappa, ",
            cov.model,
            "model requires two values for the argument kappa"
          )
        )
    }
    else{
      if (length(kappa) != ns)
        stop('wrong length for kappa')
    }
    ## settings for power model (do not reverse order of the next two lines!)
    phi[cov.model == "linear"] <- 1
    cov.model[cov.model == "linear"] <- "power"
    ## checking input for cov. models with extra parameter(s)
    if (any(cov.model == 'gneiting.matern') && ns > 1)
      stop('nested models including the gneiting.matern are not implemented')
    for (i in 1:ns) {
      if (any(
        cov.model[i] == c(
          "matern",
          "powered.exponential",
          "cauchy",
          "gneiting.matern",
          "gencauchy",
          "askey",
          "c2wendland",
          "c4wendland",
          "sinepower",
          "multiquadric"
        )
      )) {
        if (any(cov.model[i] == c("gneiting.matern", "gencauchy"))) {
          if (any(is.na(kappa)) || length(kappa) != 2 * ns)
            stop(
              paste(
                cov.model[i],
                "correlation function model requires a vector with 2 parameters in the argument kappa"
              )
            )
        }
        else{
          if (is.na(kappa[i]) | is.null(kappa[i]))
            stop(
              "for matern, powered.exponential, gencauchy, askey, c2wendland,
              c4wendland, sinepower, multiquadric and cauchy covariance
              functions the parameter kappa must be provided"
            )
        }
        if ((
          cov.model[i] == "matern" | cov.model[i] == "powered.exponential" |
          cov.model[i] == "askey" | cov.model[i] == "c2wendland" |
          cov.model[i] == "c4wendland" |
          cov.model[i] == "sinepower" |
          cov.model[i] == "multiquadric"
        ) & length(kappa) != 1 * ns)
          stop("kappa must have 1 parameter for this correlation function")
        if (cov.model[i] == "matern" &
            kappa[i] == 0.5)
          cov.model[i] == "exponential"
        }
      if (cov.model[i] == "power")
        if (any(phi[i] >= 2) | any(phi[i] <= 0))
          stop("for power model the phi parameters must be in the interval ]0,2[")
      }
    if (!is.null(env)) {
      assign("sigmasq", sigmasq, envir = env)
      assign("phi", phi, envir = env)
      assign("kappa", kappa, envir = env)
      assign("ns", ns, envir = env)
      assign("cov.model", cov.model, envir = env)
    }
    if (output)
      return(list(
        cov.model = cov.model,
        sigmasq = sigmasq,
        phi = phi,
        kappa = kappa,
        ns = ns
      ))
    else
      return(invisible())
  }
