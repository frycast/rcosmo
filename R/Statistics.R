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
#' An object of the class CMBcovariance that is a modification of \code{\link[geoR]{variog}}
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

  covs <- covCMB_internal2(cmbdf[, c("x", "y", "z", "I")], breaks)

  # Reverse order since cosine is decreasing
  covs[2:nrow(covs), 1] <- rev(covs[2:nrow(covs), 1])
  covs[2:nrow(covs), 2] <- rev(covs[2:nrow(covs), 2])
  v <- c(0, rev(acos(breaks)))
  # Drop the throw-away bin (distances greater than max.dist)
  covs <- covs[-nrow(covs), ]

  # Find the bin centers (break_{i+1} - break_i)/2 with break_0 = 0.
  centers <-
    c(0, v[1:(length(v) - 1)] + (v[2:length(v)] - v[1:(length(v) - 1)]) / 2)

  pairs.min <- 2
  n <- as.integer(c(covs[, 2][1], covs[, 2][-1] / 2))
  indp <- (n >= pairs.min)
  sd1 <- n
  result <- list(
    u = centers,
    v = covs[, 1],
    n = n,
    sd1,
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
  oldClass(result) <- "CMBcovariance"
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
#' An object of the class CMBcorrelation that is a modification of \code{\link[geoR]{variog}}
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
                         calc.max.dist = FALSE)
{corrCMB<- covCMB(cmbdf, num.bins,
                 sample.size,
                 max.dist,
                 breaks,
                 equiareal ,
                 calc.max.dist)
corrCMB$v <- corrCMB$v/corrCMB$v[1]
oldClass(corrCMB) <- "CMBcorrelation"
return(corrCMB)}

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
                   calc.max.dist = FALSE)
{varCMB<- covCMB(cmbdf, num.bins ,
                     sample.size,
                     max.dist ,
                     breaks,
                     equiareal ,
                     calc.max.dist )
varCMB$v <- varCMB$v[1]-varCMB$v
oldClass(varCMB) <- "variogram"
return(varCMB)}

#'Plot variogram
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
#' ## Download the map first
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # cmbdf <- sampleCMB(df, sample.size = 100000)
#' # varcmb <- variogramCMB(cmbdf, max.dist = 0.1, num.bins = 30, sample.size=1000)
#' # plot(varcmb)
#'
#'@name plot.variogram
#'
#' @aliases plot.variogram Plot variogram
#'
NULL

#'Plot CMBcovariance
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
plot.CMBcovariance <-  function (x, ...) {
      x0 <- x
      attributes(x0)$class <- "variogram"
      graphics::plot(x0, ylab = "sample covariance", ...)
  }

#'Plot CMBcorrelation
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
plot.CMBcorrelation <-  function (x, ...) {
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



#'Plot angular scaterplots and means
#'
#'For specified measurements from \code{\link{CMBDataFrame}} this function
#'produces scaterplots and binned means versus theta and phi angles.
#'
#'@param cmbdf A  full \code{\link{CMBDataFrame}} or a windowed
#'\code{\link{CMBDataFrame}}
#'
#'@param intensities  A CMBDataFrame column with measured values
#'
#'@return
#' 2x2 plot. The first row shows scaterplots. The second row gives
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
#' This funcion is a modification of standard \link{qqplot} functions to work
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
#' This funcion is a modification of standard \link{qqnorm} functions to work
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
#'region. The functions employes the function \link{entropy} and uses histogram
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
#'\code{\link{CMBDataFrame}}. The functions employes the function \link{chi2.empirical} and uses histogram
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
#'\code{\link{CMBWindow}}s.  For smal sample sizes and many zero counts
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

