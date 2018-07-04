#' Nested Search
#'
#' Finds the closest HEALPix pixel center to a given \code{target} point,
#' specified in Cartesian coordinates, using an efficient nested search
#' algorithm. HEALPix indices are all assumed to be in the "nested"
#' ordering scheme.
#'
#' @param target is a vector of Cartesian coordinates for the target
#' point on S^2
#' @param nside is the nside for which the HEALPix points are searched
#' @param demo.plot If TRUE then a plot will be produced with
#' target pixel in yellow and closest pixel at each step in red
#'
#' @return if \code{index.only = TRUE} then the output will be a HEALPix index.
#' If \code{index.only} FALSE then the output is the list containing the HEALPix index
#' and Cartesian coordinate vector of the HEALPix point closest to \code{target}.
#'
#' @examples
#' # Find the pix index and Cartesian coordinates of the HEALPix point
#' # at nside closest to the target point c(0,0,1)
#' h <- nestSearch(c(0,0,1), nside=1024)
#' cat("Closest HEALPix point to (0,0,1) at nside = 1024 is (",h$xyz,")")
#'
#' @export
nestSearch <- function(target, nside,
                       index.only = FALSE,
                       j = 0:log2(nside),
                       demo.plot = FALSE)
{
  j <- c(0,j)
  h <- list(pix = 0)
  for ( i in 2:length(j) )
  {
    h <- rcosmo::nestSearch_step( target, j2=j[i],
                                  j1=j[i-1], pix.j1 = h$pix,
                                  demo.plot = demo.plot)
  }

  if (index.only==TRUE) {
    return(h$pix)
  }
  else {
    return(list( xyz=h$xyz, pix=h$pix ))
  }
}




#' nestSearch_step
#'
#' Search for the closest HEALPix pixel to a \code{target} point,
#' where the search is restricted to within HEALPix pixel,
#' \code{pix.j1}, at resolution j1. The returned value is a
#' HEALPix pixel (and, optionally, the cartesian coordinates of
#' its center) at resolution j2, where j2 > j1. All pixels are
#' assumed to be in nested ordering scheme.
#'
#' j1 and j2 are HEALPix resolution parameters, i.e., \eqn{nside = 2^j}.
#'
#' \code{nestSearch_step(target, j2, j1, pix.j1)} searches within the subregion
#' pix.j1, where pix.j1 is a HEALPix pixel index at resolution j1.
#' The return value is the HEALPix point closest to \code{target}, at
#' resolution j2.
#'
#' Setting \code{pix.j1 = 0} (the default) searches for the HEALPix point closest
#' to \code{target} at resolution j2, among all HEALPix points at resolution j1.
#'
#' @param target is the target point on S^2 in spherical coordinates.
#' @param j1 is the lower resolution, with j1 < j2.
#' @param j2 is the upper resolution.
#' @param pix.j1 is the initial pix index at resolution j1,
#' i.e., the j1-level pixel to search in. If \code{pix.j1 = 0} then
#' all pixels will be searched (slow).
#' @param demo.plot If TRUE then a plot will be produced with
#' target pixel in yellow and closest pixel in red
#'
#' @return A list containing the Cartesian coordinates, \code{xyz},
#' and the HEALPix pixel index, \code{pix}, of the closest HEALPix
#' pixel center to the target point, \code{target}, at resolution j2
#'
#' @examples
#' # search for the HEALPix pixel center closest to North pole
#' # (0,0,1) at level 3
#' nestSearch_step(target = c(0,0,1), j2 = 3, j1 = -1, demo.plot = TRUE )
#'
#' @export
nestSearch_step <- function(target, j1 = j2, j2, pix.j1 = 0,
                            demo.plot = FALSE) {

  if ( j1 > j2 ) stop("j1 must be less than or equal to j2")

  # nside at level j2
  nside.j2 <- 2^j2

  spix.j2 <- pixelWindow(j1, j2, pix.j1)

  # Convert spix.j2 to Cartesian coordinates
  xyz.j2 <- rcosmo::pix2coords_internal(nside = nside.j2,
                                        nested = TRUE,
                                        cartesian = TRUE,
                                        spix = spix.j2)[,1:3]

  dists <-  apply(xyz.j2, MARGIN = 1,
                  function(point) {
                    sum((point - target)^2)
                  } )
  index.min <- which.min(dists)


  if ( demo.plot ){
    ### -------- JUST FOR VISUALISATION ---------- ###
    ###      PLOT TO VISUALISE THE RESULTS         ###
    # HEALPix points at level j2
    hp.j2 <- rcosmo::pix2coords(nside = nside.j2,
                                order = "nested",
                                coords = "cartesian")
    rgl::plot3d(hp.j2, col="black", size = 2,
                type = 'p', pch = 2,  add = TRUE)

    # HEALPix points in window
    pixelWin <- rcosmo::pixelWindow(j1, j2, pix.j1)
    pixelWin <- rcosmo::pix2coords(nside = nside.j2,
                                   order = "nested",
                                   coords = "cartesian",
                                   spix = pixelWin)
    rgl::plot3d(pixelWin, col="green", size = 12,
                type = 'p', pch = 2,  add = TRUE)

    # # add the target point target
    rgl::plot3d(target[1], target[2], target[3],
                col="yellow", size = 12,
                type = 'p', pch = 2,  add = TRUE)
    #
    # # add the closest HEALPix point to target
    xyz <- as.numeric(xyz.j2[index.min,])
    rgl::plot3d(xyz[1], xyz[2], xyz[3],
                col="red", size = 13,
                type = 'p', pch = 2,  add = TRUE)
    ### ----------------------------------------- ###
  }

  return( list(xyz = as.numeric( xyz.j2[index.min,] ),
               pix = spix.j2[index.min] ) )
}


#' Pixel window
#'
#' All pixels are assumed to be in nested ordering
#'
#'@param j1 is the lower resolution, with j1 < j2
#'@param j2 the upper resolution
#'@param pix.j1 the pixel index at resolution j1 within which
#'all pixels from resolution j2 will be returned. \code{pix.j1} can
#'also be a vector of non-zero pixel indices.
#'
#'@return All pixels in resolution j2 that fall within the pixel
#'pix.j1 specified at resolution j1
#'
#'@export
pixelWindow <- function(j1, j2, pix.j1)
{
  if ( length(pix.j1) == 1 && pix.j1 == 0 ) {
    # pix indices at level j2
    spix.j2 <- 1:(12*2^(2*j2)) #= 12*nside.j2^2

  } else {
  # Number of pixels at level j2 in each pixel from level j1
    lev.diff <- 2^((j2-j1)*2)
    # pix indices at level j2
    spix.j2 <- unlist(mapply(seq, from = (lev.diff*(pix.j1-1)+1),
                             to = (lev.diff*pix.j1),
                             SIMPLIFY = FALSE))
  }
  return(spix.j2)
}
