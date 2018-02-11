#' nest_search
#'
#' Search for the closest HEALPix pixel to a \code{target} point,
#' where the search is restricted to within HEALPix pixel,
#' \code{pix.j1}, at resolution j1. The returned value is a
#' HEALPix pixel (and, optionally, the cartesian coordinates of
#' its center) at resolution j2, where j2 > j1.
#'
#' j1 and j2 are HEALPix resolution parameters, i.e., \eqn{Nside = 2^j}.
#'
#' \code{nest_search(target, j2, j1, pix.j1)} searches within the subregion
#' pix.j1, where pix.j1 is a HEALPix pixel index at resolution j1.
#' The return value is the HEALPix point closest to \code{target}, at
#' resolution j2.
#'
#' \code{nest_search(target, j2)} searches for the HEALPix point closest to
#' \code{target} at resolution j2, among all HEALPix points at resolution j2.
#'
#' @param target is the target point on S^2 in spherical coordinates.
#' @param j2 is the upper resolution.
#' @param j1 is the lower resolution, with j1 < j2.
#' @param pix.j1 is the initial pix index at resolution j1,
#' i.e., the j1-level pixel to search in. If \code{pix.j1 = 0} then
#' all pixels will be searched (slow).
#'
#' @return A list containing the Cartesian coordinates, \code{xyz},
#' and the HEALPix pixel index, \code{pix}, of the closest HEALPix
#' pixel center to the target point, \code{target}, at resolution j2
#'
#' @examples
#' # search for the HEALPix pixel center closest to North pole
#' # (0,0,1) at level 3
#' nest_search(target = c(0,0,1), j2 = 3, j1 = -1, plot_points=TRUE )
#'
#' @export
nest_search <- function(target, j2, j1 = 0, pix.j1 = 0) {

  ## TEST VALUES
  # library(sphereplot)
  # target <- sph2car(0.2, 0.4, deg = FALSE)
  # j2 <- log2(16)
  # j1 <- 0
  # pix.j1 <- 0

  if ( j1 > j2 ) stop("j1 must be less than or equal to j2")

  # Nside at level j2
  nside.j2 <- 2^j2

  if ( j1 <= 0 || pix.j1 == 0 ) {
    pix.j2 <- 1:(12*nside.j2^2)
  } else {
    # number of points at level j2 in each box at level j1
    lev.diff <- 2^((j2-j1)*2)
    # pix indices at level j2
    pix.j2 <- (lev.diff*(pix.j1-1)+1):(lev.diff*pix.j1)
  }

  # Convert pix.j2 to Cartesian coordinates
  xyz.j2 <- pix2vec(Nside = nside.j2, order = "ring", spix = pix.j2)

  ### DOES THE SAME THING AS THE ABOVE, BUT 15 TIMES FASTER
  dists <-  apply(xyz.j2, MARGIN = 1,
                  function(point) {
                    sum((point - target)^2)
                  } )
  index.min <- which.min(dists)

  ### -------- JUST FOR VISUALISATION ---------- ###
  ###      PLOT TO VISUALISE THE RESULTS         ###
  # HEALPix points at level j2
  # hp.j2 <- pix2vec(Nside = nside.j2, order = "ring")
  # rgl::plot3d(hp.j2, col="blue", size = 3,
  #             type = 'p', pch = 2,  add = TRUE)
  #
  # # add the target point target
  # rgl::plot3d(target[1], target[2], target[3],
  #             col="yellow", size = 5,
  #             type = 'p', pch = 2,  add = TRUE)
  #
  # # add the closest HEALPix point to target
  # xyz <- as.numeric(xyz.j2[index.min,])
  # rgl::plot3d(xyz[1], xyz[2], xyz[3],
  #             col="red", size = 6,
  #             type = 'p', pch = 2,  add = TRUE)
  ### ----------------------------------------- ###

  return( list(xyz = as.numeric( xyz.j2[index.min,] ),
               pix = pix.j2[index.min] ) )
}

