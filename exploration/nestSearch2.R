
# Parallel nestSearch ----------------------------------------------------

target1 <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
target2 <- c(1,0,0)
target3 <- matrix(c(1,0,0,0,1,0,0,0,1), byrow = TRUE, nrow = 3)
ns <- 16


nestSearch2(target1, nside = ns)
nestSearch2(target2, nside = ns)
nestSearch2(target3, nside = ns)



world <- read.csv("../worldcities.csv")
df <- data.frame(phi = pi/180*world$lng, theta = pi/2 - pi/180*world$lat,
                   I=rep(1,length(world$lng)))
k <- 5000
dfs <- df[sample(nrow(df), k), ]
separatingNside(dfs)
a <- nestSearch2(dfs, nside = 2^12)



# This version always assumes that j2 = j1 + 1
# Takes advantage of children being resolution independent
nestSearch_step2 <- function(target, j2, pix.j1 = 0) {


  # Get the 4 kids
  if (is.list(pix.j1)) { kids <-  lapply(pix.j1, children)
    } else { kids <-  list(children(pix.j1)) }


  nside.j2 <- 2^j2
  # Coordinates of the 4 kids
  xyz.j2 <- lapply(kids, function(x) {
    rcosmo:::pix2coords_internal(nside = nside.j2,
                                 nested = TRUE,
                                 cartesian = TRUE,
                                 spix = x) })


  tlen <- length(target)
  result <- list()
  for (i in 1:tlen) {
    dots <- xyz.j2[[i]] %*% target[[i]]
    minkid <- kids[[i]][max.col(t(dots), ties.method = "first")]
    result[[i]] <- rcosmo::neighbours(minkid, j2)
  }


  return( result )
}



## This version breaks target into a list of vectors
# target may be a data.frame, matrix, or a single vector
nestSearch2 <- function(target, nside,
                       index.only = FALSE) {

  # # Convert the target to a list where elements are the row vectors
  # if ( is.numeric(target) && !is.matrix(target) ) { target <- list(target)
  # } else if ( is.matrix(target) ) { target <- as.list(as.data.frame(t(target))) }

  if ( is.matrix(target) ) {

    target <- as.list(as.data.frame(t(target)))
  } else if ( is.data.frame(target) ) {

    if ( all(c("theta","phi") %in% names(target)) ) {

      coords(target) <- "cartesian"
      target <- target[,c("x","y","z")]
    }
    target <- as.list(as.data.frame(t(target)))
  } else if ( is.numeric(target) ) {

    target <- list(target)
  } else {

    stop("Target must be data.frame, matrix or numeric vector")
  }

  j = 0:(log2(nside)+1)
  jlen <- length(j)
  tlen <- length(target)
  h <- rep(list(0), tlen)
  for ( i in 1:jlen ) {

    h <- nestSearch_step2( target, j2 = j[i], pix.j1 = h)
  }

  # Note h is now one level deeper than the target resolution.
  # Convert h to Cartesian coordinates.
  h.xyz <- lapply(h, function(x) {
  rcosmo:::pix2coords_internal(nside = 2^j[jlen],
                               nested = TRUE,
                               cartesian = TRUE,
                               spix = x) })

  # Get the parent of the closest
  result.h <- list()
  for (i in 1:tlen) {
    dots <- h.xyz[[i]] %*% target[[i]]
    min.h <- h[[i]][max.col(t(dots), ties.method = "first")]
    result.h[[i]] <- parent(min.h)
  }

  h <- unlist(result.h)

  if ( !index.only ) {
    xyz <- rcosmo:::pix2coords_internal(nside = nside,
                                 nested = TRUE,
                                 cartesian = TRUE,
                                 spix = h)
    return(list(xyz = xyz, pix = h))
  }

  return(h)
}









# PREVIOUS ----------------------------------------------------------------



#' parent
#'
#' Gives the pixel at resolution j - 1 that contains p,
#' where p is specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
parent <- function(p)
{
  (p - p %% 4 + (p %% 4 != 0)*4)/4
}

#' children
#'
#' Gives the pixels at resolution j + 1 that are contained in p,
#' where p is specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
children <- function(p)
{
  if ( any(p > 0) ) {
    1:4 + rep((p-1)*4, each = 4)
  } else { 1:12 }
}

#' siblings
#'
#' The siblings of pixel p are defined as the
#' children of the parent of p. Note this is resolution independent.
#'
#' @p Pixel index in nested order.
#'
siblings <- function(p) {
  h <- (p - p %% 4 + (p %% 4 != 0)*4)/4
  1:4 + (h-1)*4
}

#' displayPixels
#'
#' Display the pixels spix at resolution j by colouring
#' in the grandchildren of spix at resolution plot.j
#'
#' @param j The resolution that spix are specified at.
#' @param boundary.j The resolution to display boundaries at. If
#' this is missing then boundaries will not be plotted.
#' @param plot.j The resolution to plot grandchildren at
#' @param spix Integer vector. The pixel indices to display.
#' These must be in nested order.
#' @param incl.labels Integer vector of pixel indices to label at
#' resolution j.
#' @param boundary.col The boundary colour.
#' @param boundary.lwd The boundary line width.
#' @param col The colour to make the grandchildren.
#' @param size The size to make the grandchildren.
#'
#'
#'@examples
#'
#' demoNeighbours <- function(p,j) {
#'   neighbours(p, j)
#'   displayPixels(boundary.j = j, j = j, plot.j = 5,
#'                 spix = neighbours(p, j),
#'                 boundary.col = "gray",
#'                 boundary.lwd = 1,
#'                 incl.labels = neighbours(p, j),
#'                 col = "blue",
#'                 size = 3)
#'   rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
#' }
#'
#'
#'
#'
displayPixels <- function(boundary.j, j, plot.j = 5, spix,
                          boundary.col = "gray",
                          boundary.lwd = 1,
                          incl.labels = 1:(12*4^boundary.j),
                          col = "blue",
                          size = 3)
{
  if ( !missing(boundary.j) ) {
    rcosmo::displayPixelBoundaries(nside = 2^boundary.j,
                             ordering = "nested",
                             nums.col = "red",
                             col = boundary.col,
                             lwd = boundary.lwd,
                             incl.labels = incl.labels)
  }

  # We do this by plotting grandchildren of the siblings
  gchild <- rcosmo::pixelWindow(j1 = j,
                                j2 = plot.j,
                                pix.j1 = spix)

  hp <- rcosmo::HPDataFrame(nside = 2^plot.j,
                            spix = gchild)
  plot(hp, add = TRUE, col = col, size = size)
}




# We test siblings by colouring the siblings of some pixels at level j
displaySiblings <- function(p, j, boundary.j = j,
                            plot.j = 5, col = "blue",
                            size = 3,
                            label.siblings = TRUE)
{
  spix <- siblings(p)
  if (label.siblings)
  {
    labels <- siblings(p)
  }
  else
  {
    labels <- 0
  }
  displayPixels(boundary.j = boundary.j, j = j,
                plot.j = plot.j, spix = spix, col = col, size = size,
                incl.labels = labels)
}
# Looks good:
# displaySiblings(4, 1)
# displaySiblings(4, 2)
# displaySiblings(27, 1)
# displaySiblings(28, 1)
# displaySiblings(28, 2)
# displaySiblings(127, 3)
# displaySiblings(4, 1)


#' baseSiblings
#'
#' A map from the base pixel index bp to the vector of base pixels
#' that are neighbours of bp, in counterclockwise order of
#' direction: S,SE,E,NE,N,NW,W,SW. The presence of -1 indicates
#' that the corresponding direction is empty.
#'
#' @param bp The base pixel index
#'
baseSiblings <- function(bp)
{
  # order: S,SE,E,NE,N,NW,W,SW
  # corners: S,E,N,W
  switch(bp,
         # north pole
         c(9 ,6,-1,2,3,4,-1,5),
         c(10,7,-1,3,4,1,-1,6),
         c(11,8,-1,4,1,2,-1,7),
         c(12,5,-1,1,2,3,-1,8),
         # equatorial
         c(-1,9 ,6,1,-1,4,8,12),
         c(-1,10,7,2,-1,1,5,9),
         c(-1,11,8,3,-1,2,6,10),
         c(-1,12,5,4,-1,3,7,11),
         # south pole
         c(11,10,-1,6,1,5,-1,12),
         c(12,11,-1,7,2,6,-1,9),
         c(9 ,12,-1,8,3,7,-1,10),
         c(10,9 ,-1,5,4,8,-1,11))
}


#' p2ibp
#'
#' Convert a pixel index p to its index within
#' the base pixel to which p belongs
#'
#' @param p The pixel index at resolution j, in nested order.
#' @param j The resolution parameter nside = 2^j
#'
p2ibp <- function(p, j) #indexInBP
{
  (p-1) %% 4^j + 1
}


#' p2bp
#'
#' The base pixel to which pixel p belongs at resolution j
#'
#' @param p The pixel index at resolution j, in nested order.
#' @param j The resolution parameter nside = 2^j
#'
p2bp <- function(p, j)
{
  floor((p-1) / (4^j)) + 1
}

# Find the pixel index p of a given pixel at ibp in bp
ibp2p <- function(ibp, bp, j)
{
  (bp - 1)*4^j + ibp
}

# # Note that for all p we have
# alwaysTrue <- function(p, j) {
#   all((p2bp(p, j) - 1)*4^j + p2ibp(p, j) == p)
# }
# # Check:
# # j <- 3
# # all(alwaysTrue(1:(12*2^j^2), j))

# Convert binary to decimal
bin2dec <- function(x, digits)
{
  pow <- 2^(0:31)[1:digits]
  sum(pow[as.logical(x)])
}

# Convert decimal to binary
dec2bin = function(number, digits) {
  as.numeric(intToBits(number))[1:digits]
}

# Separate a binary number (e.g., output of dec2bin)
# into its even and odd digits
bin2f <- function(bin, j)
{
  even.bits <- bin[seq(2,2*j  , by = 2)]
  odd.bits  <- bin[seq(1,2*j-1, by = 2)]

  return(list(even = even.bits, odd = odd.bits))
}

# Recombine even and odd digits into a binary number
f2bin <- function(f, j)
{
  bin <- vector(mode = "integer", length = 2*j)
  bin[seq(2,2*j  , by = 2)] <- f$even
  bin[seq(1,2*j-1, by = 2)] <- f$odd
  return(bin)
}


#' onBPBoundary
#'
#' Check if the pixel p at resolution j is on the boundary of a base pixel,
#' note \code{nside = 2^j}. If \code{p} is on the boundary then a number
#' is returned corresponding to the boundary position,
#' in order from 1 to 8 this is S, SE, E, NE, N, NW, W, SW,
#' where S,E,N,W are corners and SE, NE, NW, SW are edges.
#' If \code{p} is not on the boundary then 0 is returned.
#'
#' Pixel p is on the boundary of a base pixel if and only if
#' the binary representation of \eqn{p - 1} has either
#' its even bits or its odd bits (or both) all zeros or all ones.
#'
#' @param p Integer vector. The pixel index (or indices) at resolution j,
#' in nested order.
#' @param j The resolution within which to specify p.
#'
#' @return
#' If \code{p} is on the boundary then a number
#' is returned corresponding to the boundary position,
#' in order from 1 to 8 this is S, SE, E, NE, N, NW, W, SW,
#' where S,E,N,W are corners and SE, NE, NW, SW are edges.
#' If \code{p} is not on the boundary then 0 is returned.
#'
#'@examples
#'
#' # Visualise to see which are on boundary
#' j <- 2
#' spix <- 1:16
#' rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
#' displayPixels(boundary.j = j, j = j, spix = spix,
#'               boundary.col = "gray", incl.labels = spix)
#'
#' # Corner pixels counter-clockwise S, E, N, W
#' onBPBoundary(c(1,6,16,11),j)
#'
#' # Edge pixels counter-clockwise SE, NE, NW, SW
#' onBPBoundary(c(2,8,16,9),j)
#'
#' # Non-boundary pixels
#' onBPBoundary(c(4,7,10,13),j)
#'
onBPBoundary_old <- function(p,j)
{
  result <- integer(1)
  # The number of bits in the binary representation is 2*j (giving 4^j pixels).
  vapply(p, function(x) {
    bin <- dec2bin(p2ibp(x, j) - 1, digits = 2*j)
    f <- bin2f(bin, j = j)
    se <- sum( f$even )
    so <- sum( f$odd )

    b <- 0L
    # south corner
    if ( se == 0 && so == 0 )
    {
      b <- 1L
    }
    # east corner
    else if ( se == 0 && so == j )
    {
      b <- 3L
    }
    # north corner
    else if ( se == j && so == j )
    {
      b <- 5L
    }
    # west corner
    else if ( se == j && so == 0 )
    {
      b <- 7L
    }
    # south east wall
    else if ( se == 0 )
    {
      b <- 2L
    }
    # north east wall
    else if ( so == j )
    {
      b <- 4L
    }
    # north west wall
    else if ( se == j )
    {
      b <- 6L
    }
    # south west wall
    else if ( so == 0 )
    {
      b <- 8L
    }
    return(b)
  }, FUN.VALUE = result)
}



#' a version of onBPBoundary to use with neighbours
#'
#' @param se The sum of even bits, e.g. sum( f$even )
#' @param so The sum of odd bits, e.g. sum( f$odd )
#' @param j The resolution parameter nside = 2^j
#'
onBPBoundary <- function(se, so, j)
{
    b <- 0L
    # south corner
    if ( se == 0 && so == 0 )
    {
      b <- 1L
    }
    # east corner
    else if ( se == 0 && so == j )
    {
      b <- 3L
    }
    # north corner
    else if ( se == j && so == j )
    {
      b <- 5L
    }
    # west corner
    else if ( se == j && so == 0 )
    {
      b <- 7L
    }
    # south east wall
    else if ( se == 0 )
    {
      b <- 2L
    }
    # north east wall
    else if ( so == j )
    {
      b <- 4L
    }
    # north west wall
    else if ( se == j )
    {
      b <- 6L
    }
    # south west wall
    else if ( so == 0 )
    {
      b <- 8L
    }
    return(b)
}



#' Border pattern is resolution independent
#'
#' @param pype is the output of onBPBoundary
#'
#' @return the output is useful as an index
#' to the output of baseSiblings. It will
#' return the correct neighbouring base
#' pixel in each direction. The
#' 0 indicates to stay in current BP.
#'
borderPattern <- function(ptype)
{
  switch(ptype + 1,
         # S SE  E NE  N NW  W SW
         c(0, 0, 0, 0, 0, 0, 0, 0),
         c(1, 2, 2, 0, 0, 0, 8, 8),
         c(2, 2, 2, 0, 0, 0, 0, 0),
         c(2, 2, 3, 4, 4, 0, 0, 0),
         c(0, 0, 4, 4, 4, 0, 0, 0),
         c(0, 0, 4, 4, 5, 6, 6, 0),
         c(0, 0, 0, 0, 6, 6, 6, 0),
         c(8, 0, 0, 0, 6, 6, 7, 8),
         c(8, 0, 0, 0, 0, 0, 8, 8))
}


neighbours <- function(p, j)
{
  if ( j == 0 )
  {
    bs <- baseSiblings(p)
    return(bs[bs > 0])
  }

  # Get the index in BP and the BP
  ibp <- p2ibp(p, j)
  bp <- p2bp(p, j)

  # Get even and odd binary digit information
  f <- bin2f(dec2bin(ibp-1, digits = 2*j), j = j)
  se <- sum( f$even )
  so <- sum( f$odd )

  # Get the BP border crossing info: target border pixels
  # bdr is the direction of crossing:
  # 0,1,2,3,4,5,6,7,8 = none, S, SE, E, NE, N, NW, W, SW
  bdr <- borderPattern(onBPBoundary(se, so, j))
  target.bp <- rep(bp, 9)
  target.bp[which(bdr != 0)] <- baseSiblings(bp)[bdr]

  # Separately increment/decrement the odd/even binary reps
  even.dec <- bin2dec(f$even, digits = j)
  odd.dec <- bin2dec(f$odd, digits = j)
  # Order: S,SE, E,NE, N,NW, W,SW,self
  ei <- c(-1,-1,-1, 0, 1, 1, 1, 0,   0)
  oi <- c(-1, 0, 1, 1, 1, 0,-1,-1,   0)
  nbrs <- data.frame(even = I(lapply(rep(even.dec,9) + ei,
                                     dec2bin, digits = j)),
                     odd  = I(lapply(rep(odd.dec, 9) + oi,
                                     dec2bin, digits = j)))

  # Swap some nbrs N<->S depending on region of bp
  if ( bp %in% c(1,2,3,4) ) {
    # North pole, swap S->N in directions 4, 5 and 6 (NE, N, NW)
    i4 <- which(bdr == 4)
    i5 <- which(bdr == 5)
    i6 <- which(bdr == 6)

    #SW->NW
    if (length(i4)!=0) {
      nbrs[i4,]$odd  <- nbrs[i4,]$even
      nbrs[i4,]$even <- list(rep(1,j))
    }
    #S->N
    if (length(i5) != 0) {
      nbrs[i5,]$even <- list(rep(1,j))
      nbrs[i5,]$odd  <- list(rep(1,j))
    }
    #SE->NE
    if (length(i6)!=0) {
      nbrs[i6,]$even <- nbrs[i6,]$odd
      nbrs[i6,]$odd  <- list(rep(1,j))
    }

  } else if ( bp %in% c(9,10,11,12) ) {
    # South pole, swap N->S in directions 1, 2 and 8 (S, SE, SW)
    i1 <- which(bdr == 1)
    i2 <- which(bdr == 2)
    i8 <- which(bdr == 8)

    #N->S
    if (length(i1) != 0) {
      nbrs[i1,]$even <- list(rep(0,j))
      nbrs[i1,]$odd <- list(rep(0,j))
    }
    #NW->SW
    if (length(i2) != 0) {
      nbrs[i2,]$even <- nbrs[i2,]$odd
      nbrs[i2,]$odd <- list(rep(0,j))
    }
    #NE->SE
    if (length(i8) != 0) {
      nbrs[i8,]$odd <- nbrs[i8,]$even
      nbrs[i8,]$even <- list(rep(0,j))
    }
  }

  # Recombine even and odd
  nbrs <- recombineEvenOdd(nbrs, j)

  # Convert indices in BP to actual indices p
  np <- ibp2p(nbrs, target.bp, j)
  return(np[np > 0])
}

# Takes a data.frame with column for even binary and column for odd.
# Returns the recombined digits in decimal as vector. This is
# a helper function for neighbours.
recombineEvenOdd <- function(nbrs, j)
{
  unlist(
    mapply(FUN = function(even,odd) {
      bin <- f2bin(list(even = even, odd = odd), j = j)
      p <- bin2dec(bin, digits = 2*j) + 1
      return(p)
    },
    nbrs$even, nbrs$odd, SIMPLIFY = FALSE, USE.NAMES = FALSE))
}


demoNeighbours <- function(p,j) {
  neighbours(p, j)
  displayPixels(boundary.j = j, j = j, plot.j = 5,
                spix = neighbours(p, j),
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = neighbours(p, j),
                col = "blue",
                size = 3)
  rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
}


# demoNeighbours(1, 2)
# demoNeighbours(6, 2)
# demoNeighbours(16, 2)
# demoNeighbours(11, 2)
# demoNeighbours(172, 3)


nestSearch2_step <- function(target, j1 = j2, j2, pix.j1 = 0) {

  if ( j1 > j2 ) stop("j1 must be less than or equal to j2")

  # nside at level j2
  nside.j2 <- 2^j2

  spix.j2 <- rcosmo:::pixelWindow(j1, j2, pix.j1)

  # Convert spix.j2 to Cartesian coordinates
  xyz.j2 <- rcosmo:::pix2coords_internal(nside = nside.j2,
                                         nested = TRUE,
                                         cartesian = TRUE,
                                         spix = spix.j2)[,1:3]

  dists <-  apply(xyz.j2, MARGIN = 1,
                  function(point) {
                    sum((point - target)^2)
                  } )
  index.min <- which.min(dists)

  return( neighbours(spix.j2[index.min], j2) )
}

# It is advised that j should start at 0
nestSearch2 <- function(target, nside,
                        index.only = FALSE,
                        j = 0:log2(nside))
{
  j <- c(j[1],j)
  h <- 0
  for ( i in 2:length(j) )
  {
    h <- nestSearch2_step( target, j2 = j[i],
                           j1 = j[i-1], pix.j1 = h)
  }

  return(h)
}

# j <- 3
# spix <- nestSearch2(c(0,0,1), nside = 2^3)
# displayPixels(boundary.j = j, j = j, plot.j = 5,
#               spix = spix,
#               boundary.col = "gray",
#               boundary.lwd = 1,
#               incl.labels = spix,
#               col = "blue",
#               size = 3)
# rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)



########################################################################
########################## DEMO OF RESULTS #############################
########################################################################

nestSearch_step_old <- function(target, j1 = j2, j2, pix.j1 = 0,
                            demo.plot = FALSE) {

  if ( j1 > j2 ) stop("j1 must be less than or equal to j2")

  # nside at level j2
  nside.j2 <- 2^j2

  spix.j2 <- rcosmo:::pixelWindow(j1, j2, pix.j1)

  # Convert spix.j2 to Cartesian coordinates
  xyz.j2 <- rcosmo:::pix2coords_internal(nside = nside.j2,
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
    hp.j2 <- rcosmo:::pix2coords(nside = nside.j2,
                                 order = "nested",
                                 coords = "cartesian")
    rgl::plot3d(hp.j2, col="black", size = 2,
                type = 'p', pch = 2,  add = TRUE)

    # HEALPix points in window
    pixelWin <- rcosmo:::pixelWindow(j1, j2, pix.j1)
    pixelWin <- rcosmo:::pix2coords(nside = nside.j2,
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


nestSearch_old <- function(target, nside, index.only = FALSE,
                          j = 0:log2(nside), demo.plot = FALSE)
{
  j <- c(0,j)
  h <- list(pix = 0)
  for ( i in 2:length(j) )
  {
    h <- nestSearch_step_old( target, j2=j[i],
                              j1 = j[i-1], pix.j1 = h$pix,
                              demo.plot = demo.plot)
  }

  if (index.only==TRUE) {
    return(h$pix)
  }
  else {
    return(list( xyz=h$xyz, pix=h$pix ))
  }
}


# A demo of the new nestSearch versus the old nestSearch
demoNestSearch <- function()
{
  nside <- 64
  # The problematic pixel is located at (theta, phi) = (0.95, 0)
  theta <- 0.95
  # Conduct nested search for the nearest neighbour to
  # the problematic pixel, starting at base pixel resolution
  t <- c(sin(theta), 0, cos(theta))
  j <- log2(nside)
  spix1 <- nestSearch_old(t, nside = 2^j, demo.plot = TRUE)$pix
  spix2 <- rcosmo::nestSearch(t, nside = 2^j)$pix
  rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
  displayPixels(j = j, plot.j = 8,
                spix = spix1,
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = spix1,
                col = "blue",
                size = 3)
  displayPixels(j = j, plot.j = 8,
                spix = spix2,
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = spix1,
                col = "red",
                size = 3)
}
demoNestSearch()
