
# Map01 DATA --------------------------------------------------------------

# QUESTION: Why are there different elevations in the same exact position?
world <- read.csv("../Map01_000Ma.csv")
world$elev[world$lat == 90] # Here lon is irrelevant because lat == 90


colnames(world)[1] <- "lon"
str(world)
summary(world)

usdf <- data.frame(theta  = pi / 2 - pi / 180 * world$lat,
                   phi    = pi / 180 * world$lon,
                   I      = world$elev)

world$lat <- deg2rad(world$lat)
world$lon <- deg2rad(world$lon)
world <- geo2sph(world)
usdf2 <- data.frame(theta = world$theta,
                    phi   = world$phi,
                    I     = world$elev)


################## FIXED
## The problem was with theta = 0 and theta = pi
a <- usdf[usdf$theta == pi, ]
b  <- usdf2[usdf$theta == pi, ]

rcosmo:::sph2car(a)
rcosmo:::sph2car(b)

rcosmo:::car2sph(rcosmo:::sph2car(a))
rcosmo:::car2sph(rcosmo:::sph2car(b))

a2 <- usdf[usdf$theta == 0, ]
b2  <- usdf2[usdf$theta == 0, ]

rcosmo:::sph2car(a2)
rcosmo:::sph2car(b2)

rcosmo:::car2sph(rcosmo:::sph2car(a2))
rcosmo:::car2sph(rcosmo:::sph2car(b2))
########################



# The difference here is because usdf has phi in [-pi,pi]
all.equal(usdf, usdf2)
summary(usdf)
summary(usdf2)

rcosmo::separatingNside(usdf)
rcosmo::separatingNside(usdf2)

# usdf doesn't have duplicates (but it should)
# the problem is in the domain of phi
nrow(dup(usdf))
nrow(dup(usdf2))

# Coords will produce phi in [0, 2*phi]
rcosmo::coords(usdf) <- "cartesian"
rcosmo::coords(usdf) <- "spherical"

nrow(dup(usdf))
nrow(dup(usdf2))

# Now they are both equal
all.equal(usdf, usdf2)

# Coords will produce phi to [0, 2*phi]
rcosmo::coords(usdf2) <- "cartesian"
rcosmo::coords(usdf2) <- "spherical"

# Still both equal
all.equal(usdf, usdf2)

usdf <- ndup(usdf)
dup(usdf)

separatingNside(usdf) #Large
hpdf <- HPDataFrame(usdf, auto.spix = TRUE)

map2color <- function(x, pal, limits = range(x)){
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1),
                   all.inside=TRUE)]
}

plot(hpdf, col = map2color(hpdf$I, terrain.colors(10000)), size = 5)


hpdf2 <- HPDataFrame(usdf, nside = 64, auto.spix = TRUE,
                    delete.duplicates = TRUE)

plot(hpdf2, col = map2color(hpdf2$I, terrain.colors(10000)), size = 10)



# ARTIFICIAL DATA ---------------------------------------------------------

# # Has no duplicates
# world <- expand.grid(lon = seq(-180, 160, by = 20),
#                     lat = seq(-90, 80, by = 10))

# Has duplicates
world <- expand.grid(lon = seq(-180, 180, by = 20),
                     lat = seq(-90, 90, by = 10))

world$elev <- rep(1, nrow(world))

  colnames(world)[1] <- "lon"
  plot(world$lat, world$lon)
  str(world)
  summary(world)
  usdf <- data.frame(theta  = pi / 2 - pi / 180 * world$lat,
                     phi    = pi / 180 * world$lon,
                     I      = world$elev)


  # geo2sph will standardise phi when theta in [0,pi]
  geo2sph(lat = pi/2, lon = 2)
  geo2sph(lat = pi/2, lon = 1)

  world$lat <- deg2rad(world$lat)
  world$lon <- deg2rad(world$lon)
  world <- geo2sph(world)
  usdf2 <- data.frame(theta = world$theta,
                      phi   = world$phi,
                      I     = world$elev)

# The difference is because usdf has phi in [-pi,pi]
all.equal(usdf, usdf2)
summary(usdf)
summary(usdf2)

rcosmo::separatingNside(usdf)
rcosmo::separatingNside(usdf2)


# Both should have duplicates
nrow(dup(usdf))
nrow(dup(usdf2))


# Jitter plots to check for overlap
X11()

plot(usdf2$theta, usdf2$phi)
plot(usdf2$theta, jitter(usdf2$ph))
plot(jitter(usdf2$theta), usdf2$phi)

plot(usdf$theta, usdf$phi)
plot(jitter(usdf$theta), usdf$phi)
plot(usdf$theta, jitter(usdf$phi))


# Coords will produce phi in [0, 2*phi]
rcosmo::coords(usdf) <- "cartesian"
rcosmo::coords(usdf) <- "spherical"

# Should share same number of duplicates
nrow(dup(usdf))
nrow(dup(usdf2))

# The positions agree now
X11()
plot(usdf2$theta, usdf2$phi)
plot(usdf2$theta, jitter(usdf2$ph))
plot(jitter(usdf2$theta), usdf2$phi)

plot(usdf$theta, usdf$phi)
plot(jitter(usdf$theta), usdf$phi)
plot(usdf$theta, jitter(usdf$phi))

# Now they are both equal
all.equal(usdf, usdf2)


rcosmo::coords(usdf2) <- "cartesian"
rcosmo::coords(usdf2) <- "spherical"

# Still both equal
all.equal(usdf, usdf2)





# HELPER FUNCTIONS --------------------------------------------------------

dup <- function(usdf, digits = 4) {
  usdf[duplicated(round(usdf[,c(1,2)], digits = digits)), ]
}

ndup <- function(usdf, digits = 4){
  usdf[!duplicated(round(usdf[,c(1,2)], digits = digits)), ]
}


deg2rad <- function(deg) {
  pi/180*deg
}


#' geo2sph
#'
#' Convert latitude (lat) and longitude (lon) to spherical
#' coordinates (theta, phi) with theta in [0,pi] and
#' phi in [0,2*pi).
#' All values are assumed to be in radians.
#'
#' @param ... A data.frame with columns lat and lon,
#' or named vectors of lat and lon.
#'
#'
#' @export
geo2sph <- function(...) {

  df <- data.frame(...)

  if ( all(c("lat","lon") %in% names(df)) ) {

    lat <- df$lat
    lon <- df$lon

    theta <- pi/2 - lat
    phi <- lon

    negs <- theta < -1e-13
    theta[negs] <- -theta[negs]

    high <- theta > pi+1e-13
    theta[high] <- 2*pi - theta[high]

    # These are now in e.g., [0,1e-13]
    zeros <- theta <= 0
    pies <- theta >= pi

    theta[zeros] <- 0
    theta[pies] <- pi
    phi[zeros] <- 0
    phi[pies] <- 0

    negs <- phi < 0
    while ( any(negs) ) {
      phi[negs] <- phi[negs] + 2*pi
      negs <- phi < 0
    }

    high <- phi > 2*pi
    while ( any(high) ) {
      phi[high] <- phi[high] - 2*pi
      high <- phi > 2*pi
    }

    df$lat <- theta
    df$lon <- phi
    names(df)[names(df) == "lat"] <- "theta"
    names(df)[names(df) == "lon"] <- "phi"

  }

  return(df)
}







#############


library(rcosmo)
worldcities <- read.csv("../worldcities.csv")
#uscities <- worldcities[worldcities$country == "United States",]
uscities <- worldcities
usdf <- data.frame(phi = pi/180*uscities$lng, theta = pi/2 - pi/180*uscities$lat,
                   I=rep(1,length(uscities$lng)))

# k <- 1000
# usdf <- usdf[sample(nrow(usdf), k), ]
plot(usdf$phi, usdf$theta)

usdf[duplicated(usdf), ]
usdf<- usdf[!duplicated(usdf), ]
usdf[duplicated(usdf), ]

ushp <- HPDataFrame(usdf, auto.spix = TRUE)
plot(ushp)

usdf <- coords(usdf, new.coords = "cartesian")
ushp <- HPDataFrame(usdf, auto.spix = TRUE)
plot(ushp, size = 3)



## SOLVED
library(rcosmo)
hp1 <- HPDataFrame(I=rnorm(5), nside = 1, spix = c(1,1,2,2,3))
hp1 <- coords(hp1, new.coords = "cartesian")
separatingNside(hp1) #Inf
hpdf <- HPDataFrame(hp1, auto.spix = TRUE) #Error

df <- data.frame(theta = c(0,0.0001), phi = c(pi/2, pi/2))
coords(df) <- "cartesian"
separatingNside(df) #8192 (large)
HPDataFrame(df, auto.spix = TRUE) # Error
