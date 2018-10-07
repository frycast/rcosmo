
# Map01 DATA --------------------------------------------------------------

# QUESTION: Why are there different elevations in the same exact position?
world <- read.csv("../Map01_000Ma.csv")
world$elev[world$lat == 90] # Here lon is irrelevant because lat == 90


colnames(world)[1] <- "lon"
plot(world$lat, world$lon)
str(world)
summary(world)
usdf <- data.frame(theta  = pi / 2 - pi / 180 * world$lat,
                   phi    = pi / 180 * world$lon,
                   I      = world$elev)
usdf2 <- data.frame(theta = lat2theta( deg2rad(world$lat) ),
                    phi   = lon2phi( deg2rad(world$lon) ),
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

dup(usdf)
dup(usdf2)


# Now they are both equal
all.equal(usdf, usdf2)

# Coords will produce phi to [0, 2*phi]
rcosmo::coords(usdf2) <- "cartesian"
rcosmo::coords(usdf2) <- "spherical"

# Still both equal
all.equal(usdf, usdf2)

plot(usdf2$theta, usdf2$phi)
plot(jitter(usdf2$theta), usdf2$phi)
plot(usdf2$theta, jitter(usdf2$phi))

# Difference is only floating point errors in coord conversion
all.equal(usdf, usdf2)
summary(usdf)
summary(usdf2)
sum(usdf$phi - usdf2$phi)
all.equal(round(usdf, 10), round(usdf2, 10))


summary(round(usdf ,1))
summary(round(usdf2,1))








# ARTIFICIAL DATA ---------------------------------------------------------

world <- expand.grid(lon = seq(-180, 160, by = 20),
                    lat = seq(-90, 80, by = 10))
world$elev <- rep(1, nrow(world))

  colnames(world)[1] <- "lon"
  plot(world$lat, world$lon)
  str(world)
  summary(world)
  usdf <- data.frame(theta  = pi / 2 - pi / 180 * world$lat,
                     phi    = pi / 180 * world$lon,
                     I      = world$elev)
  usdf2 <- data.frame(theta = lat2theta( deg2rad(world$lat) ),
                      phi   = lon2phi( deg2rad(world$lon) ),
                      I     = world$elev)

# The difference is because usdf has phi in [-pi,pi]
all.equal(usdf, usdf2)
summary(usdf)
summary(usdf2)

rcosmo::separatingNside(usdf)
rcosmo::separatingNside(usdf2)

# Jitter plots to check for overlap
plot(usdf$theta, usdf$phi)
plot(usdf2$theta, usdf2$phi)
plot(jitter(usdf2$theta), usdf2$phi)
plot(jitter(usdf$theta), usdf$phi)
plot(usdf2$theta, jitter(usdf2$ph))
plot(usdf$theta, jitter(usdf$phi))

# Neither has duplicates
nrow(dup(usdf))
nrow(dup(usdf2))


# Coords will produce phi in [0, 2*phi]
rcosmo::coords(usdf) <- "cartesian"
rcosmo::coords(usdf) <- "spherical"

dup(usdf)
dup(usdf2)

plot(usdf$theta, usdf$phi)
plot(jitter(usdf$theta), usdf$phi)
plot(usdf$theta, jitter(usdf$phi))

# Now they are both equal
all.equal(usdf, usdf2)

# Coords will produce phi to [0, 2*phi]
rcosmo::coords(usdf2) <- "cartesian"
rcosmo::coords(usdf2) <- "spherical"

# Still both equal
all.equal(usdf, usdf2)

plot(usdf2$theta, usdf2$phi)
plot(jitter(usdf2$theta), usdf2$phi)
plot(usdf2$theta, jitter(usdf2$phi))




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

lon2phi <- function(lon) {
  phi <- lon
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

  return(phi)
}

lat2theta <- function(lat) {
  theta <- pi/2 - lat

  negs <- theta < 0
  theta[negs] <- -theta[negs]

  high <- theta > pi
  theta[high] <- 2*pi - theta[high]

  return(theta)
}





# Duplicates (4 digit tolerance)
# dup(usdf)
all.equal(ndup(usdf), ndup(usdf2))


# Remove duplicates (4 digit tolerance)
usdf <- ndup(usdf)

plot(usdf$theta, usdf$phi)








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
