library(FITSio)
library(Rcpp)
sourceCpp("pix2ang.cpp")

cmbdat <- readFITS("CMB_map_smica1024.fits")

# HAVING A LOOK AT THE METADATA: -------------------------------------------------------

cmbdat$colNames                 # Clearly this map is polarised (because it has 5 planes).
cmbdat$colUnits                 # All units are Kcmb
cmbdat$hdr[72]; cmbdat$hdr[73]  # Pixtype:         healpix
cmbdat$hdr[76]; cmbdat$hdr[77]  # Coordsys:        galactic
cmbdat$hdr[78]; cmbdat$hdr[79]  # Ordering:        nested
cmbdat$hdr[80]; cmbdat$hdr[81]  # Nside:           1024
cmbdat$hdr[88]; cmbdat$hdr[89]  # Bad_data:        -1.63750E+30
cmbdat$hdr[98]; cmbdat$hdr[99]  # Further details: "http://wiki.cosmos.esa.int/planckpla2015"   


# CONSTANTS ----------------------------------------------------------------------------
Nside <- as.numeric(cmbdat$hdr[81])        # HEALPix resolution parameter.
Nrings <- 4*Nside - 1                      # Number of isolatitude rings.
Npix <- 12*Nside^2                         # Total number of pixels.
NpixPC <- sum(4*seq(1,Nside-1))            # Numer of pixels in each polar cap region
NringsPC <- Nside - 1                      # Number of rings in each polar cap region
pixArea <- 4*pi/numPixels                  # Area of each pixel (unit sphere).

# There are 2Nside - 1 rings in the equatorial zone.
# There are 4Nside - 1 rings in total.
# Hence there are Nside rings in each polar zone.
# Eq-zone rings have 4Nside pixels each.
# Num polar zone rings = Nside closest to eq-zone and decreases by 4 pix per ring.
# Hence:
pixPerRing <- c(4*seq(1,Nside),          # Stores number of pixels in ring i.
                rep(4*Nside,2*Nside - 1),  
                4*rev(seq(1,Nside)))     # Check: length(pixPerRing) == numRings
cPixPerRing <- cumsum(pixPerRing)        # Cumulative sum of pix per ring.

# Isolatitude-ring indices for North Polar, North Equatorial, South Equatorial and South Polar regions.
NPindex <- c(1,Nside - 1)
NEindex <- c(Nside,2*Nside)
SEindex <- c(2*Nside+1,3*Nside)
SPindex <- c(3*Nside+1,4*Nside - 1)

# Pixel indices for regional boundaries with RING ordering scheme
bpiRingNP <- c(0,cPixPerRing[NPindex[2]]-1)
bpiRingNE <- c(cPixPerRing[NPindex[2]],cPixPerRing[NEindex[2]] - 1) 
bpiRingSE <- c(cPixPerRing[NEindex[2]],cPixPerRing[SEindex[2]] - 1)
bpiRingSP <- c(cPixPerRing[SEindex[2]],cPixPerRing[SPindex[2]] - 1)

##CHECK: e.g. for Nside = 3:
#        bpiRingNP  # = 0 3
#        bpiRingNE  # = 12 59
#        bpiRingSE  # = 60 95
#        bpiRingSP  # = 96 107


# CREATING A BASIC CMB DATA FRAME -----------------------------------------

CMB <- data.frame(I = cmbdat$col[[1]], 
                  Q = cmbdat$col[[2]],
                  U = cmbdat$col[[3]],
                  TMASK = cmbdat$col[[4]],
                  PMASK = cmbdat$col[[5]])

str(CMB)


# CREATING INTENSITY DATA FRAME WITH SPHERICAL COORDS ---------------------

angSphereC <- pix2angC(Nside, FALSE)

CMBI <- data.frame(long = angSphereC[,2],          #long
                   lat = pi/2 - angSphereC[,1],    #lat
                   I = cmbdat$col[[1]])


# CREATING INTENSITY DATA FRAME WITH CARTESIAN COORDS ---------------------

library(sphereplot)
m <- matrix(c(CMBI$long, CMBI$lat, rep(1,Npix)), nrow = Npix)
m2 <- matrix(c(CMBI$long, CMBI$lat, CMBI$I*1e5), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
mx2 <- sph2car(m2, deg = FALSE)

# cool plot
plot3d(mx2, col = rainbow(1000), cex = 30, pch = 3, add = TRUE)




# NESTED ORDERING PIX2ANG -------------------------------------------------

#helper function converts int to bit vector.
intToBitVect <- function(x){
  tmp <- rev(as.numeric(intToBits(x)))
  id <- seq_len(match(1,tmp,length(tmp))-1)
  tmp[-id]
}
#helper function converts bit to int
bitVectToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

f <- 1  # f in [0,11] is the base resolution pixel index.
p <- 11 # p in [0,Npix-1] is the pixel index/order.

F1 <- floor(f/4) + 2                     # south vertex of f in latitude
F2 <- 2(f %% 4) - floor(f/4) %% 2 + 1    # south vertex of f in longitude 

pp <- p %% Nside^2                       # nested pixel index within f

# The following binary numbers are in reverse.
ppBin <- rev(intToBitVect(pp))
xBin <- intToBitVect(pp)[seq(0,length(ppBin),2)] 
yBin <- intToBitVect(pp)[seq(1,length(ppBin),2)]
# Convert back to int
x <- bitVectToInt(rev(xBin))
y <- bitVectToInt(rev(yBin))
# Vertical and horizontal coords
v <- x + y
h <- x - y

# Ring index
i <- F1*Nside - v - 1
j <- (F2*Nside + h + s)/2




