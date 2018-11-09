library(rcosmo)
## With earth data.
## Download World Cities Database from
## https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.4.zip
## unpack the file worldcities.csv


#################################################################
#              NON-UNIQUE PIXEL HPDataFrame World
#################################################################

worldcities <- read.csv("worldcities.csv")

## Prepare a data frame with cities' coordinates
sph <- geo2sph(data.frame(lon = pi/180*worldcities$lng,
                          lat = pi/180*worldcities$lat))
df <- data.frame(phi = sph$phi,
                 theta = sph$theta,
                 I = rep(1,nrow(sph)))

## Create and plot the corresponding HPDataFrame with
## pixel indices that are not necessarily unique
## by choosing your desired resolution (nside)
hp <- HPDataFrame(df, auto.spix = TRUE, nside = 1024)
plot(hp, size = 3, col = "darkgreen", back.col = "white")
## Add some pixels to visualise the sphere
plot(CMBDataFrame(nside = 64), add = TRUE, col = "gray")


#################################################################
#              UNIQUE PIXEL HPDataFrame Australia
#################################################################

worldcities <- read.csv("worldcities.csv")
au.cities <- worldcities[worldcities$country == "Australia",]

# Prepare a data frame with cities' coordinates
sph <- geo2sph(data.frame(lon = pi/180*au.cities$lng,
                          lat = pi/180*au.cities$lat))
df <- data.frame(phi = sph$phi,
                 theta = sph$theta,
                 I = rep(1,nrow(sph)))

# Create and plot the corresponding HPDataFrame with unique pixels
hp <- HPDataFrame(df, auto.spix = TRUE)
nside(hp)
assumedUniquePix(hp)
plot(hp, size = 5, col = "darkgreen", back.col = "white")
## Add some pixels to visualise the sphere
plot(CMBDataFrame(nside = 64), add = TRUE, col = "gray")


#################################################################
#              BENCHMARK nest2ring versus ring2nest
#################################################################

df_n <-CMBDataFrame(nside = 16, ordering = "nested")
df_r <-CMBDataFrame(nside = 16, ordering = "ring")

microbenchmark::microbenchmark(
  nest2ring  = ordering(df_n, new.ordering = "ring"),
  ring2nest  = ordering(df_r, new.ordering = "nested")
)

