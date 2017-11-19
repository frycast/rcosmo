# Mollview Map for CMB Data
# 20/08/2017
# Ming Li
rm(list=ls())
library(mapproj)
library(R.matlab)
library(sphereplot)
library(FITSio)
source("pix2vec.R")  # from Yuguang's nest search folder
source("readFITScmb.R")
cmbdat <- readFITScmb("CMB_map_smica1024.fits")
nside<- 2048

Npix<-12*nside^2
pix<- 1:Npix
type="nest"
Cardi_xyz <- pix2vec(nside,type,pix)   #note that this differs from pix2angC
Cardi_xyz<-data.frame(x=Cardi_xyz[1,],y=Cardi_xyz[2,],z=Cardi_xyz[3,])
Sph_long_lat<-car2sph(Cardi_xyz)  # obtain the correct longitudes and latitudes 
CMBI<- data.frame(long=Sph_long_lat[,1], lat = Sph_long_lat[,2], I = cmbdat$col[[1]])

# use the mapproject function in mapproj package
b<-mapproject(Sph_long_lat[,1], Sph_long_lat[,2], projection="mollweide")
mat <- readMat("cmbmap.mat")
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])

#Plot
cols <- colmap[cut(CMBI$I,length(colmap))] # partition the colour map
## it is expensive!!!
plot(b$x,b$y,col=cols,type = "p")


#####further notes#####
## Once the CMBDataFrame is fixed, I can rewrite this scrpit as a function


####
b_dataframe<- data.frame(b$x,b$y)
sNside <- 256
sNpix <- 12*sNside^2
sample <- sample(seq(1,Npix),sNpix)
b_sample<- b_dataframe[sample,]
CMBI_sample<- CMBI[sample,]
cols_sample <- colmap[cut(CMBI_sample$I,length(colmap))]
plot(b_sample[,1],b_sample[,2],cols_sample,type = "p")