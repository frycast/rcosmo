# library(R.matlab)
# # CMBcolmap.mat was provided by Yu Guang to suit nside 2048
# mat <- readMat("data-raw/CMBcolmap.mat")
# colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])

# CMBcolmap.txt was intended for nside = 2048, found here:
# https://github.com/zonca/paperplots/blob/master/data/Planck_Parchment_RGB.txt
planck <- read.table("data-raw/CMBcolmap.txt")
colmap <- rgb(planck$V1, planck$V2, planck$V3, maxColorValue = 255)
#devtools::use_data(colmap, internal = TRUE)

# Generate colours for nside = 1024 though we should be doing nside = 2048
library(rcosmo)
cmbdf <- CMBDataFrame("../CMB_map_smica1024.fits")
cols <- colmap[cut(cmbdf$I, length(colmap))]
CMBcols1024 <- factor(cols)
devtools::use_data(CMBcols1024, internal = TRUE, overwrite = TRUE)

