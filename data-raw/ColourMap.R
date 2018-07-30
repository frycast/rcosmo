library(R.matlab)


# INTERNAL DATA -----------------------------------------------------------

# CMBcolmap.mat was provided by Yu Guang to suit nside 2048
# mat <- readMat("data-raw/CMBcolmap.mat")
# colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])

# CMBcolmap.txt was intended for nside = 2048, found here:
# https://github.com/zonca/paperplots/blob/master/data/Planck_Parchment_RGB.txt
planck <- read.table("data-raw/CMBcolmap.txt")
colmap <- rgb(planck$V1, planck$V2, planck$V3, maxColorValue = 255)

map <- CMBReadFITS("../CMB_map_smica1024.fits", mmap = TRUE)
sky <- CMBDataFrame(map)
rx <- range(sky$I, na.rm = TRUE)
breaks1024 <- seq.int(rx[1L], rx[2L], length.out = 257L)

devtools::use_data(colmap, breaks1024, internal = TRUE, overwrite = TRUE)


# NON-INTERNAL DATA -------------------------------------------------------
# powspec <- read.table(paste0("http://pla.esac.esa.int/pla/aio/",
#                              "product-action?COSMOLOGY.FILE_ID=",
# "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt"),
#                       quote="\"",col.names = c("L","TT","TE","EE","BB","PP"))
# devtools::use_data(powspec, internal = FALSE)
