filename <- "../CMB_map_smica1024.fits"
filename <- "C:/Users/danie/Downloads/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits"
zz <- file(filename, "rb")

# num <- 36
# cols <- 80
# maxLines <- 5000
# maxHdrs <- maxLines/num
# header <- "dummy start line"
# image <- character(num)
# start <- seq(1, 2880, by = 80)


## In CMBReadFITS we run readFITSHeader twice.
## Each time, the following commented code (roughly) is run.
## Based on the counts here, we can just run
## readChar 3 times, though this probably won't
## generalise well to other FITS files, especially
## ones that aren't for CMB data.
# count <- 1
# for (i in 1:maxHdrs)
# {
#   ## ALL THE READING WORK IS HERE IN readChar:
#   inpString <- readChar(zz, 2880)
#   nchar(inpString) # should be 2880
#   for (j in 1:num) {
#     image[j] <- substr(inpString, start[j], start[j] +
#                          79)
#   }
#   idx <- grep("^ *END +", image, ignore.case = TRUE)
#   if (length(idx) > 0) break;
#   count <- count + 1
# }

readChar(zz, 2880)
readChar(zz, 2880)
readChar(zz, 2880)

swap <- "big" != .Platform$endian

# For 2048 the following works:
.Internal(readBin(zz, "double", 1, 4, TRUE, swap))
.Internal(readBin(zz, "integer", 1, 1, FALSE, swap))

# For 1024 the following works:
.Internal(readBin(zz, "double", 1, 4, TRUE, swap))
.Internal(readBin(zz, "double", 1, 4, TRUE, swap))
.Internal(readBin(zz, "double", 1, 4, TRUE, swap))
.Internal(readBin(zz, "integer", 1, 1, FALSE, swap))
.Internal(readBin(zz, "integer", 1, 1, FALSE, swap))

#-9.201023e-05

# for (j in 1:num) {
#   image[j] <- substr(inpString, start[j], start[j] +
#                        79)
# }
# idx <- grep("^ *END +", image, ignore.case = TRUE)
#
#
# if (length(idx) > 0) {
#   image <- image[-(idx:num)]
#   header <- c(header, image)
#   header <- header[-1]
#   #return(header)
# }

close(zz)
