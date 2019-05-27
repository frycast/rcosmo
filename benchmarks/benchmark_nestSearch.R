library(microbenchmark)

target <- CMBDataFrame(nside = 2, coords = "cartesian")
target$I <- NULL

microbenchmark(nestSearch(target, nside = 4))

### 27/05/2019 on Desktop Doony
# Unit: milliseconds
# min      lq       mean    median   uq      max neval
# 13.60496 13.88428 14.7277 14.08012 14.4223 23.63274   100
