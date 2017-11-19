library(R.matlab)
mat <- readMat("cmbmap.mat") # cmbmap.mat was provided by Yu Guang.
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])
devtools::use_data(colmap, internal = TRUE)
