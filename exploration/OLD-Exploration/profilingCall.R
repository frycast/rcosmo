library(lineprof)
source("exploration/profiling.R")
l <- lineprof(HPDataFrame(sky.s, auto.spix = TRUE))
shine(l)




# TRYING TO VECTORISE nestSearch_step -------------------------------

m <- matrix(1:15, ncol = 3, byrow = TRUE)
tav <- c(1,1,1)
ta <- matrix(c(1,2,3,4), ncol = 3, nrow = 4)
tad <- data.frame(ta)

targets <- matrix(c(0,0,1,0,1,0,1,0,0,0.6,0.8,0),
                  byrow = TRUE, nrow = 4)
targets
h <- nestSearch_step(targets, j2 = 2)
h
class(h)

neighbours(c(1,2,3), 2)


nestSearch(c(0,0,1), nside = 4, index.only = TRUE)
nestSearch(c(0,1,0), nside = 4, index.only = TRUE)
nestSearch(targets, nside = 4, index.only = TRUE)

# TESTS AND CHECKS ---------------------------------------

library(microbenchmark)
remove <- 1:30
microbenchmark(
  remove[seq_along(remove) %% 2 > 0],
  remove[c(TRUE, FALSE)],
  remove[seq(1,length(remove),2)])


m <- matrix(1:1000, ncol = 10)
d <- 1:10

microbenchmark(
  m %*% diag(d),
  sweep(m, 2, d, FUN = "*"))

microbenchmark(
  sweep(m,2, d, FUN = function(x,d){(d-x)^2}),
  matrix((as.numeric(t(m)) - d)^2, nrow = nrow(m), byrow = TRUE),
  matrix(as.numeric(t(m))*d, nrow = nrow(m), byrow = TRUE),
  m %*% diag(d))

m1 <- rowSums(matrix(as.numeric(t(m))*d, nrow = nrow(m), byrow = TRUE))
m2 <- .rowSums(m %*% diag(d), m = nrow(m), n = ncol(m))
m3 <- as.numeric(m %*% d)
all.equal(m1,m2)
all.equal(m2,m3)

microbenchmark(
  rowSums(m),
  .rowSums(m, m = nrow(m), n = ncol(m)))


neighbours(1, 0)
demoNeighbours <- function(p,j) {
  neighbours(p, j)
  displayPixels(boundary.j = j, j = j, plot.j = j + 3,
                spix = neighbours(p, j),
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = neighbours(p, j),
                col = "blue",
                size = 3)
  rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
}

demoNeighbours(1,1)


point <- c(0.6,0.8,0)
j <- 2
cpoint <- nestSearch(point, nside = 2^j)

## Plot the closest pixel center in blue and the point (0.6,0.8,0) in red
displayPixels(j, j, plot.j=j, spix=c(cpoint$pix), size=5, incl.labels =FALSE)
rgl::plot3d(point[1], point[2], point[3], col="red", size = 5, add = TRUE)




points <- matrix(c(0,0,1,0,1,0,1,0,0,0.6,0.8,0),
                 byrow = TRUE, nrow = 4)
points
j <- 2
cpoints <- nestSearch(points, nside = 2^j)

## Plot the closest pixel center in blue and the point (0.6,0.8,0) in red
displayPixels(j, j, plot.j=j, spix=c(cpoints$pix), size=5, incl.labels =FALSE)
rgl::plot3d(points[,1], points[,2], points[,3], col="red", size = 5, add = TRUE)


