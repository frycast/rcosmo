y <- 2


f <- function() {
  cat("y = ", y, "\n")
  y <- 3
  g()
}


g <- function() {
  cat("y = ", y, "\n")
}


g()
f()


ff <- function(x) gg(x)
gg <- function(y) sys.status()
str(ff(1))
ff(1)
gg(1)



gg <- function(y) {
  ggg <- function() {
    cat("current frame is", sys.nframe(), "\n")
    cat("parents are", sys.parents(), "\n")
    print(sys.function(0)) # ggg
    print(sys.function(2)) # gg
  }
  if(y > 0) gg(y-1) else ggg()
}
gg(3)


#############


library(rcosmo)

worldcities <- read.csv("../worldcities.csv")

df <- data.frame(phi=worldcities$lng, theta=worldcities$lat,
                 I=rep(1,length(worldcities$lng)))
k <- 1000
df1 <- df[sample(nrow(df), k), ]
plot(df1$phi, df1$theta)
df1 <- coords(df1, new.coords = "cartesian")
df1[duplicated(df1), ]
df2 <- df1[!duplicated(df1), ]
df2[duplicated(df2), ]
hpdf <- HPDataFrame(df2, auto.spix = TRUE)
hpdf
str(hpdf)
plot(hpdf,size = 2)


library(rcosmo)
worldcities <- read.csv("../worldcities.csv")
uscities <- worldcities[worldcities$country == "United States",]
usdf <- data.frame(phi = pi/180*uscities$lng, theta = pi/2 - pi/180*uscities$lat,
                   I=rep(1,length(uscities$lng)))
k <- 1000
usdf <- usdf[sample(nrow(usdf), k), ]
plot(usdf$phi, usdf$theta)

usdf[duplicated(usdf), ]
usdf<- usdf[!duplicated(usdf), ]
usdf[duplicated(usdf), ]

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
