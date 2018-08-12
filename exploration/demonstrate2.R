




##### Examples

## Specify locations as vectors and use auto.spix
hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
                   nside = 1, auto.spix = TRUE)
class(hp)
pix(hp)
plot(hp, size = 5, hp.boundaries = 1)
plotHPBoundaries(nside = 1, ordering = "nested", col = "gray")

## Specify locations in data.frame and specify spix manually
d <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
hp2 <- HPDataFrame(d, nside = 1, spix = c(1,2,1))
plot(hp2, size = 5, hp.boundaries = 1)

## Do not specify locations (get all pixels at nside)
hp3 <- HPDataFrame(I = rep(0,12), nside = 1)
plot(hp3, size = 5, hp.boundaries = 1)
