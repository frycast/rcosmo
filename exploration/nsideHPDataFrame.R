

f <- function(...)
{
  return(data.frame(list(...)))
}

f(data.frame(x = c(1,2,3), y = c(1,2,3)))




r <- apply(sky.s[, c("x","y","z")], MARGIN = 1, nestSearch, nside = 32, index.only = TRUE)



nestSearch(target = sky[,c("x","y","z")], nside = 1, index.only = TRUE)
