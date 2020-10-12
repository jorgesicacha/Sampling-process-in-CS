
library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(spatstat)
library(maptools)
library(gdata)
library(ggplot2)
library(gridExtra)
library(PCDSpline)
library(foreach)
library(doParallel)
library(viridis)
library(RandomFields)


rotate <- function(x) (apply(t(x), 2, rev))
anti_t <- function (m){
  p <- nrow(m)
  j <- matrix(ncol = p, nrow = p, data = 0)
  for (i in 1:p) {
    j[i, p - i + 1] <- 1
  }
  return(j %*% t(m) %*% j)
}
