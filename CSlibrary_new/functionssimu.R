
library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(maptools)
library(spatstat)
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

im2rast <- function(im){
  gridlocs <- expand.grid(im$xcol,im$yrow)
  cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(im$v))))))
  r <- raster(cov1.sp)
  xstp <- im$xstep
  ystp <- im$ystep
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}

sppixels2raster <- function(sppixels){
  cov1.sp <- SpatialPointsDataFrame(coords = sppixels@coords,data = sppixels@data,proj4string = crs(sppixels))
  r <- raster(cov1.sp)
  xs <- unique(sppixels@coords[,1])
  ys <- unique(sppixels@coords[,2])
  xstp <- xs[2] - xs[1]
  ystp <- ys[2] - ys[1]
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}


