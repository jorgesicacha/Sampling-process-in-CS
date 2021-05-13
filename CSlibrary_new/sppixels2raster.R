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

