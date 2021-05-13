im2rast <- function(im){
  gridlocs <- expand.grid(im$xcol,im$yrow)
  cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(im$v))))))
  r <- raster(cov1.sp)
  xstp <- im$xstep
  ystp <- im$ystep
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}
