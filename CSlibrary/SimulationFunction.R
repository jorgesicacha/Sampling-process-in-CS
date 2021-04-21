## A function for simulating CS data ##
###########
## INPUT ##
###########

nspecies = 4

## For the species we're interested in (Species 1)
input <- data.frame(
  sigma2 = c(0.2, 1.2, 2, 0.1),
  range = c(1.2, 2.5, 3.2, 0.22),
  beta0 = c(0.8, 2.5, -1.5, 1.2),
  betacov = c(1.5, -0.12, 2,-0.4),
  alpha0 = c(2,-0.3, 5, 1.2),
  alphacov = c(-2, -0.5, -2.5, 2)
)

thin_input <- list(
  mean_thin = 1.3,
  betacov_thin = -1.5,
  sigma2x_thin = 0.2,
  range_thin = 2.5
)

seed= 1036610602

## Covariates ##

x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = 100)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = 100)
gridlocs <- expand.grid(x0,y0)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
covariate_thin.im <- im(gridcov_thin, x0, y0)

gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_detect.im <- im(gridcov_det, x0, y0)




csdata <- function(nspecies,input,thin_input,cov,domain=NULL,seed){
  
  source("functionssimu.R")
  set.seed(seed)
  RFoptions(seed=seed)
  
  
  
  ## Where the simulation is performed ##
  if(is.null(domain)){
  coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
  aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
  win <- as.owin(aa)
  }
  
  ## Mesh for models ##
  ## Maximum distance in the extent of study area ##
  ext <- extent(aa)
  pts <- SpatialPoints(coords = matrix(ext,ncol=2),proj4string = crs(aa))
  max.dist <- spDists(pts)[1,2]
  mesh <- inla.mesh.2d(loc.domain = coordsmat,max.edge = c(0.02*max.dist,0.10*max.dist),offset = c(0.3, 1),cutoff = 0.2*max.dist)
  
  if(is.list(cov)){
    
    classes <- unique(lapply(cov, class))
    if(length(classes)==1){
      
      if(classes=="im"){
        covs.im <- cov
        covs.raster <- lapply(cov,im2rast)
        covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
      }
      else{if(classes=="RasterLayer"){
        covs.im <- lapply(cov,as.im)
        covs.raster <- cov
        covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
      }
        else{if(classes=="SpatialPixelsDataFrame"){
          covs.im <- lapply(cov,as.im)
          covs.raster <- lapply(cov, sppixels2raster)
          covs.sppixels <- cov
          }
          else{
            stop("Covariates must be of 'im', 'RasterLayer' or 'SpatialPixelsDataFrame'")
          }
          }}
      
    }
    else{
      stop("All the covariate must be in the same format")
    }

    

  }
  else{
    stop("Covariates input should be a list")
  }
}
