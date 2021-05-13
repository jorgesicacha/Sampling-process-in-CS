## A function for simulating CS data ##
###########
## INPUT ##
###########
#setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds")
#source("functionssimu.R")

#Number of species
nspecies = 4

## For the species we're interested in (Species 1)
# input <- data.frame(
#   sigma2 = c(0.2, 1.2, 2, 0.1),
#   range = c(1.2, 2.5, 3.2, 0.22),
#   beta0 = c(0.8, 2.5, -1.5, 1.2),
#   betacov = c(1.5, -0.12, 2,-0.4),
#   alpha0 = c(2,-0.3, 5, 1.2),
#   alphacov = c(-2, -0.5, -2.5, 2)
# )


#Input for the simulations
input <- list(
  ecological = list(
    fixed.effect=list(
      intercept = c(0.8, 0.5, 0.5, 1.2),
      betacov = c(1.5, -2.12, 4,-0.4)
    ),
    hyperparameters = list(
      sigma2 = c(2.2, 1.2, 2, 0.1),
      range = c(3.2, 2.5, 3.2, 0.22)
    )
  ),
  sampling = list(
    fixed.effect = list(
      intercept = c(1.3),
      betacov = c(-1.5)
    ),
    hyperparameters=list(
      sigma2 = c(0.2),
      range = c(2.5)
    )
  ),
  detection = list(
    
    fixed.effect = list(
      intercept=c(2,-0.3, 5, 1.2),
      betacov = c(-2, -0.5, -2.5, 2)
    )
  ),
  
  misclassification = list(
    
    class_prob <- matrix(c(0.9, 0.02, 0.04, 0.04,
                           0.05, 0.89, 0.04, 0.02,
                           0.1,0.1, 0.8, 0,
                           0, 0.05, 0.25, 0.7),
                         nrow=4, ncol=4, byrow = TRUE)
    
    
  )
)


seed= 1036610620

## Simulating the Covariates ##

#Grid for simulation
x0 <- seq(-1, 4, length = 100)
y0 <- seq(-1,4, length = 100)
gridlocs <- expand.grid(x0,y0)

# Covariates for true ecological state
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

#Covariate for the sampling process
gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
#gridcov_thin <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate_thin.im <- im(gridcov_thin, x0, y0)

#Covariate for the detection
gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_detect.im <- im(gridcov_det, x0, y0)

#Make a list for the covariates
cov <- list(covariate.im,covariate_thin.im,covariate_detect.im)

#Get an index for each covariate 
idxs <- list(eco=c(1),sampling=c(2),detection=c(3))


## Plotting configuration ##
plot <- list(ecological=TRUE,detection=FALSE,sampling=TRUE,all=TRUE,classification=FALSE)

par(mfrow=c(2,2))
## Citizen Science data generation function ##
csdata <- function(nspecies,input,cov,idxs,domain=NULL,seed,plot=list(all=TRUE),colmatrix=NULL){
  
  #Check if the number of species provided is a number
  if(class(nspecies)!="numeric") stop("Number of species must be a number")
  
  #Setting seed for the simulation
  set.seed(seed)
  RFoptions(seed=seed)
  
  # The domain where the simulations are done
  # The default is provided here
  if(is.null(domain)){
    coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
    aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
    win <- as.owin(aa) ##Need maptools
  }
  
  ## Mesh for models ##
  ## Maximum distance in the extent of study area ##
  ext <- extent(aa)
  pts <- SpatialPoints(coords = matrix(ext,ncol=2),proj4string = crs(aa))
  max.dist <- spDists(pts)[1,2]
  mesh <- inla.mesh.2d(loc.domain = coordsmat,max.edge = c(0.02*max.dist,0.10*max.dist),offset = c(0.3, 1),cutoff = 0.2*max.dist)
  
  ## Converting the covariates provided into a raster, image and pixels ##
  if(is.list(cov)){
    
    classes <- unique(sapply(cov, class))
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
    
    ## Indexes of the covariates##  
    eco_idxs <- idxs$eco
    sampling_idxs <- idxs$sampling
    detection_idxs <- idxs$detection
    
    eco_covs.im <- covs.im[eco_idxs]
    eco_covs.raster <- covs.raster[eco_idxs]
    eco_covs.sppixels <- covs.sppixels[eco_idxs]
    
    sampling_covs.im <- covs.im[sampling_idxs]
    sampling_covs.raster <- covs.raster[sampling_idxs]
    sampling_covs.sppixels <- covs.sppixels[sampling_idxs]
    
    detection_covs.im <- covs.im[detection_idxs]
    detection_covs.raster <- covs.raster[detection_idxs]
    detection_covs.sppixels <- covs.sppixels[detection_idxs]
    
    ## Generate ecological process ##
    
    #Ecological Process
    Eco_PP <- list()
    
    ## Bringing ecological covariates ##
    p_eco <- length(eco_idxs)
    eco_form <- c()
    for(j in 1:p_eco){
      eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+1,"]][i]*eco_covs.im[[",j,"]]")
    }
    
    eco_linpred <- paste(c("input$ecological$fixed.effect[[1]][i]",eco_form),collapse="+")
    
    for(i in 1:nspecies){
      x0 <- eco_covs.im[[1]]$xcol
      y0 <- eco_covs.im[[1]]$yrow
      Eco_PP[[i]] <- rLGCP(model="matern",mu=eval(parse(text=eco_linpred)),
                           var=input$ecological$hyperparameters[[1]][i],scale=input$ecological$hyperparameters[[2]][i]/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))
      
    }
    
    # Storing the raster of the true intensity for each species
    species_rast <- list()
    for(i in 1:nspecies){
      Lam <- attr(Eco_PP[[i]], 'Lambda')
      
      Eco_GRF  <- log(Lam$v)
      gridlocs <- expand.grid(x0,y0) ##  Matching resolutions between gridlocs and covariates
      df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v)))))))
      r <- raster(df.sp)
      r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
      w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
      w1.rastaa <- crop(w1.rast,aa)
      species_rast[[i]] <- mask(w1.rastaa,aa)
    }
    
    #environment_list <- as.list(environment())
    #### First thinning stage ##
    #print(length(environment_list))
    firststage <- firstthinning(input)
    
    ## Second thinning stage ##
    environment_list <- as.list(environment())
    secondstage <- secondthinning(input,environment_list)
    
    ## Second thinning stage ##
    environment_list <- as.list(environment())
    thirdstage <- thirdthinning(input,environment_list)
    
    if(plot$all==TRUE){
      ##True ecological ##
      lapply(1:nspecies,function(x){plot(species_rast[[x]],axes=F,box=F,main=paste0("True Ecological State for species",x))
        points(Eco_PP[[x]],pch=19,cex=0.5) })
      
      ## Detection##
      
      lapply(1:nspecies,function(x){plot(secondstage$detectionprobraster[[x]],axes=F,box=F,main=paste0("Detection probability for species ",x))
        points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.5)
        points(secondstage$Eco_PPFinal_detect[[x]],pch=19,cex=0.3,col="red") })
        
      ## Sampling ##
      lapply(1:nspecies,function(x){plot(firststage$retainprobraster,axes=F,box=F,main=paste0("Retaining probability for species ",x))
        points(Eco_PP[[x]],pch=19,cex=0.5)
        points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.3,col="red")
      })
      
      ## Classification ##
      if(is.null(colmatrix)){
        colmatrix <- matrix(NA,nrow=nspecies,ncol=2)
        colmatrix[,1] <- sample(colors(),size = nspecies,replace = FALSE)
        colmatrix[,2] <- as.numeric(sample(0:25,size = nspecies,replace = FALSE))
        colmatrix <- data.frame(colmatrix)
        names(colmatrix) <- c("color","pch")
        colmatrix$pch <- as.numeric(colmatrix$pch)
      }
      
      for(i in 1:nspecies){
        plot(aa,axes=F,main=paste0("CS reports for species ",i))
        points(thirdstage$classifications[which(thirdstage$classifications$true_species==i),],pch=colmatrix[i,2],cex=1,col=colmatrix[i,1])
        indexes0 <- 1:nspecies
        indexes <- indexes0[-indexes0[which(indexes0==i)]]
        for(j in indexes){
          try(points(thirdstage$classifications[which(thirdstage$classifications$true_species==i & thirdstage$classifications$error==j),],pch=colmatrix[j,2],cex=0.6,col=colmatrix[j,1]))
        }
        legendtext <- paste0("Species ",1:nspecies)
        legend("right", legend=legendtext,col=colmatrix$color,pch=colmatrix$pch)
      }
      
    }
    else{
      plotnew <- within(plot,rm(all))
      which.plot <- names(plotnew[sapply(plotnew,isTRUE)])
      
      if("ecological"%in%which.plot){
        lapply(1:nspecies,function(x){plot(species_rast[[x]],axes=F,box=F,main=paste0("True Ecological State for species",x))
          points(Eco_PP[[x]],pch=19,cex=0.5)
        })
      }
      
      if("detection"%in%which.plot){
        lapply(1:nspecies,function(x){plot(secondstage$detectionprobraster[[x]],axes=F,box=F,main=paste0("Detection probability for species ",x))
          points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.5)
          points(secondstage$Eco_PPFinal_detect[[x]],pch=19,cex=0.3,col="red")
          
        })
        
      }
      
      if("sampling"%in%which.plot){
        lapply(1:nspecies,function(x){plot(firststage$retainprobraster,axes=F,box=F,main=paste0("Retaining probability for species ",x))
          points(Eco_PP[[x]],pch=19,cex=0.5)
          points(firststage$Eco_PPFinal[[x]],pch=19,cex=0.3,col="red")
        })
      }
      
      if("classification"%in%which.plot){
        if(is.null(colmatrix)){
        colmatrix <- matrix(NA,nrow=nspecies,ncol=2)
        colmatrix[,1] <- sample(colors(),size = nspecies,replace = FALSE)
        colmatrix[,2] <- as.numeric(sample(0:25,size = nspecies,replace = FALSE))
        colmatrix <- data.frame(colmatrix)
        names(colmatrix) <- c("color","pch")
        colmatrix$pch <- as.numeric(colmatrix$pch)
        }
        
        for(i in 1:nspecies){
          plot(aa,axes=F,main=paste0("CS reports for species ",i))
          points(thirdstage$classifications[which(thirdstage$classifications$true_species==i),],pch=colmatrix[i,2],cex=1,col=colmatrix[i,1])
          indexes0 <- 1:nspecies
          indexes <- indexes0[-indexes0[which(indexes0==i)]]
          for(j in indexes){
          try(points(thirdstage$classifications[which(thirdstage$classifications$true_species==i & thirdstage$classifications$error==j),],pch=colmatrix[j,2],cex=0.6,col=colmatrix[j,1]))
          }
          legendtext <- paste0("Species ",1:nspecies)
          legend("right", legend=legendtext,col=colmatrix$color,pch=colmatrix$pch)
          }
        
      }
      
    }
    
    
  }
  else{
    stop("Covariates input should be a list")
  }
}
#par(mfrow=c(2,2))
csdata(nspecies=nspecies,input=input,cov=cov,idxs=idxs,seed=seed)

#Out put an object and plot late on.
#We can do much more with the object returned.

# The negative should give the probability less than 0.

# Identifiability of the model
#Which parameter values would make the model identifiable and those that won't .

