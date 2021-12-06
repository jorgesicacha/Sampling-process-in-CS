source("functionssimu.R")
source("Thinstage1.R")
source("Thinstage2.R")
source("Thinstage3.R")
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
    
    if(!is.null(input$ecological$fixed.effect$scale)){eco_scale=input$ecological$fixed.effect$scale} ##Test
    else{eco_scale=1} ##Test
    if(!is.null(input$sampling$fixed.effect$scale)){sampling_scale=input$sampling$fixed.effect$scale}
    else{sampling_scale=1}
    if(!is.null(input$detection$fixed.effect$scale)){detection_scale=input$detection$fixed.effect$scale}
    else{detection_scale=1}
    
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
    Scalelambda_test_rast <- list()
    for(i in 1:nspecies){
      Lam <- attr(Eco_PP[[i]], 'Lambda')
      Eco_GRF  <- log(Lam$v)
      Lambda_test <- eco_scale*Lam$v ##Test 
      gridlocs <- expand.grid(x0,y0) ##  Matching resolutions between gridlocs and covariates
      df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                                                          scalelambda =c(anti_t(rotate(rotate(Lambda_test))))
      )) ##Test
      r <- raster(df.sp)
      #r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
      xres <-  eco_covs.im[[1]]$xstep;yres <- eco_covs.im[[1]]$ystep## Raster resolution
      r1<-disaggregate(r, fact=res(r)/c(xres,yres))
      
      
      w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
      w1.rastaa <- crop(w1.rast,aa)
      scalelambda.rast <- rasterize(df.sp@coords,r1,df.sp$scalelambda, fun=mean,na.rm=T) ##Test
      scalelambda.rastaa <- crop(scalelambda.rast,aa) ##Test
      
      species_rast[[i]] <- mask(w1.rastaa,aa)
      Scalelambda_test_rast[[i]] <- mask(scalelambda.rastaa,aa) ##Test
      
    }
    
    #environment_list <- as.list(environment())
    #### First thinning stage ##
    #print(length(environment_list))
    
    firststage <- firstthinning(input,scale=sampling_scale)
    
    ## Second thinning stage ##
    environment_list <- as.list(environment())
    secondstage <- secondthinning(input,environment_list,scale=detection_scale)
    
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
  
  ## Generating lambda_obs ##
  lambda_obs_raster <- list() ##Test
  for(i in 1:nspecies){ ##Test
    lambda_obs_raster[[i]] <- Scalelambda_test_rast[[i]]*(firststage$retainprobraster)*(secondstage$detectionprobraster[[i]])
  }
  
  
  return(list(trueecological=Eco_PP,firststage=firststage,secondstage=secondstage,thirdstage=thirdstage,
              species_raster=species_rast,lambda_obs_raster=lambda_obs_raster))
  
}
