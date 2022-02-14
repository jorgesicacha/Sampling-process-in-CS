### Example of data generation and model fitting with Omega assumed deterministic ###

## Simulating the Covariates ##


library(spatstat)
#Grid for simulation
x0 <- seq(-1.3, 4.3, length = 100)
y0 <- seq(-1.3,4.3, length = 100)
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



## Generating the data ##

source("DataGeneration.R")

nspecies <- 4#nspecies: Number of species we want to simulate 
#input: List with the parameters of the model that generates CS data
input <-{list(
  ecological = list(
    fixed.effect=list(
      intercept = c(0.8, 2.5, -1.5, 1.2),
      betacov = c(1.5, -0.12, 2,-0.4)
    ),
    hyperparameters = list(
      sigma2 = c(0.2, 1.2, 2, 0.1),
      range = c(1.2, 2.5, 3.2, 0.22)
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
)}
# cov: a list with the covariates needed by both the data generating process and the model we want to fit.
#idxs: A very simple list which tells which covariates belong to each stage of the data generating process
idxs <- list(eco=c(1),sampling=c(2),detection=c(3))
#seed: Random seed for replicating the datasets generated
seed <- 1036610602
#plot: Do we want to plot the results? Which results? ISSUE
plot = list(all=FALSE,none=TRUE)
#colmatrix: Color scale for the plots.If plot !=NULL


simulateddata <- csdata(nspecies=nspecies,input=input,cov=cov,idxs=idxs,seed=seed,plot=plot)


## Covariates need to be in SpatialPixels format ##
#Covariates for true intensity
cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate.im$v))))))
r <- raster(cov1.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")

#Covariates for first thinning
cov2.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_thin.im$v))))))
r <- raster(cov2.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov2.rast <- rasterize(cov2.sp@coords,r1,cov2.sp$cov, fun=mean,na.rm=T)
cov2.spix <- as(cov2.rast,"SpatialPixelsDataFrame")

#Covariate for second thinning
cov3.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_detect.im$v))))))
r <- raster(cov3.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov3.rast <- rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
cov3.spix <- as(cov3.rast,"SpatialPixelsDataFrame")

## Extra information on species detection ##

### Sampling the detections from a survey
rndpts_x0 <- runif(50, 0,3)
rndpts_y0 <- runif(50, 0,3)
rndpts <- data.frame(rndpts_x0, rndpts_y0)

det_prob <- list()
for(i in 1:nspecies){
  rndpts_lin <- rndpts %>%
    mutate(linpred = input$detection$fixed.effect$intercept[i] + input$detection$fixed.effect$betacov[i]* extract(cov3.rast,rndpts))
  det_prob[[i]] <- plogis(rndpts_lin$linpred)
}
#data_det <- rbinom(length(det_prob),1,det_prob)

detection_data <- list()
for(i in 1:nspecies){
  data_det <- vector("numeric", length(det_prob[[i]]))
  for(j in 1:length(det_prob[[i]])){
    data_det[j]<- rbinom(1,1, det_prob[[i]][j])
    detection_data[[i]] <- data_det
  }
}
#Organising as spatial dataframe
data_det_spframe <- list()
for(i in 1:nspecies){
  data_det_spframe[[i]] <- SpatialPointsDataFrame(rndpts, data = data.frame(detection_data[[i]]))
  names(data_det_spframe[[i]])<-paste0("detdata",i)
}


## Fit the model using inlabru ##

## the borders of the study region 
coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
poly <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))

## the mesh
mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
                     max.edge = c(0.1, 0.5), cutoff = 0.2)

## SPDEs definition
spdes <- list()
for(i in 1: nspecies){
  spdes[[i]] <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0.05) = 0.01
                                    prior.range = c(input$ecological$hyperparameters$range[i], 0.5),
                                    # PC-prior on sigma: P(sigma > 1) = 0.01
                                    prior.sigma = c(sqrt(input$ecological$hyperparameters$sigma2[i]), 0.5))
}

#SPDEs for the thinning
spde2 <- inla.spde2.pcmatern(mesh = mesh,
                             # PC-prior on range: P(practic.range < 0.05) = 0.01
                             prior.range = c(input$sampling$hyperparameters$range, 0.5),
                             # PC-prior on sigma: P(sigma > 1) = 0.01
                             prior.sigma = c(sqrt(input$sampling$hyperparameters$sigma2), 0.5))

csdata = simulateddata$thirdstage
cssampdata = simulateddata$firststage$Samp_PPFinal
detdata = data_det_spframe
covslist <- list(cov1.spix,cov2.spix,cov3.spix)
spdeslist <- list(spdes=spdes,spde2=spde2)
covs = covslist
region=poly
mesh=mesh

data_df <- data.frame(
Y = csdata$classifications$error,
C = csdata$classifications$true_species,
eco_cov = extract(cov1.rast,csdata$classifications),
samp_cov= extract(cov2.rast,csdata$classifications),
det_cov = extract(cov3.rast,csdata$classifications))
                    
                     
source("estpar.R")
fit2 <- est_par(omega=class_prob)

## prediction and results ##

predpoints <- expand.grid(x=seq(0,3,length.out = 128),y=seq(0,3,length.out = 128))
cov1.pred <- cos(predpoints$x) - sin(predpoints$y - 2)
cov2.pred <- cos(2*predpoints$x) - sin(2*predpoints$y-4)
cov3.pred <- (predpoints$x/2)^2+(predpoints$y/2)^2

pred.median.eco <- list()
pred.sd.eco <- list()
pred.median.samp <- list()
pred.sd.samp <- list()
pred.median.det <- list()
pred.sd.det <- list()

for(i in 1:nspecies){
  predpoints <- SpatialPointsDataFrame(predpoints,data=data.frame(cov1.pred,cov2.pred,cov3.pred))
  names(predpoints) <- c(paste0("cov1",i),"cov2",paste0("cov3",i))
  
  ## Predict ecological process ##
  form_eco <- as.formula(paste0("~ beta0",i,"+cov1",i,"+w1",i))
  pr_eco = predict(fit2, predpoints,form_eco)
  r <- raster(pr_eco)
  r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
  pred.median.eco[[i]] <- rasterize(pr_eco@coords,r1,pr_eco$median, fun=mean,na.rm=T)
  pred.sd.eco[[i]] <- rasterize(pr_eco@coords,r1,pr_eco$sd, fun=mean,na.rm=T)
  
  ## Predict sampling process ##
  pr_samp = predict(fit2, predpoints,~beta0thin + cov2 + w2)
  r <- raster(pr_samp)
  r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
  pred.median.samp[[i]] <- rasterize(pr_samp@coords,r1,pr_samp$median, fun=mean,na.rm=T)
  pred.sd.samp[[i]] <- rasterize(pr_samp@coords,r1,pr_samp$sd, fun=mean,na.rm=T)
  
  ## Predict detection ##
  form_det <- as.formula(paste0("~ beta0det",i,"+cov3",i))
  pr_det = predict(fit2, predpoints,form_det)
  r <- raster(pr_det)
  r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
  pred.median.det[[i]] <- rasterize(pr_det@coords,r1,pr_det$median, fun=mean,na.rm=T)
  pred.sd.det[[i]] <- rasterize(pr_det@coords,r1,pr_det$sd, fun=mean,na.rm=T)
  
}

## Posteriors for the ecological process ##
par(mfrow=c(2,2))
for(i in 1:nspecies){
  ## Plotting ecological process
  
  ##Intercept
  plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0",i))),type="l")
  abline(v=input$ecological$fixed.effect$intercept[i],col="blue",lwd=2,lty=2)
  ##Covariate
  plot(eval(parse(text=paste0("fit2$marginals.fixed$cov1",i))),type="l")
  abline(v=input$ecological$fixed.effect$betacov[i],col="blue",lwd=2,lty=2)
  ##Range
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$",paste0("'Range for w1",i,"'")))),type="l",xlim=c(0,5))
  abline(v=input$ecological$hyperparameters$range[i],col="blue",lwd=2,lty=2)
  ##SD
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$",paste0("'Stdev for w1",i,"'")))),type="l",xlim=c(0,2))
  abline(v=sqrt(input$ecological$hyperparameters$sigma2[i]),col="blue",lwd=2,lty=2)
  
}

## Posteriors for the sampling process ##
par(mfrow=c(2,2))
#for(i in 1:nspecies){
##Intercept
plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0thin"))),type="l")
abline(v=input$sampling$fixed.effect$intercept,col="blue",lwd=2,lty=2)
##Covariate
plot(eval(parse(text=paste0("fit2$marginals.fixed$cov2"))),type="l")
abline(v=input$sampling$fixed.effect$betacov,col="blue",lwd=2,lty=2)
##Range
plot(eval(parse(text=paste0("fit2$marginals.hyperpar$`Range for w2`"))),type="l",xlim=c(0,5))
abline(v=input$sampling$hyperparameters$range,col="blue",lwd=2,lty=2)
##SD
plot(eval(parse(text=paste0("fit2$marginals.hyperpar$`Stdev for w2`"))),type="l",xlim=c(0,2))
abline(v=sqrt(input$sampling$hyperparameters$sigma2),col="blue",lwd=2,lty=2)

#}


## Posteriors for the detection process ##
par(mfrow=c(1,2))
for(i in 1:nspecies){
  ##Intercept
  plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0det",i))),type="l")
  abline(v=input$detection$fixed.effect$intercept[i],col="blue",lwd=2,lty=2)
  ##Covariate
  plot(eval(parse(text=paste0("fit2$marginals.fixed$cov3",i))),type="l")
  abline(v=input$detection$fixed.effect$betacov[i],col="blue",lwd=2,lty=2)
  
}


## Let's plot the predictions ##
par(mfrow=c(4,2),mar=c(0,0,0,0))
for(i in 1:nspecies){
  ## Ecological process ##
  plot(simulateddata$species_raster[[i]])
  plot(pred.median.eco[[i]])
}
par(mfrow=c(1,2),mar=c(0,0,0,0))
for(i in 1:nspecies){
  ## Sampling process ##
  plot(log(simulateddata$firststage$retainprobraster/(1-simulateddata$firststage$retainprobraster)))
  plot(pred.median.samp[[i]])
}

for(i in 1:nspecies){
  ## Sampling process ##
  plot(input$detection$fixed.effect$intercept[i] + input$detection$fixed.effect$betacov[i]*cov3.rast )
  plot(pred.median.det[[i]])
}
