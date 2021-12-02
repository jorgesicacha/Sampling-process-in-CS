### Making copy effect between bernoulli response and PP ###
#setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds")
## Simulating the data ##
library(devtools)
#devtools::install_github("fbachl/inlabru", ref="devel", dependencies = TRUE)
library(dplyr)

setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds")

source("functionssimu.R")
###################
## Simulating w1 ##
###################
#RFoptions(seed=seed)
## Defining a window for the model ##
coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
win <- as.owin(aa)

## Gaussian Random Field for the ecological process
seed=1036610602
RFoptions(seed=seed)

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

# Mesh for the GRF
mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
                     max.edge = c(0.1, 0.5), cutoff = 0.2)

## Mesh for the Gaussian RF
spdes <- list()
for(i in 1: nspecies){
  spdes[[i]] <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0.05) = 0.01
                                    prior.range = c(input$range[i], 0.5),
                                    # PC-prior on sigma: P(sigma > 1) = 0.01
                                    prior.sigma = c(sqrt(input$sigma2[i]), 0.5))
}

#SPDEs for the thinning
spde2 <- inla.spde2.pcmatern(mesh = mesh,
                             # PC-prior on range: P(practic.range < 0.05) = 0.01
                             prior.range = c(thin_input$range_thin, 0.5),
                             # PC-prior on sigma: P(sigma > 1) = 0.01
                             prior.sigma = c(sqrt(thin_input$sigma2x_thin), 0.5))

x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = 100)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = 100)
gridlocs <- expand.grid(x0,y0)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
covariate_thin.im <- im(gridcov_thin, x0, y0)

gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_detect.im <- im(gridcov_det, x0, y0)

par(mfrow=c(2,2))
plot(covariate.im)
plot(covariate_thin.im)
plot(covariate_detect.im)


#EcologicL Process
Eco_PP <- list()
for(i in 1:nspecies){
  Eco_PP[[i]] <- rLGCP(model="matern",mu=im(input$beta0[i] + input$betacov[i]*gridcov ,xcol=x0,yrow=y0),
                       var=input$sigma2[i],scale=input$range[i]/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))
}

# Storing the raster of the true intensity for each species
species_rast <- list()
for(i in 1:nspecies){
  Lam <- attr(Eco_PP[[i]], 'Lambda')
  
  Eco_GRF  <- log(Lam$v)
  
  df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v)))))))
  r <- raster(df.sp)
  r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
  w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
  w1.rastaa <- crop(w1.rast,aa)
  species_rast[[i]] <- mask(w1.rastaa,aa)
}

par(mfrow=c(2,2),mar=c(1,1,1,1))
lapply(1:4,function(x){plot(species_rast[[x]],axes=F,box=F,main=paste0("True Ecological State for species",x))
  points(Eco_PP[[x]],pch=19,cex=0.5)
  })


# First stage of thinning (sampling process)
Samp_PP <- rLGCP(model="matern",mu=im(thin_input$mean_thin + thin_input$betacov_thin*gridcov_thin ,xcol=x0,yrow=y0),
                 var=thin_input$sigma2x_thin,scale=thin_input$range_thin/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))



Lam <- attr(Samp_PP, 'Lambda')

Samp_GRF  <- log(Lam$v)
#Retprobs <- 1-exp(-Lam$v)

#Retprobs.im <- im(Retprobs,xcol = Lam$xcol,yrow = Lam$yrow)

df.sp2 <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                                                     ew=c(anti_t(rotate(rotate((Lam$v)))))))

#df.sp2$retprob <- 1-exp(-df.sp2$ew)
df.sp2$retprob <- psych::logistic(df.sp2$w)
r <- raster(df.sp2)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
w2.rast <- rasterize(df.sp2@coords,r1,df.sp2$w, fun=mean,na.rm=T)
ew2.rast <- rasterize(df.sp2@coords,r1,df.sp2$ew, fun=mean,na.rm=T)
prob.rast <- rasterize(df.sp2@coords,r1,df.sp2$retprob, fun=mean,na.rm=T)
par(mfrow=c(1,3))
plot(w2.rast)
plot(ew2.rast)
plot(prob.rast)
#points(Eco_PP,pch=19,cex=0.2)

# cropping the thinning probability to the study region
w2.rastaa <- crop(w2.rast,aa)
w2.rastaa <- mask(w2.rastaa,aa)
prob.rastaa <- crop(prob.rast,aa)
prob.rastaa <- mask(prob.rastaa,aa)

## Getting probability of retaining samples ##
Eco_PPFinal <- list()
#retprobslik1 <- interp.im(Retprobs.im, x = Eco_PP[[i]]$x,y = Eco_PP$y)
for(i in 1:nspecies){
  retprobslik <- extract(prob.rast,cbind(Eco_PP[[i]]$x,Eco_PP[[i]]$y))
  Eco_PP1 <- SpatialPointsDataFrame(coords = cbind(Eco_PP[[i]]$x,Eco_PP[[i]]$y),data=data.frame(retprob=retprobslik))
  Eco_PP1$retain <- apply(Eco_PP1@data, 1,function(x){rbinom(1,1,p=x[1])})
  
  Eco_PPFinal[[i]] <- Eco_PP1[which(Eco_PP1$retain==1),]
}

# data from where CS sample
Samp_PPFinal <- SpatialPoints(coords=cbind(Samp_PP$x,Samp_PP$y))

par(mfrow=c(2,2),mar=c(1,1,1,1))
lapply(1:4,function(x){plot(prob.rastaa,axes=F,box=F,main=paste0("Retaining probability for species ",x))
  points(Eco_PP[[x]],pch=19,cex=0.5)
  points(Eco_PPFinal[[x]],pch=19,cex=0.3,col="red")
})


##### Detection Thinning##############  #
df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(etadetect=c(anti_t(rotate(rotate(covariate_detect.im$v))))))
r <- raster(df.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
etadetect.rast <- rasterize(df.sp@coords,r1,df.sp$etadetect, 
                            fun=mean,na.rm=T)
plot(etadetect.rast)

# Computing detection probabilities
det_PP <- list()
for(i in 1:nspecies){
  det_PP[[i]] <-  psych::logistic(input$alpha0[i]+input$alphacov[i]*extract(etadetect.rast,Eco_PPFinal[[i]]))
}
#if(i==1) assign(paste0("detectprob",i),
#               psych::logistic(get(paste0("alpha0_detectsp",i))+
#                                 get(paste0("alpha1_detectsp",i))*extract(etadetect.rast,get("Eco_PPFinal"))))
#else assign(paste0("detectprob",i),
#           psych::logistic(get(paste0("alpha0_detectsp",i))+
#                            get(paste0("alpha1_detectsp",i))*extract(etadetect.rast,get(paste0("Eco_PPFinal_sp",i)))))
#}

#Point pattern from thinning due to detection
Eco_PPFinal_detect <- list()
for(i in 1:nspecies){
  Eco_PPFinal[[i]]$retain_detect <- sapply(1:length(Eco_PPFinal[[i]]),
                                           function(x){rbinom(1,1,p=det_PP[[i]][x])})
  Eco_PPFinal_detect[[i]] <- Eco_PPFinal[[i]][which(Eco_PPFinal[[i]]$retain_detect==1),]
}

## Making rasters of detection probability ##

etadetect.rastaa <- crop(etadetect.rast,aa)
detrasts <- list()


for(i in 1:nspecies){
  detrasts[[i]] <- psych::logistic(input$alpha0[i]+input$alphacov[i]*etadetect.rastaa)
}

par(mfrow=c(2,2),mar=c(1,1,1,1))
lapply(1:4,function(x){plot(detrasts[[x]],axes=F,box=F,main=paste0("Detection probability for species ",x))
  points(Eco_PPFinal[[x]],pch=19,cex=0.5)
  points(Eco_PPFinal_detect[[x]],pch=19,cex=0.3,col="red")
  
})


### Sampling the detections from a survey
rndpts_x0 <- runif(50, 0,3)
rndpts_y0 <- runif(50, 0,3)
rndpts <- data.frame(rndpts_x0, rndpts_y0)

det_prob <- list()
for(i in 1:nspecies){
  rndpts_lin <- rndpts %>%
    mutate(linpred = input$alpha0[i] + input$alphacov[i]* extract(etadetect.rast,rndpts))
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

# Third stage of Thinning
#Chassification probability
class_prob <- matrix(c(0.9, 0.02, 0.04, 0.04,
                       0.05, 0.89, 0.04, 0.02,
                       0.1,0.1, 0.8, 0,
                       0, 0.05, 0.25, 0.7),
                     nrow=4, ncol=4, byrow = TRUE
)

for(i in 1:nspecies){
  Eco_PPFinal_detect[[i]]$error <- apply(Eco_PPFinal_detect[[i]]@data, 1,function(x){which(rmultinom(1,1,p=(class_prob[i,]))==1)})
  Eco_PPFinal_detect[[i]]$true_species <- i
}

#Putting all the data together
cit_data <- list()
for(i in 1:nspecies){
  cit_data[[i]] <- do.call("rbind",lapply(Eco_PPFinal_detect, function(x) x[which(x@data$error ==i),]))
}

allcit_data <- do.call("rbind",cit_data)
par(mfrow=c(2,2),mar=c(1,1,1,1))
lapply(1:4,function(x){
  if(x==1){
    plot(aa,axes=F,box=F,main=paste0("A:CS reports for species ",x))
    points(allcit_data[which(allcit_data$true_species==x),],pch=19,cex=1,col="blue")
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==2),],pch=15,cex=0.6,col="green"))
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==3),],pch=17,cex=0.6,col="yellow"))
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==4),],pch=18,cex=0.6,col="red"))
    }
  else{if(x==2){
    plot(aa,axes=F,box=F,main=paste0("B:CS reports for species ",x))
    points(allcit_data[which(allcit_data$true_species==x),],pch=15,cex=1,col="green")
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==1),],pch=19,cex=0.6,col="blue"))
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==3),],pch=17,cex=0.6,col="yellow"))
    try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==4),],pch=18,cex=0.6,col="red"))
    }
    else{
      if(x==3){
        plot(aa,axes=F,box=F,main=paste0("C:CS reports for species ",x))
        points(allcit_data[which(allcit_data$true_species==x),],pch=17,cex=1,col="yellow")
        try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==1),],pch=19,cex=0.6,col="blue"))
        try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==2),],pch=15,cex=0.6,col="green"))
        try( points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==4),],pch=18,cex=0.6,col="red"))
        }
      else{
        plot(aa,axes=F,box=F,main=paste0("D:CS reports for species ",x))
        points(allcit_data[which(allcit_data$true_species==x),],pch=18,cex=1,col="red")
        try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==2),],pch=15,cex=0.6,col="green"))
        try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==3),],pch=17,cex=0.6,col="yellow"))
        try(points(allcit_data[which(allcit_data$true_species==x & allcit_data$error==4),],pch=18,cex=0.6,col="red"))
        }
    }}
  
  # points(cit_data[[x]][which(cit_data[[x]]$true_species==1),],pch=19,cex=0.6,col="blue")
  # try(points(cit_data[[x]][which(cit_data[[x]]$true_species==2),],pch=15,cex=0.6,col="green"))
  # try(points(cit_data[[x]][which(cit_data[[x]]$true_species==3),],pch=17,cex=0.6,col="yellow"))
  # try(points(cit_data[[x]][which(cit_data[[x]]$true_species==4),],pch=18,cex=0.6,col="red"))
  # 
})







####### Fitting the model ##########
#Covariates for true intensity
cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate.im$v))))))
r <- raster(cov1.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
plot(cov1.rast)
cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")

#Covariates for first thinning
cov2.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_thin.im$v))))))
r <- raster(cov2.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov2.rast <- rasterize(cov2.sp@coords,r1,cov2.sp$cov, fun=mean,na.rm=T)
plot(cov2.rast)
cov2.spix <- as(cov2.rast,"SpatialPixelsDataFrame")

#Covariate for second thinning
cov3.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_detect.im$v))))))
r <- raster(cov3.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov3.rast <- rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
plot(cov3.rast)
cov3.spix <- as(cov3.rast,"SpatialPixelsDataFrame")



est_par <- function(omega){
#Defining components of the model
cmp1 <- list()
for(i in 1:nspecies){
  cmp1[[i]] <- (paste0("+ beta0",i,"(1)", 
                       "+ beta0thin(1)",
                       "+ beta0det" ,i, "(1)",
                       "+w1",i,"(map = coordinates, model =","spdes[[",i, "]])", 
                       "+ w2(map = coordinates, model = spde2)+",
                       "cov1",i, "(map=cov1.spix,model='linear') +",
                       "cov2(map=cov2.spix,model='linear')+",
                       "cov3",i, "(map=cov3.spix,model='linear')"))
}
cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))

fun <- function(x,y,z){
  -log(1+exp((x+y+z)))
}

fun1 <- function(x,y){
  -log(1+exp((x+y)))
}



fun21 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
  ret <- log(omega[1,1]*plogis(a+b+c) + omega[2,1]*plogis(d+e+f) + 
               omega[3,1]*plogis(g+h+i) + omega[4,1]*plogis(j+k+l))
  #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
  return(ret)
}

fun22 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
  ret <- log(omega[1,2]*plogis(a+b+c) + omega[2,2]*plogis(d+e+f) + 
               omega[3,2]*plogis(g+h+i) + omega[4,2]*plogis(j+k+l))
  return(ret)
}

fun23 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
  ret <- log(omega[1,3]*plogis(a+b+c) + omega[2,3]*plogis(d+e+f) + 
               omega[3,3]*plogis(g+h+i) + omega[4,3]*plogis(j+k+l))
  #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
  return(ret)
}

fun24 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
  ret <- log(omega[1,4]*plogis(a+b+c) + omega[2,4]*plogis(d+e+f) + 
               omega[3,4]*plogis(g+h+i) + omega[4,4]*plogis(j+k+l))
  #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
  return(ret)
}

library(inlabru)


lik1 <- lik2 <- lik3 <- list()

for(i in 1:nspecies){
  lik1[[i]] <- like("cp",
                    formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",i,"+cov3",i,"+fun1(beta0det",i,", cov3",i,")",
                                                "+fun2",i,"(beta01, cov11, w11, beta02, cov12, w12,beta03, cov13, w13,beta04, cov14, w14)")),
                    data = Eco_PPFinal_detect[[i]],
                    #components = cmp,
                    domain = list(coordinates = mesh),
                    samplers = aa)
  lik2[[i]] <- like("cp",
                    formula = coordinates ~ beta0thin + cov2 + w2,
                    data = Samp_PPFinal,
                    #components = cmp,
                    domain = list(coordinates = mesh),
                    samplers = aa)
  lik3[[i]] <- like("binomial",
                    formula = as.formula(paste0("detdata",i," ~ beta0det",i," + cov3",i)),
                    data = data_det_spframe[[i]],
                    #components = cmp,
                    domain = list(coordinates = mesh),
                    samplers = aa)
}


#str(lik1[[1]])


# 
# 

predpoints <- expand.grid(x=seq(0,3,length.out = 128),y=seq(0,3,length.out = 128))
cov1.pred <- cos(predpoints$x) - sin(predpoints$y - 2)
cov2.pred <- cos(2*predpoints$x) - sin(2*predpoints$y-4)
cov3.pred <- (predpoints$x/2)^2+(predpoints$y/2)^2


# par(mfrow=c(1,3))
# plot(pred.rast.median.bru)#,zlim=c(0.9,4.1))
# plot(pred.rast.sd.bru)
# plot(w1.rast)

inlabru:::iinla.setOption("iinla.verbose", TRUE)
#fit2 <- list()
pred.median.eco <- list()
pred.sd.eco <- list()
pred.median.samp <- list()
pred.sd.samp <- list()
pred.median.det <- list()
pred.sd.det <- list()

cov1.spix$cov11 <- cov1.spix$cov12 <- cov1.spix$cov13 <- cov1.spix$cov14 <- cov1.spix$layer
cov3.spix$cov31 <- cov3.spix$cov32 <- cov3.spix$cov33 <- cov3.spix$cov34 <- cov3.spix$layer
#for(i in 1:nspecies){
# names(cov1.spix)<- paste0("cov1",i)
#nemes(cov1.spix) <- "cov11"
names(cov2.spix)<- "cov2"
#names(cov3.spix) <- paste0("cov3",i)
fit2 <- bru(cmp, lik1[[1]], lik1[[2]],lik1[[3]],lik1[[4]],
            lik2[[1]],
            lik3[[1]],lik3[[2]],lik3[[3]],lik3[[4]],
            options = list(control.inla = list(strategy = "gaussian",
                                               int.strategy = "eb"),
                           max.iter=50))
return(fit2)
}
fit2 <- est_par(class_prob)

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
save.image("output1903.Rdata")
load("output1903.Rdata")



## Posteriors for the ecological process ##
par(mfrow=c(2,2))
for(i in 1:nspecies){
  ## Plotting ecological process
  
  ##Intercept
  plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0",i))),type="l")
  abline(v=input$beta0[i],col="blue",lwd=2,lty=2)
  ##Covariate
  plot(eval(parse(text=paste0("fit2$marginals.fixed$cov1",i))),type="l")
  abline(v=input$betacov[i],col="blue",lwd=2,lty=2)
  ##Range
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$",paste0("'Range for w1",i,"'")))),type="l",xlim=c(0,5))
  abline(v=input$range[i],col="blue",lwd=2,lty=2)
  ##SD
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$",paste0("'Stdev for w1",i,"'")))),type="l",xlim=c(0,2))
  abline(v=sqrt(input$sigma2[i]),col="blue",lwd=2,lty=2)
  
}

## Posteriors for the sampling process ##
par(mfrow=c(2,2))
#for(i in 1:nspecies){
  ##Intercept
  plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0thin"))),type="l")
  abline(v=thin_input$mean_thin,col="blue",lwd=2,lty=2)
  ##Covariate
  plot(eval(parse(text=paste0("fit2$marginals.fixed$cov2"))),type="l")
  abline(v=thin_input$betacov_thin,col="blue",lwd=2,lty=2)
  ##Range
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$`Range for w2`"))),type="l",xlim=c(0,5))
  abline(v=thin_input$range_thin,col="blue",lwd=2,lty=2)
  ##SD
  plot(eval(parse(text=paste0("fit2$marginals.hyperpar$`Stdev for w2`"))),type="l",xlim=c(0,2))
  abline(v=sqrt(thin_input$sigma2x_thin),col="blue",lwd=2,lty=2)
  
#}


## Posteriors for the detection process ##
par(mfrow=c(1,2))
for(i in 1:nspecies){
  ##Intercept
  plot(eval(parse(text=paste0("fit2$marginals.fixed$beta0det",i))),type="l")
  abline(v=input$alpha0[i],col="blue",lwd=2,lty=2)
  ##Covariate
  plot(eval(parse(text=paste0("fit2$marginals.fixed$cov3",i))),type="l")
  abline(v=input$alphacov[i],col="blue",lwd=2,lty=2)
  
}


## Let's plot the predictions ##
par(mfrow=c(4,2),mar=c(0,0,0,0))
for(i in 1:nspecies){
  ## Ecological process ##
  plot(species_rast[[i]])
  plot(pred.median.eco[[i]])
}
par(mfrow=c(1,2),mar=c(0,0,0,0))
for(i in 1:nspecies){
  ## Sampling process ##
  plot(w2.rastaa)
  plot(pred.median.samp[[i]])
}

for(i in 1:nspecies){
  ## Sampling process ##
  plot(input$alpha0[i] + input$alphacov[i]*etadetect.rast )
  plot(pred.median.det[[i]])
}



