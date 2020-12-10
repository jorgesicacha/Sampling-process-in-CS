### Making copy effect between bernoulli response and PP ###

## Simulating the data ##

setwd("/Users/jorgespa/Documents/Research/DataIntegration/DeadBirds")
source("functionssimu.R")
#################################
## Simulating w1 for 4 species ##
#################################
#RFoptions(seed=seed)
## Defining a window for the model ##
coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
win <- as.owin(aa)
## Gaussian Random Field for the ecological process
seed=1036610602
RFoptions(seed=seed)
## For the species we're interested in ##
sigma2x <- 0.2
range <- 1.2
meanf <- 3
betacov <- -0.5
mean_thin <- 0.3
betacov_thin <- -1
sigma2x_thin <- 2.5
range_thin <- 0.35
## For the second species ##
sigma2xsp2 <- 1.2
rangesp2 <- 2.5
meanfsp2 <- 2.5
betacovsp2 <- -0.12
## For the third species ##
sigma2xsp3 <- 2
rangesp3 <- 3.2
meanfsp3 <- -0.1
betacovsp3 <- 0.89
## For the fourth species ##
sigma2xsp4 <- 0.1
rangesp4 <- 0.22
meanfsp4 <- 1.2
betacovsp4 <- -0.4

mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
                     max.edge = c(0.1, 0.5), cutoff = 0.2)

x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = 1000)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = 1000)
gridlocs <- expand.grid(x0,y0)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

gridcov_thin <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_thin.im <- im(gridcov_thin, x0, y0)

Eco_PP <- rLGCP(model="matern",mu=im(meanf + betacov*gridcov ,xcol=x0,yrow=y0),
                var=sigma2x,scale=range/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))

Lam <- attr(Eco_PP, 'Lambda')

Eco_GRF  <- log(Lam$v)


df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v)))))))
r <- raster(df.sp)
r1<-disaggregate(r, fact=res(r)/c(0.0056,0.0056))
w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
plot(w1.rast)

Eco_PP_sp2 <- rLGCP(model="matern",mu=im(meanfsp2 + betacovsp2*gridcov ,xcol=x0,yrow=y0),
                var=sigma2xsp2,scale=rangesp2/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))

Eco_PP_sp3 <- rLGCP(model="matern",mu=im(meanfsp3 + betacovsp3*gridcov ,xcol=x0,yrow=y0),
                    var=sigma2xsp3,scale=rangesp3/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))

Eco_PP_sp4 <- rLGCP(model="matern",mu=im(meanfsp4 + betacovsp4*gridcov ,xcol=x0,yrow=y0),
                    var=sigma2xsp4,scale=rangesp4/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))


plot(aa)
points(Eco_PP,pch=19,cex=1,col="blue")
points(Eco_PP_sp2,pch=15,cex=1,col="green")
points(Eco_PP_sp3,pch=17,cex=1,col="yellow")
points(Eco_PP_sp4,pch=18,cex=1,col="red")

Samp_PP <- rLGCP(model="matern",mu=im(mean_thin + betacov_thin*gridcov_thin ,xcol=x0,yrow=y0),
                 var=sigma2x_thin,scale=range_thin/sqrt(8),nu=1,win = win,xy=list(x=x0,y=y0))

Lam <- attr(Samp_PP, 'Lambda')

Samp_GRF  <- log(Lam$v)

df.sp2 <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                                                ew=c(anti_t(rotate(rotate((Lam$v)))))))

df.sp2$retprob <- psych::logistic(df.sp2$w)
r <- raster(df.sp2)
r1<-disaggregate(r, fact=res(r)/c(0.0056,0.0056))
w2.rast <- rasterize(df.sp2@coords,r1,df.sp2$w, fun=mean,na.rm=T)
ew2.rast <- rasterize(df.sp2@coords,r1,df.sp2$ew, fun=mean,na.rm=T)
prob.rast <- rasterize(df.sp2@coords,r1,df.sp2$retprob, fun=mean,na.rm=T)
#par(mfrow=c(1,3))
#plot(w2.rast)
#plot(ew2.rast)
#plot(prob.rast)
#points(Eco_PP,pch=19,cex=0.2)

## Getting probability of retaining samples ##
retprobslik1 <- extract(prob.rast,cbind(Eco_PP$x,Eco_PP$y))
Eco_PP <- SpatialPointsDataFrame(coords = cbind(Eco_PP$x,Eco_PP$y),data=data.frame(retprob=retprobslik1))
Eco_PP$retain <- apply(Eco_PP@data, 1,function(x){rbinom(1,1,p=x[1])})
Eco_PPFinal <- Eco_PP[which(Eco_PP$retain==1),]

retprobslik1_sp2 <- extract(prob.rast,cbind(Eco_PP_sp2$x,Eco_PP_sp2$y))
Eco_PP_sp2 <- SpatialPointsDataFrame(coords = cbind(Eco_PP_sp2$x,Eco_PP_sp2$y),data=data.frame(retprob=retprobslik1_sp2))
Eco_PP_sp2$retain <- apply(Eco_PP_sp2@data, 1,function(x){rbinom(1,1,p=x[1])})
Eco_PPFinal_sp2 <- Eco_PP_sp2[which(Eco_PP_sp2$retain==1),]

retprobslik1_sp3 <- extract(prob.rast,cbind(Eco_PP_sp3$x,Eco_PP_sp3$y))
Eco_PP_sp3 <- SpatialPointsDataFrame(coords = cbind(Eco_PP_sp3$x,Eco_PP_sp3$y),data=data.frame(retprob=retprobslik1_sp3))
Eco_PP_sp3$retain <- apply(Eco_PP_sp3@data, 1,function(x){rbinom(1,1,p=x[1])})
Eco_PPFinal_sp3 <- Eco_PP_sp3[which(Eco_PP_sp3$retain==1),]

retprobslik1_sp4 <- extract(prob.rast,cbind(Eco_PP_sp4$x,Eco_PP_sp4$y))
Eco_PP_sp4 <- SpatialPointsDataFrame(coords = cbind(Eco_PP_sp4$x,Eco_PP_sp4$y),data=data.frame(retprob=retprobslik1_sp4))
Eco_PP_sp4$retain <- apply(Eco_PP_sp4@data, 1,function(x){rbinom(1,1,p=x[1])})
Eco_PPFinal_sp4 <- Eco_PP_sp4[which(Eco_PP_sp4$retain==1),]

plot(aa)
points(Eco_PPFinal,pch=19,cex=1,col="blue")
points(Eco_PPFinal_sp2,pch=15,cex=1,col="green")
points(Eco_PPFinal_sp3,pch=17,cex=1,col="yellow")
points(Eco_PPFinal_sp4,pch=18,cex=1,col="red")

legend("right",legend = c("Species 1","Species 2", "Species 3","Species 4"),
       col = c("blue","green","yellow","red"),pch = c(19,15,17,18)
              )

