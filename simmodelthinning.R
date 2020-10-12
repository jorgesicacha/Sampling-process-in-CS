### Making copy effect between bernoulli response and PP ###

## Simulating the data ##

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
sigma2x <- 0.2
range <- 1.2
meanf <- 3
betacov <- -0.5
mean_thin <- 0.3
betacov_thin <- -1
sigma2x_thin <- 2.5
range_thin <- 0.35


mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
                     max.edge = c(0.1, 0.5), cutoff = 0.2)

spde1 <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < 0.05) = 0.01
                            prior.range = c(1.2, 0.5),
                            # PC-prior on sigma: P(sigma > 1) = 0.01
                            prior.sigma = c(0.45, 0.5))

spde2 <- inla.spde2.pcmatern(mesh = mesh,
                             # PC-prior on range: P(practic.range < 0.05) = 0.01
                             prior.range = c(2.5, 0.5),
                             # PC-prior on sigma: P(sigma > 1) = 0.01
                             prior.sigma = c(1.58, 0.5))

x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = 100)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = 100)
gridlocs <- expand.grid(x0,y0)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

gridcov_thin <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_thin.im <- im(gridcov_thin, x0, y0)

par(mfrow=c(1,2))
plot(covariate.im)
plot(covariate_thin.im)

Eco_PP <- rLGCP(model="matern",mu=im(meanf + betacov*gridcov ,xcol=x0,yrow=y0),
                 var=sigma2x,scale=range/sqrt(8),nu=1,win = win)

Lam <- attr(Eco_PP, 'Lambda')

Eco_GRF  <- log(Lam$v)

xs <- Lam$xcol
ys <- Lam$yrow
xys <- expand.grid(xs,ys)

df.sp <- SpatialPointsDataFrame(coords = xys,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v)))))))
r <- raster(df.sp)
r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
plot(w1.rast)

Samp_PP <- rLGCP(model="matern",mu=im(mean_thin + betacov_thin*gridcov_thin ,xcol=x0,yrow=y0),
                var=sigma2x_thin,scale=range_thin/sqrt(8),nu=1,win = win)



Lam <- attr(Samp_PP, 'Lambda')

Samp_GRF  <- log(Lam$v)
Retprobs <- 1-exp(-Lam$v)

Retprobs.im <- im(Retprobs,xcol = Lam$xcol,yrow = Lam$yrow)

df.sp2 <- SpatialPointsDataFrame(coords = xys,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                                                ew=c(anti_t(rotate(rotate((Lam$v)))))))

#df.sp2$retprob <- 1-exp(-df.sp2$ew)
df.sp2$retprob <- psych::logistic(df.sp2$w)
r <- raster(df.sp2)
r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
w2.rast <- rasterize(df.sp2@coords,r1,df.sp2$w, fun=mean,na.rm=T)
ew2.rast <- rasterize(df.sp2@coords,r1,df.sp2$ew, fun=mean,na.rm=T)
prob.rast <- rasterize(df.sp2@coords,r1,df.sp2$retprob, fun=mean,na.rm=T)
par(mfrow=c(1,3))
plot(w2.rast)
plot(ew2.rast)
plot(prob.rast)
points(Eco_PP,pch=19,cex=0.2)

## Getting probability of retaining samples ##
retprobslik1 <- extract(prob.rast,cbind(Eco_PP$x,Eco_PP$y))
retprobslik2 <- interp.im(Retprobs.im, 
                          x = Eco_PP$x,
                          y = Eco_PP$y)

Eco_PP <- SpatialPointsDataFrame(coords = cbind(Eco_PP$x,Eco_PP$y),data=data.frame(retprob=retprobslik2))
Eco_PP$retain <- apply(Eco_PP@data, 1,function(x){rbinom(1,1,p=x[1])})

Eco_PPFinal <- Eco_PP[which(Eco_PP$retain==1),]

Samp_PPFinal <- SpatialPoints(coords=cbind(Samp_PP$x,Samp_PP$y))

####### Fitting the model ##########
gridlocsmesh.grid <- expand.grid(x0,y0)

cov1.sp <- SpatialPointsDataFrame(coords = gridlocsmesh.grid,data = data.frame(cov=c(anti_t(rotate(rotate(covariate.im$v))))))
r <- raster(cov1.sp)
r1<-disaggregate(r, fact=res(r)/c(0.075,0.075))
cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
plot(cov1.rast)
cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")

cov2.sp <- SpatialPointsDataFrame(coords = gridlocsmesh.grid,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_thin.im$v))))))
r <- raster(cov2.sp)
r1<-disaggregate(r, fact=res(r)/c(0.075,0.075))
cov2.rast <- rasterize(cov2.sp@coords,r1,cov2.sp$cov, fun=mean,na.rm=T)
plot(cov2.rast)
cov2.spix <- as(cov1.rast,"SpatialPixelsDataFrame")


cmp<- ~ -1 + beta0 + beta0thin +
        w1(map = coordinates, model = spde1) + w2(map = coordinates, model = spde2) +
        cov1(map=cov1.spix,model="linear") + cov2(map=cov2.spix,model="linear")

fun <- function(x,y,z){
  -log(1+exp(-(x+y+z)))
}

lik1 <- like("cp",
             formula = coordinates ~ beta0  + cov1 + w1 + fun(beta0thin,cov2,w2),
             data = Eco_PPFinal,
             components = cmp,
             domain = list(coordinates = mesh),
             samplers = aa)

lik2 <- like("cp",
             formula = coordinates ~ beta0thin + cov2 + w2,
             data = Samp_PPFinal,
             components = cmp,
             domain = list(coordinates = mesh),
             samplers = aa)
# 
# 
inlabru:::iinla.setOption("iinla.verbose", TRUE)
fit2 <- bru(cmp, lik1, lik2,options = list(control.inla = list(strategy = "gaussian",
                                                               int.strategy = "eb"),
                                           max.iter=20))

fit2$summary.fixed
fit2$summary.hyperpar


par(mfrow=c(2,4))
plot(fit2$marginals.hyperpar$`Range for w1`,type="l",col="red",xlim=c(0,20))
abline(v=range,col="blue",lty=2,lwd=2)
plot(fit2$marginals.hyperpar$`Stdev for w1`,type="l",col="red",xlim=c(0,2))
abline(v=sqrt(sigma2x),col="blue",lty=2,lwd=2)
plot(fit2$marginals.fixed$beta0,type="l",col="red")
abline(v=meanf,col="blue",lty=2,lwd=2)
plot(fit2$marginals.fixed$cov1,type="l",col="red")
abline(v=betacov,col="blue",lty=2,lwd=2)

plot(fit2$marginals.hyperpar$`Range for w2`,type="l",col="red",xlim=c(0,20))
abline(v=range_thin,col="blue",lty=2,lwd=2)
plot(fit2$marginals.hyperpar$`Stdev for w2`,type="l",col="red",xlim=c(0,2))
abline(v=sqrt(sigma2x_thin),col="blue",lty=2,lwd=2)
plot(fit2$marginals.fixed$beta0thin,type="l",col="red")
abline(v=mean_thin,col="blue",lty=2,lwd=2)
plot(fit2$marginals.fixed$cov2,type="l",col="red")
abline(v=betacov_thin,col="blue",lty=2,lwd=2)


predpoints <- expand.grid(x=seq(0,3,length.out = 128),y=seq(0,3,length.out = 128))
cov1.pred <- interp.im(covariate.im, 
                      x = predpoints[,1],
                      y = predpoints[,2])
predpoints <- SpatialPointsDataFrame(predpoints,data=data.frame(cov1=cov1.pred))

pr = predict(fit2, predpoints, ~ beta0 + cov1 + w1)
r <- raster(pr)
r1<-disaggregate(r, fact=res(r)/c(0.025,0.025))
pred.rast.median.bru <- rasterize(pr@coords,r1,pr$median, fun=mean,na.rm=T)
pred.rast.sd.bru <- rasterize(pr@coords,r1,pr$sd, fun=mean,na.rm=T)
par(mfrow=c(1,3))
plot(pred.rast.median.bru)#,zlim=c(0.9,4.1))
plot(pred.rast.sd.bru)
plot(w1.rast)


