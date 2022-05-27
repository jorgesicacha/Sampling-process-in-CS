## Data generating process plots ##
library(spatstat)

## Simulating the Covariates ##

#Grid for simulation
x0 <- seq(-1.3, 4.3, length = 100)
y0 <- seq(-1.3,4.3, length = 100)
gridlocs <- expand.grid(x0,y0)

# Covariates for true ecological state
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

#Covariate for the sampling process
gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
covariate_thin.im <- im(gridcov_thin, x0, y0)

#Covariate for the detection
gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_detect.im <- im(gridcov_det, x0, y0)

# cov: a list with the covariates needed by both the data generating process and the model we want to fit.
cov <- list(covariate.im,covariate_thin.im,covariate_detect.im)

## Generating the data ##

source("DataGeneration.R")

nspecies <- 2#nspecies: Number of species we want to simulate 

#input: List with the parameters of the model that generates CS data
input <-{list(
  ecological = list(
    fixed.effect=list(
      intercept = c(0.8, 2.5),
      betacov = c(1.5, -0.12)
    ),
    hyperparameters = list(
      sigma2 = c(0.2, 1.2),
      range = c(1.2, 2.5)
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
      intercept=c(2,-0.3),
      betacov = c(-2, -0.5)
    )
  ),
  
  misclassification = list(
    
    class_prob <- matrix(c(0.9, 0.1, 
                           0.05, 0.95),
                         nrow=nspecies, ncol=nspecies, byrow = TRUE)
    
    
  )
)}

#idxs: A very simple list which tells which covariates belong to each stage of the data generating process
idxs <- list(eco=c(1),sampling=c(2),detection=c(3))

#seed: Random seed for replicating the datasets generated
seed <- 1036610602

#plot: Do we want to plot the results? Which results? ISSUE
plot = list(all=TRUE)

#Data simulation
simulateddata <-csdata(nspecies=nspecies,input=input,cov=cov,idxs=idxs,seed=seed,plot=plot)


## New figure for the manuscript
coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))


library(ggplot2)
library(raster)
library(rasterVis)
l1 <- levelplot(simulateddata$species_raster[[1]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[1]==109),
                scales=list(draw=FALSE),            
                col.regions=viridis,                   
                at=seq(-5, 5, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[1]]$x,simulateddata$trueecological[[1]]$y)),pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[1]]$x,simulateddata$trueecological[[1]]$y)),pch=19,cex=0.5,col="black",ann=FALSE,lwd=2))

l2 <- levelplot(simulateddata$species_raster[[2]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[2]==108),
                scales=list(draw=FALSE),            
                col.regions=viridis,                   
                at=seq(-0.1, 5.5, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[2]]$x,simulateddata$trueecological[[2]]$y)),pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[2]]$x,simulateddata$trueecological[[2]]$y)),pch=19,cex=0.5,col="black",ann=FALSE,lwd=2))

library(RColorBrewer)
colr <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))
l3 <- levelplot(simulateddata$firststage$retainprobraster,
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=0:1, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                scales=list(draw=FALSE),   
                main=bquote(n[1]==67),
                col.regions=colr,                   
                at=seq(0, 1, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[1]]$x,simulateddata$trueecological[[1]]$y)),pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[1]]$x,simulateddata$trueecological[[1]]$y)),pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[1]],pch=21,cex=1,col="black",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[1]],pch=19,cex=0.5,col="LimeGreen",ann=FALSE,lwd=2))


l4 <- levelplot(simulateddata$firststage$retainprobraster,
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[2]==79),
                scales=list(draw=FALSE),            
                col.regions=colr,                   
                at=seq(0, 1, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[2]]$x,simulateddata$trueecological[[2]]$y)),pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(SpatialPoints(cbind(simulateddata$trueecological[[2]]$x,simulateddata$trueecological[[2]]$y)),pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[2]],pch=21,cex=1,col="black",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[2]],pch=19,cex=0.5,col="LimeGreen",ann=FALSE,lwd=2))

l5 <- levelplot(simulateddata$secondstage$detectionprobraster[[1]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[1]==47),
                scales=list(draw=FALSE),            
                col.regions=colr,                   
                at=seq(0, 1, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[1]],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[1]],pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[1]],pch=21,cex=1,col="black",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[1]],pch=19,cex=0.5,col="LimeGreen",ann=FALSE,lwd=2))

l6 <- levelplot(simulateddata$secondstage$detectionprobraster[[2]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[2]==10),
                scales=list(draw=FALSE),            
                col.regions=colr,                   
                at=seq(0, 1, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[2]],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$firststage$Eco_PPFinal[[2]],pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[2]],pch=21,cex=1,col="black",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[2]],pch=19,cex=0.5,col="LimeGreen",ann=FALSE,lwd=2))

l7 <- levelplot(simulateddata$lambda_obs_raster[[1]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                main=bquote(n[1]==42),
                scales=list(draw=FALSE),            
                col.regions=viridis,                   
                at=seq(-5, 5, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[1]],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[1]],pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==1 & simulateddata$thirdstage$classifications$true_species==1),],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==1 & simulateddata$thirdstage$classifications$true_species==1),],pch=19,cex=0.5,col="blue",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==1 & simulateddata$thirdstage$classifications$true_species==2),],pch=23,cex=1,col="orange",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==2 & simulateddata$thirdstage$classifications$true_species==1),],pch=22,cex=1,col="red",ann=FALSE,lwd=2))



l8 <- levelplot(simulateddata$lambda_obs_raster[[2]],
                margin=FALSE,                       
                colorkey=list(
                  space='bottom',                   
                  labels=list(at=-5:5, font=4),
                  axis.line=list(col='black')       
                ),    
                par.settings=list(
                  axis.line=list(col='transparent') 
                ),
                scales=list(draw=FALSE), 
                main=bquote(n[2]==15),
                col.regions=viridis,                   
                at=seq(-5, 5, len=101)) +           
  layer(sp.polygons(aa, lwd=3))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[2]],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$secondstage$Eco_PPFinal_detect[[2]],pch=19,cex=0.5,col="grey",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==2 & simulateddata$thirdstage$classifications$true_species==2) ,],pch=21,cex=1,col="white",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==2 & simulateddata$thirdstage$classifications$true_species==2 ),],pch=19,cex=0.5,col="blue",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==2 & simulateddata$thirdstage$classifications$true_species==1),],pch=23,cex=1,col="orange",ann=FALSE,lwd=2))+
  layer(sp.points(simulateddata$thirdstage$classifications[which(simulateddata$thirdstage$classifications$error==1 & simulateddata$thirdstage$classifications$true_species==2),],pch=22,cex=1,col="red",ann=FALSE,lwd=2))

library(grid)
legd_1 <- legendGrob("True occurrences", pch=21,
                     gp=gpar(col = "white", fill = "black",cex=1))

legd_2 <- legendGrob(c("Retained observations","Deleted observations"), pch=21,
                     gp=gpar(col = c("black","white"), fill = c("LimeGreen","grey"),lwd=c(2,2),cex=c(1,1)))

legd_3 <- legendGrob(c("True Positives","False positives","False negatives"), pch=c(21,23,22),
                     gp=gpar(col = c("white","orange","red"), fill = c("blue",NA,NA),lwd=c(2,2),cex=c(1,1)))


library(cowplot)
library(grid)
library(gridExtra)
true.grob <- textGrob("True Ecological State", 
                      gp=gpar(fontface="bold", col="Black", fontsize=15), rot=0)
top_row <- plot_grid(true.grob,l1, l2, ncol=3,rel_widths=c(1,2,2))


first_stage.grob <- textGrob("First thinning stage\n(Sampling process)", 
                             gp=gpar(fontface="bold", col="Black", fontsize=15), rot=0)
second_row <- plot_grid(first_stage.grob,l3, l4, ncol=3,rel_widths=c(1,2,2))

second_stage.grob <- textGrob("Second thinning stage\n(Detectability)", 
                              gp=gpar(fontface="bold", col="Black", fontsize=15), rot=0)
third_row <- plot_grid(second_stage.grob,l5, l6, ncol=3,rel_widths=c(1,2,2))

third_stage.grob <- textGrob("Third thinning stage\n(Misclassification)", 
                             gp=gpar(fontface="bold", col="Black", fontsize=15), rot=0)
fourth_row <- plot_grid(third_stage.grob,l7, l8, ncol=3,rel_widths=c(1,2,2))

plot_grid(top_row,second_row,ncol=1)

plot_grid(true.grob,first_stage.grob,second_stage.grob,third_stage.grob,l1,l3,l5,l7,l2,l4,l6,l8,legd_1,legd_2,legd_2,legd_3,ncol=4,rel_heights = c(.5,2,2,.5))


