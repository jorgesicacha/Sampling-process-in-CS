options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(ggplot2)
library(gridExtra)

## Initial data manipulation ##
library(readr)
datagbif <- read_delim("~/Documents/GitHub/Sampling-process-in-CS/gullData/gbifData.csv",
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)

#datagbif <- readRDS("moosegbif.rds")
newcrs <- CRS("+proj=robin +datum=WGS84 +units=km")
norwaybu <-raster::getData("GADM",country="NOR",level=1)
norwaybu <- spTransform(norwaybu,newcrs)
norway <- SpatialPolygons(list(norwaybu@polygons[[6]])) ## Hedmark polygon
proj4string(norway) <- CRS(proj4string(norwaybu))
datagbif.no <- datagbif
datagbif.no <- datagbif.no[which(datagbif.no$basisOfRecord=="HUMAN_OBSERVATION"),] #Only human observation records
dates <- paste(datagbif.no$day,datagbif.no$month,datagbif.no$year,sep="/")
dates <- as.Date(dates,format = "%d/%m/%Y")
datagbif.no <- datagbif.no[which(dates > "2000-01-01"),]
coordsgbif.no <- datagbif.no[,c("decimalLongitude","decimalLatitude")]
coordsgbif.no<- distinct(coordsgbif.no)
coordsgbif.no <- coordsgbif.no[complete.cases(coordsgbif.no),]
sppointscoords.no <- SpatialPoints(coordsgbif.no,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
sppointscoords.no <- spTransform(sppointscoords.no,newcrs)
sppoints.gbif <- gIntersection(norway,sppointscoords.no)


## Making the models ##

source('spde-book-functions.R')
premesh<- inla.sp2segment(norway)
mesh <- inla.mesh.2d(boundary = premesh, loc = sppoints.gbif, offset = c(10, -0.2),
                     max.edge = c(7,50),cutoff=4,min.angle = 30)

nv <- mesh$n
meshpoints <- mesh$loc[,1:2]
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = '')
points(sppoints.gbif, col = 2, pch = 19)

#### Defining PC priors ####
spde <- inla.spde2.pcmatern(mesh,
                            prior.sigma = c(1, 0.05),
                            prior.range = c(15, 0.05))


#### Constructing the dual mesh #### Partly taken from Krainski book
dmesh <- book.mesh.dual(mesh)

#### Converting the domain polygon into a SpatialPolygons class ####
domain.polys <- norway@polygons
domainSP0 <- SpatialPolygons(domain.polys)

w1 <- pblapply(1:length(dmesh), function(i) {
  gIntersects(dmesh[i, ], domainSP0)})
w1.1 <- do.call(rbind,w1)
table(w1.1)
w1.1.1 <- which(w1.1)
length(w1.1.1)
## Let's focus on the largest polygon
areas.no <- c()
for(i in 1:length(norway@polygons[[1]]@Polygons)){
  areas.no[i]<- norway@polygons[[1]]@Polygons[[i]]@area
  print(i)
}
summary(areas.no)
which.max(areas.no)
sum(areas.no)
domain.polys1 <- norway@polygons[[1]]@Polygons[[which.max(areas.no)]]
domainSP1 <- SpatialPolygons(list(Polygons(list(domain.polys1),1)))

w2 <- pblapply(1:length(dmesh), function(i) {
  gIntersects(dmesh[i, ], domainSP1)})
w2.1 <- do.call(rbind,w2)
table(w2.1)
w2.1.1 <- which(w2.1)
length(w2.1.1)
w2.2 <- pblapply(w2.1.1, function(i) {
  gArea(gIntersection(dmesh[i, ], domainSP1))})
w2.2.1 <- do.call(rbind,w2.2) ##Needed areas

no.int <- setdiff(1:length(dmesh),w1.1.1)
df.val0 <- data.frame(poly=no.int,area=0)
int.1 <- w2.1.1
df.val1 <- data.frame(poly=int.1,area=w2.2.1)
df.val <- rbind(df.val0,df.val1)
melt.areas <- melt(df.val,id="poly")
areas.final <- cast(melt.areas,poly~variable,sum)
sum(areas.final$area) ##For validation
sum(areas.no)
w <- areas.final$area
#### Summary of these weights ####
sum(w)
table(w > 0)

##################################
par(mar = c(2, 2, 1, 1), mgp = 2:0)
plot(mesh$loc, asp = 1,  pch = 19, xlab = '', ylab = '',col=(w==0)+1)
plot(dmesh, add = TRUE)
points(sppoints.gbif, col = "blue", pch = 19,cex=0.1)
##################################

### Getting the covariates ###

ext <- extent(c(min(mesh$loc[,1])-10,max(mesh$loc[,1])+10,min(mesh$loc[,2])-10,max(mesh$loc[,2])+10))
ext2_0 <- spTransform(SpatialPoints(matrix(ext,ncol=2,byrow=F),proj4string =CRS(proj4string(norway)) ),CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ext2 <- extent(c(min(ext2_0@coords[,1]),max(ext2_0@coords[,1]),min(ext2_0@coords[,2]),max(ext2_0@coords[,2])))

#Calculating distance to road
