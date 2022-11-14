if(!require("tidyverse")) install.packages("tidyverse")
if(!require("osmdata")) install.packages("osmdata")
if(!require("sf")) install.packages("sf")
if(!require("dodgr")) remotes::install_git("https://git.sr.ht/~mpadge/dodgr")
if(!require("geosphere")) install.packages("geosphere")
if(!require("classInt")) install.packages("classInt")
if(!require("extrafont")) install.packages("extrafont")
if(!require("ggmap")) install.packages("ggmap")
if(!require("dplyr")) install.packages("dplyr")
library(osmar)



library(tidyverse, quietly=T) # data processing
library(osmdata, quietly=T) # load osm data
library(sf, quietly=T) # use spatial vector data
library(dodgr, quietly=T) # driving distance
library(geosphere, quietly=T) # aerial distance
library(classInt, quietly=T) # legend
library(extrafont, quietly=T) # font
library(ggmap)
library(dplyr)

#register_google(key = "AIzaSyDf3SR3rv2R5YfTiBX7dOepvNAa8P5lP6c")

#distGBIFdata <- function(gbif_data, bioclim = FALSE){
gbif_data <- datagbif.no
  message("Getting the countries from the GBIF data")
  countries <- unique(gbif_data$countryCode)

  if(bioclim == TRUE){
    message("Extracting the Bioclim data")
    alt<- raster::getData('worldclim', var='alt', res=10)
    bio <- raster::getData('worldclim', var='bio', res=10)
  }

  data_df <- list()
  for(country.tag in 1: length(countries)){
    message(paste("Extracting states and bioclim data for", countries[country.tag]))
    all_data <- gbif_data%>%
      dplyr::filter(countryCode == countries[country.tag])
    States <- raster::getData("GADM", country = countries[country.tag], level = 2)
    # USborder <- rgeos::gUnaryUnion(States, id = States$ISO)
    #boundary_points <- fortify(USborder)
    sp_df <- sp::SpatialPointsDataFrame(cbind(as.numeric(all_data$decimalLongitude),
                                              as.numeric(all_data$decimalLatitude)),
                                        all_data,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

    df_data <- raster::intersect(sp_df, States)
    data_with_admin <- df_data@data

    if(bioclim == TRUE){
      altitude <- raster::extract(alt, df_data, df=T)
      bio_data <- raster::extract(bio, df_data, df = T)
      data_df[[country.tag]] <- cbind(data_with_admin, altitude, bio_data)%>%
        data.frame()
    }else{
      data_df[[country.tag]] <- cbind(data_with_admin)%>%
        data.frame()
    }
  }

  data_df <- do.call("rbind", data_df)

  #save(gull_df, file = "gull_df.RData")

  #######################
  # Calculate distance to road
  #load filtered data

  distance_estimate <- function(i, data){

    # message("Obtaining the states from the data")
    states <- unique(data$NAME_2)
    states <- gsub("\\s*\\([^\\)]+\\)","",as.character(states))
    #states <- states[states != "Rhein-Neckar-Kreis"]
    # define Belgrade's bounding box based on ggmap's bounding box
    message(paste("Extracting maps for", states[i]))
    bg_map <- ggmap::get_map(getbb(states[i]),
                             maptype = "toner-lite",
                             source = "stamen",
                             color="bw",
                             force=T)

    bg_bbox <- attr(bg_map, 'bb')

    message(paste("map coordinates for", states[i]))
    bbox <- c(xmin=bg_bbox[,2],
              ymin= bg_bbox[,1],
              xmax= bg_bbox[,4],
              ymax=bg_bbox[,3])


    #get states's paved roads
    message(paste("Getting roads for", states[i]))
    bg_way <- osmdata::opq(bbox = bbox, timeout = 1240, memsize = 10004857600) %>%
      osmdata::add_osm_feature(
        key = 'highway') %>%
      osmdata::osmdata_sf(quiet = F)

    message(paste("Extracting the lines", states[i]))

    bg_r <- bg_way$osm_lines

    message(paste("Formatting data from", states[i], " to Polygons"))
    pk <- data %>%
      dplyr::filter(NAME_2 == states[i]) %>%
      sf::st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"),
                   crs = st_crs(bg_r))

    message(paste("Calculating distance for ", states[i]))
    dists_0 <- sf::st_distance(pk,bg_r)
    dist2road <- apply(dists_0,1,min)

    message(paste("Returning distance for ", states[i]))
    ret <- data %>%
      dplyr::filter(NAME_2 == states[i]) %>%
      dplyr::mutate(distance = dist2road)

    return(ret)

  }

  states_length <- length(unique(data_df$NAME_2))
  dataWithDistance1 <- lapply(as.list(seq(31,40)), function(x){
    ret <- distance_estimate(x, data_df)
  })

  #return(data_with_distance)
save(dataWithDistance1, file = "dataWithDistance1.RData")
#}

distanceData <- distGBIFdata(datagbif.no)

# api_list <- c('http://overpass-api.de/api/interpreter',
#               'https://lz4.overpass-api.de/api/interpreter',
#               'https://z.overpass-api.de/api/interpreter',
#               'https://overpass.kumi.systems/api/interpreter')
#
# api_to_use <- sample(1:length(api_list), 1)
#
# set_overpass_url(api_list[api_to_use])
