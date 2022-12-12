library(tmap)
library(sf)
library(raster)
library(sp)
library(dplyr)
library(stringr)
library(mapview)
library(RColorBrewer)
sf::sf_use_s2(FALSE)

proj1     <- "+proj=aeqd +lat_0=78 +lon_0=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km"
proj.aezd <- CRS("+proj=aeqd  +lat_0=73  +lon_0=33 +units=km")

BS.ext     <- extent(-60,130,40,89.5)
Npland_ext <- extent(-518, 920, -155, 1116)

lons <- c(5,52)
lats <- c(73,83)
pts = matrix(c(lons[1]  ,seq(lons[1],lons[2],0.1), lons[2] ,rev(seq(lons[1],lons[2],0.1)),
               lats[1], rep(lats[2],length(seq(lons[1],lons[2],0.1))) ,lats[1], rep(lats[1],length(seq(lons[1],lons[2],0.1)))),
             ncol=2, byrow=F)
study.extent = st_as_sf(spPolygons((pts)))
st_crs(study.extent) <- 4326
study.extent <- st_transform(study.extent, CRS(proj1))
bbox  <- st_bbox(study.extent)
ratio <- (bbox[4]-bbox[2])/(bbox[3]-bbox[1])


Npland.area <- matrix(c(-518,920,920,-518,-518, 
                        -155,-155,1116,1116,-155),5,byrow = F)
Npland.area <- st_sfc(st_polygon(list(Npland.area)))
st_crs(Npland.area) <- proj.aezd

NPland <- st_read("data/map data/Coastline_GSHHS_NPI.shp")
NPland <- st_transform(NPland, proj.aezd)
NPland <- st_crop(NPland, Npland_ext)
NPland <- st_cast(NPland, "POLYGON")

land      <- st_read("data/map data","ne_10m_land")
land      <- st_crop(land, BS.ext)
land      <- st_transform(land, proj.aezd)
land      <- st_difference(land, Npland.area)


load("data/Full GLS location dataset.RData")
path <- st_as_sf(path2[!is.na(path2$lat),], coords = c("lon", "lat"), crs = 4326)
  
load("data/Seasonal GLS centroids.RData")
m2$bear_id <- str_split_fixed(m2$bear.gls, " ",2)[,1]
m2$year    <- as.numeric(str_split_fixed(m2$bys, " ",3)[,2])
cap        <- st_as_sf(m2[!duplicated(m2$bear.gls),], coords = c("cap.lon", "cap.lat"), crs = 4326)
gls        <- st_as_sf(m2, coords = c("gls.lon", "gls.lat"), crs = 4326)
table(gls$bear_id)

autumn_c <- readRDS("data/Autumn sea ice concentration 15% contour 2012-2021.RDS")
spring_c <- readRDS("data/Spring sea ice concentration 15% contour 2012-2021.RDS")

points <- cap.points <- lines <- vector(mode="list")

ids <- c(26205, 23881, 26088, 26236)#         23689, 23992)
id.chosen <- 23689 #26205 # 23992
for(j in 1:length(ids)){
  points[[j]] <- gls[gls$bear_id == ids[j] & gls$del.ratio > 0.5 & gls$lat.diff > -4,]
  points[[j]] <- points[[j]][order(points[[j]]$year),]
  cap.points[[j]] <- cap[cap$bear_id == ids[j],]
  # points <- path[path$bear_id == id.chosen,]
  for(i in 1:nrow(cap.points[[j]])) {
    x <- points[[j]][points[[j]]$year >= cap.points[[j]]$year[i],]
    if(i < nrow(cap.points[[j]])) x <- x[x$year < cap.points[[j]]$year[i+1],]
    pp <- rbind(cap.points[[j]][i,c("bys","bear.gls","gls.n","del.ratio")], x[,c("bys","bear.gls","gls.n","del.ratio")])
    if(i == 1) pp2 <- pp else pp2 <- rbind(pp2,pp)
  } 
  lines[[j]]  <- pp2 %>% st_combine() %>% st_cast("LINESTRING")
}

symbol.size <- 0.5
season.col <- brewer.pal(4, "Set2")[c(2,3)]
mainmap <- tm_shape(study.extent) +
  tm_polygons(col=grey(1), border.col="transparent") +
  
  tm_shape(autumn_c) +
  tm_lines(col=season.col[1], lwd = 0.5) +
  
  tm_shape(spring_c) +
  tm_lines(col=season.col[2], lwd = 0.5) +
  
  tm_shape(NPland) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
  
  tm_shape(land) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
  
  tm_graticules(col=grey(0.7), lwd=0.5, labels.size=1) +
  
  tm_shape(lines[[1]]) +
  tm_lines(col="black") +
  tm_shape(cap.points[[1]]) +
  tm_symbols(col = 'white', shape = 25, size = symbol.size) +
  tm_shape(points[[1]]) +
  tm_symbols(col = 'season', shape = 25, size = symbol.size, legend.col.show = F, palette = season.col) +
  tm_shape(cap.points[[1]][nrow(cap.points[[1]]),]) +
  tm_symbols(col = 'white', shape = 25, size = symbol.size) +
  
  tm_shape(lines[[2]]) +
  tm_lines(col="black") +
  tm_shape(cap.points[[2]]) +
  tm_symbols(col = 'white', shape = 23, size = symbol.size) +
  tm_shape(points[[2]]) +
  tm_symbols(col = 'season', shape = 23, size = symbol.size, legend.col.show = F, palette = season.col) +
  tm_shape(cap.points[[2]][nrow(cap.points[[2]]),]) +
  tm_symbols(col = 'white', shape = 23, size = symbol.size) +
  
  tm_shape(lines[[3]]) +
  tm_lines(col="black") +
  tm_shape(cap.points[[3]]) +
  tm_symbols(col = 'white', shape = 21, size = symbol.size) +
  tm_shape(points[[3]]) +
  tm_symbols(col = 'season', shape = 21, size = symbol.size, legend.col.show = F, palette = season.col) +
  tm_shape(cap.points[[3]][nrow(cap.points[[3]]),]) +
  tm_symbols(col = 'white', shape = 21, size = symbol.size) +
  
  tm_shape(lines[[4]]) +
  tm_lines(col="black") +
  tm_shape(cap.points[[4]]) +
  tm_symbols(col = 'white', shape = 24, size = symbol.size) +
  tm_shape(points[[4]]) +
  tm_symbols(col = 'season', shape = 24, size = symbol.size, legend.col.show = F, palette = season.col) +
  tm_shape(cap.points[[4]][nrow(cap.points[[4]]),]) +
  tm_symbols(col = 'white', shape = 24, size = symbol.size) +
  
  tm_add_legend(type="symbol", col = c(season.col, "white"), labels = c("autumn","spring","capture"), shape = 22, size = 1.5) +
  
  tm_layout(legend.text.size = 1.1,
            legend.position = c("left","bottom"),
            legend.outside = F,
            legend.frame = F)

png(paste0("figures/example seasonal centroid tracks with sea ice 5.png"), res = 800, width=20, height = 20*ratio, units="cm")
tmap_options (bg.color = grey(1), basemaps.alpha = 1)
print(mainmap)
dev.off()

