library(tmap)
library(sf)
library(raster)
library(sp)
library(dplyr)
library(mapview)
sf::sf_use_s2(FALSE)

mini_world <- read_sf('data/map data/ne_110m_land.shp')

globe_plot <- function(lat, lon, land.col = grey(0.5), grat.col = grey(0.8)) {
  # Define the orthographic projection
  # Choose lat_0 with -90 <= lat_0 <= 90 and lon_0 with -180 <= lon_0 <= 180
  ortho <- paste0('+proj=ortho +lat_0=',
                  lat,
                  ' +lon_0=',
                  lon,
                  ' +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs'
  )
  
  circle <- st_point(x = c(lon, lat)) %>% st_buffer(dist = 6371000) %>% st_sfc(crs = ortho)
  
  # Project this polygon in lat-lon
  circle_longlat <- circle %>% st_transform(crs = 4326)
  
  # circle_longlat cannot be used as it is
  # You must decompose it into a string with ordered longitudes
  # Then complete the polygon definition to cover the hemisphere
  if(lat != 0) {
    circle_longlat <- st_boundary(circle_longlat)
    
    circle_coords <- st_coordinates(circle_longlat)[, c(1,2)]
    circle_coords <- circle_coords[order(circle_coords[, 1]),]
    circle_coords <- circle_coords[!duplicated(circle_coords),]
    
    # Rebuild line
    circle_longlat <- st_linestring(circle_coords) %>% st_sfc(crs = 4326)
    
    if(lat > 0) {
      rectangle <- list(rbind(circle_coords,
                              c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                              c(X = 180, Y = 90),
                              c(X = -180, Y = 90),
                              c(X = -180, circle_coords[1, 'Y']),
                              circle_coords[1, c('X','Y')])) %>% 
        st_polygon() %>% st_sfc(crs = 4326)
    } else {
      rectangle <- list(rbind(circle_coords,
                              c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                              c(X = 180, Y = -90),
                              c(X = -180, Y = -90),
                              c(X = -180, circle_coords[1, 'Y']),
                              circle_coords[1, c('X','Y')])) %>% 
        st_polygon() %>% st_sfc(crs = 4326)
    }
    
    circle_longlat <- st_union(st_make_valid(circle_longlat), st_make_valid(rectangle))
  }
  
  # A small negative buffer is necessary to avoid polygons still disappearing in a few pathological cases
  # I should not change the shapes too much
  visible <- st_intersection(st_make_valid(mini_world), st_buffer(circle_longlat, -0.09)) %>%
    st_transform(crs = ortho)
  
  # DISCLAIMER: This section is the outcome of trial-and-error and I don't claim it is the best approach 
  # Resulting polygons are often broken and they need to be fixed
  # Get reason why they're broken
  broken_reason <- st_is_valid(visible, reason = TRUE)
  
  # First fix NA's by decomposing them
  # Remove them from visible for now
  na_visible <- visible[is.na(broken_reason),]
  visible <- visible[!is.na(broken_reason),]
  
  # Open and close polygons
  na_visible <- st_cast(na_visible, 'MULTILINESTRING') %>% 
    st_cast('LINESTRING', do_split=TRUE)
  na_visible <- na_visible %>% mutate(npts = npts(geometry, by_feature = TRUE))
  # Exclude polygons with less than 4 points
  na_visible <- na_visible %>%
    filter(npts >=4) %>%
    select(-npts) %>%
    st_cast('POLYGON')
  
  # Fix other broken polygons
  broken <- which(!st_is_valid(visible))
  for(land in broken) {
    result = tryCatch({
      # visible[land,] <- st_buffer(visible[land,], 0) # Sometimes useful sometimes not
      visible[land,] <- st_make_valid(visible[land,]) %>% 
        st_collection_extract()  
    }, error = function(e) {
      visible[land,] <<- st_buffer(visible[land,], 0)
    })
  }
  
  # Bind together the two tables
  visible <- rbind(visible, na_visible)
  
  
  #qtm(plot_point, projection = ortho)
  # Plot
  tm_shape(circle, projection = ortho) +
    tm_polygons(col = "white") +
    tm_graticules(labels.show = F, col = grat.col, lwd = 0.9) +
    tm_shape(visible) +
    tm_polygons(border.col = land.col, col = land.col) +
    tm_shape(st_cast(circle, "LINESTRING")) +
    tm_lines(col = land.col) +
    tm_layout(frame = F, bg.color="transparent") %>%
    return()
}

proj1       <- "+proj=aeqd +lat_0=78 +lon_0=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km"
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

# load("data/Seasonal GLS centroids 2022.RData")
# m2$bear_id <- str_split_fixed(m2$bear.gls, " ",2)[,1]
# m2$year    <- as.numeric(str_split_fixed(m2$bys, " ",3)[,2])
# cap        <- st_as_sf(m2[!duplicated(m2$bear.gls),], coords = c("cap.lon", "cap.lat"), crs = 4326)
load("data/Full GLS meta data.RData")
columns <- c("bear.gls","Capture.GLS.on.lon", "Capture.GLS.on.lat","Capture.GLS.off.lon", "Capture.GLS.off.lat")
m2 <- st_as_sf(meta2[,columns], coords = c("Capture.GLS.on.lon", "Capture.GLS.on.lat"), crs = 4326)
m3 <- st_as_sf(meta2[,columns], coords = c("Capture.GLS.off.lon", "Capture.GLS.off.lat"), crs = 4326)
colnames(m2) <- colnames(m3) <- c("bear.gls","lon","lat","geometry")
cap <- rbind(m2,m3)
cap <- cap[!duplicated(cap$lat, cap$lon),]

# Greenland capture locations
GL.cap <- data.frame(gls.id = c("Q178", "Q181"),
                   date   = as.Date(c("2016-04-10", "2016-05-14")),
                   lon    = c(-43.169, -32.551),   
                   lat    = c( 60.515,  68.520))
GL.cap <- st_as_sf(GL.cap, coords=c("lon","lat"), crs = 4326)



mainmap <- tm_shape(study.extent) +
  tm_polygons(col=grey(1), border.col="transparent") +
  tm_shape(NPland) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
  
  tm_shape(land) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
    
  tm_graticules(col=grey(0.7), lwd=0.5, labels.size=1) +
  
  tm_shape(cap) +
    tm_symbols(border.col = grey(1), col="red", shape=21, size =0.5) +
  
  tm_shape(st_as_sf(data.frame(lat = 78.2166, lon = 15.6333), coords = c("lon","lat"), crs = 4326)) + 
  tm_symbols(col='gold', size =1.2, shape = 22) + 
  
  # tm_shape(gls) +
  # tm_symbols(col="season") +
  
  # tm_scale_bar(position = c("RIGHT","BOTTOM")) +
  tm_add_legend('symbol',labels = c("Capture locations",
                                    "Longyearbyen"),
                shape      = c(21,22), 
                border.col = c(grey(1), grey(0)),
                col        = c("red", "gold")) +
  
  tm_layout(legend.text.size = 1.1,
            legend.position = c("left","bottom"),
            legend.outside = F,
            main.title.fontface = "italic", 
            main.title.position = "center")


# x1 <- 5
# x2 <- 35
# y1 <- 75
# y2 <- 81
# my.df <- data.frame(Plot="study extend", lon = c(seq(x1,x2,0.1), rep(x2,length(seq(y1,y2,0.1))), 
#                                                  rev(seq(x1,x2,0.1)), rep(x1,length(seq(y1,y2,0.1)))), 
#                     lat = c(rep(y1,length(seq(x1,x2,0.1))), seq(y1,y2,0.1),
#                             rep(y2,length(seq(x1,x2,0.1))), rev(seq(y1,y2,0.1))))
# sf <- sfheaders::sf_polygon(obj = my.df, x = "lon", y = "lat", polygon_id = "Plot")
# sf::st_crs( sf ) <- 4326
# Svalbard.area <- st_transform(sf, crs= ortho)
# 
# x1 <- -48
# x2 <- -40
# y1 <- 59
# y2 <- 63
# my.df <- data.frame(Plot="study extend", lon = c(seq(x1,x2,0.1), rep(x2,length(seq(y1,y2,0.1))), 
#                                                  rev(seq(x1,x2,0.1)), rep(x1,length(seq(y1,y2,0.1)))), 
#                     lat = c(rep(y1,length(seq(x1,x2,0.1))), seq(y1,y2,0.1),
#                             rep(y2,length(seq(x1,x2,0.1))), rev(seq(y1,y2,0.1))))
# sf <- sfheaders::sf_polygon(obj = my.df, x = "lon", y = "lat", polygon_id = "Plot")
# sf::st_crs( sf ) <- 4326
# Greenland.area <- st_transform(sf, crs= ortho)


inset <- globe_plot(60, 12) + 
  tm_shape(st_cast(study.extent, "LINESTRING")) + 
  tm_lines(col='red', lwd = 2, size =0.8) +
  
  tm_shape(GL.cap) +
  tm_symbols(border.col = grey(1), col="red", shape=21, size =0.5) 
  
  # tm_shape(st_cast(Greenland.area, "LINESTRING")) + tm_lines(col='red', lwd = 2) 


png(paste0("figures/study extent figure v3.png"), res = 800, width=20, height = 20*ratio, units="cm")
tmap_options (bg.color = grey(1), basemaps.alpha = 1)
print(mainmap)
print(inset, vp = grid::viewport(0.80, 0.24, width = 0.45, height = 0.45))
dev.off()


