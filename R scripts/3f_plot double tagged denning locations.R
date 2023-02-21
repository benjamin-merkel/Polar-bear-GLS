


proj1     <- "+proj=aeqd +lat_0=78 +lon_0=30 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km"
lons <- c(14,33)
lats <- c(74.5,81)
pts = matrix(c(lons[1]  ,seq(lons[1],lons[2],0.1), lons[2] ,rev(seq(lons[1],lons[2],0.1)),
               lats[1], rep(lats[2],length(seq(lons[1],lons[2],0.1))) ,lats[1], rep(lats[1],length(seq(lons[1],lons[2],0.1)))),
             ncol=2, byrow=F)
study.extent = st_as_sf(spPolygons((pts)))
st_crs(study.extent) <- 4326
study.extent <- st_transform(study.extent, CRS(proj1))
bbox  <- st_bbox(study.extent)
ratio <- (bbox[4]-bbox[2])/(bbox[3]-bbox[1])


# lon and lat during spring for denning bears
path3 <- path2[path2$bear_year %in% mdb$bear_year[mdb$Den %in% 1] & str_split_fixed(path2$bys, " ", 3)[,3] == "spring",]
path3 <- st_as_sf(path3[!is.na(path3$lat),], coords = c("lon", "lat"), crs = 4326)
gps3  <- gps[gps$bear_year %in% unique(path3$bear_year) & as.numeric(gps$month) %in% c(3,4),]
gps3  <- st_as_sf(gps3[!is.na(gps3$x),], coords = c("x", "y"), crs = 4326)
path3 <- path3[path3$bear_year %in% unique(gps$bear_year),]


m3 <- m2[m2$bear_year %in% mdb$bear_year[mdb$Den %in% 1] & m2$season == "spring" & !is.na(m2$gps.lon),]
m3$distance.gps <- distCosine(matrix(c(m3$gls.lon, m3$gls.lat), ncol=2), matrix(c(m3$gps.lon, m3$gps.lat), ncol=2))/1000
m3.gls <- st_as_sf(m3, coords = c("gls.lon", "gls.lat"), crs = 4326)
m3.gps <- st_as_sf(m3, coords = c("gps.lon", "gps.lat"), crs = 4326)
summary(m3)


lines <- vector(mode="list")
for(j in 1:nrow(m3)){lines[[j]] <- st_linestring(rbind(c(st_coordinates(m3.gps[j,]), st_coordinates(m3.gls[j,])))}
lines <- st_sfc(st_multilinestring(lines))
st_crs(lines) <- 4326

symbol.size <- 0.5
col.pal <- brewer.pal(4, "Set1")[c(1,4)]
mainmap <- tm_shape(study.extent) +
  tm_symbols(col=grey(1)) +
  
  tm_shape(NPland) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
  
  tm_shape(land) +
  tm_polygons(col=grey(0.7), border.col="transparent", size = 0.0001) +  #tm_lines(col="transparent") +
  
  tm_graticules(col=grey(0.7), lwd=0.5, labels.size=1) +
  
  tm_shape(lines) +
  tm_lines(col="black") +
  
  tm_shape(m3.gls) +
  tm_symbols(col = col.pal[1], shape = 22, size = 1.4, border.col = "white") +
  tm_shape(m3.gps) +
  tm_symbols(col = col.pal[2], shape = 21, size = 1.4, border.col = "white") +
  
  tm_add_legend(type="symbol", 
                col = rev(col.pal), border.col = "white", 
                labels = c("ST collar den location", "GLS den location"), 
                shape = c(21,22), size = c(2)) +
  
  tm_layout(legend.text.size = 1.1,
            legend.position = c("right","bottom"),
            legend.outside = F,
            legend.frame = F)

png(paste0("figures/seasonal centroid during spring after denning with GLS and GPS.png"), res = 800, width=20, height = 20*ratio, units="cm")
tmap_options (bg.color = grey(1), basemaps.alpha = 1)
print(mainmap)
dev.off()


