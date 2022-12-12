library(readxl)
library(sf)

load("data/GLS files.RData")

gpsdata <- readRDS("data/Bear collar data provided by Clement Dec 2020.rds")
gpsdata$bear.date <- paste(gpsdata$ID.NR, gpsdata$date)
gps <- gpsdata[,c("x","y","bear.date")]
gps <- gps[!duplicated(gps$bear.date),]
colnames(gps) <- c("gps.lon","gps.lat","bear.date")

save(gps, file="data/Daily gps location.RData")

meta  <- data.frame(read_excel("data/Metadata_Positions_10_02_22.xlsx"))
meta2 <- data.frame(meta[!is.na(meta$GLS_ID),])
meta2 <- meta2[meta2$GLS_ID %in% gls.ids$gls,]
meta2$overlap <- 0
for(m in 1:nrow(meta2)){
  meta3 <- meta2[m,]
  gps2  <- gpsdata[gpsdata$ID.NR %in% meta3$Bear_ID & 
                   gpsdata$acquisition_time>=meta3$GLS_on_date & 
                   gpsdata$acquisition_time<=meta3$GLS_off_date,]
  if(nrow(gps2)>0) meta2$overlap[m] <- 1
}

on <- st_as_sf(meta2, coords=c("Capture.GLS.on.lon","Capture.GLS.on.lat"), crs=4326)
off <- st_as_sf(meta2, coords=c("Capture.GLS.off.lon","Capture.GLS.off.lat"), crs=4326)
meta2$capture.distance <- as.numeric(diag(st_distance(on,off)))/1000
meta2$bear.gls <- paste(meta2$Bear_ID, meta2$GLS_ID)
meta2$period <- as.numeric(difftime(meta2$GLS_off_date, meta2$GLS_on_date, units ="days"))
save(meta2 ,file="data/Full GLS meta data.RData")


meta3 <- data.frame(meta[!is.na(meta$Collar_ID),])
meta3 <- meta3[meta3$Bear_ID %in% meta2$Bear_ID,]
meta3$period <- as.numeric(difftime(meta3$Collar_off_date, meta3$Collar_on_date, units ="days"))
meta3 <- meta3[meta3$period > 0,]
save(meta3 ,file="data/Full ST collar meta data for bear that also have GLS meta data.RData")

summary(meta3)
length(unique(meta3$Bear_ID))
