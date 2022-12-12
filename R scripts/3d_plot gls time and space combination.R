library(readxl)
library(GeoLight)
library(TwGeos)
library(sf)
library(SGAT)
library(MASS)
library(lubridate)
library(geosphere)
library(trip)
library(stringr)
library(geodist)
library(scales)
library(raster)
library(scales)
library(RColorBrewer)
source("R scripts/FUNCTIONS.R")

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

ct      <- c("guess","guess",rep("date",4),"guess",rep("date",4),rep("guess",11))
meta    <- read_excel("data/Metadata_Positions_21_01_21.xlsx", col_types = ct)
gpsdata <- readRDS("data/Bear collar data provided by Clement Dec 2020.rds")

land    <- read_sf("data/map data/Coastline_GSHHS_NPI.shp")
land    <- land[land$SOURCE=="NPI 1MILL",]
land    <- st_cast(land,"POLYGON")



load("data/Full twilight times dataset GLS.RData")
load("data/Full deg dataset GLS.RData")
load("data/Full GLS meta data.RData")
load("data/GLS files.RData")
load("data/Seasonal GLS centroids.RData")
load("data/Full GLS location dataset.RData")





sort(table(path2$bear_id[path2$bear_id %in% unique(gpsdata$ID.NR)]))


id.chosen <- "23980"
id.chosen <- "26088"

# load metadata
mt <- as.data.frame(meta[meta$Bear_ID == id.chosen,])

capture_dates <- unique(as.Date(c(mt$Collar_capture_on_date,mt$Collar_capture_off_date,mt$Collar_on_date,mt$Collar_off_date,mt$GLS_on_date,mt$GLS_off_date)))
capture_dates <- sort(capture_dates[!is.na(capture_dates)])

cap <- data.frame(date=capture_dates,
                  year=as.numeric(strftime(capture_dates,"%Y")),
                  month=as.numeric(strftime(capture_dates,"%m")))
cap$year2 <- cap$year
cap$year2[cap$month > 6] <- cap$year2[cap$month > 6] + 1
cap$year2[cap$year==2012] <- 2013
cap$year2[cap$year==2020] <- 2020

# load GLS data
p4 <- path2[path2$bear_id == id.chosen,]
p4$date <- as.Date(p4$datetime)
p4$season <- str_split_fixed(p4$bys," ",3)[,3]

dates <- data.frame(date= as.Date(c(as.Date(mt$GLS_on_date): as.Date(mt$GLS_off_date))))
p4 <- merge(p4,dates,all = T)
p4$s5lon <- rollapply(p4$lon, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
p4$s5lat <- rollapply(p4$lat, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
p4$year <- as.numeric(strftime(p4$date,"%Y"))
p4$month <- as.numeric(strftime(p4$date,"%m"))
p4$year2 <- p4$year
p4$year2[p4$month > 6] <- p4$year2[p4$month > 6] + 1

p5   <- p4[!is.na(p4$s5lat),]
p5$timediff.after  <- -difftime(p5$date,c(p5$date[2:nrow(p5)],NA),units = 'days')
p6 <- p5[p5$timediff.after > 7,]
p6 <- p6[!is.na(p6$s5lat),]


# load GLS centroid
m2$bear_id <- str_split_fixed(m2$bear.gls, " ",2)[,1]
m3 <- m2[m2$bear_id == id.chosen,]
m3$year <-as.numeric(str_split_fixed(m3$bys," ",3)[,2])
m3$year2<- m3$year
m3$year2[m3$season=="fall"]<- m3$year2[m3$season=="fall"]+1
m3 <- m3[m3$del.ratio > 0.5,]

# load GPS data
g4 <- gpsdata[gpsdata$ID.NR == id.chosen,]
g4$month <- as.numeric(strftime(g4$date,"%m"))
g4$year2 <- g4$year
g4$year2[g4$month > 6] <- g4$year2[g4$month > 6] + 1






file_select <- gls.ids[gls.ids$gls %in% mt$GLS_ID,]
file_select <- file_select[grepl(mt$Bear_ID[1],file_select$file),]
file_select <- file_select[!duplicated(file_select$gls),]

for(f in 1:nrow(file_select)){
  meta3 <- mt[mt$GLS_ID==file_select$gls[f],]
  meta3 <- meta3[!is.na(meta3$GLS_ID),]
  
  glsa          <- read.table(paste("data",file_select$file[f],sep="/"), skip=19, header = T)
  glsa$gls_id   <- file_select$gls[f]
  glsa$datetime <- as.POSIXct(strptime(paste(glsa[,1], glsa[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
  glsa$date     <- as.Date(glsa$datetime)
  glsa          <- glsa[as.Date(glsa$datetime) > as.Date(meta3$GLS_on_date) & 
                        as.Date(glsa$datetime) < as.Date(meta3$GLS_off_date),]
  
  if(f==1) gls <- glsa else gls <- rbind(gls,glsa)
}
gls$Light    <- (gls$light.lux.)
gls$year     <- as.numeric(strftime(gls$datetime,"%Y")) 
gls$month    <- as.numeric(strftime(gls$datetime,"%m")) 
gls$year2    <- gls$year
gls$year2[gls$month>6] <- gls$year2[gls$month>6]+1


twl3b <- twl6b[twl6b$bear.id==id.chosen & twl6b$Deleted==F,]


season.col <- brewer.pal(4, "Set2")[c(2,3)]
tw.colour <- brewer.pal(8, "Set2")[6]
Sys.setlocale("LC_TIME", "English")
png(paste0("figures/space and timeline figure ",id.chosen," e.png"),res=500,units="cm",width=15,height=30)
# all.years <- 2013:2017
all.years <- 2013:2020
# par(mfrow=c(length(all.years), 2), mar=c(0.1,0.1,0.1,0.1), oma=c(2,4,2,2))
layout(matrix(c(1,2,2, 3,4,4, 5,6,6, 7,8,8, 9,10,10, 11,12,12, 13,14,14, 15,16,16), length(all.years), 3, byrow = TRUE),
       widths=1, heights=c(1))
par(oma=c(2,0,0,0))

meta0 <- mt[mt$Bear_ID == id.chosen,]
meta0 <- meta0[!is.na(meta0$Bear_ID),]

for(yr in all.years){
  
  glsb <- gls[gls$year2==yr,]
  years <- unique(glsb$year)
  gls2 <- glsb[,c("datetime","Light")]
  colnames(gls2) <- c("Date","Light")
  
  meta3 <- mt[mt$GLS_ID==glsb$gls_id[1],]
  meta3 <- meta3[!is.na(meta3$GLS_ID),]
  
  cap2 <- cap[cap$year2==yr,]
  mt2 <- mt[as.Date(mt$Collar_off_date) == cap2$date | as.Date(mt$Collar_on_date) == cap2$date |
            as.Date(mt$GLS_on_date) == cap2$date | as.Date(mt$GLS_off_date) == cap2$date,]
  
  m4  <- m3[m3$year2==yr,]
  p44 <- p4[p4$year2==yr,]
  
  gps3 <- g4[g4$year2==yr,]
  
  ###### plot map
  par(mar=c(0,1,0,0))
  plot(land$geometry, col=grey(0.7), border=grey(0.7), ylim=c(75,81), xlim=c(3,35),lwd=0.5)
  lines(gps3$x,gps3$y,col=grey(0.1))
  lines(gps3$x[gps3$month %in% c(8:11)],gps3$y[gps3$month %in% c(8:11)],col=season.col[1])
  lines(gps3$x[gps3$month %in% c(2:4)],gps3$y[gps3$month %in% c(2:4)],col=season.col[2])
  
  # points(m4$gls.lon, m4$gls.lat)
  lines(rep(median(p44$lon[p44$season=="fall"], na.rm=T),2),c(75,81), col=season.col[1],lty=3)
  if(median(p44$lon[p44$season=="spring"], na.rm=T) > 5) lines(rep(median(p44$lon[p44$season=="spring"], na.rm=T),2),c(75,81), col=season.col[2],lty=3)
  
  m4.la <- m4.ls <- T
  if(nrow(m4[m4$season=="fall",])==0)   m4.la <- F else {if(is.na(m4$gls.lat[m4$season=="fall"]))   m4.la <- F }
  if(nrow(m4[m4$season=="spring",])==0) m4.ls <- F else {if(is.na(m4$gls.lat[m4$season=="spring"])) m4.ls <- F }
  
  if(m4.la) lines(c(5,35),rep(median(p44$lat[p44$season=="fall"], na.rm=T),2), col=season.col[1],lty=3)
  if(m4.ls) lines(c(5,35),rep(median(p44$lat[p44$season=="spring"], na.rm=T),2), col=season.col[2],lty=3)
  
  points(m4$gls.lon[m4$season=="fall"],   m4$gls.lat[m4$season=="fall"], pch=21, bg=season.col[1],cex=2)
  points(m4$gls.lon[m4$season=="spring"], m4$gls.lat[m4$season=="spring"], pch=21, bg=season.col[2],cex=2)
  
  if(length(p44$lon[p44$season=="fall" & !is.na(p44$lon)])>1){
    den.lon <- density(p44$lon[p44$season=="fall"], na.rm=T, from = 5, to =35)
    den.lon$x <- c(5,den.lon$x,35)
    den.lon$y <- c(0,den.lon$y,0)
    polygon(den.lon$x,den.lon$y/max(den.lon$y)*0.5+75,border="transparent",col=alpha(season.col[1],0.2))
  }
  if(length(p44$lon[p44$season=="spring" & !is.na(p44$lon)])>1){
    den.lon <- density(p44$lon[p44$season=="spring"], na.rm=T, from = 5, to =35)
    den.lon$x <- c(5,den.lon$x,35)
    den.lon$y <- c(0,den.lon$y,0)
    polygon(den.lon$x,den.lon$y/max(den.lon$y)*0.5+75,border="transparent",col=alpha(season.col[2],0.2))
  }
  if(m4.la & length(p44$lat[p44$season=="fall" & !is.na(p44$lat)])>1){
    den.lat <- density(p44$lat[p44$season=="fall"], na.rm=T, from = 75, to =81)
    den.lat$x <- c(81,den.lat$x,75)
    den.lat$y <- c(0,den.lat$y,0)
    polygon(den.lat$y/max(den.lat$y)*2+5,den.lat$x,border="transparent",col=alpha(season.col[1],0.2))
  }
  if(m4.ls & length(p44$lat[p44$season=="spring" & !is.na(p44$lat)])>1){
    den.lat <- density(p44$lat[p44$season=="spring"], na.rm=T, from = 75, to =81)
    den.lat$x <- c(81,den.lat$x,75)
    den.lat$y <- c(0,den.lat$y,0)
    polygon(den.lat$y/max(den.lat$y)*2+5,den.lat$x,border="transparent",col=alpha(season.col[2],0.2))
  }
  
  yr2 <- yr
  if(yr==2013)  yr2 <- 2012
  # if(strftime(meta3$GLS_on_date,"%Y") %in% yr2)  text(as.numeric(meta3$Capture.GLS.on.lon), as.numeric(meta3$Capture.GLS.on.lat),  "C", col="purple", cex=1.5)
  # if(strftime(meta3$GLS_off_date,"%Y") %in% yr2) text(as.numeric(meta3$Capture.GLS.off.lon), as.numeric(meta3$Capture.GLS.off.lat), "C", col="purple", cex=1.5)
  if(strftime(meta3$GLS_on_date,"%Y") %in% yr2)  points(as.numeric(meta3$Capture.GLS.on.lon), as.numeric(meta3$Capture.GLS.on.lat), pch=23, cex=2,bg=grey(1))
  if(strftime(meta3$GLS_off_date,"%Y") %in% yr2) points(as.numeric(meta3$Capture.GLS.off.lon), as.numeric(meta3$Capture.GLS.off.lat), pch=23, cex=2,bg=grey(1))
  
  axis(1,at=seq(5,35,10),labels = paste0(seq(5,35,10),'°N'),pos=75)
  axis(2,at=seq(75,81,2),labels = paste0(seq(75,81,2),'°E'),pos=5,las=2)
  # text(meta3$Capture.GLS.on.lon, meta3$Capture.GLS.on.lat, labels = "S",col=4,cex=2)
  # text(mt$Capture.GLS.off.lon, mt$Capture.GLS.off.lat, labels = "E",col=2,cex=2)

  # mtext(paste0(years[1],"/", years[2]), side = 2, line = 0.5, las=1)
  if(yr==2013) mtext("Space", side = 3, line = 0.5, cex=1.3)
  
  
  ##### plot timeline
  par(mar=c(0.5,4,0,0.5))
  offset <- -2 # adjusts the y-axis to put night (dark shades) in the middle
  lightImage( tagdata = gls2, offset = offset, zlim = c(0, 20), ylab = "hour [UTC]", main = "", xlim =as.POSIXct(c(paste0(yr-1,"-08-01"),paste0(yr,"-05-15"))))
  # MytsimagePoints(twl3b$Twilight, offset = offset, pch = 16, cex = 1.2, col = ifelse(twl3b$Rise, "firebrick", "cornflowerblue"))
  tsimageDeploymentLines(gls2$Date, lon = as.numeric(meta3$Capture.GLS.on.lon), lat = as.numeric(meta3$Capture.GLS.on.lat),
                         offset = offset, lwd = 3, col = adjustcolor(tw.colour, alpha.f = 0.8))
  
  years <- unique(glsb$year)
  if(length(years)==1) years <- rep(years,2)
  #mtext(paste0(years[1],"/", years[2]), side = 4, line = 0.5)
  
  # equinox
  fall.doy <- quantile(p4$doy[is.na(p4$lat) & p4$season=="fall"], c(0.05,0.95),na.rm=T)
  tm.calib <- as.POSIXct(as.Date(fall.doy, paste0(years[1],"-01-01")), tz = "UTC")
  rect(tm.calib[1],-10,tm.calib[2],40, border="transparent",density=30,col=alpha(season.col[1],0.5))
  
  spring.doy <- quantile(p4$doy[is.na(p4$lat) & p4$season=="spring"], c(0.05,0.95),na.rm=T)
  tm.calib <- as.POSIXct(as.Date(spring.doy, paste0(years[2],"-01-01")), tz = "UTC")
  rect(tm.calib[1],-10,tm.calib[2],40, border="transparent",density=30,col=alpha(season.col[2],0.5))
  # abline(v = tm.calib, lwd = 2, lty = 3, col = "orange")
  # tm.calib <- as.POSIXct(c(paste0(years,"-03-21"), paste0(years,"-09-21")), tz = "UTC")
  # abline(v = tm.calib, lwd = 2, lty = 3, col = "orange")
  
  # ticks for position estimations
  rug(side=1, x=p44$datetime[p44$season=="fall"],col=season.col[1],lwd=1.5)
  rug(side=1, x=p44$datetime[p44$season=="spring"],col=season.col[2],lwd=1.5)
  
  # capture timing
  gls.dates <- as.POSIXct(unique(c(meta0$GLS_on_date, meta0$GLS_off_date)))
  gls.dates <- gls.dates[!is.na(gls.dates)]
  gls.dates[gls.dates==min(gls.dates)] <- as.Date("2012-09-01")
  gls.dates[gls.dates==max(gls.dates)] <- as.Date("2020-04-15")
  
  if(yr==2013) arrows(as.POSIXct("2012-09-01"),12+offset,as.POSIXct("2012-08-02"),12+offset, lwd=3,length=0.1)
  if(yr==2020) arrows(as.POSIXct("2020-04-15"),12+offset,as.POSIXct("2020-05-14"),12+offset, lwd=3,length=0.1)
  
  points(gls.dates, rep(12 + offset, length(gls.dates)), pch=23, bg=grey(1), cex=4.5)
  text(gls.dates, 12 + offset, lwd=3,lty=1, "C", col=1, cex=1.9)
  
  if(yr==2013) mtext("Time", side = 3, line = 0.5, cex=1.3)
  
  
}
dev.off()







plot(land$geometry, col=grey(0.7), ylim=c(74,83), xlim=c(0,50))
points(m3$gls.lon,m3$gls.lat)
points(p4$s5lon, p4$s5lat,type="o",pch=19)
lines(c(mt$Capture.GLS.on.lon, p5$s5lon[1]), c(mt$Capture.GLS.on.lat, p5$s5lat[1]), lty=3)
lines(c(mt$Capture.GLS.off.lon, p5$s5lon[nrow(p5)]), c(mt$Capture.GLS.off.lat, p5$s5lat[nrow(p5)]), lty=3)
if(nrow(p6)>0){
  for(pp in 1:nrow(p6)){
    p7 <- p5[p5$date>=p6$date[pp],][1:2,]
    lines(p7$s5lon, p7$s5lat,lty=2)
  }
}
text(mt$Capture.GLS.on.lon, mt$Capture.GLS.on.lat, labels = "S",col=4,cex=2)
text(mt$Capture.GLS.off.lon, mt$Capture.GLS.off.lat, labels = "E",col=2,cex=2)
legend("topright",legend=c("track","gap (>1wk)","to start/end","start capture","end capture"),lty=c(1,2,3,NA,NA),pch=c(rep("",3),"S","E"),col=c(rep(1,3),4,2))




plot(p4$date, p4$lat,type="o",pch=19,col=grey(0.5),xlim=range(capture_dates))
points(g4$date, g4$y,type="l" ,col=2)
points(p4$date, p4$lat,type="o",pch=19,col=grey(0.5))
points(p4$date, p4$s5lat,type="o",pch=19)
abline(v=capture_dates)

plot(p4$date, p4$lon,type="o",pch=19,col=grey(0.5),xlim=range(capture_dates))
points(g4$date, g4$x,type="l" ,col=2)
points(p4$date, p4$lon,type="o",pch=19,col=grey(0.5))
points(p4$date, p4$s5lon,type="o",pch=19)
abline(v=capture_dates)


