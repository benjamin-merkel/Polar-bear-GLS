library(GeoLight)
library(SGAT)
library(stringr)
library(geosphere)
library(sf)

load("data/Full twilight times dataset GLS dusk adjusted.RData")


cap.cal <- read.table("data/solar angle calibration with capture locations.txt")


gpsdata <- readRDS("data/Bear collar data provided by Clement Dec 2020.rds")
gpsdata$bear.date <- paste(gpsdata$ID.NR, gpsdata$date)
# gps <- gpsdata[,c("x","y","bear.date")]
# gps <- gps[!duplicated(gps$bear.date),]
# colnames(gps) <- c("gps.lon","gps.lat","bear.date")


twl6c <- twl6b[twl6b$month %in% c(2:4,8:11),]
# twl6c$Twilight[twl6c$Rise==F] <- twl6c$Twilight[twl6c$Rise==F] - (9 * 60)
bys <- unique(twl6c$bys)
m2 = path2 = NULL
for(tt in 1:length(bys)){
  cat("\r", tt, "-", length(bys))
  
  t0        <- twl6c[twl6c$bys == bys[tt],]
  del.ratio <- table(t0$Deleted)[1]/sum(table(t0$Deleted))
  # t0        <- t0[t0$Deleted==F & t0$month %in% c(2:4,8:11),]
  # t1 <- export2GeoLight(t0)
  
  if(nrow(t0)>4){
    path <- thresholdLocation (t0$Twilight, t0$Rise, zenith = median(cap.cal$a1), tol=0.18)
    path <- data.frame(lon=path[[2]][,1], lat=path[[2]][,2],datetime=path[[1]])
    path$doy <- as.numeric(strftime(path$datetime,"%j"))
    path$bear.gls <- t0$bear.gls[1]
    path$lat[path$lat <74] <- 100
    path$lon[path$lat == 100] <- NA
    path$gls.n    = nrow(path)
    path$del.ratio= del.ratio
    path$bys = bys[tt]
    path <- path[!is.na(path$lon),]
    path <- path[!duplicated(paste(path[,1],path[,2],path[,3],path[,4])),]
    path$date <- as.Date(path$datetime)
    
    if(nrow(path)>0){
      gps0=NULL
      gps0 <- gpsdata[gpsdata$ID.NR == t0$bear.id[1] & gpsdata$date >= min(t0$date) & gpsdata$date <= max(t0$date),]
      gps0 <- gps0[!is.na(gps0$x) & !is.na(gps0$y),]
      if(nrow(gps0)>0){
        
        gps0 <- gps0[,c("x", "y", "date", "acquisition_time")]
        gpsx <- aggregate(x ~ date, data = gps0, mean)
        gpsy <- aggregate(y ~ date, data = gps0, mean)
        gps1 <- data.frame(date = gpsx$date,
                           gps.lon = gpsx$x,
                           gps.lat = gpsy$y)
        
        path <- merge(path,gps1,all.x=T,by="date")
      } else {
        path$gps.lon <- NA
        path$gps.lat <- NA
      }
      m1 <- data.frame(bys = bys[tt],
                       bear.gls = t0$bear.gls[1],
                       gls.lon  = median(path$lon, na.rm=T),
                       gls.lat  = median(path$lat, na.rm=T),
                       gls.n    = nrow(path),
                       del.ratio= del.ratio,
                       cap.lon  = t0$Capture.GLS.on.lon[1],
                       cap.lat  = t0$Capture.GLS.on.lat[1],
                       gps.lon  = NA,
                       gps.lat  = NA,
                       months = paste(unique(t0$month), collapse = " - "),
                       id=tt)
      if(nrow(gps0)>0){
        m1$gps.lon <- median(gps0$x)
        m1$gps.lat <- median(gps0$y)
      }
      
      if(is.null(m2)) m2 <- m1 else m2 <- rbind(m2,m1)
      if(is.null(path2)) path2 <- path else path2 <- rbind(path2,path)
    }
  }
}
path2$bear_id <- str_split_fixed(path2$bear.gls, " ", 2)[,1]
m2$lon.diff <- m2$gls.lon - m2$cap.lon
m2$lat.diff <- m2$gls.lat - m2$cap.lat
m2$lon.diff.gps <- m2$gls.lon - m2$gps.lon
m2$lat.diff.gps <- m2$gls.lat - m2$gps.lat
m2 <- m2[!is.na(m2$gls.lat),]
m2$distance <- distCosine(matrix(c(m2$gls.lon, m2$gls.lat), ncol=2), matrix(c(m2$cap.lon, m2$cap.lat), ncol=2))/1000
m2$season <- str_split_fixed(m2$bys," ",3)[,3]

save(m2 ,file="data/Seasonal GLS centroids.RData")
save(path2 ,file="data/Full GLS location dataset.RData")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



load("data/Seasonal GLS centroids.RData")
load("data/Full GLS location dataset.RData")

land    <- read_sf("data/map data/Coastline_GSHHS_NPI.shp")
land    <- land[land$SOURCE=="NPI 1MILL",]
land    <- st_cast(land,"POLYGON")


m2$year <- str_split_fixed(m2$bys, " ", 3)[,2]

pathg <- path2[!is.na(path2$gps.lat) & !is.na(path2$gps.lon) & !is.na(path2$lat),]
pathg$distance <- distCosine(matrix(c(pathg$lon, pathg$lat), ncol=2), matrix(c(pathg$gps.lon, pathg$gps.lat), ncol=2))/1000
pathg$doy <- as.numeric(strftime(pathg$date, "%j"))
pathg$doy <- factor(pathg$doy,levels = 1:365)
pathg$month <- as.numeric(strftime(pathg$date, "%m"))
boxplot(pathg$distance~pathg$doy)



png("figures/Daily gps vs gls distance.png",res=500,units="cm",width=20,height=20)
opar<-par(mar=c(4.1,4.1,0.1,2.1))
boxplot(pathg$distance~ pathg$bear_id, horizontal = T,las=1, ylab="",xlab="distance between GLS and GPS locations [km]",col="white")
points(tapply(pathg$distance, pathg$bear_id, mean), 1:length(unique(pathg$bear_id)), pch=19, col=2)
axis(4, at= 1:length(unique(pathg$bear_id)), labels = table(pathg$bear_id), las = 1)
abline(v=mean(tapply(pathg$distance, pathg$bear_id, mean)),lty=3, col=2)
abline(v=mean(tapply(pathg$distance, pathg$bear_id, median)),lty=3, col=1)
par(opar)
dev.off()

# lm1 <- lm(distance ~ gls.n*del.ratio, data = m2)
# summary(lm1)
# opar <- par(mfrow=c(2,2));plot(lm1);par(opar)

m3 <- m2[m2$gls.n >9 & m2$del.ratio<0.6,]
boxplot(m3$del.ratio ~m3$year)
table(m3$year)


m3 <- m2[m2$gls.n >9 & m2$del.ratio>0.6,]
summary(m3)


png("figures/Seasonal Capture distance vs gls n.png",res=500,units="cm",width=20,height=20)
opar<-par(mar=c(4.1,4.1,0.1,0.1))
plot(m2$gls.n, m2$distance, bg = grey(1-m2$del.ratio),pch=21, col=grey(1),cex=1.1,las=1, lwd = 0.5, ylim=c(0,800), xlim = c(0,250),
     ylab="distance between centroid and capture location [km]", xlab="# of estimated GLS locations per season")
legend("topright", legend = paste(seq(0.9,0,-0.2)*100, "%"), pch = 19, pt.cex=1, col =grey(seq(0.9,0,-0.2)), title = "Ratio deleted", box.col = "transparent")
box()
par(opar)
dev.off()


m2g <- m2[!is.na(m2$gps.lat),]
m2g$distance.gps <- distCosine(matrix(c(m2g$gls.lon, m2g$gls.lat), ncol=2), matrix(c(m2g$gps.lon, m2g$gps.lat), ncol=2))/1000
summary(m2g)
summary(m2g[m2g$del.ratio>0.6,])

png("figures/Seasonal GPS distance vs gls n.png",res=500,units="cm",width=20,height=20)
opar<-par(mar=c(4.1,4.1,0.1,0.1))
plot(m2g$gls.n, m2g$distance.gps, bg = grey(1-m2g$del.ratio),pch=21, col=grey(1),cex=1.1,las=1, lwd = 0.5, ylim=c(0,800), xlim = c(0,250),
     ylab="distance between centroid and capture location [km]", xlab="# of estimated GLS locations per season")
legend("topright", legend = paste(seq(0.9,0,-0.2)*100, "%"), pch = 19, pt.cex=1, col =grey(seq(0.9,0,-0.2)), title = "Ratio deleted", box.col = "transparent")
box()
par(opar)
dev.off()

m3 <- m2[m2$gls.n >9 & m2$del.ratio>0.6,]
png("figures/ALL seasonal GLS median vs capture locations.png",res=500,units="cm",width=20,height=25)
opar<-par(mar=c(0,0,0,0))
plot(land$geometry, col=grey(0.7), border=grey(0.7), ylim=c(74,83)) #ylim=c(75,82.5)
for(tt in unique(m3$bys)) lines(c(m3$gls.lon[m3$bys==tt], m3$cap.lon[m3$bys==tt]), c(m3$gls.lat[m3$bys==tt], m3$cap.lat[m3$bys==tt]), col=grey(0.3))
points(m3$cap.lon[m3$season=="fall"],   m3$cap.lat[m3$season=="fall"], pch=23, col=2, bg=grey(1))
points(m3$cap.lon[m3$season=="spring"], m3$cap.lat[m3$season=="spring"], pch=21, col=4, bg=grey(1))

points(m3$gls.lon[m3$season=="fall"],   m3$gls.lat[m3$season=="fall"], pch=23, bg=2)
points(m3$gls.lon[m3$season=="spring"], m3$gls.lat[m3$season=="spring"], pch=21, bg=4)
legend("topright",legend=c("spring centroid","fall centroid","GLS","Initial capture"),pch=c(21,23,22,22), pt.bg= c(4,2,1,grey(1)))
mtext(paste("Solar angle:",round(median(cap.cal$a1),3)),3,line = -1)
par(opar)
dev.off()


m3 <- m2[m2$gls.n >9 & !is.na(m2$gps.lat),]
m3$distance.gps <- distCosine(matrix(c(m3$gls.lon, m3$gls.lat), ncol=2), matrix(c(m3$gps.lon, m3$gps.lat), ncol=2))/1000
summary(m3)

plot(m3$gls.n, m3$distance.gps)
plot(m3$del.ratio, m3$distance.gps)

png("figures/ALL seasonal GLS vs GPS median locations.png",res=500,units="cm",width=25,height=25)
opar<-par(mar=c(0,0,0,0))
plot(land$geometry, col=grey(0.7), border=grey(0.7), ylim=c(74,83)) #ylim=c(75,82.5)
for(tt in unique(m3$bys)) lines(c(m3$gls.lon[m3$bys==tt], m3$gps.lon[m3$bys==tt]), c(m3$gls.lat[m3$bys==tt], m3$gps.lat[m3$bys==tt]), col=grey(0.3))
points(m3$gps.lon[m3$season=="fall"],   m3$gps.lat[m3$season=="fall"], pch=23, bg=grey(1), col=2)
points(m3$gps.lon[m3$season=="spring"], m3$gps.lat[m3$season=="spring"], pch=21, bg=grey(1),col=4)

points(m3$gls.lon[m3$season=="fall"],   m3$gls.lat[m3$season=="fall"], pch=23, bg=2,cex=2)
points(m3$gls.lon[m3$season=="spring"], m3$gls.lat[m3$season=="spring"], pch=21, bg=4,cex=2)
# points(m3$gls.lon[m3$season=="fall"],   m3$gls.lat[m3$season=="fall"], pch=23, bg=2, cex=m3$gls.n[m3$season=="fall"]/max(m3$gls.n)*2)
# points(m3$gls.lon[m3$season=="spring"], m3$gls.lat[m3$season=="spring"], pch=21, bg=4, cex=m3$gls.n[m3$season=="spring"]/max(m3$gls.n)*2)
legend("topright",legend=c("spring centroid","fall centroid","GLS","GPS"),pch=c(21,23,22,22), pt.bg= c(4,2,1,grey(1)))
mtext(paste("Solar angle:",round(median(cap.cal$a1),3)),3,line = -1)
par(opar)
dev.off()

path3              <- path2[!is.na(path2$gps.lat) & !is.na(path2$lat),]
path3$distance.gps <- distCosine(matrix(c(path3$lon, path3$lat), ncol=2), matrix(c(path3$gps.lon, path3$gps.lat), ncol=2))/1000
path3$season       <- "autumn"
path3$season[path3$doy < 150] <- "spring"

boxplot(path3$distance.gps ~ path3$season)


path3 <- path2[path2$gls.n >9 & path2$del.ratio>0.6,]


plot(path2$doy,path2$lat,col=grey(0.7))
points(path3$doy,path3$lat)

plot(path2$doy,path2$lon,col=grey(0.7))
points(path3$doy,path3$lon)


tt="26018 B128"
tt="26205 J113"
tt="23811 Q162"

daily.color <- wheel("darkblue", num = 12)
daily.color <- wheel("tomato", num = 12, init.angle = -100)


path3 <- path3[order(path3$datetime),]
path3$month <- as.numeric(strftime(path3$datetime,"%m"))
for(tt in unique(path3$bear.gls)) {
  mt <- meta2[meta2$bear.gls==tt,]
  mt$Capture.GLS.on.lon  <- as.numeric(mt$Capture.GLS.on.lon)
  mt$Capture.GLS.on.lat  <- as.numeric(mt$Capture.GLS.on.lat)
  mt$Capture.GLS.off.lon <- as.numeric(mt$Capture.GLS.off.lon)
  mt$Capture.GLS.off.lat <- as.numeric(mt$Capture.GLS.off.lat)
  
  p4       <- path2[path2$bear.gls==tt,]
  p4$date  <- as.Date(p4$datetime)
  p4$month <- as.numeric(strftime(p4$date, "%m"))
  
  dates <- data.frame(date= as.Date(c(as.Date(mt$GLS_on_date): as.Date(mt$GLS_off_date))))
  p4 <- merge(p4,dates,all = T)
  
  # p4$timediff.after  <- -difftime(p4$datetime,c(p4$datetime[2:nrow(p4)],NA),units = 'days')
  # p4$distance.after  <- distCosine(matrix(c(p4$lon, p4$lat), ncol=2), matrix(c(c(p4$lon[2:nrow(p4)],NA), c(p4$lat[2:nrow(p4)],NA)), ncol=2))/1000
  p4$s5lon <- rollapply(p4$lon, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
  p4$s5lat <- rollapply(p4$lat, width = 5, FUN = mean, na.rm = TRUE, align = "center", partial = T)
  
  # p5       <- p4[!is.na(p4$s5lat),]
  p5       <- p4[!is.na(p4$lat),]
  p5$timediff.after  <- -difftime(p5$date,c(p5$date[2:nrow(p5)],NA),units = 'days')
  # p5$timediff.before <- p5$datetime - c(p5$datetime[1:(nrow(p5)-1)], NA)
  p6 <- p5[p5$timediff.after > 7,]
  # p6 <- p6[!is.na(p6$s5lat),]
  # 
  # plot(p4$date, p4$lat,type="o",pch=19,col=grey(0.5))
  # points(p4$date, p4$s5lat,type="o",pch=19)
  # 
  # plot(p4$date, p4$lon,type="o",pch=19,col=grey(0.5))
  # points(p4$date, p4$s5lon,type="o",pch=19)
  
  
  png(paste0("figures/ind maps/2022_",tt,".png"),res=500,units="cm",width=20,height=20)
  opar<-par(mar=c(0,0,0,0))
  plot(land$geometry, col=grey(0.7), ylim=c(74,83), xlim=c(-10,70))
  # points(p4$s5lon, p4$s5lat,type="o",pch=19)
  # lines(c(mt$Capture.GLS.on.lon, p5$s5lon[1]), c(mt$Capture.GLS.on.lat, p5$s5lat[1]), lty=3)
  # lines(c(mt$Capture.GLS.off.lon, p5$s5lon[nrow(p5)]), c(mt$Capture.GLS.off.lat, p5$s5lat[nrow(p5)]), lty=3)
  # if(nrow(p6)>0){
  #   for(pp in 1:nrow(p6)){
  #     p7 <- p5[p5$date>=p6$date[pp],][1:2,]
  #     lines(p7$s5lon, p7$s5lat,lty=2)
  #   }
  # }
  points(p4$lon, p4$lat,type="o",pch=19)
  points(p4$lon, p4$lat,type="p",pch=21, bg=daily.color[p4$month], cex=1.1)
  lines(c(mt$Capture.GLS.on.lon, p5$lon[1]), c(mt$Capture.GLS.on.lat, p5$lat[1]), lty=3)
  lines(c(mt$Capture.GLS.off.lon, p5$lon[nrow(p5)]), c(mt$Capture.GLS.off.lat, p5$lat[nrow(p5)]), lty=3)
  if(nrow(p6)>0){
    for(pp in 1:nrow(p6)){
      p7 <- p5[p5$date>=p6$date[pp],][1:2,]
      lines(p7$lon, p7$lat,lty=2)
    }
  }
  text(mt$Capture.GLS.on.lon, mt$Capture.GLS.on.lat, labels = "S",col=4,cex=2)
  text(mt$Capture.GLS.off.lon, mt$Capture.GLS.off.lat, labels = "E",col=2,cex=2)
  legend("topright",legend=c("track","gap (>1wk)","to start/end","start capture","end capture"),lty=c(1,2,3,NA,NA),pch=c(rep("",3),"S","E"),col=c(rep(1,3),4,2))
  legend("right", legend=1:12, pch = 21, pt.cex=2, pt.bg =daily.color, title = "month")
  par(opar)
  dev.off()
}


