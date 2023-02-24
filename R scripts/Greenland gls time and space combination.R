library(sf)
library(readxl)
library(stringr)
library(GeoLight)
library(RColorBrewer)
library(SGAT)
library(geosphere)
library(raster)
library(scales)
library(TwGeos)


land2 <- read_sf('E://map data/ne_50m_land.shp')

cap.cal <- read.table("data/solar angle calibration with capture locations.txt")


# capture/recapture for Q181 is infered from ARGOS locations and GLS longitude data
meta <- data.frame(gls.id = c("Q178", "Q178", "Q181", "Q181"),
                   type   = c("capture", "recapture","capture","recapture"),
                   date   = as.Date(c("2016-04-10", "2018-05-01","2016-05-14","2017-03-20")),
                   lon    = c(-43.169, -43.329, -32.551, -32.074),   
                   lat    = c(60.515, 60.43, 68.520, 68.372))
# Q 178 shot May 1 2018 in Kuugaarmiut which is 43´19´39.6 W 60´25´45.3 N



rerun.GLS <- F
if(rerun.GLS){
  gls_files_GL         <- list.files(path="data/Migratech_GLS_data/Greenland", pattern = "driftadj.lux$", recursive = F)
  for(i in 1:length(gls_files_GL)){
    # gls <- as.character(read.table(paste("data",gls_files_notdriftadj[i],sep="/"),skip=2,nrow=1,header = F)[1,3])
    gls <- as.character(read.table(paste("data/Migratech_GLS_data/Greenland",gls_files_GL[i],sep="/"),skip=2,nrow=1,header = F)[1,3])
    gls <- data.frame(gls=gls, driftadj=T)
    gls$driftadj[!grepl("driftadj", gls_files_GL[i])]<-F
    gls$file <- gls_files_GL[i]
    
    if(i==1) gls.ids <- gls else gls.ids <- rbind(gls.ids, gls)
  }
  gls.ids <- gls.ids[order(gls.ids$driftadj, decreasing = T),]
  save(gls.ids ,file="data/GLS files GL.RData")
} else {
  load("data/GLS files GL.RData")
}


rerun.twl <- F
if(rerun.twl){
  twl6 <- NULL
  for(m in 1:nrow(gls.ids)){
    
    cat("\r",m)
    
    
    file_select <- gls.ids$file[m]
    gls          <- read.table(paste("data/Migratech_GLS_data/Greenland",file_select[1],sep="/"), skip=19, header = T)
    gls$datetime <- as.POSIXct(strptime(paste(gls[,1], gls[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
    gls$date     <- as.Date(gls$datetime)
    # gls          <- gls[as.Date(gls$datetime) > as.Date(meta3$GLS_on_date) & 
    #                       as.Date(gls$datetime) < as.Date(meta3$GLS_off_date),]
    
    twl <- twilightCalc(datetime = gls$datetime, light = log(gls$light.lux.), ask = F, maxLight = 10)
    twl$Rise    <- T
    twl$Rise[twl$type==2] <- F
    twl3b <- NULL
    
    if(nrow(twl)>6){
      
      # split into chunks divided by longer gaps of two weeks or more
      twl$diff  <- as.numeric(difftime(twl$tSecond, twl$tFirst, units = "hour"))
      twl$chunk <- 1
      for(i in 1:(nrow(twl)-1)) if(twl$diff[i]> 24 * 14) twl$chunk[(i+1):nrow(twl)] <- twl$chunk[(i+1):nrow(twl)] + 1
      twl <- twl[twl$chunk %in% names(table(twl$chunk))[table(twl$chunk)>14],]
      twl <- twl[twl$diff < 24 *14,]
      
      twl.switch <- twl
      twl.switch$tFirst <- twl.switch$tSecond
      twl.switch$Rise[twl.switch$type==2] <- T
      twl.switch$Rise[twl.switch$type==1] <- F
      twl2 <- rbind(twl, twl.switch)
      twl2 <- twl2[!duplicated(twl2[,1]),]
      twl2 <- twl2[,c(1,4,6)]
      colnames(twl2) <- c("Twilight","Rise","chunk")
      
      # remove nights or days less than 3 hours long
      twl2$diff  <- as.numeric(difftime(c(twl2$Twilight[2:nrow(twl2)],NA), twl2$Twilight, units = "hour"))
      twl2 <- twl2[!(twl2$diff<3 & twl2$Rise==F),]
      twl2 <- twl2[!is.na(twl2$Twilight),]
      
      offset <- -6 # adjusts the y-axis to put night (dark shades) in the middle
      par(mar=c(4,4,1,1))
      for(chunks in unique(twl2$chunk)){
        
        
        try(twl3b <- MytwilightEdit(twilights = twl2[twl2$chunk==chunks,1:2],
                                    offset = offset,
                                    window = 4,           # two days before and two days after
                                    outlier.mins = 30,    # difference in mins
                                    stationary.mins = 15,  # are the other surrounding twilights within x mins of one another
                                    plot = T), silent = T)
        if(is.null(twl3b)) twl3b <- cbind(twl2[twl2$chunk==chunks,1:2], Deleted = F, Edited = F, Twilight0 = twl2[twl2$chunk==chunks,1])
        twl3b$chunk <- chunks
        if(chunks == unique(twl2$chunk)[1]) twl3 <- twl3b else twl3 <- rbind(twl3, twl3b)
      }
      
      
      twl3$year  <- as.numeric(strftime(twl3$Twilight,'%Y'))
      twl3$month <- as.numeric(strftime(twl3$Twilight,'%m'))
      
      twl4 <- twl3#[twl3$Deleted==F,]
      
      twl4$Twilight <- with_tz(twl4$Twilight, "UTC")
      twl4$date     <- as.Date(twl4$Twilight)
      
      twl4$gls.id  <- gls.ids$gls[m]
      twl4$doy     <- as.numeric(strftime(twl4$Twilight, "%j"))
      
      if(is.null(twl6)) twl6 <- twl4 else twl6 <- rbind(twl6, twl4)
    }
  }
  twl6 <- twl6[!is.na(twl6$year),]
  twl6$season <- "spring"
  twl6$season[twl6$month > 4] <- "summer"
  twl6$season[twl6$month > 7] <- "fall"
  twl6$season[twl6$month  %in% c(12,1)] <- "winter"
  twl6$year2    <- twl6$year
  twl6$year2[twl6$month > 7] <- twl6$year2[twl6$month>7]+1
  twl6$gys <- paste(twl6$gls.id, twl6$year2, twl6$season)
  save(twl6 ,file="data/Full twilight times dataset GLS GL.RData")
} else {
  load("data/Full twilight times dataset GLS GL.RData")
}




raw <- read_excel("data/Greenland Argos/Raw_24054.xlsx") 
raw$difftime <- (as.numeric(raw$Date) - c(NA,head(raw$Date, -1))) / 3600 / 24
summary(raw$difftime)

raw <- read_excel("data/Greenland Argos/Raw_24042.xlsx") 
raw$difftime <- (as.numeric(raw$Date) - c(NA,head(raw$Date, -1))) / 3600 / 24
summary(raw$difftime)




ids=unique(meta$gls.id)
for(i in 1:length(ids)){
tt <- twl6[twl6$date > meta$date[meta$gls.id == ids[i] & meta$type=="capture"] & 
           twl6$date < meta$date[meta$gls.id == ids[i] & meta$type=="recapture"] & 
           twl6$gls.id == ids[i],]
if(i == 1) twl6b <- tt else twl6b <- rbind(twl6b, tt)
}
twl6 <- twl6b

tapply(twl6$Twilight, twl6$gls.id, summary)

# Argos for Q178
argos <- read_excel("data/Greenland Argos/OPD_24054.xlsx") 
argos$year     <- as.numeric(strftime(argos$Date,"%Y")) 
argos$month    <- as.numeric(strftime(argos$Date,"%m")) 
argos$year2    <- argos$year
argos$year2[argos$month > 7] <- argos$year2[argos$month>7]+1
argos$season <- "spring"
argos$season[argos$month > 4] <- "summer"
argos$season[argos$month > 7] <- "fall"
argos$season[argos$month  %in% c(12,1)] <- "winter"
argos$gls.id <- "Q178"
# argos <- st_as_sf(argos, coords = c("Longitude","Latitude"),crs=st_crs())
argos$gys <- paste(argos$gls.id,argos$year2,argos$season) 
argos$difftime <- (as.numeric(argos$Date) - c(NA,head(argos$Date, -1))) / 3600 / 24

# Argos for Q181
argos2 <- read_excel("data/Greenland Argos/OPD_24042.xlsx") 
argos2$year     <- as.numeric(strftime(argos2$Date,"%Y")) 
argos2$month    <- as.numeric(strftime(argos2$Date,"%m")) 
argos2$year2    <- argos2$year
argos2$year2[argos2$month > 7] <- argos2$year2[argos2$month>7]+1
argos2$season <- "spring"
argos2$season[argos2$month > 4] <- "summer"
argos2$season[argos2$month > 7] <- "fall"
argos2$season[argos2$month  %in% c(12,1)] <- "winter"
argos2$gls.id <- "Q181"
# argos2 <- st_as_sf(argos2, coords = c("Longitude","Latitude"),crs=st_crs())
argos2$gys <- paste(argos2$gls.id,argos2$year2,argos2$season) 
argos2$difftime <- (as.numeric(argos2$Date) - c(NA,head(argos2$Date, -1))) / 3600 / 24

arg2 <- rbind(argos[,names(argos2)],argos2)




gys <- unique(twl6$gys)
m2 = path2 = NULL
for(tt in 1:length(gys)){
  t0 <- twl6[twl6$gys == gys[tt],]
  t0$Deleted <- factor(t0$Deleted, levels = c("FALSE",'TRUE'))
  del.ratio  <- table(t0$Deleted)[1]/sum(table(t0$Deleted))
  t0 <- t0[t0$Deleted==F,]
  # t1 <- export2GeoLight(t0)
  a0 <- arg2[arg2$gys == gys[tt],]
  
  if(nrow(t0)>4){
    path <- thresholdLocation (t0$Twilight, t0$Rise, zenith = median(cap.cal$a1), tol=0.18) # 93.833
    path <- data.frame(lon=path[[2]][,1], lat=path[[2]][,2],datetime=path[[1]])
    
    if(nrow(path) > 0){
      path$doy <- as.numeric(strftime(path$datetime,"%j"))
      path$date <- as.Date(path$datetime)
      path$bear.gls <- t0$gls.id[1]
      path$lat[path$lat < 40] <- 100
      # path$lon[path$lon < -150] <- NA
      # path$lon[path$lon > 5] <- NA
      path$lon[path$lat == 100] <- NA
      path$gls.n    = nrow(path)
      path$del.ratio= del.ratio
      # path <- path[!is.na(path$lon),]
      path <- path[!duplicated(paste(path[,1],path[,2],path[,3],path[,4])),]
      
      gps0=NULL
      gps0 <- arg2[arg2$gls.id == t0$gls.id[1] & as.Date(arg2$Date) >= min(t0$date) & as.Date(arg2$Date) <= max(t0$date),]
      gps0 <- gps0[!is.na(gps0$Longitude) & !is.na(gps0$Latitude),]
      gps0$x <- gps0$Longitude
      gps0$y <- gps0$Latitude
      if(nrow(gps0)>0){
        
        gps0 <- gps0[,c("x", "y", "Date", "Date2")]
        gpsx <- aggregate(x ~ Date2, data = gps0, mean)
        gpsy <- aggregate(y ~ Date2, data = gps0, mean)
        gps1 <- data.frame(date = as.Date(gpsx$Date2),
                           gps.lon = gpsx$x,
                           gps.lat = gpsy$y)
        
        path <- merge(path, gps1, all.x=T, by="date")
      } else {
        path$gps.lon <- NA
        path$gps.lat <- NA
      }
      
      m1 <- data.frame(gys = gys[tt],
                       gls.id   = t0$gls.id[1],
                       gls.lon  = median(path$lon, na.rm=T),
                       gls.lat  = median(path$lat, na.rm=T),
                       gls.n    = nrow(path),
                       argos.lon  = median(a0$Longitude, na.rm=T),
                       argos.lat  = median(a0$Latitude, na.rm=T),
                       argos.n    = nrow(a0),
                       del.ratio= del.ratio,
                       months = paste(unique(t0$month), collapse = " - "),
                       id=tt)
      
      if(is.null(m2)) m2 <- m1 else m2 <- rbind(m2,m1)
      if(is.null(path2)) path2 <- path else path2 <- rbind(path2,path)
    }
  }
}
m2$season   <- str_split_fixed(m2$gys," ",3)[,3]
# m2 <- m2[!is.na(m2$gls.lat) & m2$season!="winter" & !is.na(m2$argos.lat),]
m2$lon.diff <- m2$gls.lon - m2$argos.lon
m2$lat.diff <- m2$gls.lat - m2$argos.lat
m2$distance <- NA
for(i in 1:nrow(m2)){
  if(!is.na(m2$lat.diff[i])){
    m2$distance[i] <- distCosine(matrix(c(m2$gls.lon[i], m2$gls.lat[i]), ncol=2), matrix(c(m2$argos.lon[i], m2$argos.lat[i]), ncol=2))/1000
  }
}
path2$distance <- NA
for(i in 1:nrow(path2)){
  if(!is.na(path2$gps.lat[i]) & !is.na(path2$lat[i])){
    path2$distance[i] <- distCosine(matrix(c(path2$lon[i], path2$lat[i]), ncol=2), matrix(c(path2$gps.lon[i], path2$gps.lat[i]), ncol=2))/1000
  }
}


# % location coverage throughout the year
path_year <- path2[!is.na(path2$lat) & !duplicated(paste(path2$bear.gls, as.Date(path2$datetime))),]
path_year$year <- as.numeric(strftime(path_year$datetime, "%Y"))
summary(as.numeric(table(paste(path_year$bear.gls,path_year$year))))/720
table(paste(path_year$bear.gls))/c(751,310)
#% using just longitudes
path_year <- path2[!duplicated(paste(path2$bear.gls, as.Date(path2$datetime))),]
path_year$year <- as.numeric(strftime(path_year$datetime, "%Y"))
summary(as.numeric(table(paste(path_year$bear.gls,path_year$year))))/720
table(paste(path_year$bear.gls))/c(751,310)



m2$distance[is.na(m2$lat.diff)] <- NA
m2$year     <- str_split_fixed(m2$gys," ",3)[,2]
summary(m2)
summary(path2)

nrow(path2[!is.na(path2$distance),])
nrow(m2[!is.na(m2$distance),])

m2[,c(2:8,16,12:15)]




plot(st_as_sf(arg2, coords = c("Longitude","Latitude"),crs=4326)$geometry, col=4)
plot(land2$geometry,add=T)
points(m2$gls.lon, m2$gls.lat,pch=19,col=2)




GLS <- "Q178" 
GLS <- "Q181" 
for(GLS in c("Q181", "Q178")){
  
  capture <- meta[meta$gls.id == GLS,]
  
  file_select <- gls.ids[gls.ids$gls %in% GLS,]
  for(f in 1:nrow(file_select)){
    
    glsa          <- read.table(paste("data/Migratech_GLS_data/Greenland",file_select$file[f],sep="/"), skip=19, header = T)
    glsa$gls_id   <- file_select$gls[f]
    glsa$datetime <- as.POSIXct(strptime(paste(glsa[,1], glsa[,2]), "%d/%m/%Y %H:%M:%S"), tz = "UTC")
    glsa$date     <- as.Date(glsa$datetime)
    
    if(f==1) gls <- glsa else gls <- rbind(gls,glsa)
  }
  gls$Light    <- (gls$light.lux.)
  gls$year     <- as.numeric(strftime(gls$datetime,"%Y")) 
  gls$month    <- as.numeric(strftime(gls$datetime,"%m")) 
  gls$year2    <- gls$year
  gls$year2[gls$month>7] <- gls$year2[gls$month>7]+1
  gls <- gls[gls$date > capture$date[capture$type == "capture"],]
  gls <- gls[gls$date < capture$date[capture$type == "recapture"],]
  
  p4 <- path2[path2$bear.gls == GLS,]
  p4$year     <- as.numeric(strftime(p4$datetime,"%Y")) 
  p4$month    <- as.numeric(strftime(p4$datetime,"%m")) 
  p4$year2    <- p4$year
  p4$year2[p4$month>7] <- p4$year2[p4$month>7]+1
  p4$season <- "spring"
  p4$season[p4$month > 4] <- "summer"
  p4$season[p4$month > 7] <- "fall"
  p4$season[p4$month  %in% c(12,1)] <- "winter"
  
  
  all.years <- unique(gls$year)
  if(GLS == "Q181") all.years <- 2017
  
  
  
  Sys.setlocale("LC_TIME", "English")
  png(paste0("figures/timeline Greenland ",GLS,".png"),res=500,units="cm",width=25,height=7*length(all.years))
  
  cols <- c(brewer.pal(4, "Set1")[c(3)], brewer.pal(4, "Set2")[c(2,3)],8)
  names(cols) <- c("summer","fall","spring","winter")
  tw.colour <- brewer.pal(8, "Set2")[6]
  
  layout(matrix(c(1,2,2, 3,4,4, 5,6,6, 7,8,8, 9,10,10), length(all.years), 3, byrow = TRUE),
         widths=1, heights=c(1))
  par(oma=c(2,11,3,0))
  for(yr in all.years){
    
    glsb <- gls[gls$year2==yr,]
    years <- unique(glsb$year)
    gls2 <- glsb[,c("datetime","Light")]
    colnames(gls2) <- c("Date","Light")
    
    p44 <- p4[p4$year2==yr,]
    
    argos3 <- arg2[arg2$year2==yr & arg2$gls.id == GLS,]
    
    if(GLS=="Q178") ext <- extent(-48,-40,59,63)
    if(GLS=="Q181") ext <- extent(-36,-29,67.5,69.5)
    
    ###### plot map
    par(mar=c(0,0,0,0))
    plot(st_crop(land2$geometry, ext), col=grey(0.85), border=grey(0.85), xlim=c(ext[1],ext[2]),ylim=c(ext[3],ext[4]), lwd=0.5) #ylim=c(59,61), 
    lines(argos3$Longitude,argos3$Latitude,col=grey(0.1))
    
    for(c in 1:length(cols)) lines(argos3$Longitude[argos3$season==names(cols[c])],argos3$Latitude[argos3$season==names(cols[c])],col=cols[c])
    for(c in 1:length(cols)) lines(rep(median(p44$lon[p44$season==names(cols[c])], na.rm=T),2),c(ext[3],ext[4]), col=cols[c],lty=3)
    for(c in 1:length(cols)) lines(c(ext[1],ext[2]),rep(median(p44$lat[p44$season==names(cols[c])], na.rm=T),2), col=cols[c],lty=3)
    
    if(GLS == "Q181") points(capture$lon[1], capture$lat[1], pch=23, cex=2,bg=grey(1))
    if(yr == as.numeric(strftime(capture$date[1], "%Y"))) points(capture$lon[1], capture$lat[1], pch=23, cex=2,bg=grey(1))
    
    for(c in 1:length(cols)) points(median(p44$lon[p44$season==names(cols[c])],na.rm=T),   median(p44$lat[p44$season==names(cols[c])],na.rm=T), pch=21, bg=cols[c], cex=2)
    
    if(yr == as.numeric(strftime(capture$date[2], "%Y"))) points(capture$lon[2], capture$lat[2], pch=23, cex=2,bg=grey(1))
    
    sess <- c("summer","fall","spring","winter")
    if(GLS=='Q178') sess <- c("summer","fall","spring")
    
    for(ses in sess){
      if(length(p44$lon[p44$season==ses & !is.na(p44$lon)])>1){
        den.lon <- density(p44$lon[p44$season==ses], na.rm=T, from = ext[1], to = ext[2])
        den.lon$x <- c(ext[1],den.lon$x,ext[2])
        den.lon$y <- c(0,den.lon$y,0)
        polygon(den.lon$x,den.lon$y/max(den.lon$y)*0.5+ext[3],border="transparent",col=alpha(cols[ses],0.2))
      }
    }
    
    for(ses in sess){
      if(length(p44$lon[p44$season==ses & !is.na(p44$lon)])>1){
        den.lat <- density(p44$lat[p44$season==ses], na.rm=T, from = ext[3], to = ext[4])
        den.lat$x <- c(ext[3],den.lat$x,ext[4])
        den.lat$y <- c(0,den.lat$y,0)
        polygon(den.lat$y/max(den.lat$y)+ext[1],den.lat$x,border="transparent",col=alpha(cols[ses],0.2))
      }
    }
    
    yr2 <- yr
    axis(1,at=seq(ext[1],ext[2],2),labels = paste0(seq(ext[1],ext[2],2),'°'),pos=ext[3])
    axis(2,at=seq(ext[3],ext[4],1),labels = paste0(seq(ext[3],ext[4],1),'°'),pos=ext[1],las=2)
    
    if(yr==min(all.years)) mtext(paste0(years[1]-1,"/", years[1]), side = 2, line = 2.5, las=1, cex=0.8)
    if(yr>min(all.years)) mtext(paste0(years[1],"/", years[2]), side = 2, line = 2.5, las=1,cex=0.8)
    if(yr==min(all.years)) mtext("Space", side = 3, line = 0.5, cex=1.3)
    
    
    ##### plot timeline
    par(mar=c(0.5,2,0,0.5))
    offset <- 3 # adjusts the y-axis to put night (dark shades) in the middle
    if(yr<max(all.years)) lightImage( tagdata = gls2, offset = offset, zlim = c(0, 20), ylab = "hour [UTC]", main = "", 
                                      xlim =as.POSIXct(c(paste0(yr-1,"-07-31"),paste0(yr,"-08-01"))),xaxt="n")
    if(yr==max(all.years)) lightImage( tagdata = gls2, offset = offset, zlim = c(0, 20), ylab = "hour [UTC]", main = "", 
                                       xlim =as.POSIXct(c(paste0(yr-1,"-07-31"),paste0(yr,"-08-01"))))
    tsimageDeploymentLines(gls2$Date, lon = capture$lon, lat = capture$lat, offset = offset, lwd = 3, col = adjustcolor(tw.colour, alpha.f = 0.8))
    
    years <- unique(glsb$year)
    if(length(years)==1) years <- rep(years,2)
    
    # equinox
    fall.doy <- quantile(p4$doy[is.na(p4$lat) & p4$season=="fall"], c(0.05,0.95),na.rm=T)
    tm.calib <- as.POSIXct(as.Date(fall.doy, paste0(years[1],"-01-01")), tz = "UTC")
    rect(tm.calib[1],-10,tm.calib[2],40, border="transparent",density=30,col=alpha(cols[2],0.5))
    
    spring.doy <- quantile(p4$doy[is.na(p4$lat) & p4$season=="spring"], c(0.05,0.95),na.rm=T)
    tm.calib <- as.POSIXct(as.Date(spring.doy, paste0(years[2],"-01-01")), tz = "UTC")
    if(GLS == "Q178") tm.calib <- tm.calib[tm.calib>= as.POSIXct(capture$date[1], tz = "UTC")] 
    rect(tm.calib[1],-10,tm.calib[2],40, border="transparent",density=30,col=alpha(cols[3],0.5))
    
    # ticks for position estimations
    for(c in 1:length(cols)) rug(side=1, x=p44$datetime[p44$season==names(cols[c])],col=cols[c],lwd=1.5)
    
    # capture timing
    if(GLS=="Q181") {arrows(as.POSIXct("2016-09-01"),12+offset,as.POSIXct("2016-08-02"),12+offset, lwd=3,length=0.1)
      points(as.POSIXct("2016-09-01"), 12 + offset, pch=23, bg=grey(1), cex=4.5)
      text(as.POSIXct("2016-09-01"), 12 + offset, lwd=3,lty=1, "C", col=1, cex=1.9)
    }
    
    gls.dates <- as.POSIXct(capture$date)
    gls.dates <- gls.dates[!is.na(gls.dates)]
    
    points(gls.dates, rep(12 + offset, length(gls.dates)), pch=23, bg=grey(1), cex=4.5)
    text(gls.dates, 12 + offset, lwd=3,lty=1, "C", col=1, cex=1.9)
    # 
    if(yr==min(all.years)) mtext("Time", side = 3, line = 0.5, cex=1.3)
    
    
  }
  
  dev.off()
  
}
